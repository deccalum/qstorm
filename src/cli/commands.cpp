#include "commands.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <expected>
#include <filesystem>
#include <limits>
#include <optional>
#include <print>
#include <span>
#include <string>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "algorithm/dsp_spectrum.h"
#include "algorithm/render_waveform.h"
#include "atom/atom_constant.h"
#include "atom/atom_extractor.h"
#include "io/bridge.h"
#include "io/database.h"
#include "io/interpret_file.h"
#include "io/player.h"

namespace qstorm::cli {

// ********************************************************
// *** [ LIBRARY ] ****************************************
// ********************************************************
namespace library {
void ingest_source(const char* directory) {
  const std::filesystem::path folder(directory);
  const std::string           origin = folder.parent_path().filename().string();
  const std::string           domain = folder.filename().string();
  vec<std::filesystem::path>  file(std::filesystem::directory_iterator(folder), {});
  std::ranges::sort(file);
  for (const std::filesystem::path& each : file) {
    if (!io::audio(each)) continue;
    const std::string stem = each.stem().string();
    const int64_t     id   = std::stoll(stem.substr(0, stem.find(env::described)));

    const std::expected<io::Probe, int> probed = io::probe(each);
    if (!probed) {
      std::fprintf(stderr, "%s: %s exit %d\n", each.c_str(),
                   env::require(env::var::ffprobe).c_str(), probed.error());
      return;
    }
    const std::expected<pcmvec, int> transcoded = io::transcode(each);
    if (!transcoded) {
      std::fprintf(stderr, "%s: %s exit %d\n", each.c_str(), env::require(env::var::ffmpeg).c_str(),
                   transcoded.error());
      return;
    }
    const vec<std::byte> audio = io::flac_encode(*transcoded, probed->channels, probed->hz);

    table::write_audio(id, domain, origin, *probed, env::ingest_bit, stem, audio);
    decode_source(id);
    std::println("{} -> {} bytes flac, {} ch @ {} kHz", id, audio.size(), probed->channels,
                 probed->hz);
  }
}

void decode_source(int64_t source_id) {
  const table::source source = table::find<table::source>(source_id);
  const pcmvec        xc     = io::flac_decode(*source.audio, source.bit_depth);
  const size_t        nc     = source.channels;
  double              input  = 0.0;
  for (env::pcm::shifted sample : xc) input = std::max(input, std::fabs(sample));
  const size_t nx = xc.size() / nc;

  vvec ch(nc);
  for (size_t i = 0; i < nx; ++i)
    for (size_t j = 0; j < nc; ++j) {
      const double amp = xc[i * nc + j];
      ch[j].push_back(algorithm::db<double>::ratio(std::fabs(amp) / input));
    }

  const Qtz<double> x{nx, source.hz, 0.0};
  const Qtz<double> y = algorithm::ydb(ch, algorithm::depth_floor<double>(source.bit_depth));
  const Qtz<double> z{1, algorithm::ratio<double>(2, source.hz), 0.0};

  const Spectrum<double> sample = algorithm::as_spectrum(ch, x, y, z);
  vvec                   log_c;
  const Spectrum<double> log_mono = algorithm::log_mono(sample, log_c);
  vvec                   lin_c;
  const Spectrum<double> lin_mono = algorithm::lin_mono(sample, lin_c);
  vvvec<float>           lvl_ch;

  vvec<float>      amp;
  const Level      envelope = algorithm::xc_spectrum(xc, nc, source.hz, 0.0, y.as<float>(), z, amp);
  const vec<Level> waveform = algorithm::waveform_level(envelope, lvl_ch, 2);

  table::write_decoded(source_id, sample, log_mono, lin_mono, waveform);
}

vvec<std::byte>
encode_component(const atom::Atom& found, const Spectrum<double>& env, const zvec<double>& stft,
                 const table::source& source, const matrix::Matrix<int, double>& power) {
  vvec         spectrum;
  vvec         activation;
  vvec<size_t> member;
  for (size_t k = 0; k < found.atom.size(); ++k) {
    spectrum.push_back(std::get<atom::Raw>(found.atom[k].component[0]).val);
    activation.push_back(std::get<atom::Raw>(found.atom[k].component[1]).val);
    member.push_back({k});
  }
  const size_t      nc = source.channels;
  const vec<pcmvec> masked
    = matrix::component_audio<env::pcm::shifted>(env, stft, spectrum, activation, power, member);
  vvec<std::byte> audio;
  for (const pcmvec& each : masked) audio.push_back(io::flac_encode(each, nc, source.hz));
  return audio;
}
}

// ********************************************************
// *** [ EXTRACT ] ****************************************
// ********************************************************

namespace extract {
void root_atom(int64_t id, int64_t dna_id) {
  const table::decoded   dcode  = table::find<table::decoded>(id);
  const table::dna       dna    = table::find<table::dna>(dna_id);
  const Spectrum<double> window = *dcode.log_mono;
  const vec<atom::Atom>  atom   = atom::extract_root_atom(window, dna.factor);
  const table::source    source = table::find<table::source>(id);
  const size_t           nc     = source.channels;
  const pcmvec           xc     = io::flac_decode(*source.audio, source.bit_depth);
  vvec<std::byte>        audio;
  for (const atom::Atom& each : atom) {
    const atom::ADSR& adsr = std::get<atom::ADSR>(each.component.front());
    const pcmspan     part = algorithm::xc_window(xc, nc, source.hz, adsr.x0, adsr.x1);
    audio.push_back(io::flac_encode(part, nc, source.hz));
  }
  table::write_root_extraction(dcode, dna, window, atom, audio);
}

void atomize(int64_t root_id, int64_t index_id, int64_t ard_id, int ch) {
  using namespace qstorm::matrix;
  const table::root    root = table::latest<table::root>(root_id, index_id);
  const table::stft    stft = table::find<table::stft>(root.extraction, root.index, root.revision);
  const table::ard     ard  = table::find<table::ard>(ard_id);
  const table::source    source = table::follow<table::source>(root);
  const table::decoded   dcode  = table::follow<table::decoded>(root);
  const zvec<double>     zc     = algorithm::select_ch(*stft.stft, ch);
  const pcmvec           xc     = io::flac_decode(*root.audio, source.bit_depth);
  vvec                   amp;
  const Spectrum<double> env  = algorithm::xc_spectrum(xc, source.channels, source.hz, root.adsr.x0,
                                                       dcode.sample.qy, dcode.sample.qz, amp);
  const Matrix<int, double> V = power_matrix(zc);
  auto [W, H]                 = init_factorization(V, ard.factor);
  const vec             wgt   = band_row_weight(zc);
  const Relevance       rel   = non_neg_matrix_factorization(V, W, H, ard.factor, wgt);
  const vec             resid = band_residual(zc, rel.residual);
  const atom::Atom      atom  = atom::through_nmf(zc, W, H, rel);
  const vvec<std::byte> audio = library::encode_component(atom, env, *stft.stft, source, V);
  table::write_atomization(root, ard, ch, rel, resid, atom, audio);
}

void group_component(int64_t atomization_id) {
  const table::atomization    a    = table::find<table::atomization>(atomization_id);
  const vec<table::component> live = table::select_live_component(atomization_id);
  const auto [z, activation]       = table::split_spectrum_activation(live);
  const vvec cos                   = matrix::cos_matrix(z);
  const vvec corr                  = matrix::corr_matrix(activation);
  vvec       density;
  for (const table::component& each : live) density.push_back(each.density);
  const vvec<bool>  related  = matrix::shift_related(activation);
  const vec<size_t> spectrum = matrix::complete_linkage(cos, matrix::overlap_related(density));
  const vec<size_t> event    = matrix::complete_linkage(corr, related);
  table::write_group(a, live, cos, corr, related, spectrum, event);
}

void link_component(const char* domain, const char* source, const char* atomization) {
  const table::linking             scope  = table::link_scope(domain, source, atomization);
  const vec<table::atomization>    window = table::select_link_window(scope);
  vec<std::pair<int64_t, int64_t>> member;
  vec<Qtz<double>>                 axis;
  vvec                             density;
  for (const table::atomization& a : window)
    for (const table::component& each : table::select_live_component(a.id)) {
      member.emplace_back(a.id, each.index);
      axis.push_back(each.octave);
      density.push_back(each.density);
    }

  const Qtz<double> octave  = algorithm::octave_union(axis);
  const vvec        share   = algorithm::octave_resample(axis, density, octave);
  const vvec        cos     = matrix::cos_matrix(share);
  const vvec<bool>  related = matrix::overlap_related(share);
  const vec<size_t> set     = matrix::complete_linkage(cos, related);
  table::write_link(scope, member, set);
}

void subtract_reconstruction(int64_t atomization_id) {
  const table::atomization    a    = table::find<table::atomization>(atomization_id);
  const table::stft           stft = table::find<table::stft>(a.extraction, a.index, a.revision);
  const vec<table::component> part = table::select<table::component>(atomization_id);
  const zvec<double>          zc   = algorithm::select_ch(*stft.stft, a.channel);
  const auto [z, activation]       = table::split_spectrum_activation(part);
  vvvec ch;

  const zvec<double> layer = matrix::residual_layer(zc, z, activation, ch);
  table::write_residual(a, layer);
}
}

// ********************************************************
// *** [ RENDER ] *****************************************
// ********************************************************
namespace render {
void play_component_set(int64_t atomization_id, const char* component) {
  const table::atomization    a      = table::find<table::atomization>(atomization_id);
  const table::root           root   = table::find<table::root>(a.extraction, a.index, a.revision);
  const table::stft           stft   = table::find<table::stft>(a.extraction, a.index, a.revision);
  const table::source         source = table::follow<table::source>(root);
  const table::decoded        dcode  = table::follow<table::decoded>(root);
  const vec<table::component> part   = table::select<table::component>(atomization_id);
  const vec<int64_t>          comp   = table::comma_number(component);
  const auto [z, activation]         = table::split_spectrum_activation(part);
  vec<size_t> member;
  for (int64_t index : comp)
    member.push_back(std::ranges::find(part, index, &table::component::index) - part.begin());
  const zvec<double>                zc = algorithm::select_ch(*stft.stft, a.channel);
  const matrix::Matrix<int, double> V  = matrix::power_matrix(zc);
  const pcmvec                      xc = io::flac_decode(*root.audio, source.bit_depth);
  vvec                              amp;
  const Spectrum<double> env = algorithm::xc_spectrum(xc, source.channels, source.hz, root.adsr.x0,
                                                      dcode.sample.qy, dcode.sample.qz, amp);
  const vec<pcmvec>      mask
    = matrix::component_audio<env::pcm::shifted>(env, *stft.stft, z, activation, V, {member});
  io::player::stream(mask.front(), source.channels, source.hz, dcode.sample.qy.as<float>());
}

void play_residual(int64_t atomization_id) {
  const table::atomization    a      = table::find<table::atomization>(atomization_id);
  const table::root           root   = table::find<table::root>(a.extraction, a.index, a.revision);
  const table::stft           stft   = table::find<table::stft>(a.extraction, a.index, a.revision);
  const table::source         source = table::follow<table::source>(root);
  const table::decoded        dcode  = table::follow<table::decoded>(root);
  const vec<table::component> part   = table::select<table::component>(atomization_id);
  const auto [z, activation]           = table::split_spectrum_activation(part);
  const zvec<double>                zc = algorithm::select_ch(*stft.stft, a.channel);
  const matrix::Matrix<int, double> V  = matrix::power_matrix(zc);
  const pcmvec                      xc = io::flac_decode(*root.audio, source.bit_depth);
  vvec                              amp;
  const Spectrum<double> env = algorithm::xc_spectrum(xc, source.channels, source.hz, root.adsr.x0,
                                                      dcode.sample.qy, dcode.sample.qz, amp);
  const pcmvec mask = matrix::residual_audio<env::pcm::shifted>(env, *stft.stft, z, activation, V);
  io::player::stream(mask, source.channels, source.hz, dcode.sample.qy.as<float>());
}

void play_audio(int64_t audio, double pos) {
  const qr::oid                                 id{audio};
  const table::source                           source = table::audio_source(id);
  const table::decoded                          dcode  = table::find<table::decoded>(source.id);
  const table::blob<std::span<const std::byte>> flac{id, {}};
  const pcmvec                                  xc = io::flac_decode(*flac, source.bit_depth);
  const pcmspan range = algorithm::xc_window(xc, source.channels, source.hz, pos, k::INF<>);
  io::player::stream(range, source.channels, source.hz, dcode.sample.qy.as<float>());
}

void render_component_list(int64_t atomization_id, int64_t width) {
  const table::atomization    a     = table::find<table::atomization>(atomization_id);
  const table::stft           stft  = table::find<table::stft>(a.extraction, a.index, a.revision);
  const Spectrum<double>&     band  = stft.stft->front();
  const vec<table::component> part  = table::select<table::component>(atomization_id);
  const size_t                col   = width;
  const double                ncomp = part.size();
  const uint64_t              shape = col;
  const std::array<double, 4> header{band.x0(), band.x1(), a.floor, ncomp};

  std::fwrite(header.data(), sizeof(double), header.size(), stdout);
  std::fwrite(&shape, sizeof(uint64_t), 1, stdout);
  const double absent = std::numeric_limits<double>::quiet_NaN();
  for (const table::component& each : part) {
    const double peak  = each.spectral ? each.spectral->peak : absent;
    const double x0    = each.adsr ? each.adsr->x0 : absent;
    const double decay = each.adsr && each.adsr->decay ? *each.adsr->decay : absent;
    const double index = each.index;

    const std::array<double, 5> key{index, each.lambda, peak, x0, decay};
    std::fwrite(key.data(), sizeof(double), key.size(), stdout);
  }

  const double      band_x0 = band.x0();
  const double      band_x1 = band.x1();
  const Qtz<double> row{1, 1.0, 0.0};
  for (const table::component& each : part) {
    vvec                   extreme{each.activation, each.activation};
    vvec                   reduced;
    const Spectrum<double> activation = algorithm::as_spectrum(extreme, band.qx, band.qy, row);
    algorithm::column_peak(activation, band_x0, band_x1, col, reduced, 2);
    const vec<float> peak(reduced[1].begin(), reduced[1].end());
    std::fwrite(peak.data(), sizeof(float), col, stdout);
  }
}

void render_waveform(int64_t id, double x0, double x1, int64_t width) {
  vvec<float>          ch;
  const table::decoded dcode  = table::find<table::decoded>(id);
  const vec<Level>&    lvl    = *dcode.waveform;
  const Level&         sample = lvl.front();
  const double         src_x0 = sample.x0();
  const double         src_x1 = sample.x1();
  const double         x_0    = std::max(x0, src_x0);
  const double         x_1    = std::min(x1, src_x1);
  const size_t         col    = width;
  const Level          viz    = algorithm::view_wave(lvl, x_0, x_1, col, ch, 2);
  const double         dx     = sample.x0().time_fraction();
  const size_t         nx     = sample.window(x_0, x_1).qx.n;
  const double         fpp    = algorithm::ratio<double>(nx, col);

  const std::array<double, 5> header{src_x0, src_x1, dx, sample.qy.pos, fpp};
  vec                         xs;
  for (ms<double> x = viz.x0(); x < viz.x1(); ++x) xs.push_back(x);
  const std::array<uint64_t, 2> shape{viz.ch.size() / 2, xs.size()};
  std::fwrite(header.data(), sizeof(double), header.size(), stdout);
  std::fwrite(shape.data(), sizeof(uint64_t), shape.size(), stdout);
  std::fwrite(xs.data(), sizeof(double), xs.size(), stdout);
  for (const float* row : viz.ch) std::fwrite(row, sizeof(float), viz.qx.n, stdout);
}

void write_stft(int64_t root_id, int64_t index_id) {
  vvvec                  ch;
  const table::root      root   = table::latest<table::root>(root_id, index_id);
  const table::decoded   dcode  = table::follow<table::decoded>(root);
  const table::source    source = table::follow<table::source>(root);
  const pcmvec           xc     = io::flac_decode(*root.audio, source.bit_depth);
  vvec                   amp;
  const Spectrum<double> env  = algorithm::xc_spectrum(xc, source.channels, source.hz, root.adsr.x0,
                                                       dcode.sample.qy, dcode.sample.qz, amp);
  const vec<Qtz<double>> band = algorithm::stft_fft_ms(env);
  const Qtz<double>      col  = algorithm::stft_hop_ms(env, band.back());
  const zvec<double>     stft = algorithm::stft(env, col, band, ch);
  table::write_stft(root, stft);
}

void render_band(const zvec<double>& stored) {
  const zvec<double> band = algorithm::select_ch(stored, 0);
  const Qtz<double>  qx   = band.front().qx;
  const uint64_t     nb   = band.size();
  std::fwrite(&qx, sizeof qx, 1, stdout);
  std::fwrite(&nb, sizeof nb, 1, stdout);
  for (const Spectrum<double>& each : band) std::fwrite(&each.qz, sizeof each.qz, 1, stdout);
  for (const Spectrum<double>& each : band) {
    const size_t     cell = each.qx.n * each.qz.n;
    const double*    row  = each.ch.front();
    const vec<float> db(row, row + cell);
    std::fwrite(db.data(), sizeof(float), cell, stdout);
  }
}

void render_stft(int64_t root_id, int64_t index_id) {
  const table::root root = table::latest<table::root>(root_id, index_id);
  const table::stft stft = table::find<table::stft>(root.extraction, root.index, root.revision);
  render_band(*stft.stft);
}

void render_residual(int64_t atomization_id) {
  const table::residual residual = table::find<table::residual>(atomization_id);
  render_band(*residual.layer);
}

void render_reconstruction(int64_t atomization_id) {
  const table::atomization    a     = table::find<table::atomization>(atomization_id);
  const table::stft           stft  = table::find<table::stft>(a.extraction, a.index, a.revision);
  const vec<table::component> part  = table::select<table::component>(atomization_id);
  const zvec<double>          zc    = algorithm::select_ch(*stft.stft, a.channel);
  const auto [spectrum, activation] = table::split_spectrum_activation(part);
  vvvec ch;

  const zvec<double> layer = matrix::reconstruct_layer(zc, spectrum, activation, ch);
  render_band(layer);
}

void render_component(int64_t atomization_id, int64_t component_index) {
  const table::atomization a    = table::find<table::atomization>(atomization_id);
  const table::stft        stft = table::find<table::stft>(a.extraction, a.index, a.revision);
  const table::component   part = table::find<table::component>(atomization_id, component_index);
  const zvec<double>       zc   = algorithm::select_ch(*stft.stft, a.channel);
  const vvec<>             spectrum{part.spectrum};
  const vvec<>             activation{part.activation};
  vvvec                    ch;

  const zvec<double> layer = matrix::reconstruct_layer(zc, spectrum, activation, ch);
  render_band(layer);
}
}

// ********************************************************
// *** [ QUERIES  ] ***************************************
// ********************************************************
namespace query {
void read_table(const char* table, const char* parent) { table::read(table, parent); }
void edit_root(int64_t run_id, int64_t index, const char* column, double value) {
  table::root          t      = table::edit_root(run_id, index, column, value);
  const table::source  source = table::follow<table::source>(t);
  const size_t         nc     = source.channels;
  const pcmvec         xc     = io::flac_decode(*source.audio, source.bit_depth);
  const pcmspan        part   = algorithm::xc_window(xc, nc, source.hz, t.adsr.x0, t.adsr.x1);
  const vec<std::byte> audio  = io::flac_encode(part, nc, source.hz);
  t.audio                     = table::write_blob(std::span<const std::byte>(audio));
  table::write(t);
}
void restore_root(int64_t run_id, int64_t index) { table::restore_root(run_id, index); }
void edit_dna(int64_t dna_id, const char* column, double value) {
  table::edit_dna(dna_id, column, value);
}
void copy_dna(int64_t dna_id) { std::println("{}", table::copy_dna(dna_id)); }
void clear_cache(const char* table, int64_t parent) { table::clear_cache(table, parent); }
void clear_root_extraction() {
  const size_t removed = table::clear_root_extraction();
  std::println("{} extractions removed", removed);
}
}
}
