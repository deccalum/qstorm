#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <expected>
#include <filesystem>
#include <fstream>
#include <meta>
#include <qr/qr.hpp>
#include <span>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "algorithm/dsp_pyramid.h"
#include "algorithm/dsp_spectrum.h"
#include "atom/atom_constant.h"
#include "io/bridge.h"
#include "io/interpret_file.h"
#include "io/database.h"
#include "io/player.h"

namespace qstorm::cli {

template <typename T>
using Qtz = algorithm::Quantization<T>;
template <typename X, typename Y = X, typename Z = X>
using Spectrum = algorithm::Spectrogram<X, Y, Z>;
using Level = algorithm::Spectrogram<double, float, double>;

inline void decode(int64_t source_id) {
  const auto [id, audio] = io::table::source(source_id);
  io::write_bytes(io::path(id, io::env::ext::flac), audio);
  io::decode_flac(audio);
  double input = 0.0;
  for (io::env::pcm::shifted sample : io::env::pcm::interleaved)
    input = std::max(input, std::fabs(static_cast<double>(sample)));
  const size_t frame = io::env::pcm::interleaved.size() / io::env::pcm::count;

  io::env::pcm::pcm_db<double>.assign(io::env::pcm::count, std::vector<double>());
  for (size_t i = 0; i < frame; ++i)
    for (size_t j = 0; j < io::env::pcm::count; ++j) {
      const double amp = static_cast<double>(io::env::pcm::interleaved[i*io::env::pcm::count +j]);
      io::env::pcm::pcm_db<double>[j].push_back(algorithm::db<double>::ratio(std::fabs(amp) / input));
    }

  const Qtz<double> x{frame, io::env::pcm::ms, 0.0};
  const Qtz<double> y = algorithm::y_db(io::env::pcm::pcm_db<double>, algorithm::depth_floor<double>(io::env::bit_depth::decoded));
  const Qtz<double> z{1, algorithm::ratio<double>(2, io::env::pcm::sample), 0.0};
  const Spectrum<double> signal = algorithm::as_spectrogram(io::env::pcm::pcm_db<double>, x, z, y);
  std::ofstream file(io::path(id, io::env::ext::bin), std::ios::binary);
  io::write_spectrogram(file, signal);
}

inline void run_atomextract(int64_t id, int64_t dna_id, double start, double end) {
  std::vector<std::vector<double>> channel;
  const atom::Domain::Factor dna = io::table::factor(dna_id);
  const Spectrum<double>  signal = io::map_spectrogram<double, double, double>(io::path(id, io::env::ext::bin)).front();
  const Spectrum<double>    mono = algorithm::make_mono(signal, channel);
  const Spectrum<double>  window = mono.window(start, end);
  const auto               found = atom::extract_root_atom(window, dna);
  io::table::record(id, dna, window, found);
}

inline void ingest(const char* directory) {
  const std::filesystem::path folder(directory);
  const std::string origin = folder.parent_path().filename().string();
  const std::string domain = folder.filename().string();
  std::vector<std::filesystem::path> file(std::filesystem::directory_iterator(folder), {});
  std::ranges::sort(file);
  for (const std::filesystem::path& each : file) {
    if (!io::audio(each)) continue;
    const std::string stem = each.stem().string();
    const size_t   first = stem.find(io::env::described);
    const size_t  second = stem.find(io::env::described, first + io::env::described.size());
    const bool described = second != std::string::npos;
    const int64_t     id = std::stoll(stem.substr(0, first));
    const std::expected<std::pair<int, int>, int> probed = io::probe(each);
    if (!probed) {
      std::fprintf(stderr, "%s: %s exit %d\n", each.c_str(), io::env::value.ffprobe.c_str(), probed.error());
      return;
    }
    const auto [hz, channels] = *probed;
    const std::expected<void, int> transcoded = io::transcode(each);
    if (!transcoded) {
      std::fprintf(stderr, "%s: %s exit %d\n", each.c_str(), io::env::value.ffmpeg.c_str(), transcoded.error());
      return;
    }
    const std::vector<std::byte> audio = io::encode_flac(io::env::pcm::interleaved, static_cast<size_t>(channels), hz);

    io::table::source_info row{};
    row.id = id;
    row.domain = domain;
    row.origin = origin;
    row.hz = hz;
    row.channels = channels;
    row.url = std::string(io::env::freesound) + std::to_string(id) + "/";
    if (described) {
      row.recorded_by = stem.substr(first + io::env::described.size(), second - first - io::env::described.size());
      row.title = stem.substr(second + io::env::described.size());
    }
    io::table::store(row, audio);
    std::printf("%lld -> %zu bytes flac, %d ch @ %d Hz\n", static_cast<long long>(id), audio.size(), channels, hz);
  }
}

inline void play(int64_t id, double start, double end) {
  if (!std::filesystem::exists(io::path(id, io::env::ext::bin))) decode(id);
  const Spectrum<double> db = io::map_spectrogram<double, double, double>(io::path(id, io::env::ext::bin)).front();
  const Qtz<float> y = db.y_db.as<float>();
  const io::player::Sink<io::env::pcm::shifted> sink{io::env::playback_passthrough};
  io::decode_flac_range(io::path(id, io::env::ext::flac), start, end);
  using PlaybackControl = io::player::PlaybackControl<io::player::Sink<io::env::pcm::shifted>, float>;
  const size_t full = io::env::gain_step - 1;
  PlaybackControl control{sink, y, full, std::vector<size_t>(io::env::pcm::count, full), false};
  const bool finished = io::player::play(sink, control);
  if (!finished) sink.abort();
  io::player::close(sink);
}

inline void pcm(int64_t id, double start, double end) {
  io::decode_flac_range(io::path(id, io::env::ext::flac), start, end);
  const size_t frame = io::env::pcm::interleaved.size() / io::env::pcm::count;
  const Qtz<double> x_ms{frame, io::env::pcm::ms, start};
  std::fwrite(&x_ms, sizeof x_ms, 1, stdout);
  std::fwrite(&io::env::pcm::count, sizeof io::env::pcm::count, 1, stdout);
  const std::vector<float> amplitude = algorithm::amplitude<float>(io::env::pcm::interleaved);
  std::fwrite(amplitude.data(), sizeof(float), amplitude.size(), stdout);
}

inline void pyramid(int64_t id) {
  if (!std::filesystem::exists(io::path(id, io::env::ext::bin))) decode(id);
  const Spectrum<double> db = io::map_spectrogram<double, double, double>(io::path(id, io::env::ext::bin)).front();
  const double start = db.begin();
  const double   end = db.end();
  io::decode_flac_range(io::path(id, io::env::ext::flac), start, end);
  const size_t frame = io::env::pcm::interleaved.size() / io::env::pcm::count;
  const size_t count = io::env::pcm::count;

  std::vector<float> amplitude = algorithm::amplitude<float>(io::env::pcm::interleaved);
  std::vector<io::env::pcm::shifted>().swap(io::env::pcm::interleaved);
  std::vector<std::vector<float>> sample = algorithm::deinterleave<float>(amplitude, count);
  std::vector<float>().swap(amplitude);
  const Qtz<double> x{frame, io::env::pcm::ms, start};
  const Qtz<float> y = db.y_db.as<float>();
  const Level signal = algorithm::as_spectrogram(sample, x, db.z_hz, y);
  std::vector<std::vector<std::vector<float>>> channel;
  const std::vector<Level> level = algorithm::pyramid(signal, channel, 2);

  std::ofstream file(io::path(id, io::env::ext::pyramid), std::ios::binary);
  for (const Level& l : level) io::write_spectrogram(file, l);
}

inline void clear_cache(int64_t id) {
  std::filesystem::remove(io::path(id, io::env::ext::pyramid));
}

inline void view(int64_t id, double t0, double t1, int64_t width) {
  std::vector<std::vector<float>> channel;
  if (!std::filesystem::exists(io::path(id, io::env::ext::pyramid))) pyramid(id);
  const std::vector<Level> level = io::map_spectrogram<double, float, double>(io::path(id, io::env::ext::pyramid));
  const Level&   sample = level.front();
  const double    start = sample.begin();
  const double      end = sample.end();
  const double     from = std::max(t0, start);
  const double    until = std::min(t1, end);
  const size_t   column = static_cast<size_t>(width);
  const Level     shown = algorithm::view(level, from, until, column, channel, 2);
  const double frame_ms = sample.begin().time_fraction();
  const size_t    frame = sample.window(from, until).x_ms.extent;
  const double      fpp = algorithm::ratio<double>(frame, column);
  const std::array<double, 5> header{start, end, frame_ms, sample.y_db.offset, fpp};
  std::vector<double> time;
  for (algorithm::ms<double> x = shown.begin(); x < shown.end(); ++x) time.push_back(x);
  const std::array<uint64_t, 2> shape{shown.channel.size() / 2, time.size()};
  std::fwrite(header.data(), sizeof(double), header.size(), stdout);
  std::fwrite(shape.data(), sizeof(uint64_t), shape.size(), stdout);
  std::fwrite(time.data(), sizeof(double), time.size(), stdout);
  for (const float* row : shown.channel) std::fwrite(row, sizeof(float), shown.x_ms.extent, stdout);
}

inline void edit(int64_t run_id, int64_t index, const char* column, double value) {
  io::table::edit(run_id, index, column, value);
}

template <typename T>
void read_rows(const char* parent) {
  const std::vector<T> rows = parent
                            ? io::table::select<T>(std::strtoll(parent, nullptr, 10))
                            : io::table::select<T>();
  qr::write_rows(rows, stdout);
}

inline void read(const char* table, const char* parent) {
  template for (constexpr std::meta::info type : qr::detail::members_of<^^io::table>()) {
    if constexpr (std::meta::is_type(type) && std::meta::is_class_type(type)) {
      if (std::meta::identifier_of(type) == table) read_rows<typename [:type:]>(parent);
    }
  }
}

inline int studio() { return io::studio(); }

}
