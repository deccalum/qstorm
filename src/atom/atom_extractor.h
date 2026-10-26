#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <utility>
#include <variant>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "algorithm/dsp_envelope.h"
#include "algorithm/dsp_spectrum.h"
#include "algorithm/math_scalar.h"
#include "algorithm/sim_matrix.h"
#include "atom/atom_constant.h"
#include "atom/atom_thunder.h"

namespace qstorm::atom {
using Domain = std::variant<thunder::Root>;
using Factor = std::variant<thunder::Factor>;

struct Atom {
  std::optional<Domain> domain;
  vec<Component>        component;
  vec<Atom>             atom;
};

template <size_t I = 0>
vec<Atom> extract_root_atom(const Spectrum<double>& window, const Factor& factor) {
  if constexpr (I + 1 < std::variant_size_v<Factor>) {
    if (factor.index() != I) return extract_root_atom<I + 1>(window, factor);
  }
  vec<ADSR>                                        adsr;
  const vec<std::variant_alternative_t<I, Domain>> dom = root(window, std::get<I>(factor), adsr);
  vec<Atom>                                        atom;
  for (size_t j = 0; j < dom.size(); ++j) atom.push_back(Atom{dom[j], {adsr[j]}, {}});
  return atom;
}

template <typename X, typename Y, typename Z>
std::optional<Spectral> spectral_shape(const zvec<X, Y, Z>& band) {
  const std::optional<double> lvl = algorithm::spectral_level(band);
  std::optional<Spectral>     spectral;
  if (lvl) {
    const double          center = algorithm::spectral_centroid(band);
    const std::pair<Z, Z> edge   = algorithm::spectral_edge(band, *lvl);

    spectral = Spectral{
      .peak     = algorithm::spectral_peak(band),
      .center   = center,
      .spread   = algorithm::spectral_spread(band, center),
      .low      = edge.first,
      .high     = edge.second,
      .tilt     = algorithm::spectral_tilt(band),
      .flatness = algorithm::spectral_flatness(band),
    };
  }
  return spectral;
}

template <typename X, typename Y, typename Z>
std::optional<ADSR> adsr(const Spectrum<X, Y, Z>& env, double energy) {
  const std::optional<double> lvl = algorithm::envelope_level(env);
  std::optional<ADSR>         adsr;
  if (lvl) {
    const Y*     rms     = env.ch.front();
    const size_t crest   = algorithm::argmax(std::span<const Y>(rms, env.qx.n));
    const double y0      = *lvl - rms[crest];
    const double peak    = algorithm::grid(env.qx.pos, env.qx, static_cast<int64_t>(crest));
    const double attack  = algorithm::attack_time(env, y0, 0.0);
    const double release = algorithm::release_time(env, y0, 0.0);
    const double sustain = algorithm::sustain_level(env, y0);
    const double decay   = algorithm::release_time(env, sustain, 0.0);
    const double center  = algorithm::temporal_centroid(env);

    adsr = ADSR{
      .x0      = { peak - attack, &env.qx},
      .attack  = {        attack, &env.qx},
      .peak    = {          peak, &env.qx},
      .decay   = ms<double>{         decay, &env.qx},
      .sustain = sustain,
      .release = {       release, &env.qx},
      .x1      = {peak + release, &env.qx},
      .energy  = energy,
      .slope   = algorithm::decay_rate(env, y0),
      .center  = {        center, &env.qx}
    };
  }
  return adsr;
}

template <typename X, typename Y, typename Z>
Spectrum<X, Y, Z> activation_envelope(const zvec<X, Y, Z>& band, const vec<>& spectrum,
                                      const vec<>& activation, vvec<Y>& ch) {
  const double z1 = *std::ranges::max_element(spectrum);
  const Qtz<Z> single{1, 1, 0};
  ch.assign(1, {});
  for (double h : activation) ch[0].push_back(algorithm::db<Y>::ratio(std::sqrt(h * z1)));
  const Qtz<Y> qy = algorithm::ydb(ch, band.front().qy);
  return algorithm::as_spectrum(ch, band.front().qx, qy, single);
}

template <typename X, typename Y, typename Z, typename R, typename K, typename C>
Atom through_nmf(const zvec<X, Y, Z>& band, const matrix::Matrix<R, K>& spectrum,
                 const matrix::Matrix<K, C>& activation, const matrix::Relevance& relevance) {
  const size_t      row     = spectrum.row.n;
  const size_t      rank    = spectrum.col.n;
  const size_t      nx      = activation.col.n;
  const double      row_pho = spectrum.row.rho;
  const double      row_pos = spectrum.row.pos;
  const double      col_pho = activation.col.rho;
  const double      col_pos = activation.col.pos;
  const Qtz<double> x{nx, col_pho, col_pos};
  const Qtz<double> z{row, row_pho, row_pos};
  const vec<Qtz<Z>> layout = algorithm::band_axis(band);
  const Qtz<double> octave = algorithm::octave_axis(layout);
  Atom              atom;

  for (size_t i = 0; i < rank; ++i) {
    Raw zraw;
    zraw.val.resize(row);
    for (size_t j = 0; j < row; ++j) zraw.val[j] = spectrum.cell[j * rank + i];
    zraw.axis = z;

    Raw araw;
    araw.val.resize(nx);
    for (size_t j = 0; j < nx; ++j) araw.val[j] = activation.cell[i * nx + j];
    araw.axis         = x;
    const double act1 = *std::ranges::max_element(araw.val);

    Raw draw;
    draw.val  = algorithm::octave_density(layout, octave, zraw.val);
    draw.axis = octave;

    vvvec<Y>      col_lvl(band.size());
    zvec<X, Y, Z> col;
    size_t        i0 = 0;
    for (size_t b = 0; b < band.size(); ++b) {
      vec<Y> level;
      for (size_t j = i0; j < i0 + band[b].qz.n; ++j)
        level.push_back(algorithm::db<Y>::ratio(std::sqrt(zraw.val[j] * act1)));
      col_lvl[b].push_back(std::move(level));
      const Qtz<X> once{1, band[b].qx.rho, band[b].qx.pos};
      const Qtz<Y> qy = algorithm::ydb(col_lvl[b], band[b].qy);
      col.push_back(algorithm::as_spectrum(col_lvl[b], once, qy, band[b].qz));
      i0 += band[b].qz.n;
    }

    vvec<Y>                       ch;
    const Spectrum<X, Y, Z>       aenv  = activation_envelope(band, zraw.val, araw.val, ch);
    const std::optional<Spectral> shape = spectral_shape(col);
    const std::optional<ADSR>     env   = adsr(aenv, relevance.share[i]);

    Atom a;
    a.component.push_back(std::move(zraw));
    a.component.push_back(std::move(araw));
    a.component.push_back(std::move(draw));
    if (shape) a.component.push_back(*shape);
    if (env) a.component.push_back(*env);
    atom.atom.push_back(std::move(a));
  }
  return atom;
}
}
