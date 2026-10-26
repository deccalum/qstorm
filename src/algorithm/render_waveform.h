#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <span>

#include "algorithm/axis_quantizer.h"
#include "algorithm/math_scalar.h"

namespace qstorm::algorithm {
inline constexpr size_t wave_plane = 3;

template <typename X, typename Y, typename Z>
Spectrum<X, Y, Z>
sample_peak(const Spectrum<X, Y, Z>& env, vvec<Y>& magnitude, size_t plane = wave_plane) {
  Spectrum<X, Y, Z> triple{{}, env.qx, env.qy, env.qz};
  if (plane == wave_plane) magnitude.assign(env.ch.size(), vec<Y>());
  for (size_t c = 0; c < env.ch.size(); ++c) {
    triple.ch.push_back(env.ch[c]);
    triple.ch.push_back(env.ch[c]);
    if (plane < wave_plane) continue;
    magnitude[c].reserve(env.qx.n);
    for (Y sample : ch(env, c)) magnitude[c].push_back(std::abs(sample));
    triple.ch.push_back(magnitude[c].data());
  }
  return triple;
}

template <typename Y>
Y rms(std::span<const Y> level) {
  Y pwr = 0;
  for (Y each : level) pwr += each * each;
  return std::sqrt(pwr / static_cast<Y>(level.size()));
}

template <typename X, typename Y, typename Z>
void span_peak(const Spectrum<X, Y, Z>& env, vvec<Y>& ch, size_t plane = wave_plane) {
  for (size_t c = 0; c < ch.size(); c += plane) {
    ch[c].push_back(std::ranges::min(algorithm::ch(env, c)));
    ch[c + 1].push_back(std::ranges::max(algorithm::ch(env, c + 1)));
    if (plane == wave_plane) ch[c + 2].push_back(rms<Y>(algorithm::ch(env, c + 2)));
  }
}

template <typename X, typename Y, typename Z>
Spectrum<X, Y, Z> pair_peak(const Spectrum<X, Y, Z>& env, vvec<Y>& ch, size_t plane = wave_plane) {
  const X      dx = 2 * env.x0().time_fraction();
  const size_t nx = (env.qx.n + 1) / 2;
  ch.assign(env.ch.size(), vec<Y>());
  for (vec<Y>& plane : ch) plane.reserve(nx);
  for (ms<X> x = env.x0(); x < env.x1();) {
    ms<X> next = x;
    ++next;
    ++next;
    span_peak(env.window(x, std::min<X>(next, env.x1())), ch, plane);
    x = next;
  }
  const Qtz<X> qx = Qtz<X>::stepped(ch.front().size(), dx, env.x0());
  return as_spectrum(ch, qx, env.qy, env.qz);
}

template <typename X, typename Y, typename Z>
zvec<X, Y, Z>
waveform_level(const Spectrum<X, Y, Z>& env, vvvec<Y>& ch, size_t plane = wave_plane) {
  ch.emplace_back();
  zvec<X, Y, Z> level{sample_peak(env, ch.back(), plane)};
  while (level.back().qx.n > 1) {
    ch.emplace_back();
    level.push_back(pair_peak(level.back(), ch.back(), plane));
  }
  return level;
}

template <typename X, typename Y, typename Z>
const Spectrum<X, Y, Z>& level_for(const zvec<X, Y, Z>& level, X x0, X x1, size_t nx) {
  for (size_t i = level.size(); i-- > 0;)
    if (level[i].window(x0, x1).qx.n >= nx) return level[i];
  return level.front();
}

template <typename X, typename Y, typename Z>
Spectrum<X, Y, Z> column_peak(const Spectrum<X, Y, Z>& env, X x0, X x1, size_t nx, vvec<Y>& ch,
                              size_t plane = wave_plane) {
  const X      dx   = ratio<X>(x1 - x0, nx);
  const Qtz<X> qx   = Qtz<X>::stepped(nx, dx, x0);
  const X      xn   = grid(x0, qx, nx);
  ms<X>        xe   = env.x1();
  --xe;
  ch.assign(env.ch.size(), vec<Y>());
  for (ms<X> x{x0, &qx}; x < xn; ++x) {
    ms<X> next = x;
    ++next;
    const Spectrum<X, Y, Z> slice = env.window(x, next);
    if (slice.qx.n)
      span_peak(slice, ch, plane);
    else
      span_peak(env.frame(std::min<X>(x, xe)), ch, plane);
  }
  return as_spectrum(ch, qx, env.qy, env.qz);
}

template <typename X, typename Y, typename Z>
Spectrum<X, Y, Z> view_wave(const zvec<X, Y, Z>& level, X x0, X x1, size_t nx, vvec<Y>& ch,
                            size_t plane = wave_plane) {
  const Spectrum<X, Y, Z>& sample = level.front();
  if (sample.window(x0, x1).qx.n >= nx)
    return column_peak(level_for(level, x0, x1, nx), x0, x1, nx, ch, plane);
  ms<X> lo{x0, &sample.qx};
  ms<X> hi{x1, &sample.qx};
  --lo;
  ++hi;
  return sample.window(std::max<X>(lo, sample.x0()), std::min<X>(hi, sample.x1()));
}
}
