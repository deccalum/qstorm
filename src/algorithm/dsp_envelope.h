#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <optional>
#include <span>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "algorithm/dsp_window.h"
#include "algorithm/math_scalar.h"

namespace qstorm::algorithm {
template <typename X, typename Y, typename Z>
double attack_time(const Spectrum<X, Y, Z>& env, double dy0, double dy1) {
  const Qtz<X>& qx   = env.qx;
  const Y*      rms  = env.ch.front();
  const size_t  peak = argmax(std::span<const Y>(rms, qx.n));
  const double  y1   = rms[peak];

  size_t i0 = peak;
  while (i0 > 0 && rms[i0 - 1] >= dy0 + y1) i0--;

  size_t i1 = peak;
  for (size_t i = i0; i < peak; i++)
    if (rms[i] >= dy1 + y1) {
      i1 = i;
      break;
    }

  const ms<X>  x = env.x0();
  const double n = i1 - i0;
  return n * x.time_fraction();
}

template <typename X, typename Y, typename Z>
double release_time(const Spectrum<X, Y, Z>& env, double dy0, double dy1) {
  const Qtz<X>& qx   = env.qx;
  const Y*      rms  = env.ch.front();
  const size_t  peak = argmax(std::span<const Y>(rms, qx.n));
  const double  y1   = rms[peak];

  size_t i0 = peak;
  for (size_t i = peak; i < qx.n; i++)
    if (rms[i] <= dy1 + y1) {
      i0 = i;
      break;
    }

  size_t i1 = qx.n - 1;
  for (size_t i = i0; i < qx.n; i++)
    if (rms[i] <= dy0 + y1) {
      i1 = i;
      break;
    }
  const ms<X>  x = env.x0();
  const double n = i1 - i0;
  return n * x.time_fraction();
}

template <typename X, typename Y, typename Z>
double sustain_level(const Spectrum<X, Y, Z>& env, double dy1) {
  const Qtz<X>& qx   = env.qx;
  const Y*      rms  = env.ch.front();
  const size_t  peak = argmax(std::span<const Y>(rms, qx.n));
  const double  y1   = rms[peak];

  size_t i0 = peak;
  while (i0 > 0 && rms[i0 - 1] >= dy1 + y1) i0--;

  size_t i1 = peak;
  for (size_t i = peak; i < qx.n; i++)
    if (rms[i] <= dy1 + y1) {
      i1 = i;
      break;
    }

  double tot = 0.0;
  size_t n   = 0;
  for (size_t i = i0; i < i1; i++) {
    tot += rms[i];
    n++;
  }
  return tot / n - y1;
}

template <typename X, typename Y, typename Z>
double noise_floor(const Spectrum<X, Y, Z>& env, double dy0) {
  const Qtz<X>& qx   = env.qx;
  const Y*      rms  = env.ch.front();
  const size_t  peak = argmax(std::span<const Y>(rms, qx.n));
  const double  y1   = rms[peak];

  vec level;
  for (size_t i = 0; i < qx.n; i++) level.push_back(on_axis(rms[i], env.qy));
  return average_below(level, std::identity{}, dy0 + y1) - y1;
}

template <typename X, typename Y, typename Z>
std::optional<double> envelope_level(const Spectrum<X, Y, Z>& env) {
  const double y0  = env.qy.pos;
  const Y*     rms = env.ch.front();
  vec          lvl;
  for (size_t i = 0; i < env.qx.n; i++)
    if (rms[i] > y0) lvl.push_back(rms[i]);
  const vec weight(lvl.size(), 1.0);
  return otsus_method(lvl, weight, k::Lin);
}

template <typename X, typename Y, typename Z>
double decay_rate(const Spectrum<X, Y, Z>& env, double dy0) {
  const Qtz<X>& qx   = env.qx;
  const Y*      rms  = env.ch.front();
  const size_t  peak = argmax(std::span<const Y>(rms, qx.n));
  const double  y1   = rms[peak];

  vec lvl;
  for (size_t i = peak; i < qx.n && rms[i] > dy0 + y1; i++) lvl.push_back(rms[i]);
  const ms<X> x = env.x0();
  return slope(lvl, x.time_fraction());
}

template <typename X, typename Y, typename Z>
double temporal_centroid(const Spectrum<X, Y, Z>& env) {
  double weight = 0.0;
  double moment = 0.0;
  for (ms<X> x = env.x0(); x < env.x1(); ++x) {
    const double power = db<Y>::lin(env[x]);
    weight += power;
    moment += x * power;
  }
  return moment / weight;
}

template <typename X, typename Y, typename Z>
vec<X> macro_point(const Spectrum<X, Y, Z>& env, double smooth_ms, double tolerance) {
  vvec<Y>                 smooth_ch;
  const Spectrum<X, Y, Z> smooth = moving_average(env, smooth_ms, smooth_ch);
  vec<X>                  xs;
  vec                     lvl;
  for (ms<X> x = smooth.x0(); x < smooth.x1(); ++x) {
    xs.push_back(x);
    lvl.push_back(smooth[x]);
  }
  const double y1 = *std::ranges::max_element(lvl);
  const double x0 = xs.front();
  const double dx = xs.back() - x0;

  vec x_norm;
  vec y_norm;
  for (size_t i = 0; i < xs.size(); ++i) {
    x_norm.push_back((xs[i] - x0) / dx);
    y_norm.push_back(lvl[i] - y1);
  }

  const vec<size_t> vertex = rdp(x_norm, y_norm, tolerance);
  vec<X>            x_vertex;
  for (size_t i : vertex) x_vertex.push_back(xs[i]);
  return x_vertex;
}

inline double lag_one_autocorr(const vec<>& deviation) {
  const double mean = average(deviation, std::identity{});
  double       num  = 0.0;
  double       den  = 0.0;
  for (size_t i = 0; i + 1 < deviation.size(); i++)
    num += (deviation[i] - mean) * (deviation[i + 1] - mean);
  for (double v : deviation) {
    const double d = v - mean;
    den += d * d;
  }
  const double r1 = num / den;
  return r1;
}

template <typename X, typename Y, typename Z>
vec<> spike_interval(const Spectrum<X, Y, Z>& env, double smooth_ms, double deviation) {
  vvec<Y>                 smooth_ch;
  const Spectrum<X, Y, Z> smooth = moving_average(env, smooth_ms, smooth_ch);

  const double y0 = env.qy.pos;
  vec          ratio;
  vec<X>       ratio_time;
  for (ms<X> x = env.x0(); x < env.x1(); ++x) {
    ratio_time.push_back(x);
    if (smooth[x] > y0)
      ratio.push_back(env[x] - smooth[x]);
    else
      ratio.push_back(db<Y>::ratio(1.0));
  }
  const size_t n = ratio.size();

  const double mean      = average(ratio, std::identity{});
  const double sd        = std::sqrt(variance(ratio, std::identity{}, 0));
  const double threshold = mean + deviation * sd;

  vec xs;
  for (size_t i = 1; i + 1 < n; i++)
    if (ratio[i] > threshold && ratio[i] >= ratio[i - 1] && ratio[i] >= ratio[i + 1])
      xs.push_back(ratio_time[i]);

  return difference(xs, std::identity{});
}
}
