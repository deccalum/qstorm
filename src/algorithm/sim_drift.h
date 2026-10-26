#pragma once

#include <optional>
#include <span>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "algorithm/dsp_envelope.h"
#include "algorithm/dsp_window.h"
#include "algorithm/math_constant.h"
#include "algorithm/math_scalar.h"
#include "algorithm/sim_distribution.h"

namespace qstorm::algorithm {
struct Drift {
  vec<k::Point<>> point;
  Distribution    deviation;
  double          wander;
};

template <typename X, typename Y, typename Z>
Drift drift(const Spectrum<X, Y, Z>& env, double smooth_ms, double tolerance) {
  vvec<Y>                 smooth_channel;
  const Spectrum<X, Y, Z> smooth = moving_average(env, smooth_ms, smooth_channel);
  const double            y0     = env.qy.pos;
  vec                     deviation;
  for (ms<X> x = env.x0(); x < env.x1(); ++x)
    if (smooth[x] > y0 && env[x] > y0) deviation.push_back(env[x] - smooth[x]);
  const vec                     deviation_sorted = sort(deviation);
  const std::span<const double> span{deviation_sorted};
  vec<Distribution>             ranked   = rank(estimate_all(span), fidelity, span);
  const double                  wander   = lag_one_autocorr(deviation);
  const vec<X>                  x_vertex = macro_point(env, smooth_ms, tolerance);
  vec<k::Point<>>               point;
  for (X x : x_vertex) point.push_back(k::Point<>{x, env[x]});

  return Drift{point, ranked.front(), wander};
}

template <typename X, typename Y, typename Z>
std::optional<Distribution>
recurrence(const Spectrum<X, Y, Z>& env, double smooth_ms, double deviation, double significance) {
  const vec                   gap = spike_interval(env, smooth_ms, deviation);
  std::optional<Distribution> distribution;
  if (gap.size() >= 2) {
    const vec         gap_sorted = sort(gap);
    vec<Distribution> dist  = estimate_all(gap_sorted);
    vec<Distribution> ranked
      = rank(std::move(dist), fidelity, std::span<const double>(gap_sorted));
    if (fidelity(ranked.front(), std::span<const double>(gap_sorted)) >= significance)
      distribution = ranked.front();
  }
  return distribution;
}
}
