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
  std::vector<k::Point<>> point;
  Distribution residual;
  double wander;
};

template <typename X, typename Y, typename Z>
Drift drift(const Spectrogram<X, Y, Z>& env, double smooth_ms, double tolerance) {
  std::vector<std::vector<Y>> smooth_channel;
  const Spectrogram<X, Y, Z> smooth = moving_average(env, smooth_ms, smooth_channel);

  const double axis = env.y_db.offset;
  std::vector<double> residual;
  for (ms<X> x = env.begin(); x < env.end(); ++x)
    if (smooth[x] > axis && env[x] > axis)
      residual.push_back(env[x] - smooth[x]);
  const std::vector<double> residual_sorted = sort(residual);
  const std::span<const double> span{residual_sorted};
  std::vector<Distribution> ranked = rank(estimate_all(span), fidelity, span);
  const double wander = residual_autocorr(env, smooth_ms);
  const std::vector<X> x_vertex = macro_point(env, smooth_ms, tolerance);
  std::vector<k::Point<>> point;
  for (X x : x_vertex) point.push_back(k::Point<>{x, env[x]});

  return Drift{point, ranked.front(), wander};
}

template <typename X, typename Y, typename Z>
std::optional<Distribution> recurrence(const Spectrogram<X, Y, Z>& envelope, double smooth_ms, double z,
                                       double significance) {
  const std::vector<double>   gap = spike_interval(envelope, smooth_ms, z);
  std::optional<Distribution> result;
  if (gap.size() >= 2) {
    const std::vector<double> gap_sorted = sort(gap);
    std::vector<Distribution> candidate = estimate_all(gap_sorted);
    std::vector<Distribution> ranked
      = rank(std::move(candidate), fidelity, std::span<const double>(gap_sorted));
    if (fidelity(ranked.front(), std::span<const double>(gap_sorted)) >= significance)
      result = ranked.front();
  }
  return result;
}

}
