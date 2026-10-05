#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <ranges>
#include <span>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "algorithm/math_constant.h"
#include "algorithm/math_scalar.h"

namespace qstorm::algorithm {

template <typename T>
size_t centroid_index(std::span<const T> t) {
  const double  mean = average(t, std::identity{});
  double       total = 0.0;
  size_t           n = 0;
  for (size_t i = 0; i < t.size(); ++i) {
    if (t[i] < mean) continue;
    total += static_cast<double>(i);
    ++n;
  }
  return static_cast<size_t>(total / static_cast<double>(n));
}

template <std::ranges::forward_range R, typename F>
std::vector<double> window(R&& r, size_t w, F f) {
  std::vector<double> reduced;
  for (auto&& slice : r | std::views::slide(w)) reduced.push_back(f(slice));
  return reduced;
}

template <std::ranges::random_access_range R>
std::vector<double> trailing_percentile(R&& r, size_t w, double fraction) {
  const size_t n = static_cast<size_t>(std::ranges::size(r));
  w = std::min(w, n);
  std::vector<double> series(n);
  for (size_t i = w - 1; i < n; ++i) {
    auto tail = r | std::views::drop(i + 1 - w) | std::views::take(w);
    series[i] = percentile(tail, fraction);
  }
  for (size_t i = 0; i + 1 < w; ++i) series[i] = series[w - 1];
  return series;
}

template <std::ranges::forward_range L, std::ranges::forward_range R>
double side_ratio(L&& l, R&& r) {
  double  energy_mid = 0.0;
  double energy_side = 0.0;
  for (auto&& [a, b] : std::views::zip(l, r)) {
    const double  mid = 0.5 * (static_cast<double>(a) +static_cast<double>(b));
    const double side = 0.5 * (static_cast<double>(a) -static_cast<double>(b));
    energy_mid  += mid*mid;
    energy_side += side*side;
  }
  return energy_side / (energy_mid + energy_side);
}

template <std::ranges::forward_range L, std::ranges::forward_range R>
double energy_balance(L&& l, R&& r) {
  double energy_l = 0.0;
  double energy_r = 0.0;
  for (auto&& [a, b] : std::views::zip(l, r)) {
    energy_l += std::abs(static_cast<double>(a));
    energy_r += std::abs(static_cast<double>(b));
  }
  return energy_r / (energy_l + energy_r);
}

template <std::ranges::input_range R, typename P, typename M>
std::optional<double> otsus_method(R&& r, P p, M m, double min_significance) {
  std::vector<double> v;
  for (auto&& e : r) {
    const double x = m.fwd(std::invoke(p, e));
    if (!std::isfinite(x)) continue;
    v.push_back(x);
  }
  std::ranges::sort(v);
  const size_t          n = v.size();
  std::optional<double> result;
  if (n >= 4) {
    const std::vector<double> cumulative = qstorm::algorithm::running_total(v, std::identity{});
    const double total = cumulative.back();

    size_t best_i = 1;
    double best_between = -1.0;
    for (size_t i = 1; i < n; ++i) {
      const double      w0 = static_cast<double>(i) / static_cast<double>(n);
      const double      m0 = cumulative[i] / static_cast<double>(i);
      const double      m1 = (total - cumulative[i]) / static_cast<double>(n - i);
      const double   delta = m0 - m1;
      const double  weight = w0 * (1.0 - w0);
      const double between = weight*delta*delta;
      if (between > best_between) {
        best_between = between;
        best_i = i;
      }
    };

    const double v0 = qstorm::algorithm::variance(v | std::views::take(best_i), std::identity{}, 0);
    const double v1 = qstorm::algorithm::variance(v | std::views::drop(best_i), std::identity{}, 0);
    const double w0 = v0*static_cast<double>(best_i);
    const double w1 = v1*static_cast<double>(n - best_i);
    const double       within = (w0 + w1) / static_cast<double>(n);
    const double significance = within > k::Epsilon
                              ? best_between / within
                              : std::numeric_limits<double>::infinity();
    if (significance >= min_significance) {
      const double mid = std::midpoint(v.at(best_i - 1), v.at(best_i));
      result = m.inv(mid);
    }
  }
  return result;
}

template <std::ranges::input_range R, typename P, typename Typical, typename M> std::optional<double>
widest_gap_threshold(R&& r, P p, Typical typical_of, M m, double factor, double search_fraction) {
  std::vector<double> v;
  for (auto&& e : r) {
    const double x = m.fwd(std::invoke(p, e));
    if (!std::isfinite(x)) continue;
    v.push_back(x);
  }
  std::ranges::sort(v);
  const size_t          n = v.size();
  std::optional<double> result;
  if (n >= 4) {
    const std::vector<double> gap = qstorm::algorithm::difference(v, std::identity{});
    const size_t w
      = std::max<size_t>(1, static_cast<size_t>(search_fraction*static_cast<double>(gap.size())));
    const size_t peak = qstorm::algorithm::argmax(gap | std::views::take(w));

    const double typical = typical_of(gap);
    if (gap[peak] >= factor*typical) result = m.inv(0.5 * (v.at(peak) + v.at(peak + 1)));
  }
  return result;
}

template <std::ranges::random_access_range Mask>
std::vector<size_t> sustain_edge(Mask&& mask, size_t hold, size_t refractory) {
  const auto   m = std::ranges::begin(mask);
  const size_t n = static_cast<size_t>(std::ranges::size(mask));
  std::vector<size_t> hit;
  long last = -static_cast<long>(refractory);

  size_t i = 0;
  while (i < n) {
    const long since = static_cast<long>(i) - last;
    const bool ready = since >= static_cast<long>(refractory);
    if (m[i] && ready) {
      const size_t hold_end = std::min(n, i + hold);
      bool         held = true;
      for (size_t j = i; j < hold_end; ++j) {
        if (!m[j]) {
          held = false;
          break;
        }
      }
      if (held) {
        hit.push_back(i);
        last = static_cast<long>(i);
        i += refractory;
        continue;
      }
    }
    ++i;
  }
  return hit;
}

constexpr double iou(double a0, double a1, double b0, double b1) {
  const double inter = std::min(a1, b1) - std::max(a0, b0);
  const double   uni = std::max(a1, b1) - std::min(a0, b0);
  return (inter <= 0.0 || uni <= 0.0) ? 0.0 : inter / uni;
}

template <typename X, typename Y, typename Z> Spectrogram<X, Y, Z>
moving_average(const Spectrogram<X, Y, Z>& env, X half, std::vector<std::vector<Y>>& channel) {
  channel.assign(1, std::vector<Y>());
  std::vector<Y>&           smooth = channel.front();
  const Spectrogram<X, Y, Z> reach = env.window(env.begin(), env.begin() + half);
  const X                     snap = reach.end() - reach.begin();
  for (ms<X> x = env.begin(); x < env.end(); ++x) {
    ms<X> next{x + snap, &env.x_ms};
    ++next;
    const X lo = std::max<X>(x - snap, env.begin());
    const X hi = std::min<X>(next, env.end());
    const Spectrogram<X, Y, Z> slice = env.window(lo, hi);
    double total = 0;
    size_t     n = 0;
    for (ms<X> slice_x = slice.begin(); slice_x < slice.end(); ++slice_x) {
      total += slice[slice_x];
      n++;
    }
    smooth.push_back(static_cast<Y>(total / static_cast<double>(n)));
  }
  return as_spectrogram(channel, env.x_ms, env.z_hz, env.y_db);
}

template <std::ranges::random_access_range X, std::ranges::random_access_range Y, typename P>
long peak_lag(X&& x, Y&& y, long max_lag, P p) {
  const auto xs = std::ranges::begin(x);
  const auto ys = std::ranges::begin(y);
  const long  n = static_cast<long>(std::ranges::size(x));

  double energy_x = 0.0;
  for (auto&& e : x) {
    const double v = static_cast<double>(std::invoke(p, e));
    energy_x += v*v;
  }

  double energy_y = 0.0;
  for (auto&& e : y) {
    const double v = static_cast<double>(std::invoke(p, e));
    energy_y += v*v;
  }
  const double norm = std::sqrt(energy_x*energy_y);

  long  peak = 0;
  double top = -2.0;
  for (long lag = -max_lag; lag <= max_lag; ++lag) {
    const long lo = std::max<long>(0, -lag);
    const long hi = std::min(n, n - lag);
    double  total = 0.0;
    for (long i = lo; i < hi; ++i) {
      const double a = static_cast<double>(std::invoke(p, xs[i]));
      const double b = static_cast<double>(std::invoke(p, ys[i + lag]));
      total += a*b;
    }
    const double value = total / norm;
    if (value > top) {
      top  = value;
      peak = lag;
    }
  }
  return peak;
}

template <std::ranges::random_access_range X, std::ranges::random_access_range Y>
size_t chord_peak(X&& x, Y&& y, size_t lo, size_t hi) {
  const double  dx = x[hi] - x[lo];
  const double  dy = y[hi] - y[lo];
  const double len = std::hypot(dx, dy);
  size_t      peak = lo + 1;
  double       top = 0.0;
  for (size_t i = lo + 1; i < hi; i++) {
    const double distance = std::abs(dx * (y[i] - y[lo]) - dy * (x[i] - x[lo])) / len;
    if (distance > top) {
      top = distance;
      peak = i;
    }
  }
  return peak;
}

template <std::ranges::random_access_range X, std::ranges::random_access_range Y>
void rdp_split(X&& x, Y&& y, double tolerance, std::vector<size_t>& vertex, size_t lo, size_t hi) {
  if (hi > lo + 1) {
    const size_t peak = chord_peak(x, y, lo, hi);
    const double  dx = x[hi] - x[lo];
    const double  dy = y[hi] - y[lo];
    const double len = std::hypot(dx, dy);
    const double distance = std::abs(dx * (y[peak] - y[lo]) - dy * (x[peak] - x[lo])) / len;
    if (distance > tolerance) {
      vertex.push_back(peak);
      rdp_split(x, y, tolerance, vertex, lo, peak);
      rdp_split(x, y, tolerance, vertex, peak, hi);
    }
  }
}

template <std::ranges::random_access_range X, std::ranges::random_access_range Y>
std::vector<size_t> rdp(X&& x, Y&& y, double tolerance) {
  std::vector<size_t> vertex{0, std::ranges::size(x) - 1};
  rdp_split(x, y, tolerance, vertex, size_t{0}, vertex.back());
  std::ranges::sort(vertex);
  return vertex;
}

}
