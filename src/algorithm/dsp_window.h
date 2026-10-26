#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <ranges>
#include <span>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "algorithm/math_scalar.h"

namespace qstorm::algorithm {
template <std::ranges::forward_range R, typename F>
vec<> window(R&& r, size_t w, F fn) {
  vec reduced;
  for (auto&& slice : r | std::views::slide(w)) reduced.push_back(fn(slice));
  return reduced;
}

template <std::ranges::random_access_range R>
vec<> trailing_percentile(R&& r, size_t w, double fraction) {
  const size_t n = std::ranges::size(r);
  w              = std::min(w, n);
  vec serie(n);
  for (size_t i = w - 1; i < n; ++i) {
    auto tail = r | std::views::drop(i + 1 - w) | std::views::take(w);
    serie[i]  = percentile(tail, fraction);
  }
  for (size_t i = 0; i + 1 < w; ++i) serie[i] = serie[w - 1];
  return serie;
}

template <std::ranges::forward_range L, std::ranges::forward_range R>
double side_ratio(L&& left, R&& right) {
  double nrg_mid  = 0.0;
  double nrg_side = 0.0;
  for (auto&& [a, b] : std::views::zip(left, right)) {
    const double l    = a;
    const double r    = b;
    const double mid  = (l + r) / 2;
    const double side = (l - r) / 2;
    nrg_mid += mid * mid;
    nrg_side += side * side;
  }
  return nrg_side / (nrg_mid + nrg_side);
}

template <std::ranges::forward_range L, std::ranges::forward_range R>
double energy_balance(L&& left, R&& right) {
  double nrg_l = 0.0;
  double nrg_r = 0.0;
  for (auto&& [a, b] : std::views::zip(left, right)) {
    const double left  = a;
    const double right = b;
    nrg_l += std::abs(left);
    nrg_r += std::abs(right);
  }
  return nrg_r / (nrg_l + nrg_r);
}

template <std::ranges::random_access_range M>
vec<size_t> sustain_edge(M&& mask, size_t hold, size_t refractory) {
  const auto   m = std::ranges::begin(mask);
  const size_t n = std::ranges::size(mask);
  vec<size_t>  hit;
  const long   space = refractory;
  long         i0    = -space;

  size_t i = 0;
  while (i < n) {
    const long at    = i;
    const long since = at - i0;
    const bool ready = since >= space;
    if (m[i] && ready) {
      const size_t hold_end = std::min(n, i + hold);
      bool         held     = true;
      for (size_t j = i; j < hold_end; ++j) {
        if (!m[j]) {
          held = false;
          break;
        }
      }
      if (held) {
        hit.push_back(i);
        i0 = i;
        i += refractory;
        continue;
      }
    }
    ++i;
  }
  return hit;
}

template <typename X, typename Y, typename Z>
Spectrum<X, Y, Z> moving_average(const Spectrum<X, Y, Z>& env, X dx, vvec<Y>& ch) {
  const size_t nx     = env.qx.n;
  const size_t radius = env.qx.frame(dx);
  const vec    tot    = run_sum(std::span<const Y>(env.ch.front(), nx), std::identity{});
  ch.assign(1, vec<Y>(nx));
  vec<Y>& smooth = ch.front();
  for (size_t i = 0; i < nx; ++i) {
    const size_t i0 = i >= radius ? i - radius : 0;
    const size_t i1 = std::min(nx, i + radius + 1);
    const double n  = i1 - i0;
    smooth[i]       = (tot[i1] - tot[i0]) / n;
  }
  return as_spectrum(ch, env.qx, env.qy, env.qz);
}

inline vec<> periodic_hann(size_t n) {
  vec hann;
  for (size_t i = 0; i < n; ++i) {
    const double theta = k::tau * i / n;
    hann.push_back((1.0 - std::cos(theta) / 2));
  }
  return hann;
}

template <std::ranges::random_access_range X, std::ranges::random_access_range Y, typename P>
long peak_lag(X&& x, Y&& y, long max_lag, P p) {
  const auto xs = std::ranges::begin(x);
  const auto ys = std::ranges::begin(y);
  const long n  = std::ranges::size(x);

  double nrg_x = 0.0;
  for (auto&& e : x) {
    const double v = std::invoke(p, e);
    nrg_x += v * v;
  }

  double nrg_y = 0.0;
  for (auto&& e : y) {
    const double v = std::invoke(p, e);
    nrg_y += v * v;
  }
  const double norm = std::sqrt(nrg_x * nrg_y);

  long   arg = 0;
  double y1  = -2.0;
  for (long lag = -max_lag; lag <= max_lag; ++lag) {
    const long i0    = std::max<long>(0, -lag);
    const long i1    = std::min(n, n - lag);
    double     tot   = 0.0;
    for (long i = i0; i < i1; ++i) {
      const double a = std::invoke(p, xs[i]);
      const double b = std::invoke(p, ys[i + lag]);
      tot += a * b;
    }
    const double value = tot / norm;
    if (value > y1) {
      y1  = value;
      arg = lag;
    }
  }
  return arg;
}
}
