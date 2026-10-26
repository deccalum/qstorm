#pragma once
#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <numeric>
#include <optional>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

#include "algorithm/math_constant.h"

namespace qstorm::algorithm {
template <std::ranges::input_range R, typename P>
constexpr double sum(R&& r, P p) {
  double tot = 0.0;
  for (auto&& e : r) {
    const double x = std::invoke(p, e);
    tot += x;
  }
  return tot;
}

template <std::ranges::forward_range R>
constexpr size_t argmax(R&& r) {
  auto distance = std::ranges::distance(std::ranges::begin(r), std::ranges::max_element(r));
  return distance;
}

template <typename T, typename F, typename... R>
vec<T> rank(vec<T> input, F fn, const R&... ref) {
  vec<std::pair<double, size_t>> order;
  order.reserve(input.size());
  for (size_t i = 0; i < input.size(); ++i) order.emplace_back(fn(input[i], ref...), i);
  std::ranges::sort(order, std::ranges::greater{});
  vec<T> ranked;
  ranked.reserve(input.size());
  for (const auto& [fidelity, i] : order) ranked.push_back(std::move(input[i]));
  return ranked;
}

template <std::ranges::forward_range R, typename P>
vec<> difference(R&& r, P p) {
  vec delta;
  for (auto&& [prev, next] : r | std::views::adjacent<2>) {
    const double x = std::invoke(p, prev);
    const double y = std::invoke(p, next);
    delta.push_back(y - x);
  }
  return delta;
}

template <std::ranges::input_range R, typename P>
vec<> run_sum(R&& r, P p) {
  vec cumulative{0.0};
  for (auto&& e : r) {
    const double x = std::invoke(p, e);
    cumulative.push_back(cumulative.back() + x);
  }
  return cumulative;
}

template <std::ranges::input_range R>
constexpr double average(R&& r) {
  double tot = 0.0;
  size_t n   = 0;
  for (const auto& value : r) {
    tot += value;
    ++n;
  }
  return tot / n;
}

template <std::ranges::input_range R, typename P>
  requires std::ranges::range<std::ranges::range_value_t<R>>
constexpr double average_nested(R&& r, P p) {
  return average(r | std::views::join, p);
}

template <std::ranges::input_range R, typename P>
constexpr double average(R&& r, P p) {
  double tot = 0.0;
  size_t n   = 0;
  for (auto&& e : r) {
    const double x = std::invoke(p, e);
    tot += x;
    ++n;
  }
  return tot / n;
}

template <typename T, typename A, typename B>
constexpr T ratio(A numerator, B denominator) {
  const T num = numerator;
  const T den = denominator;
  return num / den;
}

template <std::ranges::input_range R>
vec<std::ranges::range_value_t<R>> sort(R&& r) {
  vec<std::ranges::range_value_t<R>> v = std::forward<R>(r) | std::ranges::to<std::vector>();
  std::ranges::sort(v);
  return v;
}

template <std::ranges::input_range R>
double percentile(R&& r, double fraction) {
  const std::vector v     = sort(r);
  const size_t      n     = v.size();
  const size_t      rank  = fraction * n;
  const size_t      index = std::min(n - 1, rank);
  return v.at(index);
}

template <std::ranges::input_range R>
double median(R&& r) {
  const std::vector v     = sort(r);
  const size_t      mid   = v.size() / 2;
  const double      m0    = v.at((v.size() - 1) / 2);
  const double      m1    = v.at(mid);
  return (m0 + m1) / 2.0;
}

template <typename R>
double median(R&& r) {
  return median(r | std::views::join);
}

template <std::ranges::forward_range R, typename P>
double variance(R&& r, P p, size_t d) {
  const double mean = average(r, p);
  double       tot  = 0.0;
  size_t       n    = 0;
  for (auto&& e : r) {
    const double x     = std::invoke(p, e);
    const double delta = x - mean;
    tot += delta * delta;
    ++n;
  }
  return tot / (n - d);
}

template <std::ranges::forward_range X, std::ranges::forward_range Y>
double covariance(X&& x, Y&& y) {
  const double xbar  = average(x, std::identity{});
  const double ybar  = average(y, std::identity{});
  double       tot   = 0.0;
  size_t       n     = 0;
  for (auto&& [xi, yi] : std::views::zip(x, y)) {
    const double dx = xi - xbar;
    const double dy = yi - ybar;
    tot += dx * dy;
    ++n;
  }
  return tot / n;
}

template <std::ranges::forward_range X, std::ranges::forward_range Y>
double slope(X&& x, Y&& y) {
  const double cov   = covariance(x, y);
  const double var_x = variance(x, std::identity{}, 0);
  return cov / var_x;
}

template <std::ranges::forward_range Y>
double slope(Y&& y, double dx) {
  const size_t n = std::ranges::distance(y);
  vec          x(n);
  for (size_t i = 0; i < n; ++i) x[i] = i * dx;
  return slope(x, y);
}

template <std::ranges::forward_range X, std::ranges::forward_range Y>
double correlation(X&& x, Y&& y) {
  const double cov  = covariance(x, y);
  const double sd_x = std::sqrt(variance(x, std::identity{}, 0));
  const double sd_y = std::sqrt(variance(y, std::identity{}, 0));
  return cov / (sd_x * sd_y);
}

template <std::ranges::forward_range X, std::ranges::forward_range Y>
double cosine(X&& x, Y&& y) {
  double dot      = 0.0;
  double nrg_x = 0.0;
  double nrg_y = 0.0;
  for (auto&& [a, b] : std::views::zip(x, y)) {
    const double u = a;
    const double v = b;
    dot += u * v;
    nrg_x += u * u;
    nrg_y += v * v;
  }
  return dot / std::sqrt(nrg_x * nrg_y);
}

template <std::ranges::input_range R, typename P>
double average_below(R&& r, P p, double threshold) {
  double tot = 0.0;
  size_t n   = 0;
  for (auto&& e : r) {
    const double x = std::invoke(p, e);
    if (x >= threshold) continue;
    tot += x;
    ++n;
  }
  return tot / n;
}

template <std::invocable<int> F>
double sigma(F fn, int i0, int i1) {
  return std::ranges::fold_left(std::views::iota(i0, i1 + 1) | std::views::transform(std::move(fn)),
                                0.0, std::plus<>{});
}

inline double gamma_series(double a, double x) {
  double term = 1.0 / a;
  double tot  = term;
  for (double ap = a + 1.0; std::abs(term) > std::abs(tot) * k::Epsilon; ap += 1.0) {
    term *= x / ap;
    tot += term;
  }
  return tot;
}

inline double gamma_fraction(double a, double x) {
  double b = x + 1.0 - a;
  double c = 1.0 / k::Tiny;
  double d = 1.0 / b;
  double h = d;
  for (int i = 1;; ++i) {
    const double an = -i * (i - a);
    b += 2.0;
    d = an * d + b;
    if (std::abs(d) < k::Tiny) d = k::Tiny;
    c = b + an / c;
    if (std::abs(c) < k::Tiny) c = k::Tiny;
    d                = 1.0 / d;
    const double del = d * c;
    h *= del;
    if (std::abs(del - 1.0) <= k::Epsilon) return h;
  }
}

inline double gamma_p(double a, double x) {
  const double prefactor = std::exp(a * std::log(x) - x - std::lgamma(a));
  return x < a + 1.0 ? prefactor * gamma_series(a, x) : 1.0 - prefactor * gamma_fraction(a, x);
}

inline double beta_fraction(double a, double b, double x) {
  double c = 1.0;
  double d = 1.0 - (a + b) * x / (a + 1.0);
  if (std::abs(d) < k::Tiny) d = k::Tiny;
  d        = 1.0 / d;
  double h = d;
  for (int m = 1;; ++m) {
    const double even = m * (b - m) * x / ((a + 2.0 * m - 1.0) * (a + 2.0 * m));
    d                 = 1.0 + even * d;
    if (std::abs(d) < k::Tiny) d = k::Tiny;
    d = 1.0 / d;
    c = 1.0 + even / c;
    if (std::abs(c) < k::Tiny) c = k::Tiny;
    h *= d * c;
    const double odd = -(a + m) * (a + b + m) * x / ((a + 2.0 * m) * (a + 2.0 * m + 1.0));
    d                = 1.0 + odd * d;
    if (std::abs(d) < k::Tiny) d = k::Tiny;
    d = 1.0 / d;
    c = 1.0 + odd / c;
    if (std::abs(c) < k::Tiny) c = k::Tiny;
    const double del = d * c;
    h *= del;
    if (std::abs(del - 1.0) <= k::Epsilon) return h;
  }
}

inline double beta_p(double a, double b, double x) {
  const double prefactor = std::exp(std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b)
                                    + a * std::log(x) + b * std::log(1.0 - x));

  const bool flip = x > (a + 1.0) / (a + b + 2.0);
  return flip ? 1.0 - prefactor * beta_fraction(b, a, 1.0 - x) / b
              : prefactor * beta_fraction(a, b, x) / a;
}

template <std::invocable<double> F>
double bisect(F fn, double x0, double x1, double target, int iterations) {
  const bool x0_above = fn(x0) > target;
  for (int i = 0; i < iterations; ++i) {
    const double mid = std::midpoint(x0, x1);

    ((fn(mid) > target) == x0_above ? x0 : x1) = mid;
  }
  return std::midpoint(x0, x1);
}

template <std::ranges::forward_range R, typename P>
double kurtosis(R&& r, P p) {
  const double mean = average(r, p);
  double       mu2  = 0.0;
  double       mu4  = 0.0;
  size_t       n    = 0;
  for (auto&& e : r) {
    const double d  = std::invoke(p, e) - mean;
    const double d2 = d * d;
    mu2 += d2;
    mu4 += d2 * d2;
    ++n;
  }
  mu2 /= n;
  mu4 /= n;
  return mu4 / (mu2 * mu2) - 3.0;
}

template <typename T>
size_t centroid_index(std::span<const T> t) {
  const double mean = average(t, std::identity{});
  double       tot  = 0.0;
  size_t       n    = 0;
  for (size_t i = 0; i < t.size(); ++i) {
    if (t[i] < mean) continue;
    tot += i;
    ++n;
  }
  return tot / n;
}

template <typename M>
std::optional<double>
otsus_method(std::span<const double> value, std::span<const double> weight, M m) {
  vec<std::pair<double, double>> v;
  for (size_t i = 0; i < value.size(); ++i) {
    const double x = m.fwd(value[i]);
    if (!std::isfinite(x)) continue;
    v.emplace_back(x, weight[i]);
  }
  std::ranges::sort(v);
  const size_t          n = v.size();
  std::optional<double> result;
  if (n >= 4) {
    vec omega{0.0};
    vec mu{0.0};
    for (const auto& [x, w] : v) {
      omega.push_back(omega.back() + w);
      mu.push_back(mu.back() + w * x);
    }
    const double omega_t = omega.back();
    const double mu_t    = mu.back();

    size_t arg = 1;
    double y1  = -1.0;
    for (size_t i = 1; i < n; ++i) {
      const double w0     = omega[i] / omega_t;
      const double m0     = mu[i] / omega[i];
      const double m1     = (mu_t - mu[i]) / (omega_t - omega[i]);
      const double delta  = m0 - m1;
      const double weight = w0 * (1.0 - w0);
      const double var_b  = weight * delta * delta;
      if (var_b > y1) {
        y1  = var_b;
        arg = i;
      }
    }
    result = m.inv(std::midpoint(v.at(arg - 1).first, v.at(arg).first));
  }
  return result;
}

constexpr double iou(double a0, double a1, double b0, double b1) {
  const double inter = std::min(a1, b1) - std::max(a0, b0);
  const double uni   = std::max(a1, b1) - std::min(a0, b0);
  return (inter <= 0.0 || uni <= 0.0) ? 0.0 : inter / uni;
}

template <std::ranges::random_access_range X, std::ranges::random_access_range Y>
double chord_distance(X&& x, Y&& y, size_t i0, size_t i1, size_t i) {
  const double dx  = x[i1] - x[i0];
  const double dy  = y[i1] - y[i0];
  const double len = std::hypot(dx, dy);
  return std::abs(dx * (y[i] - y[i0]) - dy * (x[i] - x[i0])) / len;
}

template <std::ranges::random_access_range X, std::ranges::random_access_range Y>
size_t chord_peak(X&& x, Y&& y, size_t i0, size_t i1) {
  size_t arg = i0 + 1;
  double y1  = 0.0;
  for (size_t i = i0 + 1; i < i1; i++) {
    const double distance = chord_distance(x, y, i0, i1, i);
    if (distance > y1) {
      y1  = distance;
      arg = i;
    }
  }
  return arg;
}

template <std::ranges::random_access_range X, std::ranges::random_access_range Y>
void rdp_split(X&& x, Y&& y, double tol, vec<size_t>& vertex, size_t i0, size_t i1) {
  if (i1 > i0 + 1) {
    const size_t peak     = chord_peak(x, y, i0, i1);
    const double distance = chord_distance(x, y, i0, i1, peak);
    if (distance > tol) {
      vertex.push_back(peak);
      rdp_split(x, y, tol, vertex, i0, peak);
      rdp_split(x, y, tol, vertex, peak, i1);
    }
  }
}

template <std::ranges::random_access_range X, std::ranges::random_access_range Y>
vec<size_t> rdp(X&& x, Y&& y, double tol) {
  vec<size_t> vertex{0, std::ranges::size(x) - 1};
  rdp_split(x, y, tol, vertex, size_t{0}, vertex.back());
  std::ranges::sort(vertex);
  return vertex;
}
}
