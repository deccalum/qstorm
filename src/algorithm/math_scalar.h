#pragma once
#include <algorithm>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <numeric>
#include <ranges>
#include <utility>
#include <vector>

#include "algorithm/math_constant.h"

namespace qstorm::algorithm {

template <std::ranges::input_range R, typename P>
constexpr double sum(R&& r, P p) {
  double total = 0.0;
  for (auto&& e : r) {
    const double x = static_cast<double>(std::invoke(p, e));
    total += x;
  }
  return total;
}

template <std::ranges::forward_range R>
constexpr size_t argmax(R&& r) {
  auto distance = std::ranges::distance(std::ranges::begin(r), std::ranges::max_element(r));
  return static_cast<size_t>(distance);
}

template <typename T, typename F, typename... R>
std::vector<T> rank(std::vector<T> candidate, F f, const R&... ref) {
  std::vector<double> fidelity;
  fidelity.reserve(candidate.size());
  for (const T& e : candidate) fidelity.push_back(f(e, ref...));
  std::vector<size_t> order(candidate.size());
  std::iota(order.begin(), order.end(), size_t{0});
  std::ranges::sort(order, std::ranges::greater{}, [&](size_t i) { return fidelity[i]; });
  std::vector<T> ranked;
  ranked.reserve(candidate.size());
  for (size_t i : order) ranked.push_back(std::move(candidate[i]));
  return ranked;
}

template <std::ranges::forward_range R, typename P>
std::vector<double> difference(R&& r, P p) {
  std::vector<double> delta;
  for (auto&& [prev, next] : r | std::views::adjacent<2>) {
    const double x = static_cast<double>(std::invoke(p, prev));
    const double y = static_cast<double>(std::invoke(p, next));
    delta.push_back(y - x);
  }
  return delta;
}

template <std::ranges::input_range R, typename P>
std::vector<double> running_total(R&& r, P p) {
  std::vector<double> cumulative{0.0};
  for (auto&& e : r) {
    const double x = static_cast<double>(std::invoke(p, e));
    cumulative.push_back(cumulative.back() + x);
  }
  return cumulative;
}

template <std::ranges::input_range R>
constexpr double average(R&& r) {
  double total = 0.0;
  size_t n = 0;
  for (const auto& value : r) {
    total += static_cast<double>(value);
    ++n;
  }
  return total / static_cast<double>(n);
}

template <std::ranges::input_range R, typename P>
  requires std::ranges::range<std::ranges::range_value_t<R>>
constexpr double average_nested(R&& r, P p) {
  return average(r | std::views::join, p);
}

template <std::ranges::input_range R, typename P>
constexpr double average(R&& r, P p) {
  double total = 0.0;
  size_t n = 0;
  for (auto&& e : r) {
    const double x = static_cast<double>(std::invoke(p, e));
    total += x;
    ++n;
  }
  return total / static_cast<double>(n);
}

template <typename R, typename P>
constexpr double average(R&& r, P p) {
  return average(r | std::views::join, p);
}

template <std::ranges::input_range R, std::indirect_unary_predicate<std::ranges::iterator_t<R>> P>
constexpr double ratio(R&& r, P p) {
  return average(r, p);
}

template <typename T, typename A, typename B>
constexpr T ratio(A numerator, B denominator) {
  const T top = static_cast<T>(numerator);
  const T bot = static_cast<T>(denominator);
  return top / bot;
}

template <std::ranges::input_range R>
std::vector<std::ranges::range_value_t<R>> sort(R&& r) {
  std::vector<std::ranges::range_value_t<R>> v = std::forward<R>(r) | std::ranges::to<std::vector>();
  std::ranges::sort(v);
  return v;
}

template <std::ranges::input_range R>
double percentile(R&& r, double fraction) {
  const std::vector v = sort(r);
  const size_t      n = v.size();
  const size_t   rank = static_cast<size_t>(fraction*static_cast<double>(n));
  const size_t  index = std::min(n - 1, rank);
  return static_cast<double>(v.at(index));
}

template <std::ranges::input_range R>
double median(R&& r) {
  const std::vector v = sort(r);
  const size_t mid = v.size() / 2;
  return (v.size() % 2 == 0)
         ? (static_cast<double>(v.at(mid - 1)) + static_cast<double>(v.at(mid))) / 2.0
         : static_cast<double>(v.at(mid));
}

template <typename R>
double median(R&& r) {
  return median(r | std::views::join);
}

template <std::ranges::forward_range R, typename P>
double variance(R&& r, P p, size_t d) {
  const double mean = average(r, p);
  double      total = 0.0;
  size_t          n = 0;
  for (auto&& e : r) {
    const double x = static_cast<double>(std::invoke(p, e));
    const double delta = x - mean;
    total += delta*delta;
    ++n;
  }
  return total / static_cast<double>(n - d);
}

template <std::ranges::forward_range X, std::ranges::forward_range Y>
double covariance(X&& x, Y&& y) {
  const double xbar = average(x, std::identity{});
  const double ybar = average(y, std::identity{});
  double      total = 0.0;
  size_t          n = 0;
  for (auto&& [xi, yi] : std::views::zip(x, y)) {
    const double dx = static_cast<double>(xi) -xbar;
    const double dy = static_cast<double>(yi) -ybar;
    total += dx*dy;
    ++n;
  }
  return total / static_cast<double>(n);
}

template <std::ranges::forward_range X, std::ranges::forward_range Y>
double slope(X&& x, Y&& y) {
  const double   cov = covariance(x, y);
  const double var_x = variance(x, std::identity{}, 0);
  return cov / var_x;
}

template <std::ranges::forward_range Y>
double slope(Y&& y, double dx) {
  const size_t n = static_cast<size_t>(std::ranges::distance(y));
  std::vector<double> x(n);
  for (size_t i = 0; i < n; ++i) x[i] = static_cast<double>(i) * dx;
  return slope(x, y);
}

template <std::ranges::forward_range X, std::ranges::forward_range Y>
double correlation(X&& x, Y&& y) {
  const double  cov = covariance(x, y);
  const double sd_x = std::sqrt(variance(x, std::identity{}, 0));
  const double sd_y = std::sqrt(variance(y, std::identity{}, 0));
  return cov / (sd_x*sd_y);
}

template <std::ranges::forward_range X, std::ranges::forward_range Y>
double cosine(X&& x, Y&& y) {
  double      dot = 0.0;
  double energy_x = 0.0;
  double energy_y = 0.0;
  for (auto&& [a, b] : std::views::zip(x, y)) {
    const double u = static_cast<double>(a);
    const double v = static_cast<double>(b);
    dot      += u*v;
    energy_x += u*u;
    energy_y += v*v;
  }
  return dot / std::sqrt(energy_x*energy_y);
}

template <std::ranges::input_range R, typename P>
double average_below(R&& r, P p, double threshold) {
  double total = 0.0;
  size_t n = 0;
  for (auto&& e : r) {
    const double x = static_cast<double>(std::invoke(p, e));
    if (x >= threshold) continue;
    total += x;
    ++n;
  }
  return total / static_cast<double>(n);
}

template <std::invocable<int> F>
double sigma(F f, int lo, int hi) {
  return std::ranges::fold_left(std::views::iota(lo, hi + 1) | std::views::transform(std::move(f)),
                                0.0, std::plus<>{});
}

inline double gamma_series(double a, double x) {
  double  term = 1.0 / a;
  double total = term;
  for (double ap = a + 1.0; std::abs(term) > std::abs(total)*k::Epsilon; ap += 1.0) {
    term *= x / ap;
    total += term;
  }
  return total;
}

inline double gamma_fraction(double a, double x) {
  double b = x + 1.0 - a;
  double c = 1.0 / k::Tiny;
  double d = 1.0 / b;
  double h = d;
  for (int i = 1;; ++i) {
    const double an = -i * (i - a);
    b += 2.0;
    d = an*d + b;
    if (std::abs(d) < k::Tiny) d = k::Tiny;
    c = b + an / c;
    if (std::abs(c) < k::Tiny) c = k::Tiny;
    d = 1.0 / d;
    const double del = d*c;
    h *= del;
    if (std::abs(del - 1.0) <= k::Epsilon) return h;
  }
}

inline double gamma_p(double a, double x) {
  const double prefactor = std::exp(a*std::log(x) - x - std::lgamma(a));
  return x < a + 1.0 ? prefactor*gamma_series(a, x) : 1.0 - prefactor*gamma_fraction(a, x);
}

inline double beta_fraction(double a, double b, double x) {
  double c = 1.0;
  double d = 1.0 - (a + b) * x / (a + 1.0);
  if (std::abs(d) < k::Tiny) d = k::Tiny;
  d = 1.0 / d;
  double h = d;
  for (int m = 1;; ++m) {
    const double even = m * (b - m) * x / ((a + 2.0*m - 1.0) * (a + 2.0*m));
    d = 1.0 + even*d;
    if (std::abs(d) < k::Tiny) d = k::Tiny;
    d = 1.0 / d;
    c = 1.0 + even / c;
    if (std::abs(c) < k::Tiny) c = k::Tiny;
    h *= d*c;
    const double odd = -(a + m) * (a + b + m) * x / ((a + 2.0*m) * (a + 2.0*m + 1.0));
    d = 1.0 + odd*d;
    if (std::abs(d) < k::Tiny) d = k::Tiny;
    d = 1.0 / d;
    c = 1.0 + odd / c;
    if (std::abs(c) < k::Tiny) c = k::Tiny;
    const double del = d*c;
    h *= del;
    if (std::abs(del - 1.0) <= k::Epsilon) return h;
  }
}

inline double beta_p(double a, double b, double x) {
  const double prefactor =
    std::exp(std::lgamma(a +b) -std::lgamma(a) -std::lgamma(b) +a*std::log(x) +b*std::log(1.0 -x));
  const bool flip = x > (a + 1.0) / (a + b + 2.0);
  return flip ? 1.0 - prefactor*beta_fraction(b, a, 1.0 - x) / b
              : prefactor*beta_fraction(a, b, x) / a;
}

template <std::invocable<double> F>
double bisect(F f, double lo, double hi, double target, int iterations) {
  const bool lo_above = f(lo) > target;
  for (int i = 0; i < iterations; ++i) {
    const double mid = std::midpoint(lo, hi);
    ((f(mid) > target) == lo_above ? lo : hi) = mid;
  }
  return std::midpoint(lo, hi);
}

template <std::ranges::forward_range R, typename P>
double kurtosis(R&& r, P p) {
  const double mean = average(r, p);
  double        mu2 = 0.0;
  double        mu4 = 0.0;
  size_t          n = 0;
  for (auto&& e : r) {
    const double d = std::invoke(p, e) - mean;
    const double d2 = d*d;
    mu2 += d2;
    mu4 += d2*d2;
    ++n;
  }
  mu2 /= static_cast<double>(n);
  mu4 /= static_cast<double>(n);
  return mu4 / (mu2*mu2) - 3.0;
}

}
