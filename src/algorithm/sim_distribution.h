#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iterator>
#include <numbers>
#include <optional>
#include <span>
#include <type_traits>
#include <variant>
#include <vector>

#include "algorithm/dsp_window.h"
#include "algorithm/math_constant.h"
#include "algorithm/math_scalar.h"

namespace qstorm::algorithm {

struct Normal          {double mean;     double sigma;};
struct Lognormal       {double mu;       double sigma;};
struct Gamma           {double shape;    double scale;};
struct Weibull         {double shape;    double scale;};
struct Exponential     {double rate;                  };
struct ChiSquared      {double dof;                   };
struct ExtremeValue    {double location; double scale;};
struct Cauchy          {double location; double scale;};
struct StudentT        {double dof;      double location; double scale;};
struct PiecewiseLinear {std::vector<k::Point<double, double>> point;};
using Distribution = std::variant<Normal, Lognormal, Gamma, Weibull, Exponential, ChiSquared, ExtremeValue, Cauchy, StudentT, PiecewiseLinear>;

constexpr double kWeibullShapeLo     = 0.01;
constexpr double kWeibullShapeHi     = 50.0;
constexpr int    kWeibullBisectIter  = 100;
constexpr double kPiecewiseTolerance = 0.01;
constexpr int    kKsTerm             = 100;

inline double cdf(const       Normal& dist, double x) {return 0.5*std::erfc(-(x - dist.mean) / (dist.sigma*std::numbers::sqrt2));}
inline double cdf(const    Lognormal& dist, double x) {return x <= 0.0 ? 0.0 : 0.5*std::erfc(-(std::log(x) - dist.mu) / (dist.sigma*std::numbers::sqrt2));}
inline double cdf(const        Gamma& dist, double x) {return x <= 0.0 ? 0.0 : gamma_p(dist.shape, x / dist.scale);}
inline double cdf(const      Weibull& dist, double x) {return x <= 0.0 ? 0.0 : 1.0 - std::exp(-std::pow(x / dist.scale, dist.shape));}
inline double cdf(const  Exponential& dist, double x) {return x <= 0.0 ? 0.0 : 1.0 - std::exp(-dist.rate*x);}
inline double cdf(const   ChiSquared& dist, double x) {return x <= 0.0 ? 0.0 : gamma_p(dist.dof / 2.0, x / 2.0);}
inline double cdf(const ExtremeValue& dist, double x) {return std::exp(-std::exp(-(x - dist.location) / dist.scale));}
inline double cdf(const       Cauchy& dist, double x) {return 0.5 + std::atan((x - dist.location) / dist.scale) / k::pi;}

inline double cdf(const StudentT& dist, double x) {
  const double t = (x - dist.location) / dist.scale;
  const double ib = beta_p(0.5, dist.dof / 2.0, t*t / (t*t + dist.dof));
  return 0.5 + std::copysign(0.5*ib, t);
}

inline double cdf(const PiecewiseLinear& dist, double x) {
  const auto hi = std::ranges::upper_bound(dist.point, x, {}, &k::Point<>::x);
  double result = 0.0;
  if (hi == dist.point.begin()) result = dist.point.front().y;
  else if (hi == dist.point.end()) result = dist.point.back().y;
  else {
    const k::Point<>& a = *std::prev(hi);
    const k::Point<>& b = *hi;
    result = std::lerp(a.y, b.y, (x - a.x) / (b.x - a.x));
  }
  return result;
}

inline double cdf(const Distribution& dist, double x) {
  return std::visit([x](const auto& family) { return cdf(family, x); }, dist);
}

inline std::optional<Normal>
estimate(std::type_identity<Normal>, std::span<const double> v) {
  const double mean = average(v, std::identity{});
  const double sd   = std::sqrt(variance(v, std::identity{}, 0));
  return Normal{mean, sd};
}

inline std::optional<Lognormal>
estimate(std::type_identity<Lognormal>, std::span<const double> v) {
  std::optional<Lognormal> result;
  if (v.front() > 0.0) {
    std::vector<double> logged;
    for (const double x : v) logged.push_back(std::log(x));
    const double mu = average(logged, std::identity{});
    const double sd = std::sqrt(variance(logged, std::identity{}, 0));
    result = Lognormal{mu, sd};
  }
  return result;
}

inline std::optional<Gamma>
estimate(std::type_identity<Gamma>, std::span<const double> v) {
  std::optional<Gamma> result;
  if (v.front() > 0.0) {
    const double mean = average(v, std::identity{});
    const double var  = variance(v, std::identity{}, 0);
    result = Gamma{mean*mean / var, var / mean};
  }
  return result;
}

inline double weibull_cv(double shape) {
  const double g1 = std::tgamma(1.0 + 1.0 / shape);
  const double g2 = std::tgamma(1.0 + 2.0 / shape);
  return std::sqrt(g2 / (g1*g1) - 1.0);
}


inline std::optional<Weibull>
estimate(std::type_identity<Weibull>, std::span<const double> v) {
  std::optional<Weibull> result;
  if (v.front() > 0.0) {
    const double mean  = average(v, std::identity{});
    const double sd    = std::sqrt(variance(v, std::identity{}, 0));
    const double cv    = sd / mean;
    const double shape = bisect(weibull_cv, kWeibullShapeLo, kWeibullShapeHi, cv, kWeibullBisectIter);
    const double scale = mean / std::tgamma(1.0 + 1.0 / shape);
    result = Weibull{shape, scale};
  }
  return result;
}

inline std::optional<Exponential>
estimate(std::type_identity<Exponential>, std::span<const double> v) {
  return v.front() <= 0.0 ? std::nullopt
                          : std::optional{Exponential{1.0 / average(v, std::identity{})}};
}

inline std::optional<ChiSquared>
estimate(std::type_identity<ChiSquared>, std::span<const double> v) {
  return v.front() <= 0.0 ? std::nullopt
    : std::optional{ChiSquared{average(v, std::identity{})}};
}

inline std::optional<ExtremeValue>
estimate(std::type_identity<ExtremeValue>, std::span<const double> v) {
  const double mean  = average(v, std::identity{});
  const double sd    = std::sqrt(variance(v, std::identity{}, 0));
  const double scale = sd*std::sqrt(6.0) / k::pi;
  return ExtremeValue{mean - std::numbers::egamma*scale, scale};
}

inline std::optional<Cauchy>
estimate(std::type_identity<Cauchy>, std::span<const double> v) {
  return Cauchy{median(v), (percentile(v, 0.75) - percentile(v, 0.25)) / 2.0};
}

inline std::optional<StudentT>
estimate(std::type_identity<StudentT>, std::span<const double> v) {
  const double kurt = kurtosis(v, std::identity{});
  std::optional<StudentT> result;
  if (kurt > 0.0) {
    const double dof   = 4.0 + 6.0 / kurt;
    const double var   = variance(v, std::identity{}, 0);
    const double scale = std::sqrt(var * (dof - 2.0) / dof);
    result = StudentT{dof, median(v), scale};
  }
  return result;
}

inline std::optional<PiecewiseLinear>
estimate(std::type_identity<PiecewiseLinear>, std::span<const double> v) {
  const double n     = static_cast<double>(v.size());
  const double lo    = v.front();
  const double width = v.back() - lo;
  std::optional<PiecewiseLinear> result;
  if (width > 0.0) {
    std::vector<double> x(v.size());
    std::vector<double> y(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
      x[i] = (v[i] - lo) / width;
      y[i] = (static_cast<double>(i) + 0.5) / n;
    }
    const std::vector<size_t> vertex = rdp(x, y, kPiecewiseTolerance);
    std::vector<k::Point<>>   point(vertex.size());
    for (size_t i = 0; i < vertex.size(); ++i) point[i] = {v[vertex[i]], y[vertex[i]]};
    result = PiecewiseLinear{point};
  }
  return result;
}

inline std::vector<Distribution> estimate_all(std::span<const double> v) {
  std::vector<Distribution> candidate;
  if (auto dist = estimate(std::type_identity<Normal>{}, v))          candidate.push_back(*dist);
  if (auto dist = estimate(std::type_identity<Lognormal>{}, v))       candidate.push_back(*dist);
  if (auto dist = estimate(std::type_identity<Gamma>{}, v))           candidate.push_back(*dist);
  if (auto dist = estimate(std::type_identity<Weibull>{}, v))         candidate.push_back(*dist);
  if (auto dist = estimate(std::type_identity<Exponential>{}, v))     candidate.push_back(*dist);
  if (auto dist = estimate(std::type_identity<ChiSquared>{}, v))      candidate.push_back(*dist);
  if (auto dist = estimate(std::type_identity<ExtremeValue>{}, v))    candidate.push_back(*dist);
  if (auto dist = estimate(std::type_identity<Cauchy>{}, v))          candidate.push_back(*dist);
  if (auto dist = estimate(std::type_identity<StudentT>{}, v))        candidate.push_back(*dist);
  if (auto dist = estimate(std::type_identity<PiecewiseLinear>{}, v)) candidate.push_back(*dist);
  return candidate;
}

inline double ks_d(std::span<const double> v, const Distribution& dist) {
  const double n = static_cast<double>(v.size());
  double top = 0.0;
  for (size_t i = 0; i < v.size(); ++i) {
    const double F     = cdf(dist, v[i]);
    const double Fn_lo = static_cast<double>(i) / n;
    const double Fn_hi = static_cast<double>(i + 1) / n;
    top = std::max({top, std::abs(Fn_hi - F), std::abs(Fn_lo - F)});
  }
  return top;
}

inline double ks_p(size_t n, double d) {
  const double t = std::sqrt(static_cast<double>(n)) * d;
  double s = 0.0;
  const bool small = t < 1.0;
  for (int j = 1; j <= kKsTerm; ++j)
    s += small
       ? std::exp(-(2*j - 1.0) * (2*j - 1.0) * k::pi*k::pi / (8.0*t*t))
       : (j % 2 ? 1.0 : -1.0) * std::exp(-2.0 * j*j * t*t);
  return small
       ? 1.0 - std::sqrt(k::tau) / t*s
       : 2.0*s;
}

inline double fidelity(const Distribution& dist, std::span<const double> v) {
  return ks_p(v.size(), ks_d(v, dist));
}

}
