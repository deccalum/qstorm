#pragma once

#include <limits>
#include <numbers>
#include <vector>

namespace qstorm {
template <typename T = double>
using vec = std::vector<T>;
template <typename T = double>
using vvec = std::vector<std::vector<T>>;
template <typename T = double>
using vvvec = std::vector<std::vector<std::vector<T>>>;
}

namespace qstorm::k {
inline constexpr double pi /*         π                            */ = std::numbers::pi;
inline constexpr double tau /*        τ = 2π                       */ = std::numbers::pi * 2.0;
inline constexpr double hfpi /*       π / 2                        */ = std::numbers::pi / 2;
inline constexpr double sqrt2 /*                                   */ = std::numbers::sqrt2;
inline constexpr double F /*          20·log₁₀(DBL_MIN)            */ = -6153.053111371775;
inline constexpr double DbPerBit /*   20·log₁₀(2)                  */ = 6.020599913279624;
inline constexpr double QuantNoise /* 10·log₁₀(3/2)                */ = 1.760912590556813;
inline constexpr double Epsilon /*    ε, convergence tolerance     */ = 1e-12;
inline constexpr double Tiny /*       ε₀, floor for Lentz divisors */ = 1e-300;
inline constexpr double WeibullShapeLo                                = 0.01;
inline constexpr double WeibullShapeHi                                = 50.0;
inline constexpr int    WeibullBisectIter                             = 100;
inline constexpr double PiecewiseTolerance                            = 0.01;
inline constexpr int    KsTerm                                        = 100;

template <typename T = double>
inline constexpr T INF = std::numeric_limits<T>::infinity();

template <typename X = double, typename Y = double>
struct Point {
  X x;
  Y y;
};

struct Modifier {
  double (*fwd /* threshold search */)(double);
  double (*inv /* threshold back   */)(double);
};

inline constexpr Modifier Lin{[](double x) { return x; }, [](double x) { return x; }};

}  // namespace k, global constants
