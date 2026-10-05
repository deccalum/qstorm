#pragma once

#include <numbers>

namespace qstorm::k {
inline double pi = std::numbers::pi;                     // π
inline double tau = std::numbers::pi*2.0;                // τ = 2π
inline constexpr double F = -6153.053111371775;          // 20·log₁₀(DBL_MIN)
inline constexpr double DbPerBit = 6.020599913279624;    // 20·log₁₀(2)
inline constexpr double QuantNoise = 1.760912590556813;  // 10·log₁₀(3/2)
inline constexpr double Epsilon = 1e-12;                 // ε, convergence tolerance
inline constexpr double Tiny = 1e-300;                   // ε₀, floor for Lentz divisors

template <typename X = double, typename Y = double>
struct Point {X x; Y y;};
} // namespace k, global constants
