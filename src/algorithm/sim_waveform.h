#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <functional>
#include <ranges>
#include <span>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "algorithm/dsp_envelope.h"
#include "algorithm/math_constant.h"
#include "algorithm/math_scalar.h"

namespace qstorm::algorithm {
struct FourierSeries {
  int  stride;
  int  power;
  bool alternative;
};

inline constexpr FourierSeries kSine{1, 0, false};
inline constexpr FourierSeries kSaw{1, 1, false};
inline constexpr FourierSeries kSquare{2, 1, false};
inline constexpr FourierSeries kTriangle{2, 2, true};

struct Arc {
  FourierSeries fourier;
  double        ms;
  int           low;
  int           high;
  double        amp;
  double        pos;
  double        phase;
};

struct Locus {
  vec<Arc> path;
};

struct Sigmoid {
  double k;
  double x0;
};

struct Term {
  int    k;
  double weight;
  double phase;
};

struct Square {
  static constexpr Term term(int n) {
    return {
      2 * n - 1,
      1.0 / (2 * n - 1),
      0.0,
    };
  }
};

struct Saw {
  static constexpr Term term(int n) {
    return {
      n,
      1.0 / n,
      0.0,
    };
  }
};

struct Triangle {
  static constexpr Term term(int n) {
    return {
      2 * n - 1,
      1.0 / ((2 * n - 1) * (2 * n - 1)),
      k::pi * ((n + 1) % 2),
    };
  }
};

struct Sine {
  static constexpr Term term(int n) { return {n, 1.0, 0.0}; }
};

inline int harmonic(const FourierSeries& fourier, int i) { return 1 + fourier.stride * (i - 1); }

inline double weight(const FourierSeries& fourier, int i) {
  const int    j    = harmonic(fourier, i);
  const double sign = fourier.alternative && i % 2 == 0 ? -1.0 : 1.0;
  return sign / std::pow(j, fourier.power);
}

inline double wave(const Arc& a, double x) {
  double tot = 0.0;
  for (int i = a.low; i <= a.high; ++i)
    tot += weight(a.fourier, i) * std::sin(k::tau * harmonic(a.fourier, i) * (x + a.phase));
  return a.pos + a.amp * tot;
}

inline double sigmoid(const Sigmoid& s, double x) {
  return 1.0 / (1.0 + std::exp(-s.k * (x - s.x0)));
}

template <typename X, typename Y, typename Z>
vec<std::complex<double>> coefficient(const Spectrum<X, Y, Z>& env) {
  const Qtz<X>&            qx = env.qx;
  const std::span<const Y> r(env.ch.front(), qx.n);
  const double             mean = average(r, std::identity{});
  const double             n    = qx.n;

  vec<std::complex<double>> c(qx.n / 2 + 1);
  for (size_t i = 0; i < c.size(); i++) {
    std::complex<double> tot = 0.0;
    for (size_t j = 0; j < qx.n; j++) {
      const double phase = j / n;
      const double angle = k::tau * i * phase;
      tot += (r[j] - mean) * std::polar(1.0, angle);
    }
    c[i] = 2.0 * tot / n;
  }
  return c;
}

inline std::array<double, 2>
projection(std::span<const std::complex<double>> c, const FourierSeries& fourier, int i0, int i1,
           double phase) {
  double dot = 0.0;
  double nrg = 0.0;
  for (int i = i0; i <= i1; ++i) {
    const int    j = harmonic(fourier, i);
    const double w = weight(fourier, i);
    dot += w * std::imag(c[j] * std::polar(1.0, k::tau * j * phase));
    nrg += w * w;
  }
  return {dot / nrg, dot * dot / nrg};
}

template <typename X>
vec<> arc_time(const Arc& a, const Qtz<X>& qx) {
  vec   value(qx.n);
  ms<X> x{qx.pos, &qx};
  for (size_t i = 0; i < qx.n; i++, ++x) {
    const double x0 = qx.pos;
    const double dx = x - x0;
    value[i]        = wave(a, dx / a.ms);
  }
  return value;
}

template <typename X, typename Y, typename Z>
Arc arc_spectrum(const Spectrum<X, Y, Z>& env, std::span<const FourierSeries> fourier) {
  const Qtz<X>& qx  = env.qx;
  const int     nx  = qx.n;
  const double  n   = nx;
  const double  dx  = n / qx.rho;
  const double  pos = average(std::span<const Y>(env.ch.front(), qx.n), std::identity{});
  const vec<std::complex<double>> c = coefficient(env);

  Arc    arc{fourier.front(), dx, 1, 1, 0.0, pos, 0.0};
  double y1 = -1.0;
  for (const FourierSeries each : fourier)
    for (int i = 1; harmonic(each, i) <= nx / 2; i++)
      for (int j = i; harmonic(each, j) <= nx / 2; j++)
        for (int m = 0; m < nx; m++) {
          const double                phase  = m / n;
          const std::array<double, 2> proj   = projection(c, each, i, j, phase);
          const double                amp    = proj[0];
          const double                energy = proj[1];
          if (energy > y1) {
            y1  = energy;
            arc = {each, dx, i, j, amp, pos, phase};
          }
        }
  return arc;
}

template <typename X, typename Y, typename Z>
Locus locus(const Spectrum<X, Y, Z>& env, std::span<const FourierSeries> fourier, double smooth_ms,
            double tolerance) {
  const vec<X> x_vertex = macro_point(env, smooth_ms, tolerance);
  Locus        l;
  for (auto [prev, next] : x_vertex | std::views::adjacent<2>)
    l.path.push_back(arc_spectrum(env.window(prev, next), fourier));
  return l;
}

double residual(auto wave, std::span<const double> t, int i0, int i1) {
  vec synth(t.size());
  arc(wave, synth, 1.0, i0, i1);
  double tot = 0.0;
  for (size_t i = 0; i < t.size(); ++i) {
    const double d = t[i] - synth[i];
    tot += d * d;
  }
  return tot;
}

template <typename Wave>
void sweep(Wave wave, std::span<const double> t, size_t peak, double& best_err, int& best_lo,
           int& best_hi) {
  const int highest = peak;
  for (int i = 1; i <= highest; ++i)
    for (int j = i; j <= highest; ++j) {
      const double err = residual(wave, t, i, j);
      if (err < best_err) {
        best_err = err;
        best_lo  = i;
        best_hi  = j;
      }
    }
}

template <typename T>
Arc wave_solution(std::span<const T> t) {
  const size_t peak     = centroid_index(t);
  double       best_err = k::INF<>;
  int          best_lo = 1, best_hi = 1;

  sweep(Square{}, t, peak, best_err, best_lo, best_hi);
  sweep(Saw{}, t, peak, best_err, best_lo, best_hi);
  sweep(Triangle{}, t, peak, best_err, best_lo, best_hi);

  return {0.0, *std::ranges::max_element(t), best_lo, best_hi};
}
}
