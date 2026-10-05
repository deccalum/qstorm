#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <functional>
#include <limits>
#include <ranges>
#include <span>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "algorithm/dsp_envelope.h"
#include "algorithm/dsp_window.h"
#include "algorithm/math_constant.h"
#include "algorithm/math_scalar.h"

namespace qstorm::algorithm {

struct FourierSeries {
  int stride;
  int power;
  bool alternative;
};

inline constexpr FourierSeries     kSine{1, 0, false};
inline constexpr FourierSeries      kSaw{1, 1, false};
inline constexpr FourierSeries   kSquare{2, 1, false};
inline constexpr FourierSeries kTriangle{2, 2,  true};

struct Arc {
  FourierSeries fourier;  //
  double ms;              //
  int low;                // first harmonic boundary
  int high;               // last  harmonic boundary
  double amp;             //
  double offset;          //
  double phase;           //
};

struct Locus {std::vector<Arc> path;};

struct Sigmoid {
  double k;
  double x0;
};

struct Term {
  int k;
  double weight;
  double phase;
};

struct Square {
  static constexpr Term term(int n) {
    return {
      2*n - 1,
      1.0 / (2*n - 1),
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
      2*n - 1,
      1.0 / ((2*n - 1) * (2*n - 1)),
      k::pi * ((n + 1) % 2),
    };
  }
};

struct Sine {
  static constexpr Term term(int n) {
    return {
      n,
      1.0,
      0.0
    };
  }
};

inline int harmonic(const FourierSeries& f, int i) {
  return 1 + f.stride * (i - 1);
}

inline double weight(const FourierSeries& f, int i) {
  const int j = harmonic(f, i);
  const double sign = f.alternative && i % 2 == 0 ? -1.0 : 1.0;
  return sign / std::pow(j, f.power);
}

inline double wave(const Arc& a, double x) {
  double total = 0.0;
  for (int i = a.low; i <= a.high; ++i)
    total += weight(a.fourier, i) * std::sin(k::tau * harmonic(a.fourier, i) * (x + a.phase));
  return a.offset + a.amp*total;
}

inline double sigmoid(const Sigmoid& s, double x) {
  return 1.0 / (1.0 + std::exp(-s.k * (x - s.x0)));
}

template <typename X, typename Y, typename Z>
std::vector<std::complex<double>> coefficient(const Spectrogram<X, Y, Z>& env) {
  const Quantization<X>& x_ms = env.x_ms;
  const std::span<const Y> r(env.channel.front(), x_ms.extent);
  const double mean = average(r, std::identity{});
  const double n = static_cast<double>(x_ms.extent);
  std::vector<std::complex<double>> c(x_ms.extent / 2 + 1);
  for (size_t i = 0; i < c.size(); i++) {
    std::complex<double> total = 0.0;
    for (size_t j = 0; j < x_ms.extent; j++) {
      const double x = static_cast<double>(j) / n;
      const double angle = k::tau*static_cast<double>(i)*x;
      total += (r[j] - mean) * std::polar(1.0, angle);
    }
    c[i] = 2.0*total / n;
  }
  return c;
}

inline std::array<double, 2>
projection(std::span<const std::complex<double>> c, const FourierSeries& f, int lo, int hi, double phase) {
  double dot = 0.0;
  double energy = 0.0;
  for (int i = lo; i <= hi; ++i) {
    const int j = harmonic(f, i);
    const double w = weight(f, i);
    dot += w*std::imag(c[j] * std::polar(1.0, k::tau*j*phase));
    energy += w*w;
  }
  return {dot / energy, dot*dot / energy};
}

template <typename X>
std::vector<double> arc_time(const Arc& a, const Quantization<X>& x_ms) {
  std::vector<double> value(x_ms.extent);
  ms<X> x{x_ms.offset, &x_ms};
  for (size_t i = 0; i < x_ms.extent; i++, ++x) {
    const double elapsed = static_cast<double>(x) - static_cast<double>(x_ms.offset);
    value[i] = wave(a, elapsed / a.ms);
  }
  return value;
}

template <typename X, typename Y, typename Z>
Arc arc_spectrum(const Spectrogram<X, Y, Z>& env, std::span<const FourierSeries> f) {
  const Quantization<X>& x_ms = env.x_ms;
  const int n = static_cast<int>(x_ms.extent);
  const double x = static_cast<double>(n) / x_ms.density;
  const double offset = average(std::span<const Y>(env.channel.front(), x_ms.extent), std::identity{});
  const std::vector<std::complex<double>> c = coefficient(env);

  Arc best{f.front(), x, 1, 1, 0.0, offset, 0.0};
  double best_energy = -1.0;
  for (const FourierSeries fs : f)
    for (int i = 1; harmonic(fs, i) <= n / 2; i++)
      for (int j = i; harmonic(fs, j) <= n / 2; j++)
        for (int k = 0; k < n; k++) {
          const double phase = static_cast<double>(k) / n;
          const std::array<double, 2> proj = projection(c, fs, i, j, phase);
          const double amp = proj[0];
          const double energy = proj[1];
          if (energy > best_energy) {
            best_energy = energy;
            best = {fs, x, i, j, amp, offset, phase};
          }
        }
  return best;
}

template <typename X, typename Y, typename Z>
Locus locus(const Spectrogram<X, Y, Z>& env, std::span<const FourierSeries> f, double smooth_ms,
            double tolerance) {
  const std::vector<X> x_vertex = macro_point(env, smooth_ms, tolerance);
  Locus l;
  for (auto [prev, next] : x_vertex | std::views::adjacent<2>)
    l.path.push_back(arc_spectrum(env.window(prev, next), f));
  return l;
}

double residual(auto wave, std::span<const double> t, int lo, int hi) {
  std::vector<double> synth(t.size());
  arc(wave, synth, 1.0, lo, hi);
  double total = 0.0;
  for (size_t i = 0; i < t.size(); ++i) {
    const double d = t[i] - synth[i];
    total += d*d;
  }
  return total;
}

template <typename Wave>
void sweep(Wave wave, std::span<const double> t, size_t peak, double& best_err, int& best_lo,
           int& best_hi) {
  for (int i = 1; i <= static_cast<int>(peak); ++i)
    for (int j = i; j <= static_cast<int>(peak); ++j) {
      const double err = residual(wave, t, i, j);
      if (err < best_err) {
        best_err = err;
        best_lo = i;
        best_hi = j;
      }
    }
}

template <typename T>
Arc wave_solution(std::span<const T> t) {
  const size_t peak = centroid_index(t);
  double best_err = std::numeric_limits<double>::max();
  int best_lo = 1, best_hi = 1;

  sweep(Square{}, t, peak, best_err, best_lo, best_hi);
  sweep(Saw{}, t, peak, best_err, best_lo, best_hi);
  sweep(Triangle{}, t, peak, best_err, best_lo, best_hi);

  return {0.0, *std::ranges::max_element(t), best_lo, best_hi};
}

}
