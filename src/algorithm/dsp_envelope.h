#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <span>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "algorithm/dsp_window.h"
#include "algorithm/math_scalar.h"

namespace qstorm::algorithm {

template <typename X, typename Y, typename Z>
double attack_time(const Spectrogram<X, Y, Z>& env, double lo, double hi) {
  const Quantization<X>& x_ms = env.x_ms;
  const Y*      rms = env.channel.front();
  const size_t peak = argmax(std::span<const Y>(rms, x_ms.extent));
  const double  top = rms[peak];

  size_t a_lo = 0;
  for (size_t i = 0; i < peak; i++)
    if (rms[i] >= lo + top) {
      a_lo = i;
      break;
    }

  size_t a_hi = peak;
  for (size_t i = a_lo; i < peak; i++)
    if (rms[i] >= hi + top) {
      a_hi = i;
      break;
    }

  const ms<X> x = env.begin();
  return static_cast<double>(a_hi - a_lo) * x.time_fraction();
}

template <typename X, typename Y, typename Z>
double release_time(const Spectrogram<X, Y, Z>& env, double lo, double hi) {
  const Quantization<X>& x_ms = env.x_ms;
  const Y*      rms = env.channel.front();
  const size_t peak = argmax(std::span<const Y>(rms, x_ms.extent));
  const double  top = rms[peak];

  size_t r_hi = peak;
  for (size_t i = peak; i < x_ms.extent; i++)
    if (rms[i] <= hi + top) {
      r_hi = i;
      break;
    }

  size_t r_lo = x_ms.extent - 1;
  for (size_t i = r_hi; i < x_ms.extent; i++)
    if (rms[i] <= lo + top) {
      r_lo = i;
      break;
    }
  const ms<X> x = env.begin();
  return static_cast<double>(r_lo - r_hi) * x.time_fraction();
}

template <typename X, typename Y, typename Z>
double sustain_level(const Spectrogram<X, Y, Z>& env, double hi) {
  const Quantization<X>& x_ms = env.x_ms;
  const Y*      rms = env.channel.front();
  const size_t peak = argmax(std::span<const Y>(rms, x_ms.extent));
  const double  top = rms[peak];

  size_t a_hi = peak;
  for (size_t i = 0; i < peak; i++)
    if (rms[i] >= hi + top) {
      a_hi = i;
      break;
    }

  size_t r_hi = peak;
  for (size_t i = peak; i < x_ms.extent; i++)
    if (rms[i] <= hi + top) {
      r_hi = i;
      break;
    }

  double total = 0.0;
  size_t n = 0;
  for (size_t i = a_hi; i < r_hi; i++) {
    total += rms[i];
    n++;
  }
  return total / static_cast<double>(n) - top;
}

template <typename X, typename Y, typename Z>
double noise_floor(const Spectrogram<X, Y, Z>& env, double lo, double hi) {
  const Quantization<X>& x_ms = env.x_ms;
  const Y*      rms = env.channel.front();
  const size_t peak = argmax(std::span<const Y>(rms, x_ms.extent));
  const double  top = rms[peak];

  size_t a_lo = 0;
  for (size_t i = 0; i < peak; i++)
    if (rms[i] >= lo + top) {
      a_lo = i;
      break;
    }
  size_t r_hi = peak;
  for (size_t i = peak; i < x_ms.extent; i++)
    if (rms[i] <= hi + top) {
      r_hi = i;
      break;
    }
  size_t r_lo = x_ms.extent - 1;
  for (size_t i = r_hi; i < x_ms.extent; i++)
    if (rms[i] <= lo + top) {
      r_lo = i;
      break;
    }

  double total = 0.0;
  size_t n = 0;
  for (size_t i = 0; i < a_lo; i++) {
    total += on_axis(rms[i], env.y_db);
    n++;
  }
  for (size_t i = r_lo + 1; i < x_ms.extent; i++) {
    total += on_axis(rms[i], env.y_db);
    n++;
  }
  return total / static_cast<double>(n) - top;
}

template <typename X, typename Y, typename Z>
std::vector<X> macro_point(const Spectrogram<X, Y, Z>& env, double smooth_ms, double tolerance) {
  std::vector<std::vector<Y>> smooth_channel;
  const Spectrogram<X, Y, Z>  smooth = moving_average(env, smooth_ms, smooth_channel);

  std::vector<X> time;
  std::vector<double> level;
  for (ms<X> x = smooth.begin(); x < smooth.end(); ++x) {
    time.push_back(x);
    level.push_back(smooth[x]);
  }
  const double      top = *std::ranges::max_element(level);
  const double    first = time.front();
  const double duration = time.back() - first;

  std::vector<double> x_norm;
  std::vector<double> y_norm;
  for (size_t i = 0; i < time.size(); ++i) {
    x_norm.push_back((time[i] - first) / duration);
    y_norm.push_back(level[i] - top);
  }

  const std::vector<size_t> vertex = rdp(x_norm, y_norm, tolerance);
  std::vector<X> x_vertex;
  for (size_t i : vertex) x_vertex.push_back(time[i]);
  return x_vertex;
}

template <typename X, typename Y, typename Z>
double residual_autocorr(const Spectrogram<X, Y, Z>& env, double smooth_ms) {
  std::vector<std::vector<Y>> smooth_channel;
  const Spectrogram<X, Y, Z>  smooth = moving_average(env, smooth_ms, smooth_channel);

  const double axis = env.y_db.offset;
  std::vector<double> residual;
  for (ms<X> x = env.begin(); x < env.end(); ++x)
    if (smooth[x] > axis && env[x] > axis) residual.push_back(env[x] - smooth[x]);

  const double mean = average(residual, std::identity{});
  double        num = 0.0;
  double        den = 0.0;
  for (size_t i = 0; i + 1 < residual.size(); i++)
    num += (residual[i] - mean) * (residual[i + 1] - mean);
  for (double v : residual) {
    const double d = v - mean;
    den += d*d;
  }
  const double r1 = num / den;
  return r1;
}

template <typename X, typename Y, typename Z>
std::vector<double> spike_interval(const Spectrogram<X, Y, Z>& env, double smooth_ms, double z) {
  std::vector<std::vector<Y>> smooth_channel;
  const Spectrogram<X, Y, Z>  smooth = moving_average(env, smooth_ms, smooth_channel);

  const double axis = env.y_db.offset;
  std::vector<double> ratio;
  std::vector<X> ratio_time;
  for (ms<X> x = env.begin(); x < env.end(); ++x) {
    ratio_time.push_back(x);
    if (smooth[x] > axis)
      ratio.push_back(env[x] - smooth[x]);
    else
      ratio.push_back(db<Y>::ratio(1.0));
  }
  const size_t n = ratio.size();

  const double mean = average(ratio, std::identity{});
  const double sd = std::sqrt(variance(ratio, std::identity{}, 0));
  const double threshold = mean + z*sd;

  std::vector<double> time;
  for (size_t i = 1; i + 1 < n; i++)
    if (ratio[i] > threshold && ratio[i] >= ratio[i - 1] && ratio[i] >= ratio[i + 1])
      time.push_back(ratio_time[i]);

  return difference(time, std::identity{});
}

}
