#pragma once

#include <complex>
#include <cstddef>

#include "algorithm/axis_quantizer.h"
#include "algorithm/math_scalar.h"

namespace qstorm::algorithm {

template <typename X, typename Y, typename Z> Spectrogram<X, Y, Z>
make_mono(const Spectrogram<X, Y, Z>& spectrum, std::vector<std::vector<Y>>& channel) {
  const size_t count = spectrum.channel.size();
  const size_t sample = spectrum.x_ms.extent*spectrum.z_hz.extent;
  channel.assign(1, std::vector<Y>(sample));
  for (size_t i = 0; i < sample; i++) {
    Y total = 0;
    for (size_t j = 0; j < count; j++) total += spectrum.channel[j][i];
    channel.front()[i] = total / static_cast<Y>(count);
  }
  return as_spectrogram(channel, spectrum.x_ms, spectrum.z_hz, spectrum.y_db);
}

template <typename X, typename Y, typename Z> std::vector<double>
hann_smoothing(const Spectrogram<X, Y, Z>& spectrogram, bool periodic = true) {
  ms<X> last = spectrogram.end();
  --last;
  const X denominator = periodic ? spectrogram.end() - spectrogram.begin()
                                 : last - spectrogram.begin();

  std::vector<double> shape;
  for (ms<X> x = spectrogram.begin(); x < spectrogram.end(); ++x) {
    const double elapsed = x - spectrogram.begin();
    const double theta   = k::tau*elapsed / denominator;
    const double hann    = 0.5 * (1.0 - std::cos(theta));
    shape.push_back(hann);
  }
  return shape;
}

template <typename X>
void fft(std::span<std::complex<X>> x) {
  const size_t n = x.size();
  for (size_t i = 1, j = 0; i < n; i++) {
    size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(x[i], x[j]);
  }
  for (size_t w = 2; w <= n; w <<= 1) {
    const size_t half = w / 2;
    const X angle = -2*std::numbers::pi_v<X> / static_cast<X>(w);
    const std::complex<X> step = std::polar(X{1}, angle);
    for (size_t i = 0; i < n; i += w) {
      std::complex<X> twiddle = 1;
      for (size_t j = 0; j < half; j++) {
        const std::complex<X> even = x[i + j];
        const std::complex<X> odd = x[i + j + half] * twiddle;
        x[i + j] = even + odd;
        x[i + j + half] = even - odd;
        twiddle *= step;
      }
    }
  }
}

template <typename X, typename Y, typename Z> Spectrogram<X, Y, Z>
stft(const Spectrogram<X, Y, Z> env, X fft_ms, X hop_ms, std::vector<std::vector<Y>>& channel) {
  const Spectrogram<X, Y, Z> reach_fft = env.window(env.begin(), env.begin() + fft_ms);
  const size_t                fft_size = reach_fft.x_ms.extent;
  const X                     fft_snap = reach_fft.end() - reach_fft.begin();
  const Spectrogram<X, Y, Z> reach_hop = env.window(env.begin(), env.begin() + hop_ms);
  const size_t                hop_size = reach_hop.x_ms.extent;
  const X                     hop_snap = reach_hop.end() - reach_hop.begin();

  const size_t pad = std::bit_ceil(fft_size);
  const size_t hop_count = env.x_ms.extent < fft_size ? 0 : (env.x_ms.extent - fft_size) / hop_size + 1;
  const Quantization<X> x_ms = Quantization<X>::stepped(hop_count, hop_snap, env.begin());
  const size_t hertz = env.x_ms.frame(X{1000});
  const Quantization<Z> z_hz{pad / 2 + 1, ratio<Z>(pad, hertz), Z{0}};

  const std::vector<double> shape = hann_smoothing(reach_fft);
  const double        total = sum(shape, std::identity{});
  const double norm_general = total / 2.0;
  const double    norm_edge = total;

  std::vector<double> magnitude;
  double top = 0.0;
  ms<X> x0 = env.begin();
  for (size_t i = 0; i < x_ms.extent; ++i) {
    const Spectrogram<X, Y, Z> slice = env.window(x0, x0 + fft_snap);
    std::vector<std::complex<double>> dft;
    for (ms<X> x = slice.begin(); x < slice.end(); ++x)
      dft.push_back(static_cast<double>(slice[x]));
    for (size_t j = 0; j < fft_size; ++j) dft[j] *= shape[j];
    dft.resize(pad);
    fft(std::span<std::complex<double>>(dft));

    for (size_t k = 0; k < z_hz.extent; k++) {
      const double norm = (k == 0 || k == pad / 2) ? norm_edge : norm_general;
      magnitude.push_back(std::abs(dft[k]) / norm);
      top = std::max(top, magnitude.back());
    }
    x0 = x0 + hop_snap;
  }

  channel.assign(1, std::vector<Y>());
  for (double e : magnitude) channel.front().push_back(db<Y>::ratio(std::fabs(e) / top));

  return as_spectrogram(channel, x_ms, z_hz, y_db(channel, env.y_db));
}

template <typename X, typename Y, typename Z>
double spectral_centroid(const Spectrogram<X, Y, Z>& spectrum) {
  const double axis = spectrum.y_db.offset;
  const ms<X> x = spectrum.begin();
  double weight = 0.0;
  double moment = 0.0;
  for (hz<Z> z = spectrum.z_begin(); z < spectrum.z_end(); ++z) {
    const double height = spectrum[x, z] - axis;
    weight += height;
    moment += z*height;
  }
  return moment / weight;
}

template <typename X, typename Y, typename Z>
double spectral_spread(const Spectrogram<X, Y, Z>& spectrum, double centroid) {
  const double axis = spectrum.y_db.offset;
  const ms<X>     x = spectrum.begin();
  double     weight = 0.0;
  double dispersion = 0.0;
  for (hz<Z> z = spectrum.z_begin(); z < spectrum.z_end(); ++z) {
    const double height = spectrum[x, z] -axis;
    const double d = z -centroid;
    weight += height;
    dispersion += d*d*height;
  }
  return std::sqrt(dispersion / weight);
}

template <typename X, typename Y, typename Z>
double spectral_tilt(const Spectrogram<X, Y, Z>& spectrum) {
  const Quantization<Z>& z_hz = spectrum.z_hz;
  const Y* mag = spectrum.channel.front();
  std::vector<double> octave;
  std::vector<double> level;
  hz<Z> z = spectrum.z_begin();
  for (size_t i = 1; i < z_hz.extent; i++) {
    ++z;
    octave.push_back(std::log2(static_cast<double>(z)));
    level.push_back(mag[i]);
  }
  return slope(octave, level);
}

template <typename X, typename Y, typename Z>
double spectral_flatness(const Spectrogram<X, Y, Z>& spectrum) {
  const double   axis = spectrum.y_db.offset;
  const        Y* mag = spectrum.channel.front();
  const size_t extent = spectrum.z_hz.extent;
  double    log_total = 0.0;
  double        total = 0.0;
  size_t            n = 0;
  for (size_t i = 0; i < extent; ++i) {
    if (mag[i] <= axis) continue;
    const double height = mag[i] - axis;
    log_total += std::log(height);
    total += height;
    n++;
  }
  const double geometric  = std::exp(log_total / static_cast<double>(n));
  const double arithmetic = total / static_cast<double>(n);
  return geometric / arithmetic;
}

template <typename X, typename Y, typename Z>
double spectral_rolloff(const Spectrogram<X, Y, Z>& spectrum, double fraction) {
  const Quantization<Z>& z_hz = spectrum.z_hz;
  const double axis = spectrum.y_db.offset;
  const Y* mag = spectrum.channel.front();
  std::vector<double> cumulative(z_hz.extent + 1);
  cumulative[0] = 0.0;
  for (size_t i = 0; i < z_hz.extent; ++i)
    cumulative[i + 1] = cumulative[i] + (mag[i] - axis);
  const double target = fraction*cumulative.back();
  size_t rolloff_bin = z_hz.extent - 1;
  for (size_t i = 0; i < z_hz.extent; ++i)
    if (cumulative[i + 1] >= target) {
      rolloff_bin = i;
      break;
    }
  return spectrum.z_grid(spectrum.first() + rolloff_bin).value;
}

}
