#pragma once

#include <algorithm>
#include <bit>
#include <cmath>
#include <complex>
#include <cstddef>
#include <functional>
#include <optional>
#include <utility>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "algorithm/dsp_window.h"
#include "algorithm/math_constant.h"
#include "algorithm/math_scalar.h"

namespace qstorm::algorithm {
template <typename X, typename Y, typename Z>
Spectrum<X, Y, Z> log_mono(const Spectrum<X, Y, Z>& env, vvec<Y>& ch) {
  const size_t nc = env.ch.size();
  const size_t n  = env.qx.n * env.qz.n;
  ch.assign(1, vec<Y>(n));
  for (size_t i = 0; i < n; i++) {
    Y tot = 0;
    for (size_t c = 0; c < nc; c++) tot += env.ch[c][i];
    ch.front()[i] = tot / nc;
  }
  return as_spectrum(ch, env.qx, env.qy, env.qz);
}

template <typename X, typename Y, typename Z>
Spectrum<X, Y, Z> lin_mono(const Spectrum<X, Y, Z>& env, vvec<Y>& ch) {
  const size_t nc = env.ch.size();
  const size_t n  = env.qx.n * env.qz.n;
  ch.assign(1, vec<Y>(n));
  for (size_t i = 0; i < n; i++) {
    Y tot = 0;
    for (size_t c = 0; c < nc; c++) tot += db<Y>::lin(env.ch[c][i]);
    const Y pwr   = tot / nc;
    ch.front()[i] = db<Y>::ratio(std::sqrt(pwr));
  }
  return as_spectrum(ch, env.qx, env.qy, env.qz);
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
    const size_t          hf    = w / 2;
    const X               angle = -2 * std::numbers::pi_v<X> / w;
    const std::complex<X> step  = std::polar(X{1}, angle);
    for (size_t i = 0; i < n; i += w) {
      std::complex<X> twiddle = 1;
      for (size_t j = 0; j < hf; j++) {
        const std::complex<X> even = x[i + j];
        const std::complex<X> odd  = x[i + j + hf] * twiddle;
        x[i + j]                   = even + odd;
        x[i + j + hf]              = even - odd;
        twiddle *= step;
      }
    }
  }
}

template <typename X, typename Z>
size_t fft_size(const Qtz<X>& qx, const Qtz<Z>& qz) {
  return std::round(qz.rho * qx.rho);
}

template <typename X, typename Y, typename Z>
vvec<std::complex<double>>
frame_spectrum(const Spectrum<X, Y, Z>& env, const Qtz<X>& qx, const Qtz<Z>& qz) {
  const size_t n    = fft_size(env.qx, qz);
  const vec    hann = periodic_hann(n);
  const size_t hf   = n / 2;
  const size_t nz   = hf + 1;
  const size_t hop  = env.qx.frame(X{1} / qx.rho);
  const size_t nc   = env.ch.size();

  vvec<std::complex<double>> spectrum(nc, vec<std::complex<double>>(qx.n * nz));
  for (size_t i = 0; i < qx.n; ++i) {
    const size_t center = i * hop;
    for (size_t c = 0; c < nc; ++c) {
      vec<std::complex<double>> dft(n);
      for (size_t j = 0; j < n; ++j) {
        const size_t at     = center + j;
        const bool   in     = at >= hf && at - hf < env.qx.n;
        const double sample = in ? env.ch[c][at - hf] : 0.0;
        dft[j]              = sample * hann[j];
      }
      fft(std::span<std::complex<double>>(dft));
      for (size_t k = 0; k < nz; ++k) spectrum[c][i * nz + k] = dft[k];
    }
  }
  return spectrum;
}

template <typename X, typename Y, typename Z>
vvec<> frame_power(const Spectrum<X, Y, Z>& env, const Qtz<X>& qx, const Qtz<Z>& qz) {
  const size_t                     n            = fft_size(env.qx, qz);
  const vec                        hann         = periodic_hann(n);
  const double                     tot          = sum(hann, std::identity{});
  const double                     norm_general = tot / 2.0;
  const double                     norm_edge    = tot;
  const size_t                     half         = n / 2;
  const size_t                     nz           = half + 1;
  const vvec<std::complex<double>> z            = frame_spectrum(env, qx, qz);
  vvec                             power(z.size(), vec(qx.n * nz));
  for (size_t c = 0; c < z.size(); ++c)
    for (size_t i = 0; i < qx.n; ++i)
      for (size_t k = 0; k < nz; ++k) {
        const double norm    = (k == 0 || k == half) ? norm_edge : norm_general;
        const double amp     = std::abs(z[c][i * nz + k]) / norm;
        power[c][i * nz + k] = amp * amp;
      }
  return power;
}

template <typename X, typename Y, typename Z>
Qtz<X> band_frame(const Spectrum<X, Y, Z>& env, const Qtz<Z>& qz) {
  const size_t n  = fft_size(env.qx, qz);
  const size_t hf = n / 2;
  const X      dx = hf / env.qx.rho;
  const size_t nx = (env.qx.n + hf - 1) / hf + 1;
  return Qtz<X>::stepped(nx, dx, env.x0());
}

template <typename Z>
double crossover_step(const vec<Qtz<Z>>& qz, size_t boundary, double z) {
  const Qtz<Z>& next     = qz[boundary + 1];
  const double  mainlobe = 4.0;
  const double  edge     = next.pos;
  const double  dz       = mainlobe / next.rho;
  const double  t        = std::clamp((z - edge) / dz + 0.5, 0.0, 1.0);
  return 0.5 - 0.5 * std::cos(k::pi * t);
}

template <typename Z>
double band_weight(const vec<Qtz<Z>>& qz, size_t b, double z) {
  const double below = b == 0 ? 1.0 : crossover_step(qz, b - 1, z);
  const double above = b + 1 == qz.size() ? 0.0 : crossover_step(qz, b, z);
  return below - above;
}

template <typename X, typename Y, typename Z>
vvvec<Y> band_signal(const Spectrum<X, Y, Z>& env, const vec<Qtz<Z>>& qz) {
  const size_t nx = env.qx.n;
  const size_t nc = env.ch.size();
  const size_t n  = std::bit_ceil(2 * nx);

  vvvec<Y> part(qz.size(), vvec<Y>(nc, vec<Y>(nx)));
  for (size_t c = 0; c < nc; ++c) {
    vec<std::complex<double>> whole(n);
    for (size_t i = 0; i < nx; ++i) whole[i] = env.ch[c][i];
    fft(std::span<std::complex<double>>(whole));
    for (size_t b = 0; b < qz.size(); ++b) {
      vec<std::complex<double>> filtered(n);
      for (size_t k = 0; k < n; ++k) {
        const size_t mirror = k <= n / 2 ? k : n - k;
        const double z      = ratio<double>(mirror, n) * env.qx.rho;
        filtered[k]         = std::conj(whole[k] * band_weight(qz, b, z));
      }
      fft(std::span<std::complex<double>>(filtered));
      for (size_t i = 0; i < nx; ++i) part[b][c][i] = filtered[i].real() / n;
    }
  }
  return part;
}

template <typename X, typename Y, typename Z>
vvvec<std::complex<double>> band_spectrum(const Spectrum<X, Y, Z>& env, const vec<Qtz<Z>>& qz) {
  vvvec<Y>                    part = band_signal(env, qz);
  vvvec<std::complex<double>> spectrum;
  for (size_t b = 0; b < qz.size(); ++b) {
    const Spectrum<X, Y, Z> band = as_spectrum(part[b], env.qx, env.qy, env.qz);
    spectrum.push_back(frame_spectrum(band, band_frame(env, qz[b]), qz[b]));
  }
  return spectrum;
}

template <typename Z>
vec<size_t> mask_row(const vec<Qtz<Z>>& qz, size_t b, size_t nz) {
  vec<size_t> i0{0};
  for (const Qtz<Z>& each : qz) i0.push_back(i0.back() + each.n);
  vec<size_t> row(nz);
  for (size_t k = 0; k < nz; ++k) {
    const Qtz<Z>& band  = qz[b];
    const double  rho   = band.rho;
    const double  z     = k / rho;
    size_t        bz    = 0;
    while (bz + 1 < qz.size() && qz[bz + 1].pos <= z) ++bz;
    const Qtz<Z>& qbz     = qz[bz];
    const double  u       = std::round((z - qbz.pos) * qbz.rho);
    const double  i1      = qbz.n - 1;
    const size_t  band_in = std::clamp(u, 0.0, i1);
    row[k]                = i0[bz] + band_in;
  }
  return row;
}

template <typename X, typename Y, typename Z>
vvec<> inverse_stft_masked(const Spectrum<X, Y, Z>& env, const Qtz<X>& qx, const vec<Qtz<Z>>& qz,
                           const vvvec<std::complex<double>>& spectrum, const vvec<>& mask) {
  const size_t nc  = env.ch.size();
  const size_t hop = env.qx.frame(X{1} / qx.rho);
  vvec         tot(nc, vec(env.qx.n));
  for (size_t b = 0; b < qz.size(); ++b) {
    const size_t      n     = fft_size(env.qx, qz[b]);
    const size_t      half  = n / 2;
    const size_t      nz    = half + 1;
    const size_t      nx    = band_frame(env, qz[b]).n;
    const vec         hann  = periodic_hann(n);
    const vec<size_t> bz    = mask_row(qz, b, nz);
    vvec              part(nc, vec(env.qx.n));
    vec               coverage(env.qx.n);

    for (size_t i = 0; i < nx; ++i) {
      const double u  = ratio<double>(i * half, hop);
      const size_t i0 = std::min<size_t>(u, qx.n - 1);
      const size_t i1 = std::min(i0 + 1, qx.n - 1);
      const double t  = std::min(u - i0, 1.0);
      vec          gain(nz);
      for (size_t k = 0; k < nz; ++k) gain[k] = (1.0 - t) * mask[bz[k]][i0] + t * mask[bz[k]][i1];
      const bool silent = *std::ranges::max_element(gain) == 0.0;
      for (size_t c = 0; c < nc && !silent; ++c) {
        vec<std::complex<double>> dft(n);
        for (size_t k = 0; k < nz; ++k) dft[k] = std::conj(spectrum[b][c][i * nz + k] * gain[k]);
        for (size_t k = 1; k < half; ++k) dft[n - k] = std::conj(dft[k]);
        fft(std::span<std::complex<double>>(dft));
        for (size_t j = 0; j < n; ++j) {
          const size_t at = i * half + j;
          if (at < half || at - half >= env.qx.n) continue;
          part[c][at - half] += dft[j].real() / n;
        }
      }
      for (size_t j = 0; j < n; ++j) {
        const size_t at = i * half + j;
        if (at < half || at - half >= env.qx.n) continue;
        coverage[at - half] += hann[j];
      }
    }
    for (size_t c = 0; c < nc; ++c)
      for (size_t i = 0; i < env.qx.n; ++i) tot[c][i] += part[c][i] / coverage[i];
  }
  return tot;
}

template <typename X, typename Y, typename Z>
vec<Qtz<Z>> stft_fft_ms(const Spectrum<X, Y, Z>& env) {
  const size_t nmax = std::bit_floor(env.qx.n);
  const size_t nz   = nmax / 2 + 1;
  const Z      zmin = ratio<Z>(env.qx.rho, nmax);
  const Z      zmax = ratio<Z>(nz * env.qx.rho, nmax);
  const double dx   = ratio<double>(env.qx.n, env.qx.rho);
  const double dz   = zmax - zmin;
  const double x    = std::sqrt(dx / dz);
  const size_t nx   = x * env.qx.rho;
  const size_t nmin = std::bit_floor(nx);
  vec<Qtz<Z>>  band;
  Z            z1 = zmax;
  for (size_t n = nmin; z1 > zmin; n = std::min(2 * n, nmax)) {
    const Z      z0  = n == nmax ? zmin : std::max(zmin, z1 / 2);
    const size_t k0  = std::ceil(z0 * n / env.qx.rho);
    const size_t k1  = std::ceil(z1 * n / env.qx.rho);
    const Z      rho = ratio<Z>(n, env.qx.rho);
    const Z      pos = k0 / rho;
    band.push_back(Qtz<Z>{k1 - k0, rho, pos});
    z1 = z0;
  }
  std::ranges::reverse(band);
  return band;
}

template <typename X, typename Y, typename Z>
Qtz<X> stft_hop_ms(const Spectrum<X, Y, Z>& env, const Qtz<Z>& qz) {
  const size_t n   = fft_size(env.qx, qz);
  const size_t hop = n / 2;
  const size_t nx  = (env.qx.n + hop - 1) / hop;
  const X      dx  = hop / env.qx.rho;
  return Qtz<X>::stepped(nx, dx, env.x0());
}

template <typename X, typename Y, typename Z>
zvec<X, Y, Z>
stft(const Spectrum<X, Y, Z>& env, const Qtz<X>& qx, const vec<Qtz<Z>>& qz, vvvec<Y>& ch) {
  const size_t hop = env.qx.frame(X{1} / qx.rho);
  const size_t nc  = env.ch.size();
  ch.assign(qz.size(), vvec<Y>(nc + 1));

  zvec<X, Y, Z> band;
  for (size_t b = 0; b < qz.size(); ++b) {
    const Qtz<Z>& qzb = qz[b];
    const size_t  n   = fft_size(env.qx, qz[b]);
    const size_t  hf  = n / 2;
    const size_t  nz  = hf + 1;
    const size_t  k0  = std::round(qzb.pos * qzb.rho);
    const Qtz<X>  qxb = band_frame(env, qz[b]);
    const size_t  nx  = qxb.n;
    const vvec    pwr = frame_power(env, qxb, qz[b]);

    for (size_t i = 0; i < qx.n; ++i) {
      const double u  = ratio<double>(i * hop, hf);
      const size_t i0 = u;
      const size_t i1 = std::min(i0 + 1, nx - 1);
      const double t  = u - i0;
      for (size_t k = 0; k < qz[b].n; ++k) {
        double mono = 0.0;
        for (size_t c = 0; c < nc; ++c) {
          const double y0 = pwr[c][i0 * nz + k0 + k];
          const double y1 = pwr[c][i1 * nz + k0 + k];
          const double y  = (1.0 - t) * y0 + t * y1;
          ch[b][c + 1].push_back(db<Y>::ratio(std::sqrt(y)));
          mono += y / nc;
        }
        ch[b].front().push_back(db<Y>::ratio(std::sqrt(mono)));
      }
    }
    const Qtz<Y> qy = ydb(ch[b], env.qy);
    band.push_back(as_spectrum(ch[b], qx, qy, qz[b]));
  }
  return band;
}

template <typename X, typename Y, typename Z>
double spectral_peak(const zvec<X, Y, Z>& band) {
  double arg = 0.0;
  double y1  = -k::INF<>;
  for (const Spectrum<X, Y, Z>& each : band) {
    const ms<X>  x  = each.x0();
    const double dy = 10.0 * std::log10(each.qz.rho);
    for (hz<Z> z = each.z0(); z < each.z1(); ++z) {
      const double y = each[x, z] + dy;
      if (y <= y1) continue;
      y1  = y;
      arg = z;
    }
  }
  return arg;
}

template <typename X, typename Y, typename Z>
double spectral_centroid(const zvec<X, Y, Z>& band) {
  double weight = 0.0;
  double moment = 0.0;
  for (const Spectrum<X, Y, Z>& each : band) {
    const ms<X> x = each.x0();
    for (hz<Z> z = each.z0(); z < each.z1(); ++z) {
      const double pwr = db<Y>::lin(each[x, z]);
      weight += pwr;
      moment += z * pwr;
    }
  }
  return moment / weight;
}

template <typename X, typename Y, typename Z>
double spectral_spread(const zvec<X, Y, Z>& band, double centroid) {
  double weight     = 0.0;
  double dispersion = 0.0;
  for (const Spectrum<X, Y, Z>& each : band) {
    const ms<X> x = each.x0();
    for (hz<Z> z = each.z0(); z < each.z1(); ++z) {
      const double pwr = db<Y>::lin(each[x, z]);
      const double dz  = z - centroid;
      weight += pwr;
      dispersion += dz * dz * pwr;
    }
  }
  return std::sqrt(dispersion / weight);
}

template <typename X, typename Y, typename Z>
double spectral_tilt(const zvec<X, Y, Z>& band) {
  vec octave;
  vec lvl;
  for (const Spectrum<X, Y, Z>& each : band) {
    const ms<X>  x  = each.x0();
    const double y0 = each.qy.pos;
    const double dy = 10.0 * std::log10(each.qz.rho);
    for (hz<Z> z = each.z0(); z < each.z1(); ++z) {
      if (each[x, z] <= y0) continue;
      octave.push_back(std::log2(z));
      lvl.push_back(each[x, z] + dy);
    }
  }
  return slope(octave, lvl);
}

template <typename X, typename Y, typename Z>
double spectral_flatness(const zvec<X, Y, Z>& band) {
  double dz_tot  = 0.0;
  double log_tot = 0.0;
  double tot     = 0.0;
  for (const Spectrum<X, Y, Z>& each : band) {
    const ms<X>  x  = each.x0();
    const double y0 = each.qy.pos;
    const double dz = 1.0 / each.qz.rho;
    for (hz<Z> z = each.z0(); z < each.z1(); ++z) {
      if (each[x, z] <= y0) continue;
      const double power = db<Y>::lin(each[x, z]);
      dz_tot += dz;
      log_tot += dz * std::log(power / dz);
      tot += power;
    }
  }
  const double geometric  = std::exp(log_tot / dz_tot);
  const double arithmetic = tot / dz_tot;
  return geometric / arithmetic;
}

template <typename X, typename Y, typename Z>
std::optional<double> spectral_level(const zvec<X, Y, Z>& band) {
  vec lvl;
  vec octave;
  for (const Spectrum<X, Y, Z>& each : band) {
    const ms<X>  x  = each.x0();
    const double y0 = each.qy.pos;
    const double dz = 1.0 / each.qz.rho;
    const double dy = 10.0 * std::log10(each.qz.rho);
    for (hz<Z> z = each.z0(); z < each.z1(); ++z) {
      if (each[x, z] <= y0) continue;
      lvl.push_back(each[x, z] + dy);
      octave.push_back(std::log2((z + dz / 2.0) / (z - dz / 2.0)));
    }
  }
  return otsus_method(lvl, octave, k::Lin);
}

template <typename X, typename Y, typename Z>
std::pair<Z, Z> spectral_edge(const zvec<X, Y, Z>& band, double y) {
  Z z0 = k::INF<Z>;
  Z z1 = 0;
  for (const Spectrum<X, Y, Z>& each : band) {
    const ms<X>  x  = each.x0();
    const double dy = 10.0 * std::log10(each.qz.rho);
    for (hz<Z> z = each.z0(); z < each.z1(); ++z) {
      if (each[x, z] + dy < y) continue;
      z0 = std::min<Z>(z0, z);
      z1 = std::max<Z>(z1, z);
    }
  }
  return {z0, z1};
}

template <typename Z>
Qtz<double> octave_axis(const vec<Qtz<Z>>& qz) {
  double zmin = k::INF<>;
  double zmax = -k::INF<>;
  double rho  = 0.0;
  for (const Qtz<Z>& each : qz) {
    const double dz = 1.0 / each.rho;
    const double z0 = std::log2(each.pos - dz / 2.0);
    const double z1 = std::log2(each.pos + (each.n - 0.5) * dz);

    zmin = std::min(zmin, z0);
    zmax = std::max(zmax, z1);
    rho  = std::max(rho, each.n / (z1 - z0));
  }
  const size_t n = std::ceil((zmax - zmin) * rho);
  return {n, rho, zmin};
}

inline void
octave_spread(double z0, double z1, double level, const Qtz<double>& octave, vec<>& density) {
  const double dz    = 1.0 / octave.rho;
  const size_t i0    = (z0 - octave.pos) * octave.rho;
  const size_t ceil1 = std::ceil((z1 - octave.pos) * octave.rho);
  const size_t i1    = std::min(octave.n, ceil1);
  for (size_t i = i0; i < i1; ++i) {
    const double z       = octave.pos + i * dz;
    const double overlap = std::min(z1, z + dz) - std::max(z0, z);
    if (overlap > 0.0) density[i] += overlap * level * octave.rho;
  }
}

template <typename Z>
vec<> octave_density(const vec<Qtz<Z>>& qz, const Qtz<double>& octave, const vec<>& spectrum) {
  vec    density(octave.n);
  size_t row = 0;
  for (const Qtz<Z>& each : qz) {
    const double dz = 1.0 / each.rho;
    for (size_t k = 0; k < each.n; ++k) {
      const double z   = each.pos + k * dz;
      const double z0  = std::log2(z - dz / 2.0);
      const double z1  = std::log2(z + dz / 2.0);
      const double lvl = spectrum[row] / (z1 - z0);
      octave_spread(z0, z1, lvl, octave, density);
      ++row;
    }
  }
  return density;
}

inline Qtz<double> octave_union(const vec<Qtz<double>>& qz) {
  double zmin = k::INF<>;
  double zmax = -k::INF<>;
  double rho  = 0.0;
  for (const Qtz<double>& each : qz) {
    zmin = std::min(zmin, each.pos);
    zmax = std::max(zmax, each.pos + each.n / each.rho);
    rho  = std::max(rho, each.rho);
  }
  const size_t n = std::ceil((zmax - zmin) * rho);
  return {n, rho, zmin};
}

inline vvec<>
octave_resample(const vec<Qtz<double>>& qz, const vvec<>& density, const Qtz<double>& octave) {
  vvec resample(density.size(), vec(octave.n));
  for (size_t s = 0; s < density.size(); ++s) {
    const double dz = 1.0 / qz[s].rho;
    for (size_t i = 0; i < density[s].size(); ++i) {
      const double z0 = qz[s].pos + i * dz;
      octave_spread(z0, z0 + dz, density[s][i], octave, resample[s]);
    }
  }
  return resample;
}

template <typename X, typename Y, typename Z>
double spectral_rolloff(const Spectrum<X, Y, Z>& env, double fraction) {
  const Qtz<Z>& qz  = env.qz;
  const double  y0  = env.qy.pos;
  const Y*      mag = env.ch.front();
  vec           cumulative(qz.n + 1);
  cumulative[0] = 0.0;
  for (size_t k = 0; k < qz.n; ++k) cumulative[k + 1] = cumulative[k] + (mag[k] - y0);
  const double target      = fraction * cumulative.back();
  size_t       rolloff_bin = qz.n - 1;
  for (size_t k = 0; k < qz.n; ++k)
    if (cumulative[k + 1] >= target) {
      rolloff_bin = k;
      break;
    }
  return env.z_grid(env.first() + rolloff_bin).value;
}
}
