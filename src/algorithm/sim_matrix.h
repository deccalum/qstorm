#pragma once

#include <algorithm>
#include <bit>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <random>
#include <span>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "algorithm/dsp_spectrum.h"
#include "algorithm/exec_threads.h"
#include "algorithm/math_constant.h"
#include "algorithm/math_scalar.h"

namespace qstorm::matrix {
using namespace algorithm;

inline constexpr size_t dot_lane  = 16;
inline constexpr float  negligible = std::numeric_limits<float>::min();

template <typename R, typename C>
struct Matrix {
  Qtz<R>     row;
  Qtz<C>     col;
  vec<float> cell;
};

// **ARD** - Automatic Relevance Determination
//
// Defalt values for the non-negative matrix factorization.
//
// - `seed`: sets the random starting entries of W and H. IS-NMF has
//           local minima (each start can converge to a different stationary
//           point) meaning a different seed would give different results.
// - `l1`: chooses which prior W and H get as well as sets shape anf formula
//         for `b`. defaults to true because ℓ1-ARD recovers the true number
//         of components more consistently across runs than ℓ2-ARD (§6.1).
//    - true (ℓ1): exponential prior resulting in sparser spectrums and
//                 activations.
//                 Component mass is Σw + Σh, update penalty is 1/λ_k.
//    - false (ℓ2): half-normal prior resulting in smoother components by
//                  shrinking everything proportionally.
//                  Component mass is ½Σw² + ½Σh², update penalty is w/λ_k.
// - `K`: initial number of components for the factorization, which sets
//        the width of W and the height of H. defaults to `64` because ARD
//        only prunes, so K must exceed the expected number of components.
//        The paper uses 2 to 3 times the true count (§6).
// - `tol`: loop limit in relative change of λ. Iteration stops when
//          max_k |λ_k − λ_k,prev| / λ_k,prev < tol (τ in §4.4).
// - `a`: inverse-gamma prior shape affecting λ_k.
//        A small `a` gives a weak, heavy-tailed prior, where the data mostly
//        decides each λ_k. A large `a` pull every λ_k toward the centre.
//        defaults to 5 because the paper recommends a small value w.r.t.
//        F + N (§4.5.2) and uses a = 5 on a piano spectrogram (§6.3).
//
// Prior: λ_k ~ IG(a, b), shared by column k of W and row k of H (eq. 14).
//        b is set from the sample mean of V (eq. 38).
// IS: Itakura-Saito divergence, the β = 0 case of the β-divergence, used
//     on power spectrograms (Févotte, Bertin & Durrieu 2009).
struct Ard {
  uint32_t seed /* mt19937 seed for the starting W and H       */ = 1;
  bool     l1 /*   true: ℓ1 exponential, false: ℓ2 half-normal */ = true;
  int      K /*    starting components, W width and H height   */ = 64;
  float    tol /*  stop when max relative λ change is below    */ = 1e-5;
  int      a /*    inverse-gamma shape of the λ_k prior        */ = 5;
};

struct Relevance {
  vec<float> lambda;      // λ_k per component
  double     floor;       // λ of a pruned component, b / shape
  double     divergence;  // IS divergence at convergence
  int        iteration;   // update passes until λ settled
  vec<>      trace;       // IS divergence before each pass, then at convergence
  vec<>      residual;    // IS divergence of each V row at convergence
  vec<>      removed;     // IS divergence gained when w_k h_k leaves WH
  vec<>      share;       // sum(W_k)·sum(H_k) / sum(V) per component
  vec<bool>  relevant;    // 2·removed > (F + N)·ln(measured cells), BIC
};

template <typename X, typename Y, typename Z>
Matrix<int, X> power_matrix(const zvec<X, Y, Z>& band) {
  const size_t nx = band.front().qx.n;
  size_t       nz = 0;
  for (const Spectrum<X, Y, Z>& each : band) nz += each.qz.n;

  const Qtz<int> qz{nz, 1, 0};
  const Qtz<X>   qx = band.front().qx;
  Matrix<int, X> grid;
  grid.row = qz;
  grid.col = qx;
  grid.cell.resize(nz * nx);

  size_t i0 = 0;
  for (const Spectrum<X, Y, Z>& each : band) {
    for (ms<X> x = each.x0(); x < each.x1(); x++)
      for (hz<Z> z = each.z0(); z < each.z1(); z++) {
        const size_t idx = i0 * nx + cell<Z, X>(z, x, each.qz, each.qx);
        const db<Y>  y   = each[x, z];
        grid.cell[idx]   = db<Y>::lin(y);
      }
    i0 += each.qz.n;
  }
  return grid;
}

template <typename X, typename Y, typename Z>
vec<> band_row_weight(const zvec<X, Y, Z>& band) {
  const double rho1 = band.back().qz.rho;
  vec          weight;
  for (const Spectrum<X, Y, Z>& each : band)
    for (size_t f = 0; f < each.qz.n; ++f) weight.push_back(rho1 / each.qz.rho);
  return weight;
}

template <typename X, typename Y, typename Z>
vec<> band_residual(const zvec<X, Y, Z>& band, const vec<>& residual) {
  const size_t nx = band.front().qx.n;
  vec          mean;
  size_t       i0 = 0;
  for (const Spectrum<X, Y, Z>& each : band) {
    const std::span<const double> rows = std::span<const double>(residual).subspan(i0, each.qz.n);
    const double                  tot  = sum(rows, std::identity{});
    const size_t                  n    = each.qz.n * nx;
    mean.push_back(tot / n);
    i0 += each.qz.n;
  }
  return mean;
}

inline vec<> reconstruct_row(const vvec<>& spectrum, const vvec<>& activation, size_t f) {
  vec restruct(activation.front().size());
  for (size_t k = 0; k < spectrum.size(); ++k)
    for (size_t n = 0; n < restruct.size(); ++n) restruct[n] += spectrum[k][f] * activation[k][n];
  return restruct;
}

template <typename X, typename Y, typename Z>
zvec<X, Y, Z> reconstruct_layer(const zvec<X, Y, Z>& band, const vvec<>& spectrum,
                                const vvec<>& activation, vvvec<Y>& channel) {
  const size_t nx = activation.front().size();
  channel.assign(band.size(), vvec<Y>(1));
  zvec<X, Y, Z> layer;
  size_t        i0 = 0;
  for (size_t b = 0; b < band.size(); ++b) {
    const size_t nz = band[b].qz.n;
    channel[b][0].resize(nx * nz);
    for (size_t j = 0; j < nz; ++j) {
      const vec restruct = reconstruct_row(spectrum, activation, i0 + j);
      for (size_t n = 0; n < nx; ++n)
        channel[b][0][n * nz + j] = db<Y>::ratio(std::sqrt(restruct[n]));
    }

    const Qtz<Y> qy = ydb(channel[b], band[b].qy);
    layer.push_back(as_spectrum(channel[b], band[b].qx, qy, band[b].qz));
    i0 += nz;
  }

  return layer;
}

template <typename X, typename Y, typename Z>
zvec<X, Y, Z> residual_layer(const zvec<X, Y, Z>& band, const vvec<>& spectrum,
                             const vvec<>& activation, vvvec<Y>& channel) {
  const Matrix<int, X> V  = power_matrix(band);
  const size_t         nx = V.col.n;
  channel.assign(band.size(), vvec<Y>(1));
  zvec<X, Y, Z> layer;
  size_t        i0 = 0;
  for (size_t b = 0; b < band.size(); ++b) {
    const size_t nz = band[b].qz.n;
    channel[b][0].resize(nx * nz);
    for (size_t j = 0; j < nz; ++j) {
      const size_t f        = i0 + j;
      const vec    restruct = reconstruct_row(spectrum, activation, f);
      for (size_t n = 0; n < nx; ++n) {
        const double left         = std::max(0.0, V.cell[f * nx + n] - restruct[n]);
        channel[b][0][n * nz + j] = db<Y>::ratio(std::sqrt(left));
      }
    }
    const Qtz<Y> qy = ydb(channel[b], band[b].qy);
    layer.push_back(as_spectrum(channel[b], band[b].qx, qy, band[b].qz));
    i0 += nz;
  }
  return layer;
}

inline vvec<> mean_centered(const vvec<>& row) {
  vvec centered;
  for (const vec<>& each : row) {
    const double mean = average(each, std::identity{});
    vec          deviation;
    for (double value : each) deviation.push_back(value - mean);
    centered.push_back(deviation);
  }
  return centered;
}

inline vvec<> unit_length(const vvec<>& row) {
  vvec unit;
  for (const vec<>& each : row) {
    const double dot = std::inner_product(each.begin(), each.end(), each.begin(), 0.0);
    const double len = std::sqrt(dot);
    vec          scaled;
    for (double value : each) scaled.push_back(value / len);
    unit.push_back(scaled);
  }
  return unit;
}

inline vec<size_t> complete_linkage(const vvec<>& similarity, const vvec<bool>& related) {
  const size_t n = similarity.size();
  vvec<size_t> member;
  for (size_t i = 0; i < n; ++i) member.push_back({i});

  bool merged = true;
  while (merged) {
    merged      = false;
    double best = -k::INF<>;
    size_t into = 0;
    size_t from = 0;
    for (size_t a = 0; a < member.size(); ++a)
      for (size_t b = a + 1; b < member.size(); ++b) {
        bool   every = true;
        double least = k::INF<>;
        for (size_t i : member[a])
          for (size_t j : member[b]) {
            every = every && related[i][j];
            least = std::min(least, similarity[i][j]);
          }
        if (every && least > best) {
          best   = least;
          into   = a;
          from   = b;
          merged = true;
        }
      }
    if (merged) {
      member[into].insert(member[into].end(), member[from].begin(), member[from].end());
      member.erase(member.begin() + from);
    }
  }

  vec<size_t> set(n);
  for (const vec<size_t>& each : member) {
    const size_t first = *std::ranges::min_element(each);
    for (size_t i : each) set[i] = first;
  }

  return set;
}

inline vec<std::complex<double>> pad_transform(const vec<>& row) {
  vec<std::complex<double>> x(std::bit_ceil(2 * row.size()));
  for (size_t n = 0; n < row.size(); ++n) x[n] = row[n];
  fft(std::span<std::complex<double>>(x));
  return x;
}

inline vec<> linear_corr(const vec<std::complex<double>>& a, const vec<std::complex<double>>& b) {
  const size_t              n = a.size();
  vec<std::complex<double>> sab(n);
  for (size_t i = 0; i < n; ++i) sab[i] = std::conj(std::conj(a[i]) * b[i]);
  fft(std::span<std::complex<double>>(sab));
  vec lag(n);
  for (size_t i = 0; i < n; ++i) lag[i] = sab[i].real() / n;
  return lag;
}

inline vec<>
circular_corr(const vec<std::complex<double>>& a, const vec<std::complex<double>>& b, size_t n) {
  const vec    lag = linear_corr(a, b);
  const size_t pad = lag.size();
  vec          circular(n);
  for (size_t i = 0; i < n; ++i) circular[i] = lag[i] + lag[pad + i - n];
  return circular;
}

inline size_t autocorr_width(const vec<>& autocorr) {
  size_t lag = 1;
  while (lag < autocorr.size() && autocorr[lag] > 0.0) ++lag;
  return lag;
}

inline vvec<bool> shift_related(const vvec<>& row) {
  const vvec                 unit = unit_length(mean_centered(row));
  const size_t               nk   = unit.size();
  const size_t               nx   = unit.front().size();
  vvec<std::complex<double>> transform;
  vec<size_t>                width;
  for (const vec<>& each : unit) {
    transform.push_back(pad_transform(each));
    width.push_back(autocorr_width(circular_corr(transform.back(), transform.back(), nx)));
  }
  vvec<bool> related(nk, vec<bool>(nk));
  for (size_t i = 0; i < nk; ++i)
    for (size_t j = i + 1; j < nk; ++j) {
      const vec    lag    = circular_corr(transform[i], transform[j], nx);
      const size_t near   = std::max(width[i], width[j]);
      bool         tested = false;
      bool         beaten = true;
      for (size_t shift = near; shift + near <= nx; ++shift) {
        tested = true;
        beaten = beaten && lag[0] > lag[shift];
      }
      related[i][j] = tested && beaten;
      related[j][i] = related[i][j];
    }
  return related;
}

inline vvec<bool> overlap_related(const vvec<>& row) {
  const vvec                 unit = unit_length(mean_centered(row));
  const size_t               nk   = row.size();
  const size_t               nx   = row.front().size();
  vvec<std::complex<double>> transform;
  vvec                       tot;
  vvec                       square;
  vec<size_t>                width;
  for (size_t i = 0; i < nk; ++i) {
    const vec<std::complex<double>> centered = pad_transform(unit[i]);
    width.push_back(autocorr_width(circular_corr(centered, centered, nx)));
    transform.push_back(pad_transform(row[i]));
    vec running{0.0};
    vec running_square{0.0};

    for (double value : row[i]) {
      running.push_back(running.back() + value);
      running_square.push_back(running_square.back() + value * value);
    }

    tot.push_back(running);
    square.push_back(running_square);
  }

  vvec<bool> related(nk, vec<bool>(nk));
  for (size_t i = 0; i < nk; ++i)
    for (size_t j = i + 1; j < nk; ++j) {
      const vec    lag  = linear_corr(transform[i], transform[j]);
      const size_t pad  = lag.size();
      const size_t near = std::max(width[i], width[j]);
      vec          corr;
      for (size_t gap = 0; gap + near <= nx; gap = gap == 0 ? near : gap + 1)
        for (const bool forward : {true, false}) {
          const size_t over    = nx - gap;
          const size_t a       = forward ? 0 : gap;
          const size_t b       = forward ? gap : 0;
          const double sab     = lag[forward ? gap : (pad - gap) % pad];
          const double sa      = tot[i][a + over] - tot[i][a];
          const double sb      = tot[j][b + over] - tot[j][b];
          const double saa     = square[i][a + over] - square[i][a];
          const double sbb     = square[j][b + over] - square[j][b];
          const double cov     = sab - sa * sb / over;
          const double var     = (saa - sa * sa / over) * (sbb - sb * sb / over);
          corr.push_back(cov / std::sqrt(var));
        }
      const double aligned = corr.front();
      bool         beaten  = corr.size() > 2;
      for (size_t m = 2; m < corr.size(); ++m) beaten = beaten && aligned > corr[m];
      related[i][j] = beaten;
      related[j][i] = related[i][j];
    }
  return related;
}

inline vvec<> cos_matrix(const vvec<>& row) {
  const size_t n = row.size();
  vvec         similarity(n, vec(n));
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j) similarity[i][j] = cosine(row[i], row[j]);
  return similarity;
}
inline vvec<> corr_matrix(const vvec<>& row) { return cos_matrix(mean_centered(row)); }

inline vvec<> reconstruct_power(const vvec<>& spectrum, const vvec<>& activation) {
  vvec restruct;
  for (size_t f = 0; f < spectrum.front().size(); ++f)
    restruct.push_back(reconstruct_row(spectrum, activation, f));
  return restruct;
}

inline vvec<>
component_mask(const vvec<>& spectrum, const vvec<>& activation, const vvec<>& restruct,
               const Matrix<int, double>& power, const vec<size_t>& member) {
  const size_t nx = activation.front().size();
  vvec         mask(restruct.size(), vec(nx));
  for (size_t f = 0; f < restruct.size(); ++f)
    for (size_t n = 0; n < nx; ++n) {
      double share = 0.0;
      for (size_t k : member) share += spectrum[k][f] * activation[k][n];
      const double whole = std::max<double>(power.cell[f * nx + n], restruct[f][n]);
      mask[f][n]         = whole > 0.0 ? share / whole : 0.0;
    }
  return mask;
}

inline vvec<> residual_mask(const vvec<>& restruct, const Matrix<int, double>& power) {
  const size_t nx = power.col.n;
  vvec         mask(restruct.size(), vec(nx));
  for (size_t f = 0; f < restruct.size(); ++f)
    for (size_t n = 0; n < nx; ++n) {
      const double v     = power.cell[f * nx + n];
      const double whole = std::max(v, restruct[f][n]);
      mask[f][n]         = whole > 0.0 ? std::max(0.0, v - restruct[f][n]) / whole : 0.0;
    }
  return mask;
}

template <typename S>
vvec<S> component_audio(const Spectrum<double>& env, const zvec<double>& stft,
                        const vvec<>& spectrum, const vvec<>& activation,
                        const Matrix<int, double>& power, const vvec<size_t>& member) {
  const vec<Qtz<double>> qz      = band_axis(stft);
  const Qtz<double>      qx      = stft.front().qx;
  const auto             framed  = band_spectrum(env, qz);
  const vvec             restruc = reconstruct_power(spectrum, activation);
  vvec<S>                audio;
  for (const vec<size_t>& each : member) {
    const vvec mask = component_mask(spectrum, activation, restruc, power, each);
    const vvec part = inverse_stft_masked(env, qx, qz, framed, mask);
    audio.push_back(xc_qtzd<S>(part));
  }
  return audio;
}

template <typename S>
vec<S> residual_audio(const Spectrum<double>& env, const zvec<double>& stft, const vvec<>& spectrum,
                      const vvec<>& activation, const Matrix<int, double>& power) {
  const vec<Qtz<double>> qz      = band_axis(stft);
  const Qtz<double>      qx      = stft.front().qx;
  const auto             framed  = band_spectrum(env, qz);
  const vvec             restruc = reconstruct_power(spectrum, activation);
  const vvec             mask    = residual_mask(restruc, power);
  const vvec             part    = inverse_stft_masked(env, qx, qz, framed, mask);
  return xc_qtzd<S>(part);
}

inline vec<int> min_cost_match(const vvec<>& a) {
  const int n = a.size();
  vec       u(n + 1, 0.0), v(n + 1, 0.0);
  vec<int>  p(n + 1, 0), way(n + 1, 0);
  for (int i = 1; i <= n; ++i) {
    p[0]         = i;
    int       j0 = 0;
    vec       minv(n + 1, k::INF<>);
    vec<bool> used(n + 1, false);
    do {
      used[j0]        = true;
      const int i0    = p[j0];
      double    delta = k::INF<>;
      int       j1    = -1;
      for (int j = 1; j <= n; ++j)
        if (!used[j]) {
          const double cur = a[i0 - 1][j - 1] - u[i0] - v[j];
          if (cur < minv[j]) {
            minv[j] = cur;
            way[j]  = j0;
          }
          if (minv[j] < delta) {
            delta = minv[j];
            j1    = j;
          }
        }
      for (int j = 0; j <= n; ++j) {
        if (used[j]) {
          u[p[j]] += delta;
          v[j] -= delta;
        } else
          minv[j] -= delta;
      }
      j0 = j1;
    } while (p[j0] != 0);
    do {
      const int j1 = way[j0];
      p[j0]        = p[j1];
      j0           = j1;
    } while (j0);
  }
  vec<int> match(n, -1);
  for (int j = 1; j <= n; ++j) match[p[j] - 1] = j - 1;
  return match;
}

template <typename R, typename C>
std::pair<Matrix<R, int>, Matrix<int, C>>
init_factorization(const Matrix<R, C>& V, const Ard& ard) {
  const size_t nz = V.row.n;
  const size_t nx = V.col.n;
  const size_t nk = ard.K;

  const Qtz<R>   qz = V.row;
  const Qtz<C>   qx = V.col;
  const Qtz<int> qk{nk, 1, 0};

  Matrix<R, int> W;
  W.row = qz;
  W.col = qk;
  W.cell.resize(nz * nk);

  Matrix<int, C> H;
  H.row = qk;
  H.col = qx;
  H.cell.resize(nk * nx);

  std::mt19937                          rng(ard.seed);
  std::uniform_real_distribution<float> uniform(0.1f, 1.0f);
  for (float& e : W.cell) e = uniform(rng);
  for (float& e : H.cell) e = uniform(rng);
  return {W, H};
}

template <typename R, typename P, typename C>
void matmul_range(size_t i0, size_t i1, const Matrix<R, P>& A, const Matrix<P, C>& B,
                  Matrix<R, C>& result) {
  const size_t pair = A.col.n;
  const size_t col  = B.col.n;
  for (size_t i = i0; i < i1; ++i) {
    float*       out = result.cell.data() + i * col;
    const float* a   = A.cell.data() + i * pair;
    std::fill(out, out + col, 0.0f);
    size_t k = 0;
    for (; k + 4 <= pair; k += 4) {
      const float  a0 = a[k];
      const float  a1 = a[k + 1];
      const float  a2 = a[k + 2];
      const float  a3 = a[k + 3];
      const float* f0 = B.cell.data() + k * col;
      const float* f1 = f0 + col;
      const float* f2 = f1 + col;
      const float* f3 = f2 + col;
      for (size_t j = 0; j < col; ++j)
        out[j] = out[j] + a0 * f0[j] + a1 * f1[j] + a2 * f2[j] + a3 * f3[j];
    }
    for (; k < pair; ++k) {
      const float  a0   = a[k];
      const float* from = B.cell.data() + k * col;
      for (size_t j = 0; j < col; ++j) out[j] += a0 * from[j];
    }
  }
}

template <typename R, typename P, typename C>
void matmul(const Matrix<R, P>& A, const Matrix<P, C>& B, Matrix<R, C>& result) {
  const size_t row  = A.row.n;
  const size_t pair = A.col.n;
  const size_t col  = B.col.n;
  in_parallel(row, row * pair * col, matmul_range<R, P, C>, A, B, result);
}

template <typename P, typename R, typename C>
void matmul_row_range(size_t i0, size_t i1, const Matrix<P, R>& A, const Matrix<P, C>& B,
                      Matrix<R, C>& result) {
  const size_t pair = A.row.n;
  const size_t row  = A.col.n;
  const size_t col  = B.col.n;
  std::fill(result.cell.data() + i0 * col, result.cell.data() + i1 * col, 0.0f);
  size_t k = 0;
  for (; k + 4 <= pair; k += 4) {
    const float* f0 = B.cell.data() + k * col;
    const float* f1 = f0 + col;
    const float* f2 = f1 + col;
    const float* f3 = f2 + col;
    for (size_t i = i0; i < i1; ++i) {
      const float a0  = A.cell[k * row + i];
      const float a1  = A.cell[(k + 1) * row + i];
      const float a2  = A.cell[(k + 2) * row + i];
      const float a3  = A.cell[(k + 3) * row + i];
      float*      out = result.cell.data() + i * col;
      for (size_t j = 0; j < col; ++j)
        out[j] = out[j] + a0 * f0[j] + a1 * f1[j] + a2 * f2[j] + a3 * f3[j];
    }
  }
  for (; k < pair; ++k) {
    const float* from = B.cell.data() + k * col;
    for (size_t i = i0; i < i1; ++i) {
      const float a0  = A.cell[k * row + i];
      float*      out = result.cell.data() + i * col;
      for (size_t j = 0; j < col; ++j) out[j] += a0 * from[j];
    }
  }
}

template <typename R, typename P, typename C>
void matmul_col_range(size_t i0, size_t i1, const Matrix<R, P>& A, const Matrix<C, P>& B,
                      Matrix<R, C>& result) {
  const size_t col  = B.row.n;
  const size_t pair = A.col.n;
  for (size_t i = i0; i < i1; ++i) {
    const float* a = A.cell.data() + i * pair;
    for (size_t j = 0; j < col; ++j) {
      const float* b = B.cell.data() + j * pair;
      float        lane[dot_lane]{};
      size_t       k = 0;
      for (; k + dot_lane <= pair; k += dot_lane)
        for (size_t m = 0; m < dot_lane; ++m) lane[m] += a[k + m] * b[k + m];
      float tot = 0.0f;
      for (float each : lane) tot += each;
      for (; k < pair; ++k) tot += a[k] * b[k];
      result.cell[i * col + j] = tot;
    }
  }
}

template <typename P, typename R, typename C>
void matmul_row(const Matrix<P, R>& A, const Matrix<P, C>& B, Matrix<R, C>& result) {
  const size_t pair = A.row.n;
  const size_t row  = A.col.n;
  const size_t col  = B.col.n;
  in_parallel(row, pair * row * col, matmul_row_range<P, R, C>, A, B, result);
}

template <typename R, typename P, typename C>
void matmul_col(const Matrix<R, P>& A, const Matrix<C, P>& B, Matrix<R, C>& result) {
  const size_t row  = A.row.n;
  const size_t col  = B.row.n;
  const size_t pair = A.col.n;
  in_parallel(row, row * col * pair, matmul_col_range<R, P, C>, A, B, result);
}

template <typename R, typename K, typename C>
void relevance_weight_range(size_t i0, size_t i1, const Matrix<R, K>& W, const Matrix<K, C>& H,
                            const Ard& ard, double b, double shape, vec<float>& lambda) {
  const size_t nz = W.row.n;
  const size_t nk = W.col.n;
  const size_t nx = H.col.n;
  for (size_t k = i0; k < i1; k++) {
    double mass = 0.0;
    for (size_t f = 0; f < nz; ++f) {
      const double w = W.cell[f * nk + k];
      mass += ard.l1 ? w : 0.5 * w * w;
    }
    for (size_t n = 0; n < nx; ++n) {
      const double h = H.cell[k * nx + n];
      mass += ard.l1 ? h : 0.5 * h * h;
    }
    lambda[k] = (mass + b) / shape;
  }
}

template <typename R, typename K, typename C>
void relevance_weight(const Matrix<R, K>& W, const Matrix<K, C>& H, const Ard& ard, double b,
                      double shape, vec<float>& lambda) {
  const size_t nz   = W.row.n;
  const size_t nk   = W.col.n;
  const size_t nx   = H.col.n;
  const size_t work = nk * (nz + nx);
  in_parallel(nk, work, relevance_weight_range<R, K, C>, W, H, ard, b, shape, lambda);
}

template <typename R, typename C>
void cell_ratio(size_t i0, size_t i1, const Matrix<R, C>& V, const Matrix<R, C>& WH,
                Matrix<R, C>& P, Matrix<R, C>& Q) {
  for (size_t i = i0; i < i1; ++i) {
    const float ratio = V.cell[i] / WH.cell[i];
    P.cell[i]         = ratio / WH.cell[i];
    Q.cell[i]         = 1.0f / WH.cell[i];
  }
}

template <typename K, typename C>
void cell_activation(size_t i0, size_t i1, const Matrix<K, C>& WtP, const Matrix<K, C>& WtQ,
                     const Ard& ard, const vec<float>& lambda, Matrix<K, C>& H) {
  const size_t nx = H.col.n;
  for (size_t i = i0; i < i1; ++i) {
    const size_t k       = i / nx;
    const float  penalty = ard.l1 ? 1.0f / lambda[k] : H.cell[i] / lambda[k];
    const float  updated = H.cell[i] * std::sqrt(WtP.cell[i] / (WtQ.cell[i] + penalty));
    H.cell[i]            = updated < negligible ? 0.0f : updated;
  }
}

template <typename R, typename C>
void row_divergence(size_t i0, size_t i1, const Matrix<R, C>& V, const Matrix<R, C>& WH,
                    vec<>& residual) {
  const size_t nx = V.col.n;
  for (size_t f = i0; f < i1; ++f) {
    double row_tot = 0.0;
    for (size_t n = 0; n < nx; ++n) {
      const double observed = V.cell[f * nx + n];
      const double ratio    = observed / WH.cell[f * nx + n];
      row_tot += ratio - std::log(ratio) - 1.0;
    }
    residual[f] = row_tot;
  }
}

template <typename R, typename K, typename C>
void row_removal(size_t i0, size_t i1, const Matrix<R, C>& V, const Matrix<R, K>& W,
                 const Matrix<K, C>& H, vec<>& residual, vec<>& removed_row) {
  const size_t nx = V.col.n;
  const size_t nk = W.col.n;
  vec          part(nk);
  for (size_t f = i0; f < i1; ++f) {
    double row_tot = 0.0;
    for (size_t n = 0; n < nx; ++n) {
      double restruct = 0.0;
      for (size_t k = 0; k < nk; ++k) {
        const double w = W.cell[f * nk + k];
        part[k]        = w * H.cell[k * nx + n];
        restruct += part[k];
      }
      const double value = V.cell[f * nx + n];
      const double kept  = value / restruct - std::log(value / restruct) - 1.0;
      row_tot += kept;
      for (size_t k = 0; k < nk; ++k) {
        const double rest    = restruct - part[k];
        const double without = rest > 0.0 ? value / rest - std::log(value / rest) - 1.0 : k::INF<>;
        removed_row[f * nk + k] += without - kept;
      }
    }
    residual[f] = row_tot;
  }
}

template <typename R, typename K>
void cell_spectrum(size_t i0, size_t i1, const Matrix<R, K>& PHt, const Matrix<R, K>& QHt,
                   const Ard& ard, const vec<float>& lambda, Matrix<R, K>& W) {
  const size_t nk = W.col.n;
  for (size_t i = i0; i < i1; ++i) {
    const size_t k       = i % nk;
    const float  penalty = ard.l1 ? 1.0f / lambda[k] : W.cell[i] / lambda[k];
    const float  updated = W.cell[i] * std::sqrt(PHt.cell[i] / (QHt.cell[i] + penalty));
    W.cell[i]            = updated < negligible ? 0.0f : updated;
  }
}

// non-negative matrix factorization
//
// V ≈ W·H with every entry ≥ 0, which forces a parts-based representation.
//
// - V: power spectrogram (nz × nx)
// - W: dictionary, one spectrum per component (nz × K)
// - H: activations, one timeline per component (K × nx)
//
// Minimizes the IS divergence D(V | WH) plus the ARD penalty on W and H.
// W and H start from random positive values (`init_factorization`).
//
// Each pass updates:
// - H ← H ∘ (Wᵀ(V / (WH)²) / (Wᵀ(1 / WH) + penalty))^½
// - W ← W ∘ ((V / (WH)²)Hᵀ / ((1 / WH)Hᵀ + penalty))^½
// - λ_k ← (mass_k + b) / shape
//
// Ratios of positive numbers keep W and H ≥ 0.
// Stops when the largest relative change in λ falls below `ard.tol`.
template <typename R, typename K, typename C>
Relevance non_neg_matrix_factorization(const Matrix<R, C>& V, Matrix<R, K>& W, Matrix<K, C>& H,
                                       const Ard& ard, const vec<>& row_weight) {
  const size_t nz    = V.row.n;
  const size_t nx    = V.col.n;
  const size_t nk    = W.col.n;
  const size_t plane = nz * nx;

  const double tot   = sum(V.cell, std::identity{});
  const double mean  = tot / plane;
  const double comp  = nk;
  const double edge  = nz + nx;

  const double l1_scale = (ard.a - 1.0) * (ard.a - 2.0);
  const double l1_shape = edge + ard.a + 1.0;
  const double l2_shape = edge / 2.0 + ard.a + 1.0;

  const double b     = ard.l1 ? std::sqrt(l1_scale * mean / comp)
                              : k::hfpi * (ard.a - 1.0) * mean / comp;
  const double shape = ard.l1 ? l1_shape : l2_shape;

  Matrix<R, C> WH;
  WH.row = V.row;
  WH.col = V.col;
  WH.cell.resize(plane);

  Matrix<R, C> P = WH;
  Matrix<R, C> Q = WH;

  Matrix<K, C> WtP;
  WtP.row = W.col;
  WtP.col = V.col;
  WtP.cell.resize(nk * nx);

  Matrix<K, C> WtQ = WtP;

  Matrix<R, K> PHt;
  PHt.row = V.row;
  PHt.col = H.row;
  PHt.cell.resize(nz * nk);

  Matrix<R, K> QHt = PHt;

  const size_t nkx = nk * nx;
  const size_t nzk = nz * nk;

  vec<float> lambda(nk);
  vec<float> prev(nk);
  relevance_weight(W, H, ard, b, shape, lambda);

  vec residual(nz);
  vec removed_row(nzk);
  vec trace;
  int iteration = 0;

  float change = ard.tol;
  while (change >= ard.tol) {
    prev = lambda;

    matmul(W, H, WH);
    in_parallel(nz, plane, row_divergence<R, C>, V, WH, residual);
    trace.push_back(sum(residual, std::identity{}));
    in_parallel(plane, plane, cell_ratio<R, C>, V, WH, P, Q);
    matmul_row(W, P, WtP);
    matmul_row(W, Q, WtQ);
    in_parallel(nkx, nkx, cell_activation<K, C>, WtP, WtQ, ard, lambda, H);

    matmul(W, H, WH);
    in_parallel(plane, plane, cell_ratio<R, C>, V, WH, P, Q);
    matmul_col(P, H, PHt);
    matmul_col(Q, H, QHt);
    in_parallel(nzk, nzk, cell_spectrum<R, K>, PHt, QHt, ard, lambda, W);

    relevance_weight(W, H, ard, b, shape, lambda);
    change = 0.0f;
    for (size_t k = 0; k < nk; ++k) {
      const float step = std::fabs(lambda[k] - prev[k]) / prev[k];
      change           = std::max(change, step);
    }
    iteration++;
  }
  in_parallel(nz, plane * nk, row_removal<R, K, C>, V, W, H, residual, removed_row);
  const double divergence = sum(residual, std::identity{});
  trace.push_back(divergence);

  vec removed(nk);
  vec share(nk);
  for (size_t k = 0; k < nk; ++k) {
    double spectrum = 0.0;
    for (size_t f = 0; f < nz; ++f) {
      spectrum += W.cell[f * nk + k];
      removed[k] += removed_row[f * nk + k] * row_weight[f];
    }
    double activation = 0.0;
    for (size_t n = 0; n < nx; ++n) activation += H.cell[k * nx + n];
    share[k] = spectrum * activation / tot;
  }

  const double floor   = b / shape;
  double       measure = 0.0;
  for (double weight : row_weight) measure += weight * nx;
  const double penalty = edge * std::log(measure);
  vec<bool>    relevant(nk);
  for (size_t k = 0; k < nk; ++k) relevant[k] = 2.0 * removed[k] > penalty;

  return {lambda, floor, divergence, iteration, trace, residual, removed, share, relevant};
}
}
