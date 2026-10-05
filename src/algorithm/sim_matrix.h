#pragma once

#include <cstddef>
#include <random>
#include <vector>
#include "algorithm/axis_quantizer.h"

namespace qstorm::algorithm {

template <typename R, typename C>
struct Matrix {
  Quantization<R>    row;
  Quantization<C>    col;
  std::vector<float> cell;
};

template <typename R, typename C>
Matrix<R, C> matrix(const Spectrogram<float>& spectrogram) {
  const size_t hz_extent = spectrogram.z_hz.extent;
  const size_t     frame = spectrogram.x_ms.extent;
  const size_t channel   = spectrogram.channel.size();
  const size_t plane     = hz_extent*frame;

  const Quantization<R> Row = spectrogram.z_hz.as<R>();
  const Quantization<C> Cell = spectrogram.x_ms.as<C>();
  Matrix<R, C> grid;
  grid.row = Row;
  grid.col = Cell;
  grid.cell.resize(channel*plane);

  for (size_t i = 0; i < channel; i++)
    for (size_t j = 0; j < frame; j++)
      for (size_t k = 0; k < hz_extent; k++)
        grid.cell[i*plane + k*frame +j] = spectrogram.channel[i][j*hz_extent +k];
  return grid;
}

template <typename R, typename C> std::pair<Matrix<R, int>, Matrix<int, C>>
init_factorization(const Spectrogram<float>& spectrogram, int rank, uint32_t seed) {
  const size_t hz_extent = spectrogram.z_hz.extent;
  const size_t     frame = spectrogram.x_ms.extent;
  const size_t   channel = spectrogram.channel.size();
  const size_t component = static_cast<size_t>(rank);

  const Quantization<R>  Hertz = spectrogram.z_hz.as<R>();
  const Quantization<C>   Time = spectrogram.x_ms.as<C>();
  const Quantization<int> Rank{component, 1, 0};

  Matrix<R, int> W;
  W.row = Hertz;
  W.col = Rank;
  W.cell.resize(channel*hz_extent*component);

  Matrix<int, C> H;
  H.row = Rank;
  H.col = Time;
  H.cell.resize(component*frame);

  std::mt19937 rng(seed);
  std::uniform_real_distribution<float> uniform(0.1f, 1.0f);
  for (float& e : W.cell) e = uniform(rng);
  for (float& e : H.cell) e = uniform(rng);
  return {W, H};
}

template <typename R, typename P, typename C>
void matmul(const Matrix<R, P>& A, const Matrix<P, C>& B, Matrix<R, C>& result) {
  const size_t  row = A.row.extent;
  const size_t pair = A.col.extent;
  const size_t  col = B.col.extent;

  std::ranges::fill(result.cell, 0.0f);
  for (size_t i = 0; i < row; ++i)
    for (size_t k = 0; k < pair; ++k) {
      const float a = A.cell[i*pair + k];
      for (size_t j = 0; j < col; ++j) result.cell[i*col +j] += a*B.cell[k*col +j];
    }
}

template <typename P, typename R, typename C>
void matmul_row(const Matrix<P, R>& A, const Matrix<P, C>& B, Matrix<R, C>& result) {
  const size_t pair = A.row.extent;
  const size_t  row = A.col.extent;
  const size_t  col = B.col.extent;

  std::ranges::fill(result.cell, 0.0f);
  for (size_t k = 0; k < pair; ++k)
    for (size_t i = 0; i < row; ++i) {
      const float a = A.cell[k*row + i];
      for (size_t j = 0; j < col; ++j)
        result.cell[i*col +j] += a*B.cell[k*col +j];
    }
}

template <typename R, typename P, typename C>
void matmul_col(const Matrix<R, P>& A, const Matrix<C, P>& B, Matrix<R, C>& result) {
  const size_t  row = A.row.extent;
  const size_t  col = B.row.extent;
  const size_t pair = A.col.extent;

  for (size_t i = 0; i < row; ++i)
    for (size_t j = 0; j < col; ++j) {
      float total = 0.0f;
      for (size_t k = 0; k < pair; ++k) total += A.cell[i*pair +k] * B.cell[j*pair +k];
      result.cell[i*col + j] = total;
    }
}

template <typename R, typename K, typename C>
void normalize_col(Matrix<R, K>& W, Matrix<K, C>& H, float eps) {
  const size_t  row = W.row.extent;
  const size_t rank = W.col.extent;
  const size_t  col = H.col.extent;

  for (size_t i = 0; i < rank; ++i) {
    float energy = 0.0f;
    for (size_t j = 0; j < row; ++j) energy += W.cell[j*rank +i] * W.cell[j*rank +i];
    const float norm = std::sqrt(energy) + eps;
    for (size_t j = 0; j < row; ++j) W.cell[j*rank +i] /= norm;
    for (size_t j = 0; j < col; ++j) H.cell[i*col  +j] *= norm;
  }
}

template <typename R, typename K, typename C>
float non_neg_matrix_factorization(const Matrix<R, C>& V, Matrix<R, K>& W, Matrix<K, C>& H,
                                   int iteration, float tol, float eps) {
  const size_t hz_extent = V.row.extent;
  const size_t     frame = V.col.extent;
  const size_t      rank = W.col.extent;
  const size_t     plane = hz_extent*frame;

  double energy = 0.0;
  for (size_t j = 0; j < plane; ++j) {
    const float e = V.cell[j];
    energy += static_cast<double>(e)*e;
  }

  float prev = 0.0f;

  Matrix<K, C> WtV;
  WtV.row = W.col;
  WtV.col = V.col;
  WtV.cell.resize(rank*frame);

  Matrix<K, K> WtW;
  WtW.row = W.col;
  WtW.col = W.col;
  WtW.cell.resize(rank*rank);

  Matrix<K, C> WtWH;
  WtWH.row = W.col;
  WtWH.col = H.col;
  WtWH.cell.resize(rank*frame);

  Matrix<R, K> VHt;
  VHt.row = V.row;
  VHt.col = H.row;
  VHt.cell.resize(hz_extent*rank);

  Matrix<K, K> HHt;
  HHt.row = H.row;
  HHt.col = H.row;
  HHt.cell.resize(rank*rank);

  Matrix<R, K> WHHt;
  WHHt.row = W.row;
  WHHt.col = H.row;
  WHHt.cell.resize(hz_extent*rank);

  const size_t rank_frame = WtV.cell.size();
  const size_t    hz_rank = hz_extent*rank;
  const size_t  rank_rank = WtW.cell.size();

  for (int i = 0; i < iteration; ++i) {
    matmul_row(W, V, WtV);
    matmul_row(W, W, WtW);
    matmul(WtW, H, WtWH);
    for (size_t j = 0; j < rank_frame; ++j)
      H.cell[j] *= WtV.cell[j] / (WtWH.cell[j] +eps);

    matmul_col(V, H, VHt);
    matmul_col(H, H, HHt);
    matmul(W, HHt, WHHt);
    for (size_t j = 0; j < hz_rank; ++j)
      W.cell[j] *= VHt.cell[j] / (WHHt.cell[j] +eps);

    normalize_col(W, H, eps);
    matmul_row(W, V, WtV);
    matmul_row(W, W, WtW);
    matmul_col(H, H, HHt);
    double cross = 0.0;
    for (size_t j = 0; j < rank_frame; ++j)
      cross += static_cast<double>(WtV.cell[j]) * H.cell[j];
    double quad = 0.0;
    for (size_t j = 0; j < rank_rank; ++j)
      quad += static_cast<double>(WtW.cell[j]) * HHt.cell[j];
    double raw = energy - 2.0*cross +quad;
    float err = static_cast<float>(std::sqrt(std::max(raw, 0.0) / energy));
    if (i > 0 && prev - err < tol*prev) {
      prev = err;
      break;
    }
    prev = err;
  }
  return prev;
}
}
