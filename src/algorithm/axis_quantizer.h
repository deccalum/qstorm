#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <span>
#include <utility>
#include <vector>

#include "algorithm/math_constant.h"

namespace qstorm::algorithm {

template <typename X, typename Y, typename Z>
struct Spectrogram;

template <typename T>
struct Quantization {
  size_t extent;
  T      density;
  T      offset;

  static Quantization stepped(size_t extent, T step, T offset) {
    const T density = T{1} / step;
    return {extent, density, offset};
  }

  template <typename U>
  Quantization<U> as() const {
    const U density_as = static_cast<U>(density);
    const U offset_as = static_cast<U>(offset);
    return {extent, density_as, offset_as};
  }

  size_t frame(T duration) const {
    const T scaled = duration*density;
    return static_cast<size_t>(std::round(scaled));
  }

 private:
  template <typename Discrete = int64_t>
  constexpr Discrete quantize(T continuous) const {
    return static_cast<Discrete>(std::round((continuous - offset) * density));
  }

  template <typename Discrete>
  constexpr T dequantize(Discrete discrete) const {
    return static_cast<T>(discrete) / density + offset;
  }

  template <typename U>
  friend U grid(U value, const Quantization<U>& a, int64_t n);
  template <typename U>
  friend U fraction(const Quantization<U>& a);
  template <typename X, typename Z>
  friend size_t cell(X x, Z z, const Quantization<X>& x_ms, const Quantization<Z>& z_hz);
  template <typename X>
  friend Quantization<X> axis(const Quantization<X>& x, X axis_a, X axis_b);
  template <typename X, typename Y, typename Z>
  friend struct Spectrogram;
};

template <typename T>
T grid(T value, const Quantization<T>& a, int64_t n) {
  const int64_t i = a.quantize(value);
  return a.dequantize(i + n);
}

template <typename T>
T on_axis(T value, const Quantization<T>& a) {
  const bool silent = value == -std::numeric_limits<T>::infinity();
  return silent ? a.offset : value;
}

template <typename T>
T fraction(const Quantization<T>& a) {
  const T  first = a.dequantize(0);
  const T second = a.dequantize(1);
  return second - first;
}

template <typename X, typename Z>
size_t cell(X x, Z z, const Quantization<X>& x_ms, const Quantization<Z>& z_hz) {
  const size_t   f = x_ms.template quantize<size_t>(x);
  const size_t   k = z_hz.template quantize<size_t>(z);
  const size_t bin = z_hz.extent;
  return f*bin + k;
}

template <typename X>
Quantization<X> axis(const Quantization<X>& x, X axis_a, X axis_b) {
  const int64_t lo = x.quantize(axis_a);
  const int64_t hi = x.quantize(axis_b);
  const size_t  extent = static_cast<size_t>(hi - lo);
  const X       offset = x.dequantize(lo);
  return {extent, x.density, offset};
}

template <typename Y>
std::vector<Y*> rows(const std::vector<Y*>& channel, size_t first) {
  std::vector<Y*> moved;
  for (Y* row : channel) moved.push_back(row + first);
  return moved;
}

template <typename Derived, typename T>
struct grid_step {
  Derived& next() {
    Derived& self = static_cast<Derived&>(*this);
    T& value = self.*Derived::value_member;
    const Quantization<T>& axis = *(self.*Derived::axis_member);
    value = grid(value, axis, +1);
    return self;
  }
  Derived& prev() {
    Derived& self = static_cast<Derived&>(*this);
    T& value = self.*Derived::value_member;
    const Quantization<T>& axis = *(self.*Derived::axis_member);
    value = grid(value, axis, -1);
    return self;
  }

  Derived& operator++() { return next(); }
  Derived& operator--() { return prev(); }
  Derived operator++(int) {
    const Derived before = static_cast<Derived&>(*this);
    next();
    return before;
  }
  Derived operator--(int) {
    const Derived before = static_cast<Derived&>(*this);
    prev();
    return before;
  }
};

template <typename X>
struct ms : grid_step<ms<X>, X> {
  X x;
  const Quantization<X>* x_ms;

  static constexpr X ms::* value_member = &ms::x;
  static constexpr const Quantization<X>* ms::* axis_member = &ms::x_ms;

  ms(X value, const Quantization<X>* axis) : x(value), x_ms(axis) {}
  operator X() const { return x; }
  X time_fraction() const { return fraction(*x_ms); }

  ms& operator=(X value) {
    x = value;
    return *this;
  }
};

template <typename Y>
struct db : grid_step<db<Y>, Y> {
  Y y;
  const Quantization<Y>* y_db;

  static constexpr Y db::* value_member = &db::y;
  static constexpr const Quantization<Y>* db::* axis_member = &db::y_db;

  static Y ratio(double linear) { return static_cast<Y>(20.0*std::log10(linear)); }

  db(Y value, const Quantization<Y>* axis) : y(on_axis(value, *axis)), y_db(axis) {}

  operator Y() const { return y; }

  db& operator=(Y value) {
    y = on_axis(value, *y_db);
    return *this;
  }
};

template <typename Y>
Y vol_db(size_t i, size_t count) {
  const double fraction = static_cast<double>(i) / static_cast<double>(count - 1);
  return db<Y>::ratio(fraction);
}

template <typename Y, size_t... i> std::array<db<Y>, sizeof...(i)>
gain_step(const Quantization<Y>& y, std::index_sequence<i...>) {
  return {db<Y>(vol_db<Y>(i, sizeof...(i)), &y)...};
}

template <typename Z>
struct hz : grid_step<hz<Z>, Z> {
  Z z;
  const Quantization<Z>* z_hz;

  static constexpr Z hz::* value_member = &hz::z;
  static constexpr const Quantization<Z>* hz::* axis_member = &hz::z_hz;

  hz(Z value, const Quantization<Z>* axis) : z(value), z_hz(axis) {}

  operator Z() const { return z; }

  hz& operator=(Z value) {
    z = value;
    return *this;
  }
};

template <typename X, typename Y = X, typename Z = X>
struct Spectrogram {
  std::vector<Y*> channel;
  Quantization<X> x_ms;
  Quantization<Y> y_db;
  Quantization<Z> z_hz;

  ms<X> begin() const { return {x_ms.offset, &x_ms}; }
  ms<X> end() const { return {x_ms.dequantize(x_ms.extent), &x_ms}; }

  db<Y> operator[](X x) const {
    const size_t i = cell(x, z_hz.offset, x_ms, z_hz);
    return {channel.front()[i], &y_db};
  }
  db<Y> operator[](X x, Z z) const {
    const size_t i = cell(x, z, x_ms, z_hz);
    return {channel.front()[i], &y_db};
  }

  hz<Z> z_begin() const { return {z_hz.offset, &z_hz}; }
  hz<Z> z_end() const { return {z_hz.dequantize(z_hz.extent), &z_hz}; }

  Spectrogram window(X axis_a, X axis_b) const {
    const size_t first = cell(axis_a, z_hz.offset, x_ms, z_hz);
    return {rows(channel, first), axis(x_ms, axis_a, axis_b), y_db, z_hz};
  }

  Spectrogram frame(X x) const {
    const X next = grid(x, x_ms, 1);
    return window(x, next);
  }
};

template <typename T>
std::vector<std::vector<T>> deinterleave(std::span<const T> interleaved, size_t count) {
  const size_t extent = interleaved.size() / count;
  std::vector<std::vector<T>> channel(count, std::vector<T>(extent));
  for (size_t i = 0; i < extent; i++)
    for (size_t j = 0; j < count; j++) channel[j][i] = interleaved[i*count + j];
  return channel;
}

template <typename S, typename Y>
std::vector<S> mix(std::span<const S> interleaved, size_t count, std::span<const Y> gain, bool mono) {
  const size_t extent = interleaved.size() / count;
  std::vector<S> mix(interleaved.size());
  for (size_t i = 0; i < extent; i++) {
    Y sum = 0;
    for (size_t j = 0; j < count; j++) sum += static_cast<Y>(interleaved[i*count + j]);
    const Y mean = sum / static_cast<Y>(count);
    for (size_t j = 0; j < count; j++) {
      const Y source = mono ? mean : static_cast<Y>(interleaved[i*count + j]);
      mix[i*count + j] = static_cast<S>(source*gain[j]);
    }
  }
  return mix;
}

template <typename S>
int shift(int bits) {
  return std::numeric_limits<S>::digits + 1 - bits;
}

template <typename Y>
Y depth_floor(int bits) {
  return static_cast<Y>(-static_cast<double>(bits) * k::DbPerBit);
}

template <typename Y, typename S>
std::vector<Y> amplitude(const std::vector<S>& interleaved) {
  const Y full_scale = std::ldexp(Y{1}, std::numeric_limits<S>::digits);
  std::vector<Y> scaled(interleaved.size());
  for (size_t i = 0; i < interleaved.size(); i++) scaled[i] = static_cast<Y>(interleaved[i]) / full_scale;
  return scaled;
}

template <typename Y>
Quantization<Y> y_db(const std::vector<std::vector<Y>>& channel, Y floor) {
  double min = 0.0;
  for (const std::vector<Y>& row : channel)
    for (Y value : row)
      if (std::isfinite(value)) min = std::min(min, static_cast<double>(value));
  const double density = 1.0;
  const double span = std::ceil(-min*density);
  return {static_cast<size_t>(span), static_cast<Y>(density), floor};
}

template <typename Y>
Quantization<Y> y_db(const std::vector<std::vector<Y>>& channel, const Quantization<Y>& parent) {
  return y_db(channel, parent.offset);
}

template <typename X, typename Y, typename Z>
std::span<Y> channel(const Spectrogram<X, Y, Z>& spectrum, size_t c) {
  return std::span<Y>(spectrum.channel[c], spectrum.x_ms.extent*spectrum.z_hz.extent);
}

template <typename X, typename Y, typename Z>
Spectrogram<X, Y, Z> as_spectrogram(std::vector<std::vector<Y>>& channel, Quantization<X> x_ms,
                                    Quantization<Z> z_hz, Quantization<Y> y_db) {
  Spectrogram<X, Y, Z> spectrum;
  for (std::vector<Y>& row : channel) spectrum.channel.push_back(row.data());
  spectrum.x_ms = x_ms;
  spectrum.z_hz = z_hz;
  spectrum.y_db = y_db;
  return spectrum;
}

}
