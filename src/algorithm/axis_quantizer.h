#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <span>
#include <type_traits>
#include <utility>

#include "algorithm/math_constant.h"

namespace qstorm::algorithm {
template <typename X, typename Y, typename Z>
struct Spectrum;

template <typename T>
struct Qtz {
  size_t n;    // Number of grid points
  T      rho;  // Density Per Axis - grid points per axis unit
  T      pos;  // Axis coordinate at grid index 0

  static Qtz stepped(size_t n, T step, T pos) {
    const T density = T{1} / step;
    return {n, density, pos};
  }

  template <typename U>
  Qtz<U> as() const {
    const U rho_as = rho;
    const U pos_as = pos;
    return {n, rho_as, pos_as};
  }

  size_t frame(T dx) const {
    const T scaled = dx * rho;
    return std::round(scaled);
  }

 private:
  template <typename Discrete = int64_t>
  constexpr Discrete qtz(T continuous) const {
    return std::round((continuous - pos) * rho);
  }
  template <typename Discrete>
  constexpr T dqtz(Discrete discrete) const {
    return discrete / rho + pos;
  }

  template <typename U>
  friend U grid(U val, const Qtz<U>& q, int64_t n);
  template <typename U>
  friend U fraction(const Qtz<U>& q);
  template <typename X, typename Z>
  friend size_t cell(X x, Z z, const Qtz<X>& qx, const Qtz<Z>& qz);
  template <typename X>
  friend Qtz<X> axis(const Qtz<X>& qx, X x0, X x1);
  template <typename X, typename Y, typename Z>
  friend struct Spectrum;
};

template <typename T>
T grid(T val, const Qtz<T>& q, int64_t n) {
  const int64_t i = q.qtz(val);
  return q.dqtz(i + n);
}

template <typename T>
T on_axis(T val, const Qtz<T>& q) {
  const bool silent = val == -k::INF<T>;
  return silent ? q.pos : val;
}

template <typename T>
T fraction(const Qtz<T>& q) {
  const T x0 = q.dqtz(0);
  const T x1 = q.dqtz(1);
  return x1 - x0;
}

template <typename X, typename Z>
size_t cell(X x, Z z, const Qtz<X>& qx, const Qtz<Z>& qz) {
  const size_t i  = qx.template qtz<size_t>(x);
  const size_t k  = qz.template qtz<size_t>(z);
  const size_t nz = qz.n;
  return i * nz + k;
}

template <typename X>
Qtz<X> axis(const Qtz<X>& qx, X x0, X x1) {
  const int64_t i0  = qx.qtz(x0);
  const int64_t i1  = qx.qtz(x1);
  const size_t  n   = i1 - i0;
  const X       pos = qx.dqtz(i0);
  return {n, qx.rho, pos};
}

template <typename Y>
vec<Y*> rows(const vec<Y*>& ch, size_t i0) {
  vec<Y*> moved;
  for (Y* row : ch) moved.push_back(row + i0);
  return moved;
}

template <typename Derived, typename T>
struct grid_step {
  Derived& next() {
    Derived&      self = static_cast<Derived&>(*this);
    T&            val  = self.*Derived::val_member;
    const Qtz<T>& q    = *(self.*Derived::axis_member);
    val                = grid(val, q, +1);
    return self;
  }
  Derived& prev() {
    Derived&      self = static_cast<Derived&>(*this);
    T&            val  = self.*Derived::val_member;
    const Qtz<T>& q    = *(self.*Derived::axis_member);
    val                = grid(val, q, -1);
    return self;
  }

  Derived& operator++() { return next(); }
  Derived& operator--() { return prev(); }
  Derived  operator++(int) {
    const Derived v0 = static_cast<Derived&>(*this);
    next();
    return v0;
  }
  Derived operator--(int) {
    const Derived v0 = static_cast<Derived&>(*this);
    prev();
    return v0;
  }
};

template <typename X>
struct ms : grid_step<ms<X>, X> {
  X                  x;
  const Qtz<X>*      qx;
  static constexpr X ms::*       val_member        = &ms::x;
  static constexpr const Qtz<X>* ms::* axis_member = &ms::qx;

  ms() = default;
  ms(X val, const Qtz<X>* q) : x(val), qx(q) {}
      operator X() const { return x; }
  X   time_fraction() const { return fraction(*qx); }
  ms& operator=(X val) {
    x = val;
    return *this;
  }
};

template <typename Y>
struct db : grid_step<db<Y>, Y> {
  Y                  y;
  const Qtz<Y>*      qy;
  static constexpr Y db::*       val_member        = &db::y;
  static constexpr const Qtz<Y>* db::* axis_member = &db::qy;

  static Y ratio(double lin) { return 20.0 * std::log10(lin); }
  static Y lin(Y db) { return std::pow(10.0, db / 10.0); }
  db(Y val, const Qtz<Y>* q) : y(on_axis(val, *q)), qy(q) {}
      operator Y() const { return y; }
  db& operator=(Y val) {
    y = on_axis(val, *qy);
    return *this;
  }
};

template <typename Y>
Y vol_db(size_t i, size_t n) {
  const double i1 = n - 1;
  const double t  = i / i1;
  return db<Y>::ratio(t);
}

template <typename Y, size_t... i>
std::array<db<Y>, sizeof...(i)> gain_step(const Qtz<Y>& qy, std::index_sequence<i...>) {
  return {db<Y>(vol_db<Y>(i, sizeof...(i)), &qy)...};
}

template <typename Z>
struct hz : grid_step<hz<Z>, Z> {
  Z                  z;
  const Qtz<Z>*      qz;
  static constexpr Z hz::*       val_member        = &hz::z;
  static constexpr const Qtz<Z>* hz::* axis_member = &hz::qz;

  hz(Z val, const Qtz<Z>* q) : z(val), qz(q) {}
      operator Z() const { return z; }
  hz& operator=(Z val) {
    z = val;
    return *this;
  }
};

template <typename X, typename Y = X, typename Z = X>
struct Spectrum {
  vec<Y*> ch;
  Qtz<X>  qx;
  Qtz<Y>  qy;
  Qtz<Z>  qz;
  ms<X>   x0() const { return {qx.pos, &qx}; }
  ms<X>   x1() const { return {qx.dqtz(qx.n), &qx}; }

  db<Y> operator[](X x, size_t c = 0) const {
    const size_t i = cell(x, qz.pos, qx, qz);
    return {ch[c][i], &qy};
  }
  db<Y> operator[](X x, Z z, size_t c = 0) const {
    const size_t i = cell(x, z, qx, qz);
    return {ch[c][i], &qy};
  }

  hz<Z> z0() const { return {qz.pos, &qz}; }
  hz<Z> z1() const { return {qz.dqtz(qz.n), &qz}; }

  Spectrum window(X x0, X x1) const {
    const size_t i0 = cell(x0, qz.pos, qx, qz);
    return {rows(ch, i0), axis(qx, x0, x1), qy, qz};
  }
  Spectrum frame(X x0) const {
    const X x1 = grid(x0, qx, 1);
    return window(x0, x1);
  }

  Spectrum select(size_t c) const { return {{ch[c]}, qx, qy, qz}; }
};

template <typename X, typename Y = X, typename Z = X>
using zvec = vec<Spectrum<X, Y, Z>>;

template <typename X, typename Y, typename Z>
zvec<X, Y, Z> select_ch(const zvec<X, Y, Z>& band, size_t ch) {
  zvec<X, Y, Z> c;
  for (const Spectrum<X, Y, Z>& each : band) c.push_back(each.select(ch));
  return c;
}

template <typename X, typename Y, typename Z>
vec<Qtz<Z>> band_axis(const zvec<X, Y, Z>& band) {
  vec<Qtz<Z>> qz;
  for (const Spectrum<X, Y, Z>& each : band) qz.push_back(each.qz);
  return qz;
}

template <typename T>
vvec<T> xc_to_cx(std::span<const T> xc, size_t nc) {
  const size_t nx = xc.size() / nc;

  vvec<T> cx(nc, vec<T>(nx));
  for (size_t i = 0; i < nx; i++)
    for (size_t c = 0; c < nc; c++) cx[c][i] = xc[i * nc + c];
  return cx;
}

template <typename S, typename Y>
vec<S> xc_mix(std::span<const S> xc, size_t nc, std::span<const Y> gain, bool mono) {
  const size_t nx = xc.size() / nc;
  vec<S>       mix(xc.size());
  for (size_t i = 0; i < nx; i++) {
    Y sum = 0;
    for (size_t c = 0; c < nc; c++) sum += static_cast<Y>(xc[i * nc + c]);
    const Y mean = sum / static_cast<Y>(nc);
    for (size_t c = 0; c < nc; c++) {
      const Y source  = mono ? mean : static_cast<Y>(xc[i * nc + c]);
      mix[i * nc + c] = static_cast<S>(source * gain[c]);
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
  return -bits * k::DbPerBit;
}

template <typename Y, typename S>
Y full_scale() {
  return std::ldexp(Y{1}, std::numeric_limits<S>::digits);
}

template <typename Y, typename S>
vec<Y> amp(const vec<S>& xc) {
  const Y full = full_scale<Y, S>();
  vec<Y>  scale(xc.size());
  for (size_t i = 0; i < xc.size(); i++) scale[i] = static_cast<Y>(xc[i]) / full;
  return scale;
}

template <typename Y, typename S>
vvec<Y> ch_amp(const vec<S>& xc, size_t nc) {
  const vec<Y> scale = amp<Y>(xc);
  return xc_to_cx<Y>(scale, nc);
}

template <typename S, typename X>
std::span<const S> xc_window(const vec<S>& xc, size_t nc, X hz, std::type_identity_t<X> x0,
                             std::type_identity_t<X> x1) {
  const size_t nx   = xc.size() / nc;
  const X      xmax = nx / hz;
  const Qtz<X> qx{nx, hz, X{0}};
  const size_t i0 = qx.frame(x0);
  const size_t i1 = qx.frame(std::min(x1, xmax));
  return std::span<const S>(xc).subspan(i0 * nc, (i1 - i0) * nc);
}

template <typename S, typename Y>
vec<S> xc_qtzd(const vvec<Y>& cx) {
  const Y      full = full_scale<Y, S>();
  const Y      y1   = static_cast<Y>(std::numeric_limits<S>::max());
  const Y      y0   = static_cast<Y>(std::numeric_limits<S>::min());
  const size_t nc   = cx.size();
  const size_t nx   = cx.front().size();
  vec<S>       xc(nc * nx);
  for (size_t i = 0; i < nx; i++)
    for (size_t c = 0; c < nc; c++)
      xc[i * nc + c] = static_cast<S>(std::clamp(std::round(cx[c][i] * full), y0, y1));
  return xc;
}

template <typename Y>
Qtz<Y> ydb(const vvec<Y>& channel, Y floor) {
  double y0 = 0.0;
  for (const vec<Y>& row : channel)
    for (Y val : row)
      if (std::isfinite(val)) y0 = std::min<double>(y0, val);
  const Y      density = 1.0;
  const size_t n       = std::ceil(-y0 * density);
  return {n, density, floor};
}

template <typename Y>
Qtz<Y> ydb(const vvec<Y>& ch, const Qtz<Y>& parent) {
  return ydb(ch, parent.pos);
}

template <typename X, typename Y, typename Z>
std::span<Y> ch(const Spectrum<X, Y, Z>& env, size_t c) {
  return std::span<Y>(env.ch[c], env.qx.n * env.qz.n);
}

template <typename X, typename Y, typename Z>
Spectrum<X, Y, Z> as_spectrum(vvec<Y>& ch, Qtz<X> qx, Qtz<Y> qy, Qtz<Z> qz) {
  Spectrum<X, Y, Z> spectrum;
  for (vec<Y>& row : ch) spectrum.ch.push_back(row.data());
  spectrum.qx = qx;
  spectrum.qy = qy;
  spectrum.qz = qz;
  return spectrum;
}

template <typename X, typename Y, typename Z, typename S>
Spectrum<X, Y, Z> xc_spectrum(const vec<S>& xc, size_t nc, X hz, std::type_identity_t<X> x0,
                              Qtz<Y> qy, Qtz<Z> qz, vvec<Y>& amp) {
  const size_t nx = xc.size() / nc;
  const Qtz<X> qx{nx, hz, x0};
  amp = ch_amp<Y>(xc, nc);
  return as_spectrum(amp, qx, qy, qz);
}
}

namespace qstorm {
using algorithm::ms;
using algorithm::Qtz;
using algorithm::Spectrum;
using algorithm::zvec;
using Level = Spectrum<double, float, double>;
}
