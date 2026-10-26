#pragma once

#include <optional>

#include "algorithm/axis_quantizer.h"
#include "atom/atom_constant.h"

namespace qstorm::atom::thunder {
using ms = algorithm::ms<double>;
using db = algorithm::db<double>;
inline thread_local Spectrum<double> window;

namespace x {
inline thread_local const Qtz<double>& axis = window.qx;
inline thread_local ms                 start{0, &axis};
inline thread_local ms                 cross{0, &axis};
inline thread_local ms                 onset{0, &axis};
inline thread_local ms                 peak{0, &axis};
inline thread_local ms                 sustain{0, &axis};
inline thread_local ms                 end{0, &axis};

constexpr struct {
  operator double() const { return window.x1(); }
} extent{};
constexpr struct {
  operator double() const { return peak - onset; }
} rise{};
}  // namespace x

namespace y {
inline thread_local const Qtz<double>& axis = window.qy;
inline thread_local db                 floor{0, &axis};
inline thread_local db                 cross{0, &axis};
inline thread_local db                 onset{0, &axis};
inline thread_local db                 peak{0, &axis};
inline thread_local db                 sustain{0, &axis};
inline thread_local db                 noise{0, &axis};

constexpr struct {
  db operator[](ms i) const { return window[i]; }
} pos{};
constexpr struct {
  operator double() const { return peak - onset; }
} rise{};
constexpr struct {
  operator double() const { return peak - floor; }
} range{};
}  // namespace y

namespace xy {
inline thread_local k::Point<const ms&, const db&> cross{x::cross, y::cross};
inline thread_local k::Point<const ms&, const db&> onset{x::onset, y::onset};
inline thread_local k::Point<const ms&, const db&> peak{x::peak, y::peak};
inline thread_local k::Point<const ms&, const db&> sustain{x::sustain, y::sustain};
inline thread_local k::Point<vec<ms>, vec<db>>     peak_run;
inline thread_local k::Point<vec<ms>, vec<db>>     rise_run;
}  // namespace xy

struct Factor {
  double decline;
  double tolerance;
  double split;
  double scale;
  double decay;
  double exponent;
  double range;
  double patience;
  double release;
  double cap;
  double end;
  double fall;

  std::optional<double> silence;
  std::optional<double> settle;
};

struct Root {
  ms               walk  = x::start;
  db               floor = y::floor;
  k::Point<ms, db> cross{x::cross, y::cross};
  k::Point<ms, db> onset{x::onset, y::onset};
  k::Point<ms, db> peak{x::peak, y::peak};
  k::Point<ms, db> sustain{x::sustain, y::sustain};
};

double hold(double span_x, db peak, db floor, double scale, double decay);
double reach(double climb_x, db peak, db floor, double exponent);
bool   rising(db level, db floor, double split);
db     silence(db floor, double factor);
void   absolute_attack(const Factor& factor);
ms     fractional_threshold(Factor& factor, double patience);
void   ascent(Factor& factor);
bool   descent(const Factor& factor);
void   sustain(const Factor& factor);
void   fractional_ascent(const Factor& factor);
bool   noise(ms from, ms to, const Factor& factor);

vec<Root> root(const Spectrum<double>& spectrum, Factor factor, vec<ADSR>& adsr);
}
