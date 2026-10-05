#pragma once

#include <vector>

#include "algorithm/axis_quantizer.h"

namespace qstorm::atom::thunder {

using ms = algorithm::ms<double>;
using db = algorithm::db<double>;
inline thread_local algorithm::Spectrogram<double> window;

namespace x {
inline thread_local const algorithm::Quantization<double>& axis = window.x_ms;
inline thread_local ms start{0,   &axis};
inline thread_local ms cross{0,   &axis};
inline thread_local ms onset{0,   &axis};
inline thread_local ms peak{0,    &axis};
inline thread_local ms sustain{0, &axis};
inline thread_local ms end{0,     &axis};
constexpr struct {operator double() const { return window.end(); }} extent{};
constexpr struct {operator double() const { return peak - onset; }} rise{};
}  // namespace x

namespace y {
inline thread_local const algorithm::Quantization<double>& axis = window.y_db;
inline thread_local db floor{0,   &axis};
inline thread_local db cross{0,   &axis};
inline thread_local db onset{0,   &axis};
inline thread_local db peak{0,    &axis};
inline thread_local db sustain{0, &axis};
constexpr struct {db operator[](ms i)  const { return window[i];    }}   pos{};
constexpr struct {   operator double() const { return peak - onset; }}  rise{};
constexpr struct {   operator double() const { return peak - floor; }} range{};
}  // namespace y

namespace xy {
inline thread_local k::Point<const ms&, const db&> cross{x::cross, y::cross};
inline thread_local k::Point<const ms&, const db&> onset{x::onset, y::onset};
inline thread_local k::Point<const ms&, const db&> peak{x::peak, y::peak};
inline thread_local k::Point<const ms&, const db&> sustain{x::sustain, y::sustain};
inline thread_local k::Point<std::vector<ms>, std::vector<db>> peak_run;
inline thread_local k::Point<std::vector<ms>, std::vector<db>> rise_run;
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
};

struct Root {
  ms start = x::start;
  ms   end = x::end;
  db floor = y::floor;
  k::Point<ms, db> cross{x::cross, y::cross};
  k::Point<ms, db> onset{x::onset, y::onset};
  k::Point<ms, db> peak{x::peak, y::peak};
  k::Point<ms, db> sustain{x::sustain, y::sustain};
};

double hold(double span_x, db peak, db floor, double scale, double decay);
double reach(double climb_x, db peak, db floor, double exponent);
bool   rising(db level, db floor, double split);

void absolute_attack(const Factor& factor);
ms   fractional_threshold(Factor& factor, double patience);
void ascent(Factor& factor);
bool descent(const Factor& factor);
void sustain(const Factor& factor);
void fractional_ascent(const Factor& factor);
bool noise(ms from, ms to, const Factor& factor);

std::vector<Root> root(Factor factor);
}
