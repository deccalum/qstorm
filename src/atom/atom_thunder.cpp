#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <vector>

#include "algorithm/math_scalar.h"
#include "atom/atom_thunder.h"

namespace qstorm::atom::thunder {

double hold(double span_x, db peak, db floor, double scale, double decay) {
  const double y_rise = peak - floor;
  double rise_scl = y_rise*y_rise*scale;
  return rise_scl * (1.0 +decay / std::log(span_x +1.0));
}

double reach(double climb_x, db peak, db floor, double exponent) {
  const double   axis = peak.y_db->offset;
  const double height = peak - axis;
  const double  depth = floor - axis;
  return height == 0.0 ? climb_x : climb_x * std::pow(height / depth, -exponent);
}

bool rising(db level, db floor, double split) {
  const double   axis = floor.y_db->offset;
  const double height = level - axis;
  const double  depth = floor - axis;
  return height / depth >= split;
}

void fractional_ascent(const Factor& factor) {
  double total = 0;
  size_t n = 0;
  double rise_total = 0;
  xy::rise_run.x.clear();
  xy::rise_run.y.clear();

  for (ms i = x::start; i < x::extent; ++i) {
    const db level = y::pos[i];
    total += level;
    ++n;
    y::floor = total / n;

    xy::rise_run.x.push_back(i);
    xy::rise_run.y.push_back(level);
    rise_total += level;
    y::cross = rise_total / xy::rise_run.y.size();

    const double elapsed = i - x::start;
    const double held = i - xy::rise_run.x.front();
    const double frac = hold(elapsed, y::cross, y::floor, factor.scale, factor.decay);
    if (held < frac) continue;

    if (rising(y::cross, y::floor, factor.split)
        && !noise(xy::rise_run.x.front(), i, factor)) {
      x::cross = i;
      return;
    }
    xy::rise_run.x.clear();
    xy::rise_run.y.clear();
    rise_total = 0;
  }

  x::cross = x::extent;
}

bool noise(ms from, ms to, const Factor& factor) {
  const double stretch = to - from;
  double before_total = 0;
  size_t before_n = 0;
  ms i = from;
  for (--i; i >= x::axis.offset && from - i <= stretch; --i) {
    before_total += y::pos[i];
    ++before_n;
  }
  if (before_n == 0) return false;

  const double    before = before_total / before_n;
  const double     range = y::cross - y::floor;
  const double tolerance = factor.range*range;
  return std::abs(y::cross - before) <= tolerance;
}

ms fractional_threshold(Factor& factor, double patience) {
  const double stretch = x::cross - xy::rise_run.x.front();
  double block_total = 0;
  size_t     block_n = 0;
  ms         block_x = x::peak;
  double   run_total = 0;
  size_t       run_n = 0;
  ms           run_x = x::peak;
  bool         above = false;
  bool         below = false;
  xy::peak_run.x.assign(1, x::peak);
  xy::peak_run.y.assign(1, y::peak);
  for (ms i = x::peak; i < x::extent; ++i) {
    block_total += y::pos[i];
    ++block_n;
    if (i - block_x < stretch) continue;

    const double mean = block_total / block_n;
    if (!above && !below) {
      const double y_decline = y::peak - mean;
      const double x_decline = i - x::peak;
      above = mean >= y::peak;
      below = !above && y_decline / y::range > patience / x_decline;
      run_x = block_x;
      run_total = 0;
      run_n = 0;
    }
    if (above || below) {
      run_total += mean;
      ++run_n;
      const db level{run_total / run_n, y::peak.y_db};
      const double length
        = hold(i - x::start, above ? level : y::peak, y::floor, factor.scale, factor.decay);
      if (i - run_x >= length) {
        if (above && level >= y::peak) {
          y::peak = level;
          x::peak = run_x;
          xy::peak_run.x.push_back(run_x);
          xy::peak_run.y.push_back(level);
        }
        if (below) {
          const double y_decline = y::peak - level;
          const double x_decline = i - x::peak;
          factor.decline = y_decline / y::range;
          factor.tolerance = patience / x_decline;
          if (factor.decline > factor.tolerance) return i;
        }
        above = false;
        below = false;
      }
    }
    block_total = 0;
    block_n = 0;
    block_x = i;
  }
  return {x::extent, &x::axis};
}

// Rewind to find the physical start of the transient below the floor, bounded by
// `hold` so onset cannot drift back further than the attack's own climb warrants.
void absolute_attack(const Factor& factor) {
  const double stretch = x::cross - xy::rise_run.x.front();
  const double climb = x::peak - x::cross;
  const double cap = reach(climb, y::peak, y::floor, factor.exponent);

  ms limit{x::axis.offset, &x::axis};
  if (x::cross - cap > limit) limit = x::cross - cap;
  ms     i = x::cross;
  ms     edge = x::cross;
  ms     run_x = x::cross;
  double run_total = 0;
  size_t run_n = 0;
  bool   quiet = false;
  bool   held = false;
  double block_total = 0;
  size_t block_n = 0;
  double mean = y::cross;
  while (i > limit) {
    block_total += y::pos[i];
    ++block_n;
    if (edge - i >= stretch) {
      mean = block_total / block_n;
      if (!quiet && mean <= y::floor) {
        quiet = true;
        run_x = edge;
        run_total = 0;
        run_n = 0;
      }
      if (quiet) {
        run_total += mean;
        ++run_n;
        const double length = hold(x::cross - i, y::peak, y::floor, factor.scale, factor.decay);
        if (run_x - i >= length) {
          held = run_total / run_n <= y::floor;
          if (held) break;
          quiet = false;
        }
      }
      block_total = 0;
      block_n = 0;
      edge = i;
    }
    --i;
  }
  x::onset = held ? run_x : edge;
  y::onset = held ? run_total / run_n : mean;
}

void ascent(Factor& factor) {
  x::peak = x::cross;
  y::peak = y::cross;

  const double patience_cross = (x::cross - x::start) * factor.patience;
  fractional_threshold(factor, patience_cross);
  absolute_attack(factor);

  const double patience_onset = x::rise*factor.patience;
  fractional_threshold(factor, patience_onset);
}

bool descent(const Factor& factor) {
  const double stretch = x::cross - xy::rise_run.x.front();
  const double length = hold(x::rise + stretch, y::peak, y::floor, factor.fall, factor.decay);
  const double window_end = x::peak + length;
  while (true) {
    const double held
      = std::max(stretch, hold(x::peak - x::start, y::peak, y::floor, factor.scale, factor.decay));
    std::vector<double> time;
    std::vector<double> level;
    double block_total = 0;
    size_t   block_n = 0;
    ms       block_x = x::peak;
    double run_total = 0;
    size_t     run_n = 0;
    ms         run_x = x::peak;
    double       top = y::peak;
    ms         top_x = x::peak;
    for (ms i = x::peak; i < x::extent; ++i) {
      if (i > window_end) break;
      block_total += y::pos[i];
      ++block_n;
      if (i - block_x < stretch) continue;

      const double mean = block_total / block_n;
      time.push_back(0.5 * (block_x + i));
      level.push_back(mean);
      run_total += mean;
      ++run_n;
      if (i - run_x >= held) {
        const double run_mean = run_total / run_n;
        if (run_mean > top) {
          top = run_mean;
          top_x = run_x;
        }
        run_total = 0;
        run_n = 0;
        run_x = i;
      }
      block_total = 0;
      block_n = 0;
      block_x = i;
    }
    if (time.size() < 2) return false;

    const double     n = static_cast<double>(time.size());
    const double  xbar = algorithm::average(time, std::identity{});
    const double  ybar = algorithm::average(level, std::identity{});
    const double slope = algorithm::slope(time, level);

    if (slope >= 0) {
      if (top <= y::peak) return false;
      y::peak = top;
      x::peak = top_x;
      continue;
    }

    double residual = 0;
    for (size_t i = 0; i < time.size(); ++i) {
      const double yhat = ybar + slope * (time[i] - xbar);
      residual += (level[i] - yhat) * (level[i] - yhat);
    }
    residual = std::sqrt(residual / n);

    const double yhat_drop = -slope * (time.back() - time.front());
    const double at_peak = ybar + slope * (x::peak - xbar);
    const double to_floor = at_peak - y::floor;
    if (yhat_drop <= residual || to_floor <= 0) return false;

    x::end = x::peak + to_floor / -slope;
    return true;
  }
}

void sustain(const Factor& factor) {
  const double stretch = x::cross - xy::rise_run.x.front();
  const double width = factor.release * (x::end - x::peak);
  const double limit = x::end + width;
  double block_total = 0;
  size_t block_n = 0;
  ms block_x = x::peak;
  ms run_x = x::peak;
  double run_total = 0;
  size_t run_n = 0;
  bool quiet = false;
  for (ms i = x::peak; i < x::extent; ++i) {
    if (i > limit) {
      x::sustain = i;
      return;
    }
    block_total += y::pos[i];
    ++block_n;
    if (i - block_x < stretch) continue;

    const double mean = block_total / block_n;
    if (!quiet && mean <= y::floor) {
      quiet = true;
      run_x = block_x;
      run_total = 0;
      run_n = 0;
    }
    if (quiet) {
      run_total += mean;
      ++run_n;
      const double length = hold(i - x::onset, y::peak, y::floor, factor.scale, factor.decay);
      if (i - run_x >= length) {
        y::sustain = run_total / run_n;
        if (y::sustain <= y::floor) {
          x::sustain = run_x;
          return;
        }
        quiet = false;
      }
    }
    block_total = 0;
    block_n = 0;
    block_x = i;
  }
  x::sustain = x::extent;
}

std::vector<Root> root(Factor factor) {
  x::start = x::axis.offset;

  std::vector<Root> out;
  while (x::start < x::extent) {
    fractional_ascent(factor);
    if (x::cross >= x::extent) break;

    ascent(factor);
    bool keep = descent(factor);
    if (keep) {
      absolute_attack(factor);
      sustain(factor);
      const double width = factor.release * (x::end - x::peak);
      keep = std::abs(x::sustain - x::end) <= width;
      x::end = x::sustain;
    }
    if (keep) {
      const double least = hold(x::end - x::onset, y::peak, y::floor, factor.end, factor.decay);
      if (x::end - x::peak < least) x::end = x::peak + least;
      out.push_back(Root{});
    }

    const ms next = keep ? x::end : x::cross;
    if (next > x::start)
      x::start = next;
    else
      ++x::start;
  }
  return out;
}

}
