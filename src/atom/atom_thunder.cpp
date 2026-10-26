#include "atom/atom_thunder.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <span>
#include <vector>

#include "algorithm/dsp_envelope.h"
#include "algorithm/math_scalar.h"

namespace qstorm::atom::thunder {
db silence(db floor, double factor) {
  const double axis  = floor.qy->pos;
  const double depth = y::noise - axis;
  const double bound = axis + factor * depth;
  return {std::max<double>(floor, bound), floor.qy};
}

double hold(double span_x, db peak, db floor, double scale, double decay) {
  const double y_rise   = peak - floor;
  double       rise_scl = y_rise * y_rise * scale;
  return rise_scl * (1.0 + decay / std::log(span_x + 1.0));
}

double reach(double climb_x, db peak, db floor, double exponent) {
  const double axis   = peak.qy->pos;
  const double height = peak - axis;
  const double depth  = floor - axis;
  return height == 0.0 ? climb_x : climb_x * std::pow(height / depth, -exponent);
}

bool rising(db level, db floor, double split) {
  const double axis   = floor.qy->pos;
  const double height = level - axis;
  const double depth  = floor - axis;
  return height / depth >= split;
}

void fractional_ascent(const Factor& factor) {
  double tot      = 0;
  size_t n        = 0;
  double rise_tot = 0;
  xy::rise_run.x.clear();
  xy::rise_run.y.clear();

  for (ms i = x::start; i < x::extent; ++i) {
    const db level = y::pos[i];
    tot += level;
    ++n;
    y::floor = tot / n;
    y::floor = silence(y::floor, factor.silence.value());
    y::floor = std::min<double>(y::floor, y::noise);

    xy::rise_run.x.push_back(i);
    xy::rise_run.y.push_back(level);
    rise_tot += level;
    y::cross = rise_tot / xy::rise_run.y.size();

    const double elapsed = i - x::start;
    const double held    = i - xy::rise_run.x.front();
    const double frac    = hold(elapsed, y::cross, y::floor, factor.scale, factor.decay);
    if (held < frac) continue;

    if (rising(y::cross, y::floor, factor.split) && !noise(xy::rise_run.x.front(), i, factor)) {
      x::cross = xy::rise_run.x.front();
      return;
    }
    xy::rise_run.x.clear();
    xy::rise_run.y.clear();
    rise_tot = 0;
  }

  x::cross = x::extent;
}

bool noise(ms from, ms to, const Factor& factor) {
  const double stretch      = to - from;
  double       before_tot   = 0;
  size_t       before_n     = 0;
  ms           i            = from;
  for (--i; i >= x::axis.pos && from - i <= stretch; --i) {
    before_tot += y::pos[i];
    ++before_n;
  }
  if (before_n == 0) return false;

  const double before    = before_tot / before_n;
  const double range     = y::cross - y::floor;
  const double tolerance = factor.range * range;
  return std::abs(y::cross - before) <= tolerance;
}

ms fractional_threshold(Factor& factor, double patience) {
  const double stretch     = xy::rise_run.x.back() - xy::rise_run.x.front();
  double       blk_tot     = 0;
  size_t       block_n     = 0;
  ms           block_x     = x::peak;
  double       run_tot     = 0;
  size_t       run_n       = 0;
  ms           run_x       = x::peak;
  bool         above       = false;
  bool         below       = false;
  xy::peak_run.x.assign(1, x::peak);
  xy::peak_run.y.assign(1, y::peak);
  for (ms i = x::peak; i < x::extent; ++i) {
    blk_tot += y::pos[i];
    ++block_n;
    if (i - block_x < stretch) continue;

    const double mean = blk_tot / block_n;
    if (!above && !below) {
      const double y_decline = y::peak - mean;
      const double x_decline = i - x::peak;
      above                  = mean >= y::peak;
      below                  = !above && y_decline / y::range > patience / x_decline;
      run_x                  = block_x;
      run_tot                = 0;
      run_n                  = 0;
    }
    if (above || below) {
      run_tot += mean;
      ++run_n;
      const db     level{run_tot / run_n, y::peak.qy};
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
          factor.decline         = y_decline / y::range;
          factor.tolerance       = patience / x_decline;
          if (factor.decline > factor.tolerance) return i;
        }
        above = false;
        below = false;
      }
    }
    blk_tot     = 0;
    block_n     = 0;
    block_x     = i;
  }
  return {x::extent, &x::axis};
}

// Rewind to find the physical start of the transient below the floor, bounded by
// `hold` so onset cannot drift back further than the attack's own climb warrants.
void absolute_attack(const Factor& factor) {
  const double stretch = xy::rise_run.x.back() - xy::rise_run.x.front();
  const double climb   = x::peak - x::cross;
  const double cap     = reach(climb, y::peak, y::floor, factor.exponent);

  ms limit{x::axis.pos, &x::axis};
  if (x::cross - cap > limit) limit = x::cross - cap;
  ms     i           = x::cross;
  ms     edge        = x::cross;
  ms     run_x       = x::cross;
  size_t run_n       = 0;
  bool   quiet       = false;
  bool   held        = false;
  double run_tot     = 0;
  double blk_tot     = 0;
  size_t block_n     = 0;
  double mean        = y::cross;
  while (i > limit) {
    blk_tot += y::pos[i];
    ++block_n;
    if (edge - i >= stretch) {
      mean = blk_tot / block_n;
      if (!quiet && mean <= y::floor) {
        quiet     = true;
        run_x     = edge;
        run_tot   = 0;
        run_n     = 0;
      }
      if (quiet) {
        run_tot += mean;
        ++run_n;
        const double length = hold(x::cross - i, y::peak, y::floor, factor.scale, factor.decay);
        if (run_x - i >= length) {
          held = run_tot / run_n <= y::floor;
          if (held) break;
          quiet = false;
        }
      }
      blk_tot     = 0;
      block_n     = 0;
      edge        = i;
    }
    --i;
  }
  x::onset = held ? run_x : edge;
  y::onset = held ? run_tot / run_n : mean;
}

void ascent(Factor& factor) {
  x::peak = x::cross;
  y::peak = y::cross;

  const double patience_cross = (x::cross - x::start) * factor.patience;
  fractional_threshold(factor, patience_cross);
  absolute_attack(factor);

  const double patience_onset = x::rise * factor.patience;
  fractional_threshold(factor, patience_onset);
}

bool descent(const Factor& factor) {
  const double stretch    = xy::rise_run.x.back() - xy::rise_run.x.front();
  const double length     = hold(x::rise + stretch, y::peak, y::floor, factor.fall, factor.decay);
  const double window_end = x::peak + length;
  while (true) {
    const double held
      = std::max(stretch, hold(x::peak - x::start, y::peak, y::floor, factor.scale, factor.decay));
    vec    xs;
    vec    ys;
    double blk_tot = 0;
    size_t block_n = 0;
    ms     block_x = x::peak;
    double run_tot = 0;
    size_t run_n   = 0;
    ms     run_x   = x::peak;
    double y1      = y::peak;
    ms     arg     = x::peak;
    for (ms i = x::peak; i < x::extent; ++i) {
      if (i > window_end) break;
      blk_tot += y::pos[i];
      ++block_n;
      if (i - block_x < stretch) continue;

      const double mean = blk_tot / block_n;
      xs.push_back((block_x + i) / 2);
      ys.push_back(mean);
      run_tot += mean;
      ++run_n;
      if (i - run_x >= held) {
        const double run_mean = run_tot / run_n;
        if (run_mean > y1) {
          y1  = run_mean;
          arg = run_x;
        }
        run_tot   = 0;
        run_n     = 0;
        run_x     = i;
      }
      blk_tot     = 0;
      block_n     = 0;
      block_x     = i;
    }
    if (xs.size() < 2) return false;

    const double n     = static_cast<double>(xs.size());
    const double xbar  = algorithm::average(xs, std::identity{});
    const double ybar  = algorithm::average(ys, std::identity{});
    const double slope = algorithm::slope(xs, ys);

    if (slope >= 0) {
      if (y1 <= y::peak) return false;
      y::peak = y1;
      x::peak = arg;
      continue;
    }

    double residual = 0;
    for (size_t i = 0; i < xs.size(); ++i) {
      const double yhat = ybar + slope * (xs[i] - xbar);
      residual += (ys[i] - yhat) * (ys[i] - yhat);
    }
    residual = std::sqrt(residual / n);

    const double yhat_drop = -slope * (xs.back() - xs.front());
    const double at_peak   = ybar + slope * (x::peak - xbar);
    const double to_floor  = at_peak - y::floor;
    if (yhat_drop <= residual || to_floor <= 0) return false;

    x::end = x::peak + to_floor / -slope;
    return true;
  }
}

void sustain(const Factor& factor) {
  const double stretch     = xy::rise_run.x.back() - xy::rise_run.x.front();
  const double width       = factor.release * (x::end - x::peak);
  const double limit       = x::end + width;
  double       blk_tot     = 0;
  size_t       block_n     = 0;
  ms           block_x     = x::peak;
  ms           run_x       = x::peak;
  double       run_tot     = 0;
  size_t       run_n       = 0;
  bool         quiet       = false;
  for (ms i = x::peak; i < x::extent; ++i) {
    if (i >= limit) {
      x::sustain = x::end - width;
      return;
    }
    blk_tot += y::pos[i];
    ++block_n;
    if (i - block_x < stretch) continue;

    const db mean{blk_tot / block_n, &y::axis};
    if (!quiet && !rising(mean, y::floor, factor.settle.value())) {
      quiet     = true;
      run_x     = block_x;
      run_tot   = 0;
      run_n     = 0;
    }
    if (quiet) {
      run_tot += mean;
      ++run_n;
      const double length = hold(i - x::onset, y::peak, y::floor, factor.scale, factor.decay);
      if (i - run_x >= length) {
        y::sustain = run_tot / run_n;
        if (!rising(y::sustain, y::floor, factor.settle.value())) {
          x::sustain = run_x;
          return;
        }
        quiet = false;
      }
    }
    blk_tot     = 0;
    block_n     = 0;
    block_x     = i;
  }
  x::sustain = x::extent;
}

vec<Root> root(const Spectrum<double>& spectrum, Factor factor, vec<ADSR>& adsr) {
  window              = spectrum;
  x::start            = x::axis.pos;
  const double* rms   = window.ch.front();
  const size_t  arg   = algorithm::argmax(std::span<const double>(rms, window.qx.n));
  double        tot   = 0;
  double        pwr   = 0;
  size_t        n     = 0;
  for (ms i = x::start; i < x::extent; ++i) {
    tot += y::pos[i];
    pwr += db::lin(y::pos[i]);
    ++n;
  }
  const double plain = tot / n - rms[arg];
  y::noise           = rms[arg] + algorithm::noise_floor(window, plain);

  vec<Root> out;
  while (x::start < x::extent) {
    fractional_ascent(factor);
    if (x::cross >= x::extent) break;

    ascent(factor);
    bool keep = descent(factor);
    if (keep) {
      absolute_attack(factor);
      sustain(factor);
      const double width = factor.release * (x::end - x::peak);
      keep               = x::sustain >= x::end - width && x::sustain <= x::end + width;
      x::end             = x::sustain;
    }
    if (keep) {
      const double least = hold(x::end - x::onset, y::peak, y::floor, factor.end, factor.decay);
      if (x::end - x::peak < least) x::end = x::peak + least;
      if (x::end > x::extent) x::end = x::extent;
      double nrg = 0;
      for (ms i = x::onset; i < x::end; ++i) nrg += db::lin(y::pos[i]);
      const Spectrum<double> span   = window.window(x::onset, x::end);
      const double           center = algorithm::temporal_centroid(span);
      const double           slope  = (y::sustain - y::peak) / (x::sustain - x::peak);

      out.push_back(Root{});
      adsr.push_back(ADSR{
        .x0      = x::onset,
        .attack  = {  x::peak - x::onset, &x::axis},
        .peak    = x::peak,
        .decay   = ms{x::sustain - x::peak, &x::axis},
        .sustain = y::sustain,
        .release = { x::end - x::sustain, &x::axis},
        .x1      = x::end,
        .energy  = nrg / pwr,
        .slope   = slope,
        .center  = {              center, &x::axis}
      });
    }

    const ms next = keep ? x::end : xy::rise_run.x.back();
    if (next > x::start)
      x::start = next;
    else
      ++x::start;
  }
  return out;
}
}
