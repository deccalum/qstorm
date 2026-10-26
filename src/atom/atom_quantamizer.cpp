#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "algorithm/dsp_spectrum.h"
#include "algorithm/sim_drift.h"
#include "algorithm/sim_matrix.h"
#include "atom/atom_extractor.h"

namespace qstorm::atom {
struct Quantamization {
  double ms;
  double overlap;
  double link;
  size_t gap;
  double spike_z;
  double recurrence_floor;
  double drift;
  double tol;
};

vec<Atom> atom_track(Atom atom, const vec<int>& track) {
  vec<Atom> serie;
  for (size_t i = 0; i < track.size(); ++i)
    if (track[i] != -1) serie.push_back(atom.atom[i].atom[track[i]]);
  return serie;
}

void open_track(vvec<int>& track, size_t window_count, size_t i, int atom) {
  track.emplace_back(window_count, -1);
  track.back()[i] = atom;
}

vvec<int> hungarian_link(Atom atom, double threshold) {
  const size_t window_count = atom.atom.size();
  vvec<int>    track;
  for (int i = 0; i < static_cast<int>(atom.atom[0].atom.size()); ++i)
    open_track(track, window_count, 0, i);

  for (size_t i = 1; i < window_count; ++i) {
    vec<size_t> alive;
    for (size_t j = 0; j < track.size(); ++j)
      if (track[j][i - 1] != -1) alive.push_back(j);
    const int alive_count = static_cast<int>(alive.size());
    const int atom_count  = static_cast<int>(atom.atom[i].atom.size());
    const int dimension   = std::max(alive_count, atom_count);

    vvec cost(dimension, vec(dimension, 0.0));
    for (int j = 0; j < alive_count; ++j) {
      const int last = track[alive[j]][i - 1];
      for (int k = 0; k < atom_count; ++k) {
        const vec<>& z   = std::get<Raw>(atom.atom[i - 1].atom[last].component[0]).val;
        const vec<>& zy  = std::get<Raw>(atom.atom[i].atom[k].component[0]).val;
        const double sim = algorithm::cosine(z, zy);
        if (sim >= threshold) cost[j][k] = -sim;
      }
    }

    const vec<int> match = matrix::min_cost_match(cost);
    vec<bool>      atom_taken(atom_count, false);
    for (int j = 0; j < alive_count; ++j) {
      const int k = match[j];
      if (k < atom_count && cost[j][k] < 0.0) {
        track[alive[j]][i] = k;
        atom_taken[k]      = true;
      }
    }
    for (int k = 0; k < atom_count; ++k)
      if (!atom_taken[k]) open_track(track, window_count, i, k);
  }
  return track;
}

vvec<int> threshold_link(Atom atom, double threshold, size_t gap) {
  const size_t window_count = atom.atom.size();
  vvec<int>    track;
  vec<size_t>  last_w;
  for (int i = 0; i < static_cast<int>(atom.atom[0].atom.size()); ++i) {
    open_track(track, window_count, 0, i);
    last_w.push_back(0);
  }

  for (size_t i = 1; i < window_count; ++i) {
    const int atom_count = static_cast<int>(atom.atom[i].atom.size());
    vvec      cost(track.size(), vec(atom_count, 0.0));
    for (size_t j = 0; j < track.size(); ++j) {
      if (i - last_w[j] > gap) continue;
      const int    last     = track[j][last_w[j]];
      const vec<>& spectrum = std::get<Raw>(atom.atom[last_w[j]].atom[last].component[0]).val;
      for (int k = 0; k < atom_count; ++k) {
        const vec<>& spectrum_b = std::get<Raw>(atom.atom[i].atom[k].component[0]).val;
        const double sim        = algorithm::cosine(spectrum, spectrum_b);
        if (sim >= threshold) cost[j][k] = -sim;
      }
    }

    vec<bool> atom_taken(atom_count, false);
    vec<bool> track_taken(track.size(), false);
    for (;;) {
      double best   = 0.0;
      size_t best_t = 0;
      int    best_b = -1;
      for (size_t j = 0; j < track.size(); ++j) {
        if (track_taken[j]) continue;
        for (int k = 0; k < atom_count; ++k) {
          if (atom_taken[k] || cost[j][k] >= best) continue;
          best   = cost[j][k];
          best_t = j;
          best_b = k;
        }
      }
      if (best_b == -1) break;
      track_taken[best_t] = true;
      atom_taken[best_b]  = true;
      track[best_t][i]    = best_b;
      last_w[best_t]      = i;
    }
    for (int k = 0; k < atom_count; ++k)
      if (!atom_taken[k]) {
        open_track(track, window_count, i, k);
        last_w.push_back(i);
      }
  }
  return track;
}

double reconstruction_fidelity(const vvec<int>& track, Atom atom, const Spectrum<float>& spectrum,
                               double hop_ms, double window_ms) {
  const Qtz<float>& qx     = spectrum.qx;
  const Qtz<float>& qz     = spectrum.qz;
  const size_t      hop    = std::round(hop_ms * qx.rho);
  const size_t      wframe = std::round(window_ms * qx.rho);
  vec               restruct(qz.n * qx.n, 0.0);
  vec               cover(qx.n, 0.0);
  for (const vec<int>& t : track) {
    vec rep(qz.n, 0.0);
    for (size_t i = 0; i < t.size(); ++i) {
      const int b = t[i];
      if (b == -1) continue;
      const vec<>& spectrum = std::get<Raw>(atom.atom[i].atom[b].component[0]).val;
      for (size_t j = 0; j < qz.n; ++j) rep[j] += spectrum[j];
    }
    double norm = 0.0;
    for (double e : rep) norm += e * e;
    norm = std::sqrt(norm);
    for (double& e : rep) e /= norm;
    for (size_t i = 0; i < t.size(); ++i) {
      const int b = t[i];
      if (b == -1) continue;
      const vec<>& activation = std::get<Raw>(atom.atom[i].atom[b].component[1]).val;
      const size_t t0         = i * hop;
      for (size_t j = 0; j < activation.size(); ++j)
        for (size_t k = 0; k < qz.n; ++k) restruct[k * qx.n + (t0 + j)] += rep[k] * activation[j];
    }
  }

  for (size_t i = 0; i < atom.atom.size(); ++i) {
    const size_t t0 = i * hop;
    for (size_t j = 0; j < wframe && t0 + j < qx.n; ++j) cover[t0 + j] += 1.0;
  }

  double err    = 0.0;
  double energy = 0.0;
  for (size_t i = 0; i < qz.n; ++i)
    for (size_t j = 0; j < qx.n; ++j) {
      const double rebuilt = cover[j] > 0.0 ? restruct[i * qx.n + j] / cover[j] : 0.0;
      const double truth   = spectrum.ch.front()[j * qz.n + i];
      const double d       = truth - rebuilt;
      err += d * d;
      energy += truth * truth;
    }
  return 1.0 - std::sqrt(err / energy);
}

vvec<int> link_track(Atom atom, const Spectrum<float>& spectrum, double hop_ms, double window_ms,
                     double threshold, size_t gap) {
  vvvec<int> input;
  input.push_back(threshold_link(atom, threshold, 1));
  input.push_back(hungarian_link(atom, threshold));
  input.push_back(threshold_link(atom, threshold, gap));
  vvvec<int> ranked
    = algorithm::rank(std::move(input), reconstruction_fidelity, atom, spectrum, hop_ms, window_ms);
  return std::move(ranked.front());
}

Atom assemble_atom(const vvec<int>& best, Atom window, const Qtz<double>& qx,
                   const Qtz<double>& parent_db, const Quantamization& spec) {
  Atom assembly;
  for (const vec<int>& track : best) {
    vec<Atom> series = atom_track(window, track);
    vec       center;
    vec       width;
    vec       energy;
    for (Atom& e : series) {
      Raw&         z       = std::get<Raw>(e.component[0]);
      const Raw&   act_win = std::get<Raw>(e.component[1]);
      const double z_ref   = *std::ranges::max_element(z.val);
      vvec         ch(1);
      for (double v : z.val) ch[0].push_back(algorithm::db<double>::ratio(std::fabs(v) / z_ref));

      const Qtz<double>      qx{1, act_win.axis.rho / act_win.axis.n, act_win.axis.pos};
      const Qtz<double>      qy = algorithm::ydb(ch, parent_db);
      const Spectrum<double> qz = algorithm::as_spectrum(ch, qx, qy, z.axis);

      const zvec<double> band{qz};
      center.push_back(algorithm::spectral_centroid(band));
      width.push_back(algorithm::spectral_spread(band, center.back()));
      const vec<>& activation = act_win.val;
      double       tot        = 0.0;
      for (double x : activation) tot += x * x;
      energy.push_back(tot);
    }

    size_t i0 = 0;
    while (track[i0] == -1) ++i0;
    ms<double> x0{qx.pos, &qx};
    for (size_t i = 0; i < i0; ++i) ++x0;
    const Qtz<double>  axis{center.size(), qx.rho, x0};
    const Qtz<double>& freq = std::get<Raw>(series.front().component[0]).axis;
    const Qtz<double>  full_band{1, freq.rho / freq.n, freq.pos};

    const double center_ref = *std::ranges::max_element(center);
    const double width_ref  = *std::ranges::max_element(width);
    const double energy_ref = *std::ranges::max_element(energy);
    vvec         center_ch(1);
    vvec         width_ch(1);
    vvec         energy_ch(1);
    for (size_t i = 0; i < center.size(); ++i) {
      center_ch[0].push_back(algorithm::db<double>::ratio(std::fabs(center[i]) / center_ref));
      width_ch[0].push_back(algorithm::db<double>::ratio(std::fabs(width[i]) / width_ref));
      energy_ch[0].push_back(algorithm::db<double>::ratio(std::fabs(energy[i]) / energy_ref));
    }

    Qtz cy = algorithm::ydb(center_ch, parent_db);
    Qtz wy = algorithm::ydb(width_ch, parent_db);
    Qtz ey = algorithm::ydb(energy_ch, parent_db);

    const Spectrum<double> center_env = algorithm::as_spectrum(center_ch, axis, cy, full_band);
    const Spectrum<double> width_env  = algorithm::as_spectrum(width_ch, axis, wy, full_band);
    const Spectrum<double> energy_env = algorithm::as_spectrum(energy_ch, axis, ey, full_band);

    Atom child;
    child.component.push_back(algorithm::drift(center_env, spec.drift, spec.tol));
    child.component.push_back(algorithm::drift(width_env, spec.drift, spec.tol));
    child.component.push_back(algorithm::drift(energy_env, spec.drift, spec.tol));
    assembly.atom.push_back(std::move(child));
  }
  return assembly;
}
}
