#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "algorithm/dsp_spectrum.h"
#include "algorithm/sim_drift.h"

#include "atom/atom_constant.h"

namespace qstorm::atom {

struct Quantamization {
  double ms;
  double overlap;
  double link;
  size_t gap;
  double spike_z;
  double recurrence_floor;
  double drift_ms;
  double drift_tolerance;
};

std::vector<Atom> atom_track(Atom window, const std::vector<int>& track) {
  std::vector<Atom> serie;
  for (size_t i = 0; i < track.size(); ++i)
    if (track[i] != -1) serie.push_back(window.atom[i].atom[track[i]]);
  return serie;
}

std::vector<int> min_cost_match(const std::vector<std::vector<double>>& a) {
  const int n = static_cast<int>(a.size());
  const double kInf = 1e18;
  std::vector<double> u(n + 1, 0.0), v(n + 1, 0.0);
  std::vector<int> p(n + 1, 0), way(n + 1, 0);
  for (int i = 1; i <= n; ++i) {
    p[0] = i;
    int j0 = 0;
    std::vector<double> minv(n + 1, kInf);
    std::vector<bool> used(n + 1, false);
    do {
      used[j0] = true;
      const int i0 = p[j0];
      double delta = kInf;
      int j1 = -1;
      for (int j = 1; j <= n; ++j)
        if (!used[j]) {
          const double cur = a[i0 - 1][j - 1] - u[i0] - v[j];
          if (cur < minv[j]) {
            minv[j] = cur;
            way[j] = j0;
          }
          if (minv[j] < delta) {
            delta = minv[j];
            j1 = j;
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
      p[j0] = p[j1];
      j0 = j1;
    } while (j0);
  }
  std::vector<int> match(n, -1);
  for (int j = 1; j <= n; ++j) match[p[j] - 1] = j - 1;
  return match;
}

void open_track(std::vector<std::vector<int>>& track, size_t window_count, size_t i, int atom) {
  track.emplace_back(window_count, -1);
  track.back()[i] = atom;
}

std::vector<std::vector<int>> hungarian_link(Atom window, double threshold) {
  const size_t window_count = window.atom.size();
  std::vector<std::vector<int>> track;
  for (int i = 0; i < static_cast<int>(window.atom[0].atom.size()); ++i)
    open_track(track, window_count, 0, i);

  for (size_t i = 1; i < window_count; ++i) {
    std::vector<size_t> alive;
    for (size_t j = 0; j < track.size(); ++j)
      if (track[j][i - 1] != -1) alive.push_back(j);
    const int alive_count = static_cast<int>(alive.size());
    const int atom_count = static_cast<int>(window.atom[i].atom.size());
    const int dimension = std::max(alive_count, atom_count);
    std::vector<std::vector<double>> cost(dimension, std::vector<double>(dimension, 0.0));
    for (int j = 0; j < alive_count; ++j) {
      const int last = track[alive[j]][i - 1];
      for (int k = 0; k < atom_count; ++k) {
        const std::vector<double>& spectrum =
          std::get<Raw>(window.atom[i - 1].atom[last].component[0]).value;
        const std::vector<double>& spectrum_b =
          std::get<Raw>(window.atom[i].atom[k].component[0]).value;
        const double sim = algorithm::cosine(spectrum, spectrum_b);
        if (sim >= threshold) cost[j][k] = -sim;
      }
    }

    const std::vector<int> match = min_cost_match(cost);
    std::vector<bool> atom_taken(atom_count, false);
    for (int j = 0; j < alive_count; ++j) {
      const int k = match[j];
      if (k < atom_count && cost[j][k] < 0.0) {
        track[alive[j]][i] = k;
        atom_taken[k] = true;
      }
    }
    for (int k = 0; k < atom_count; ++k)
      if (!atom_taken[k]) open_track(track, window_count, i, k);
  }
  return track;
}

std::vector<std::vector<int>> greedy_link(Atom window, double threshold) {
  const size_t window_count = window.atom.size();
  std::vector<std::vector<int>> track;
  for (int i = 0; i < static_cast<int>(window.atom[0].atom.size()); ++i)
    open_track(track, window_count, 0, i);
  for (size_t i = 1; i < window_count; ++i) {
    const int atom_count = static_cast<int>(window.atom[i].atom.size());
    std::vector<std::vector<double>> cost(track.size(), std::vector<double>(atom_count, 0.0));
    for (size_t j = 0; j < track.size(); ++j) {
      const int last = track[j][i - 1];
      if (last == -1) continue;
      const std::vector<double>& spectrum =
        std::get<Raw>(window.atom[i - 1].atom[last].component[0]).value;
      for (int k = 0; k < atom_count; ++k) {
        const std::vector<double>& spectrum_b =
          std::get<Raw>(window.atom[i].atom[k].component[0]).value;
        const double sim = algorithm::cosine(spectrum, spectrum_b);
        if (sim >= threshold) cost[j][k] = -sim;
      }
    }

    std::vector<bool> atom_taken(atom_count, false);
    std::vector<bool> track_taken(track.size(), false);
    for (;;) {
      double best = 0.0;
      size_t best_t = 0;
      int best_b = -1;
      for (size_t j = 0; j < track.size(); ++j) {
        if (track_taken[j]) continue;
        for (int k = 0; k < atom_count; ++k) {
          if (atom_taken[k] || cost[j][k] >= best) continue;
          best = cost[j][k];
          best_t = j;
          best_b = k;
        }
      }
      if (best_b == -1) break;
      track_taken[best_t] = true;
      atom_taken[best_b] = true;
      track[best_t][i] = best_b;
    }
    for (int k = 0; k < atom_count; ++k)
      if (!atom_taken[k]) open_track(track, window_count, i, k);
  }
  return track;
}

std::vector<std::vector<int>> threshold_link(Atom window, double threshold, size_t gap) {
  const size_t window_count = window.atom.size();
  std::vector<std::vector<int>> track;
  std::vector<size_t> last_w;
  for (int i = 0; i < static_cast<int>(window.atom[0].atom.size()); ++i) {
    open_track(track, window_count, 0, i);
    last_w.push_back(0);
  }

  for (size_t i = 1; i < window_count; ++i) {
    const int atom_count = static_cast<int>(window.atom[i].atom.size());
    std::vector<std::vector<double>> cost(track.size(), std::vector<double>(atom_count, 0.0));
    for (size_t j = 0; j < track.size(); ++j) {
      if (i - last_w[j] > gap) continue;
      const int last = track[j][last_w[j]];
      const std::vector<double>& spectrum =
        std::get<Raw>(window.atom[last_w[j]].atom[last].component[0]).value;
      for (int k = 0; k < atom_count; ++k) {
        const std::vector<double>& spectrum_b =
          std::get<Raw>(window.atom[i].atom[k].component[0]).value;
        const double sim = algorithm::cosine(spectrum, spectrum_b);
        if (sim >= threshold) cost[j][k] = -sim;
      }
    }

    std::vector<bool> atom_taken(atom_count, false);
    std::vector<bool> track_taken(track.size(), false);
    for (;;) {
      double best = 0.0;
      size_t best_t = 0;
      int best_b = -1;
      for (size_t j = 0; j < track.size(); ++j) {
        if (track_taken[j]) continue;
        for (int k = 0; k < atom_count; ++k) {
          if (atom_taken[k] || cost[j][k] >= best) continue;
          best = cost[j][k];
          best_t = j;
          best_b = k;
        }
      }
      if (best_b == -1) break;
      track_taken[best_t] = true;
      atom_taken[best_b] = true;
      track[best_t][i] = best_b;
      last_w[best_t] = i;
    }
    for (int k = 0; k < atom_count; ++k)
      if (!atom_taken[k]) {
        open_track(track, window_count, i, k);
        last_w.push_back(i);
      }
  }
  return track;
}

double reconstruction_fidelity(const std::vector<std::vector<int>>& track, Atom window,
                               const algorithm::Spectrogram<float>& spectrogram, double hop_ms,
                               double window_ms) {
  const algorithm::Quantization<float>& x_ms = spectrogram.x_ms;
  const algorithm::Quantization<float>& z_hz = spectrogram.z_hz;
  const size_t    hop = static_cast<size_t>(std::round(hop_ms*x_ms.density));
  const size_t wframe = static_cast<size_t>(std::round(window_ms*x_ms.density));
  std::vector<double> reconstruct(z_hz.extent*x_ms.extent, 0.0);
  std::vector<double> cover(x_ms.extent, 0.0);
  for (const std::vector<int>& t : track) {
    std::vector<double> rep(z_hz.extent, 0.0);
    for (size_t i = 0; i < t.size(); ++i) {
      const int b = t[i];
      if (b == -1) continue;
      const std::vector<double>& spectrum =
        std::get<Raw>(window.atom[i].atom[b].component[0]).value;
      for (size_t j = 0; j < z_hz.extent; ++j) rep[j] += spectrum[j];
    }
    double norm = 0.0;
    for (double e : rep) norm += e*e;
    norm = std::sqrt(norm);
    for (double& e : rep) e /= norm;
    for (size_t i = 0; i < t.size(); ++i) {
      const int b = t[i];
      if (b == -1) continue;
      const std::vector<double>& activation =
        std::get<Raw>(window.atom[i].atom[b].component[1]).value;
      const size_t t0 = i*hop;
      for (size_t j = 0; j < activation.size(); ++j)
        for (size_t k = 0; k < z_hz.extent; ++k)
          reconstruct[k*x_ms.extent + (t0 + j)] += rep[k] * activation[j];
    }
  }

  for (size_t i = 0; i < window.atom.size(); ++i) {
    const size_t t0 = i*hop;
    for (size_t j = 0; j < wframe && t0 + j < x_ms.extent; ++j) cover[t0 + j] += 1.0;
  }

  double err = 0.0;
  double energy = 0.0;
  for (size_t i = 0; i < z_hz.extent; ++i)
    for (size_t j = 0; j < x_ms.extent; ++j) {
      const double model = cover[j] > 0.0 ? reconstruct[i*x_ms.extent +j] / cover[j] : 0.0;
      const double truth = spectrogram.channel.front()[j*z_hz.extent +i];
      const double     d = truth - model;
      err += d*d;
      energy += truth*truth;
    }
  return 1.0 - std::sqrt(err / energy);
}

std::vector<std::vector<int>> link_track(Atom window,
                                         const algorithm::Spectrogram<float>& spectrogram,
                                         double hop_ms, double window_ms, double threshold,
                                         size_t gap) {
  std::vector<std::vector<std::vector<int>>> candidate;
  candidate.push_back(greedy_link(window, threshold));
  candidate.push_back(hungarian_link(window, threshold));
  candidate.push_back(threshold_link(window, threshold, gap));
  std::vector<std::vector<std::vector<int>>> ranked = algorithm::rank(
      std::move(candidate), reconstruction_fidelity, window, spectrogram, hop_ms, window_ms);
  return std::move(ranked.front());
}

Atom assemble_atom(const std::vector<std::vector<int>>& best, Atom window,
                   const algorithm::Quantization<double>& x_ms,
                   const algorithm::Quantization<double>& parent_db, const Quantamization& spec) {
  Atom root;
  for (const std::vector<int>& track : best) {
    std::vector<Atom> series = atom_track(window, track);
    std::vector<double> center;
    std::vector<double> width;
    std::vector<double> energy;
    for (Atom& e : series) {
      Raw&                   spectrum = std::get<Raw>(e.component[0]);
      const Raw&    activation_window = std::get<Raw>(e.component[1]);
      const double spectrum_reference = *std::ranges::max_element(spectrum.value);
      std::vector<std::vector<double>> spectrum_channel(1);
      for (double v : spectrum.value)
        spectrum_channel[0].push_back(
            algorithm::db<double>::ratio(std::fabs(v) / spectrum_reference));

      const algorithm::Quantization<double> spectrum_time
      {
        1,
        activation_window.axis.density / static_cast<double>(activation_window.axis.extent),
        activation_window.axis.offset,
      };

      const algorithm::Spectrogram<double> qz =
          algorithm::as_spectrogram(spectrum_channel, spectrum_time, spectrum.axis,
                                    algorithm::y_db(spectrum_channel, parent_db));

      center.push_back(algorithm::spectral_centroid(qz));
      width.push_back(algorithm::spectral_spread(qz, center.back()));
      const std::vector<double>& activation = activation_window.value;
      double total = 0.0;
      for (double x : activation) total += x*x;
      energy.push_back(total);
    }

    size_t first = 0;
    while (track[first] == -1) ++first;
    algorithm::ms<double> start{x_ms.offset, &x_ms};
    for (size_t i = 0; i < first; ++i) ++start;
    const algorithm::Quantization<double> axis{center.size(), x_ms.density, start};
    const algorithm::Quantization<double>& frequency =
        std::get<Raw>(series.front().component[0]).axis;
    const algorithm::Quantization<double> full_band
    {
      1,
      frequency.density / static_cast<double>(frequency.extent),
      frequency.offset,
    };

    const double center_reference = *std::ranges::max_element(center);
    const double  width_reference = *std::ranges::max_element(width);
    const double energy_reference = *std::ranges::max_element(energy);
    std::vector<std::vector<double>> center_channel(1);
    std::vector<std::vector<double>> width_channel(1);
    std::vector<std::vector<double>> energy_channel(1);
    for (size_t i = 0; i < center.size(); ++i) {
      center_channel[0].push_back(
          algorithm::db<double>::ratio(std::fabs(center[i]) / center_reference));
      width_channel[0].push_back(
          algorithm::db<double>::ratio(std::fabs(width[i]) / width_reference));
      energy_channel[0].push_back(
          algorithm::db<double>::ratio(std::fabs(energy[i]) / energy_reference));
    }
    const algorithm::Spectrogram<double> center_env = algorithm::as_spectrogram(
        center_channel, axis, full_band, algorithm::y_db(center_channel, parent_db));
    const algorithm::Spectrogram<double> width_env = algorithm::as_spectrogram(
        width_channel, axis, full_band, algorithm::y_db(width_channel, parent_db));
    const algorithm::Spectrogram<double> energy_env = algorithm::as_spectrogram(
        energy_channel, axis, full_band, algorithm::y_db(energy_channel, parent_db));

    Atom child;
    child.component.push_back(algorithm::drift(center_env, spec.drift_ms, spec.drift_tolerance));
    child.component.push_back(algorithm::drift(width_env,  spec.drift_ms, spec.drift_tolerance));
    child.component.push_back(algorithm::drift(energy_env, spec.drift_ms, spec.drift_tolerance));
    root.atom.push_back(std::move(child));
  }
  return root;
}
}
