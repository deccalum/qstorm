#pragma once
#include <cstddef>
#include <optional>
#include <variant>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "algorithm/sim_drift.h"
#include "algorithm/sim_waveform.h"
#include "atom/atom_thunder.h"

namespace qstorm::atom {

struct Domain {
  using Factor = thunder::Factor;
  using Root = thunder::Root;
};

struct Spectral {
  algorithm::Drift center;
  algorithm::Drift spread;
  algorithm::Drift shape;
  algorithm::Drift flatness;
};

struct ADSR {
  size_t time;
  double attack_ms;
  size_t peak;
  std::optional<double> decay_ms;
  std::optional<double> sustain;
  double release_ms;
  size_t end;
  int energy;
};

struct Transient {
  double prominence_db;
  double ms;
  ADSR adsr;
  Spectral spectral;
  algorithm::Locus locus;
};

template <typename T>
struct Series {
  std::vector<T> step;
  algorithm::Locus locus;
  algorithm::Drift envelope;
};

struct Raw {
  std::vector<double> value;
  algorithm::Quantization<double> axis;
};

using Component = std::variant<
Raw, Spectral, algorithm::Drift, ADSR, algorithm::Locus, Series<Spectral>, Series<algorithm::Drift>
>;

struct Atom {
  std::optional<Domain> domain;
  std::vector<Component> component;
  std::vector<Atom> atom;
};

template <typename D>
auto extract_root_atom(const algorithm::Spectrogram<double>& window, const D& dna) {
  thunder::window = window;
  return root(dna);
}

}
