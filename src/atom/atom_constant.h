#pragma once
#include <optional>
#include <variant>

#include "algorithm/axis_quantizer.h"
#include "algorithm/sim_drift.h"
#include "algorithm/sim_waveform.h"

namespace qstorm::atom {

struct Spectral {
  double peak;      // highest power density
  double center;    // power-weighted centroid
  double spread;    // power-weighted standard deviation around center
  double low;       // lowest bin at the spectrum's level
  double high;      // highest bin at the spectrum's level
  double tilt;      // power density dB per octave
  double flatness;  // geometric over arithmetic mean of the power density
};

struct ADSR {
  ms<double>                x0;
  ms<double>                attack;
  ms<double>                peak;
  std::optional<ms<double>> decay;
  std::optional<double>     sustain;
  ms<double>                release;
  ms<double>                x1;
  double                    energy;
  double                    slope;
  ms<double>                center;
};

struct Transient {
  double           prominence;
  double           ms;
  ADSR             adsr;
  Spectral         spectral;
  algorithm::Locus locus;
};

template <typename T>
struct Series {
  vec<T>           step;
  algorithm::Locus locus;
  algorithm::Drift envelope;
};

struct Raw {
  vec<>       val;
  Qtz<double> axis;
};

using Component = std::variant<Raw, Spectral, algorithm::Drift, ADSR, algorithm::Locus,
                               Series<Spectral>, Series<algorithm::Drift>>;
}
