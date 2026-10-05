#include <cstddef>
#include <utility>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "algorithm/sim_matrix.h"
#include "atom/atom_constant.h"

namespace qstorm::atom {

template <typename R, typename K, typename C>
Atom through_nmf(const algorithm::Matrix<R, K>& spectrum,
                 const algorithm::Matrix<K, C>& activation) {
  const size_t   row = spectrum.row.extent;
  const size_t  rank = spectrum.col.extent;
  const size_t frame = activation.col.extent;
  const size_t row_density = static_cast<double>(spectrum.row.density);
  const size_t row_offset  = static_cast<double>(spectrum.row.offset);
  const size_t col_density = static_cast<double>(activation.col.density);
  const size_t col_offset  = static_cast<double>(activation.col.offset);
  const algorithm::Quantization<double> freq{row,   row_density, row_offset};
  const algorithm::Quantization<double> time{frame, col_density, col_offset};
  Atom atom;

  for (size_t i = 0; i < rank; ++i) {
    Raw spectrum_raw;
    spectrum_raw.value.resize(row);
    for (size_t j = 0; j < row; ++j) spectrum_raw.value[j] = spectrum.cell[j*rank +i];
    spectrum_raw.axis = freq;

    Raw activation_raw;
    activation_raw.value.resize(frame);
    for (size_t j = 0; j < frame; ++j) activation_raw.value[j] = activation.cell[i*frame +j];
    activation_raw.axis = time;

    Atom component;
    component.component.push_back(std::move(spectrum_raw));
    component.component.push_back(std::move(activation_raw));
    atom.atom.push_back(std::move(component));
  }
  return atom;
}

}
