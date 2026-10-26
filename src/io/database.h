#pragma once

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <optional>
#include <qr/qr.hpp>
#include <span>
#include <string>
#include <type_traits>
#include <utility>

#include "algorithm/axis_quantizer.h"
#include "algorithm/sim_matrix.h"
#include "atom/atom_extractor.h"
#include "io/interpret_file.h"

namespace qstorm::table {
inline std::optional<qr::Connection> db;

template <typename T>
struct blob {
  struct held {
    vec<std::byte> bytes;
    T              value;

    explicit held(vec<std::byte> read)
      : bytes(std::move(read)), value([&] {
        if constexpr (std::is_same_v<T, std::span<const std::byte>>)
          return std::span<const std::byte>(bytes);
        else
          return io::from_bytes<T>(std::span<const std::byte>(bytes));
      }()) {}
  };

  qr::oid id;

  mutable std::shared_ptr<const held> loaded;

  const T& operator*() const {
    if (!loaded) loaded = std::make_shared<const held>(db->read_blob(id));
    return loaded->value;
  }
  const T* operator->() const { return &**this; }
};

template <typename T>
qr::oid to_column(const blob<T>& column) {
  return column.id;
}

template <typename T>
blob<T> from_column(qr::oid id, std::type_identity<blob<T>>) {
  return {id, {}};
}

template <typename X, typename Y = X, typename Z = X>
struct spectrum {
  Qtz<X>                  qx;
  Qtz<Y>                  qy;
  Qtz<Z>                  qz;
  int64_t                 count;
  blob<Spectrum<X, Y, Z>> data;

  const Spectrum<X, Y, Z>& operator*() const { return *data; }
  const Spectrum<X, Y, Z>* operator->() const { return data.operator->(); }
};

struct source {
  QR_KEY int64_t                   id;
  std::optional<std::string>       domain;
  std::string                      origin;
  double                           hz;
  int                              channels;
  int                              bit_depth;
  blob<std::span<const std::byte>> audio;
  std::optional<std::string>       title;
  std::optional<std::string>       uploader;
  std::optional<std::string>       url;
  std::optional<io::Original>      original;
};

struct decoded {
  QR_KEY           QR_REFERENCES(table::source) int64_t source_id;
  spectrum<double> sample;
  spectrum<double> log_mono;
  spectrum<double> linear_mono;
  blob<vec<Level>> waveform;
};

struct atom_version {
  QR_KEY QR_SERIAL int64_t id;
  std::string              commit;
  std::string              header_blob;
  std::string              source_blob;
  std::string              header_text;
  std::string              source_text;
};

struct dna {
  QR_KEY QR_SERIAL int64_t   id;
  std::optional<std::string> name;
  atom::Factor               factor;
  QR_REFERENCES(table::atom_version) std::optional<int64_t> version;
};

struct ard {
  QR_KEY QR_SERIAL int64_t   id;
  std::optional<std::string> name;
  matrix::Ard                factor;
};

struct extraction {
  QR_KEY QR_SERIAL int64_t id;
  QR_REFERENCES(table::decoded) QR_CASCADE int64_t source_id;
  std::chrono::system_clock::time_point time_point;
  atom::Factor                          factor;
  ms<double>                            x0;
  ms<double>                            x1;
  QR_REFERENCES(table::atom_version) std::optional<int64_t> version;
};

struct root {
  QR_KEY                           QR_REFERENCES(table::extraction) QR_CASCADE int64_t extraction;
  QR_KEY int64_t                   index;
  QR_KEY int64_t                   revision;
  atom::Domain                     domain;
  atom::ADSR                       adsr;
  blob<std::span<const std::byte>> audio;
  bool                             rejected;
};

struct stft {
  QR_KEY             QR_REFERENCES(table::extraction) QR_CASCADE int64_t extraction;
  QR_KEY int64_t     index;
  QR_KEY int64_t     revision;
  blob<zvec<double>> stft;
};

struct atomization {
  QR_KEY QR_SERIAL int64_t id;
  QR_REFERENCES(table::extraction) QR_CASCADE int64_t extraction;
  int64_t                               index;
  int64_t                               revision;
  std::chrono::system_clock::time_point time_point;
  int                                   channel;
  matrix::Ard                           factor;
  double                                floor;
  double                                divergence;
  int                                   iteration;
  vec<>                                 trace;
  vec<>                                 residual;
};

struct component {
  QR_KEY                           QR_REFERENCES(table::atomization) QR_CASCADE int64_t atomization;
  QR_KEY int64_t                   index;
  double                           lambda;
  double                           divergence;
  bool                             relevant;
  std::optional<atom::Spectral>    spectral;
  std::optional<atom::ADSR>        adsr;
  vec<>                            spectrum;
  vec<>                            activation;
  Qtz<double>                      octave;
  vec<>                            density;
  blob<std::span<const std::byte>> audio;
};

struct residual {
  QR_KEY             QR_REFERENCES(table::atomization) QR_CASCADE int64_t atomization;
  blob<zvec<double>> layer;
};

struct grouping {
  QR_KEY QR_SERIAL int64_t id;
  QR_REFERENCES(table::atomization) QR_CASCADE int64_t atomization;
  std::chrono::system_clock::time_point time_point;
};

struct component_pair {
  QR_KEY         QR_REFERENCES(table::grouping) QR_CASCADE int64_t grouping;
  QR_KEY int64_t component;
  QR_KEY int64_t other;
  double         cosine;
  double         correlation;
  bool           related;
};

struct component_group {
  QR_KEY         QR_REFERENCES(table::grouping) QR_CASCADE int64_t grouping;
  QR_KEY int64_t component;
  int64_t        spectrum;
  int64_t        event;
};

struct linking {
  QR_KEY QR_SERIAL int64_t              id;
  std::chrono::system_clock::time_point time_point;
  std::string                           domain;
  vec<int64_t>                          source;
  vec<int64_t>                          atomization;
};

struct component_link {
  QR_KEY         QR_REFERENCES(table::linking) QR_CASCADE int64_t linking;
  QR_KEY         QR_REFERENCES(table::atomization) QR_CASCADE int64_t atomization;
  QR_KEY int64_t component;
  int64_t        set_atomization;
  int64_t        set_component;
};

struct column_setting {
  QR_KEY std::string table;
  QR_KEY std::string column;
  vec<>              code;
  std::string        token;
};

void create();

template <typename T>
vec<T> select();
template <typename T>
vec<T> select(int64_t parent);
template <typename T, typename... Key>
T find(Key... key);
template <typename T, typename... Key>
T latest(Key... key);
template <typename Target, typename T>
Target follow(const T& row);
source audio_source(qr::oid audio);
void   remove_unreferenced_blob();
template <typename T>
void remove(int64_t parent);
void read(const char* table, const char* parent);
void clear_cache(const char* table, int64_t parent);

template <typename T>
void write(const T& row) {
  db->upsert(row);
}
template <typename T>
void write(const vec<T>& rows) {
  db->insert(rows);
}

template <typename T>
int64_t write_serial(const T& row) {
  return db->insert_returning(row);
}

template <typename T>
blob<T> write_blob(const T& value) {
  if constexpr (std::is_same_v<T, std::span<const std::byte>>) {
    const qr::oid id = db->write_blob(value);
    return {id, {}};
  } else {
    const vec<std::byte> bytes = io::to_bytes(value);
    const qr::oid        id    = db->write_blob(bytes);
    return {id, {}};
  }
}

template <typename X, typename Y, typename Z>
spectrum<X, Y, Z> write_spectrum(const Spectrum<X, Y, Z>& value) {
  const int64_t                 count = static_cast<int64_t>(value.ch.size());
  const blob<Spectrum<X, Y, Z>> data  = write_blob(value);
  return {value.qx, value.qy, value.qz, count, data};
}

void write_audio(int64_t id, const std::string& domain, const std::string& origin,
                 const io::Probe& probed, int bit_depth, const std::string& stem,
                 std::span<const std::byte> audio);
void write_decoded(int64_t source_id, const Spectrum<double>& sample,
                   const Spectrum<double>& log_mono, const Spectrum<double>& linear_mono,
                   const vec<Level>& waveform);
void write_root_extraction(const decoded& source, const dna& factor, const Spectrum<double>& window,
                           std::span<const atom::Atom>     found,
                           std::span<const vec<std::byte>> audio);
void write_atomization(const root& row, const ard& factor, int channel,
                       const matrix::Relevance& relevance, const vec<>& residual,
                       const atom::Atom& atom, std::span<const vec<std::byte>> audio);
void write_stft(const root& row, const zvec<double>& stft);
void write_residual(const atomization& source, const zvec<double>& layer);
void write_group(const atomization& source, const vec<component>& live, const vvec<>& cosine,
                 const vvec<>& correlation, const vvec<bool>& related, const vec<size_t>& spectrum,
                 const vec<size_t>& event);
vec<component>            select_live_component(int64_t atomization_id);
std::pair<vvec<>, vvec<>> split_spectrum_activation(const vec<component>& part);
int64_t                   thunder_version_id();
root                      edit_root(int64_t id, int64_t index, const char* column, double value);
void                      restore_root(int64_t id, int64_t index);
void                      edit_dna(int64_t id, const char* column, double value);
int64_t                   copy_dna(int64_t id);
size_t                    clear_root_extraction();
vec<int64_t>              comma_number(const char* text);
linking          link_scope(const char* domain, const char* source, const char* atomization);
vec<atomization> select_link_window(const linking& scope);
void             write_link(const linking& scope, const vec<std::pair<int64_t, int64_t>>& member,
                            const vec<size_t>& set);
}

namespace qstorm::setting {
struct part {};
inline constexpr char   token[5]{0, 'd', 'm', 'H', 't'};
inline constexpr size_t D = 5;
inline constexpr size_t M = 9;
inline constexpr size_t T = 17;
inline constexpr size_t H[16]{12, 13, 14, 15, 0, 0, 0, 0, 0, 0, 0, 0, 12, 13, 14, 15};

using source       = part[1][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2][2];
using root         = part[2][2][2][0][2][2][2][2][2][2][1][2][2][1][0][0][0][0][0][1][2][2][2][2][0];
using Factor       = part[0][0][2][2][2][2][2][2][2][0][2][2][2][2];
using atom_version = part[0][1][1][1][0][0];
using dna          = part[0][0][0][0][1][1][1][1][1][1][1][0][1][1][1][1][1];
using time_point   = part[D][D][M][M][2 [H][H]];
}
