#pragma once

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <meta>
#include <optional>
#include <qr/qr.hpp>
#include <span>
#include <string>
#include <tuple>
#include <vector>

#include "atom/atom_constant.h"

namespace qstorm::io::table {

inline std::optional<qr::Connection> db;

struct source_info {
  QR_KEY int64_t id;
  std::optional<std::string> domain;
  std::string origin;
  std::optional<double> start_ms;
  std::optional<double> end_ms;
  int hz;
  int channels;
  qr::oid audio;
  std::optional<std::string> title;
  std::optional<std::string> recorded_by;
  std::optional<std::string> url;
  std::optional<int64_t> windows;
  std::optional<int64_t> atoms;
};

struct window {
  QR_KEY QR_REFERENCES(table::source_info) int64_t id;
  QR_KEY double start;
  double end_ms;
  std::optional<std::vector<double>> inner;
};

struct atom_version {
  QR_KEY QR_SERIAL int64_t id;
  std::string commit;
  std::string header_blob;
  std::string source_blob;
  std::string header_text;
  std::string source_text;
};

struct dna {
  QR_KEY QR_SERIAL int64_t id;
  std::optional<std::string> name;
  atom::Domain::Factor factor;
  QR_REFERENCES(table::atom_version) std::optional<int64_t> version;
};

struct run {
  QR_KEY QR_REFERENCES(table::source_info) int64_t source_id;
  std::chrono::system_clock::time_point time_point;
  atom::Domain::Factor factor;
  double start;
  double end;
  QR_REFERENCES(table::atom_version) int64_t version;
};

struct root : atom::Domain::Root {
  QR_KEY QR_REFERENCES(table::run) int64_t run;
  QR_KEY int64_t index;
  QR_KEY bool edit;
  double length;
  bool rejected;
};

struct column_setting {
  QR_KEY std::string table;
  QR_KEY std::string column;
  std::vector<double> code;
  std::string token;
};

void create();

template <typename T>
std::vector<T> select() {
  return db->select<T>();
}

template <typename T>
std::vector<T> select(int64_t parent) {
  std::vector<T> found;
  template for (constexpr std::meta::info member : qr::detail::nonstatic_members_of<T>()) {
    if constexpr (qr::detail::referenced(member) != std::meta::info{})
      found = db->select<T>(qr::column<T, &[:member:]> == parent);
  }
  return found;
}

template <>
std::vector<root> select<root>(int64_t run_id);

std::tuple<int64_t, std::vector<std::byte>> source(int64_t id);
void store(source_info row, std::span<const std::byte> audio);
std::vector<std::tuple<int64_t, std::vector<std::byte>>> sources(atom::Domain domain);
std::vector<window> windows();
atom::Domain::Factor factor(int64_t dna_id);
int64_t version_id();

int64_t record(int64_t source_id, const atom::Domain::Factor& factor, const algorithm::Spectrogram<double>& window, std::span<const atom::Domain::Root> found);
void edit(int64_t run_id, int64_t index, const char* column, double value);

std::size_t clear(atom::Domain domain);

}

namespace setting {
struct part {};
inline constexpr char token[5]{0, 'd', 'm', 'H', 't'};
inline constexpr std::size_t D = 5;
inline constexpr std::size_t M = 9;
inline constexpr std::size_t T = 17;
inline constexpr std::size_t H[16]{12, 13, 14, 15, 0, 0, 0, 0, 0, 0, 0, 0, 12, 13, 14, 15};

using source_info = part[1][2][2][0][0][2][2][2][2][2][2][1][1];
using        root = part[0][1][2][2][2][2][2][2][1][2][2][2][2][2][T][0];
using      factor = part[0][0][2][2][2][2][2][2][2][0][2][2][2][2];
using atom_version = part[0][1][1][1][0][0];
using         dna = part[0][0][0][0][1][1][1][1][1][1][1][0][1][1][1][1][1];
using  time_point = part[D][D][M][M][2[H][H]];
}
