#include <chrono>
#include <cstddef>
#include <cstdint>
#include <meta>
#include <qr/qr.hpp>
#include <span>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#include "atom_version.h"
#include "io/bridge.h"
#include "io/database.h"

namespace qstorm::io::table {

consteval std::meta::info lookup(std::meta::info name) {
  for (std::meta::info each : qr::detail::members_of<^^setting>())
    if (std::meta::is_type_alias(each) && std::meta::identifier_of(each) == std::meta::identifier_of(name))
      return each;
  return std::meta::info{};
}

template <typename Shape>
void extents(std::vector<std::size_t>& bounds) {
  if constexpr (std::rank_v<Shape> > 0) {
    bounds.push_back(std::extent_v<Shape>);
    extents<std::remove_extent_t<Shape>>(bounds);
  }
}

template <typename Shape>
void split(std::vector<std::vector<std::size_t>>& columns) {
  std::vector<std::size_t> whole;
  extents<Shape>(whole);
  for (std::size_t bound : whole) columns.push_back({bound});
}

template <typename M, std::meta::info found>
void collect_member(std::vector<std::vector<std::size_t>>& columns) {
  const std::vector<std::size_t> plain{1};
  if constexpr (qr::is_column_value<M>()) {
    if constexpr (found == std::meta::info{}) {
      columns.push_back(plain);
    } else {
      std::vector<std::size_t> bounds;
      extents<typename [:found:]>(bounds);
      columns.push_back(bounds);
    }
  } else if constexpr (found == std::meta::info{}) {
    columns.insert(columns.end(), qr::detail::leaf_count<M>(), plain);
  } else {
    static_assert(std::rank_v<typename [:found:]> == qr::detail::leaf_count<M>());
    split<typename [:found:]>(columns);
  }
}

template <typename T>
void collect(std::vector<std::vector<std::size_t>>& columns) {
  template for (constexpr std::meta::info base : qr::detail::bases_of<T>()) {
    collect<typename [:std::meta::type_of(base):]>(columns);
  }
  template for (constexpr std::meta::info member : qr::detail::nonstatic_members_of<T>()) {
    collect_member<typename [:std::meta::type_of(member):], lookup(member)>(columns);
  }
}

template <typename T>
std::vector<std::vector<std::size_t>> bounds_of() {
  std::vector<std::vector<std::size_t>> columns;
  constexpr std::meta::info whole = lookup(^^T);
  if constexpr (whole == std::meta::info{}) {
    collect<T>(columns);
  } else {
    static_assert(std::rank_v<typename [:whole:]> == qr::detail::leaf_count<T>());
    split<typename [:whole:]>(columns);
  }
  return columns;
}

template <typename T>
std::vector<column_setting> settings_of() {
  const std::vector<qr::ColumnSpec>& columns = qr::columns_of<T>();
  const std::vector<std::vector<std::size_t>> bounds = bounds_of<T>();
  std::vector<column_setting> rows;
  for (std::size_t i = 0; i < columns.size(); ++i) {
    column_setting row{std::string(qr::table_name_of<T>()), columns[i].name, {}, {}};
    for (std::size_t bound : bounds[i]) {
      const char each = setting::token[bound / 4];
      row.code.push_back(static_cast<double>(bound % 4));
      if (each) row.token += each;
    }
    rows.push_back(row);
  }
  return rows;
}

void create() {
  db.emplace(env::value.database_url);
  db->create<source_info>();
  db->create<window>();
  db->create<atom_version>();
  db->create<dna>();
  db->create<run>();
  db->create<root>();
  db->create<column_setting>();
  db->insert_missing(settings_of<source_info>());
  db->insert_missing(settings_of<window>());
  db->insert_missing(settings_of<atom_version>());
  db->insert_missing(settings_of<dna>());
  db->insert_missing(settings_of<run>());
  db->insert_missing(settings_of<root>());
}

std::tuple<int64_t, std::vector<std::byte>> source(int64_t id) {
  const std::vector<source_info> rows = db->select<source_info>(qr::column<source_info, &source_info::id> == id);
  const source_info& row = rows[0];
  return std::tuple{row.id, db->read_blob(row.audio)};
}

void store(source_info row, std::span<const std::byte> audio) {
  const std::vector<source_info> stored = db->select<source_info>(qr::column<source_info, &source_info::id> == row.id);
  for (const source_info& each : stored) {
    row.start_ms = each.start_ms;
    row.end_ms = each.end_ms;
    row.windows = each.windows;
    row.atoms = each.atoms;
  }
  row.audio = db->write_blob(audio);
  db->upsert(row);
  for (const source_info& each : stored) db->remove_blob(each.audio);
}

std::vector<window> windows() { return db->select<window>(); }

template <>
std::vector<root> select<root>(int64_t run_id) {
  std::vector<root> found;
  for (const root& each : db->select<root>(qr::column<root, &root::run> == run_id)) {
    if (!found.empty() && found.back().index == each.index) found.back() = each;
    else found.push_back(each);
  }
  return found;
}

atom::Domain::Factor factor(int64_t dna_id) {
  const std::vector<dna> rows = db->select<dna>(qr::column<dna, &dna::id> == dna_id);
  return rows[0].factor;
}

int64_t version_id() {
  const std::vector<atom_version> stored = db->select<atom_version>(qr::column<atom_version, &atom_version::header_blob> == std::string(qstorm::version::header_blob) && qr::column<atom_version, &atom_version::source_blob> == std::string(qstorm::version::source_blob));
  for (const atom_version& each : stored) return each.id;
  static constexpr char header_text[] = {
#embed "../atom/atom_thunder.h" suffix(, 0)
  };
  static constexpr char source_text[] = {
#embed "../atom/atom_thunder.cpp" suffix(, 0)
  };
  db->insert(atom_version{.id = 0, .commit = qstorm::version::commit, .header_blob = qstorm::version::header_blob, .source_blob = qstorm::version::source_blob, .header_text = header_text, .source_text = source_text});
  return version_id();
}

int64_t record(int64_t source_id, const atom::Domain::Factor& factor, const algorithm::Spectrogram<double>& window, std::span<const atom::Domain::Root> found) {
  db->upsert(run{.source_id = source_id, .time_point = std::chrono::system_clock::now(), .factor = factor, .start = window.begin(), .end = window.end(), .version = version_id()});
  db->remove(qr::delete_from<root>().where(qr::column<root, &root::run> == source_id));
  std::vector<root> rows;
  for (size_t i = 0; i < found.size(); ++i) rows.push_back(root{found[i], source_id, static_cast<int64_t>(i), false, found[i].end - found[i].onset.x});
  db->insert(rows);
  source_info source = db->select<source_info>(qr::column<source_info, &source_info::id> == source_id)[0];
  source.windows = static_cast<int64_t>(found.size());
  db->upsert(source);
  return source_id;
}

struct column_editor {
  const std::vector<qr::ColumnSpec>& columns;
  const char* column;
  double value;
  std::size_t i = 0;

  template <typename V>
  void cell(V& cell) {
    if constexpr (std::is_assignable_v<V&, double>)
      if (columns[i].name == column) cell = value;
    ++i;
  }
};

void edit(int64_t run_id, int64_t index, const char* column, double value) {
  const std::vector<root> rows = db->select<root>(qr::column<root, &root::run> == run_id && qr::column<root, &root::index> == index);
  root after = rows.back();
  const std::vector<qr::ColumnSpec>& columns = qr::columns_of<root>();
  column_editor visit{columns, column, value};
  qr::each_value(after, visit);
  after.length = after.end - after.onset.x;
  after.edit = true;
  db->upsert(after);
}

std::size_t clear(atom::Domain) {
  const std::size_t roots = db->remove(qr::delete_from<root>());
  const std::size_t runs = db->remove(qr::delete_from<run>());
  return roots + runs;
}

}
