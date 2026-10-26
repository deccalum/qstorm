#include "io/database.h"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <meta>
#include <qr/qr.hpp>
#include <span>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "atom_version.h"
#include "io/bridge.h"

namespace qstorm::table {
template <typename T>
consteval std::span<const std::meta::info> keys_of() {
  vec<std::meta::info> found;
  for (std::meta::info member : qr::detail::nonstatic_members_of<T>())
    if (qr::detail::annotated(member)) found.push_back(member);
  return std::define_static_array(found);
}

template <typename T>
consteval std::meta::info parent_of() {
  for (std::meta::info member : qr::detail::nonstatic_members_of<T>())
    if (qr::detail::referenced(member) != std::meta::info{}) return member;
  return std::meta::info{};
}

consteval bool reaches(auto from, auto target) {
  constexpr std::meta::access_context context = std::meta::access_context::unchecked();
  for (std::meta::info member : std::meta::nonstatic_data_members_of(from, context)) {
    const std::meta::info next = qr::detail::referenced(member);
    if (next == target) return true;
    if (next != std::meta::info{} && reaches(next, target)) return true;
  }
  return false;
}

template <typename T, size_t... i, typename... Key>
qr::Predicate keyed(std::index_sequence<i...>, const Key&... key) {
  constexpr std::span<const std::meta::info> key_member = keys_of<T>();
  return (... && (qr::column<T, &[:key_member[i]:]> == key));
}

template <typename T>
bool same_entry(const T& a, const T& b) {
  bool same = true;
  template for (constexpr std::meta::info member : qr::detail::nonstatic_members_of<T>()) {
    if constexpr (qr::detail::annotated(member) && std::meta::identifier_of(member) != "revision")
      same = same && a.[:member:] == b.[:member:];
  }
  return same;
}

template <typename T>
vec<T> current(const vec<T>& rows) {
  vec<T> found;
  for (const T& each : rows) {
    if (!found.empty() && same_entry(found.back(), each))
      found.back() = each;
    else
      found.push_back(each);
  }
  return found;
}

template <typename T>
vec<T> select() {
  const vec<T> found = db->select<T>();
  return current(found);
}

template <typename T>
vec<T> select(int64_t parent) {
  constexpr std::meta::info member = parent_of<T>();
  const vec<T>              found  = db->select<T>(qr::column<T, &[:member:]> == parent);
  return current(found);
}

template <typename T, typename... Key>
T find(Key... key) {
  const vec<T> found = db->select<T>(keyed<T>(std::index_sequence_for<Key...>{}, key...));
  return found.front();
}

template <typename T, typename... Key>
T latest(Key... key) {
  const vec<T> found = db->select<T>(keyed<T>(std::index_sequence_for<Key...>{}, key...));
  return found.back();
}

template <typename Target, typename T>
Target follow(const T& row) {
  Target found;
  template for (constexpr std::meta::info member : qr::detail::nonstatic_members_of<T>()) {
    constexpr std::meta::info next = qr::detail::referenced(member);
    if constexpr (next == ^^Target)
      found = find<Target>(row.[:member:]);
    else if constexpr (next != std::meta::info{} && reaches(next, ^^Target))
      found = follow<Target>(find<typename[:next:]>(row.[:member:]));
  }
  return found;
}

template <typename T>
consteval std::meta::info audio_member() {
  for (std::meta::info member : qr::detail::nonstatic_members_of<T>())
    if (std::meta::type_of(member) == ^^blob<std::span<const std::byte>>) return member;
  return std::meta::info{};
}

source audio_source(qr::oid audio) {
  const blob<std::span<const std::byte>> held{audio, {}};
  source                                 found;
  template for (constexpr std::meta::info type : qr::detail::members_of<^^qstorm::table>()) {
    if constexpr (std::meta::is_type(type) && !std::meta::is_type_alias(type)
                  && std::meta::is_class_type(type)) {
      if constexpr (audio_member<typename[:type:]>() != std::meta::info{}) {
        constexpr std::meta::info   member = audio_member<typename[:type:]>();
        const vec<typename[:type:]> row
          = db->select<typename[:type:]>(qr::column<typename[:type:], &[:member:]> == held);
        for (const typename[:type:]& each : row) {
          if constexpr (type == ^^source)
            found = each;
          else
            found = follow<source>(each);
        }
      }
    }
  }
  return found;
}

template <typename T>
void remove(int64_t parent) {
  constexpr std::meta::info member = parent_of<T>();
  db->transaction([&] {
    db->remove(qr::delete_from<T>().where(qr::column<T, &[:member:]> == parent));
    remove_unreferenced_blob();
  });
}

template <typename T>
void read_rows(const char* parent) {
  if constexpr (parent_of<T>() == std::meta::info{}) {
    const vec<T> rows = select<T>();
    qr::write_rows(rows, stdout);
  } else {
    const vec<T> rows = parent ? select<T>(std::strtoll(parent, nullptr, 10)) : select<T>();
    qr::write_rows(rows, stdout);
  }
}

void read(const char* table, const char* parent) {
  template for (constexpr std::meta::info type : qr::detail::members_of<^^qstorm::table>()) {
    if constexpr (std::meta::is_type(type) && !std::meta::is_type_alias(type)
                  && std::meta::is_class_type(type))
      if (std::meta::identifier_of(type) == table) read_rows<typename[:type:]>(parent);
  }
}

void remove_unreferenced_blob() {
  vec<qr::ColumnRef> referencing;
  template for (constexpr std::meta::info type : qr::detail::members_of<^^qstorm::table>()) {
    if constexpr (std::meta::is_type(type) && !std::meta::is_type_alias(type)
                  && std::meta::is_class_type(type)) {
      for (const qr::ColumnSpec& column : qr::columns_of<typename[:type:]>()) {
        if (column.type != qr::sql_type<qr::oid>()) continue;
        referencing.push_back(qr::ColumnRef{qr::table_name_of<typename[:type:]>(), column.name});
      }
    }
  }
  db->run(qr::compile_unlink_unreferenced_blobs(referencing));
}

void clear_cache(const char* table, int64_t parent) {
  template for (constexpr std::meta::info type : qr::detail::members_of<^^qstorm::table>()) {
    if constexpr (std::meta::is_type(type) && !std::meta::is_type_alias(type)
                  && std::meta::is_class_type(type))
      if constexpr (parent_of<typename[:type:]>() != std::meta::info{})
        if (std::meta::identifier_of(type) == table) remove<typename[:type:]>(parent);
  }
}

consteval std::meta::info lookup(std::meta::info type) {
  const std::meta::info plain = std::meta::dealias(type);
  const std::meta::info named
    = std::meta::has_template_arguments(plain) ? std::meta::template_of(plain) : plain;
  if (!std::meta::has_identifier(named)) return std::meta::info{};
  for (std::meta::info each : qr::detail::members_of<^^setting>())
    if (std::meta::is_type_alias(each)
        && std::meta::identifier_of(each) == std::meta::identifier_of(named))
      return each;
  return std::meta::info{};
}

template <typename Shape>
void extents(vec<size_t>& bounds) {
  if constexpr (std::rank_v<Shape> > 0) {
    bounds.push_back(std::extent_v<Shape>);
    extents<std::remove_extent_t<Shape>>(bounds);
  }
}

template <typename Shape>
void split(vvec<size_t>& columns) {
  vec<size_t> whole;
  extents<Shape>(whole);
  for (size_t bound : whole) columns.push_back({bound});
}

template <typename M, std::meta::info found, size_t K = 0>
void collect_member(vvec<size_t>& columns) {
  const vec<size_t> plain{1};
  if constexpr (qr::is_column_value<M>()) {
    if constexpr (found == std::meta::info{})
      columns.push_back(plain);
    else {
      vec<size_t> bounds;
      extents<typename[:found:]>(bounds);
      columns.push_back(bounds);
    }
  } else if constexpr (qr::IsVariantType<M>::value) {
    using Alternative = std::variant_alternative_t<K, M>;
    collect_member<Alternative, lookup(^^Alternative)>(columns);
    if constexpr (K + 1 < std::variant_size_v<M>) collect_member<M, found, K + 1>(columns);
  } else if constexpr (found == std::meta::info{}) {
    columns.insert(columns.end(), qr::detail::member_leaf_count<M>(), plain);
  } else {
    static_assert(std::rank_v<typename[:found:]> == qr::detail::leaf_count<M>());
    split<typename[:found:]>(columns);
  }
}

template <typename T>
void collect(vvec<size_t>& columns) {
  template for (constexpr std::meta::info base : qr::detail::bases_of<T>()) {
    collect<typename[:std::meta::type_of(base):]>(columns);
  }
  template for (constexpr std::meta::info member : qr::detail::nonstatic_members_of<T>()) {
    collect_member<typename[:std::meta::type_of(member):], lookup(std::meta::type_of(member))>(
      columns);
  }
}

template <typename T>
vvec<size_t> bounds_of() {
  vvec<size_t>              columns;
  constexpr std::meta::info whole = lookup(^^T);
  if constexpr (whole == std::meta::info{}) {
    collect<T>(columns);
  } else {
    static_assert(std::rank_v<typename[:whole:]> == qr::detail::leaf_count<T>());
    split<typename[:whole:]>(columns);
  }
  return columns;
}

template <typename T>
vec<column_setting> settings_of() {
  const vec<qr::ColumnSpec>& columns = qr::columns_of<T>();
  const vvec<size_t>         bounds  = bounds_of<T>();
  vec<column_setting>        rows;
  for (size_t i = 0; i < columns.size(); ++i) {
    column_setting row{std::string(qr::table_name_of<T>()), columns[i].name, {}, {}};
    for (size_t bound : bounds[i]) {
      const char each = setting::token[bound / 4];
      row.code.push_back(static_cast<double>(bound % 4));
      if (each) row.token += each;
    }
    rows.push_back(row);
  }
  return rows;
}

template <typename T>
void migrate_table(const vec<std::string>& claimed) {
  const std::string          table(qr::table_name_of<T>());
  const vec<qr::ColumnSpec>& columns = qr::columns_of<T>();
  const vec<std::string> existing = db->run_rows<std::string>(qr::compile_select_stored_tables());
  if (std::ranges::contains(existing, table)) return;
  std::string best;
  size_t      best_shared = 0;
  size_t      best_count  = 0;
  for (const std::string& name : existing) {
    if (std::ranges::contains(claimed, name)) continue;
    const vec<qr::StoredColumn> stored
      = db->run_rows<qr::StoredColumn>(qr::compile_select_stored_columns(name));
    size_t shared = 0;
    for (const qr::StoredColumn& column : stored)
      shared += std::ranges::contains(columns, column.name, &qr::ColumnSpec::name);
    if (shared > best_shared) {
      best        = name;
      best_shared = shared;
      best_count  = 1;
    } else if (shared == best_shared)
      ++best_count;
  }
  if (best_shared == 0 || best_count != 1) return;
  db->run(qr::compile_rename_table<T>(best));
  db->update(qr::update<column_setting>()
               .set(qr::column<column_setting, &column_setting::table>, table)
               .where(qr::column<column_setting, &column_setting::table> == best));
}

template <typename T>
void migrate_columns() {
  const std::string          table(qr::table_name_of<T>());
  const vec<qr::ColumnSpec>& columns = qr::columns_of<T>();
  const vec<int64_t> wanted_type = db->run_rows<int64_t>(qr::compile_select_column_type_oids<T>());
  const vec<qr::StoredColumn> stored
    = db->run_rows<qr::StoredColumn>(qr::compile_select_stored_columns(table));
  vec<qr::StoredColumn> dropped;
  for (const qr::StoredColumn& column : stored)
    if (!std::ranges::contains(columns, column.name, &qr::ColumnSpec::name))
      dropped.push_back(column);
  vec<int64_t> added_type;
  for (size_t i = 0; i < columns.size(); ++i)
    if (!std::ranges::contains(stored, columns[i].name, &qr::StoredColumn::name))
      added_type.push_back(wanted_type[i]);
  for (size_t i = 0; i < columns.size(); ++i) {
    const qr::ColumnSpec& column  = columns[i];
    const int64_t         type    = wanted_type[i];
    const auto            kept    = std::ranges::find(stored, column.name, &qr::StoredColumn::name);
    const auto            renamed = std::ranges::find(dropped, type, &qr::StoredColumn::type);
    const bool       unique_type  = std::ranges::count(dropped, type, &qr::StoredColumn::type) == 1
                                 && std::ranges::count(added_type, type) == 1;
    qr::StoredColumn now{column.name, type, false};
    if (kept != stored.end()) {
      now = *kept;
    } else if (unique_type) {
      db->run(qr::compile_rename_column<T>(renamed->name, column.name));
      db->update(qr::update<column_setting>()
                   .set(qr::column<column_setting, &column_setting::column>, column.name)
                   .where(qr::column<column_setting, &column_setting::table> == table
                          && qr::column<column_setting, &column_setting::column> == renamed->name));
      now.not_null = renamed->not_null;
      dropped.erase(renamed);
    } else {
      db->run(qr::compile_add_column<T>(column));
    }
    if (now.type != type) db->run(qr::compile_retype_column<T>(column));
    if (now.not_null && column.nullable) db->run(qr::compile_drop_not_null<T>(column.name));
    if (!now.not_null && !column.nullable) {
      db->run(qr::compile_delete_null_rows<T>(column.name));
      db->run(qr::compile_set_not_null<T>(column.name));
    }
  }
  for (const qr::StoredColumn& column : dropped) {
    db->run(qr::compile_drop_column<T>(column.name));
    db->remove(qr::delete_from<column_setting>().where(
      qr::column<column_setting, &column_setting::table> == table
      && qr::column<column_setting, &column_setting::column> == column.name));
  }
  const vec<qr::StoredPrimaryKey> key
    = db->run_rows<qr::StoredPrimaryKey>(qr::compile_select_stored_primary_key(table));
  vec<std::string> stored_key;
  for (const qr::StoredPrimaryKey& each : key) stored_key.push_back(each.column);
  vec<std::string> wanted_key;
  for (const qr::ColumnSpec& column : columns)
    if (column.key) wanted_key.push_back(column.name);
  std::ranges::sort(stored_key);
  std::ranges::sort(wanted_key);
  if (stored_key == wanted_key) return;
  if (!key.empty()) db->run(qr::compile_drop_constraint<T>(key.front().constraint));
  db->run(qr::compile_add_primary_key<T>());
}

bool same_foreign_key(const qr::ColumnSpec& column, const qr::StoredForeignKey& key) {
  return column.name == key.column && column.references == key.references
      && column.cascade == key.cascade;
}

template <typename T>
void migrate_references() {
  const std::string               table(qr::table_name_of<T>());
  const vec<qr::ColumnSpec>&      columns = qr::columns_of<T>();
  const vec<qr::StoredForeignKey> stored
    = db->run_rows<qr::StoredForeignKey>(qr::compile_select_stored_foreign_keys(table));
  for (const qr::StoredForeignKey& key : stored) {
    bool kept = false;
    for (const qr::ColumnSpec& column : columns) kept = kept || same_foreign_key(column, key);
    if (!kept) db->run(qr::compile_drop_constraint<T>(key.constraint));
  }
  for (const qr::ColumnSpec& column : columns) {
    if (column.references.empty()) continue;
    bool present = false;
    for (const qr::StoredForeignKey& key : stored)
      present = present || same_foreign_key(column, key);
    if (!present) db->run(qr::compile_add_foreign_key<T>(column));
  }
}

template <typename... T>
void migrate() {
  const vec<std::string> claimed{std::string(qr::table_name_of<T>())...};
  (migrate_table<T>(claimed), ...);
  (db->create<T>(), ...);
  (migrate_columns<T>(), ...);
  (migrate_references<T>(), ...);
  (db->insert_missing(settings_of<T>()), ...);
}

void remove_stft_misaligned_with_root() {
  const vec<root>                                         roots = db->select<root>();
  std::map<std::tuple<int64_t, int64_t, int64_t>, double> x0;
  for (const root& each : roots)
    x0.insert_or_assign(std::tuple{each.extraction, each.index, each.revision}, each.adsr.x0);
  const vec<stft> written = db->select<stft>();
  for (const stft& each : written) {
    const std::tuple<int64_t, int64_t, int64_t> window{each.extraction, each.index, each.revision};
    if (!x0.contains(window)) continue;
    const vec<std::byte> head = db->read_blob_head(each.stft.id, sizeof(Qtz<double>));
    Qtz<double>          qx;
    std::memcpy(&qx, head.data(), sizeof qx);
    if (qx.pos == x0.at(window)) continue;
    db->remove(qr::delete_from<stft>().where(
      keyed<stft>(std::index_sequence_for<int64_t, int64_t, int64_t>{}, each.extraction, each.index,
                  each.revision)));
    db->remove(qr::delete_from<atomization>().where(
      qr::column<atomization, &atomization::extraction> == each.extraction
      && qr::column<atomization, &atomization::index> == each.index
      && qr::column<atomization, &atomization::revision> == each.revision));
  }
}

void create() {
  db.emplace(env::require(env::var::db_url));
  db->transaction([] {
    migrate<source, decoded, atom_version, dna, ard, extraction, root, stft, atomization, component,
            residual, grouping, component_pair, component_group, linking, component_link,
            column_setting>();
    remove_stft_misaligned_with_root();
    remove_unreferenced_blob();
    if (db->select<ard>().empty())
      db->insert(ard{.id = 0, .name = "default", .factor = matrix::Ard{}});
  });
}

void write_audio(int64_t id, const std::string& domain, const std::string& origin,
                 const io::Probe& probed, int bit_depth, const std::string& stem,
                 std::span<const std::byte> audio) {
  const size_t               i0        = stem.find(env::described);
  const size_t               i1        = stem.find(env::described, i0 + env::described.size());
  const bool                 described = i1 != std::string::npos;
  std::optional<std::string> uploader;
  std::optional<std::string> title;
  if (described) {
    uploader = stem.substr(i0 + env::described.size(), i1 - i0 - env::described.size());
    title    = stem.substr(i1 + env::described.size());
  }
  db->transaction([&] {
    const source row{.id        = id,
                     .domain    = domain,
                     .origin    = origin,
                     .hz        = probed.hz,
                     .channels  = probed.channels,
                     .bit_depth = bit_depth,
                     .audio     = write_blob(audio),
                     .title     = title,
                     .uploader  = uploader,
                     .url       = std::string(env::freesound) + std::to_string(id) + "/",
                     .original  = probed.original};
    db->upsert(row);
    remove_unreferenced_blob();
  });
}

void write_decoded(int64_t source_id, const Spectrum<double>& sample,
                   const Spectrum<double>& log_mono, const Spectrum<double>& linear_mono,
                   const vec<Level>& waveform) {
  db->transaction([&] {
    const spectrum<double> sample_row      = write_spectrum(sample);
    const spectrum<double> log_mono_row    = write_spectrum(log_mono);
    const spectrum<double> linear_mono_row = write_spectrum(linear_mono);
    const blob<vec<Level>> waveform_row    = write_blob(waveform);
    db->upsert(decoded{source_id, sample_row, log_mono_row, linear_mono_row, waveform_row});
    remove_unreferenced_blob();
  });
}

void write_root_extraction(const decoded& source, const dna& factor, const Spectrum<double>& window,
                           std::span<const atom::Atom>     found,
                           std::span<const vec<std::byte>> audio) {
  db->transaction([&] {
    const extraction extraction{.id         = 0,
                                .source_id  = source.source_id,
                                .time_point = std::chrono::system_clock::now(),
                                .factor     = factor.factor,
                                .x0         = window.x0(),
                                .x1         = window.x1(),
                                .version    = thunder_version_id()};
    const int64_t    id = db->insert_returning(extraction);
    vec<root>        rows;
    for (size_t i = 0; i < found.size(); ++i) {
      const blob<std::span<const std::byte>> window_audio
        = write_blob(std::span<const std::byte>(audio[i]));
      const atom::ADSR& adsr = std::get<atom::ADSR>(found[i].component.front());
      rows.push_back(
        root{id, static_cast<int64_t>(i), 0, *found[i].domain, adsr, window_audio, false});
    }
    db->insert(rows);
  });
}

void write_atomization(const root& row, const ard& factor, int channel,
                       const matrix::Relevance& relevance, const vec<>& residual,
                       const atom::Atom& atom, std::span<const vec<std::byte>> audio) {
  db->transaction([&] {
    const atomization a{.id         = 0,
                        .extraction = row.extraction,
                        .index      = row.index,
                        .revision   = row.revision,
                        .time_point = std::chrono::system_clock::now(),
                        .channel    = channel,
                        .factor     = factor.factor,
                        .floor      = relevance.floor,
                        .divergence = relevance.divergence,
                        .iteration  = relevance.iteration,
                        .trace      = relevance.trace,
                        .residual   = residual};

    const int64_t  atomization_id = db->insert_returning(a);
    vec<component> rows;
    for (size_t k = 0; k < atom.atom.size(); ++k) {
      component each{.atomization = atomization_id,
                     .index       = static_cast<int64_t>(k),
                     .lambda      = relevance.lambda[k],
                     .divergence  = relevance.removed[k],
                     .relevant    = relevance.relevant[k],
                     .spectral    = {},
                     .adsr        = {},
                     .spectrum    = std::get<atom::Raw>(atom.atom[k].component[0]).val,
                     .activation  = std::get<atom::Raw>(atom.atom[k].component[1]).val,
                     .octave      = std::get<atom::Raw>(atom.atom[k].component[2]).axis,
                     .density     = std::get<atom::Raw>(atom.atom[k].component[2]).val,
                     .audio       = write_blob(std::span<const std::byte>(audio[k]))};
      for (const atom::Component& part : atom.atom[k].component) {
        if (const atom::Spectral* found = std::get_if<atom::Spectral>(&part))
          each.spectral = *found;
        if (const atom::ADSR* found = std::get_if<atom::ADSR>(&part)) each.adsr = *found;
      }
      rows.push_back(each);
    }
    db->insert(rows);
  });
}

void write_group(const atomization& source, const vec<component>& live, const vvec<>& cosine,
                 const vvec<>& correlation, const vvec<bool>& related, const vec<size_t>& spectrum,
                 const vec<size_t>& event) {
  db->transaction([&] {
    const grouping run{
      .id = 0, .atomization = source.id, .time_point = std::chrono::system_clock::now()};
    const int64_t        id = db->insert_returning(run);
    vec<component_pair>  pair;
    vec<component_group> member;
    for (size_t i = 0; i < live.size(); ++i) {
      for (size_t j = i + 1; j < live.size(); ++j)
        pair.push_back(component_pair{id, live[i].index, live[j].index, cosine[i][j],
                                      correlation[i][j], related[i][j]});
      member.push_back(
        component_group{id, live[i].index, live[spectrum[i]].index, live[event[i]].index});
    }
    db->insert(pair);
    db->insert(member);
  });
}

void write_stft(const root& row, const zvec<double>& stft) {
  db->transaction([&] {
    const blob<zvec<double>> stft_row = write_blob(stft);
    db->upsert(table::stft{row.extraction, row.index, row.revision, stft_row});
    remove_unreferenced_blob();
  });
}

void write_residual(const atomization& source, const zvec<double>& layer) {
  db->transaction([&] {
    const blob<zvec<double>> layer_row = write_blob(layer);
    db->upsert(table::residual{source.id, layer_row});
    remove_unreferenced_blob();
  });
}

vec<component> select_live_component(int64_t atomization_id) {
  const vec<component> part = select<component>(atomization_id);
  vec<component>       live;
  for (const component& each : part)
    if (each.relevant) live.push_back(each);
  return live;
}

std::pair<vvec<>, vvec<>> split_spectrum_activation(const vec<component>& part) {
  vvec spectrum;
  vvec activation;
  for (const component& each : part) {
    spectrum.push_back(each.spectrum);
    activation.push_back(each.activation);
  }
  return {spectrum, activation};
}

int64_t thunder_version_id() {
  const vec<atom_version> stored
    = db->select<atom_version>(qr::column<atom_version, &atom_version::header_blob>
                                 == std::string(qstorm::version::header_blob)
                               && qr::column<atom_version, &atom_version::source_blob>
                                    == std::string(qstorm::version::source_blob));
  for (const atom_version& each : stored) return each.id;
  static constexpr char header_text[] = {
#embed "../atom/atom_thunder.h" suffix(, 0)
  };
  static constexpr char source_text[] = {
#embed "../atom/atom_thunder.cpp" suffix(, 0)
  };
  db->insert(atom_version{.id          = 0,
                          .commit      = qstorm::version::commit,
                          .header_blob = qstorm::version::header_blob,
                          .source_blob = qstorm::version::source_blob,
                          .header_text = header_text,
                          .source_text = source_text});
  return thunder_version_id();
}

struct column_editor {
  const vec<qr::ColumnSpec>& columns;
  const char*                column;
  double                     value;
  size_t                     i = 0;

  template <typename V>
  void cell(V& cell) {
    if constexpr (std::is_assignable_v<V&, double>)
      if (columns[i].name == column) cell = value;
    ++i;
  }

  template <size_t K = 0, typename... S>
  void group(std::variant<S...>& target) {
    if (target.index() == K)
      qr::each_value(std::get<K>(target), *this);
    else
      i += qr::detail::leaf_count<std::variant_alternative_t<K, std::variant<S...>>>();
    if constexpr (K + 1 < sizeof...(S)) group<K + 1>(target);
  }
};

root edit_root(int64_t id, int64_t index, const char* column, double value) {
  const root                 r0      = latest<root>(id, index);
  root                       r1      = r0;
  const vec<qr::ColumnSpec>& columns = qr::columns_of<root>();
  column_editor              visit{columns, column, value};
  qr::each_value(r1, visit);
  r1.revision = r0.revision + 1;
  return r1;
}

void restore_root(int64_t id, int64_t index) {
  const root original = find<root>(id, index, int64_t{0});
  const root r0       = latest<root>(id, index);
  root       r1       = r0;
  r1.domain           = original.domain;
  r1.adsr             = original.adsr;
  r1.audio            = original.audio;
  r1.revision         = r0.revision + 1;
  db->insert(r1);
}

void edit_dna(int64_t id, const char* column, double value) {
  dna                        after   = find<dna>(id);
  const vec<qr::ColumnSpec>& columns = qr::columns_of<dna>();
  column_editor              visit{columns, column, value};
  qr::each_value(after, visit);
  db->upsert(after);
}

int64_t copy_dna(int64_t id) {
  dna copy  = find<dna>(id);
  copy.id   = 0;
  copy.name = std::nullopt;
  return db->insert_returning(copy);
}

size_t clear_root_extraction() {
  const size_t extracts = db->remove(qr::delete_from<extraction>());
  db->remove(qr::delete_from<linking>());
  remove_unreferenced_blob();
  return extracts;
}

vec<std::string> comma_list(const char* text) {
  const std::string_view all(text);
  vec<std::string>       item;
  size_t                 i0 = 0;
  while (all != "all" && i0 <= all.size()) {
    const size_t comma = all.find(',', i0);
    const size_t i1    = comma == std::string_view::npos ? all.size() : comma;
    item.emplace_back(all.substr(i0, i1 - i0));
    i0 = i1 + 1;
  }
  return item;
}

vec<int64_t> comma_number(const char* text) {
  vec<int64_t> number;
  for (const std::string& each : comma_list(text)) number.push_back(std::stoll(each));
  return number;
}

linking link_scope(const char* domain, const char* source, const char* atomization) {
  return {.id          = 0,
          .time_point  = std::chrono::system_clock::now(),
          .domain      = domain,
          .source      = comma_number(source),
          .atomization = comma_number(atomization)};
}

vec<atomization> select_link_window(const linking& scope) {
  const vec<table::source>      every_source     = select<table::source>();
  const vec<extraction>         every_extraction = select<extraction>();
  const vec<table::atomization> every_window     = select<table::atomization>();
  std::map<int64_t, int64_t>    newest_extraction;
  for (const extraction& each : every_extraction)
    newest_extraction.insert_or_assign(each.source_id,
                                       std::max(newest_extraction[each.source_id], each.id));
  std::map<std::pair<int64_t, int64_t>, table::atomization> newest_window;
  for (const table::atomization& each : every_window)
    newest_window.insert_or_assign(std::pair{each.extraction, each.index}, each);
  const vec<std::string>  domain = comma_list(scope.domain.c_str());
  vec<table::atomization> taken;
  for (const table::atomization& each : every_window) {
    const extraction& run = *std::ranges::find(every_extraction, each.extraction, &extraction::id);
    const table::source& src = *std::ranges::find(every_source, run.source_id, &table::source::id);
    const bool           domain_taken
      = domain.empty() || (src.domain && std::ranges::contains(domain, *src.domain));
    const bool source_taken = scope.source.empty() || std::ranges::contains(scope.source, src.id);
    const bool window_taken
      = scope.atomization.empty()
        ? run.id == newest_extraction[src.id]
            && newest_window[std::pair{each.extraction, each.index}].id == each.id
        : std::ranges::contains(scope.atomization, each.id);
    if (domain_taken && source_taken && window_taken) taken.push_back(each);
  }
  return taken;
}

void write_link(const linking& scope, const vec<std::pair<int64_t, int64_t>>& member,
                const vec<size_t>& set) {
  db->transaction([&] {
    const int64_t       id = db->insert_returning(scope);
    vec<component_link> rows;
    for (size_t i = 0; i < member.size(); ++i)
      rows.push_back(component_link{.linking         = id,
                                    .atomization     = member[i].first,
                                    .component       = member[i].second,
                                    .set_atomization = member[set[i]].first,
                                    .set_component   = member[set[i]].second});
    db->insert(rows);
  });
}

template vec<component> select<component>(int64_t);
template source         find<source, int64_t>(int64_t);
template decoded        find<decoded, int64_t>(int64_t);
template dna            find<dna, int64_t>(int64_t);
template ard            find<ard, int64_t>(int64_t);
template atomization    find<atomization, int64_t>(int64_t);
template root           find<root, int64_t, int64_t, int64_t>(int64_t, int64_t, int64_t);
template residual       find<residual, int64_t>(int64_t);
template stft           find<stft, int64_t, int64_t, int64_t>(int64_t, int64_t, int64_t);
template component      find<component, int64_t, int64_t>(int64_t, int64_t);
template root           latest<root, int64_t, int64_t>(int64_t, int64_t);
template source         follow<source, root>(const root&);
template decoded        follow<decoded, root>(const root&);
}
