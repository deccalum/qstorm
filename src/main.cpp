#include <algorithm>
#include <array>
#include <charconv>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <expected>
#include <string_view>
#include <type_traits>
#include <utility>

#include "cli/commands.h"
#include "io/bridge.h"
#include "io/database.h"

namespace qstorm {
template <typename T>
T argument(const char* text) {
  if constexpr (std::is_same_v<T, const char*>) {
    return text;
  } else {
    T value{};
    std::from_chars(text, text + std::strlen(text), value);
    return value;
  }
}

template <typename... P>
consteval size_t parameter_count(void (*)(P...)) {
  return sizeof...(P);
}

template <typename... P, size_t... i>
void call_with_argv(void (*command)(P...), char** argv, std::index_sequence<i...>) {
  command(argument<P>(argv[2 + i])...);
}

template <auto command>
void run_command(char** argv) {
  call_with_argv(command, argv, std::make_index_sequence<parameter_count(command)>{});
}

struct Command {
  std::string_view name;
  void (*run)(char**);
};

inline constexpr std::array command{
  Command{          "ingest_source",           run_command<&cli::library::ingest_source>},
  Command{          "decode_source",           run_command<&cli::library::decode_source>},
  Command{              "root_atom",               run_command<&cli::extract::root_atom>},
  Command{                "atomize",                 run_command<&cli::extract::atomize>},
  Command{        "group_component",         run_command<&cli::extract::group_component>},
  Command{         "link_component",          run_command<&cli::extract::link_component>},
  Command{"subtract_reconstruction", run_command<&cli::extract::subtract_reconstruction>},
  Command{             "play_audio",               run_command<&cli::render::play_audio>},
  Command{     "play_component_set",       run_command<&cli::render::play_component_set>},
  Command{          "play_residual",            run_command<&cli::render::play_residual>},
  Command{  "render_component_list",    run_command<&cli::render::render_component_list>},
  Command{        "render_waveform",          run_command<&cli::render::render_waveform>},
  Command{             "write_stft",               run_command<&cli::render::write_stft>},
  Command{            "render_stft",              run_command<&cli::render::render_stft>},
  Command{        "render_residual",          run_command<&cli::render::render_residual>},
  Command{  "render_reconstruction",    run_command<&cli::render::render_reconstruction>},
  Command{       "render_component",         run_command<&cli::render::render_component>},
  Command{             "read_table",                run_command<&cli::query::read_table>},
  Command{              "edit_root",                 run_command<&cli::query::edit_root>},
  Command{           "restore_root",              run_command<&cli::query::restore_root>},
  Command{               "edit_dna",                  run_command<&cli::query::edit_dna>},
  Command{               "copy_dna",                  run_command<&cli::query::copy_dna>},
  Command{            "clear_cache",               run_command<&cli::query::clear_cache>},
  Command{  "clear_root_extraction",     run_command<&cli::query::clear_root_extraction>},
};
}

int main(int argc, char** argv) {
  using namespace qstorm;
  env::load();
  table::create();
  if (argc == 1) {
    const std::expected<int, int> exited = io::studio();
    if (!exited)
      std::fprintf(stderr, "%s: %s\n", env::require(env::var::python).c_str(),
                   std::strerror(exited.error()));
  } else {
    const Command* chosen = std::ranges::find(command, std::string_view(argv[1]), &Command::name);
    chosen->run(argv);
  }
}
