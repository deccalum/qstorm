#include <charconv>
#include <cstddef>
#include <cstring>
#include <meta>
#include <type_traits>
#include <utility>

#include "cli/commands.h"
#include "io/database.h"

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

template <std::meta::info command, std::size_t... i>
void call(char** argv, std::index_sequence<i...>) {
  [:command:](argument<typename[:std::meta::type_of(std::meta::parameters_of(command)[i]):]>(
      argv[2 + i])...);
}

void dispatch(char** argv) {
  template for (constexpr std::meta::info command : qr::detail::members_of<^^qstorm::cli>()) {
    if constexpr (std::meta::is_function(command)) {
      if (std::meta::identifier_of(command) == argv[1])
        call<command>(argv, std::make_index_sequence<std::meta::parameters_of(command).size()>{});
    }
  }
}

int main(int argc, char** argv) {
  qstorm::io::env::load();
  qstorm::io::table::create();
  if (argc == 1)
    qstorm::cli::studio();
  else
    dispatch(argv);
}
