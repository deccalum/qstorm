#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <expected>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <meta>
#include <qr/qr.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "io/bridge.h"

extern char** environ;

namespace qstorm::io {

namespace env {
std::string require(const char* key) {
  const char* value = std::getenv(key);
  return value ? std::string(value)
               : throw std::runtime_error(std::string("missing: ") + key);
}

static std::filesystem::path resolve(const std::filesystem::path& dotenv) {
  if (dotenv.is_absolute()) return dotenv;
  for (std::filesystem::path directory = std::filesystem::current_path();;
       directory = directory.parent_path()) {
    if (std::filesystem::exists(directory / dotenv)) return directory / dotenv;
    if (directory == directory.parent_path()) return dotenv;
  }
}

static std::string upper(std::string_view name) {
  std::string key;
  for (char letter : name) key += static_cast<char>(std::toupper(static_cast<unsigned char>(letter)));
  return key;
}

void load(const std::filesystem::path& dotenv) {
  std::ifstream file(resolve(dotenv));
  std::string   line;
  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') continue;
    const size_t equals = line.find('=');
    if (equals == std::string::npos) continue;
    ::setenv(line.substr(0, equals).c_str(), line.substr(equals + 1).c_str(), 0);
  }
  template for (constexpr std::meta::info field : qr::detail::nonstatic_members_of<Environment>()) {
    const std::string key = upper(std::meta::identifier_of(field));
    value.[:field:] = require(key.c_str());
  }
}
}

std::vector<char*> argv(std::vector<std::string>& word) {
  std::vector<char*> pointer;
  for (std::string& each : word) pointer.push_back(each.data());
  pointer.push_back(nullptr);
  return pointer;
}

pid_t spawn(std::vector<std::string>& word, int child_end, int child_fd, int parent_end) {
  posix_spawn_file_actions_t actions;
  posix_spawn_file_actions_init(&actions);
  posix_spawn_file_actions_adddup2(&actions, STDERR_FILENO, STDOUT_FILENO);
  posix_spawn_file_actions_adddup2(&actions, child_end, child_fd);
  posix_spawn_file_actions_addclose(&actions, parent_end);
  std::vector<char*> arguments = argv(word);
  pid_t pid;
  posix_spawnp(&pid, arguments.front(), &actions, nullptr, arguments.data(), environ);
  posix_spawn_file_actions_destroy(&actions);
  ::close(child_end);
  return pid;
}

std::expected<std::vector<std::byte>, int> capture(std::vector<std::string>& word) {
  int pipe_fd[2];
  ::pipe(pipe_fd);
  const pid_t pid = spawn(word, pipe_fd[1], STDOUT_FILENO, pipe_fd[0]);
  std::vector<std::byte> output;
  std::array<std::byte, 1 << 16> chunk;
  for (ssize_t got = ::read(pipe_fd[0], chunk.data(), chunk.size()); got > 0;
       got = ::read(pipe_fd[0], chunk.data(), chunk.size()))
    output.insert(output.end(), chunk.begin(), chunk.begin() + got);
  ::close(pipe_fd[0]);
  int status;
  waitpid(pid, &status, 0);
  const int code = WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
  if (code != 0) return std::unexpected(code);
  return output;
}

std::expected<std::pair<int, int>, int> probe(const std::filesystem::path& file) {
  std::vector<std::string> word = env::probe_stream.words({file.string()});
  const std::expected<std::vector<std::byte>, int> output = capture(word);
  if (!output) return std::unexpected(output.error());
  const char* text = reinterpret_cast<const char*>(output->data());
  const char* end = text + output->size();
  int hz;
  int channels;
  const std::from_chars_result rate = std::from_chars(text, end, hz);
  std::from_chars(rate.ptr + 1, end, channels);
  return std::pair{hz, channels};
}

std::expected<void, int> transcode(const std::filesystem::path& file) {
  std::vector<std::string> word = env::ffmpeg::words<env::pcm::shifted>(env::transcode_raw, file.string());
  const std::expected<std::vector<std::byte>, int> output = capture(word);
  if (!output) return std::unexpected(output.error());
  env::pcm::interleaved.resize(output->size() / sizeof(env::pcm::shifted));
  std::memcpy(env::pcm::interleaved.data(), output->data(), output->size());
  return {};
}

bool audio(const std::filesystem::path& file) {
  std::string extension = file.extension().string();
  for (char& letter : extension) letter = static_cast<char>(std::tolower(static_cast<unsigned char>(letter)));
  return std::ranges::contains(env::audioext_passthrough, extension);
}

int python(std::string_view script, std::initializer_list<std::string_view> args) {
  std::vector<std::string> word{env::value.python, std::string(script)};
  for (std::string_view arg : args) word.emplace_back(arg);
  std::vector<char*> argument = argv(word);

  pid_t pid;
  posix_spawn(&pid, env::value.python.c_str(), nullptr, nullptr, argument.data(), environ);
  int status;
  waitpid(pid, &status, 0);
  return WEXITSTATUS(status);
}

}
