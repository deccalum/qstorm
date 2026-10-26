#include "io/bridge.h"

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
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

extern char** environ;

namespace qstorm::env {
std::string require(const char* key) {
  const char* value = std::getenv(key);
  return value ? std::string(value) : throw std::runtime_error(std::string("missing: ") + key);
}

static std::filesystem::path resolve(const std::filesystem::path& dotenv) {
  if (dotenv.is_absolute()) return dotenv;
  for (std::filesystem::path dir = std::filesystem::current_path();; dir = dir.parent_path()) {
    if (std::filesystem::exists(dir / dotenv)) return dir / dotenv;
    if (dir == dir.parent_path()) return dotenv;
  }
}

void load(const std::filesystem::path& dotenv) {
  std::ifstream file(resolve(dotenv));
  std::string   line;
  while (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') continue;
    const size_t equal = line.find('=');
    if (equal == std::string::npos) continue;
    ::setenv(line.substr(0, equal).c_str(), line.substr(equal + 1).c_str(), 0);
  }
}
}

namespace qstorm::io {
vec<char*> argv(vec<std::string>& word) {
  vec<char*> pointer;
  for (std::string& each : word) pointer.push_back(each.data());
  pointer.push_back(nullptr);
  return pointer;
}

pid_t spawn(vec<std::string>& word, int child_end, int child_fd, int parent_end) {
  posix_spawn_file_actions_t action;
  posix_spawn_file_actions_init(&action);
  posix_spawn_file_actions_adddup2(&action, STDERR_FILENO, STDOUT_FILENO);
  posix_spawn_file_actions_adddup2(&action, child_end, child_fd);
  posix_spawn_file_actions_addclose(&action, parent_end);
  vec<char*> args = argv(word);
  pid_t      pid;
  posix_spawnp(&pid, args.front(), &action, nullptr, args.data(), environ);
  posix_spawn_file_actions_destroy(&action);
  ::close(child_end);
  return pid;
}

std::expected<vec<std::byte>, int> capture(vec<std::string>& word) {
  int pipe_fd[2];
  ::pipe(pipe_fd);
  const pid_t                    pid = spawn(word, pipe_fd[1], STDOUT_FILENO, pipe_fd[0]);
  vec<std::byte>                 output;
  std::array<std::byte, 1 << 16> chunk;
  for (ssize_t got = ::read(pipe_fd[0], chunk.data(), chunk.size()); got > 0;
       got         = ::read(pipe_fd[0], chunk.data(), chunk.size()))
    output.insert(output.end(), chunk.begin(), chunk.begin() + got);
  ::close(pipe_fd[0]);
  int status;
  waitpid(pid, &status, 0);
  const int code = WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
  if (code != 0) return std::unexpected(code);
  return output;
}

std::map<std::string, std::string> probe_field(std::string_view text) {
  std::map<std::string, std::string> field;
  size_t                             i0 = 0;
  while (i0 < text.size()) {
    const size_t           newline = text.find('\n', i0);
    const size_t           i1      = newline == std::string_view::npos ? text.size() : newline;
    const std::string_view line    = text.substr(i0, i1 - i0);
    const size_t           equal   = line.find('=');
    std::string_view       value   = line.substr(equal + 1);
    if (value.starts_with('"')) value = value.substr(1, value.size() - 2);
    field.emplace(std::string(line.substr(0, equal)), std::string(value));
    i0 = i1 + 1;
  }
  return field;
}

template <typename T>
std::optional<T>
probe_num(const std::map<std::string, std::string>& field, const std::string& key) {
  const std::string&           text = field.at(key);
  T                            value{};
  const std::from_chars_result read
    = std::from_chars(text.data(), text.data() + text.size(), value);
  std::optional<T> found;
  if (read.ec == std::errc{} && value != T{}) found = value;
  return found;
}

std::string joined_entries(std::span<const std::string_view> entry) {
  std::string joined(entry.front());
  for (std::string_view each : entry.subspan(1)) {
    joined += env::ffmpeg::entry_sep;
    joined += each;
  }
  return joined;
}

std::expected<Probe, int> probe(const std::filesystem::path& file) {
  using namespace env::ffmpeg;
  const std::string stream_list = joined_entries(env::stream_entry);
  const std::string fmt_list    = joined_entries(env::fmt_entry);
  const std::string entry
    = std::string(stream_section) + stream_list + std::string(fmt_section) + fmt_list;
  vec<std::string> word = env::probe_stream.words({show_entries, entry, file.string()});
  const std::expected<vec<std::byte>, int> output = capture(word);
  if (!output) return std::unexpected(output.error());
  const std::string_view text(reinterpret_cast<const char*>(output->data()), output->size());
  const std::map<std::string, std::string> field = probe_field(text);
  const std::string                        stream(stream_prefix);
  const std::string                        format(fmt_prefix);

  const std::optional<int>     raw_depth    = probe_num<int>(field, stream + std::string(bprs));
  const std::optional<int>     sample_depth = probe_num<int>(field, stream + std::string(bps));
  const std::optional<int64_t> stream_rate
    = probe_num<int64_t>(field, stream + std::string(bit_rate));
  const std::optional<int64_t> fmt_rate = probe_num<int64_t>(field, format + std::string(bit_rate));
  const std::string&           layout   = field.at(stream + std::string(ch_layout));
  const double                 sec      = *probe_num<double>(field, format + std::string(duration));
  const Original original{.container  = field.at(format + std::string(fmt_name)),
                          .codec      = field.at(stream + std::string(codec_name)),
                          .sample_fmt = field.at(stream + std::string(sample_fmt)),
                          .bit_depth  = raw_depth ? raw_depth : sample_depth,
                          .bit_rate   = stream_rate ? stream_rate : fmt_rate,
                          .ch_layout  = layout == absent ? std::nullopt : std::optional{layout},
                          .length     = sec * 1000.0,
                          .size       = *probe_num<int64_t>(field, format + std::string(size))};

  const int hz = *probe_num<int>(field, stream + std::string(rate));
  const int c  = *probe_num<int>(field, stream + std::string(channels));
  return Probe{hz / 1000.0, c, original};
}

std::expected<pcmvec, int> transcode(const std::filesystem::path& file) {
  vec<std::string> word = env::ffmpeg::words<env::pcm::shifted>(env::transcode_raw, file.string());
  const std::expected<vec<std::byte>, int> output = capture(word);
  if (!output) return std::unexpected(output.error());
  pcmvec xc(output->size() / sizeof(env::pcm::shifted));
  std::memcpy(xc.data(), output->data(), output->size());
  return xc;
}

bool audio(const std::filesystem::path& file) {
  std::string extension = file.extension().string();
  for (char& letter : extension)
    letter = static_cast<char>(std::tolower(static_cast<unsigned char>(letter)));
  return std::ranges::contains(env::audioext_passthrough, extension);
}

std::expected<int, int>
python(std::string_view script, std::initializer_list<std::string_view> args) {
  vec<std::string> word{env::require(env::var::python), std::string(script)};
  for (std::string_view arg : args) word.emplace_back(arg);
  vec<char*> arg = argv(word);

  pid_t     pid;
  const int spawned
    = posix_spawn(&pid, word.front().c_str(), nullptr, nullptr, arg.data(), environ);
  if (spawned != 0) return std::unexpected(spawned);
  int status;
  waitpid(pid, &status, 0);
  return WEXITSTATUS(status);
}
}
