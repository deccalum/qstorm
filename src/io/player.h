#pragma once
#include <poll.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <format>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "io/bridge.h"

namespace qstorm::io::player {

inline int fd = -1;
inline size_t last_count = 0;
inline size_t last_sample = 0;

template <typename T>
struct Sink {
  const env::ProcessArg& profile;
  mutable pid_t pid = 0;

  int open() const {
    int pipe_fd[2];
    ::pipe(pipe_fd);
    std::vector<std::string> word = env::pw_cat::words<T>(profile);
    pid = spawn(word, pipe_fd[0], STDIN_FILENO, pipe_fd[1]);
    return pipe_fd[1];
  }

  void write(int fd, std::span<const T> block, size_t) const {
    const std::span<const std::byte> byte = std::as_bytes(block);
    for (size_t sent = 0; sent < byte.size();) {
      const ssize_t written = ::write(fd, byte.data() + sent, byte.size() - sent);
      sent += static_cast<size_t>(written);
    }
  }

  void close(int fd) const {
    ::close(fd);
    waitpid(pid, nullptr, 0);
  }

  void pause()  const { ::kill(pid, SIGSTOP); }
  void resume() const { ::kill(pid, SIGCONT); }
  void abort()  const { ::kill(pid, SIGTERM); }
  void gain()   const;
};

template <typename T>
struct Source {
  const env::ProcessArg& profile;
  mutable pid_t pid = 0;

  int open() const {
    int pipe_fd[2];
    ::pipe(pipe_fd);
    std::vector<std::string> word = env::pw_cat::words<T>(profile);
    pid = spawn(word, pipe_fd[1], STDOUT_FILENO, pipe_fd[0]);
    return pipe_fd[0];
  }

  size_t read(int fd, std::span<T> block) const {
    const std::span<std::byte> byte = std::as_writable_bytes(block);
    size_t received = 0;
    while (received < byte.size()) {
      const ssize_t got = ::read(fd, byte.data() + received, byte.size() - received);
      if (got == 0) break;
      received += static_cast<size_t>(got);
    }
    return received / sizeof(T);
  }

  void close(int fd) const {
    ::close(fd);
    waitpid(pid, nullptr, 0);
  }

  void abort() const { ::kill(pid, SIGTERM); }
};

template <typename D, typename T>
concept Device = requires(const D& d, int fd, std::span<const T> block, size_t count) {
  { d.open() } -> std::same_as<int>;
  d.write(fd, block, count);
  d.close(fd);
  d.pause();
  d.resume();
  d.abort();
};

template <Device<env::pcm::shifted> D>
void close(const D& device) {
  if (fd >= 0) device.close(fd);
  fd = -1;
}

template <typename Y>
void report(const std::vector<algorithm::db<Y>>& level) {
  std::string text{env::control::level};
  for (const algorithm::db<Y>& each : level) text += std::format(" {:.1f}", each.y);
  text += '\n';
  ::write(STDOUT_FILENO, text.data(), text.size());
}

inline bool ready() {
  pollfd watch{STDIN_FILENO, POLLIN, 0};
  return ::poll(&watch, 1, 0) > 0;
}

inline std::string line() {
  std::string text;
  char each;
  while (::read(STDIN_FILENO, &each, 1) == 1 && each != '\n') text += each;
  return text;
}

inline std::vector<std::string> words(const std::string& order) {
  std::vector<std::string> found;
  size_t first = 0;
  while (first < order.size()) {
    const size_t space = order.find(' ', first);
    const size_t last = space == std::string::npos ? order.size() : space;
    if (last > first) found.push_back(order.substr(first, last - first));
    first = last + 1;
  }
  return found;
}

template <Device<env::pcm::shifted> D, typename Y>
struct PlaybackControl {
  const D& device;
  algorithm::Quantization<Y> y;
  size_t master;
  std::vector<size_t> step;
  bool mono;
  std::array<algorithm::db<Y>, env::gain_step> volume =
      algorithm::gain_step(y, std::make_index_sequence<env::gain_step>{});

  Y linear(size_t at) const {
    if (at == 0) return 0;
    const Y db = volume[at].y;
    return static_cast<Y>(std::pow(10.0, db / 20.0));
  }

  Y gain(size_t c) const {
    const Y whole = linear(master);
    const Y channel = linear(step[c]);
    return whole*channel;
  }

  std::vector<Y> gains() const {
    std::vector<Y> each;
    for (size_t c = 0; c < step.size(); ++c) each.push_back(gain(c));
    return each;
  }

  bool apply(const std::string& order) {
    const std::vector<std::string> word = words(order);
    if (word.empty()) return false;
    if (word[0] == env::control::master) {
      master = std::min<size_t>(std::stoul(word[1]), env::gain_step - 1);
      return true;
    }
    if (word[0] == env::control::volume) {
      const size_t c = std::stoul(word[1]);
      step[c] = std::min<size_t>(std::stoul(word[2]), env::gain_step - 1);
      return true;
    }
    if (word[0] == env::control::mono) {
      mono = word[1] == env::control::on;
      return true;
    }
    return false;
  }

  bool held() {
    device.pause();
    while (true) {
      const std::string order = line();
      if (order == env::control::stop) {
        device.resume();
        return true;
      }
      if (order == env::control::pause) {
        device.resume();
        return false;
      }
      if (!apply(order)) {
        device.resume();
        return true;
      }
    }
  }

  bool interrupted() {
    while (ready()) {
      const std::string order = line();
      if (order == env::control::stop) return true;
      if (order == env::control::pause) {
        if (held()) return true;
        continue;
      }
      if (!apply(order)) return true;
    }
    return false;
  }
};

template <Device<env::pcm::shifted> D, typename Y>
bool play(const D& device, PlaybackControl<D, Y>& control) {
  const algorithm::Quantization<double> x_ms{env::pcm::interleaved.size() / env::pcm::count, env::pcm::ms, 0.0};
  if (fd < 0 || env::pcm::count != last_count || env::pcm::sample != last_sample) {
    close(device);
    fd = device.open();
    last_count = env::pcm::count;
    last_sample = env::pcm::sample;
  }
  const size_t block = x_ms.frame(env::play_block_ms);
  const std::span<const env::pcm::shifted> interleaved(env::pcm::interleaved);
  for (size_t t0 = 0; t0 < x_ms.extent; t0 += block) {
    if (control.interrupted()) return false;
    const size_t frames = std::min(block, x_ms.extent - t0);
    const std::span<const env::pcm::shifted> source = interleaved.subspan(t0*env::pcm::count, frames*env::pcm::count);
    const std::vector<Y> gain = control.gains();
    const std::vector<env::pcm::shifted> mixed = algorithm::mix<env::pcm::shifted, Y>(source, env::pcm::count, gain, control.mono);
    const std::vector<Y> scale = algorithm::amplitude<Y>(mixed);
    const std::vector<std::vector<Y>> channel = algorithm::deinterleave<Y>(scale, env::pcm::count);
    std::vector<algorithm::db<Y>> level;
    for (const std::vector<Y>& row : channel) {
      const auto [lo, hi] = std::ranges::minmax(row);
      const Y peak = std::max(-lo, hi);
      level.push_back({algorithm::db<Y>::ratio(peak), &control.y});
    }
    report(level);
    device.write(fd, mixed, env::pcm::count);
  }
  return true;
}

}
