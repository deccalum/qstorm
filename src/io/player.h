#pragma once
#include <poll.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <format>
#include <span>
#include <string>
#include <utility>

#include "algorithm/axis_quantizer.h"
#include "io/bridge.h"

namespace qstorm::io::player {

inline bool   opened      = false;
inline size_t last_nc     = 0;
inline double last_sample = 0.0;

struct Sink {
  const env::ProcessArg& profile;
  mutable pid_t          pid = 0;
  mutable int            fd  = -1;

  void open(double hz, size_t nc) const {
    int pipe_fd[2];
    ::pipe(pipe_fd);
    vec<std::string> word = env::ffplay::words<env::pcm::shifted>(profile, hz, nc);
    pid                   = spawn(word, pipe_fd[0], STDIN_FILENO, pipe_fd[1]);
    fd                    = pipe_fd[1];
  }

  void write(pcmspan block) const {
    const std::span<const std::byte> byte = std::as_bytes(block);
    for (size_t sent = 0; sent < byte.size();) {
      const ssize_t written = ::write(fd, byte.data() + sent, byte.size() - sent);
      sent += static_cast<size_t>(written);
    }
  }

  void close() const {
    ::close(fd);
    waitpid(pid, nullptr, 0);
  }

  void pause() const { ::kill(pid, SIGSTOP); }
  void resume() const { ::kill(pid, SIGCONT); }

  void abort() const {
    ::kill(pid, SIGTERM);
    ::kill(pid, SIGCONT);
  }
};

template <typename D, typename T>
concept Device = requires(const D& d, std::span<const T> block, double hz, size_t nc) {
  d.open(hz, nc);
  d.write(block);
  d.close();
  d.pause();
  d.resume();
  d.abort();
};

template <Device<env::pcm::shifted> D>
void close(const D& device) {
  if (opened) device.close();
  opened = false;
}

template <typename Y>
void report(const vec<algorithm::db<Y>>& level) {
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
  char        each;
  while (::read(STDIN_FILENO, &each, 1) == 1 && each != '\n') text += each;
  return text;
}

inline vec<std::string> words(const std::string& order) {
  vec<std::string> found;
  size_t           i0 = 0;
  while (i0 < order.size()) {
    const size_t space = order.find(' ', i0);
    const size_t i1    = space == std::string::npos ? order.size() : space;
    if (i1 > i0) found.push_back(order.substr(i0, i1 - i0));
    i0 = i1 + 1;
  }
  return found;
}

template <Device<env::pcm::shifted> D, typename Y>
struct PlaybackControl {
  const D&    device;
  Qtz<Y>      qy;
  size_t      master;
  vec<size_t> step;
  bool        mono;

  std::array<algorithm::db<Y>, env::gain_step> volume
    = algorithm::gain_step(qy, std::make_index_sequence<env::gain_step>{});

  Y linear(size_t at) const {
    if (at == 0) return 0;
    const Y db = volume[at].y;
    return static_cast<Y>(std::pow(10.0, db / 20.0));
  }

  Y gain(size_t c) const {
    const Y whole   = linear(master);
    const Y channel = linear(step[c]);
    return whole * channel;
  }

  vec<Y> gains() const {
    vec<Y> each;
    for (size_t c = 0; c < step.size(); ++c) each.push_back(gain(c));
    return each;
  }

  bool apply(const std::string& order) {
    const vec<std::string> word = words(order);
    if (word.empty()) return false;
    if (word[0] == env::control::master) {
      master = std::min<size_t>(std::stoul(word[1]), env::gain_step - 1);
      return true;
    }
    if (word[0] == env::control::volume) {
      const size_t c = std::stoul(word[1]);
      step[c]        = std::min<size_t>(std::stoul(word[2]), env::gain_step - 1);
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
      if (order == env::control::stop) return true;
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
bool play(const D& device, pcmspan xc, size_t nc, double hz, PlaybackControl<D, Y>& control) {
  const Qtz<double> qx{xc.size() / nc, hz, 0.0};
  if (!opened || nc != last_nc || hz != last_sample) {
    close(device);
    device.open(hz, nc);
    opened      = true;
    last_nc     = nc;
    last_sample = hz;
  }
  const size_t block = qx.frame(env::play_block);
  for (size_t t0 = 0; t0 < qx.n; t0 += block) {
    if (control.interrupted()) return false;
    const size_t  nx      = std::min(block, qx.n - t0);
    const pcmspan source  = xc.subspan(t0 * nc, nx * nc);
    const vec<Y>  gain    = control.gains();
    const pcmvec  mixed   = algorithm::xc_mix<env::pcm::shifted, Y>(source, nc, gain, control.mono);
    const vvec<Y> channel = algorithm::ch_amp<Y>(mixed, nc);
    vec<algorithm::db<Y>> level;
    for (const vec<Y>& row : channel) {
      const auto [y0, y1] = std::ranges::minmax(row);
      const Y peak        = std::max(-y0, y1);
      level.push_back({algorithm::db<Y>::ratio(peak), &control.qy});
    }
    report(level);
    device.write(mixed);
  }
  return true;
}

inline void stream(pcmspan xc, size_t nc, double hz, const Qtz<float>& qy) {
  const Sink                   sink{env::playback_raw};
  const size_t                 full = env::gain_step - 1;
  PlaybackControl<Sink, float> control{sink, qy, full, vec<size_t>(nc, full), false};
  const bool                   finished = play(sink, xc, nc, hz, control);
  if (!finished) sink.abort();
  close(sink);
}

}
