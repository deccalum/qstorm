#pragma once

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <functional>
#include <optional>
#include <semaphore>
#include <thread>
#include <tuple>
#include <type_traits>

#include "algorithm/math_constant.h"

namespace qstorm::algorithm {
inline constexpr size_t thread_work  = size_t{1} << 15;
inline constexpr size_t thread_piece = 8;

struct Workers {
  Workers() {
    const size_t helper = std::max<size_t>(1, std::thread::hardware_concurrency()) - 1;
    for (size_t t = 0; t < helper; ++t) thread.emplace_back(&Workers::serve, this);
  }

  ~Workers() {
    stopping = true;
    ready.release(thread.size());
  }

  size_t size() const { return thread.size() + 1; }

  void run(size_t nt, size_t npiece, std::function_ref<void(size_t) const> fn) {
    job       = fn;
    piece     = npiece;
    next      = 0;
    remaining = nt - 1;
    ready.release(nt - 1);
    take();
    for (size_t left = remaining.load(); left != 0; left = remaining.load()) remaining.wait(left);
  }

 private:
  void serve() {
    while (true) {
      ready.acquire();
      if (stopping) return;
      take();
      if (remaining.fetch_sub(1) == 1) remaining.notify_one();
    }
  }

  void take() {
    for (size_t t = next.fetch_add(1); t < piece; t = next.fetch_add(1)) (*job)(t);
  }

  std::counting_semaphore<> ready{0};
  std::atomic<size_t>       next{0};
  std::atomic<size_t>       remaining{0};
  std::atomic<bool>         stopping{false};
  size_t                    piece = 0;

  std::optional<std::function_ref<void(size_t) const>> job;

  vec<std::jthread> thread;
};

inline Workers workers;

template <typename... Arg>
struct PieceCall {
  void (*range)(size_t, size_t, Arg...);
  size_t             n;
  size_t             piece;
  std::tuple<Arg...> arg;

  void operator()(size_t t) const {
    const size_t                     i0 = t * n / piece;
    const size_t                     i1 = (t + 1) * n / piece;
    const std::tuple<size_t, size_t> bound{i0, i1};
    std::apply(range, std::tuple_cat(bound, arg));
  }
};

template <typename... Arg>
void in_parallel(size_t n, size_t work, void (*range)(size_t, size_t, Arg...),
                 std::type_identity_t<Arg>... arg) {
  const size_t            needed = std::max<size_t>(1, work / thread_work);
  const size_t            nt     = std::min({workers.size(), needed, n});
  const size_t            piece  = std::min(n, nt * thread_piece);
  const PieceCall<Arg...> call{range, n, piece, {arg...}};
  workers.run(nt, piece, call);
}
}
