#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <span>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "algorithm/math_scalar.h"

namespace qstorm::algorithm {

inline constexpr size_t waveform_plane = 3;

template <typename X, typename Y, typename Z> Spectrogram<X, Y, Z>
extremum(const Spectrogram<X, Y, Z>& signal, std::vector<std::vector<Y>>& magnitude, size_t plane = waveform_plane) {
  Spectrogram<X, Y, Z> triple{{}, signal.x_ms, signal.y_db, signal.z_hz};
  if (plane == waveform_plane) magnitude.assign(signal.channel.size(), std::vector<Y>());
  for (size_t c = 0; c < signal.channel.size(); ++c) {
    triple.channel.push_back(signal.channel[c]);
    triple.channel.push_back(signal.channel[c]);
    if (plane < waveform_plane) continue;
    magnitude[c].reserve(signal.x_ms.extent);
    for (Y sample : algorithm::channel(signal, c)) magnitude[c].push_back(std::abs(sample));
      triple.channel.push_back(magnitude[c].data());
  }
  return triple;
}

template <typename Y>
Y rms(std::span<const Y> level) {
  Y power = 0;
  for (Y each : level) power += each*each;
  return std::sqrt(power / static_cast<Y>(level.size()));
}

template <typename X, typename Y, typename Z>
void reduce(const Spectrogram<X, Y, Z>& slice, std::vector<std::vector<Y>>& channel, size_t plane = waveform_plane) {
  for (size_t c = 0; c < channel.size(); c += plane) {
    channel[c].push_back(std::ranges::min(algorithm::channel(slice, c)));
    channel[c + 1].push_back(std::ranges::max(algorithm::channel(slice, c + 1)));
    if (plane == waveform_plane)
      channel[c + 2].push_back(rms<Y>(algorithm::channel(slice, c + 2)));
  }
}

template <typename X, typename Y, typename Z> Spectrogram<X, Y, Z>
coarser(const Spectrogram<X, Y, Z>& level, std::vector<std::vector<Y>>& channel, size_t plane = waveform_plane) {
  const X pair = 2*level.begin().time_fraction();
  const size_t n = (level.x_ms.extent + 1) / 2;
  channel.assign(level.channel.size(), std::vector<Y>());
  for (std::vector<Y>& plane : channel) plane.reserve(n);
  for (ms<X> x = level.begin(); x < level.end();) {
    ms<X> next = x;
    ++next;
    ++next;
    reduce(level.window(x, std::min<X>(next, level.end())), channel, plane);
    x = next;
  }
  const Quantization<X> x_ms = Quantization<X>::stepped(channel.front().size(), pair, level.begin());
  return as_spectrogram(channel, x_ms, level.z_hz, level.y_db);
}

template <typename X, typename Y, typename Z> std::vector<Spectrogram<X, Y, Z>>
pyramid(const Spectrogram<X, Y, Z>& signal, std::vector<std::vector<std::vector<Y>>>& channel, size_t plane = waveform_plane) {
  channel.emplace_back();
  std::vector<Spectrogram<X, Y, Z>> level{extremum(signal, channel.back(), plane)};
  while (level.back().x_ms.extent > 1) {
    channel.emplace_back();
    level.push_back(coarser(level.back(), channel.back(), plane));
  }
  return level;
}

template <typename X, typename Y, typename Z> const Spectrogram<X, Y, Z>&
level_for(const std::vector<Spectrogram<X, Y, Z>>& level, X t0, X t1, size_t width) {
  for (size_t k = level.size(); k-- > 0;)
    if (level[k].window(t0, t1).x_ms.extent >= width) return level[k];
  return level.front();
}

template <typename X, typename Y, typename Z> Spectrogram<X, Y, Z>
column(const Spectrogram<X, Y, Z>& level, X t0, X t1, size_t width,
       std::vector<std::vector<Y>>& channel, size_t plane = waveform_plane) {
  const X col = ratio<X>(t1 - t0, width);
  const Quantization<X> x_ms = Quantization<X>::stepped(width, col, t0);
  const X stop = grid(t0, x_ms, width);
  ms<X> last = level.end();
  --last;
  channel.assign(level.channel.size(), std::vector<Y>());
  for (ms<X> from{t0, &x_ms}; from < stop; ++from) {
    ms<X> until = from;
    ++until;
    const Spectrogram<X, Y, Z> slice = level.window(from, until);
    if (slice.x_ms.extent)
      reduce(slice, channel, plane);
    else
      reduce(level.frame(std::min<X>(from, last)), channel, plane);
  }
  return as_spectrogram(channel, x_ms, level.z_hz, level.y_db);
}

template <typename X, typename Y, typename Z> Spectrogram<X, Y, Z>
view(const std::vector<Spectrogram<X, Y, Z>>& levels, X t0, X t1, size_t width,
     std::vector<std::vector<Y>>& channel, size_t plane = waveform_plane) {
  const Spectrogram<X, Y, Z>& sample = levels.front();
  if (sample.window(t0, t1).x_ms.extent >= width)
    return column(level_for(levels, t0, t1, width), t0, t1, width, channel, plane);
  ms<X> lo{t0, &sample.x_ms};
  ms<X> hi{t1, &sample.x_ms};
  --lo;
  ++hi;
  return sample.window(std::max<X>(lo, sample.begin()), std::min<X>(hi, sample.end()));
}

}
