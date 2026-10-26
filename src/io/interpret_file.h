#pragma once

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>
#include <type_traits>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "io/bridge.h"

namespace qstorm::io {
pcmvec         flac_decode(std::span<const std::byte> bytes, int bit_depth);
vec<std::byte> flac_encode(pcmspan xc, size_t nc, double hz);

template <typename T>
void append_value(vec<std::byte>& bytes, const T& value) {
  const std::byte* first = reinterpret_cast<const std::byte*>(&value);
  bytes.insert(bytes.end(), first, first + sizeof value);
}

template <typename X, typename Y, typename Z>
void append_spectrum(vec<std::byte>& bytes, const Spectrum<X, Y, Z>& spectrum) {
  const uint64_t nc = spectrum.ch.size();
  append_value(bytes, spectrum.qx);
  append_value(bytes, spectrum.qy);
  append_value(bytes, spectrum.qz);
  append_value(bytes, nc);
  for (size_t c = 0; c < nc; ++c) {
    const std::span<const std::byte> plane = std::as_bytes(algorithm::ch(spectrum, c));
    bytes.insert(bytes.end(), plane.begin(), plane.end());
  }
}

template <typename S>
S read_spectrum(const std::byte*& at) {
  using Y = std::remove_pointer_t<typename decltype(S::ch)::value_type>;
  S        spectrum;
  uint64_t nc;
  std::memcpy(&spectrum.qx, at, sizeof spectrum.qx);
  at += sizeof spectrum.qx;
  std::memcpy(&spectrum.qy, at, sizeof spectrum.qy);
  at += sizeof spectrum.qy;
  std::memcpy(&spectrum.qz, at, sizeof spectrum.qz);
  at += sizeof spectrum.qz;
  std::memcpy(&nc, at, sizeof nc);
  at += sizeof nc;
  const size_t plane = spectrum.qx.n * spectrum.qz.n * sizeof(Y);
  for (uint64_t c = 0; c < nc; ++c) {
    spectrum.ch.push_back(const_cast<Y*>(reinterpret_cast<const Y*>(at)));
    at += plane;
  }
  return spectrum;
}

template <typename T>
vec<std::byte> to_bytes(const T& value) {
  vec<std::byte> bytes;
  if constexpr (requires { typename T::value_type; })
    for (const typename T::value_type& each : value) append_spectrum(bytes, each);
  else
    append_spectrum(bytes, value);
  return bytes;
}

template <typename T>
T from_bytes(std::span<const std::byte> bytes) {
  const std::byte* at = bytes.data();
  if constexpr (requires { typename T::value_type; }) {
    const std::byte* end = at + bytes.size();
    T                spectra;
    while (at < end) spectra.push_back(read_spectrum<typename T::value_type>(at));
    return spectra;
  } else {
    return read_spectrum<T>(at);
  }
}
}
