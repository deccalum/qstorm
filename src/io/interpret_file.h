#pragma once

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <span>
#include <string>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "io/bridge.h"

namespace qstorm::io {

void decode_flac(std::span<const std::byte> bytes);
void decode_flac_range(const std::string& path, double start, double end);
std::vector<std::byte> encode_flac(std::span<const env::pcm::shifted> interleaved, size_t count, int hz);

inline void write_bytes(const std::string& path, std::span<const std::byte> bytes) {
  std::ofstream file(path, std::ios::binary);
  file.write(reinterpret_cast<const char*>(bytes.data()),
             static_cast<std::streamsize>(bytes.size()));
}

template <typename X, typename Y, typename Z>
void write_spectrogram(std::ofstream& file, const algorithm::Spectrogram<X, Y, Z>& spectrum) {
  const uint64_t count = spectrum.channel.size();
  file.write(reinterpret_cast<const char*>(&spectrum.x_ms), sizeof spectrum.x_ms);
  file.write(reinterpret_cast<const char*>(&spectrum.y_db), sizeof spectrum.y_db);
  file.write(reinterpret_cast<const char*>(&spectrum.z_hz), sizeof spectrum.z_hz);
  file.write(reinterpret_cast<const char*>(&count), sizeof count);
  for (size_t c = 0; c < count; ++c) {
    const std::span<Y> plane = algorithm::channel(spectrum, c);
    file.write(reinterpret_cast<const char*>(plane.data()),
               static_cast<std::streamsize>(plane.size_bytes()));
  }
}

template <typename X, typename Y, typename Z>
algorithm::Spectrogram<X, Y, Z> read_spectrogram(const std::byte*& at) {
  algorithm::Spectrogram<X, Y, Z> spectrum;
  uint64_t count;
  std::memcpy(&spectrum.x_ms, at, sizeof spectrum.x_ms);
  at += sizeof spectrum.x_ms;
  std::memcpy(&spectrum.y_db, at, sizeof spectrum.y_db);
  at += sizeof spectrum.y_db;
  std::memcpy(&spectrum.z_hz, at, sizeof spectrum.z_hz);
  at += sizeof spectrum.z_hz;
  std::memcpy(&count, at, sizeof count);
  at += sizeof count;
  const size_t plane = spectrum.x_ms.extent * spectrum.z_hz.extent * sizeof(Y);
  for (uint64_t c = 0; c < count; ++c) {
    spectrum.channel.push_back(const_cast<Y*>(reinterpret_cast<const Y*>(at)));
    at += plane;
  }
  return spectrum;
}

template <typename X, typename Y, typename Z>
std::vector<algorithm::Spectrogram<X, Y, Z>> map_spectrogram(const std::string& path) {
  const int fd = ::open(path.c_str(), O_RDONLY);
  struct stat info;
  ::fstat(fd, &info);
  const size_t size = static_cast<size_t>(info.st_size);
  void* mapped = ::mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0);
  ::close(fd);
  const std::byte* at = static_cast<const std::byte*>(mapped);
  const std::byte* end = at + size;
  std::vector<algorithm::Spectrogram<X, Y, Z>> spectra;
  while (at < end) spectra.push_back(read_spectrogram<X, Y, Z>(at));
  return spectra;
}

}
