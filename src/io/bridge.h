#pragma once

#include <sys/types.h>
#include <unistd.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <filesystem>
#include <initializer_list>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace qstorm::io {

namespace env {
struct Environment {
  std::string python;
  std::string pw_cat;
  std::string cache_dir;
  std::string ffmpeg;
  std::string ffprobe;
  std::string studio;
  std::string database_url;
};
inline Environment value;

void load(const std::filesystem::path& dotenv = ".env");
std::string require(const char* key);

inline constexpr float play_block_ms = 100;
inline constexpr size_t gain_step = 500;
inline constexpr int ingest_bit = 24;

inline constexpr std::string_view pipe_end = "-";
inline constexpr std::string_view module = "-m";
inline constexpr std::string_view freesound = "https://freesound.org/s/";
inline constexpr std::string_view described = "__";

namespace pw_cat {
inline constexpr std::string_view playback = "--playback";
inline constexpr std::string_view record = "--record";
inline constexpr std::string_view raw = "--raw";
inline constexpr std::string_view format = "--format";
inline constexpr std::string_view rate = "--rate";
inline constexpr std::string_view channels = "--channels";
inline constexpr std::string_view volume = "--volume";
inline constexpr std::string_view latency = "--latency";
inline constexpr std::string_view target = "--target";
}

namespace control {
inline constexpr std::string_view master = "master";
inline constexpr std::string_view volume = "volume";
inline constexpr std::string_view level = "level";
inline constexpr std::string_view mono = "mono";
inline constexpr std::string_view pause = "pause";
inline constexpr std::string_view stop = "stop";
inline constexpr std::string_view on = "1";
}

namespace ffmpeg {
inline constexpr std::string_view log_level = "-loglevel";
inline constexpr std::string_view error = "error";
inline constexpr std::string_view hide_banner = "-hide_banner";
inline constexpr std::string_view select_streams = "-select_streams";
inline constexpr std::string_view first_audio = "a:0";
inline constexpr std::string_view show_entries = "-show_entries";
inline constexpr std::string_view rate_channels = "stream=sample_rate,channels";
inline constexpr std::string_view output_format = "-of";
inline constexpr std::string_view bare_csv = "csv=p=0";
inline constexpr std::string_view input = "-i";
inline constexpr std::string_view force_format = "-f";
}

namespace ext {
inline constexpr std::string_view bin = ".bin";
inline constexpr std::string_view flac = ".flac";
inline constexpr std::string_view pyramid = ".pyramid";
inline constexpr std::string_view wav = ".wav";
inline constexpr std::string_view mp3 = ".mp3";
inline constexpr std::string_view m4a = ".m4a";
inline constexpr std::string_view aif = ".aif";
inline constexpr std::string_view aiff = ".aiff";
inline constexpr std::string_view ogg = ".ogg";
}

namespace bit_depth {
inline int decoded = 0;

template <typename T>consteval std::string_view of();
template <>consteval std::string_view of<int8_t>()  {return  "s8";}
template <>consteval std::string_view of<uint8_t>() {return  "u8";}
template <>consteval std::string_view of<int16_t>() {return "s16";}
template <>consteval std::string_view of<int32_t>() {return "s32";}
template <>consteval std::string_view of<float>()   {return "f32";}
template <>consteval std::string_view of<double>()  {return "f64";}

template <typename T>consteval std::string_view raw_of();
template <>consteval std::string_view raw_of<int8_t>()  {return    "s8";}
template <>consteval std::string_view raw_of<uint8_t>() {return    "u8";}
template <>consteval std::string_view raw_of<int16_t>() {return "s16le";}
template <>consteval std::string_view raw_of<int32_t>() {return "s32le";}
template <>consteval std::string_view raw_of<float>()   {return "f32le";}
template <>consteval std::string_view raw_of<double>()  {return "f64le";}
}

namespace pcm {
using shifted = int32_t;
inline std::vector<shifted> interleaved;
inline std::vector<std::byte> data;
inline size_t pos = 0;
inline size_t count = 0;
inline size_t sample = 0;
inline double ms = 0;
inline float hertz = 0;

template <typename Y>
inline std::vector<std::vector<Y>> pcm_db;
}  // namespace pcm

struct ProcessArg {
  const std::string& binary;
  std::vector<std::string_view> option;

  std::vector<std::string> words(std::initializer_list<std::string_view> tail) const {
    std::vector<std::string> word{binary};
    for (std::string_view each : option) word.emplace_back(each);
    for (std::string_view each : tail) word.emplace_back(each);
    return word;
  }
};

namespace pw_cat {
template <typename T>
std::vector<std::string> words(const ProcessArg& profile) {
  const std::string sample = std::to_string(pcm::sample);
  const std::string count  = std::to_string(pcm::count);
  return profile.words({format, bit_depth::of<T>(), rate, sample, channels, count, pipe_end});
}
}

namespace ffmpeg {
template <typename T>
std::vector<std::string> words(const ProcessArg& profile, std::string_view file) {
  return profile.words({file, force_format, bit_depth::raw_of<T>(), pipe_end});
}
}

inline const ProcessArg playback_passthrough{
  .binary = value.pw_cat,
  .option = {pw_cat::playback, pw_cat::raw}
};
inline const ProcessArg stream_passthrough{
  .binary = value.pw_cat,
  .option = {pw_cat::record, pw_cat::raw}
};
inline const ProcessArg probe_stream{
  .binary = value.ffprobe,
  .option = {ffmpeg::log_level, ffmpeg::error, ffmpeg::select_streams, ffmpeg::first_audio, ffmpeg::show_entries, ffmpeg::rate_channels, ffmpeg::output_format, ffmpeg::bare_csv}
};
inline const ProcessArg transcode_raw{
  .binary = value.ffmpeg,
  .option = {ffmpeg::hide_banner, ffmpeg::log_level, ffmpeg::error, ffmpeg::input}
};
inline constexpr std::array<std::string_view, 7> audioext_passthrough = {ext::wav, ext::flac, ext::mp3, ext::m4a, ext::aif, ext::aiff, ext::ogg};
}  // namespace env

inline std::string path(int64_t source_id, std::string_view extension) {
  return env::value.cache_dir + "/" + std::to_string(source_id) + std::string(extension);
}

std::vector<char*> argv(std::vector<std::string>& word);
pid_t spawn(std::vector<std::string>& word, int child_end, int child_fd, int parent_end);
std::expected<std::vector<std::byte>, int> capture(std::vector<std::string>& word);
std::expected<std::pair<int, int>, int> probe(const std::filesystem::path& file);
std::expected<void, int> transcode(const std::filesystem::path& file);
bool audio(const std::filesystem::path& file);
int python(std::string_view script, std::initializer_list<std::string_view> args);
inline int studio() { return python(env::module, {env::value.studio}); }

}
