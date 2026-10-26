#pragma once

#include <sys/types.h>
#include <unistd.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <filesystem>
#include <initializer_list>
#include <optional>
#include <span>
#include <string>
#include <string_view>

#include "algorithm/math_constant.h"

namespace qstorm {

namespace env {
namespace var {
inline constexpr const char* python    = "PYTHON";
inline constexpr const char* cache_dir = "CACHE_DIR";
inline constexpr const char* ffmpeg    = "FFMPEG";
inline constexpr const char* ffprobe   = "FFPROBE";
inline constexpr const char* ffplay    = "FFPLAY";
inline constexpr const char* studio    = "STUDIO";
inline constexpr const char* db_url    = "DATABASE_URL";
}

void        load(const std::filesystem::path& dotenv = ".env");
std::string require(const char* key);

inline constexpr float  play_block /*                                           */       = 100;
inline constexpr size_t gain_step /* steps inside 0-100 vol                     */       = 500;
inline constexpr int    ingest_bit /* hardcoded til qtz is tested to handle bit depth */ = 24;

inline constexpr std::string_view pipe_end  = "-";
inline constexpr std::string_view module    = "-m";
inline constexpr std::string_view freesound = "https://freesound.org/s/";
inline constexpr std::string_view described = "__";

namespace control {
inline constexpr std::string_view master = "master";
inline constexpr std::string_view volume = "volume";
inline constexpr std::string_view level  = "level";
inline constexpr std::string_view mono   = "mono";
inline constexpr std::string_view pause  = "pause";
inline constexpr std::string_view stop   = "stop";
inline constexpr std::string_view on     = "1";
}

namespace ffmpeg {
inline constexpr std::string_view log_level      = "-loglevel";
inline constexpr std::string_view error          = "error";
inline constexpr std::string_view hide_banner    = "-hide_banner";
inline constexpr std::string_view select_streams = "-select_streams";
inline constexpr std::string_view first_audio    = "a:0";
inline constexpr std::string_view show_entries   = "-show_entries";
inline constexpr std::string_view stream_section = "stream=";
inline constexpr std::string_view fmt_section    = ":format=";
inline constexpr std::string_view entry_sep      = ",";
inline constexpr std::string_view rate           = "sample_rate";
inline constexpr std::string_view channels       = "channels";
inline constexpr std::string_view codec_name     = "codec_name";
inline constexpr std::string_view sample_fmt     = "sample_fmt";
inline constexpr std::string_view bps            = "bits_per_sample";
inline constexpr std::string_view bprs           = "bits_per_raw_sample";
inline constexpr std::string_view bit_rate       = "bit_rate";
inline constexpr std::string_view ch_layout      = "channel_layout";
inline constexpr std::string_view fmt_name       = "format_name";
inline constexpr std::string_view duration       = "duration";
inline constexpr std::string_view size           = "size";
inline constexpr std::string_view output_fmt     = "-of";
inline constexpr std::string_view flat           = "flat";
inline constexpr std::string_view stream_prefix  = "streams.stream.0.";
inline constexpr std::string_view fmt_prefix     = "format.";
inline constexpr std::string_view absent         = "N/A";
inline constexpr std::string_view input          = "-i";
inline constexpr std::string_view force_fmt      = "-f";
}

namespace ffplay {
inline constexpr std::string_view nodisp    = "-nodisp";
inline constexpr std::string_view autoexit  = "-autoexit";
inline constexpr std::string_view nostats   = "-nostats";
inline constexpr std::string_view rate      = "-sample_rate";
inline constexpr std::string_view ch_layout = "-ch_layout";
inline constexpr std::string_view ch_count  = "c";
}

namespace ext {
inline constexpr std::string_view bin      = ".bin";
inline constexpr std::string_view flac     = ".flac";
inline constexpr std::string_view waveform = ".waveform";
inline constexpr std::string_view wav      = ".wav";
inline constexpr std::string_view mp3      = ".mp3";
inline constexpr std::string_view m4a      = ".m4a";
inline constexpr std::string_view aif      = ".aif";
inline constexpr std::string_view aiff     = ".aiff";
inline constexpr std::string_view ogg      = ".ogg";
}

namespace bit_depth {  // clang-format off
template <typename T>
consteval std::string_view raw_of();
template <> consteval std::string_view raw_of<int8_t>()  { return "s8";    }
template <> consteval std::string_view raw_of<uint8_t>() { return "u8";    }
template <> consteval std::string_view raw_of<int16_t>() { return "s16le"; }
template <> consteval std::string_view raw_of<int32_t>() { return "s32le"; }
template <> consteval std::string_view raw_of<float>()   { return "f32le"; }
template <> consteval std::string_view raw_of<double>()  { return "f64le"; }
}  // clang-format on

namespace pcm { using shifted = int32_t; }

struct ProcessArg {
  const char*           bin_var;
  vec<std::string_view> opt;

  vec<std::string> words(std::initializer_list<std::string_view> tail) const {
    vec<std::string> word{require(bin_var)};
    for (std::string_view each : opt) word.emplace_back(each);
    for (std::string_view each : tail) word.emplace_back(each);
    return word;
  }
};

namespace ffplay {
template <typename T>
vec<std::string> words(const ProcessArg& profile, double hz, size_t nc) {
  const std::string sample = std::to_string(std::lround(hz * 1000.0));
  const std::string layout = std::to_string(nc) + std::string(ch_count);
  return profile.words({ffmpeg::force_fmt, bit_depth::raw_of<T>(), rate, sample, ch_layout, layout,
                        ffmpeg::input, pipe_end});
}
}

namespace ffmpeg {
template <typename T>
vec<std::string> words(const ProcessArg& profile, std::string_view file) {
  return profile.words({file, force_fmt, bit_depth::raw_of<T>(), pipe_end});
}
}

inline const ProcessArg playback_raw{
  .bin_var = var::ffplay,
  .opt = {ffmpeg::hide_banner, ffmpeg::log_level, ffmpeg::error, ffplay::nostats, ffplay::nodisp,
          ffplay::autoexit}
};
inline const ProcessArg probe_stream{
  .bin_var = var::ffprobe,
  .opt     = {ffmpeg::log_level, ffmpeg::error, ffmpeg::select_streams, ffmpeg::first_audio,
              ffmpeg::output_fmt, ffmpeg::flat}
};

inline constexpr std::array<std::string_view, 8> stream_entry
  = {ffmpeg::rate, ffmpeg::channels, ffmpeg::codec_name, ffmpeg::sample_fmt,
     ffmpeg::bps,  ffmpeg::bprs,     ffmpeg::bit_rate,   ffmpeg::ch_layout};
inline constexpr std::array<std::string_view, 4> fmt_entry
  = {ffmpeg::fmt_name, ffmpeg::duration, ffmpeg::size, ffmpeg::bit_rate};
inline const ProcessArg transcode_raw{
  .bin_var = var::ffmpeg,
  .opt     = {ffmpeg::hide_banner, ffmpeg::log_level, ffmpeg::error, ffmpeg::input}
};
inline constexpr std::array<std::string_view, 7> audioext_passthrough
  = {ext::wav, ext::flac, ext::mp3, ext::m4a, ext::aif, ext::aiff, ext::ogg};
}

using pcmvec  = vec<env::pcm::shifted>;
using pcmspan = std::span<const env::pcm::shifted>;

namespace io {
struct Original {
  std::string                container;
  std::string                codec;
  std::string                sample_fmt;
  std::optional<int>         bit_depth;
  std::optional<int64_t>     bit_rate;
  std::optional<std::string> ch_layout;
  double                     length;
  int64_t                    size;
};

struct Probe {
  double   hz;
  int      channels;
  Original original;
};

vec<char*> argv(vec<std::string>& word);
pid_t      spawn(vec<std::string>& word, int child_end, int child_fd, int parent_end);
std::expected<vec<std::byte>, int> capture(vec<std::string>& word);

std::expected<Probe, int>  probe(const std::filesystem::path& file);
std::expected<pcmvec, int> transcode(const std::filesystem::path& file);
bool                       audio(const std::filesystem::path& file);

std::expected<int, int>
python(std::string_view script, std::initializer_list<std::string_view> args);

inline std::expected<int, int> studio() {
  return python(env::module, {env::require(env::var::studio)});
}

inline std::string path(int64_t source_id, std::string_view ext) {
  return env::require(env::var::cache_dir) + "/" + std::to_string(source_id) + std::string(ext);
}
}  // namespace io
}
