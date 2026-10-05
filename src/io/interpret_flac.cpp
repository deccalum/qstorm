#include <FLAC/stream_decoder.h>
#include <FLAC/stream_encoder.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>
#include <string>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "io/bridge.h"

namespace qstorm::io {

void flac_error(const FLAC__StreamDecoder*, FLAC__StreamDecoderErrorStatus, void*) {}

FLAC__StreamDecoderReadStatus
flac_read(const FLAC__StreamDecoder*, FLAC__byte buffer[], size_t* bytes, void*) {
  const size_t extent = std::min(*bytes, env::pcm::data.size() - env::pcm::pos);
  if (extent == 0) {
    *bytes = 0;
    return FLAC__STREAM_DECODER_READ_STATUS_END_OF_STREAM;
  }
  std::memcpy(buffer, env::pcm::data.data() + env::pcm::pos, extent);
  env::pcm::pos += extent;
  *bytes = extent;
  return FLAC__STREAM_DECODER_READ_STATUS_CONTINUE;
}

FLAC__StreamDecoderWriteStatus
flac_write(const FLAC__StreamDecoder*, const FLAC__Frame* frame, const FLAC__int32* const buffer[], void*) {
  const unsigned channels = frame->header.channels;
  const unsigned extent = frame->header.blocksize;
  const int shift = algorithm::shift<env::pcm::shifted>(env::bit_depth::decoded);
  for (unsigned i = 0; i < extent; ++i)
    for (unsigned c = 0; c < channels; ++c) env::pcm::interleaved.push_back(buffer[c][i] << shift);
  return FLAC__STREAM_DECODER_WRITE_STATUS_CONTINUE;
}

void flac_meta(const FLAC__StreamDecoder*, const FLAC__StreamMetadata* m, void*) {
  const int bps = static_cast<int>(m->data.stream_info.bits_per_sample);
  env::pcm::count = m->data.stream_info.channels;
  env::pcm::sample = m->data.stream_info.sample_rate;
  env::pcm::ms = static_cast<double>(env::pcm::sample) / 1000;
  env::bit_depth::decoded = bps;
  env::pcm::interleaved.reserve(m->data.stream_info.total_samples * env::pcm::count);
}

FLAC__StreamEncoderWriteStatus
flac_encoded(const FLAC__StreamEncoder*, const FLAC__byte buffer[], size_t bytes, uint32_t, uint32_t, void*) {
  const std::byte* p = reinterpret_cast<const std::byte*>(buffer);
  env::pcm::data.insert(env::pcm::data.end(), p, p + bytes);
  return FLAC__STREAM_ENCODER_WRITE_STATUS_OK;
}

void decode_flac(std::span<const std::byte> bytes) {
  env::pcm::data.assign(bytes.begin(), bytes.end());
  env::pcm::pos = 0;
  env::pcm::interleaved.clear();
  FLAC__StreamDecoder* decoder = FLAC__stream_decoder_new();
  FLAC__stream_decoder_init_stream(decoder, flac_read, nullptr, nullptr, nullptr, nullptr, flac_write, flac_meta, flac_error, nullptr);
  FLAC__stream_decoder_process_until_end_of_stream(decoder);
  FLAC__stream_decoder_delete(decoder);
  env::pcm::data.clear();
  env::pcm::data.shrink_to_fit();
}

void decode_flac_range(const std::string& path, double start, double end) {
  env::pcm::interleaved.clear();
  FLAC__StreamDecoder* decoder = FLAC__stream_decoder_new();
  FLAC__stream_decoder_init_file(decoder, path.c_str(), flac_write, flac_meta, flac_error, nullptr);
  FLAC__stream_decoder_process_until_end_of_metadata(decoder);
  const uint64_t first = static_cast<uint64_t>(std::llround(start*env::pcm::ms));
  const size_t frames = static_cast<size_t>(std::llround((end - start) * env::pcm::ms));
  FLAC__stream_decoder_seek_absolute(decoder, first);
  while (env::pcm::interleaved.size() < frames * env::pcm::count && FLAC__stream_decoder_get_state(decoder) != FLAC__STREAM_DECODER_END_OF_STREAM)
    FLAC__stream_decoder_process_single(decoder);
  env::pcm::interleaved.resize(frames * env::pcm::count);
  FLAC__stream_decoder_delete(decoder);
}

std::vector<std::byte> encode_flac(std::span<const env::pcm::shifted> interleaved, size_t count, int hz) {
  env::pcm::data.clear();
  FLAC__StreamEncoder* encoder = FLAC__stream_encoder_new();
  FLAC__stream_encoder_set_channels(encoder, static_cast<uint32_t>(count));
  FLAC__stream_encoder_set_bits_per_sample(encoder, static_cast<uint32_t>(env::ingest_bit));
  FLAC__stream_encoder_set_sample_rate(encoder, static_cast<uint32_t>(hz));
  FLAC__stream_encoder_init_stream(encoder, flac_encoded, nullptr, nullptr, nullptr, nullptr);
  const size_t extent = interleaved.size();
  const int shift = algorithm::shift<env::pcm::shifted>(env::ingest_bit);
  std::vector<FLAC__int32> sample(extent);
  for (size_t i = 0; i < extent; ++i) sample[i] = interleaved[i] >> shift;
  FLAC__stream_encoder_process_interleaved(encoder, sample.data(), static_cast<uint32_t>(extent / count));
  FLAC__stream_encoder_finish(encoder);
  FLAC__stream_encoder_delete(encoder);
  return std::move(env::pcm::data);
}

}
