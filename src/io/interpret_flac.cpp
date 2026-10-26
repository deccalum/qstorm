#include <FLAC/stream_decoder.h>
#include <FLAC/stream_encoder.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <span>
#include <utility>
#include <vector>

#include "algorithm/axis_quantizer.h"
#include "io/bridge.h"
#include "io/interpret_file.h"

// clang-format off
namespace qstorm::io {
struct Decoding {
  std::span<const std::byte> bytes;
  size_t pos;
  int bit_depth;
  pcmvec xc;
};

void flac_error(const FLAC__StreamDecoder*, FLAC__StreamDecoderErrorStatus, void*) {}

FLAC__StreamDecoderReadStatus flac_read(const FLAC__StreamDecoder*, FLAC__byte buffer[], size_t* bytes, void* client) {
  Decoding& decoding = *static_cast<Decoding*>(client);
  const size_t n = std::min(*bytes, decoding.bytes.size() - decoding.pos);
  if (n == 0) {
    *bytes = 0;
    return FLAC__STREAM_DECODER_READ_STATUS_END_OF_STREAM;
  }
  std::memcpy(buffer, decoding.bytes.data() + decoding.pos, n);
  decoding.pos += n;
  *bytes = n;
  return FLAC__STREAM_DECODER_READ_STATUS_CONTINUE;
}

FLAC__bool flac_eof(const FLAC__StreamDecoder*, void* client) {
  const Decoding& decoding = *static_cast<const Decoding*>(client);
  return decoding.pos == decoding.bytes.size();
}

FLAC__StreamDecoderWriteStatus flac_write(const FLAC__StreamDecoder*, const FLAC__Frame* frame, const FLAC__int32* const buffer[], void* client) {
  Decoding& decoding = *static_cast<Decoding*>(client);
  const unsigned nc = frame->header.channels;
  const unsigned n = frame->header.blocksize;
  const int shift = algorithm::shift<env::pcm::shifted>(decoding.bit_depth);
  for (unsigned i = 0; i < n; ++i) for (unsigned c = 0; c < nc; ++c) decoding.xc.push_back(buffer[c][i] << shift);
  return FLAC__STREAM_DECODER_WRITE_STATUS_CONTINUE;
}

void flac_meta(const FLAC__StreamDecoder*, const FLAC__StreamMetadata* m, void* client) {
  Decoding& decoding = *static_cast<Decoding*>(client);
  decoding.xc.reserve(m->data.stream_info.total_samples * m->data.stream_info.channels);
}

FLAC__StreamEncoderWriteStatus flac_encoded(const FLAC__StreamEncoder*, const FLAC__byte buffer[], size_t bytes, uint32_t, uint32_t, void* client) {
  vec<std::byte>& encoded = *static_cast<vec<std::byte>*>(client);
  const std::byte* p = reinterpret_cast<const std::byte*>(buffer);
  encoded.insert(encoded.end(), p, p + bytes);
  return FLAC__STREAM_ENCODER_WRITE_STATUS_OK;
}

pcmvec flac_decode(std::span<const std::byte> bytes, int bit_depth) {
  Decoding decoding{bytes, 0, bit_depth, {}};
  FLAC__StreamDecoder* decoder = FLAC__stream_decoder_new();
  FLAC__stream_decoder_init_stream(decoder, flac_read, nullptr, nullptr, nullptr, flac_eof, flac_write, flac_meta, flac_error, &decoding);
  FLAC__stream_decoder_process_until_end_of_stream(decoder);
  FLAC__stream_decoder_delete(decoder);
  return std::move(decoding.xc);
}

vec<std::byte> flac_encode(pcmspan xc, size_t nc, double hz) {
  vec<std::byte> encoded;
  FLAC__StreamEncoder* encoder = FLAC__stream_encoder_new();
  FLAC__stream_encoder_set_channels(encoder, static_cast<uint32_t>(nc));
  FLAC__stream_encoder_set_bits_per_sample(encoder, static_cast<uint32_t>(env::ingest_bit));
  FLAC__stream_encoder_set_sample_rate(encoder, static_cast<uint32_t>(std::lround(hz*1000.0)));
  FLAC__stream_encoder_init_stream(encoder, flac_encoded, nullptr, nullptr, nullptr, &encoded);
  const size_t n = xc.size();
  const int shift = algorithm::shift<env::pcm::shifted>(env::ingest_bit);
  vec<FLAC__int32> sample(n);
  for (size_t i = 0; i < n; ++i) sample[i] = xc[i] >> shift;
  FLAC__stream_encoder_process_interleaved(encoder, sample.data(), static_cast<uint32_t>(n / nc));
  FLAC__stream_encoder_finish(encoder);
  FLAC__stream_encoder_delete(encoder);
  return encoded;
}
}
