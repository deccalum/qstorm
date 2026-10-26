#pragma once

#include <cstdint>
#include <qr/qr.hpp>

namespace qstorm::cli {

namespace library {
void ingest_source(const char* directory);
void decode_source(int64_t source_id);
}

namespace extract {
void root_atom(int64_t id, int64_t dna_id);
void atomize(int64_t root_id, int64_t index_id, int64_t ard_id, int ch);
void group_component(int64_t atomization_id);
void link_component(const char* domain, const char* source, const char* atomization);
void subtract_reconstruction(int64_t atomization_id);
}

namespace render {
void play_audio(int64_t audio, double offset);
void play_component_set(int64_t atomization_id, const char* component);
void play_residual(int64_t atomization_id);
void render_component_list(int64_t atomization_id, int64_t width);
void render_waveform(int64_t id, double x0, double x1, int64_t width);
void write_stft(int64_t root_id, int64_t index_id);
void render_stft(int64_t root_id, int64_t index_id);
void render_residual(int64_t atomization_id);
void render_reconstruction(int64_t atomization_id);
void render_component(int64_t atomization_id, int64_t component_index);
}

namespace query {
void read_table(const char* table, const char* parent);
void edit_root(int64_t run_id, int64_t index, const char* column, double value);
void restore_root(int64_t run_id, int64_t index);
void edit_dna(int64_t dna_id, const char* column, double value);
void copy_dna(int64_t dna_id);
void clear_cache(const char* table, int64_t parent);
void clear_root_extraction();
}

}
