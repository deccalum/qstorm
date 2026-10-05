# Quantum Stormsynth
Deconstructs audio into *"atoms"* to provide the building blocks for a natural
soundscape synth.

Personal project under active development. C++26 core using static
reflection, PostgreSQL storage, and a temporary Python web front end for
inspection.

## Status

**1/3 Foundation**  
`algorithm/`: low-level DSP, statistics and simulation primitives.  
`atom/`: window extractor working from the signal envelope alone.

- [`axis_quantizer.h`](src/algorithm/axis_quantizer.h): lets the foundation's
  math be written in ms / dB / Hz while frame indexing stays in the quantizer
- [`atom_thunder.h`](src/atom/atom_thunder.h): the root atom extractor, each
  stored extraction hold the git blob hashes so results trace back to the code
  that produced it.
- [`atom_constant.h`](src/atom/atom_constant.h): attempts to pre-define the
  atoms part 2 will define.

**2/3 Quantamization** (in progress)  
Defines the characteristics of a root atom. Non-negative matrix factorization
and cross-window component tracking are implemented and not yet wired into
the pipeline.

**3/3 Synthesis** (planned)
