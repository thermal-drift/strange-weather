# Changelog

## [2.0.1] - 2024-12-14

### Fixed
- **Fixed CV clipping** — Removed 3x output gain that caused square/clipped patterns at slow rates
- **Trail scales with range** — Display now samples at 10/20/60 fps for low/med/high ranges, so slow attractors show full trails
- **Fixed module tags** — Removed incorrect "Polyphonic" tag, now correctly shows "LFO" and "Function Generator"

## [2.0.0] - 2024-12-14

### Changed
- **Replaced Lorenz attractor with Sprott B** — The Lorenz attractor had issues crossing the separatrix between its two lobes, often getting stuck on one side. Sprott B is a minimal chaotic system that produces robust, reliable chaos without separatrix problems.
- **Replaced Halvorsen attractor with Dadras** — Halvorsen was unstable and produced poor visual output. Dadras is a multi-wing attractor with complex, interesting dynamics.
- **Improved CV output scaling** — Better normalization to utilize the full voltage range.

### Fixed
- Fixed bounding box decay that caused clipping artifacts (square patterns instead of smooth curves)
- Fixed warmup phase using different equations than runtime
- Added proper blow-up guards to prevent numerical instability crashes
- Improved attractor stability with tighter integration parameters

### Technical
- Bounds now only expand, never shrink (prevents clipping when attractor swings to extremes)
- Added warmup perturbation to avoid periodic lock-in
- Tightened blow-up threshold from 1e6 to 1000 for earlier detection of runaway states
