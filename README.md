# Strange Weather

A chaotic CV generator for VCV Rack 2 based on strange attractors.

## Overview

Strange Weather generates continuously evolving control voltages using mathematical chaos. Three independent attractor banks produce deterministic but unpredictable modulation — not random, but impossible to predict.

## Features

- **3 Independent Banks** — Each runs its own strange attractor
- **4 Attractor Types** — Lorenz, Rossler, Thomas, Halvorsen
- **16 CV Outputs** — 4 per bank + 4 combined
- **Real-time Visualization** — Watch the attractors evolve
- **3 Display Modes** — Trace (lines), Lissajous (phosphor dots), Scope (waveforms)
- **3D Display Mode** — See the full three-dimensional structure with rotation
- **Adjustable Trail Length** — Control how much history is displayed
- **Extreme Rate Range** — From 20-minute cycles to sub-second modulation

## Installation

### From Release

1. Download the latest release for your platform
2. Extract the zip file
3. Copy the `StrangeWeather` folder to your VCV Rack 2 plugins directory:

| Platform | Plugins Directory |
|----------|-------------------|
| macOS (Apple Silicon) | `~/Library/Application Support/Rack2/plugins-mac-arm64/` |
| macOS (Intel) | `~/Documents/Rack2/plugins/` |
| Windows | `%USERPROFILE%\Documents\Rack2\plugins\` |
| Linux | `~/.Rack2/plugins/` |

4. Restart VCV Rack

### Building from Source

Requires [VCV Rack SDK](https://vcvrack.com/manual/Building#Building-Rack-plugins) v2.x.

```bash
git clone https://github.com/kvarnelis/strange-weather.git
cd strange-weather
make RACK_DIR=/path/to/Rack-SDK
```

To install after building:

```bash
# macOS (Apple Silicon)
mkdir -p ~/Library/Application\ Support/Rack2/plugins-mac-arm64/StrangeWeather
cp plugin.dylib plugin.json ~/Library/Application\ Support/Rack2/plugins-mac-arm64/StrangeWeather/
cp -r res ~/Library/Application\ Support/Rack2/plugins-mac-arm64/StrangeWeather/

# macOS (Intel) / Linux / Windows - adjust path accordingly
```

## Controls

### Per Bank (A, B, C)

| Control | Function |
|---------|----------|
| **RATE** | Attractor evolution speed within selected range |
| **RNG** | Range selector: Low (5-20 min), Med (1s-2min), High (0.1-10s) |
| **SHAPE** | Attractor type: Lorenz, Rossler, Thomas, Halvorsen |
| **VOLT** | Output voltage: +/-5V, +/-10V, 0-5V, 0-10V |
| **CHAOS** | Primary chaos parameter - affects attractor behavior |

### Outputs Per Bank
- **x, y, z** — Attractor coordinates (scaled per VOLT setting)
- **SUM** — x + y + z

### Combined Outputs
- **SUM** — Sum of all bank sums
- **RECT** — Rectified (absolute values)
- **INV** — Inverted sum
- **DIST** — Inverse distance from center

### Display Controls
- **CYCLE** — Cycles through views: A, B, C, Combined, All
- **MODE** — Cycles display style: Trace (lines), Lissajous (tiny dots), Scope (time-based waveforms)
- **3D** — Toggles 3D rotation view
- **TRAIL** — Adjusts trail length from ~1 second to ~34 seconds of history

## Attractor Types

| Type | Character |
|------|-----------|
| **Lorenz** | The classic butterfly. Smooth, two-lobed with fold-over behavior |
| **Rossler** | Asymmetric spiral with occasional large excursions |
| **Thomas** | Cyclically symmetric, smooth rolling motion |
| **Halvorsen** | Sculptural and aggressive with sharp transitions |

## License

GPL-3.0-or-later

## Acknowledgments

- Edward Lorenz — for discovering deterministic chaos
- Worng Electronics Vector Space — output transform inspiration
- Synthesis Technology E352 — interface reference
