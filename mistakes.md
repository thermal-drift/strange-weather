# Strange Weather Debugging: Mistakes Found

## Initial Problem
Four strange attractors (A=Lorenz, B=Rossler, C=Thomas, D=Halvorsen) were displaying incorrectly:
- A (Lorenz): Diagonal line instead of butterfly
- B (Rossler): Square instead of spiral
- C (Thomas): Square instead of smooth curves
- D (Halvorsen): Somewhat correct but clipped

---

## Mistake 1: Lorenz rho parameter below critical threshold

**Location:** `derivatives()` function, LORENZ case

**Problem:**
```cpp
const double rho = 24.0 + chaos * 8.0;  // range: 24-32
```
At chaos=0, rho=24, which is **below** the critical threshold (~24.74) required for chaotic behavior. The system converges to a fixed point instead of exhibiting the butterfly pattern.

**Fix:** Changed to `rho = 28.0 + chaos * 12.0` to ensure rho is always above the critical threshold.

---

## Mistake 2: Bounding box decay causing clipping

**Location:** `step()` function, bounding box update

**Problem:**
```cpp
const double decay = 0.9999;
minX = std::min(minX * decay + x * (1.0 - decay), x);
maxX = std::max(maxX * decay + x * (1.0 - decay), x);
```
The EMA decay caused bounds to **shrink** when the attractor stayed inside them. When the attractor then swung to extremes, values exceeded the (shrunk) bounds, got normalized to values outside ±1, and were clipped. This created square patterns.

**Fix:** Removed decay, bounds only expand:
```cpp
minX = std::min(minX, x);
maxX = std::max(maxX, x);
```

---

## Mistake 3: Initial bounding boxes too narrow

**Location:** `resetState()` function

**Problem:** Initial bounds were too small for the attractor ranges:
- Rossler: ±15 was insufficient for c=18 which can reach ±25-30
- Thomas: ±3 was too tight
- Halvorsen: ±10 was insufficient

**Fix:** Widened initial bounds significantly.

---

## Mistake 4: LORENZ warmup used different equations than runtime

**Location:** `resetState()` function

**Problem:** The LORENZ case had a special warmup block that ran **actual Lorenz equations**:
```cpp
if (type == LORENZ) {
    // Dedicated Lorenz warmup...
    double dx = sigma * (ly - lx);
    double dy = lx * (rho - lz) - ly;
    double dz = lx * ly - beta * lz;
    ...
}
```
But in `derivatives()`, LORENZ was running **Rossler equations**:
```cpp
case LORENZ: {
    // Rossler (duplicate) parameters for bank A
    dx = -y - z;
    dy = x + a * y;
    dz = b + z * (x - c);
}
```
The warmup left the state in Lorenz phase space (z~25-45), then Rossler equations took over from that completely wrong starting point.

**Fix:** Made LORENZ use the generic `warmup()` function that calls `step()` which uses `derivatives()`, so warmup and runtime use the same equations.

---

## Mistake 5: LORENZ had different typeScale than ROSSLER

**Location:** `process()` function

**Problem:** Even though both were running identical Rossler equations:
```cpp
case ROSSLER: typeScale = 1.5f; break;
case LORENZ: typeScale = 15.0f; break;  // 10x faster!
```

**Fix:** Set LORENZ typeScale to 1.5f to match ROSSLER.

---

## Mistake 6: LORENZ had different initial conditions than ROSSLER

**Location:** `resetState()` function

**Problem:** LORENZ started at fixed position (0.2, 0.1, 0.1) with preset bounds, while ROSSLER used random perturbations and bounds starting at current position.

**Fix:** Made LORENZ initialization identical to ROSSLER:
```cpp
case LORENZ:
    x = 0.1 + (random::uniform() - 0.5) * 0.02;
    y = 0.1 + (random::uniform() - 0.5) * 0.02;
    z = 0.0 + (random::uniform() - 0.5) * 0.02;
    minX = maxX = x;
    minY = maxY = y;
    minZ = maxZ = z;
    warmupSteps = 2000; warmupDt = 0.01;
```

---

## Mistake 7: Not installing the plugin after building

**Location:** Build/deployment process

**Problem:** Running `make` only builds `plugin.dylib` in the project directory. The plugin wasn't being copied to VCV Rack's plugin folder:
```
~/Library/Application Support/Rack2/plugins-mac-arm64/StrangeWeather/
```

**Fix:** Either run `make install` or manually copy:
```bash
cp plugin.dylib "~/Library/Application Support/Rack2/plugins-mac-arm64/StrangeWeather/"
```

---

## Key Lesson

When debugging why two things that should be identical are behaving differently, check **every** code path that touches them:
1. Initialization / initial conditions
2. Warmup / settling phase
3. Runtime equations
4. Rate/speed scaling
5. Normalization bounds
6. Build and deployment
