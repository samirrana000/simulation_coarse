# Scorer — Custom JSON MLP (I87)

*Source anchors: `src/scorer.js:35` `PoseScorer`, `src/ml-tier.js`.*

This note documents the **JSON MLP scorer** that the browser actually ships.

> **ONNX is NOT implemented and is out of scope.** There is no `OnnxScorer`, no
> `onnxruntime-web` dependency, and no `.onnx` loader anywhere in `src/`. The
> `src/scorer-onnx.js` stub that used to sit here has been deleted — it only
> printed `console.warn("ONNX scorer not yet implemented")` and delegated to
> `PoseScorer`, so it advertised a capability that did not exist. Do not write
> code against it. See "ONNX — absent" below.

> Fidelity note: `src/mol2.js` preserves bond order for topology only — `grep -n "bond order.*topology" src/mol2.js` hits header.

## Current: custom JSON MLP (`src/scorer.js`)

The browser ships a dependency-free feed-forward MLP:

```js
import { PoseScorer } from "./src/scorer.js";
const scorer = PoseScorer.docked();          // built-in linear docked prior
// or loaded from JSON:
const cfg = await fetch("weights.json").then(r=>r.json());
const scorer = new PoseScorer(cfg);
const score = scorer.predict(feat6);         // feat6 = [nContactsNorm, burialNorm, clashNorm, ligRmsdNorm, cvNorm, coulombNorm]
```

Weights JSON schema (`src/scorer.js:38`):

```json
{
  "layers":   [6, 4, 1],
  "W":        [[/* 4×6 row-major */], [/* 1×4 */]],
  "b":        [[/* 4 biases */], [0]],
  "act":      "silu",
  "featMean": [0,0,0,0,0,0],
  "featStd":  [1,1,1,1,1,1]
}
```

- `layers` — `[n_in, n_h1, ..., 1]` (last must be 1)
- `W[L]` row-major per layer, `b[L]` bias
- `act` ∈ `silu|relu|tanh|none` (hidden layers only)
- Features are 6 normalized scalars (see `src/scorer.js:16` `POSE_FEATURE_N=6`)

This format is tiny (<1 KB), requires no runtime beyond `src/scorer.js`, and is the **default contract** used by `src/ml-tier.js` tickbox `NN pose score (MLP)`. Any JSON that validates against `PoseScorer` constructor runs.

**Exporting a new JSON scorer from Python:**

```python
import json, numpy as np
cfg = {
  "layers": [6, 8, 1],
  "W": [np.random.randn(8*6).tolist(), np.random.randn(1*8).tolist()],
  "b": [np.random.randn(8).tolist(), [0.0]],
  "act": "silu",
  "featMean": [0]*6, "featStd": [1]*6
}
json.dump(cfg, open("scorer.json","w"), indent=2)
# Browser: PoseScorer(JSON)
```

## ONNX — absent (out of scope)

There is **no ONNX support in this project**. Concretely:

- No `OnnxScorer` class, no `scorer-onnx.js` module — the stub was deleted
  2026-09-30 (see `CHANGELOG.md` → Unreleased).
- No `onnxruntime-web` / `ort` dependency. The project is zero-dependency by
  design (`package.json`); an ONNX path would break that guarantee.
- No `.onnx` file is fetched, and `OnnxScorer.fromOnnx(...)` does not exist.

**What to use instead:** the JSON `PoseScorer` in `src/scorer.js`. It takes the
identical 6-feature vector and the identical `predict(feat) → number` contract,
so a PyTorch model exported as JSON weights is a drop-in replacement with no
code change beyond supplying the weights file:

```js
import { PoseScorer } from "./src/scorer.js";
const scorer = new PoseScorer(await fetch("weights.json").then(r => r.json()));
```

**Why the feature is not on the roadmap:** the score is a cosmetic readout, not
a force-field term. The stub was imported by zero files — it never warned on a
real page load either — which is exactly the problem: it was an unreachable file
whose only observable behaviour was a warning message, kept alive by prose in
this document. Unsurfaced work has zero perceived value (wiki P3), so it was
deleted rather than advertised. If a real ONNX runtime ever lands it must arrive
as a working implementation with a parity test against `PoseScorer`, not as a
warning message.

**If you need ONNX in your own pipeline:** export the PyTorch MLP to JSON
weights as shown above, or run the ONNX graph in Python and hand the resulting
6-parameter input to `PoseScorer` in the browser.

## Grep & measurability

```bash
grep -n "PoseScorer" src/scorer.js src/ml-tier.js   # hits the real scorer and its wiring
grep -n "bond order.*topology" src/mol2.js docs/*  # hits src/mol2.js header
```

## References

- `src/scorer.js:101` `PoseScorer.docked()` built-in linear prior
- `src/ml-tier.js:27` JSON weights loader wired to the "NN pose score" tickbox
- `src/main.js:1571` the live `state.scorer.predict(feat)` call site

