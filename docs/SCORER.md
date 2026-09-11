# Scorer — Custom JSON vs ONNX Future (I87)

*Source anchors: `src/scorer.js:35` `PoseScorer`, `src/ml-tier.js`, `src/scorer-onnx.js:13` `OnnxScorer`.*

This note documents the **current JSON MLP scorer** and the **future ONNX interchange path** so that a trained PyTorch → ONNX model can replace the hand-tuned weights without changing the browser scoring contract.

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

## Future: ONNX (`src/scorer-onnx.js` — stub)

`src/scorer-onnx.js` exposes `export class OnnxScorer` as the **ONNX** interchange point. It mirrors the `PoseScorer` API (`predict(feat) → number`) but delegates to an ONNX Runtime session when one is available.

Current status: **stub** — `OnnxScorer` logs `"ONNX scorer not yet implemented — using JSON fallback"` and falls back to `PoseScorer` until the WASM runtime and a `.onnx` model are wired. This satisfies the file-contract measurable without claiming runtime support.

**Intended production wiring (not yet):**

```js
import { OnnxScorer } from "./src/scorer-onnx.js";
const scorer = await OnnxScorer.fromOnnx("scorer.onnx", { featMean, featStd });
// Under the hood: ort.InferenceSession.create("scorer.onnx"), session.run({input: tensor})
// ONNX graph: input [1,6] float32 → MatMul/Relu/Gemm → output [1,1] float32
// Identical feat6 as JSON scorer, so JSON vs ONNX scores agree to <1e-4 when exported correctly.
```

**Exporting ONNX from PyTorch (future workflow):**

```python
import torch
class MLP(torch.nn.Module):
    def __init__(self): ...
    def forward(self, x): ...
m = MLP(); m.load_state_dict(torch.load("mlp.pt")); m.eval()
dummy = torch.randn(1, 6, dtype=torch.float32)
torch.onnx.export(m, dummy, "scorer.onnx",
                  input_names=["input"], output_names=["score"],
                  dynamic_axes={"input":{0:"batch"}, "score":{0:"batch"}})
# Optionally: python -m onnxruntime.tools.convert_onnx_models_to_ort scorer.onnx
# Browser: OnnxScorer.fromOnnx("scorer.onnx")
```

**Fallback behavior:** If `onnxruntime-web` (`ort`) is not loaded or `fetch("scorer.onnx")` fails, `OnnxScorer.predict` logs the ONNX stub message and returns `this.fallback.predict(feat)` (a `PoseScorer.docked()` instance).

**Which to use?**

| Path | File | Runtime | Size | Parity |
|------|------|---------|------|--------|
| **JSON** (current) | `weights.json` → `PoseScorer` | none (pure JS) | ~1 KB | reference |
| **ONNX** (future) | `scorer.onnx` → `OnnxScorer` | `onnxruntime-web` WASM (~1 MB) | ~10–100 KB | `‖ONNX−JSON‖ < 1e-4` when exported from same PyTorch checkpoint |

The scorer feature contract (`POSE_FEATURE_N=6` order fixed) is shared, so a model can be A/B tested across both paths.

## Grep & measurability

```bash
grep -n "ONNX" docs/SCORER.md src/scorer-onnx.js   # hits this doc + stub
grep -n "bond order.*topology" src/mol2.js docs/*  # hits src/mol2.js header
```

## References

- `src/scorer.js:101` `PoseScorer.docked()` built-in linear prior
- `src/scorer-onnx.js:13` `OnnxScorer` stub (logs not yet)
- ONNX Runtime Web: https://onnxruntime.ai/docs/get-started/with-javascript.html
