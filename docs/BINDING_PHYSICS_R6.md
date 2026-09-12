# Phase R6 — Efficient Simulation-Data Capture Port (BindLog)

Loop-1, step 6. Implementation: `src/capture/bindlog.js`. Tests: `scripts/test_bindlog.mjs` (18/18).
Additive only — no existing module touched; `tests/test_all.js` stays 32/32, gate OPEN.

## 1. Why

The recorder (`src/recorder.js`) stores raw Float32 frames (2.2 KB × 500) and nothing else.
R4 needs per-frame energy components + apo/holo frame pairs; R7 needs contact timelines,
PMF hills, state hops; downstream ML wants compact binary. BindLog adds both channels
without touching the existing hot path.

## 2. DENSE channel — quantized delta frames

- Frame 0: absolute Float32 (3n).
- Frames 1..N: Int16 per-coordinate deltas × quantStep (default 0.01 Å), each frame
  quantized against the *reconstructed* previous frame (bounded drift, no error accumulation
  beyond q/2 per coordinate).
- Measured: **1.95× compression** (30,700 B vs 60,000 B for 100 frames × 50 atoms),
  **max error 5.0e-3 Å** (≤ quantStep/2 ✓) on 0.5 Å thermal-motion synthetic data.
- Same-delta clipping at ±32767 steps (±327 Å) — unreachable in practice at 0.01 Å steps.

## 3. SPARSE channel — struct-of-arrays event log

Typed arrays, zero per-event objects, amortized ×2 growth:
`evTime Float64 | evType Uint8 | evA/evB Int32 | evX/evY Float32`

| Type id | Name | a | b | x | y |
|---|---|---|---|---|---|
| 0 | energy | term id (0-6 per R4) | −1 | value | 0 |
| 1 | contact+ | ligand atom | residue | dist Å | 0 |
| 2 | contact− | ligand atom | residue | dist Å | 0 |
| 3 | hill | 0 | −1 | CV | height |
| 4 | state | stateId 0-3 | −1 | 0 | 0 |
| 5 | pocketVol | 0 | −1 | Å³ | 0 |

Memory: 10k events = **800 KB** typed arrays (vs several MB with per-event objects at GC
pressure); helpers `pushEnergyComponents/pushContact/pushHill/pushState/pushPocketVolume`.

## 4. Binary blob format (BLG1 v1)

```
+0000  4B  magic "BLG1" (0x424c4731 LE)
+0004  4B  u32 version = 1
+0008  4B  u32 headerLen
+000C  headerLen  ASCII JSON header (offsets zero-padded to fixed width —
       length-stable rewrite-in-place): {n, nF, nE, quantStep, typeNames, offsets}
+pad   0-3B zero pad to 4-byte alignment
+0     3n×4B  Float32 frame0
+1     nF×4B  Float32 frame times (ps)
+2     3n×2B×(nF−1)  Int16 delta frames
+3     32B×nE  events: f64 time | u8 type (+3 pad) | i32 a | i32 b | f32 x | f32 y
end    4B "END1" + 4B total size
```

Round-trip verified bit-identical (max drift 0.0 on the real 4W52 run; times match to 1e-4).
JSON header offsets are relative to the sections start (after header pad).

## 5. Integration points (Loop-2 S4 wired — read-only phase superseded)

- `importRecorder(recorder)` — drains existing Recorder frames into the dense channel.
- `captureFrame(pos, timePs)` — call from any loop; mini-run integration demonstrated in
  the test with the real `ForceField` + `LangevinIntegrator` (40 frames, 120 energy events).
- **Wired (S4, `src/main.js`)**: `ff-binding.js`/`heavy.js` per-term accumulators
  (`ff.trackTerms`, default OFF) → `pushEnergyComponents` (7 terms, recorder stride,
  idempotent `_nextAt` scheduler); `funnel.js deposit()` → `Funnel.onHill` guarded
  callback → `pushHill`; contact form/break diffed from the nContacts HUD pair set
  (5.5 Å) via `state._lastContacts` Map → `pushContact`; `state.bindLog = new
  BindLog()` per Build System; UI checkbox `bindlogOn` (Recording panel, default
  OFF — hot path bit-identical). `physics/network.js` hop → `pushState` (S6+);
  `cryptic_pockets.trackPocketVolume` → `pushPocketVolume` (S6+).
- Integration test: `scripts/test_bindlog_integration.mjs` (14/14) — 4W52 CG
  trackTerms split sums exactly to bindingU, OFF bit-identical, 300-step run
  → 30 frames / 210 energy events / contact+/− events / blob 0-drift.

## 6. Space accounting (4W52, 500 frames, 10k events)

| Store | Bytes |
|---|---|
| Recorder raw Float32 (today) | 1.09 MB |
| BindLog dense (q=0.01) | **0.56 MB** (frame0 1.97 KB + deltas 488 KB + times 2 KB) |
| BindLog sparse 10k events | 0.80 MB in-RAM / 0.32 MB in blob |
| Total BindLog blob | ~0.9 MB for everything (frames + times + all event streams) |

## 7. Downstream consumers

- R4 entropy: dense frames holo+apo → covariance; energy events → ΔH decomposition.
- R7 visualization: contact+/− events → interaction timeline stripes; hills → PMF
  formation animation; states → binding-state track.
- Future ML training data: one ArrayBuffer per run, directly memmap-able via typed views.
