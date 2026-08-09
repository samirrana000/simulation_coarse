/**
 * scorer.js — small neural pose-scoring network (plain JS, no dependencies).
 *
 * GNN/MLP pose scorer for the binding demo (Step 4, ML tier): takes a fixed
 * feature vector describing the current protein–ligand pose and returns a
 * scalar "pose score" (higher ≈ more plausible docked pose). The class is a
 * generic feed-forward MLP so a real trained model (silu/relu/tanh layers,
 * weights loaded from a JSON config) can be dropped in unchanged; the bundled
 * default (`PoseScorer.docked()`) is a single-layer linear model whose
 * weights encode the intuitive physics of binding (more protein–ligand
 * contacts, more burial/desolvation gain, fewer clashes, lower ligand RMSD
 * and CV distance ⇒ better pose) — an honest stand-in until an offline-trained
 * network is exported (see ml/export_esm_contacts.py for the companion
 * contact-prior exporter).
 *
 * Feature vector (length POSE_FEATURE_N, order fixed):
 *   0  nContactsNorm = min(1, nContacts / 20)
 *   1  burialNorm    = clip(−desolvU / 3.3, 0, 1.2)     (EEF1-lite burial gain)
 *   2  clashNorm     = min(1, clashCount / 5)           (hard vdW overlaps)
 *   3  ligRmsdNorm   = min(1, ligRMSD / 4)
 *   4  cvNorm        = min(1, CV / 10)                  (ligand↔pocket COM dist)
 *   5  coulombNorm    = clip(−U_el / 2, 0, 1)            (electrostatic gain)
 *
 * All features are normalized to ~[0,1] before the network (featMean/featStd
 * are applied too — the docked defaults use 0/1 so the normalization above
 * is final), then the MLP returns ONE score.
 *
 * Performance: the network is tiny (≤ a few hundred FLOPs per call); it is
 * evaluated at most once per rendered frame — orders of magnitude below the
 * physics budget.
 */

export const POSE_FEATURE_N = 6;

export class PoseScorer {
  /**
   * @param {object} w JSON config:
   *   {
   *     layers:   [n_in, n_h1, ..., n_out],
   *     W:        [Float64Array-compatible flat arrays, out-major per layer],
   *     b:        [per-layer bias vectors],
   *     act:      "silu" | "relu" | "tanh" | "none"  (hidden-layer activation),
   *     featMean: [n_in], featStd: [n_in]  (optional; default 0/1)
   *   }
   *   W[L] has length layers[L+1]*layers[L], row-major: y[o] = Σ_i W[o*in + i] * x[i] + b[o].
   */
  constructor(w) {
    this.layers = w.layers;
    this.act = w.act ?? "silu";
    this.W = (w.W ?? []).map((a) => Float64Array.from(a));
    this.b = (w.b ?? []).map((a) => Float64Array.from(a));
    this.nIn = this.layers[0];
    const nOut = this.layers[this.layers.length - 1];
    if (nOut !== 1) throw new Error(`PoseScorer: last layer must be 1 (got ${nOut})`);
    const fm = w.featMean ?? new Array(this.nIn).fill(0);
    const fs = w.featStd ?? new Array(this.nIn).fill(1);
    this.featMean = Float64Array.from(fm);
    this.featStd = Float64Array.from(fs);
    // structural sanity
    for (let L = 0; L < this.layers.length - 1; L++) {
      const nin = this.layers[L], nout = this.layers[L + 1];
      if (this.W[L].length !== nin * nout) throw new Error(`PoseScorer: W[${L}] size mismatch (${this.W[L].length} != ${nin}*${nout})`);
      if (this.b[L].length !== nout) throw new Error(`PoseScorer: b[${L}] size mismatch`);
    }
  }

  /** activation: silu / relu / tanh / none */
  _act(x) {
    if (this.act === "silu") return x > 40 ? x : x / (1 + Math.exp(-x));
    if (this.act === "relu") return x > 0 ? x : 0;
    if (this.act === "tanh") return Math.tanh(x);
    return x;
  }

  /** @param {number[]} feat length POSE_FEATURE_N (normalized ~[0,1]) @returns {number} */
  predict(feat) {
    if (feat.length !== this.nIn) throw new Error(`PoseScorer: feature vector has ${feat.length} entries, expected ${this.nIn}`);
    let x = new Float64Array(this.nIn);
    for (let i = 0; i < this.nIn; i++) x[i] = (feat[i] - this.featMean[i]) / this.featStd[i];
    for (let L = 0; L < this.layers.length - 1; L++) {
      const nin = this.layers[L], nout = this.layers[L + 1];
      const W = this.W[L], b = this.b[L];
      const y = new Float64Array(nout);
      for (let o = 0; o < nout; o++) {
        let s = b[o];
        const row = o * nin;
        for (let i = 0; i < nin; i++) s += W[row + i] * x[i];
        y[o] = L === this.layers.length - 2 ? s : this._act(s); // no activation on output layer
      }
      x = y;
    }
    return x[0];
  }

  /**
   * Bundled default: a depth-1 LINEAR model with hand-set weights that encode
   * the intuitive physics of binding. features order per POSE_FEATURE_N above.
   * Rough scale: a well-docked benzene-in-cavity pose scores ≈ +3, an exposed
   * /clashing pose ≈ −2.
   */
  static docked() {
    return new PoseScorer({
      layers: [POSE_FEATURE_N, 1],
      act: "none",
      //   nContacts, burial, clash,  ligRMSD,  cvNorm, coulomb
      W: [[2.0,      3.0,   -4.0,   -2.0,    -1.0,   1.5]],
      b: [[-2.0]],
    });
  }
}
