/**
 * bindviz.js — Binding-physics visualization renderers (R7, Loop-1).
 *
 * Pure Canvas2D renderers consuming BindLog sparse events + dense frames.
 * All renderers follow the three-state rule (no-data text / steady / active),
 * are null-safe (missing canvas or events → draw actionable empty state),
 * anti-slop flat style: no gradients, no shadows, 1px rules, 10-12px mono.
 *
 * Interaction class colors (R1 catalog):
 *   hydrophobic/vdW  #FBBF24 (amber)
 *   H-bond           #38BDF8 (cyan)
 *   salt bridge +    #F87171 (red, acidic partner)
 *   π/aromatic       #E879F9 (magenta)
 *   unknown          #94A3B8 (slate)
 */

const CLS = {
  bg: "#0B0D0E", panel: "#121517", grid: "#2A3036",
  dim: "#4a545c", text: "#cbd5e1",
  hydro: "#FBBF24", hb: "#38BDF8", salt: "#F87171", arom: "#E879F9", unk: "#94A3B8",
};

/**
 * Effective CSS-pixel size of a canvas. The UI sizer (main.js bindvizFit)
 * sets canvas.width = cssW × dpr and tags canvas._dpr = dpr before rendering;
 * dividing device px back by _dpr keeps layout in CSS px under the dpr
 * transform. Headless/mocked canvases have no _dpr → 1 (identity).
 */
const dimsOf = (canvas) => {
  const d = canvas._dpr || 1;
  return { W: canvas.width / d, H: canvas.height / d };
};

/** Draw the empty state (three-state rule, state 1). */
function drawEmpty(canvas, msg) {
  if (!canvas || !canvas.getContext) return;
  const ctx = canvas.getContext("2d");
  if (!ctx) return;
  const { W, H } = dimsOf(canvas);
  const w = W || 300, h = H || 140;
  ctx.fillStyle = CLS.bg;
  ctx.fillRect(0, 0, w, h);
  ctx.fillStyle = CLS.dim;
  ctx.font = "11px ui-monospace, monospace";
  ctx.textAlign = "center";
  ctx.fillText(msg || "no binding data — record a trajectory with BindLog attached", w / 2, h / 2);
  ctx.textAlign = "left";
}

/**
 * (a) Interaction timeline — rows per (ligand-atom, residue) contact pair.
 * contact+ starts a colored run; contact− ends it. Row label = residue id.
 * Side list: top contacts by mean lifetime.
 *
 * @param {HTMLCanvasElement} canvas
 * @param {BindLog} bindlog
 * @param {{tFrom?:number, tTo?:number}} [range]
 */
export function renderInteractionTimeline(canvas, bindlog, range = {}) {
  if (!canvas) return;
  if (!bindlog || bindlog.nEvents === 0) return drawEmpty(canvas, "no binding data — record a trajectory with BindLog attached");
  const ctx = canvas.getContext("2d");
  if (!ctx) return;
  const { W, H } = dimsOf(canvas);
  ctx.fillStyle = CLS.bg; ctx.fillRect(0, 0, W, H);

  // gather contact events
  const contacts = []; // {atom, res, segments: [[t0, t1|null]], openSince}
  const byKey = new Map();
  for (let i = 0; i < bindlog.nEvents; i++) {
    const t = bindlog.evType[i];
    if (t !== 1 && t !== 2) continue;
    const time = bindlog.evTime[i], a = bindlog.evA[i], b = bindlog.evB[i], x = bindlog.evX[i];
    const key = a + ":" + b;
    let c = byKey.get(key);
    if (!c) { c = { atom: a, res: b, segments: [], openSince: null, dists: [] }; byKey.set(key, c); contacts.push(c); }
    if (t === 1) { c.openSince = time; }
    else if (c.openSince !== null) { c.segments.push([c.openSince, time]); c.openSince = null; }
    if (t === 1) c.dists.push(x);
  }
  // close open segments at tTo
  const tTo = range.tTo ?? bindlog.evTime[0];
  for (let i = bindlog.nEvents - 1; i >= 0; i--) { if (bindlog.evType[i] === 1 || bindlog.evType[i] === 2) { /* noop */ } }
  let tMax = 0, tMin = Infinity;
  for (const c of contacts) {
    for (const s of c.segments) { tMax = Math.max(tMax, s[1]); tMin = Math.min(tMin, s[0]); }
    if (c.openSince !== null) { c.segments.push([c.openSince, tTo === tMax && tTo === 0 ? c.openSince : (range.tTo ?? c.openSince + 2)]); }
  }
  if (!isFinite(tMin)) { tMin = 0; }
  if (tMax <= tMin) tMax = tMin + 10;
  const t0 = range.tFrom ?? tMin, t1 = range.tTo ?? tMax;
  const xOf = (t) => 46 + ((t - t0) / (t1 - t0)) * (W - 56);

  // sort rows by total lifetime desc, cap at floor((H-30)/rowH)
  const life = (c) => c.segments.reduce((s, [a, b]) => s + (b - a), 0);
  contacts.sort((p, q) => life(q) - life(p));
  const rowH = 12;
  const nRows = Math.min(contacts.length, Math.floor((H - 34) / rowH));

  // axis
  ctx.strokeStyle = CLS.grid; ctx.beginPath();
  ctx.moveTo(46, H - 18); ctx.lineTo(W - 10, H - 18); ctx.stroke();
  ctx.fillStyle = CLS.dim; ctx.font = "10px ui-monospace, monospace";
  const nTicks = 5;
  for (let k = 0; k <= nTicks; k++) {
    const tt = t0 + (k / nTicks) * (t1 - t0);
    ctx.fillText(tt.toFixed(0) + "ps", xOf(tt) - 8, H - 6);
    ctx.strokeStyle = CLS.grid; ctx.beginPath(); ctx.moveTo(xOf(tt), 22); ctx.lineTo(xOf(tt), H - 18); ctx.stroke();
  }
  ctx.fillText("res", 6, 16);

  // rows
  ctx.font = "10px ui-monospace, monospace";
  for (let r = 0; r < nRows; r++) {
    const c = contacts[r];
    const y = 22 + r * rowH;
    ctx.fillStyle = CLS.dim;
    ctx.fillText(String(c.res), 8, y + 9);
    for (const [ta, tb] of c.segments) {
      const x1 = Math.max(46, xOf(ta)), x2 = Math.min(W - 10, xOf(tb));
      if (x2 <= x1) continue;
      // class by distance: <3.5 H-bond-ish cyan, <4.6 hydrophobic amber, else unknown slate
      const meanD = c.dists.length ? c.dists.reduce((s, d) => s + d, 0) / c.dists.length : 9;
      ctx.fillStyle = meanD < 3.5 ? CLS.hb : meanD < 4.6 ? CLS.hydro : CLS.unk;
      ctx.fillRect(x1, y, x2 - x1, 8);
    }
  }
  // top-contacts side list (drawn as text rows below axis if space, else overlay corner)
  let ly = 22;
  ctx.fillStyle = CLS.text;
  const top = contacts.slice(0, 5);
  // right-side column at fixed x (W-140) with 1px rule
  const lx = Math.max(46, W - 150);
  ctx.strokeStyle = CLS.grid; ctx.beginPath(); ctx.moveTo(lx - 8, 16); ctx.lineTo(lx - 8, H - 18); ctx.stroke();
  for (const c of top) {
    ctx.fillStyle = CLS.text;
    ctx.fillText(`r${c.res} ${life(c).toFixed(1)}ps`, lx, ly + 8);
    ly += 12;
    if (ly > H - 30) break;
  }
  ctx.fillStyle = CLS.dim;
  ctx.fillText(`● hb   ● hydro   ● other`, lx, Math.min(H - 26, ly + 10));
}

/**
 * (b) Energy decomposition — stacked stripes over time + enthalpy table.
 * Energy events: type 0, a = term id, x = value.
 */
export function renderEnergyDecomposition(canvas, bindlog, range = {}) {
  if (!canvas) return;
  if (!bindlog || bindlog.nEvents === 0) return drawEmpty(canvas, "no binding data — record a trajectory with BindLog attached");
  const ctx = canvas.getContext("2d");
  if (!ctx) return;
  const { W, H } = dimsOf(canvas);
  ctx.fillStyle = CLS.bg; ctx.fillRect(0, 0, W, H);

  const TERM_COLORS = ["#38BDF8", "#F87171", "#FBBF24", "#34D399", "#E879F9", "#94A3B8", "#64748B"];
  const TERM_NAMES = ["LJ", "Coul", "HB", "desolv", "pi-X", "cat-pi", "halogen"];

  // collect per-term series
  const series = new Map();
  let tMin = Infinity, tMax = -Infinity, vMin = 0, vMax = 0;
  for (let i = 0; i < bindlog.nEvents; i++) {
    if (bindlog.evType[i] !== 0) continue;
    const t = bindlog.evTime[i], term = bindlog.evA[i], v = bindlog.evX[i];
    if (!series.has(term)) series.set(term, []);
    series.get(term).push([t, v]);
    tMin = Math.min(tMin, t); tMax = Math.max(tMax, t);
    vMin = Math.min(vMin, v); vMax = Math.max(vMax, v);
  }
  if (!series.size) return drawEmpty(canvas);
  if (tMax <= tMin) tMax = tMin + 1;
  if (vMax <= vMin) vMax = vMin + 1;
  const t0 = range.tFrom ?? tMin, t1 = range.tTo ?? tMax;
  const xOf = (t) => 8 + ((t - t0) / (t1 - t0)) * (W - 16);
  const yOf = (v) => 14 + (1 - (v - vMin) / (vMax - vMin)) * (H - 58);

  // zero line
  ctx.strokeStyle = CLS.grid; ctx.beginPath();
  ctx.moveTo(8, yOf(0)); ctx.lineTo(W - 8, yOf(0)); ctx.stroke();

  // per-term polylines (positive space = above zero)
  ctx.font = "10px ui-monospace, monospace";
  let legendX = 8;
  for (const [term, pts] of series) {
    ctx.strokeStyle = TERM_COLORS[term % 7];
    ctx.beginPath();
    for (let k = 0; k < pts.length; k++) {
      const [t, v] = pts[k];
      const x = xOf(t), y = yOf(v);
      k === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
    }
    ctx.stroke();
    // legend
    ctx.fillStyle = TERM_COLORS[term % 7];
    const name = TERM_NAMES[term % 7] ?? `t${term}`;
    ctx.fillText(name, legendX, H - 24);
    legendX += 8 + name.length * 6 + 10;
    if (legendX > W - 60) break;
  }

  // enthalpy table (top-right): mean per term + total ΔH
  let tx = W - 140, ty = 16;
  ctx.strokeStyle = CLS.grid; ctx.strokeRect(tx - 6, 8, 136, 12 + Math.min(series.size, 8) * 12);
  let total = 0;
  for (const [term, pts] of series) {
    const mean = pts.reduce((s, p) => s + p[1], 0) / pts.length;
    total += mean;
    ctx.fillStyle = TERM_COLORS[term % 7];
    ctx.fillText(`${TERM_NAMES[term % 7] ?? "t" + term} ${mean >= 0 ? "+" : ""}${mean.toFixed(2)}`, tx, ty);
    ty += 12;
    if (ty > H - 60) break;
  }
  ctx.fillStyle = CLS.text;
  ctx.fillText(`ΔH = ${total.toFixed(2)} kcal/mol`, tx, ty + 4);

  ctx.fillStyle = CLS.dim;
  ctx.fillText("time →", 8, H - 8);
  ctx.fillText(vMin.toFixed(0), 8, H - 32);
  ctx.fillText(vMax.toFixed(0), 8, 16);
}

/**
 * (c) PMF formation from hill events (well-tempered Gaussians).
 * hills: type 3 events, x = CV, y = height. Renders the PMF at the cursor time.
 */
export function renderPmfFormation(canvas, bindlog, opts = {}) {
  if (!canvas) return;
  // D1 (Loop-1 review): null bindlog → actionable empty state, never throw.
  if (!bindlog) return drawEmpty(canvas, "no binding data — record a trajectory with BindLog attached");
  const ctx = canvas.getContext && canvas.getContext("2d");
  if (!ctx) return;
  const { W, H } = dimsOf(canvas);
  ctx.fillStyle = CLS.bg; ctx.fillRect(0, 0, W, H);
  // collect hills up to cursorT
  const cursorT = opts.cursorT ?? Infinity;
  const hills = [];
  for (let i = 0; i < bindlog.nEvents; i++) {
    if (bindlog.evType[i] !== 3) continue;
    if (bindlog.evTime[i] > cursorT) break;
    hills.push([bindlog.evX[i], bindlog.evY[i]]);
  }
  if (hills.length === 0) return drawEmpty(canvas, "no metadynamics hills deposited yet");

  const binCount = opts.binCount ?? 60;
  const cvRange = opts.cvRange ?? null;
  let cvMin = Infinity, cvMax = -Infinity;
  for (const [cv] of hills) { cvMin = Math.min(cvMin, cv); cvMax = Math.max(cvMax, cv); }
  if (cvRange) { cvMin = cvRange[0]; cvMax = cvRange[1]; }
  if (cvMax <= cvMin) cvMax = cvMin + 1;

  const width = opts.hillWidth ?? 0.5; // Gaussian width in CV units
  const bins = new Float64Array(binCount);
  for (let b = 0; b < binCount; b++) {
    const s = cvMin + ((b + 0.5) / binCount) * (cvMax - cvMin);
    let u = 0;
    for (const [cv, h] of hills) {
      const d = (s - cv) / width;
      u += h * Math.exp(-0.5 * d * d);
    }
    bins[b] = u;
  }
  // PMF = -bias (deposition bias lowers free energy along visited CV)
  let uMax = 0;
  for (const u of bins) uMax = Math.max(uMax, u);
  const xOf = (b) => 8 + (b / (binCount - 1)) * (W - 16);
  const yOf = (u) => 14 + (1 - u / (uMax || 1)) * (H - 40);

  // filled curve: bias accumulated
  ctx.beginPath();
  ctx.moveTo(xOf(0), yOf(bins[0]));
  for (let b = 1; b < binCount; b++) ctx.lineTo(xOf(b), yOf(bins[b]));
  ctx.lineTo(xOf(binCount - 1), yOf(0));
  ctx.lineTo(xOf(0), yOf(0));
  ctx.closePath();
  ctx.fillStyle = "rgba(232,121,249,0.18)"; // magenta @ 18%
  ctx.fill();
  ctx.strokeStyle = CLS.arom;
  ctx.beginPath();
  ctx.moveTo(xOf(0), yOf(bins[0]));
  for (let b = 1; b < binCount; b++) ctx.lineTo(xOf(b), yOf(bins[b]));
  ctx.stroke();

  ctx.fillStyle = CLS.text; ctx.font = "10px ui-monospace, monospace";
  ctx.fillText(`hills ${hills.length}`, 10, 16);
  ctx.fillText(`CV ${cvMin.toFixed(1)}–${cvMax.toFixed(1)} Å`, 10, H - 8);
  // ΔG annotation: deepest bias point
  let best = 0;
  for (let b = 1; b < binCount; b++) if (bins[b] > bins[best]) best = b;
  ctx.fillText(`max bias −${bins[best].toFixed(2)} kcal/mol @CV ${ (cvMin + ((best + 0.5) / binCount) * (cvMax - cvMin)).toFixed(1)}`, W - 240, 16);
}

/**
 * (d) Pareto frontier scatter (from R5 CSV data).
 * tiers: [{name, msPerStep, top1 (0-1), pareto:bool}]
 */
export function renderParetoFrontier(canvas, tiers) {
  if (!canvas) return;
  const ctx = canvas.getContext("2d");
  const { W, H } = dimsOf(canvas);
  ctx.fillStyle = CLS.bg; ctx.fillRect(0, 0, W, H);
  if (!tiers || tiers.length === 0) return drawEmpty(canvas, "no tier benchmark — run scripts/pareto_bench.mjs");

  // axes: x = ms/step (log), y = pose top-1
  const xs = tiers.map(t => Math.log10(Math.max(t.msPerStep, 0.01)));
  const xMin = Math.min(...xs) - 0.3, xMax = Math.max(...xs) + 0.3;
  const yMin = 0.4, yMax = 1.05;
  const xOf = (ms) => 40 + ((Math.log10(Math.max(ms, 0.01)) - xMin) / (xMax - xMin)) * (W - 60);
  const yOf = (tp) => 14 + (1 - (tp - yMin) / (yMax - yMin)) * (H - 40);

  ctx.strokeStyle = CLS.grid;
  ctx.beginPath(); ctx.moveTo(40, 14); ctx.lineTo(40, H - 26); ctx.lineTo(W - 10, H - 26); ctx.stroke();
  ctx.fillStyle = CLS.dim; ctx.font = "10px ui-monospace, monospace";
  for (const ms of [0.01, 0.1, 1, 10, 100]) {
    if (ms < Math.pow(10, xMin) || ms > Math.pow(10, xMax)) continue;
    ctx.fillText(ms + "ms", xOf(ms) - 8, H - 12);
  }
  ctx.fillText("1.0", 14, yOf(1.0) + 4);
  ctx.fillText("0.6", 14, yOf(0.6) + 4);

  // frontier line through pareto points
  const frontier = tiers.filter(t => t.pareto).sort((a, b) => a.msPerStep - b.msPerStep);
  if (frontier.length > 1) {
    ctx.strokeStyle = CLS.arom;
    ctx.beginPath();
    frontier.forEach((t, i) => {
      const x = xOf(t.msPerStep), y = yOf(t.top1);
      i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y);
    });
    ctx.stroke();
  }
  // points + labels
  for (const t of tiers) {
    const x = xOf(t.msPerStep), y = yOf(t.top1 ?? 0.9);
    ctx.fillStyle = t.pareto ? CLS.arom : CLS.unk;
    ctx.fillRect(x - 3, y - 3, 6, 6);
    ctx.fillStyle = CLS.text;
    ctx.fillText(t.name, x + 6, y + 3);
  }
  ctx.fillStyle = CLS.dim;
  ctx.fillText("ms/step →", W - 70, H - 12);
}
