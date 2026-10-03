/**
 * controllers/live-terms.js — Revolution 1 / Issue 4: live per-term mirror.
 *
 * PUBLIC API (all pure w.r.t. the DOM; the only shared state is `state`)
 *   liveTermsWanted(ff)      -> boolean. per-term data wanted whenever a
 *                               ligand is present (no checkbox involved).
 *   applyLiveTrackTerms(ff)  -> boolean. The trackTerms POLICY: BindLog
 *                               capture OR live mirror.
 *   readLiveTerms(ff)        -> snapshot | null. Last per-term vector; missing
 *                               fields read as 0, bindingU falls back to the sum.
 *   formatLiveTermsHUD(ff)   -> string. Compact HUD fragment: scalar U_bind +
 *                               live per-term split ("" when not applicable).
 *   updateLiveTermsMirror(ff, t) -> snapshot | null. Refresh the HUD-cadence
 *                               mirror (reads last ff fields, no extra compute).
 *   liveTermsBindLogView()   -> BindLog-shaped view | null. Synthetic energy
 *                               event stream over the mirror history so the
 *                               tested renderEnergyDecomposition renderer draws
 *                               live per-term stripes without BindLog capture.
 *
 * NEEDS (imports): state from ../ui.js; bindLogWanted from
 * ./physics-tier.js. Exported by src/main.js as a thin re-export so the
 * historical `import ... from "../src/main.js"` call sites keep working.
 *
 * WHY THIS EXISTS (verbatim rationale, unchanged)
 *   lj/coul/hb/desolv/pi/cpi/xb accumulators were filled only when
 *   ff.trackTerms === true (BindLog capture on); a normal Run showed the scalar
 *   U_bind but zeros in the energy channel. Fix: keep the lightweight
 *   accumulators always on while a ligand is present (a few float adds per
 *   cross pair — U/forces bit-identical, see scripts/test_bindlog_integration.mjs)
 *   and mirror the last per-term vector at HUD cadence (10 Hz) + render it into
 *   bindvizEnergy at ≤1 Hz when the BindLog carries no energy events. BindLog
 *   event capture itself stays gated on bindLogWanted() (no extra memory), so
 *   defaults stay backward compatible.
 *
 * ZERO PHYSICS HERE: this module only READS accumulator fields the kernels
 * already wrote and decides the trackTerms flag.
 */

import { state } from "../ui.js";
import { bindLogWanted } from "./physics-tier.js";

/** Ring-buffer length of the mirror history (10 Hz x 12 s). */
export const LIVE_TERMS_N = 120;

/** Live mirror wants per-term data whenever a ligand is present (no checkbox). */
export function liveTermsWanted(ff) { return !!(ff && ff.nLigAtoms > 0); }

/**
 * Apply the live trackTerms policy: BindLog capture OR live mirror.
 * @param {object} ff live force field
 * @returns {boolean} effective trackTerms flag
 */
export function applyLiveTrackTerms(ff) {
  if (!ff) return false;
  ff.trackTerms = bindLogWanted() || liveTermsWanted(ff);
  return ff.trackTerms === true;
}

/** Read the last per-term vector (always safe; missing fields → 0). */
export function readLiveTerms(ff) {
  if (!ff) return null;
  const lj = Number(ff.bindLJU) || 0, coul = Number(ff.bindCoulU) || 0,
    hb = Number(ff.bindHBU) || 0, desolv = Number(ff.desolvU) || 0,
    pi = Number(ff.piU) || 0, cpi = Number(ff.cpiU) || 0, xb = Number(ff.xbU) || 0;
  const b = Number(ff.bindingU);
  return { lj, coul, hb, desolv, pi, cpi, xb, bindingU: Number.isFinite(b) ? b : (lj + coul + hb + desolv + pi + cpi + xb) };
}

/** Compact HUD fragment: scalar U_bind + live per-term split ("" when N/A). */
export function formatLiveTermsHUD(ff) {
  if (!ff || !(ff.nLigAtoms > 0) || !Number.isFinite(ff.bindingU)) return "";
  const t = readLiveTerms(ff);
  const f2 = (v) => (Number.isFinite(v) ? v.toFixed(2) : "—");
  let s = `U_bind = ${f2(t.bindingU)} (LJ ${f2(t.lj)} · Coul ${f2(t.coul)} · HB ${f2(t.hb)} · desolv ${f2(t.desolv)}`;
  if (Math.abs(t.pi) > 1e-9 || Math.abs(t.cpi) > 1e-9 || Math.abs(t.xb) > 1e-9)
    s += ` · π ${f2(t.pi)} · cπ ${f2(t.cpi)} · XB ${f2(t.xb)}`;
  return s + ") kcal/mol  ·  ";
}

/**
 * Refresh the HUD-cadence mirror (no extra ff.compute — reads last fields).
 * @param {object} ff live force field
 * @param {number} t sim time in ps
 */
export function updateLiveTermsMirror(ff, t) {
  if (!liveTermsWanted(ff)) return null;
  const snap = readLiveTerms(ff);
  if (!snap) return null;
  state._liveTerms = snap;
  const hist = (state._liveTermsHist ??= []);
  const lastT = hist.length ? hist[hist.length - 1].t : -Infinity;
  if (hist.length === 0 || t - lastT >= 0.09) {
    hist.push({ t, ...snap });
    if (hist.length > LIVE_TERMS_N) hist.shift();
  }
  return snap;
}

/**
 * Synthetic BindLog view over the mirror history so the tested
 * renderEnergyDecomposition renderer draws live per-term stripes without
 * requiring the BindLog checkbox (same colors/legend, ≤1 Hz).
 */
export function liveTermsBindLogView() {
  const hist = state._liveTermsHist;
  if (!hist || hist.length < 2) return null;
  const n = hist.length, nE = n * 7;
  const evTime = new Float64Array(nE), evType = new Uint8Array(nE),
    evA = new Int32Array(nE), evX = new Float32Array(nE);
  let k = 0;
  for (let i = 0; i < n; i++) {
    const h = hist[i], vals = [h.lj, h.coul, h.hb, h.desolv, h.pi, h.cpi, h.xb];
    for (let c = 0; c < 7; c++, k++) { evTime[k] = h.t; evType[k] = 0; evA[k] = c; evX[k] = vals[c]; }
  }
  return { nEvents: nE, evTime, evType, evA, evX };
}