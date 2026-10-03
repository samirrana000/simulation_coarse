/**
 * scope.js — the machine-readable honest-scope table, and the ROADMAP parser
 * that proves the table has not drifted away from the document it claims to
 * quote.
 *
 * THE PROBLEM THIS SOLVES
 * -----------------------
 * The results record has to say, per value, what the number does NOT mean. The
 * canonical wording for that already exists — ROADMAP.md §1 ("No QM/MM, no
 * explicit membrane, no PME in browser v1", plus five hard-no bullets), with
 * docs/LIMITATIONS.md giving the per-value trust boundary and docs/VALIDATION.md
 * giving the measured one ("ranking-only, NOT FEP, NOT absolute Kd").
 *
 * The obvious way to build the record is to type that wording into the record
 * builder. That is exactly the hand-copied string this repo's own history
 * warns about: ROADMAP.md:96 exists as a MEASURABLE grep because a scope claim
 * that only lives in prose is a claim nothing can falsify. A second copy of the
 * scope inside src/ would be a second thing to keep true, and nothing would
 * notice when ROADMAP.md was edited.
 *
 * WHAT IS ACTUALLY DERIVED, AND WHAT IS NOT — STATED PLAINLY
 * ---------------------------------------------------------
 *   DERIVED:  the record's `scope.outOfScope` list is BUILT from the table
 *             below at call time (`scopeBlock()`), so the record can never
 *             disagree with the table. The table's ROADMAP §1 bullets, the §1
 *             guard sentence and the docs/LIMITATIONS.md anchors are checked
 *             against the REAL files by `scopeDrift()`, called by
 *             tests/test_results_record.js on every FAST-tier run — and that
 *             test PROVES it can fail by feeding it a mutated ROADMAP.
 *   NOT DERIVED: the browser cannot read ROADMAP.md at export time (no fetch of
 *             repo-root files, no build step to inline them, and `file://`
 *             origins cannot read it at all). So the table is shipped as code
 *             and held to the documents by a test rather than by a runtime read.
 *             This is the honest limit of "no build step, zero dependencies"; a
 *             build step could inline the parse result and remove the table
 *             entirely. What it costs is one failing test if someone edits
 *             ROADMAP.md without this file — which is the cheap direction to
 *             fail in.
 *
 * THE PER-VALUE HALF
 * ------------------
 * `VALUE_TRUST_BOUNDARY` keys every exported result id to what that quantity is
 * for and what it is not. A result with no entry is a REFUSED EXPORT, not a
 * silently unlabelled number: `resultsExportGate()` in src/results-record.js
 * fails closed when an id is missing, so "we forgot to say what ΔG means" can
 * never ship as a clean-looking JSON file.
 */

/**
 * The one-sentence browser-v1 scope guard. This exact string is asserted to be
 * the text after `**Browser v1 will NOT do:**` on ROADMAP.md:3, which is the
 * line the project's own grep check (ROADMAP.md:42) hits.
 */
export const BROWSER_V1_GUARD = "No QM/MM, no explicit membrane, no PME in browser v1";

// src/dg-estimators.js owns the estimator enumeration and the reconciliation
// bridge; src/dg-chain.js owns the three nature tags and the chain builders.
// Imported (not duplicated) so the record's `scope.deltaGChain` is the live
// table. All three are leaves over src/units.js, so this adds no cycle —
// tests/test_module_size.js's acyclicity check is the thing that would say
// otherwise if that ever stopped being true.
import { DG_ESTIMATORS, DG_BRIDGE } from "./dg-estimators.js";
import { TERM_NATURE } from "./dg-chain.js";

/**
 * Documents the scope statement is sourced from. `file` is repo-relative;
 * `section` is the heading inside it. Kept as data so the record cites
 * documents rather than prose.
 */
export const SCOPE_SOURCES = [
  { file: "ROADMAP.md", section: "1", role: "out-of-scope contract — the hard no's for browser v1" },
  { file: "docs/LIMITATIONS.md", section: null, role: "per-value trust boundary (what breaks, and how far)" },
  { file: "docs/VALIDATION.md", section: null, role: "measured trust boundary (ranking-only, NOT FEP, NOT absolute Kd)" },
];

/**
 * One entry per hard "no" in ROADMAP.md §1.
 *
 * `roadmapKey` is the bullet headline with the leading "No " stripped, and is
 * normalised (`normKey`) before comparison — that is what `scopeDrift()` matches
 * on, so adding a bullet to ROADMAP.md without adding a row here fails a test.
 *
 * `limitsTerms` are phrases that must be present in docs/LIMITATIONS.md. They
 * are what makes the limitations link a check rather than a citation: if
 * LIMITATIONS.md is rewritten and drops the no-PME statement, the table's claim
 * that the record's boundary is documented there stops being true and the test
 * says so.
 *
 * @type {{roadmapKey:string, summary:string, useInstead:string, limitsTerms:string[]}[]}
 */
export const OUT_OF_SCOPE = [
  {
    roadmapKey: "QM/MM",
    summary:
      "No quantum Hamiltonian, no SCF, no bond breaking/formation, no catalysis or metal redox. " +
      "The heavy force field is LJ + screened Coulomb + GB/SA + harmonic bonds/angles/dihedrals, not a semi-empirical or DFT engine.",
    useInstead: "CP2K, ORCA, NAMD/QM-MM, OpenMM + psi4",
    limitsTerms: ["no quantum Hamiltonian"],
  },
  {
    roadmapKey: "explicit membrane",
    summary:
      "No lipid bilayer, no lateral pressure coupling, no anisotropic barostat, no CHARMM36 lipid parameters. " +
      "The model is protein-only Cα ENM or heavy-atom GB/SA in implicit solvent; a GPCR, channel or transporter " +
      "simulated with it produces qualitatively wrong energetics.",
    useInstead: "GROMACS/CHARMM-GUI + explicit lipids, NAMD, OpenMM membrane builder",
    limitsTerms: ["no lipid bilayer"],
  },
  {
    roadmapKey: "PME in browser v1",
    summary:
      "No Particle-Mesh Ewald. Long-range electrostatics beyond 8.5 Å are truncated by the 6.5→8.5 Å switching " +
      "function with no lattice sum; the 8.5→12 Å convergence test leaves ≈5% residual error.",
    useInstead: "GROMACS/AMBER/OpenMM with PME",
    limitsTerms: ["no PME"],
  },
  {
    roadmapKey: "rigorous FEP/TI/MBAR",
    summary:
      "No alchemical intermediates, no soft cores, no replica exchange. The PMF is a 1-D funnel along " +
      "r = |COM_lig − COM_pocket| and the kinetics are a 4-state toy — pedagogical, not converged.",
    useInstead: "rigorous FEP/TI/MBAR or umbrella sampling in GROMACS/AMBER/OpenMM + PLUMED",
    limitsTerms: ["not a converged MSM", "illustrative"],
  },
  {
    roadmapKey: "nucleic acids / glycosylation / PTM libraries beyond simple LJ",
    summary:
      "No nucleic-acid sugar-pucker/dihedral terms, no glycosylation library, no post-translational-modification " +
      "chemistry beyond simple LJ interactions.",
    useInstead: "AMBER OL3/bcs1 for nucleic acids",
    limitsTerms: ["no DNA/RNA backbone terms"],
  },
];

/**
 * Per-value trust boundary. Every result the record can carry MUST have an entry
 * here; `resultsExportGate()` refuses a record whose values are not all listed.
 *
 * Fields per entry:
 *   meaning     what the number legitimately is
 *   notA        what it is NOT — the sentence a citing author must not write
 *   notComputed machinery this tool does not have that the quantity would need
 *               to mean more than it means
 *   substitute  the tool that does compute the thing properly
 *   limitsTerms phrases that must appear in docs/LIMITATIONS.md
 *
 * @type {Record<string, {meaning:string, notA:string, notComputed:string[], substitute:string, limitsTerms:string[]}>}
 */
export const VALUE_TRUST_BOUNDARY = {
  dg_bind: {
    meaning: "a ranking indicator between systems simulated with identical settings, seed protocol and hill budget",
    notA: "an absolute binding free energy, a ΔG° prediction, or a K_D",
    notComputed: [
      "alchemical intermediates or soft cores",
      "replica exchange",
      "converged 1-D sampling along a single radial CV",
      "long-range electrostatics beyond the 8.5 Å truncation (no PME)",
    ],
    substitute: "alchemical FEP/TI or umbrella sampling in GROMACS/AMBER/OpenMM + PLUMED",
    limitsTerms: ["1D funnel CV only"],
  },
  occupancy_bound_fraction: {
    meaning: "the fraction of recorded frames whose ligand COM lies inside the funnel pocket radius",
    notA: "a bound-state population, an equilibrium constant, or evidence of convergence",
    notComputed: [
      "blocked sampling or replica exchange across the bound/unbound transition",
      "an orthogonal CV that separates encounter and bound states",
    ],
    substitute: "a converged MSM (PyEMMA/MSMBuilder) on explicit-solvent trajectories",
    limitsTerms: ["4-state kinetics toy"],
  },
  contact_lifetime_mean_ps: {
    meaning: "the mean longest uninterrupted streak of a protein–ligand contact below 6 Å, in simulated ps",
    notA: "a residence time or an experimental off-rate",
    notComputed: [
      "an explicit-solvent sampling rate — the streak is limited by the recording stride",
      "rare-event resampling",
    ],
    substitute: "kinetic network models on explicit-solvent MD (PyEMMA/MSMBuilder)",
    limitsTerms: ["illustrative sampling"],
  },
  rmsip_pca_enm: {
    meaning: "the overlap between the top-10 simulated Cα fluctuation modes and the 10 softest ENM modes",
    notA: "a statement that the elastic network is correct, or a validation score",
    notComputed: [
      "any comparison against a reference all-atom trajectory",
      "cross-validation of the mode selection (k is chosen by the caller)",
    ],
    substitute: "comparison against an explicit-solvent reference trajectory in MDAnalysis",
    limitsTerms: ["Coarse-grained / heavy-atom approximate force field"],
  },
  b_factor_pearson_r: {
    meaning: "the Pearson correlation between simulated and crystallographic B-factors over residues that have both",
    notA: "a model-quality metric, or evidence the force field reproduces dynamics",
    notComputed: [
      "an error bar — the code computes no confidence interval on r",
      "a time-correlation function or its Fourier transform, which is what B-factors encode physically",
    ],
    substitute: "direct comparison against a longer, all-atom explicit-solvent trajectory",
    limitsTerms: ["heuristic parameters"],
  },
  b_factor_mean_sim_ang2: {
    meaning: "the mean of (8π²/3)⟨Δr²⟩ over the protein beads over the recorded window",
    notA: "an experimental temperature factor — the prefactor is the right one but the window is ~100 ps, not the ms an experiment averages",
    notComputed: [
      "a time-correlation function, so no harmonic correction and no deconvoluted slow mode",
      "an uncertainty on the mean; the code computes none",
    ],
    substitute: "long explicit-solvent trajectories analysed with MDAnalysis RMSF",
    limitsTerms: ["coarse time scales"],
  },
  rmsf_mean_ang: {
    meaning: "the root-mean-square fluctuation of the superposed protein beads about the native structure",
    notA: "a converged fluctuation amplitude",
    notComputed: [
      "an uncertainty on the mean",
      "convergence: the value depends directly on how many frames were recorded",
    ],
    substitute: "long explicit-solvent trajectories analysed with MDAnalysis RMSF",
    limitsTerms: ["coarse time scales"],
  },
  var_top_k: {
    meaning: "the fraction of total Cα variance captured by the top k=10 principal components of the recorded window",
    notA: "the variance fraction the protein actually explores",
    notComputed: [
      "convergence — an under-sampled window inflates varMode1 by construction",
      "any error bar",
    ],
    substitute: "a converged covariance from explicit-solvent trajectories",
    limitsTerms: ["coarse time scales"],
  },
  pocket_volume_mean_ang3: {
    meaning: "the mean geometric pocket volume over the recorded frames for the selected residue set",
    notA: "a pocket volume from crystallography, or a free-energy contribution",
    notComputed: [
      "an error bar on the mean — the detector returns a sample standard deviation of the per-frame track (a dispersion), which is not one",
      "any ligand-independent or biased sampling of pocket opening",
    ],
    substitute: "fpocket/pocket analysis on an ensemble from explicit-solvent MD",
    limitsTerms: ["Coarse-grained / heavy-atom approximate force field"],
  },
  pocket_volume_sd_ang3: {
    meaning: "the sample standard deviation of the per-frame pocket-volume track",
    notA: "an uncertainty on the mean volume — a track standard deviation is a dispersion, and reporting it as an error bar would be the classic false precision",
    notComputed: [
      "any standard error, confidence interval or block analysis of the track",
    ],
    substitute: "block-averaging of an explicit-solvent ensemble",
    limitsTerms: ["Coarse-grained / heavy-atom approximate force field"],
  },
  dH_total: {
    meaning: "the per-frame mean binding enthalpy decomposed over the tracked force-field terms",
    notA: "a measured binding enthalpy, or one that includes an alchemical cycle",
    notComputed: [
      "an alchemical intermediate path, so no cycle closure and no cancellation of the intramolecular terms",
      "any term the selected ligand subset does not cover",
    ],
    substitute: "FEP/TI in GROMACS/AMBER/OpenMM",
    limitsTerms: ["Implicit solvent only"],
  },
  dS_total: {
    meaning: "a Schlitter quasi-harmonic + torsion + SASA-burial entropy difference between a recorded holo leg and an internal apo relaxation",
    notA: "a configurational or solvent entropy of the binding process",
    notComputed: [
      "a converged sampling of the apo leg (it is an internal relaxation of the same coordinates, not an independent simulation)",
      "a standard error on ΔS_pocket, which the code computes none of — replica SD dominates",
      "explicit-water solvent entropy",
    ],
    substitute: "quasi-harmonic analysis of long, independent holo and apo explicit-solvent runs",
    limitsTerms: ["no explicit water entropy"],
  },
  dg_estimate: {
    meaning: "ΔH − T·ΔS from the same decomposition; a self-consistent internal quantity",
    notA: "a binding free energy in the thermodynamic sense — the apo leg is a relaxation of the holo structure, so the entropy difference is not a basin-to-basin free-energy difference",
    notComputed: [
      "anything the two terms above do not contain",
      "a standard error; VALIDATION measures replica SD ≈ 7 kcal/mol on the −TΔS leg, which dominates",
    ],
    substitute: "alchemical FEP/TI or experimental ITC/NMR",
    limitsTerms: ["Implicit solvent only"],
  },
  dG_jarzynski: {
    meaning: "a Jarzynski estimate from a 4-pull constant-velocity steered ensemble, with a seeded nonparametric bootstrap SE",
    notA: "an equilibrium free energy — the ensemble is strongly dissipative, and Jarzynski is exponentially ill-conditioned in dissipation",
    notComputed: [
      "enough pulls for the bootstrap to mean much (n = 4 here)",
      "reverse-direction or bidirectional reweighting",
    ],
    substitute: "umbrella sampling or PMF in GROMACS + PLUMED",
    limitsTerms: ["1D funnel CV only"],
  },
  koff_surrogate_kT: {
    meaning: "a unitless ranking score in kT built from the rupture work of the same 4 pulls",
    notA: "a k_off rate, and not convertible to one — it carries no attempt frequency and no state count",
    notComputed: [
      "an MSM, so no state assignment and no rate prefactor",
      "an uncertainty",
    ],
    substitute: "MSMBuilder/PyEMMA rate estimation on explicit-solvent trajectories",
    limitsTerms: ["4-state kinetics toy"],
  },
};

/** Normalise a ROADMAP bullet headline for comparison (case/space/punctuation-free-ish). */
export function normKey(s) {
  return String(s == null ? "" : s)
    .replace(/\s+/g, " ")
    .trim()
    .replace(/[.]+$/, "")
    .toLowerCase();
}

/**
 * Extract the §1 guard sentence from ROADMAP markdown.
 *
 * The source line is ROADMAP.md:3:
 *   `> **Browser v1 will NOT do:** No QM/MM, no explicit membrane, no PME in browser v1.`
 * so the sentence is everything after that bold label, minus the trailing full
 * stop. Returns null when the line is absent or reshaped.
 * @param {string} md ROADMAP.md text
 * @returns {string|null}
 */
export function parseRoadmapGuard(md) {
  const m = /^\s*>\s*\*\*Browser v1 will NOT do:\*\*\s*(.+?)\s*$/m.exec(String(md == null ? "" : md));
  return m ? m[1].replace(/[.]+$/, "") : null;
}

/**
 * Parse the hard "no" bullets out of ROADMAP.md §1.
 *
 * §1 is delimited by its own `## 1.` heading and the next `## ` heading, so the
 * §2 "what we WILL do" list cannot leak in and §4's restatement cannot either.
 * Each bullet's bold headline is the key; the rest of the bullet (summary +
 * optional `*Use instead:*` line) is folded into one whitespace-normalised blob.
 * @param {string} md ROADMAP.md text
 * @returns {{key:string, headline:string, body:string, useInstead:string|null}[]}
 */
export function parseRoadmapOutOfScope(md) {
  const text = String(md == null ? "" : md);
  const start = text.search(/^##\s+1\.\s/m);
  if (start < 0) return [];
  const rest = text.slice(start + 1);
  const end = rest.search(/^##\s+\d/m);
  const section = end < 0 ? rest : rest.slice(0, end);
  const out = [];
  const lines = section.split(/\r?\n/);
  for (let i = 0; i < lines.length; i++) {
    const head = /^-\s+\*\*(.+?)\*\*/.exec(lines[i]);
    if (!head) continue;
    // Fold the indented continuation lines of this bullet into its body.
    let body = lines[i].slice(lines[i].indexOf("**", 2) + 2);
    for (let k = i + 1; k < lines.length && /^\s{2,}\S/.test(lines[k]); k++) {
      body += " " + lines[k].trim();
      i = k;
    }
    const use = /\*Use instead:\*\s*(.+?)[.]?\s*$/i.exec(body);
    const headline = head[1].replace(/[.]+$/, "");
    out.push({
      key: normKey(headline.replace(/^No\s+/i, "")),
      headline,
      body: body.replace(/\s+/g, " ").trim(),
      useInstead: use ? use[1].replace(/[.,]\s*$/, "").trim() : null,
    });
  }
  return out;
}

/**
 * Compare the shipped table against the real documents. Returns the list of
 * problems; empty means the table still says what the docs say.
 *
 * Called by tests/test_results_record.js against ROADMAP.md, docs/LIMITATIONS.md
 * and docs/VALIDATION.md read from disk. It is a pure function of its inputs so
 * the same test can hand it a MUTATED document and prove it reports the drift.
 *
 * @param {{roadmap?:string, limitations?:string, validation?:string}} docs
 * @returns {string[]} problems (empty = no drift)
 */
export function scopeDrift({ roadmap = "", limitations = "", validation = "" } = {}) {
  const problems = [];

  const guard = parseRoadmapGuard(roadmap);
  if (guard === null) {
    problems.push("ROADMAP.md: the `> **Browser v1 will NOT do:**` guard line was not found — update BROWSER_V1_GUARD and this parser together");
  } else if (guard !== BROWSER_V1_GUARD) {
    problems.push(`BROWSER_V1_GUARD drifted from ROADMAP.md: table "${BROWSER_V1_GUARD}" vs document "${guard}"`);
  }

  const bullets = parseRoadmapOutOfScope(roadmap);
  if (bullets.length === 0) {
    problems.push("ROADMAP.md §1 yielded no `- **No …**` bullets — the §1 parser is not reading the section it claims to");
  }
  const covered = new Map(OUT_OF_SCOPE.map((e) => [normKey(e.roadmapKey), e]));
  for (const b of bullets) {
    if (!covered.has(b.key)) {
      problems.push(`ROADMAP.md §1 bullet "No ${b.headline.replace(/^No\s+/i, "")}" is not covered by OUT_OF_SCOPE (no entry with roadmapKey "${b.key}")`);
    }
  }
  const present = new Set(bullets.map((b) => b.key));
  for (const e of OUT_OF_SCOPE) {
    if (!present.has(normKey(e.roadmapKey))) {
      problems.push(`OUT_OF_SCOPE entry "${e.roadmapKey}" no longer matches any ROADMAP.md §1 bullet — delete it or fix roadmapKey`);
    }
    for (const term of e.limitsTerms) {
      if (limitations && !limitations.includes(term)) {
        problems.push(`docs/LIMITATIONS.md no longer contains "${term}", which OUT_OF_SCOPE["${e.roadmapKey}"] cites it for`);
      }
    }
  }

  if (limitations) {
    for (const [id, v] of Object.entries(VALUE_TRUST_BOUNDARY)) {
      for (const term of v.limitsTerms) {
        if (!limitations.includes(term)) {
          problems.push(`docs/LIMITATIONS.md no longer contains "${term}", cited by VALUE_TRUST_BOUNDARY.${id}`);
        }
      }
    }
  }
  if (validation) {
    for (const phrase of ["NOT FEP", "ranking-only"]) {
      if (!validation.includes(phrase)) {
        problems.push(`docs/VALIDATION.md no longer contains "${phrase}" — the record's ranking-only claim is no longer documented there`);
      }
    }
  }
  return problems;
}

/**
 * The `scope` block as it appears in an exported record: built from the tables
 * at call time, so a record can never disagree with the table it ships with.
 *
 * `deltaGChain` is the pointer to src/dg-chain.js — the ONE enumeration of every
 * ΔG-like estimator in the app, each tagged measured / assumed / out-of-scope,
 * plus the `bridge` stating why the HUD estimator and the VALIDATION.md estimator
 * cannot be reconciled with what this tool implements. It is embedded by
 * REFERENCE, not by value: scope.js and dg-chain.js are both leaves over
 * src/units.js, and duplicating the table would create a second copy of the
 * honesty claim — which is the exact failure src/results-record.js's header
 * describes for hand-copied strings. tests/test_dg_drift.js asserts the shipped
 * record's pointer resolves to the live table.
 * @returns {{browserV1Guard:string, verdict:string, outOfScope:object[], sources:object[], derivation:object, deltaGChain:object}}
 */
export function scopeBlock() {
  return {
    browserV1Guard: BROWSER_V1_GUARD,
    verdict: "ranking-only. NOT FEP. NOT an absolute K_D. This browser tool is not a substitute for rigorous alchemical free-energy calculation.",
    outOfScope: OUT_OF_SCOPE.map((e) => ({
      key: e.roadmapKey,
      summary: e.summary,
      useInstead: e.useInstead,
      documentedIn: "ROADMAP.md §1",
    })),
    sources: SCOPE_SOURCES.map((s) => ({ ...s })),
    deltaGChain: {
      module: "src/dg-chain.js",
      estimators: DG_ESTIMATORS.map((e) => ({
        id: e.id, label: e.label, codePath: e.codePath,
        nature: e.nature, onDisplayPath: e.onDisplayPath,
        uncertaintyKind: e.uncertainty?.kind ?? null,
        noErrorBar: e.uncertainty ? null : (e.noErrorBar ?? null),
      })),
      chainNatureTagSet: [...TERM_NATURE],
      bridge: DG_BRIDGE.map((b) => ({ id: b.id, headline: b.headline, nature: b.nature, missing: b.missing ?? null })),
      rule: "every ΔG-like number the app prints is in `estimators`; a ΔG that appears in the UI and is not in this table is a bug (tests/test_dg_drift.js)",
      canReconcile: false,
      canReconcileWhy: "the HUD estimator (WTM funnel bias, one radial CV) and the VALIDATION.md estimator (ΔH − T·ΔS, holo vs internal apo) share no term; `bridge` names the missing physics term by term",
    },
    derivation: {
      builtFrom: "src/scope.js OUT_OF_SCOPE + VALUE_TRUST_BOUNDARY, plus src/dg-chain.js DG_ESTIMATORS + DG_BRIDGE",
      heldToDocsBy: "tests/test_results_record.js (scopeDrift) + tests/test_dg_drift.js (doc↔code number drift), FAST tier",
      runtimeDocRead: false,
      runtimeDocReadWhy:
        "no build step and zero dependencies: the browser cannot inline ROADMAP.md at export time and a file:// origin cannot read it at all, so the table ships as code and a failing test is what catches an edit to the document",
    },
  };
}
