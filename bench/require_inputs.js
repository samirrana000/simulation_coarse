/**
 * bench/require_inputs.js — S5 preflight guard for bench scripts.
 *
 * WHY
 * ---
 * `data/coreset/**` (1712 files, 398.79 MiB) is deliberately untracked in git —
 * see data/README.md and .gitignore. Measured with
 * `node scripts/fetch_coreset.mjs --prove-unread`: **zero** files under
 * `data/coreset/` are read by any tracked code, so no bench needs them today.
 *
 * That is a fact, not a licence to measure on a partial set. The pre-existing
 * F55 bench demonstrated exactly the failure mode this module exists to stop:
 * with `1ubq.pdb` absent, `bench/pdbbind_gb.js` swallowed the error, dropped
 * 1UBQ from the average, printed
 *
 *     RMSE = 0.96 kcal/mol over 2 PDBs
 *     PASS: RMSE <2.5 kcal/mol
 *
 * and exited 0. The RMSE had *improved* (1.53 -> 0.96) purely because a target
 * vanished. Per ROADMAP.md section 1, overselling a broken run is this
 * project's stated worst failure mode.
 *
 * CONTRACT
 * --------
 * `requireDataInputs()` resolves a bench's ENTIRE declared input set before any
 * number is produced. If anything is missing, corrupt, or lives under
 * `data/coreset/`, it prints a specific, actionable report and exits 1. A bench
 * that reaches the end of preflight has a complete, hash-checked input set and
 * its printed numbers mean what they say.
 *
 * Zero dependencies: node builtins only.
 */

import fs from "fs";
import path from "path";
import crypto from "crypto";
import { fileURLToPath } from "url";

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.resolve(__dirname, "..");
const MANIFEST = path.join(ROOT, "data", "manifest.json");

/** Same search order the benches have always used, plus the coreset layout. */
function candidatesFor(name) {
  return [
    path.resolve(process.cwd(), name),
    path.resolve(ROOT, name),
    path.resolve(ROOT, path.dirname(name), path.basename(name)),
    path.resolve(process.cwd(), "data", name),
    // coreset layout: data/coreset/<id>/<id>_<kind>
    path.resolve(ROOT, "data", "coreset", name),
  ];
}

function isCoreset(abs) {
  return abs.split(path.sep).includes("coreset");
}

function sha256File(abs) {
  return "sha256:" + crypto.createHash("sha256").update(fs.readFileSync(abs)).digest("hex");
}

function loadManifestHashes() {
  try {
    const m = JSON.parse(fs.readFileSync(MANIFEST, "utf-8"));
    const flat = {};
    for (const [k, v] of Object.entries(m)) {
      if (typeof v === "string" && v.startsWith("sha256:")) flat[k] = v;
    }
    if (m.coreset && m.coreset.files) {
      for (const [k, v] of Object.entries(m.coreset.files)) flat[`data/${k}`] = v;
    }
    return flat;
  } catch {
    return {};
  }
}

/**
 * Assert that every declared input resolves to an existing, hash-clean file.
 *
 * @param {object} spec
 * @param {string} spec.bench      bench name, used in the failure header
 * @param {Array<{names: string[], why: string}>} spec.inputs
 *        `names` is the fallback chain, most-preferred first (historical order:
 *        lowercase, UPPERCASE, then the `data/` preset path).
 * @param {string} [spec.note]     one line of provenance, printed on success
 * @returns {Map<string,string>}   primary name -> absolute path
 * @throws {never} on any gap — calls process.exit(1) after printing the report
 */
export function requireDataInputs({ bench, inputs, note }) {
  const manifest = loadManifestHashes();
  const resolved = new Map();
  const problems = [];

  for (const input of inputs) {
    let found = null;
    const tried = [];
    for (const name of input.names) {
      for (const cand of candidatesFor(name)) {
        if (tried.includes(cand)) continue;
        tried.push(cand);
        if (fs.existsSync(cand) && fs.statSync(cand).isFile()) {
          found = { name, path: cand };
          break;
        }
      }
      if (found) break;
    }

    if (!found) {
      problems.push({ kind: "missing", input, tried, coreset: input.names.some((n) => n.includes("coreset")) });
      continue;
    }

    // Integrity: if data/manifest.json knows this path, its hash must match.
    const rel = path.relative(ROOT, found.path).split(path.sep).join("/");
    const want = manifest[rel] || manifest[found.name];
    if (want) {
      const got = sha256File(found.path);
      if (got !== want) {
        problems.push({ kind: "sha256", input, found, want, got, rel });
        continue;
      }
    }

    resolved.set(input.names[0], found.path);
  }

  if (problems.length === 0) {
    if (note) console.log(`[${bench}] inputs OK (${resolved.size}/${inputs.length}): ${note}`);
    return resolved;
  }

  // ---------------------------------------------------------------- FAIL LOUD
  const coresetProblems = problems.filter((p) => p.kind === "missing" && p.coreset);
  const presetProblems = problems.filter((p) => p.kind === "missing" && !p.coreset);
  const hashProblems = problems.filter((p) => p.kind === "sha256");

  console.error("");
  console.error(`[${bench}] ABORT — incomplete or corrupt input set (${problems.length} problem(s)).`);
  console.error(`[${bench}] No number is being reported. Refusing to measure on a partial set.`);
  console.error("");

  for (const p of hashProblems) {
    console.error(`[${bench}] CORRUPT: ${p.input.names[0]}`);
    console.error(`[${bench}]   resolved to : ${p.rel}`);
    console.error(`[${bench}]   expected    : ${p.want}`);
    console.error(`[${bench}]   actual      : ${p.got}`);
    console.error(`[${bench}]   This file is tracked by data/manifest.json and has been modified.`);
    console.error(`[${bench}]   Restore it:  git checkout -- ${p.rel}`);
    console.error("");
  }

  for (const p of presetProblems) {
    console.error(`[${bench}] MISSING: ${p.input.names.join("  |  ")}   (needed for: ${p.input.why})`);
    for (const t of p.tried) console.error(`[${bench}]   tried: ${t}`);
    console.error("");
  }

  if (coresetProblems.length) {
    console.error(`[${bench}] MISSING from data/coreset (untracked by design).`);
    for (const p of coresetProblems) {
      console.error(`[${bench}]   ${p.input.names.join("  |  ")}   (needed for: ${p.input.why})`);
    }
    console.error("");
    console.error(`[${bench}] run: node scripts/fetch_coreset.mjs --restore-from <extracted-pdbbind-sybyl-tree>`);
    console.error(`[${bench}] then: node scripts/fetch_coreset.mjs   (verifies SHA256 against data/manifest.json)`);
    console.error("");
  }

  console.error(`[${bench}] This bench reads ONLY repo-root presets; it does not need data/coreset.`);
  console.error(`[${bench}] See data/README.md "Runtime data contract" for the full table.`);
  console.error("");
  process.exit(1);
}

/**
 * Fail loudly if a target failed to compute *after* preflight passed.
 * Guards the "swallow the error, shrink the average, still print PASS" bug.
 *
 * @param {string} bench
 * @param {Array<{pdbId: string, error: string}>} failures
 * @param {number} expected number of targets the bench declared
 */
export function requireAllTargetsComputed(bench, failures, expected) {
  if (failures.length === 0) return;
  console.error("");
  console.error(`[${bench}] ABORT — ${failures.length}/${expected} target(s) failed to compute.`);
  for (const f of failures) console.error(`[${bench}]   ${f.pdbId}: ${f.error}`);
  console.error(
    `[${bench}] An RMSE over the surviving subset is not a result. It is biased by which\n` +
      `[${bench}] targets happened to survive, and dropping a bad target flatters the metric.`
  );
  console.error(`[${bench}] Reporting nothing. Fix the inputs and re-run.`);
  console.error("");
  process.exit(1);
}
