#!/usr/bin/env node
/**
 * tests/test_docs_index.js — no-orphan-doc test for docs/README.md.
 *
 * THE DEFECT THIS GUARDS AGAINST
 * ------------------------------
 * docs/ held ~36 markdown files / ~4 900 lines with NO index. A reader
 * arriving from README.md could not tell which document answers "can I trust
 * this ΔG number?" — so the honest-scope work (ROADMAP.md §1's out-of-scope
 * contract, docs/APPLICABILITY.md's trust-boundary table) existed but was
 * unsurfaced. Unsurfaced work has zero perceived value.
 *
 * THE CONTRACT
 * ------------
 *   1. NO ORPHANS   — every file under docs/ is linked from docs/README.md.
 *                     Headline case: every *.md under docs/.
 *   2. NO DUPLICATES— no file under docs/ is linked more than once, so the
 *                     index stays a map and does not become a second copy.
 *   3. NO BROKEN    — every relative link target in docs/README.md resolves
 *                     to a file that exists on disk.
 *   4. TRUST BOUNDARY ROUTED — the repo-root ROADMAP.md and the out-of-tree
 *                     bench/vs_gromacs.md are each linked AT LEAST ONCE from
 *                     the index (they are the two documents most in need of
 *                     surfacing, and neither lives under docs/). Repeats are
 *                     allowed here: the no-duplicate rule is scoped to
 *                     docs/, since an out-of-tree doc may legitimately be
 *                     referenced from several questions.
 *
 * ONE DIRECTIONAL EXEMPTION
 * -------------------------
 * docs/README.md itself is exempt from (1): it is the index, so it is not
 * expected to link to itself. Nothing else is exempt.
 *
 * Only *relative* links that resolve INSIDE docs/ are counted toward (1)/(2);
 * `../ROADMAP.md`, `../bench/...`, absolute paths, and URLs are handled by (4)
 * or ignored. Links are parsed as markdown inline-link targets, so prose that
 * merely names a file (e.g. "see ROADMAP.md") does NOT count as a link — the
 * doc must actually be reachable by clicking it.
 *
 * Runnable standalone: node tests/test_docs_index.js
 * Runs as part of the FAST tier: it is registered in tests/suites.js.
 */

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const HERE = path.dirname(fileURLToPath(import.meta.url));
const ROOT = path.resolve(HERE, "..");
const DOCS = path.join(ROOT, "docs");
const INDEX = path.join(DOCS, "README.md");

/** The index may not be expected to link to itself. */
const SELF = "docs/README.md";

/**
 * Out-of-docs/ links the index MUST carry. Each is linked exactly once.
 * These are the two documents a reader most needs and least finds: the
 * repo-root out-of-scope contract and the external benchmark.
 */
const REQUIRED_EXTERNAL = ["ROADMAP.md", "bench/vs_gromacs.md"];

let passed = 0;
let failed = 0;

function assert(condition, message) {
  if (condition) {
    passed++;
    console.log(`  ✓ ${message}`);
  } else {
    failed++;
    console.error(`  ✗ FAIL: ${message}`);
  }
}

/** Repo-relative, "/" separators. */
function rel(p) {
  return path.relative(ROOT, p).split(path.sep).join("/");
}

/** Recursive relative paths of every file under dir (excluding dir itself). */
function walk(dir, out = []) {
  if (!fs.existsSync(dir)) return out;
  for (const e of fs.readdirSync(dir, { withFileTypes: true })) {
    const p = path.join(dir, e.name);
    if (e.isDirectory()) walk(p, out);
    else out.push(rel(p));
  }
  return out;
}

console.log("=== docs/ index no-orphan check ===");

// ---------------------------------------------------------------------
// 0. Precondition: the index exists and is non-trivial.
if (!fs.existsSync(INDEX)) {
  console.error(`  ✗ FAIL: docs/README.md does not exist — the docs index is the thing under test`);
  console.error(`=== test_docs_index: 0 PASSED, 1 FAILED ===`);
  process.exit(1);
}
const indexBody = fs.readFileSync(INDEX, "utf-8");
assert(indexBody.length > 2000, `docs/README.md is a real index (${indexBody.length} bytes), not a stub`);

// ---------------------------------------------------------------------
// 1. Inventory: every file under docs/ (markdown explicitly, and all files
//    for good measure — the CSV and the blueprint are discoverable too).
const allFiles = walk(DOCS).sort();
const mdFiles = allFiles.filter((f) => f.toLowerCase().endsWith(".md"));
const subjects = allFiles.filter((f) => f !== SELF);

assert(
  allFiles.length > 0,
  `docs/ inventory is non-empty (${allFiles.length} files, ${mdFiles.length} markdown)`
);

// ---------------------------------------------------------------------
// 2. Parse the index's inline-link targets and resolve them.
const targets = [];
for (const m of indexBody.matchAll(/\]\(\s*([^)\s]+)(?:\s+"[^"]*")?\s*\)/g)) {
  const t = m[1];
  if (/^[a-z][a-z0-9+.-]*:/i.test(t)) continue; // http:, mailto:, data: …
  if (t.startsWith("#")) continue; // in-page anchor
  const clean = t.split("#")[0].split("?")[0];
  if (!clean) continue;
  targets.push({ raw: t, abs: path.resolve(DOCS, clean) });
}
assert(targets.length > 0, `docs/README.md contains relative links (${targets.length} parsed)`);

/** Count of in-docs/ files linked, keyed by repo-relative path. */
const inDocs = new Map();
for (const { abs } of targets) {
  const r = rel(abs);
  if (!r.startsWith("docs/") || r === SELF) continue;
  inDocs.set(r, (inDocs.get(r) || 0) + 1);
}

/** Count of out-of-docs/ links, keyed by the path *below* the repo root. */
const outOfDocs = new Map();
for (const { abs } of targets) {
  const r = rel(abs);
  if (r.startsWith("..")) continue;
  if (r.startsWith("docs/")) continue;
  if (r === SELF) continue;
  outOfDocs.set(r, (outOfDocs.get(r) || 0) + 1);
}

// ---------------------------------------------------------------------
// 3. Contract 1 — no orphans. Every file under docs/ is linked.
const orphans = subjects.filter((f) => !inDocs.has(f));
assert(
  orphans.length === 0,
  `no orphaned docs: all ${subjects.length} files under docs/ are linked from docs/README.md` +
    (orphans.length ? ` — ORPHANS: ${orphans.join(", ")}` : "")
);

const mdOrphans = mdFiles.filter((f) => f !== SELF && !inDocs.has(f));
assert(
  mdOrphans.length === 0,
  `no orphaned markdown: all ${mdFiles.length - 1} *.md under docs/ (excluding the index) are linked` +
    (mdOrphans.length ? ` — ORPHANS: ${mdOrphans.join(", ")}` : "")
);

// ---------------------------------------------------------------------
// 4. Contract 2 — no duplicates. Exactly once, so the index stays a map.
const dupes = [...inDocs.entries()].filter(([, n]) => n > 1).map(([f, n]) => `${f} x${n}`);
assert(
  dupes.length === 0,
  `no file under docs/ is linked more than once (${inDocs.size} files linked exactly once)` +
    (dupes.length ? ` — DUPLICATES: ${dupes.join(", ")}` : "")
);

// ---------------------------------------------------------------------
// 5. Contract 3 — no broken links. Every relative target resolves.
const broken = targets.filter(({ raw, abs }) => !fs.existsSync(abs)).map(({ raw }) => raw);
assert(
  broken.length === 0,
  `every relative link in docs/README.md resolves to an existing file (${targets.length} links)` +
    (broken.length ? ` — BROKEN: ${broken.join(", ")}` : "")
);

// ---------------------------------------------------------------------
// 6. Contract 4 — the trust boundary is routed from the index.
for (const want of REQUIRED_EXTERNAL) {
  const n = outOfDocs.get(want) || 0;
  assert(
    n >= 1,
    `docs/README.md routes to ${want} (${n} link${n === 1 ? "" : "s"})`
  );
}

console.log(
  `\n=== test_docs_index (${subjects.length} docs linked, ${targets.length} links, ` +
    `${REQUIRED_EXTERNAL.length} out-of-tree links required): ${passed} PASSED, ${failed} FAILED ===`
);
if (failed) process.exit(1);
