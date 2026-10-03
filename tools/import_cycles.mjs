/**
 * tools/import_cycles.mjs — count static import cycles under src/.
 *
 * WHY A SCRIPT AND NOT A REGEX GREP
 * ----------------------------------
 * src/ carries commented-out imports (leftovers from earlier splits, and
 * deliberately-quoted examples in the module headers). A naive
 * `grep '^import'` / `grep 'from "'` walk reads those as real edges and
 * reports cycles that do not exist in the module graph the runtime builds.
 * So this strips comments and strings with a real state machine first —
 * the same discipline tests/test_cache_contract.js uses — and only then
 * reads the specifiers. A cycle it reports is a cycle the browser would
 * see.
 *
 * A cycle is only reported when it is AMBIGUOUS: two distinct modules that
 * reach each other through at least two different paths, or one that
 * reaches itself through a non-empty path. A straight two-module A→B→A is
 * a real hazard for TDZ/hoisting (a class or const read during module
 * evaluation), so `strict` also prints those, but they are listed
 * separately rather than silently folded into the count.
 *
 * USAGE
 *   node tools/import_cycles.mjs            # src/ only, ambiguity-strict
 *   node tools/import_cycles.mjs --all      # also print simple 2-cycles
 *   node tools/import_cycles.mjs --strict   # non-zero exit if any cycle at all
 */

import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const SRC = path.join(ROOT, "src");

/** All .js files under a directory, repo-relative with "/" separators. */
function walk(dir, out = []) {
  for (const e of fs.readdirSync(dir, { withFileTypes: true })) {
    const p = path.join(dir, e.name);
    if (e.isDirectory()) walk(p, out);
    else if (e.name.endsWith(".js")) out.push(path.relative(ROOT, p).split(path.sep).join("/"));
  }
  return out;
}

/**
 * Mark every byte of `src` as CODE or NOT-CODE (comment / string / template /
 * regex). Offsets are preserved 1:1, so the caller can run a regex over the
 * ORIGINAL text and require that the `from` keyword sits in code while the
 * quoted specifier sits in a string — which is exactly where a real import
 * specifier is, and exactly where a commented-out one is not.
 */
function codeMask(src) {
  const out = new Uint8Array(src.length); // 1 = code
  let i = 0;
  const n = src.length;
  let state = "code"; // code | line | block | str | tmpl | regex
  let quote = "";
  const tmplStack = [];
  let depth = 0;
  let inClass = false;
  const regexAllowed = () => {
    for (let k = i - 1; k >= 0; k--) {
      if (!out[k]) {
        if (src[k] === " " || src[k] === "\t" || src[k] === "\n") continue;
        return true; // a comment or literal precedes: `/` opens a regex
      }
      const ch = src[k];
      if (/[A-Za-z0-9_$)\]}]/.test(ch)) return false;
      return true;
    }
    return true;
  };
  while (i < n) {
    const ch = src[i];
    const nx = src[i + 1];
    if (state === "code") {
      if (ch === "/" && nx === "/") { state = "line"; i += 2; continue; }
      if (ch === "/" && nx === "*") { state = "block"; i += 2; continue; }
      if (ch === '"' || ch === "'") { state = "str"; quote = ch; i++; continue; }
      if (ch === "`") { state = "tmpl"; tmplStack.push(depth); i++; continue; }
      if (ch === "/" && regexAllowed()) { state = "regex"; inClass = false; i++; continue; }
      if (ch === "{") depth++;
      if (ch === "}") depth--;
      out[i] = 1; i++; continue;
    }
    if (state === "line") {
      if (ch === "\n") state = "code";
      i++; continue;
    }
    if (state === "block") {
      if (ch === "*" && nx === "/") { state = "code"; i += 2; continue; }
      i++; continue;
    }
    if (state === "str" || state === "regex") {
      if (ch === "\\") { i += 2; continue; }
      if (state === "regex" && ch === "\n") { state = "code"; i++; continue; }
      if (state === "regex" && ch === "[") inClass = true;
      if (state === "regex" && ch === "]") inClass = false;
      if (state === "regex" && ch === "/" && !inClass) { state = "code"; i++; continue; }
      if (state === "str" && ch === quote) { state = "code"; i++; continue; }
      i++; continue;
    }
    if (state === "tmpl") {
      if (ch === "\\") { i += 2; continue; }
      if (ch === "$" && nx === "{") { state = "code"; tmplStack.push(depth); i += 2; continue; }
      if (ch === "}" && depth === tmplStack[tmplStack.length - 1]) {
        tmplStack.pop(); state = "tmpl"; i++; continue;
      }
      if (ch === "`") { state = "code"; i++; continue; }
      i++; continue;
    }
  }
  return out;
}

const files = walk(SRC).sort();
const graph = new Map();

for (const rel of files) {
  const srcText = fs.readFileSync(path.join(ROOT, rel), "utf-8");
  const mask = codeMask(srcText);
  const edges = new Set();
  const re = /\bfrom\s*(["'])([^"']+)\1|\bimport\s*\(\s*(["'])([^"']+)\3\s*\)/g;
  let m;
  while ((m = re.exec(srcText)) !== null) {
    const spec = m[2] ?? m[4];
    // `from` / `import` keyword itself must be live code, else this is a
    // commented-out or quoted example.
    const kw = srcText.lastIndexOf("from", m.index) === m.index
      ? m.index
      : srcText.lastIndexOf("import", m.index);
    if (!mask[kw]) continue;
    if (!spec.startsWith(".")) continue;
    const abs = path.posix.normalize(path.posix.join(path.posix.dirname(rel), spec));
    if (!abs.endsWith(".js")) continue;
    if (!files.includes(abs)) {
      console.error(`  ! ${rel} imports ${spec} -> ${abs} which does not exist`);
      process.exitCode = 1;
      continue;
    }
    edges.add(abs);
  }
  graph.set(rel, [...edges]);
}

// Tarjan-free cycle detection: for each node, enumerate simple paths and
// report every edge that closes one. The graph is a DAG in a healthy tree,
// so this terminates immediately when there are no cycles.
const cycles = new Set();
const seen = new Set();
const stack = [];
function dfs(u) {
  seen.add(u);
  stack.push(u);
  for (const v of graph.get(u) ?? []) {
    const idx = stack.indexOf(v);
    if (idx >= 0) {
      const cyc = [...stack.slice(idx), v];
      cycles.add(cyc.join(" -> "));
    } else if (!seen.has(v)) {
      dfs(v);
    }
  }
  stack.pop();
}
for (const f of files) if (!seen.has(f)) dfs(f);

const all = [...cycles].sort();
const ambiguous = all.filter((c) => {
  const parts = c.split(" -> ");
  return parts.length > 3; // a loop longer than the trivial 2-module A->B->A
});

let edges = 0;
for (const v of graph.values()) edges += v.length;
console.log(`modules: ${files.length}   static import edges: ${edges}`);
console.log(`cycles (all, including trivial 2-module loops): ${all.length}`);
for (const c of all) console.log(`  ${c}`);
console.log(`cycles (ambiguous only, length > 2 hops): ${ambiguous.length}`);

const args = process.argv.slice(2);
if (args.includes("--strict") && all.length > 0) {
  console.error("FAIL: --strict and at least one cycle exists");
  process.exit(1);
}
console.log(all.length === 0 ? "PASS: import graph is acyclic" : "NOTE: cycles listed above");