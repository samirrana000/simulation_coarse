#!/usr/bin/env bash
# check.sh — A01 deterministic syntax gate.
#
# Runs `node --check` on EVERY JavaScript source file in the repo and exits
# non-zero if any file fails.
#
# HISTORY OF THIS GATE (both bugs were silent, so both are worth recording)
#
#   Bug 1 — it never descended into subdirectories. The original gate was
#   `node --check src/*.js src/physics/*.js`, which missed src/chem,
#   src/compute, src/analysis, src/capture and src/physics/{forcefield,
#   solvation,integrators}.
#
#   Bug 2 — `node --check` accepts exactly ONE file argument. Any further
#   arguments land in process.argv and are silently ignored, so the whole
#   command only ever checked src/analysis-panel.js (1 of 67).
#
#   Bug 3 (fixed here) — the recursive version only globbed `*.js`. The repo
#   also carries `.mjs` modules: 14 of them under scripts/ and evolve/,
#   including every doc-citation validator, several SLOW-tier validation
#   scripts, and evolve/evolve.mjs itself (the evolution gate). A syntax
#   error in any of those shipped green and only blew up at runtime.
#   .mjs is the extension this project uses for anything run directly by
#   node rather than imported by a browser, so it is load-bearing, not
#   incidental.
#
# SCOPE
#   Roots: every directory that actually contains JavaScript, plus the repo
#   root for root-level files (cli.js). An EXTENSION is declared in EXTS
#   below, not globbed implicitly, so a new file type cannot slip in
#   unannounced: adding one means adding it here, where the count prints.
#
# VACUITY GUARD
#   Per-extension totals are checked for non-zero. An empty `find` (wrong
#   root, renamed directory, bad glob) must never read as "all clear".
#   Declared roots are also required to exist, so a rename cannot silently
#   drop a whole subtree from coverage.
#
# Usage: npm run check  OR  bash scripts/check.sh  OR  ./scripts/check.sh
# Zero dependencies: bash + node only.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."

# Directories to scan. Must all exist (checked below).
ROOTS=(src scripts tests bench evolve)
# File extensions to syntax-check. The vacuity guard is per entry.
EXTS=(.js .mjs)
# The repo root itself, non-recursively, for root-level entry points.
ROOT_FILES=1

for r in "${ROOTS[@]}"; do
  if [ ! -d "$r" ]; then
    echo "[check] FAIL: declared root '$r' does not exist — coverage is silently reduced" >&2
    exit 1
  fi
done

declare -A EXT_COUNT=()
declare -A EXT_FAILED=()
for e in "${EXTS[@]}"; do
  EXT_COUNT["$e"]=0
  EXT_FAILED["$e"]=0
done

FAILED=0
TOTAL=0

# LC_ALL=C sort for a stable, locale-independent file order. Repo paths
# contain no spaces/newlines, so newline-delimited find output is safe.
for e in "${EXTS[@]}"; do
  for r in "${ROOTS[@]}"; do
    while IFS= read -r f; do
      TOTAL=$((TOTAL + 1))
      EXT_COUNT["$e"]=$((EXT_COUNT["$e"] + 1))
      if ! node --check "$f"; then
        echo "[check] SYNTAX ERROR: $f" >&2
        FAILED=$((FAILED + 1))
        EXT_FAILED["$e"]=$((EXT_FAILED["$e"] + 1))
      fi
    done < <(find "$r" -type f -name "*${e}" -not -path '*/node_modules/*' | LC_ALL=C sort)
  done
  if [ "$ROOT_FILES" -eq 1 ]; then
    while IFS= read -r f; do
      TOTAL=$((TOTAL + 1))
      EXT_COUNT["$e"]=$((EXT_COUNT["$e"] + 1))
      if ! node --check "$f"; then
        echo "[check] SYNTAX ERROR: $f" >&2
        FAILED=$((FAILED + 1))
        EXT_FAILED["$e"]=$((EXT_FAILED["$e"] + 1))
      fi
    done < <(find . -maxdepth 1 -type f -name "*${e}" -not -path './node_modules/*' | LC_ALL=C sort)
  fi
done

# Vacuity guard, per extension: an empty result set is a broken gate, not a
# clean bill of health.
SUMMARY=""
SEP=""
for e in "${EXTS[@]}"; do
  n="${EXT_COUNT[$e]}"
  f="${EXT_FAILED[$e]}"
  SUMMARY="${SUMMARY}${SEP}${e}: ${n} checked, ${f} failed"
  SEP=" | "
  if [ "$n" -eq 0 ]; then
    echo "[check] FAIL: found 0 ${e} files under ${ROOTS[*]} — that extension's gate would be vacuous" >&2
    exit 1
  fi
done

if [ "$FAILED" -ne 0 ]; then
  echo "[check] FAIL: $FAILED of $TOTAL file(s) failed node --check" >&2
  echo "[check]   per-extension: $SUMMARY" >&2
  exit 1
fi

# Version-drift guard: package.json "version" must equal src/version.js
# VERSION. These are the only two version carriers; a mismatch means one
# was bumped without the other.
SRC_VERSION=$(sed -n 's/^export const VERSION = "\(.*\)";$/\1/p' src/version.js)
PKG_VERSION=$(sed -n 's/^[[:space:]]*"version":[[:space:]]*"\([^"]*\)",*$/\1/p' package.json | head -1)
if [ -z "$SRC_VERSION" ] || [ "$SRC_VERSION" != "$PKG_VERSION" ]; then
  echo "[check] FAIL: version drift — src/version.js='$SRC_VERSION' package.json='$PKG_VERSION'" >&2
  exit 1
fi

# Cache-bust guard: zero hand-typed version query literals in shipped code.
#
# HISTORY (recorded because both halves of this were silent)
#   Every module edge used to carry a literal cache-bust suffix —
#   `import { ForceField } from "./forcefield.js?v=10"` — 121 of them across
#   42 src/ files plus index.html's module tag, frozen at `10` while
#   src/version.js said `1.1.0-fp7`. A module edited without also editing
#   its neighbours' literals kept the same url, the browser served the cached
#   bytes, and the symptom was "my fix did nothing". Freshness is now sw.js's
#   job (network-first, cache:"no-store"), so a literal here is at best
#   redundant and at worst a url that pins a cached copy forever.
#
#   tests/test_cache_contract.js is the authoritative check for the same
#   invariant (it also runs sw.js in a sandbox and checks the single ui.js
#   instance). This is the cheap pre-commit copy, so `npm run check` is red
#   before anyone reaches for `npm test`.
CACHE_LITS=$(grep -rnE '\.js\?|\?v=' src index.html 2>/dev/null || true)
if [ -n "$CACHE_LITS" ]; then
  echo "[check] FAIL: version query literal in shipped code — freshness is sw.js's job (see sw.js)" >&2
  echo "$CACHE_LITS" | head -20 | sed 's/^/[check]   /' >&2
  echo "[check]   (tests/test_cache_contract.js is the full check; ?v= appears in prose only after a comment strip)" >&2
  exit 1
fi
if [ ! -s sw.js ]; then
  echo "[check] FAIL: sw.js missing or empty — nothing forces a fresh module graph on reload" >&2
  exit 1
fi

echo "[check] PASS: $TOTAL files checked, 0 syntax errors"
echo "[check]   per-extension: $SUMMARY"
echo "[check]   roots: ${ROOTS[*]} (+ repo root)"
echo "[check]   VERSION=$SRC_VERSION"
