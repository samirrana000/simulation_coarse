#!/usr/bin/env bash
# check.sh — A01 deterministic syntax gate.
#
# Recursively runs `node --check` on EVERY .js module under src/ and exits
# non-zero if any file fails. Replaces the old `node --check src/*.js
# src/physics/*.js` gate, which was broken two ways:
#   1. it never descended into src/chem, src/compute, src/analysis,
#      src/capture or src/physics/{forcefield,solvation,integrators};
#   2. `node --check` accepts exactly ONE file argument — any further
#      arguments land in process.argv and are silently ignored, so the
#      whole command only ever checked src/analysis-panel.js (1 of 67).
#
# Usage: npm run check  OR  bash scripts/check.sh  OR  ./scripts/check.sh
# Zero dependencies: bash + node only.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."

COUNT=0
FAILED=0

# LC_ALL=C sort for a stable, locale-independent file order. Repo paths
# contain no spaces/newlines, so newline-delimited find output is safe.
while IFS= read -r f; do
  COUNT=$((COUNT + 1))
  if ! node --check "$f"; then
    echo "[check] SYNTAX ERROR: $f" >&2
    FAILED=$((FAILED + 1))
  fi
done < <(find src -type f -name '*.js' | LC_ALL=C sort)

if [ "$COUNT" -eq 0 ]; then
  echo "[check] FAIL: found 0 .js files under src/ — the gate would be vacuous" >&2
  exit 1
fi

if [ "$FAILED" -ne 0 ]; then
  echo "[check] FAIL: $FAILED of $COUNT file(s) failed node --check" >&2
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

echo "[check] PASS: $COUNT files checked, 0 syntax errors (VERSION=$SRC_VERSION)"
