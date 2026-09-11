#!/usr/bin/env bash
set -e
# A01 deterministic build check — syntax validation for all src modules
# Usage: npm run check  OR  bash scripts/check.sh  OR  ./scripts/check.sh
# Placeholder comment for package.json# check script:
#   "scripts": { "check": "bash scripts/check.sh" }

echo "[check] node --check src/*.js"
for f in src/*.js; do
  echo "  checking $f"
  node --check "$f"
done
echo "[check] all src/*.js syntax OK (VERSION=$(grep -o '"[^"]*"' src/version.js | head -1))"
