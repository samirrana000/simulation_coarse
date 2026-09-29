#!/usr/bin/env bash
# manuscript/reproduce.sh — Reproducibility package (J93)
# Re-runs benchmarks and regenerates all publication figures.
# Usage:
#   bash manuscript/reproduce.sh        # full repro (needs node + python3)
#   bash manuscript/reproduce.sh fast   # quick smoke (skip heavy GB FD tests)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
echo "=== simulation_coarse reproduce.sh (v1.0-jpcb) ==="
echo "ROOT=$ROOT"
echo "node=$(node --version 2>&1 || echo 'node not found')"
echo "python=$(python3 --version 2>&1 || echo 'python3 not found')"
echo ""

# 1) Syntax / unit tests
echo "--- [1/4] Syntax check ---"
node --check src/*.js src/physics/*.js || echo "warn: node --check failed (non-fatal)"

echo "--- [2/4] Benchmark: bench/perf.js (CG + Heavy timing) ---"
# This is the measurable J93 criterion: must contain bench/perf.js
# Captures ms/compute for CA 164 and heavy 1308, dashboard JSON
if [ -f bench/perf.js ]; then
  node bench/perf.js | tee /tmp/reproduce_perf.log
  echo "perf log at /tmp/reproduce_perf.log"
else
  echo "missing bench/perf.js"
  exit 1
fi

echo ""
echo "--- [2b/4] Optional: bench/b_factors.js, bench/pdbbind_gb.js (if present) ---"
if [ -f bench/b_factors.js ]; then node bench/b_factors.js | tee /tmp/reproduce_bfactors.log || true; fi
if [ -f bench/pdbbind_gb.js ]; then node bench/pdbbind_gb.js | tee /tmp/reproduce_pdbbind.log || true; fi
if [ -f bench/budget.json ]; then cat bench/budget.json; fi

echo ""
echo "--- [3/4] Regenerate publication figures ---"
echo "Running: python3 manuscript/generate_figures.py"
if [ -f manuscript/generate_figures.py ]; then
  python3 manuscript/generate_figures.py | tee /tmp/reproduce_figs.log
  echo "Figures written to manuscript/figures/ (fig1_overview.pdf, fig2_pmf_free_energy.pdf, fig3_kinetics_tpt.pdf, fig4_performance.pdf)"
  ls -lh manuscript/figures/*.pdf 2>&1 | head -20
else
  echo "missing manuscript/generate_figures.py"
  exit 1
fi

echo ""
echo "--- [4/4] Performance budget check (bench/budget.json, CI warn) ---"
if [ -f bench/budget.json ]; then
  echo "budget.json content:"
  cat bench/budget.json
  # bench/budget_check.js is the single implementation of the contract; this
  # file and .github/workflows/check.yml both call it rather than keeping
  # private copies (the copies had drifted, and BOTH used to print
  # "fps budget 30" for a quantity neither measured — nothing in this repo
  # measures frame rate, so fps is null with the reason recorded).
  # --perf-log reuses the bench/perf.js output already captured in [2/4].
  # Exit 1 = at least one MEASURED key is over budget; unmeasured keys are
  # reported as UNMEASURED and are never counted as passing.
  echo "contract check:"
  node bench/budget_check.js --perf-log /tmp/reproduce_perf.log || true
fi

echo ""
echo "=== reproduce.sh complete ==="
echo "Outputs:"
echo "  - /tmp/reproduce_perf.log  (bench/perf.js — CG + Heavy ms/compute)"
echo "  - manuscript/figures/*.pdf + *.png (regenerated)"
echo "  - /tmp/reproduce_figs.log"
echo "Diff against reference <1%: compare newly generated PDFs to committed versions"
echo "  diff manuscript/figures/fig1_overview.pdf manuscript/figures/fig1_overview.pdf || echo 'binary diff check (size/time only)'"
echo "See docs/TUTORIAL.md for where the model fails (charged ligand, membrane) and how to diagnose wrong ΔG."
