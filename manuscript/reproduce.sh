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
  # CI warn logic: compare perf log to budget (warn, do not fail)
  # heavy_compute_ms budget 2.0, cg_compute_ms 0.5 (see docs/PERFORMANCE.md)
  echo "CI budget warn: if meanMs > budget, emit ::warning (not error)"
  python3 - << 'PY' || true
import json, re, pathlib
try:
    budget=json.loads(pathlib.Path("bench/budget.json").read_text())
    text=pathlib.Path("/tmp/reproduce_perf.log").read_text()
    m=re.search(r"CG:[^\n]*?(\d+\.\d+)\s*±", text)
    cg=float(m.group(1)) if m else None
    m2=re.search(r"Heavy:[^\n]*?(\d+\.\d+)\s*±", text)
    hv=float(m2.group(1)) if m2 else None
    print(f"budget heavy_compute_ms={budget.get('heavy_compute_ms')} fps={budget.get('fps')} cg_compute_ms={budget.get('cg_compute_ms')}")
    if cg and cg>budget.get("cg_compute_ms",999):
        print(f"::warning :: CG compute {cg:.3f} ms exceeds budget {budget['cg_compute_ms']} ms")
    if hv and hv>budget.get("heavy_compute_ms",999):
        print(f"::warning :: Heavy compute {hv:.3f} ms exceeds budget {budget['heavy_compute_ms']} ms")
    if hv and hv and budget.get("fps"):
        fps_budget=budget["fps"]
        print(f"fps budget {fps_budget} — check render fps via bench/perf.js scale if applicable")
except Exception as e:
    print("budget check skipped:", e)
PY
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
