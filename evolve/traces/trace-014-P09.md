# Evolution 14 — P09 (mid)

**Goal:** Cut src LOC below the bloat knee: 23403 LOC vs the 20000 penalty start
**Verdict:** ACCEPTED
**Note:** fps 30 measured by nothing, ever; now null with _meta.measured=false, plus 9 more published fps claims the planner missed

**Acceptance:** src LOC (sum of lines over src/**/*.js as counted by evolve.mjs readTrackedStats) <= 20000, or the bloat LOC term is replaced by a term with an argued-for justification. `node evolve/evolve.mjs gate` reports bloat >= 0.8353.

**Touched:** src/, tests/, evolve/reports/
