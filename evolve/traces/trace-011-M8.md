# Evolution 11 — M8 (mid)

**Goal:** Kill the ?v=10 cache-bust sprawl by making the version real
**Verdict:** ACCEPTED
**Note:** 121 ?v= literals removed via sw.js no-store; 27.8h heuristic-freshness staleness eliminated; 64-module graph provably unchanged

**Acceptance:** Either zero hardcoded ?v= literals in src (injected from version.js), or `npm run check` fails when any ?v= mismatches src/version.js — and is wired into CI. Version bump is one edit. Verified by: edit a module, reload, see the change.

**Touched:** src/, scripts/check.sh, index.html, .github/workflows/check.yml
