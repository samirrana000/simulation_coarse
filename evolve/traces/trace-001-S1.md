# Evolution 1 — S1 (short)

**Goal:** Fix the silently-broken CI (npm test with no package.json)
**Verdict:** ACCEPTED
**Note:** syntax gate 1/67 -> 67/67, provably fails on planted error; npm test 352/0 exit 0

**Acceptance:** package.json exists with scripts.test -> tests/test_all.js; CI syntax-checks ALL src/**/*.js recursively; `npm test` exits 0; green locally.

**Touched:** package.json, .github/workflows/check.yml, scripts/check.sh
