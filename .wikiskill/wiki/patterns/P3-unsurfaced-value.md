# Pattern P3 — Unsurfaced work has zero perceived value

**Root cause:** Real, measured improvements (dashboard numbers) locked in
developer artifacts (TRANSFORMATION_PLAN_100.md, JSON) never reach the
running product, so users cannot answer "what actually improved?".

**Evidence:** TRANSFORMATION_DASHBOARD.json holds: CG 0.20 ms/step, Heavy
14.7 ms/step, forces FD maxRel 7.4e-6, NVE drift 0.05%, scale exponent 2.06,
B-factor median R 0.27 (target 0.45, honest miss). No UI surface or
root-level human changelog presents these as a before/after.

**Fix that works:** A short, honest, numbers-first changelog (WHAT_CHANGED.md)
linked from the app header; separates [measured] from [docs/infra] from
[honest misses]. Never claim aspirational stubs as wins (GPU R>0.999 is
aspirational; barriers are heuristic).

**Generalization:** For every "make it better" batch, ship a diff-of-value
next to the diff-of-code, or the work is invisible. Honesty about misses
increases trust more than the wins do.
