# Pattern P4 — Feature parity ≠ feature equality (progressive disclosure)

**Root cause:** Every feature added gets equal default visual weight
(7 numbered `<details open>` panels, 6 toolbar checkboxes, always-open
equations panel), producing the "AI slop" dashboard feel: dense, unordered,
no hierarchy of importance.

**Evidence:** index.html:24-245 — seven panels all `open`; index.html:253-261
toolbar with 6 simultaneous toggles; index.html:264-275 equations always
rendered. CSS offers no collapsed/quiet state (css/style.css:54-66).

**Fix that works:** One visible core flow (Load → Build → Run → Record),
analysis + kinetics collapsed by default, advanced physics/ML behind a
single "Advanced" disclosure, equations removed from default view. The full
feature set stays reachable — it is never deleted, only ranked.

**Generalization:** Decluttering is ranking, not removal. The anti-slop test:
a first-time user should be able to state the primary action within 5
seconds; if not, the hierarchy has failed.
