# WikiSkill Runtime — adapted for simulation_coarse

Instance of the three-layer co-evolution framework from:
Tang et al., "WikiSkill: Compiling Agent Experience into Persistent Knowledge
for Skill Evolution", arXiv:2608.27454 (2026).

    raw/   (immutable)  execution traces — write-only, never rolled back
    wiki/  (persistent) distilled patterns + skill-impact log — never rolled back
    skills ← .agents/skills/ — executable skills, gated by validation

## The loop (run once per work session)

1. **Execute** tasks with the current skills. Save the session trace to
   `raw/traces/trace-<date>-<slug>.md` (append-only; include what was
   attempted, what failed, what the tests said — verbatim numbers).
   The Inference Agent does NOT read wiki/ during execution (paper finding:
   blind execution yields diagnostic failures).
2. **Maintain** (Wiki Maintainer role): read the new trace against
   `wiki/patterns/*.md`; add new pattern pages or patch existing ones with
   root causes, evidence, and the fix that worked. Never delete history —
   supersede by linking.
3. **Propose** (Skill Proposer role): from the updated wiki, propose skill
   edits (create/edit/delete) in `.agents/skills/`. Check
   `.wikiskill/wiki/evolution/skill-impact.md` first to avoid re-proposing rejected ideas.
4. **Gate** (Gating & Rollback): run `node scripts/wikiskill_gate.js`.
   Accept the skill set only if the measured validation score beats the
   stored best `R_best`; otherwise roll back skills — but ALWAYS keep the
   wiki update and log the rejection in `.wikiskill/wiki/evolution/skill-impact.md`.

## Local conventions

- Validation gate = `node tests/test_all.js` (32 assertions) plus
  `node --check` over changed src files; extended gate runs the full
  `tests/test_*.js` suite and benches.
- A pattern page MUST cite `file:line` evidence and at least one measured
  number (test output, bench ms, R value). No unsourced claims.
- Skills are validated machine-side (`node --check`, tests) — never trust
  untested prose as a skill.
- Impact log verdicts: ACCEPTED / REJECTED / ROLLED_BACK, each with the
  gating numbers that decided it.
