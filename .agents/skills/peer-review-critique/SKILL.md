---
name: peer-review-critique
description: Two-reviewer adversarial critique protocol for peer-reviewing computational physical chemistry and biomolecular modeling manuscripts prior to journal submission.
---

# Two-Reviewer Adversarial Critique Protocol

## Reviewer Profiles

### Reviewer 1 (Theoretical Physical Chemist / Statistical Mechanics Specialist)
- **Focus Areas**:
  - Symplectic measure preservation and detailed balance in Langevin dynamics.
  - Correctness of partition function integrations and volume corrections ($\Delta G^\circ_{\rm bind}$).
  - Thermodynamic consistency of concentration-dependent chemical potentials ($\mu$).
  - Mathematical completeness of continuous SDEs and discrete master equations.
  - Detection of mathematical hand-waving, unverified claims, and dimensional inconsistencies.

### Reviewer 2 (Biophysical Chemist / Structural & Computational Chemist)
- **Focus Areas**:
  - Force field parameter choices (AMBER charges, Born radii, dielectric constants, SASA tensions).
  - Accuracy and experimental alignment of molecular benchmarks (PDB 4W52, T4 Lysozyme, small molecule binding affinities).
  - Practical usability, code availability, reproducible test suites, and hardware scaling.
  - Figure readability, axis labeling, data-to-ink ratio, and clarity of biophysical insights.
  - Elimination of redundancy, overlap, and ungrounded assertions.

## Review Loop Protocol
At each stage of manuscript refinement:
1. Reviewer 1 and Reviewer 2 independently identify 2–3 critical weaknesses, overlapping paragraphs, or imprecise formulations.
2. Authors draft direct, surgical revisions addressing each critique.
3. Both reviewers verify that the applied revision resolves the technical objection without introducing new inconsistencies.
