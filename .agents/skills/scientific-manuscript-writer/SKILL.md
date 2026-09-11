---
name: scientific-manuscript-writer
description: Standard operating procedure and stylistic rules for writing high-impact, rigorous Physical Chemistry & Biophysics manuscripts for ACS journals (JPCB/JCTC/JACS) without AI clichés, repetitive boilerplate, or superficial hand-waving.
---

# Scientific Manuscript Writing Protocol (Physical Chemistry & Biophysics)

## 1. Zero Tolerance for "AI Slop" & Cliché Tropes
- **Banned Words & Phrases**: Never use "delve into", "testament to", "tapestry", "revolutionize", "game-changer", "beacon", "in conclusion, it is evident", "serves as a crucial", "pivotal role", "paradigm shift" (unless historically justified), "seamless integration", "synergy", or repetitive summarizing transitions.
- **Direct Voice**: Write in active, precise scientific voice. State equations, state parameters, state numerical outcomes with uncertainty bounds ($\pm \sigma$), and explain physical mechanisms directly.
- **No Hand-Waving**: Every theoretical claim must be backed by an exact mathematical equation, Hamiltonian component, or explicit algorithm step.

## 2. JPCB Manuscript Structure & Density
1. **Title**: Precise, descriptive, avoiding hype.
2. **Abstract**: 150–250 words. Quantitative, self-contained: (i) Biophysical problem, (ii) Exact methodology, (iii) Key quantitative result with numerical comparison to experiment, (iv) Broader implication.
3. **Introduction**:
   - Trace the physical dilemma: Why standard atomistic MD is slow, why standard coarse-grained networks lack chemical specificity, why existing server tools create computational bottlenecks.
   - Establish the exact gap this work fills.
4. **Theoretical Models and Computational Methods**:
   - Complete Hamiltonian breakdown.
   - SDE Integrator & symplectic operator splitting.
   - Solvation continuum & pairwise descreening integrals (HCT/GB).
   - Enhanced sampling & free energy partition function integration.
   - Master equations & Transition Path Theory committors.
5. **Results and Discussion**:
   - Concrete benchmark on model systems (e.g. T4 Lysozyme L99A / Benzene).
   - Thermal fluctuations and correlation with crystallographic B-factors.
   - PMF profiles and standard-state binding free energy ($\Delta G^\circ_{\rm bind}$) vs ITC experimental data.
   - Reactive flux and committor analysis.
   - Hardware performance (WebGPU vs Workers vs single-thread).
6. **Conclusions**: Direct summary of findings and open physical questions.
7. **Code & Data Availability**: Explicit URLs, build commands, and artifact locations.

## 3. Strict Deduplication Rule
- Never repeat the same explanation in the Introduction, Theory, and Results.
- If a method is derived in Section 2, Section 3 must discuss only the observed physical data and comparative interpretation.
