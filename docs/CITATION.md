# Citation — Code as Supplement Zenodo (J95)

> **Tag:** `v1.0-jpcb` · **Zenodo DOI placeholder:** `10.5281/zenodo.XXXXXXX`
> Replace `XXXXXXX` on deposition. Cited in `manuscript/manuscript_jpcb.tex`
> (§ Code and Data Availability) and in `CITATION.cff`.

## How to cite

### Software (CITATION.cff)

`CITATION.cff:1` `cff-version: 1.2.0`, `CITATION.cff:13` `version: v1.0-jpcb`,
`CITATION.cff:24` `10.5281/zenodo.XXXXXXX` (Zenodo DOI placeholder).

```bibtex
@software{simulation_coarse_v1_0_jpcb,
  title        = {simulation\_coarse --- Client-Side Multi-Scale Biomolecular Dynamics and Free Energy Landscapes},
  version      = {v1.0-jpcb},
  author       = {Rana, Samir},
  year         = {2026},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.XXXXXXX},
  url          = {https://doi.org/10.5281/zenodo.XXXXXXX},
  note         = {MIT License. Tag v1.0-jpcb. Replace XXXXXXX on Zenodo deposition.}
}
```

### Manuscript + code supplement

In `manuscript/manuscript_jpcb.tex:359` Code and Data Availability, the
repository is cited and the Zenodo archive is referenced as the supplement:

> All source code ... are open-source under the MIT License at
> `https://github.com/samirr/simulation_coarse` and archived as
> `v1.0-jpcb` on Zenodo (`10.5281/zenodo.XXXXXXX`, placeholder until
> deposition). See `CITATION.cff` and `docs/CITATION.md`.

## Archiving checklist (do not publish before completing)

- [ ] `git tag v1.0-jpcb && git push origin v1.0-jpcb`
- [ ] GitHub release `v1.0-jpcb` triggers Zenodo-GitHub webhook (or manual upload to Zenodo)
- [ ] Replace `XXXXXXX` in `CITATION.cff:24` and this file with real Zenodo record DOI
- [ ] Update `manuscript/manuscript_jpcb.tex` DOI placeholder in the next revision
- [ ] Verify `curl -L https://doi.org/10.5281/zenodo.XXXXXXX` resolves

## Verification

```bash
grep -n "v1.0-jpcb" CITATION.cff docs/CITATION.md    # hits
grep -n "zenodo" CITATION.cff docs/CITATION.md -i   # hits
grep -n "Zenodo" manuscript/manuscript_jpcb.tex       # hits (after patch)
```

*See also `README.md:59` Code and Data Availability and `ROADMAP.md:1` for v1 scope.*

