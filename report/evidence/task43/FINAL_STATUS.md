# Task 43 final status

Complete provenance audit and standalone report, on the separate, unmerged branch `analysis/task-43-parameter-provenance`.

## Exact source attribution

- Frozen production base: `c4d4f207f702d8eee09ab94d2d90ae90552dd641`.
- Recovered Task 43 checkpoint: `590fefd5ea3c32bd52192870e9dce0095928c964`; this is the predecessor, not the final report commit.
- Verified Task 44 numerical source: `f783785a863440df459e9c5530beba4c7eb59110` on `analysis/task-44-full-system-mathematics`. The crosswalk reads this exact Git commit and records seven source-file SHA256 hashes. Later report-only Task 44 commits do not change this attribution.
- The final Task 43 commit is the commit containing this status. Its report and inventory bytes are identified below, avoiding a self-referential head SHA.

## Completed and verified

The inventory contains 137 records, 109 active flags, five inactive legacy fields, twelve constructor seeds and thirteen frozen state coordinates. Every active record includes units, equation role, selected-law use, source location, provenance, constraints, limitations, uncertainty status and manuscript dependence. All actual production dataclass fields and genotype interventions are covered. These counts are not independent fitted degrees of freedom.

Fresh reconstruction reproduces the saved inventory and numerical-control outputs exactly. Verification checks the original parameter/state hashes, 75 unchanged production source files, all required record fields and generated report tables. The crosswalk reproduces 16 stationary parameter records and 16 inverse observable rows (eight parameters, two integration observables) directly from the pinned companion commit. Its companion checks include 128 independently solved nearby WT/null equilibria and independently integrated inverse cases. Task 43 performs no new sensitivity simulations or parameter fits.

The 11-page standalone report compiles without unresolved citations, references or TeX warnings. Rendered pages were visually checked for clipping, overlap and table legibility.

Reproduce from a checkout with the companion commit available in fetched Git history:

```bash
python analysis/43_parameter_provenance/verify.py
bash analysis/43_parameter_provenance/build.sh
```

## What is justified and what remains assumed

Published coefficients document a kinetic/model lineage; they do not measure salivary carrier abundance. Resting-state, charge and NBC capacities have conditional derivations. In particular, the 70/30 loading partition remains assumed. No calibrated physiological uncertainty interval is claimed. Logical domains and numerical perturbation steps are not confidence intervals. Task 41 remains a target-selected construction.

The most consequential uncertain settings for the article are pump capacity/distribution, NBC capacity and saturation width, NKCC effective capacity, NHE carrier amount and stimulated multiplier, the finite buffer pool, and CaCC recruitment. Priorities are functional pump and ion-conditioned NKCC measurements; stimulated bicarbonate/alkalinity transport and acid-extrusion measurements; buffer/TIC/volume time series; and matched CaCC currents during acute AE4 loss. Global robustness remains unresolved without defensible parameter intervals.

## Corrections and scientific limits

- Fixed Cha Mod2 exponents, genotype settings and frozen initial values are now explicit. Nominal CCh/IPR doses are inactive input metadata. Beta occupancy is dimensionless; Task 41 requires `0 < b <= 1`; the numerical integration starts at `1e-6` s.
- The stationary AE4-null secretion deficit is 0.94617%, distinct from frozen Task 40's 600 s cumulative 3.8575% result. Exact equal-routing cancellation does not establish zero total AE4 sensitivity.
- NBC supplies 95.48% of the analysed AE4 alkalinity demand. The conditional no-NBC capacity bound covers 57.96% of that specified demand. A physiological no-NBC equilibrium nevertheless survives (pH 6.96833; secretion 0.00153598 pL/s, 3.035% lower), with AE4 loading reduced from 0.0750724 to 0.00632651 fmol/s. A broad claim of universal NBC necessity for secretion must therefore be qualified.
- Slow local modes already exist: approximately 346 s in WT and 460 s in the null. They must be interpreted with state/output projection, not dismissed as absent.
- The first crossing found in the prescribed Task 41 numerical scan needs approximately 89.51% AE4-dependent stimulated CaCC recruitment. No global lower bound or exclusion of unsampled earlier crossings is proved. The 30.3% comparison is mean minus reported SEM, not a confidence bound. Adaptive and sampled observables remain separate. Continued-stimulus null cases later exceed the pH gate at about 1485 s; their declared 600 s results are unchanged.
- Strong modelled NKCC compensation remains in tension with the isolated 2015 assay. The structural NBC pathway is not assigned to a specific isoform.

No model source, Tasks 1-42, frozen outputs, `main` or manuscript files were modified. Neither analysis branch was merged.

## Deliverables and SHA256

The inventory, audit summary, numerical controls, source manifest, verification, sensitivity crosswalk and generated TeX tables are in `analysis/43_parameter_provenance/output/`. The report files are `analysis/43_parameter_provenance/report.tex` and `analysis/43_parameter_provenance/AE4_parameter_provenance.pdf`.

- `report.tex`: `173c8855813b6ff124cae9bbe192168c6a3abc88e6ceae65cb2d8546f08100d0`
- `AE4_parameter_provenance.pdf`: `7d3018ab4b720bf998c0a46f7ffbf18650db314919be21c0af05fa776d3f9643`
- `output/parameter_inventory.json`: `1e3d5763a171dfd64b95fbb72a8a35eae4bbd3fcd1a61939c2886808296d994b`
- `output/task44_sensitivity_crosswalk.json`: `cbe0570de39fe11b64199b0d340befa8789126cf9e9d734a37f743cd834147ef`
