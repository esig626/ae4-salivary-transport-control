# Results ledger

Only validated results belong in this file.

For each result record

- result identifier
- concise statement
- status as reproduced or new
- analysis path
- code commit
- input data or parameter set
- numerical tolerances
- figure or table pointer
- robustness status
- caveats

## Phase 10 validated results

| ID | Statement | Status | Analysis / code | Input and tolerance | Pointer | Robustness and caveats |
| --- | --- | --- | --- | --- | --- | --- |
| R10-01 | At a regular steady state, `D_theta u*=-F_u^-1 F_theta` and `J_h=-h_u F_u^-1 F_theta` | New derivation for this project; standard theorem | `analysis/10_identifiability_discrimination/theory.md`; Phase 10 branch commit | `C1` maps, nonsingular `F_u` | Theory, Proposition 1 setup | Applies only on a regular local branch |
| R10-02 | A regular scalar observation of two activities has a local one-dimensional equivalence set; a rank-two panel is locally identifying | New derivation for this project; standard theorem | `analysis/10_identifiability_discrimination/theory.md`; Phase 10 branch commit | Constant/regular rank assumptions | Theory, Propositions 1--2 | No global-identifiability claim |
| R10-03 | AE2/AE4 panel rank equals the rank of their projected stoichiometric signatures when both activity factors are nonzero | New organization for this project | `analysis/10_identifiability_discrimination/theory.md`; Phase 10 branch commit | Published cycle signatures; regular `F_u` | Theory, Proposition 3 | The criterion is general; physiological numerical rank remains unknown |
| R10-04 | The affine diagnostic implicit sensitivity matches central differences | Reproduced computational check | `src/identifiability_discrimination/`; Phase 10 branch commit | Published concentration anchor; diagnostic `A,S`; step `1e-6`; acceptance `2e-10` absolute | `results/10_identifiability_discrimination/derivative_verification.json` | Maximum error `5.43e-11`; diagnostic is not a physiological reproduction |
| R10-05 | A consistent complete residual-row sign change leaves implicit sensitivity invariant | New derivation and reproduced check | `analysis/10_identifiability_discrimination/theory.md`; Phase 10 branch commit | Diagonal signs `+/-1`; acceptance `1e-14` | `derivative_verification.json` | Does not resolve physical water-flux orientation |
| R10-06 | The published 2018 numerical steady state is not independently reproducible from the citable specification | Reproduction audit result | `analysis/10_identifiability_discrimination/reduction.md`; Phase 10 branch commit | Phase 00 source map and published appendices | `results/10_identifiability_discrimination/baseline_check.json` | Baseline anchor passes by construction; `full_published_model_reproduced=false` |
