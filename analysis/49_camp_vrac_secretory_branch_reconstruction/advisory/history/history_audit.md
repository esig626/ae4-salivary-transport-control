# Advisory Task 49 repository-history audit

Pinned source: `2bf53b773c14f3c49f539957ea3f1964bc8d227f`. Read `AGENTS.md` first. No canonical modification, model run, inference, test rerun, commit or push. Scientific acceptance belongs to the orchestrator.

## Immediate literal defect: VRAC was tested under Vbeta

`analysis/49_camp_vrac_secretory_branch_reconstruction/EXCLUSION_LEDGER.md:121` is false as written. The pinned tree already contains:

- `src/modern_full_model/vbeta_diagnostic.py`: V1 conductance `g_V beta max(V_i/V_rest-1,0)`; `vbeta_effective_conductance_S` at lines 157–174; sign/conservation-consistent fixed-voltage current at lines 177–245; exact-rest anti-deadlock assessment at lines 248 onward. `V1_PRODUCTION_ELIGIBLE=False` (line 67).
- `src/modern_full_model/run_vbeta_diagnostic.py`: diagnostic runner on historical native WT roots.
- `tests/test_modern_vbeta_diagnostic.py`: exact-off nesting, swelling gate, chloride/current conservation, deadlock, and ineligible-exogenous-swelling checks.
- `analysis/13B_modern_full_model/vbeta_diagnostic.md`: title explicitly “Vbeta/VRAC pre-reveal diagnostic”; source is the same Catalán PNAS 2015 article.
- `analysis/13B_modern_full_model/evidence_ledger.md:266` onward: IPR/TMEM16A independence, blockers DCPIB/NPPB, slow 12.5 ± 0.2% swelling and rapid CCh shrinkage. This is a historical extraction, not a replacement for Task49 primary-source verification.
- `results/13B_modern_full_model/vbeta_native_exact_rest_summary.json`: final saved numerical record.

The final JSON records **389 roots / 778 solver rows, 69 roots with positive gate, 10 WT-rest-passing roots, zero WT-rest-passing activated roots**. All solvers succeeded and gate classifications agreed. Maximum relative swelling anywhere was 0.0187221128, but no WT-valid root activated. Capacity status remains `UNLICENSED_NO_WHOLE_CELL_CURRENT_MAP`; production eligibility is false. The superseded reference root shrank 1.3000→1.273730 pL and had identically zero gate for every nonnegative conductance.

The diagnostic note still describes an interim 35-root audit; use the final JSON for counts. Preserve both historical artefacts rather than silently amending frozen Task13B.

**Precise replacement claim:** the source-backed branch is absent from active production equations, but a Vbeta/VRAC swelling-gate diagnostic was already implemented and tested. The test was isolated and never inserted into production voltage closure. Its historical rejection is scoped to those roots and to absent amplitude identification, not a global exclusion on the later Task48 chassis. Task49 should explicitly carry the deadlock/identifiability result forward, without repeating the old panel or representing V1 as wholly new.

## Exclusion verification

| Task | Verified repository result and qualification | Reference |
|---|---|---|
| 12 | C1–C8/rest-compatible positive-PKA variants failed; AE4-null vector field independent of AE4 kinetics. Do not promote finite tested-family exclusion to a global impossibility theorem. | `analysis/12_ae4_mechanism_reconstruction/final_answer.md` |
| 13C | 0.10/0.25/0.50 uM calcium panel did not meet absolute WT flow requirement; calcium escalation is not an AE4-specific correction. | `analysis/13C_calcium_fast_screen/final_answer.md` |
| 19 | 50% pump cases oppose harmful AE4-pump interaction. Task classification is NUMERICALLY INCONCLUSIVE because 10% pump failed resting continuation; it does not prove no possible pump mechanism. | `analysis/19_ae4_loss_nak_pump_coupling/final_answer.md` |
| 31 | Cha NHE implemented; WT Cl 60.305 mM misses data; two exact-null resting attempts did not close within limits. This is bounded numerical failure, not proof no root exists. | `analysis/31_nhe1_mechanistic_repair/final_answer.md` |
| 36 | NBC is explicitly a stimulus-recruited increment; zero basal gate chosen to nest Task31 rest, not measured native absence. | `analysis/36_minimal_nahco3_alkalinity/design.md` |
| 37 | Valid stimulated WT; integrated positive Cl loading 79.9248% NKCC / 20.0752% AE4 / 0% AE2. | `analysis/37_wt_nbc_validation/final_answer.md` |
| 38 | Acute shared-WT-rest AE4 null retains 98.1605% WT cumulative flow; compensation remains large. | `analysis/38_ae4_perturbation_validation/final_answer.md` |
| 39 | Source-fixed Palk/Benjamin core increases null NKCC compensation to 33.1013%; null deficit 0.627320%. Core accepts only intracellular substrates. | `analysis/39_palk_nkcc1_full_validation/final_answer.md` |
| 40 | Equal routing null deficit 3.857486%; NKCC replaces about 88.51% of lost AE4 loading. | `analysis/40_ae4_equal_cation_routing/final_answer.md` |
| 41 | AE4-dependent stimulated CaCC recruitment selected against secretion target; 30.26% loss is not validation; temporal pattern and mutant Cl/pH/volume fail. | `analysis/41_ae4_loss_algebraic_design/final_answer.md` |
| 42 | Seven source classes; no held-out success under inherited scalar law. Three carbonate classes fail WT pH. Preserves source-class/kinetics qualification. | `analysis/42_catalan_2025_ae4_mechanism_classes/final_answer.md` |
| 46 | Full reference export remains feasible with AE4 absent; isolated NKCC initial-rate assay cannot justify an arbitrary sustained cap. | `analysis/46_physiology_constrained_model_reconstruction/DYNAMIC_RECONSTRUCTION_HANDOFF.md` |
| 47 | K efflux +18.0671%, pump +1.8590%; no binding K/pump ceiling on frozen trajectories; not a physiological capacity validation. | `analysis/47_dynamic_potassium_recycling_reconstruction/DYNAMIC_POTASSIUM_RECYCLING_REPORT.md` |
| 48 | AE4-null core beta invariance gives CCh=combined and zero contrast Jacobian in all 33 allowed shared parameter directions; genotype rest remains high-Cl/alkaline. | `analysis/48_joint_experimental_constraint_reconstruction/JOINT_EXPERIMENTAL_CONSTRAINT_REPORT.md` |

No literal production implementation error was established in the cited exclusions. Their inferential scope must remain visible: Task19 is numerically inconclusive beyond the valid 50% cases; Task31 local nonconvergence does not prove nonexistence; Task42 does not reject all possible kinetics of the source classes; Task47 does not validate biological channel capacity; Task48 invariance assumes a common source-faithful observation rule.

## Task36 basal NBC provenance and equation

The ledger's NBC warning is accurate. `analysis/36_minimal_nahco3_alkalinity/design.md:70` states explicitly that the implemented NBC is the stimulus-recruited increment and does not claim native basal NBC absence. At lines 74–79 it defines

`u_sec=clip((Ca-0.058)/(0.25-0.058),0,1)`.

Its reversible electrogenic rate is `J_B=G_B u_sec tanh(A_B/2)` with `A_B=log(Na_e HCO3_e²/(Na_i HCO3_i²))+V_b/V_T`. Source signature in `(Na,K,Cl,TIC,TA)` is `(1,0,0,2,2)`; cell charge source is −1 equivalent/cycle and outward basolateral conventional current is `F J_B`. The zero-rest gate was chosen to nest the accepted Task31 state exactly.

Active implementation confirms this at `src/modern_full_model/nbc_minimal.py`: `normalized_secretory_activation` (129 onward), `cycle = capacity*u*tanh(...)` (218), exact zero-activation membrane closure delegation (267 onward), and full-wrapper nesting return (607–645). Task37 reports bitwise-identical REST RHS and zero NBC flux/current.

This provides an exact gating assumption for audit. It does **not** select basal NBC as the answer. A rank-deficient fixed source basis or uncertain hidden resting coordinates must not be resolved by switching it on or fitting genotype multipliers.

## Other equation checks relevant to scope

The Palk core `src/modern_full_model/nkcc1_palk2010.py:20–39` has `X=Na_i K_i Cl_i²` and `J=alpha*m*(A1-A2 X)/(A3+A4 X)`; it has no bath-substrate arguments. Its frozen kinetics therefore do not exactly represent low-extracellular-Cl or bicarbonate-free isolated uptake experiments. This is documented and is not a Task49 licence to replace NKCC.

Task48's equation-level contradiction is explicitly documented in `JOINT_EXPERIMENTAL_CONSTRAINT_REPORT.md:43–79`: after AE4 deletion the core is `dx/dt=F_KO(x,u_CCh;theta)`, independent of beta. Its contrast sensitivity is zero for all 33 provenance-eligible shared parameter directions. The existing N1 regulator accepts calcium only and cannot break that invariance by parameter fitting. StageI's no-beta-NKCC restriction is therefore scientifically consequential, not merely organisational.

`history_audit.json` contains the findings, final historical counters, and SHA-256/size records for every referenced primary repository artefact at the pinned head.
