# Task 13C — Calcium-only WT secretion sweep

## Mission

Task 13B ended with ten retained native-source WT resting roots whose 600-second WT dynamics passed the numerical, conservation, solver, sustainment, regulatory, co-stimulation, and nearby-state gates, but every root failed the absolute one-SMG secretion-scale gate. The strict AE4-null stimulated-secretion holdout remained sealed.

Before introducing another transport mechanism or declaring a new experiment decision-critical, test the simplest outstanding possibility.

> Holding the entire Task 13B model fixed, does increasing only the stimulated intracellular calcium amplitude close the WT absolute secretion deficit?

Read `AGENTS.md` first and obey it exactly.

Do not edit `archive/`.

Do not merge to `main`.

This is deliberately a small, focused numerical experiment. Do not restart Task 13B, rebuild the model, or launch a broad multi-agent reconstruction.

---

## 1. Starting point that must remain frozen

Use the exact Task 13B branch state inherited into this branch.

The controlling Task 13B result is recorded in

- `analysis/13B_modern_full_model/final_answer.md`
- `analysis/13B_modern_full_model/generation_log.md`
- `analysis/13B_modern_full_model/dynamic_validation.md`
- `results/13B_modern_full_model/final_classification.json`
- `results/13B_modern_full_model/native_dynamic_contract_manifest.json`
- `results/13B_modern_full_model/native_dynamic_contract_group_gate.csv`
- `results/13B_modern_full_model/native_source_map.json`

Reproduce enough of the saved hashes and final gate records to establish that the ten final WT roots and their parameters are exactly the Task 13B ones before running the sweep.

The ten retained roots are the final `N_ABS_NKCC` scale-4 production roots spanning the five topology panels in the two retained routing/pump contexts. Do not re-solve, refit, rank, prune, or replace them.

### Frozen parameters

No model parameter may be optimized or changed except the stimulated calcium amplitude.

In particular do not change

- NKCC1 scale or kinetics;
- AE4 carrier amount, rate gauge, stoichiometry, routing, or regulation;
- NHE1 or AE2;
- Na/K pump capacity or membrane fraction;
- K-channel total conductance or membrane fraction;
- CaCC maximum conductance;
- paracellular conductances;
- hydraulic permeabilities;
- lumen dead volume or outflow rate;
- cell/gland geometry;
- bath composition;
- cAMP/PKA gains or time constants;
- resting roots or initial states;
- solver tolerances;
- one-SMG scale ceiling or WT 9--10 uL/min observation envelope.

Do not add a new mechanism in this task.

---

## 2. Why this test is required

The final Task 13B native dynamic contract used a stimulated calcium amplitude of `0.10 uM`. That value was a declared model sensitivity, not a direct measurement under the exact 2015 SMG CCh+IPR protocol.

Two later papers from the same modelling programme establish that calcium amplitude strongly controls mean acinar flow and that mean calcium is sufficient when mean flow is the output of interest.

### Vera-Sigüenza et al. 2019

Bull Math Biol 81:699--721, DOI `10.1007/s11538-018-0534-z`.

Relevant facts

- excised mouse parotid lobules were stimulated with `300 nM CCh`;
- resting cell-wide calcium was about `72 nM`;
- after stimulation the cell-wide calcium oscillated around a mean of about `200 nM`;
- the corresponding secretion model produced about `78 um^3/s = 0.078 pL/s` mean flow from about `22 um^3/s` at rest;
- replacing the oscillatory calcium trace with its spatio-temporal mean gave essentially the same mean flow;
- the paper therefore supports a mean-calcium input for mean secretion as a modelling approximation.

### Vera-Sigüenza et al. 2020

Bull Math Biol 82:38, DOI `10.1007/s11538-020-00712-3`.

Relevant facts

- a seven-cell parotid acinus at mean calcium around `0.18--0.20 uM` produced about `481 um^3/s` total flow;
- representative individual cells produced about `68.7--69.7 um^3/s`;
- the paper states that the earlier dual-agonist CCh+IPR experiments underlying the 2018/2015 lineage corresponded to intracellular calcium around `0.5 uM` and that raising the model calcium to `0.5 uM` reproduces the larger secretion increase;
- the same paper establishes that the total flow of an acinus is well approximated by the sum of independent representative-cell flows.

These are parotid results and must not be relabelled as direct matched-SMG calcium measurements. Their role here is to show that the Task 13B value `0.10 uM` is too narrow to be treated as the only physiologically relevant calcium amplitude and to define a source-informed sensitivity range.

---

## 3. Predeclared calcium grid

Keep resting calcium at `0.058 uM`.

For stimulation use exactly

```text
0.10
0.15
0.18
0.20
0.25
0.30
0.35
0.40
0.50
```

all in `uM`.

Do not choose extra grid points after seeing the results except for the one-dimensional threshold refinement described below.

For each calcium value report the existing model's calcium-gate activation

\[
p(C)=\frac{C^{1.46}}{0.26^{1.46}+C^{1.46}}.
\]

Do not modify the gate itself.

---

## 4. Protocol arms

### Primary scientific arm

Use the Task 13B `CCH_IPR` protocol

- `0.3 uM` CCh;
- `5 uM` IPR;
- `600 s` duration;
- beta input exactly as frozen in Task 13B;
- only `stimulated_calcium_uM` varies across the grid.

This arm determines the Task 13C scientific classification.

### Contextual control

Also run `CCH_ONLY` at the same calcium grid for the preferred regulatory member on the representative retained root, solely to display the calcium-flow relation in a beta-free arm and connect to the 2019/2020 CCh-only parotid context.

Do not use the parotid flow values as an SMG fitting target.

---

## 5. Computational design

The point of Task 13C is to avoid another 24-hour broad search.

### Stage A — calcium screen

For each of the ten retained WT roots, run the nine calcium amplitudes under

- the preferred Task 13B R1 regulatory member `R1_G125_P21_PROSE_TREFERENCE`;
- the primary `CCH_IPR` arm;
- the production Radau solver.

This is exactly `10 x 9 = 90` primary screening trajectories.

For every trajectory record

- root ID and full frozen-parameter hash;
- calcium amplitude and gate activation;
- flow at 60, 120, ..., 600 s;
- 0--600 s mean flow;
- cumulative flow;
- minimum and maximum flow after stimulation;
- sustainment ratio;
- shared WT observation-scale interval using the existing `9--10 uL/min` envelope;
- one-SMG geometry-scale pass/fail using the unchanged ceiling `3088.15386 uL min^-1 per (pL s^-1)`;
- implied effective cell count if a feasible scale exists;
- endpoint intracellular Na, K, Cl, pH, and volume;
- endpoint luminal ions/pH/volume;
- apical and basolateral potentials;
- all existing conservation/current/carbon/water residual gates;
- solver success and positivity.

Use the existing Task 13B gate implementation wherever possible rather than reimplementing the logic.

### Stage B — threshold refinement

If the one-SMG absolute WT gate changes from fail to pass between two adjacent predeclared calcium values, refine only calcium within that bracket.

Requirements

- do not change any other parameter;
- use deterministic bisection or a monotone one-dimensional search only if the relevant flow/gate statistic is empirically monotone over the bracket;
- target calcium precision `<= 0.005 uM`;
- if monotonicity fails, do not force a threshold; report the nonmonotone region.

Record a threshold separately for each retained root. Also report the minimum calcium at which any root passes and the calcium at which all ten roots pass, if such values exist.

### Stage C — confirmatory robustness

Do not rerun the full 20-member regulatory ensemble at all nine calcium values.

Confirm only

- every calcium value that first produces a Stage A one-SMG pass for at least one root;
- the immediately lower predeclared calcium value;
- any refined threshold bracket endpoints.

For those confirmatory calcium values, run across

- all ten retained roots;
- all frozen Task 13B regulatory members required by the native dynamic contract;
- both production solvers Radau and BDF.

Reuse the same numerical, conservation, sustainment, co-stimulation, and one-SMG gates from Task 13B.

A calcium value counts as a robust WT pass only if the predeclared Task 13B composite WT criteria pass without parameter changes.

If no Stage A calcium value passes, Stage C should crosscheck only `0.50 uM` with BDF and the regulatory ensemble to verify the negative endpoint. Do not launch additional mechanism searches.

---

## 6. Interpretation checks

The final analysis must answer these questions directly.

1. How strongly does per-cell secretion increase over `0.10--0.50 uM`?
2. Is the relation monotone for all ten roots?
3. At what calcium, if any, does each root first satisfy the one-SMG absolute WT scale?
4. Does a calcium in the `0.18--0.20 uM` CCh-only parotid context already close the SMG amplitude deficit?
5. Does only the high `0.35--0.50 uM` range close it?
6. Does `0.50 uM` still fail?
7. Does increasing calcium cause a different physiological gate to fail first, such as excessive ion depletion, nonphysical volume, membrane-potential pathology, loss of sustainment, or solver instability?
8. How do the final per-cell WT flows compare, as contextual benchmarks only, with the `~0.069--0.078 pL/s` parotid single-cell values from 2019/2020?
9. Is the Task 13B statement that a new per-cell secretion measurement is decision-critical still justified after this sweep?

Do not answer question 9 in advance. Let the sweep decide it.

---

## 7. Holdout firewall

The strict AE4-null stimulated secretion target remains sealed throughout Task 13C.

Do not

- open or parse `heldout_targets.csv` for the phenotype values;
- run `AE4_NULL`;
- run genotype comparison code;
- inspect the KO/WT secretion ratio;
- inspect early-versus-sustained knockout timing;
- use any knockout saliva value to choose calcium.

The calcium grid was frozen before this run from WT/parotid context and must remain independent of the held-out AE4 phenotype.

---

## 8. Implementation discipline

Prefer a new small runner such as

`src/modern_full_model/run_calcium_sweep.py`

that calls the existing frozen Task 13B model and gate code.

Do not duplicate the full model.

Add focused tests that prove

- every non-calcium frozen parameter is byte-identical or hash-identical across sweep members;
- resting roots are unchanged;
- only the stimulated calcium field changes;
- `0.10 uM` reproduces the Task 13B reference result within tolerance;
- the holdout code path is never invoked;
- calcium-gate activation is computed correctly;
- Stage B refinement changes only calcium;
- final artifacts are deterministic.

Before launching confirmatory simulations, make and push a checkpoint commit containing the runner, tests, and Stage A results so work cannot be lost to runtime/token limits.

---

## 9. Outputs

Create

### Analysis

- `analysis/13C_calcium_input_sweep/evidence_and_freeze.md`
- `analysis/13C_calcium_input_sweep/calcium_response.md`
- `analysis/13C_calcium_input_sweep/final_answer.md`

### Machine-readable results

- `results/13C_calcium_input_sweep/calcium_screen.csv`
- `results/13C_calcium_input_sweep/calcium_confirmatory.csv`
- `results/13C_calcium_input_sweep/calcium_threshold.json`
- `results/13C_calcium_input_sweep/frozen_manifest.json`
- `results/13C_calcium_input_sweep/artifact_hashes.csv`

`calcium_screen.csv` must have one row per root/calcium Stage A trajectory with explicit fields for the flow, gate, physiology, and numerical quantities listed above.

`calcium_threshold.json` must distinguish

- first provisional Stage A pass;
- refined calcium threshold if defined;
- robust Stage C pass threshold if defined;
- whether all ten roots agree qualitatively;
- whether the threshold lies inside a source-informed calcium context or only at the high edge.

---

## 10. Required figure/data summaries

Produce machine-readable data sufficient for these plots. PNG/PDF plots are optional, but the CSV data are mandatory.

1. per-cell mean flow versus calcium for all ten roots;
2. required whole-gland observation scale versus calcium with the one-SMG ceiling marked;
3. calcium-gate activation versus calcium;
4. endpoint Cl, K, Na, pH, and volume versus calcium;
5. minute-wise flow traces for the baseline `0.10 uM`, the first passing calcium if any, and `0.50 uM`.

---

## 11. Completion classifications

End Task 13C with exactly one of the following.

### Outcome 1

`WT ABSOLUTE FLOW GATE PASSES WITH CALCIUM-ONLY CHANGE`

Use when at least one source-informed non-edge calcium value gives a robust Task 13B composite WT pass with no other parameter changes.

### Outcome 2

`WT ABSOLUTE FLOW GATE PASSES ONLY AT HIGH CALCIUM EDGE — MATCHED-SMG CALCIUM REMAINS DECISION-CRITICAL`

Use when WT flow closes only near the upper edge of the `0.35--0.50 uM` sensitivity range and existing evidence cannot establish that such calcium applies to the exact SMG protocol.

### Outcome 3

`CALCIUM-ONLY CHANGE INSUFFICIENT TO CLOSE WT ABSOLUTE FLOW DEFICIT`

Use when even `0.50 uM` does not produce a robust WT absolute-scale pass.

### Outcome 4

`CALCIUM SWEEP EXPOSES PHYSIOLOGICAL OR NUMERICAL FAILURE BEFORE FLOW CLOSURE`

Use when increasing calcium raises flow but causes another predeclared physiological or numerical gate to fail before the absolute flow gate can be satisfied.

Do not reveal the AE4-null phenotype in any outcome.

---

## 12. Final repository actions

Run the focused Task 13C tests and the relevant inherited Task 13B regression tests.

Commit all Task 13C code, analyses, and results to the current branch.

Push the branch.

Do not merge to `main`.
