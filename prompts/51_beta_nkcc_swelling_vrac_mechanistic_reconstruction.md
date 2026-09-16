# Task 51 — beta-NKCC -> swelling -> VRAC mechanistic reconstruction

## Mission

There are **two** funded/scientific Codex shots remaining, including this Task 51 run.

Task 50 already provides a frozen phenotype-calibrated proof-of-sufficiency benchmark. Do not replace, rewrite or lose it. Task 51 asks one harder question:

> Can independently documented beta/cAMP salivary physiology generate the missing effective coupling mechanistically, without the Task 50 multiplier and without fitting the AE4 secretion phenotype?

Test exactly one predeclared chain:

`beta/IPR -> beta/cAMP activation of NKCC1 -> Na/K/Cl solute loading -> cell swelling -> VRAC-like apical anion conductance -> chloride/fluid secretion`

in parallel with the already represented

`beta/PKA -> AE4 activation`.

This is not a mechanism search.

## Repository and branch

Repository:

`esig626/ae4-salivary-transport-control`

Scientific parent / clean merged main:

`4bb2c89870fe3f4d4bb887cf9d40ba7de377a725`

Work only on:

`analysis/task-51-beta-nkcc-swelling-vrac-mechanistic-reconstruction`

Read `AGENTS.md` first, then execute from

`analysis/51_beta_nkcc_swelling_vrac_mechanistic_reconstruction/CODEX_START_HERE.md`.

The authoritative project memory is `docs/MANDATORY_RESEARCH_LEDGER.md`. Re-read it before **every scientific decision**, not only at startup.

Also binding:

- `docs/PHENOTYPE_TARGET_CONVENTION.md`
- `docs/TASK50_EFFECTIVE_COUPLING_BENCHMARK.md`
- `docs/TASK51_BETA_NKCC_SWELLING_VRAC_EVIDENCE.md`
- `analysis/51_beta_nkcc_swelling_vrac_mechanistic_reconstruction/NOVELTY_CHECK.md`
- `analysis/51_beta_nkcc_swelling_vrac_mechanistic_reconstruction/DECISION_LOG.md`

## Frozen Task 50 benchmark — preserve, do not inherit

The exact working Task 50 model is preserved at

`archive/task-50-working-effective-coupling-benchmark`

pinned to

`b16c30094b95f61a73d8f1977cd58e79c7bb50f6`.

Its fixed effective law is

`1 - 0.89488127156712 * a_Ca * beta * (1 - e_AE4)`.

Its frozen combined-stimulus benchmark includes:

- WT cumulative 0-600 s secretion `0.992524544081835 pL`;
- 5% AE4 cumulative deficit `23.163365263893244%`;
- AE4-null cumulative deficit `30.26115064458511%`;
- AE4-null endpoint flow deficit `27.038018389118477%`;
- substantial deficit across the saved 60-600 s window.

These are **TARGET-CALIBRATED CONSTRUCTION** outputs. They are the quantitative result Task 51 should try to explain mechanistically, not numbers Task 51 may fit.

### Absolute prohibition

Do not import, call, multiply by, reproduce algebraically, or otherwise hide the Task 50 AE4-expression-dependent conductance factor inside the Task 51 model.

A source-backed Task 51 model must generate any AE4 dependence through transport/state dynamics alone.

## Why Task 51 is new relative to Task 49

Task 49 tested only

`beta * positive swelling -> VRAC-like current`

on a parent core in which beta affected chemistry through AE4 only. In AE4-null IPR-only stimulation, exact genotype rest therefore remained chemically static: no beta-driven loading, no swelling, no VRAC activation, zero predicted uptake.

Task 51 does not change the downstream gate family. It adds the independently documented upstream fact Task 49 deliberately excluded:

`beta/cAMP -> increased NKCC1 activity`.

That new input can change Na/K/Cl loading and volume even when AE4 is absent and therefore can self-start the beta secretory branch.

Do not spend this run rediscovering the Task 49 zero-gate theorem.

## Source facts that justify the one mechanism

Use and independently verify the detailed evidence in `docs/TASK51_BETA_NKCC_SWELLING_VRAC_EVIDENCE.md`.

### Catalán et al. 2015 — mouse SMG

DOI `10.1073/pnas.1415739112`, PMID `25646474`, PMCID `PMC4343136`.

Source-backed:

- IPR/beta secretion persists after acinar TMEM16A deletion;
- it is also not abolished by CFTR or ClC-2 deletion;
- DCPIB/NPPB strongly suppress the IPR response;
- IPR causes approximately `12.5 +/- 0.2%` cell swelling by the paper's narrative value;
- IPR activates an outwardly rectifying chloride conductance with reversal near chloride equilibrium and a DCPIB-sensitive component;
- Task 49 source extraction records approximately `208.6 +/- 35.2 pA` IPR-induced whole-cell current and approximately `107.8 +/- 17.2 pA` DCPIB-blocked component at the reported voltage condition;
- molecular VRAC identity is not proven.

Retrieve and inspect the paper's supporting information before conductance inference. Reuse Task 49 digitisation/source files rather than re-digitising validated figures unless a documented defect is found.

### Beta/cAMP regulation of salivary NKCC1

Source anchors:

- Turner et al. 1992, DOI `10.1172/JCI115695`, PMID `1313447`: rat parotid, approximately 3-fold functional NKCC increase after brief IPR, `K1/2 = 21.5 nM`, cAMP/kinase dependent.
- Tanimura et al. 1995, DOI `10.1074/jbc.270.42.25252`, PMID `7559664`: IPR-regulated NKCC phosphorylation, half-maximal effect about 20 nM; strong regulated activity recruitment.
- Kurihara et al. 1999, DOI `10.1152/ajpcell.1999.277.6.C1184`, PMID `10600770`: beta stimulation increases NKCC activity/phosphorylation and bumetanide-binding sites.
- Kurihara et al. 2002, DOI `10.1152/ajpcell.00352.2001`, PMID `11880270`: cAMP reproduces regulatory phosphorylation and PKA inhibition blocks it.
- rat submandibular acini, PMID `9880083`: IPR increases NKCC-mediated NH4 influx approximately 2.5-fold; forskolin mimics the response and cAMP-dependent protein-kinase inhibition blocks it.

Do not copy a rat gain as if it were a measured mouse-SMG constant. Rat values establish sign, mechanism and plausible scale. Mouse-SMG Catalán swelling/current observations have priority for effective mouse parameter identification.

### Existing AE4 beta arm

Peña-Münzenmayer et al. 2021, DOI `10.1152/ajpgi.00145.2021`, PMID `34585968`, PMCID `PMC8887885`, establishes beta/PKA activation of AE4 and S173 dependence in mouse SMG. The model already includes this arm. Do not add it again.

## Task 51 model architecture

### A. NKCC core remains frozen

Retain:

- the source-fixed Palk/Benjamin concentration-response law;
- exact `1 Na : 1 K : 2 Cl` stoichiometry;
- thermodynamic reversal;
- existing Ca/CCh regulation at beta=0;
- resting NKCC multiplier exactly one.

Add only one beta/cAMP contribution to the existing NKCC **activity multiplier**.

### B. Minimal beta-NKCC algebra

Use the smallest no-interaction law.

Let `M_Ca(t)` be the existing inherited NKCC activity multiplier, which is exactly one at rest and contains the current Ca/CCh arm.

Let `M_beta` be the full beta-only multiplier relative to basal when `beta=1`.

The preferred predeclared combination is additive above baseline:

`M_NKCC(t) = M_Ca(t) + beta(t) * (M_beta - 1)`.

Properties:

- beta=0 -> exactly the inherited parent NKCC regulator;
- IPR-only with resting Ca -> `M_NKCC = M_beta`, so beta can start chemistry from rest;
- CCh+IPR -> Ca and beta increments add without an unmeasured synergy term;
- no change to concentration dependence or reversal.

If implementation details of the inherited regulator require an algebraically equivalent form, prove the equivalence and record it before code changes. Do not add a Ca×beta interaction, Hill coefficient, delay, phosphorylation ODE or saturation parameter unless a primary source makes that unavoidable. If so, stop and report before expanding the model.

### C. VRAC-like current

Reuse the Task 49 conservation-tested structure:

`G_V = g_V * beta * max(V_i / V_rest,genotype - 1, 0)`

with the same current/sign/source conventions and no direct AE4-expression factor.

VRAC must contribute to intracellular/luminal chloride balances and electrical closure as a real current, not as a water/secretion multiplier.

No alternative VRAC gate family is authorised in Task 51.

## Parameter identification — no AE4 phenotype allowed

There are at most two new effective quantities:

1. beta-NKCC strength `M_beta`;
2. VRAC-like conductance scale `g_V`.

No other fitted scientific parameter is permitted.

### Primary identification route

Before any AE4 phenotype comparison, determine whether the Catalán supporting information permits a defensible current-to-conductance observation map. Explicitly derive extracellular/intracellular chloride conditions, reversal potential, voltage convention and whether the DCPIB-sensitive current can be treated as the relevant conductance component.

If it does, identify/freeze `g_V` from the mouse-SMG current data and identify/check `M_beta` using the independent salivary beta-NKCC measurements plus the mouse-SMG swelling response. Use analytic equations, scalar roots or at most one two-parameter deterministic continuous solve. No grid, random search or model family.

### Predeclared fallback if direct current mapping is not identifiable

Do **not** stop the whole task merely because the current-to-conductance map remains underidentified.

Use this fallback, fixed before AE4 reveal:

1. Fix the primary beta-NKCC strength to the independently measured **2.5-fold rat submandibular** IPR response as the central cross-species effective value, clearly labelled an ASSUMPTION/SENSITIVITY because it is not mouse.
2. Use the approximately **3-fold rat parotid** functional response as a predeclared source sensitivity point, not as a value selected after seeing AE4 results.
3. With the primary `M_beta` fixed, identify the single `g_V` scalar from the independent mouse-SMG IPR swelling magnitude (approximately 12.5%) using one deterministic scalar solve.
4. Use the Catalán current/blocker data as a consistency/validation check rather than forcing an unsupported current conversion.

The approximately sixfold context in the older parotid phosphorylation literature is not an automatic Task 51 parameter and must not become a post hoc rescue value.

### Identification failure

If neither primary mapping nor the predeclared swelling fallback yields a physically admissible, independently fixed model, publish that non-identifiability/failure. Do not use the AE4 35% target or Task 50 lambda to close the system.

## Held-out data boundary

Before the immutable 51D prediction checkpoint, do not use to choose or change a parameter:

- JBC 2015 AE4-null 10-min secretion deficit `35 +/- 4.7%`;
- AE4-null CCh+IPR uptake `0.90 +/- 0.09 x10^-3 s^-1`;
- AE4-null CCh-only uptake `2.30 +/- 0.10 x10^-3 s^-1` except as a later nesting/validation comparison;
- AE4-null IPR-only uptake `0.20 +/- 0.03 x10^-3 s^-1`;
- Task 50 30.2612% null or 23.1634% 5%-AE4 deficits.

The values may be present in the ledger because repository memory is transparent. Governance means they are not objectives/calibration residuals before the freeze.

## Two simulation tracks — do not mix them

### Track 1: Task 50 apples-to-apples benchmark

Use the same Task 40/41 common rest, bath, onset and CCh+IPR protocol as the frozen Task 50 benchmark. Run WT, 5% AE4 and AE4-null with the Task 51 mechanism and **without Task 50 multiplier**.

This track answers whether the literature-backed mechanism can replace the effective Task 50 coupling on the same chassis.

### Track 2: established-genotype / Task 48 validation

Where executable without changing the Task 48 bath/observation rules, use genotype-specific rests for WT, AE4 KO and AE2 KO/control to test:

- IPR-only self-start and chloride response;
- CCh-only versus CCh+IPR distinction;
- AE2 specificity;
- resting-state limitations.

Do not substitute Track 1 output for Track 2 or vice versa. The chronic resting Cl/pH failure remains a separate unresolved limitation unless Task 51 changes it naturally.

## Required diagnostics

For every production trajectory record at least:

- total NKCC activity multiplier and its Ca versus beta contributions;
- NKCC cycle and chloride loading flux;
- AE4 and AE2 chloride loading;
- VRAC conductance, current and chloride flux;
- TMEM16A/CaCC current/flux separately;
- intracellular Na/K/Cl, TIC/HCO3, pH and volume;
- apical/basolateral/transepithelial voltages;
- pump and K fluxes;
- luminal osmolarity/volume and water outflow;
- cumulative secretion;
- conservation/current-closure residuals.

Explicitly report whether beta-NKCC actually produces positive swelling from AE4-null IPR rest and therefore removes the exact Task 49 self-start failure.

## Checkpoints and publication

Use one orchestrator and up to five advisory/read-only workers. Suggested workers:

1. repository/no-repeat/provenance audit;
2. Catalán main+supplement extraction and conductance map;
3. beta-NKCC literature reconstruction and cross-species qualification;
4. independent conservation/electrical implementation review;
5. validation/freeze/reveal audit.

Only the orchestrator may modify the canonical branch, accept scientific conclusions, run authoritative inference/trajectories, commit or push.

### 51A — source and architecture freeze

Before production edits:

- verify every source claim;
- retrieve/inspect Catalán SI;
- reuse Task 49 digitisation unless a documented defect exists;
- prove the beta-NKCC algebra and beta=0 nesting;
- prove Task 50 multiplier is absent from the planned model;
- predeclare parameter-identification route and fallback;
- update ledger if any source fact changes.

Commit, push, remotely verify.

### 51B — implementation and unit/conservation verification

Implement only the beta contribution to the inherited NKCC regulator plus production reuse/promotion of the Task 49 VRAC-like current.

Required tests:

- REST exact parent nesting;
- beta=0/CCh-only exact parent nesting;
- IPR-only beta increases NKCC without moving its thermodynamic reversal;
- exact 1Na:1K:2Cl NKCC bookkeeping;
- zero VRAC when beta=0 or swelling gate=0;
- correct VRAC reversal/current/Cl source sign;
- electrical/current closure with VRAC;
- no direct VRAC dependence on AE4 expression;
- no import/use of Task 50 effective-coupling factor;
- disabling VRAC leaves NKCC/AE4 laws unchanged;
- disabling beta-NKCC leaves the Task 49 self-start limitation reproducible as a regression fixture.

Synthetic states are allowed for software tests only and must not become scientific calibration points.

Commit, push, remotely verify.

### 51C — independent parameter identification and beta-pathway validation

Use only independent source data.

Execute the primary conductance mapping if identifiable; otherwise execute the predeclared 2.5-fold beta-NKCC + mouse-SMG swelling fallback. Freeze all parameter values, source uncertainty/sensitivity points and observation mappings before any AE4 phenotype reveal.

Validate against independent beta-pathway facts where representable:

- IPR causes positive swelling in mouse SMG direction/magnitude;
- a TMEM16A-off idealisation retains an IPR-driven secretory/anionic branch;
- VRAC blockade suppresses that branch in the model idealisation;
- predicted IPR-induced current direction/magnitude is compatible with source data to the extent the observation map permits;
- beta-NKCC source sensitivity values remain within the predeclared salivary evidence envelope.

Do not use AE4 phenotype residuals.

Commit, push, remotely verify the frozen parameter file.

### 51D — immutable pre-AE4 prediction checkpoint

With the independently fixed model, generate and publish predictions before phenotype comparison.

Track 1:

- WT CCh+IPR;
- AE4 5% CCh+IPR;
- AE4-null CCh+IPR.

Track 2 where executable:

- WT CCh, IPR, CCh+IPR;
- AE4 KO CCh, IPR, CCh+IPR;
- AE2 control/KO corresponding protocols.

Save all required diagnostics and source/parameter hashes.

Commit, push, remotely verify. This remote SHA is the immutable reveal boundary.

### 51E — held-out AE4 reveal and Task 50 comparison

Only after 51D is verified, compare frozen predictions against:

- experimental AE4-null secretion magnitude;
- CCh-only vs CCh+IPR uptake contrast;
- positive IPR-only AE4-KO response;
- AE2 specificity;
- frozen Task 50 benchmark.

No retuning.

Report explicitly:

1. Did beta-NKCC remove Task 49's IPR self-start failure?
2. Did the mechanism produce a substantial persistent AE4-null secretion deficit on the Task 50 benchmark chassis?
3. What is the 5%-AE4 deficit?
4. How close is the independently predicted magnitude to Task 50 and experiment, without calling closeness a fit?
5. Does NKCC still compensate excessively, or does the extra beta-secretory demand expose an AE4 dependence?
6. What physical/state constraints fail, if any?
7. Does AE2 remain approximately secretion-neutral?
8. What chronic rest mismatch remains?

Commit, push, remotely verify.

### 51F — final mechanistic decision

Create:

`analysis/51_beta_nkcc_swelling_vrac_mechanistic_reconstruction/BETA_NKCC_VRAC_RECONSTRUCTION_REPORT.md`

Classify the outcome precisely:

- **source-backed prospective mechanistic success** if independently constrained beta-NKCC/swelling/VRAC physiology produces the substantial AE4 phenotype without Task 50 multiplier;
- **partial mechanistic success** if the beta architecture/self-start is fixed but phenotype magnitude remains insufficient or other key constraints fail;
- **negative/exclusion result** if the source-backed chain cannot produce the required phenotype on the conservation-explicit chassis;
- **unresolved/non-identifiable** if independent source data cannot fix the required quantities.

If Task 51 fails, do not invent another mechanism. Preserve Task 50 as the quantitative proof that an additional unidentified beta-conditioned network interaction is sufficient/required by this model class.

Update `docs/MANDATORY_RESEARCH_LEDGER.md` with the controlling result in the same checkpoint.

Commit, push and remotely verify 51F. Stop. Do not start Task 52.

## Absolute prohibitions

- no Task 50 multiplier in Task 51 equations;
- no AE4-phenotype calibration before 51D;
- no post-reveal retuning;
- no broad parameter or mechanism search;
- no grid search over beta gain/VRAC conductance;
- no alternative NKCC core/stoichiometry/cap;
- no alternative VRAC gate family;
- no AE4-dependent TMEM16A/VRAC multiplier;
- no NHE/NBC/pump/K/calcium/AE4 stoichiometry/routing repair;
- no exact-time fitting;
- no new signalling ODE;
- no claim that DCPIB proves LRRC8 molecular identity;
- no claim that rat NKCC gains are exact mouse-SMG constants;
- no automatic Task 52 work.
