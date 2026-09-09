# Adversarial audit of the state-resolved AE4 reconstruction

**Status: checklist frozen before the Stage A, Stage B, or held-out results were
available.** This audit is deliberately independent of model selection. A
green numerical run is not a successful biological mechanism, and a failed
phenotype reconstruction is not automatically evidence for one particular
missing protein.

## Audit boundary and decision rule

The audit will try to falsify both possible end states:

1. If a model is called successful, it must remain successful after leakage,
   branch, solver, scaling, uncertainty, nested-model, thermodynamic, and
   source-provenance challenges.
2. If the result is called unresolved, the proposed bottleneck must follow
   from exact source geometry or independently reproduced numerical evidence.
   It must not be a euphemism for failed optimization, an undocumented state
   bound, an incomplete root search, or lack of implementation effort.

The evidence priority and claim classes in `AGENTS.md` and
`evidence_freeze.md` control this audit. In particular, the 2015 AE4-null
secretion magnitude and time course are never calibration inputs. The Stage B
ionic phenotypes may localize a missing balance only after all admissible
state-resolved families fail Stage A on the fixed chassis.

No audit item may be marked `PASS` solely because a unit test asserts the
author's intended behavior. The underlying equations, saved inputs, and
machine-readable results must agree.

## Frozen falsification checklist

### A. Holdout firewall and model ledgers

- [x] Hash or otherwise freeze the calibration and held-out ledgers before
  reading predictions; verify that AE4-null secretion fields are absent from
  every optimization objective, bound, model-selection score, stopping rule,
  protocol tuning rule, and module-selection rule.
- [x] Inspect code, not only prose, for numerical constants derived from
  `0.65`, `0.35`, `35.851%`, digitized null trajectories, late-null/WT ratios,
  or a three-minute divergence. A renamed constant still counts as leakage.
- [x] Verify every `model_registry.csv` and `model_scores.csv` row records
  whether knockout secretion was exposed before freeze. Any exposed row must
  be excluded from predictive rankings, even if exposure was accidental.
- [x] Verify Stage A did not use AE4-null Cl, pH, uptake, or any claimed volume
  target. Verify Stage B uses only the permitted ionic localization evidence;
  the primary 2015 paper contains no direct AE4-null acinar-volume target.
- [x] Check that AE2-null evidence and AE4-line WT data are not pooled across
  the distinct mouse backgrounds, and that `no detected difference` is not
  encoded as exact equality or an invented equivalence margin.
- [x] Check whether any post-reveal parameter-domain narrowing, branch choice,
  neutral-value convention, or regularization weight was selected because it
  improved the held-out secretion result. This is leakage even if the target
  is not present in the objective function.

### B. Transporter family completeness and nested equivalence

- [x] Confirm that at least one implemented family is genuinely state
  resolved: explicit carrier-state conservation, inward/outward conformations,
  reversible transitions, and separate Na/K coordination. A steady scalar
  rate wrapped in unused state variables does not qualify.
- [x] Derive the steady-state flux of every survivor independently. Test
  whether it is observationally identical to a Task 12 scalar law over the
  calibrated protocols. If so, state structure has not earned its parameters.
- [x] Confirm SR1 handles Na and K as separate species-specific parameter
  sets, not a fitted `[Na]+[K]` pool. Confirm SR2 truly shares carrier states
  rather than summing two independent scalar currents.
- [x] Challenge SR3's fixed-species `1:1:1:1` interpretation against direct
  2016 evidence that both Na and K can support transport separately. A
  thermodynamically feasible proposed cycle is not thereby an adequate sole
  mechanism.
- [x] Challenge SR4 with explicit carbonate mass balance, charge, acid-base
  equilibrium, and conversion kinetics. Carbonate may not be represented by a
  relabelled bicarbonate variable or by an unconstrained hidden reservoir.
- [x] Confirm SR5 retains direct cation transport and does not reinterpret
  measured Na/K signals as allosteric gating alone.
- [x] Confirm SR6 does not infer direct S173 phosphorylation, a physiological
  onset time, or a universal turnover multiplier from the 2021 experiments.
  Test the nested common-gate limit before crediting state-specific PKA action.
- [x] Penalize microscopic rates that transporter data do not identify. A
  model cannot win merely by fitting one effective EC50/Hill curve with many
  exchangeable rates.

### C. Stoichiometry, thermodynamics, and mutant logic

- [x] Recompute every source vector and verify zero net transported charge per
  completed cycle. Do not count charge cancellation caused by an omitted
  species or by the whole-cell membrane solver.
- [x] Recompute every cycle affinity from chemical potentials and verify local
  detailed balance around all independent graph cycles. Forward/reverse rate
  ratios must agree with the declared standard-state convention.
- [x] Verify zero flux at equilibrium and the correct sign on both sides of
  reversal over Na and K assay baths and physiological salivary ranges.
- [x] Test transporter conservation, state positivity, and probability sum
  under rest, stimulation, reversal, mutation, PKA, and knockout limits.
- [x] Verify genotype zero annihilates the complete AE4 contribution,
  including bound-state transients, regulation, and carbonate chemistry. If
  transporter states store ions, document whether deletion also removes those
  stored amounts and whether this matters only instantaneously.
- [x] Reproduce Na `EC50=49 mM`, Hill `2.0`, and K `EC50=62 mM`, Hill `1.8`
  only as effective assay constraints. Because fit uncertainty/covariance was
  not reported, challenge conclusions across a declared conservative range
  rather than inventing precise confidence bounds.
- [x] Verify the 2025 T756A-T448I qualitative discriminator: Na-dependent
  activity approaches nontransfected background while substantial K-dependent
  activity remains. Do not convert approximate single-mutant percentages into
  precise microscopic free energies.
- [x] Retain the HCO3 versus CO3 alternative as unmeasured unless a primary
  experiment discriminates it. A carbonate variant may be rejected for model
  inadequacy or missing chemistry, but not claimed biologically absent.

### D. Fixed-chassis Stage A

- [x] Reproduce all reported WT resting roots with deterministic multistart
  over the declared physical domain. Record remote/disconnected or
  bound-touching roots; do not silently select the one with favorable
  secretion.
- [x] Verify raw, scaled, and species-level residuals, positivity,
  electroneutrality, volume balance, both membrane-current closures, and all
  declared thermodynamic signs at each accepted root.
- [x] Confirm all non-AE4 chassis equations and parameters are identical
  across the state-resolved candidates. Candidate-specific WT calibration may
  vary only through predeclared allowed combinations.
- [x] Independently verify the key structural fact: with AE4 activity exactly
  zero, all candidates embedded in an unchanged chassis have the same
  AE4-null vector field. If true, state structure alone cannot repair a wrong
  null equilibrium; if false, identify the hidden candidate-dependent term.
- [x] Check AE2-null resting Cl, pH, uptake, and secretion separately. A near-
  unity dynamic flow ratio does not validate the chronic ionic phenotype.
- [x] Challenge all numerical time labels. Historical code time and absolute
  flow units remain uncertified unless Task 13 supplies a new independent
  dimensional ledger and WT protocol calibration.

### E. Stage B missing-balance localization

- [x] Verify Stage A failure is established across every transporter-level
  survivor before releasing the null ionic targets.
- [x] Recompute the AE4-null residual at the actual localization state. If Na,
  K, volume, membrane potential, or luminal variables are unmeasured, expose
  them as nuisance coordinates or an affine family; do not plug in convenient
  WT or historical values and call the correction unique.
- [x] Distinguish a vector-field correction at an observed coordinate from a
  displacement between equilibria. Convert between them only through a
  documented local Jacobian or re-solved nonlinear equilibrium.
- [x] Verify units and scaling before comparing Na, K, Cl, HCO3/H, current,
  and water residual components. Cosines in arbitrarily normalized
  coordinates are not physical evidence.
- [x] Recompute exact signature rank/span and sign-constrained one- and
  two-module projections. Report null spaces and nonunique decompositions,
  not only the smallest optimizer support.
- [x] Stress projections under uncertainty in the unobserved state and local
  Jacobian. A module is localized only if the support remains necessary over
  that admissible set.
- [x] Ensure current and water coordinates are not double-counted if they are
  algebraic consequences of ion fluxes in the chosen reduction.
- [x] Check that unchanged assay-level NKCC1 and NHE signals are treated as
  restrictions, not exact whole-cell flux equalities, and that arbitrary
  genotype-specific upregulation is prohibited.

### F. Minimal chassis reconstruction

- [x] Require an independent source for every added topology. The 2018
  Almassy abstract supports apical pump presence and an apical Ca-activated K
  current qualitatively; it does not by itself provide a measured apical
  pump/K fraction, whole-cell conductance, or genotype-specific recruitment.
- [x] Audit the membrane orientation signs of apical and basolateral Na/K
  ATPase and K conductance. Splitting one total source into two labels without
  distinct luminal/electrical consequences is not a new mechanism.
- [x] Remove each added parameter/module in turn. A claimed minimal extension
  must repair the localized residual and preserve WT/AE2 constraints; nuisance
  additions that do not alter the decisive balance must be deleted.
- [x] Fix unsupported localization fractions, turnover scales, or regulatory
  gains to neutral/source-boundary values. If the conclusion changes, report
  it as assumption-dependent rather than source-supported.
- [x] Test cation-current architecture alone, acid-base architecture alone,
  and the smallest coupled pair before accepting a larger module set.
- [x] Verify that HCO3/CO2 buffering conserves carbon and charge under the
  declared open-bath exchanges, and that pH is not repaired by an implicit
  infinite buffer or sign error.
- [x] Verify luminal outflow removes all relevant solutes consistently and
  that water/osmotic closure is not a freely scaled mapping from chloride
  current to the desired secretion ratio.

### G. Dynamics and held-out secretion

- [x] Freeze all parameters, branch rules, and input waveforms before opening
  the AE4-null flow columns. Preserve hashes or timestamps showing the order.
- [x] Reproduce WT stimulation first. A normalized KO/WT ratio cannot rescue a
  WT trajectory with wrong scale, sign, or protocol.
- [x] Test at least two stiff solvers or independent implementations, tighter
  and looser tolerances, smaller output steps, and perturbations of initial
  states. Report failures and positivity interventions.
- [x] Evaluate the 10-minute cumulative ratio, endpoint/late ratio, first
  2-3-minute comparability, and sustained late divergence separately. Matching
  one scalar integral is not success.
- [x] Do not fit a free delay, PKA time constant, calcium waveform, water
  coefficient, or pump/K recruitment law to the null curve. A delayed
  separation is predictive only if the timing is constrained independently.
- [x] Stress reasonable uncertainty in Na/K effective affinities and Hill
  slopes, PKA fold, cation localization fractions, and WT parameters. Report
  whether qualitative conclusions or classification change.
- [x] Check that digitization uncertainty is not treated as biological SEM,
  and that the exact 2015 result is `35 +/- 4.7%` less total 10-minute saliva;
  `30%` is rounded context and `24%` is a historical model output.

### H. Source and historical-claim discipline

- [x] Trace every substantive numerical bound and topology claim to a primary
  source, exact prior-project result, or explicitly labelled new assumption.
- [x] Separate heterologous CHO/HEK assay constraints from native acinar and
  whole-gland measurements; temperature and bath context travel with the
  number.
- [x] Do not promote the Catalan 2025 MD state order or proposed cycle drawing
  from hypothesis to measurement. Mutagenesis constrains functional cation
  nonidentity, not a unique transition sequence.
- [x] Do not claim the 2018 historical article was wrong or that AE4 itself is
  insufficient. Allowed exclusions are tied to an explicit model family,
  tested domain, chassis, protocol, and evidence gate.
- [x] Do not treat failure to access the full Almassy experimental Methods as
  a license to import fractions from secondary prose. Qualitative topology and
  quantitative calibration evidence must remain distinct.

### I. Irreducibility of an unresolved conclusion

- [x] Verify the alleged missing experiment corresponds to a parameter or
  flux direction that changes the model-class decision, rather than merely
  reducing uncertainty.
- [x] Demonstrate at least two source-consistent values or module
  parameterizations compatible with all released non-secretion evidence but
  giving substantively different held-out secretion predictions. Otherwise
  the missing quantity has not been shown decision-relevant.
- [x] Show that the ambiguity persists after multistart roots, independent
  solvers, and exact balance reductions. An optimizer failure is not an
  irreducible biological ambiguity.
- [x] State one executable experiment with genotype, preparation, stimulus,
  measured species/flux, timing, and discrimination criterion. “Measure pump
  activity” or “more Na/K data” is not specific enough.
- [x] Verify the proposed experiment observes the missing direction rather
  than a correlated downstream secretion phenotype already held out.

## Planned audit outputs

The completed audit will include:

- a table of every challenge, route, tolerance/domain, result, and
  disposition;
- a leakage crosswalk over all machine-readable ledgers;
- an independent calculation of source ranks/projections and the accepted or
  unresolved model's decisive numerical quantities;
- explicit discrepancies sent to the lead agent before synthesis; and
- a final `PASS`, `FAIL`, or `CONDITIONAL` disposition scoped to the exact
  Task 13 classification rather than to AE4 biology.

---

## Completed audit after Rounds 1--9

### Auditor disposition

**CONDITIONAL PASS** for the exact classification

`PHENOTYPE UNRESOLVED — SPECIFIC MISSING EXPERIMENT OR FLUX IDENTIFIED`

This is not a pass for an AE4 mechanism, a reconstructed secretion phenotype,
or a calibrated modern chassis. It is a pass for the narrower substantive
conclusion that the present evidence and implementations cannot identify the
missing whole-cell coupling, and that the remaining experimental bottleneck
has been localized precisely enough to act on. The conclusion survives the
leakage, alternate-root, optimizer, solver, nested-model, stoichiometric,
carbonate, PKA, topology-provenance, and historical-claim challenges below.

Three failures remain scientifically controlling:

1. no transporter-surviving row supplies a valid WT capacity or physiology
   gate on the fixed chassis;
2. the smallest independently supported pump/K topology extension also
   supplies no valid WT row and leaves the knockout Cl/pH mismatch large; and
3. the released knockout Cl/pH pair leaves a five-dimensional target-state
   family and cannot identify either membrane-split cation support or the
   conserved-carbon/NHE/buffer fluxes.

The complete 36-row challenge ledger is
`results/13_state_resolved_ae4/adversarial_checks.csv`.

### Most consequential falsification results

| Challenge | Attempted falsification | Result | Scientific disposition |
| --- | --- | --- | --- |
| E16-07 observable | Recompute the shared-network assay output as the selected branch current | That mapping is wrong; the measurement is total alkalinization. The code and regression test now use net `HCO3 + 2 CO3` source | Original implementation failure corrected before final inference |
| Carbonate artifact | Re-evaluate SR4 after correcting the observation map | SR4A gives finite EC50 values near `49/62 mM`, but Hill values near `1.11/1.13`; it still fails. SR4B still cannot express the direct-Na/mutant constraints | Carbonate remains a biological alternative; it is not embeddable in the bicarbonate-only chassis |
| Post-reveal chronology | Check whether the corrected transporter table silently replaced the revealed bytes | Original table hash `785d5c...a51f` is preserved; amended/current hash `e4c269...7e40`; explicit amendment records no knockout ionic or secretion use; survivors remain `5 -> 5` | Pass after auditable amendment |
| Sequential optimizer | Jointly fit SR5 gauge `0.1` rather than accepting the sequential fit | Exact dose summaries and permitted Hills are reachable, but the salivary chloride source remains `-0.0111649851` per capacity | Optimizer order does not rescue the required WT source sign |
| Wider transporter search | Differential evolution over barrier, two energy lumps, and K attempt scale | No gate-feasible positive chloride source in the declared finite domain; best is `-1.75670e-5` at an energy bound | Strong finite-domain exclusion, not a theorem over all microscopic rates |
| QSS state resolution | Ask whether explicit states add a whole-cell state or delayed memory | Whole-cell dimension remains seven; occupancies are algebraic diagnostics and the source is exactly representable by static branch functions `J_Na(x),J_K(x)` | State competition is real, but no carrier relaxation mode was tested |
| Fixed-chassis knockout field | Compare every candidate after genotype-zero deletion | The vector fields are exactly identical; bounded common root is `Cl=48.782215 mM`, `pH=7.378474` | More QSS AE4-only parameterization cannot move this knockout field |
| Global root claim | Search outside the declared physiological domain | Exact asymptotic roots occur around height `1.25e14` to `7.69e14 um` | Bounded-domain result survives; global uniqueness claim is prohibited |
| Missing-vector uniqueness | Release only knockout Cl/pH and recompute the admissible state/source family | Five state coordinates remain free, and the required-source family has numerical rank five at the mean and all one-SEM corners | There is no unique missing source vector at this evidence level |
| Local two-parameter repair | Inspect exact local pairs and their required finite changes | Exact two-output fits require log changes about `1.7--18.3`, change with the assumed target slice, and are outside a defensible local repair interpretation | Rank artifacts, not reconstructed chassis mechanisms |
| Split-cation topology | Re-solve eight predeclared pump/K topology sensitivities | Best row changes normalized Cl/pH distance only `25.6016 -> 24.8677` (`2.87%`) | M1 is a genuine omitted topology but fails as an ionic repair |
| Round-5 root robustness | Repeat seven saved M0/M1 roots with 12 starts and an alternate bounded route | One bounded root per tested topology; saved states reproduced; maximum raw residual `1.56e-12` | Finite-domain root result reproduced; remote branches are not excluded |
| Round-5 WT optimizer | Use 25 state/capacity starts on four coupled/cross topologies | Zero valid rows; height approaches `200 um`, raw closure `2.46e-5--4.79e-5`, pH `+8.72--+8.76` SEM | Negative gate is robust in this challenge, although the saved five-start minimum is not unique |
| Dynamic result | Try Radau and BDF on the only licensed knockout IVP | Numerical traces agree, but no WT-valid denominator, minute map, absolute-flow scale, or delay state exists | Numerical diagnostic only; no held-out secretion claim is licensed |

### Holdout and chronology crosswalk

The evidence and question files retained their frozen SHA-256 hashes:

- evidence freeze: `48022d5f69c4459edf6c9cceb914dd5ce092ffe73ffaffdbf9ceedc578005631`;
- question tree: `34343e076d46b00ce6eadc7a419688747b13917afce576611a2149c79aba019e`.

The evidence API fails closed if either file changes. `LOC-B` access is denied
until frozen Stage-A failure, and `STRICT` secretion IDs are denied to
calibration. The machine-ledger audit found:

| Artifact | Rows audited | Exposure result |
| --- | ---: | --- |
| `model_registry.csv` | 11 | every `stage_a_holdout_read_before_freeze=False` |
| amended transporter summary | 33 | every `heldout_target_read=False` |
| frozen survivor parameters | 5 | every `heldout_target_read=False` |
| WT capacity attempts / WT roots / common null root | 5 / 5 / 1 | every exposure field false |
| model scores / rejected models | 33 / 33 | no strict secretion used |
| CTMC frozen conditions | 17 | every `heldout_ko_data_read=False` |
| Round-5 chassis ledger | 33 | every `secretion_target_read=False` |
| dynamic status ledger | 21 | every strict-target fit and prediction-join field false |
| dynamic diagnostic trajectory | 24 | every `wt_denominator_valid=false` and `heldout_prediction=false` |
| held-out target ledger | 4 | every `read_before_model_freeze=False` |

The final dynamic table contains the strict secretion observations only as
reference rows after the numerical diagnostics. All ten candidate rows are
`blocked_before_dynamics`; endpoint and integral ratio cells are blank. No
target value is passed to the target-free validation API, which raises before
forming a ratio when the WT-capacity, physiology, or independent-chassis gates
are false.

### Transporter-level audit

The implementation contains genuinely state-resolved reversible graphs with
carrier conservation, separate Na/K coordination, and local detailed balance.
SR2 shares one carrier pool rather than adding independent scalar currents;
SR5 contains explicit catalytic-bound detours; SR6 contains explicit
unphosphorylated/phosphorylated layers and declared transition placements.
All registered branches are electroneutral and satisfy the declared affinity,
reversal, probability-conservation, and detailed-balance tests.

That structural statement must be separated from identifiability. The pure-
cation dose summaries identify effective loaded-state/attempt combinations,
not a microscopic binding order. The reported `Rmin` values are fixed
observation-map nuisances because no zero-cation point was available, and the
common assay scale is not an independently measured transporter capacity.
The broad EC50/Hill gates are new compatibility tolerances, not confidence
intervals. The T448I-T756A Na-collapse/K-retention result is imposed through a
declared mutant branch limit, so passing it demonstrates expressivity rather
than prediction. The SR6 survivors require large effective regulatory-layer
factors and have no independently measured phosphorylation or activation time.
No AIC-like or microscopic-mechanism winner is justified.

At whole-cell resolution the stationary carrier distribution is eliminated at
every RHS call. Therefore any implemented survivor is exactly a pair of
instantaneous scalar maps over whole-cell state. It can express finite-pool
competition and Na/K slip that a single pooled scalar cannot, but it adds no
carrier memory. A future dynamic-occupancy model remains untested and would
need an independent transporter relaxation time and bound-ion accounting
before it could explain the gland-level delay.

### Stage-A branch and solver audit

All five conditional transporter survivors approach the declared
`height=200 um` boundary during WT-Cl capacity calibration and retain raw
closure around `3.05e-5`, above the `1e-8` acceptance threshold. Frozen-
capacity exact roots exist in the declared physiological domain but miss WT pH
by about `+6.69` reported SEM. Failed calibration coordinates are never
promoted as parameter estimates.

After AE4 deletion the exact genotype-zero short circuit removes every
candidate contribution, including PKA and carbonate branches. The common
bounded root is reproduced independently. The unbounded positive search also
finds an ill-conditioned asymptotic ionic branch at enormous height; different
starts drift along that branch. This means the defensible statement is
candidate independence plus failure of all roots found in the declared
physiological domain, not global uniqueness over all positive states.

### Missing-balance and minimal-chassis audit

At the released means, electroneutrality fixes bicarbonate only after choosing
`(Na_l,K_l,height,Na_i,K_i)`. Those five unobserved values parameterize a
five-dimensional target manifold. The required source `-F_0(x)` rotates with
numerical local rank five across that manifold. Rankings of one- and two-module
projections change between the model-anchored and WT-reference slices. No
displayed one- or two-module cell signature is exact on both slices.

Apical and basolateral pump cycles have the same intracellular chemical
signature, as do apical and basolateral K efflux. Their split adds exact rank
two only after lumen and membrane-current coordinates are retained and rank
zero in the intracellular projection used by the released Cl/pH evidence.
This proves why the old projection cannot identify membrane fractions; it does
not prove the split is inert or sufficient.

The Round-5 M1 implementation passes the important structural tests: at zero
apical fractions it is algebraically nested in the fixed chassis; total pump
capacity and K conductance are partitioned rather than increased; the apical
pump uses local lumen K; lumen Na/K amount terms and both current closures are
explicit. The primary literature supports the topology qualitatively, but not
quantitative pump/K fractions in mouse submandibular acini. The `0.2`, `0.4`,
and `0.398` values are consequently sensitivities and may not be optimized
against secretion.

All eight exact knockout sensitivities remain more than 24 normalized Cl/pH
units from the released means. All five saved WT-only M0/M1 capacity attempts
are invalid. A stronger 25-start challenge found lower numerical minima for
some rows—for coupled `0.4/0.4`, raw error improves from `5.40e-5` to
`3.37e-5`—but every challenged topology still approaches height `200 um`,
fails exact closure, and misses WT pH by more than eight SEM. Thus the saved
optimizer coordinate is not unique, while the finite-domain gate conclusion
is robust. Pure pump and pure K continuation failures are retained as solver
failures and are not used to exclude remote branches.

For acid-base coupling, the declared NHE1, carbon hydration, and buffer
dissociation signatures have exact full rank three. Their direct `(Cl,H)`
projection has rank one and flux nullity two, even before unknown total carbon,
buffer pool, and boundary exchange are introduced. The mass-action structural
block conserves carbon, buffer sites, and charge exactly. It is not a
calibrated dynamic M2. The unchanged knockout NKCC1 and NHE assay results
restrict arbitrary genotype-specific gains but do not establish exact equality
of every native whole-cell flux.

### Dynamic and held-out audit

Round 6 did not stop at a failed integrator. A candidate-independent knockout
IVP was run with Radau and BDF at `rtol=1e-8`, `atol=1e-10` on the same
452-point code grid. The endpoint native-flow values differ by approximately
`8.3e-15`, post-step code integrals by `3.6e-11`, and the maximum shared-grid
flow difference is `1.28e-9`; conservation remains below `4.3e-14` and all
states stay positive. This verifies numerical integrability only.

No candidate has a valid WT denominator. In addition, the historical code
time has no certified minute conversion, the water output has no certified
absolute-flow scale, the Ca step is not a dose-to-input reconstruction of the
matched stimulation protocol, SR6 uses an immediate PKA input, and carrier
occupancies have no dynamic memory. It would therefore be leakage to fit a
delay or to relabel code time as the observed first two-to-three minutes. The
10-minute cumulative ratio, endpoint ratio, early comparability, sustained
late deficit, WT absolute flow, and AE2 dynamic result are correctly `N/A`.

### Why the unresolved conclusion is not merely a solver failure

The bottleneck follows from four results of different kinds:

1. an exact algebraic identity makes the fixed-chassis knockout field
   candidate independent;
2. exact target-state dimension and a reproduced rank-five source-family
   calculation show that two released observations cannot select the missing
   balance;
3. exact module-signature ranks show that the source-supported membrane split
   is invisible in the released intracellular projection, while the acid-base
   block retains a two-dimensional flux nullspace; and
4. independent CTMC, nonlinear-root, exact-rank, sparse-projection, and
   dual-solver calculations reproduce the controlling numerical results.

The frozen checklist asked for two source-consistent, non-secretion-calibrated
models with different held-out secretion predictions. That particular test is
**not licensable here**: zero model has a valid WT denominator, so manufacturing
two ratios would promote failed rows and weaken rather than strengthen the
audit. Decision relevance is established instead by the exact null spaces:
membrane-split current measurements determine two directions absent from the
cell-only projection, while carbon/NHE/buffer measurements determine at least
two acid-base directions left unresolved by Cl/pH. The missing data change
which chassis family can be calibrated; they do not merely tighten an already
accepted prediction.

### Executable discriminating experiment

Use WT and `Ae4-/-` mouse submandibular acini under the matched `0.3 uM`
carbachol plus `5 uM` isoproterenol stimulation. Measure, on a common physical
time axis:

- intracellular Na, K, Cl, pH, and cell volume;
- luminal or collected-effluent Na and K;
- membrane-sided ouabain-sensitive Na/K-pump current and Ca-activated K
  current;
- cAMP/PKA activity and AE4-associated uptake;
- total inorganic carbon plus an NHE- or buffer-flux readout; and
- simultaneously scaled gland flow for the physical time/flow map.

Use WT data to freeze membrane fractions, capacities, carbon/buffer parameters,
input kinetics, and unit mappings. Then predict the knockout trajectory without
changing a parameter. If carbon cannot be co-measured, the cation/current plus
luminal-ion protocol remains a decisive first branch test but is not a complete
cation-versus-acid-base discriminator. This wording is important: omitting
luminal Na/K would leave two coordinates of the released target manifold
unobserved, and omitting carbon flux would redirect rather than identify the
acid-base residual.

### Source and historical-claim discipline

The transporter doses and kinetics are kept in their heterologous assay
context; parotid topology evidence is not promoted to quantitative
submandibular fractions; the 2025 state order remains a hypothesis constrained
by mutagenesis rather than a directly observed unique sequence. Carbonate is
not claimed absent, and the PKA paper is not claimed to provide a direct S173
phosphorylation rate or gland-level delay.

All exclusions are scoped to the implemented QSS state graphs, declared
parameter/root domains, tested fixed or nested chassis, and available evidence
gates. The audit does not conclude that AE4 itself is insufficient, that AE4 is
biologically dispensable, or that the historical article was wrong.

### Final checklist resolution

| Frozen checklist section | Final status | Remaining caveat |
| --- | --- | --- |
| A. Holdout firewall | PASS | strict reference rows exist only after calculations and remain unjoined |
| B. State families/nesting | CONDITIONAL | genuine state graphs, but QSS whole-cell effect is exactly a static two-branch map |
| C. Stoichiometry/thermodynamics | PASS | effective dose fits do not identify microscopic order; carbonate remains open |
| D. Fixed-chassis Stage A | PASS AS NEGATIVE GATE | bounded result only; asymptotic remote branch prevents global uniqueness |
| E. Missing-balance localization | PASS AS NONIDENTIFIABILITY | source family is local numerical rank five, not a global manifold theorem |
| F. Minimal chassis reconstruction | FAIL AS PHENOTYPE REPAIR | M1 is nested and source-supported but uncalibrated and far from Cl/pH; M2 is structural only |
| G. Dynamics/held-out secretion | N/A BY VALIDATION GATE | two-solver KO diagnostic is not a prediction |
| H. Source/historical discipline | PASS | claims remain scoped to implementation and evidence context |
| I. Irreducible experiment | CONDITIONAL PASS | complete experiment needs lumen Na/K and carbon/flux readout; without carbon it is a first branch test |

The unresolved classification is therefore supported. No successful mechanism
has survived because no successful mechanism has been claimed.
