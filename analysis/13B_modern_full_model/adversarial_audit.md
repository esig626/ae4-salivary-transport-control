# Task 13B adversarial scientific audit

## PRE-REVEAL FREEZE

This checklist was frozen at `2026-08-27T20:29:45Z`, before this auditor inspected or evaluated any held-out AE4-null stimulated-secretion magnitude or time-course result. The held-out target file was treated as opaque; only its SHA-256 digest was computed. Its pre-reveal digest is `939f684ab92371f35b68fb52fbac59ebcb654fb5bf9a657feb4fd543ae1a0587`.

The repository base commit at freeze was `e831a77df803478c6731d7b4c39857761fa616d1`. The then-present `src/modern_full_model/*.py` digest ledger had aggregate SHA-256 `7c891e7415a9adf9872ffc5114faa4dd384be5dc7753e744a6c80596c4bc0acc`; calibration and non-secretion validation target digests were `bfaea270dd3c018f8e46db309695bb797b12df3aff4047bf69157d0a6589f3ff` and `3ea08c7901c4c76d5cdded089d217dcfddbe8e59e03b87e4f24c121f0acedbb5`, respectively.

The criteria below are immutable. Later work may append evidence and change only each criterion's disposition (`PASS`, `FAIL`, `BLOCKED`, or `NOT_APPLICABLE`). A criterion may not be deleted, weakened, or reworded after held-out reveal. Any post-freeze model or target change must be recorded with its reason and digest before that model is evaluated.

The criteria-only/pre-disposition checklist snapshot had SHA-256 `639472395a7ffda75e1f200fcb2bd7efe29528a6f44c791303f350290a085158`. After changing only the permitted `pre_reveal_status` field, the pre-reveal disposition snapshot has SHA-256 `26324d031fa3b28d1a1463906e8d024fbae85221ac25c6cfec75cac63a69f23c`.

## Frozen falsification criteria

| ID | Frozen criterion | Failure condition |
|---|---|---|
| H01 | No held-out AE4-null secretion value or qualitative timing may enter an objective, bound, prior, branch/root choice, input waveform, stopping rule, model-generation choice, or code constant before freeze. | Any such dependency, including an indirect target-derived choice, is found. |
| H02 | cAMP/PKA/AE4 regulatory gains, delays, and time constants must be fixed from WT/regulatory evidence independently of knockout secretion. | A regulatory parameter or family is selected using knockout secretion behavior. |
| H03 | S173A evidence must be represented only as dependence of the PKA-mediated effect unless direct phosphorylation is independently demonstrated. | Code or scientific claims treat S173 as directly phosphorylated fact without such evidence. |
| H04 | Static and dynamic regulatory families, and common-capacity versus state-specific modulation, must be compared as nested alternatives; extra states/rates require discriminating WT evidence. | A more flexible family is retained because of fit flexibility or knockout behavior rather than WT discrimination. |
| H05 | Every state, parameter, rate, current, molar flux, membrane area, compartment volume, and time conversion must have compatible physical units; current-to-moles conversion must include charge number and Faraday's constant exactly once. | A dimensional mismatch, missing/double Faraday factor, or uncertified code-time-to-physical-time label remains in an accepted result. |
| H06 | Transport source directions and apical/basolateral membrane signs must follow a declared compartment convention and reproduce independent reversal/equilibrium checks. | A sign is arbitrary, convention-dependent without a test, or gives the wrong reversal direction. |
| H07 | The acid-base block must account for intracellular and luminal inorganic carbon, buffer binding, bath exchange, and advective loss without creating/destroying carbon internally. | Carbon residual exceeds the declared numerical tolerance when boundary fluxes are included. |
| H08 | Species stoichiometry and membrane-current accounting must conserve charge, including electroneutral AE2/AE4/NKCC cycles and electrogenic pump/channel terms. | A nominally electroneutral cycle carries net charge, or charge/current residual exceeds tolerance. |
| H09 | Water flow must obey the declared osmotic/pressure law, use compatible volume units, preserve positive volumes, and close compartment water balance including outflow. | Water is created/lost internally, volume becomes nonphysical, or scaling is dimensionally inconsistent. |
| H10 | Apical and basolateral membrane current closure must each be demonstrated, rather than inferred from a single whole-cell sum. | Either membrane has a persistent unbalanced current in a state claimed as closed/steady. |
| H11 | Apical pump and K-channel fractions must be labelled unmeasured unless a primary source quantifies them; they must conserve total capacity and include the historical/baseline nesting limit. | A localization fraction is presented as measured, total capacity changes silently, or no exact nesting limit exists. |
| H12 | WT root selection must use deterministic multistart over the full declared physiological domain, retain/report all admissible branches and remote roots, and use a predeclared selection rule. | Only a convenient local root is reported, a remote physiological root is ignored, or selection changes after knockout evaluation. |
| H13 | Decisive WT and genotype trajectories must reproduce with two independent stiff solvers, tighter/looser tolerances, and nearby initial states. | Classification or qualitative behavior depends on solver, tolerance, step controls, or an unreported initial condition. |
| H14 | No accepted mechanism may depend on fitted parameters at arbitrary search bounds; profiles or perturbations must distinguish evidence-supported constraints from bound compensation. | A gate passes only at a nonphysiological/arbitrary bound or correlated parameters can compensate without disclosure. |
| H15 | Added modules must have source-matrix directions capable of correcting the documented residual; parameter fitting may not hide a wrong or rank-deficient topology. | Residual correction is achieved solely by compensating parameters in a source direction that cannot independently affect the claimed balance. |
| H16 | AE4 deletion must remove every AE4 transport and regulatory contribution while leaving unrelated parameters and inputs unchanged. | Any AE4 term remains, or another pathway/parameter is genotype-retuned. |
| H17 | AE2 deletion must remain an independent validation and must not acquire a large secretion role inconsistent with primary biology. | AE2 is used as compensatory fit freedom or gives an unsupported large output phenotype. |
| H18 | Every retained complexity increment must beat or falsifiably differ from its simpler exact nested model on allowed WT/transporter evidence. | A simpler nested model passes the same gates and the complex mechanism is nevertheless claimed as required/identified. |
| H19 | Qualitative conclusions must persist over source-supported uncertainty, including unmeasured apical fractions, weak kinetic parameters, and initial luminal states. | The classification is a fragile point prediction or reverses within allowed uncertainty without disclosure. |
| H20 | Claims about the 2018 model must remain scoped to the inherited historical-lineage implementation and tested domain; AE4 biology must not be judged from chassis failure. | The work states/implies that the paper was wrong, the old model is globally invalid, or AE4 itself is insufficient. |
| H21 | The physical time map, calcium/beta inputs, and integration horizon must be fixed using WT or independent regulatory/protocol evidence before knockout secretion reveal. | Any timing map or stimulus feature is inferred from the held-out knockout curve. |
| H22 | A valid WT resting root and WT stimulated trajectory on a physical time axis must exist before any normalized KO/WT secretion statistic is interpreted. | A normalized ratio masks an invalid denominator, failed WT gate, or uncertified time axis. |
| H23 | Unmeasured luminal states, lumen volume/outflow conditions, and apical fractions must be varied or profiled rather than silently fixed to produce a desired genotype result. | An unmeasured coordinate controls the result but is fixed without provenance or sensitivity analysis. |
| H24 | Every equation and every final mechanistic claim must carry the required provenance/claim category and distinguish measurement, derived constraint, numerical result, historical record, and new assumption. | An unsupported modeling decision is presented as measured or a numerical fit as mechanistic identification. |

The machine-readable canonical copy is `results/13B_modern_full_model/adversarial_checklist.csv`. The pre-reveal disposition for every criterion is `NOT_EVALUATED`; this records lack of exposure, not a pass.

## Evidence log (post-freeze append-only)

### 2026-08-27 pre-reveal scaffold audit

- `PYTHONPATH=src python -m unittest tests.test_modern_adversarial -v`: 10/10 focused adversarial tests passed. The tests hash the held-out file as opaque bytes, reject AE4-null secretion rows from the calibration ledger, reject sensitive-module dependencies on sealed/legacy holdout artifacts, exercise passive-flux signs, two-membrane closure, charge/carbon/water identities away from the reference point, test extreme apical partitions, and enforce assumption/root labels.
- The sealed held-out digest remained `939f684ab92371f35b68fb52fbac59ebcb654fb5bf9a657feb4fd543ae1a0587`.
- The evidence agent completed a separately logged pre-optimization evidence freeze at `2026-08-27T20:31:36Z`. Calibration and validation ledgers changed after this checklist's earlier snapshot to digests `24c86b1ed0353ab836f3526144d66a80e785da22e8ee8e0af8eb15ec75ca04a7` and `d9615c74c2c72296c2dd242483196853bee778b3eb286c3a81bbf503764287b0`; `reveal_log.json` records the reason and keeps `heldout_reveals` empty. This is a pre-optimization evidence-ledger revision, not a held-out reveal.
- The present whole-cell implementation passes local algebraic conservation tests, but it is explicitly a scaffold with placeholder parameters. These passes do not establish a physiological WT model.
- The bundled `solve_resting_root` has a conserved cell-charge left-null direction, does not parameterize the electroneutral charge manifold, uses a non-space-filling diagnostic start design, reports only one best residual point, and does not track remote branches. It is now explicitly and unconditionally marked `production_eligible=False`. H12 and H22 remain blocked until the separate production root workflow exists and passes.
- At the executable starting parameters, the diagnostic root attempt failed (`max_abs_scaled_rhs` about `4.09e-2`) and reached its broad upper bound; its returned state carried a large nonzero cell bulk charge. It is not a WT result and may not be used for model selection.
- The finite intracellular buffer pool is present in acid-base speciation but is not separately included in the osmolarity call. The model must state and test whether it is already included inside `cell_impermeant_osmoles_fmol`; otherwise osmotic particles are omitted or double-counted. H05/H09 remain blocked on this definition.
- `FullModelParameters` at this snapshot validates finiteness and the two fractions but not the physical sign/domain of many conductances, permeabilities, capacities, geometry quantities, and constants. Calibration must not be allowed to exploit negative or otherwise invalid values. H05/H14 remain blocked pending domain guards and tests.
- The Na/K-ATPase law is an ATP-replete irreversible closure without ATP/ADP/Pi or pump reversal. This is permissible only as an explicit open-energy modeling decision within the tested regime; it cannot be counted as a detailed-balance/reversal validation of the pump.
- The current Radau/BDF smoke test compares two algorithms in SciPy over 0.1 s. It is useful but is not the independent full-horizon solver reproduction required by H13.

## Pre-reveal disposition (not the final held-out audit)

Seven criteria pass at the current pre-reveal stage: H01 firewall, H03 S173 claim discipline, H06 signs/reversals, H07 carbon accounting, H08 charge accounting, H10 two-membrane current closure, and H20 historical-claim discipline. Seventeen criteria are `BLOCKED`, not failed, because the regulatory integration, production WT calibration, dynamic/time freeze, uncertainty profiles, genotype validation, and independent reproduction do not yet exist. `final_status` remains `NOT_EVALUATED` for all criteria.

### Blocking items before any held-out access

1. Produce a physiological WT resting root on an explicitly electroneutral charge manifold with non-collinear deterministic multistarts, all root clusters, rank/SVD evidence, branch/boundary distances, and a predeclared selection rule (H12, H22).
2. Integrate the smallest defensible dynamic cAMP/PKA/AE4 family and freeze its gains/time constants entirely from WT/regulatory evidence. The primary ledger states that the 2021 experiments provide no beta-onset, cAMP, PKA, phosphorylation, washout, or fitted regulatory time trace; any kinetic range is therefore a declared assumption, not an identified slow scale (H02, H04, H21).
3. Complete the dimensional/time ledger, physical-domain parameter guards, whole-gland geometry/flow scaling, and unambiguous buffer-osmole accounting (H05, H09).
4. Profile unmeasured apical pump/K fractions, luminal coordinates/outflow, regulatory kinetics, and parameter bounds; compare exact simpler nests and source-direction rank before calling any added mechanism required (H11, H14, H15, H18, H19, H23).
5. Freeze and reproduce the full WT trajectory with genuinely independent numerical work, tolerance/initial-condition sweeps, then verify exact AE4 deletion after regulatory integration and the independent AE2 gate (H13, H16, H17).
6. Complete equation/parameter/claim provenance and immutable G5 manifests before reveal (H24).

### Nonblocking warnings for the scaffold

- The ATP-replete Na/K-pump closure is irreversible and omits ATP/ADP/Pi and pump reversal. Keep it explicitly scoped as an open-energy approximation; do not cite it as a detailed-balance validation.
- The diagnostic root helper can return a non-electroneutral residual minimizer at a broad bound. It is safely labelled `production_eligible=False`, but downstream code must reject it categorically.
- A Radau/BDF smoke test within SciPy is not an independent implementation.
- Bath chloride and fixed cell anion are derived for the default tuple. If upstream bath or initial acid-base coordinates change, they must be rederived or consistency-checked. The proposed production calibration avoids upstream bath changes and reconstructs dependent chloride coordinates; that design still requires executable verification.

## Final adversarial synthesis

The final audit accepts the firewall, equation provenance, conservation,
sign/reversal, charge/current, carbon, water, physical-time, production-root,
branch, solver/tolerance, nearby-state, nesting, and historical-claim controls.
It rejects scientific promotion under H09/H22 because every native dynamic
root fails the independently frozen one-SMG absolute-flow gate. H16/H17 are
not evaluated on the replacement lineage because genotype work is prohibited
after that failure. This is correct fail-closed behavior, not missing work.

The two disclosed pre-reveal plaintext incidents remained quarantined from
model generation and did not alter objectives, bounds, topology, root choice,
stimuli, kinetics, tolerances, or stopping rules. No authorized reveal
occurred, and `heldout_reveals` remains empty.

The final Outcome 4 claim is limited to the declared source/topology/capacity
hierarchy. Its single decision-critical measurement is absolute 600-s fluid
secretion per imaged WT acinar cell under matched CCh+IPR, with 10-s
resolution. That measurement selects whether the factor-of-two absolute-scale
failure lies in the cellular flux equations or the cell-to-gland observation
map. Reduction, GSPT, bifurcation, identifiability, and manuscript gates remain
closed.
