# Task 22 WT physiology contract — declared before optimization

Task 22 continues from Task 21 final commit
`db734f52b8d71d670fcf793beab922967212bdd6`. The unchanged
`POOLED_CATION_112_NO_SLIP` mechanism retains one reversible net cycle,
donor-side cation partition and the declared 1:1:2 working stoichiometry.
No new Task 22 optimization has run or been inspected at declaration.
The machine-readable contract is
`results/22_physiological_wt_domain/precalibration_contract.json`.

| Resting quantity | Inclusive hard domain | Evidence role |
|---|---:|---|
| Cell Cl | 47.10–53.10 mM | WT 50.10 ± 1.50 mM, ±2 reported SEM screen |
| Cell pH | 6.77–7.05 | WT 6.91 ± 0.07, ±2 reported SEM screen |
| Cell Na | 2–60 mM | Inherited broad physiological/provenance envelope |
| Cell K | 60–200 mM | Inherited broad physiological/provenance envelope |
| Cell volume | 0.3–5.0 pL | Inherited broad model envelope |
| Lumen Na | 100–200 mM | Broad primary-fluid envelope |
| Lumen K | 1–30 mM | Broad primary-fluid envelope |
| Lumen Cl | 80–180 mM | Broad primary-fluid envelope |
| Lumen pH | 6.8–8.0 | Broad acid-base envelope |
| Lumen TIC | 1–80 mM | Broad carbon envelope; TIC is not HCO3 alone |
| Lumen volume | 0.02–0.50 pL | Broad geometry envelope |
| Cell/bath osmolarity | 0.80–1.20 | Production-model osmolarity |
| Lumen/bath osmolarity | 0.80–1.20 | Production-model osmolarity |
| Each adjustable positive capacity/reference | 0.01–100 | Broad anti-escape domain |
| Fractions | 0–1 | Physical domain |
| Each absolute individual solute flux / positive basolateral Cl pool | ≤100 | Broad anti-cancellation screen |

## Provenance and interpretation

The WT chloride and pH targets and their SEM values retain the Task 21
directive and `analysis/13B_modern_full_model/wt_calibration.md`. These
are screening bands, not population confidence intervals. The unchanged
Na/K/volume rationale is recorded in the Task 21 physiology contract and
`model/parameters.md`: Na references 15.5–25 mM in transferred/model contexts,
K references 120–150 mM, and a published-model cell-volume convention of
1.3 pL. None supplies a matched native mouse-SMG uncertainty interval for
these three quantities. The inherited broad limits are explicit modeling
choices, not genotype-derived targets.

[Patterson et al. (2012)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3517652/)
describes the primary acinar secretion as isotonic and plasma-like, before
ductal modification. The project's later literature lineage includes
[Su et al. (2022)](https://link.springer.com/article/10.1007/s11538-022-01041-3),
whose discussion reports primary-fluid composition ranges of Na 140–150,
K 5–15, Cl 80–130 and HCO3 30–60 mM across its cited gland/species/stimulation
contexts (Mangos 1973; Young 1971). These published contextual ranges support
broad primary-fluid screens; they are not newly verified native mouse-SMG
measurement intervals. Task 22's specified envelopes are deliberately wider.

The [2018 model](https://pmc.ncbi.nlm.nih.gov/articles/PMC5792321/) and
`model/parameters.md` carry lumen Na 118.7, K 5.6 and Cl 124.3 mM from the
Palk/Mangos lineage. The older model's bath/cell/lumen osmolarities were
292.6/296.6/297.4 mOsm. They support approximate isotonicity, but are not
substituted for the current production model's computed osmolarities.
The lumen pH/TIC endpoints are prescribed broad acid-base choices. The
volume range brackets the inherited production reference 0.10 pL and
the old model's 0.02 × 1.3 = 0.026 pL convention. Those references are
provenance, not measured volume uncertainty.

All osmotic gates use `ModernFullModel` production observables, which call
`water.compartment_osmolarity_mOsm`. Cell finite buffer particles and other
impermeant osmoles are counted exactly once. No alternative hand-written
osmolarity sum is authoritative for acceptance.

## Capacity and signed-flux definitions

The fifteen inherited adjustable coordinates are AE4, NKCC1, AE2, NHE1,
apical/basolateral Na/K pumps, apical/basolateral K conductance, CaCC,
apical/basolateral CO2 permeability, and paracellular Na/K/Cl/HCO3
conductances. Every coordinate has a positive reference in all ten roots;
the JSON records every reference. Pump/K totals and fractions are derived
from the two nonnegative components. Their totals also satisfy the fold
box. All remaining fractions stay in [0,1]. No zero-reference exception
is declared and no zero-reference pathway is activated.

For each candidate, use signed **cell chloride sources** and define
`Lplus = max(J_AE4_Cl,0) + max(2*J_NKCC1,0) + max(J_AE2,0)`.
This is the Tasks 18/20 ledger semantics; the set of positive contributors
is recomputed, not assumed to remain the legacy set. Record every individual
cycle flux and every individual Na/K/Cl/TIC/alkalinity source component,
including stoichiometric pump and exchanger sources, channels converted
from current to fmol/s, both CO2 exchanges, paracellular sources and lumen
outflow sources. Apply the ≤100 ratio to each before net cancellation.
Report water separately in pL/s because an ion-flux ratio would be
dimensionally invalid. This threshold and the capacity box are broad
anti-pathology assumptions, not fitted biological constants.

## Acceptance, ranking and checkpoints

Every WT candidate must simultaneously satisfy all above bounds, finite
physical states and parameters, positive pooled affinity and strictly
positive chloride loading, and the no-opposing-Na/K invariant. It must pass
the inherited complete production RHS, current, state-charge and all
conservation/bookkeeping tests. The scaled independent resting residual
limit remains 1e-7, omitted-row raw residual 1e-9, state charge 1e-9 fmol,
and resting current 1e-18 A. The full conservation tolerance map is copied
verbatim in the JSON, including its stricter 1e-20 A current identities.
No tolerance is relaxed. All regulatory resting rows must also close.

All three 600 s production WT trajectories at calcium 0.10, 0.25 and
0.50 uM must pass the inherited numerical, physical and conservation gates,
plus state charge ≤1e-9 fmol and valid regulatory fractions. Use the inherited
production Radau settings and 1 s output grid. Resting physiology bands are
resting screens; dynamic acceptance uses the existing dynamic gates.

Use all ten inherited WT roots and every saved Task 21 resting witness as
numerical starts. Their states do not become targets. Additional WT-only
starts and purely numerical improvements are allowed. Preserve all failed
attempts. Rank admissible solutions lexicographically by normalized squared
WT Cl/pH error, largest absolute natural-log capacity departure, then summed
squared natural-log departure. Do not penalize state movement. Retain every
distinct admissible solution; deduplicate full state and parameter vectors
at componentwise relative tolerance 1e-6 with native-unit floor 1e-12.

Bath, chemistry, buffers, geometry, water/outflow, regulation, source signs,
mechanism and pathway inventory remain unchanged. No AE4 chloride-share
target, legacy flux matching, genotype objective or post-failure bound
expansion is allowed.

This contract must be **committed, pushed and remotely verified before
running or inspecting optimization**. A later complete WT checkpoint must
contain all attempts, decisions, states, parameters, ledgers, admissible
solutions and WT dynamic evidence. It must be committed, pushed and remotely
verified by SHA and manifest before **any Task 22 genotype calculation or
historical genotype regression**. If the WT ensemble is empty, stop before
genotype evaluation and create no genotype result or comparison table.

Distinguish numerical nonconvergence, physiological-bound conflict,
capacity-bound conflict, thermodynamic conflict, conservation/dynamic-gate
failure and a proved structural contradiction. A local optimizer failure or
an active bound is not a proof of impossibility. A proved contradiction must
state the equations and domain over which it applies.
