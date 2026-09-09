# Final answer: state-resolved AE4 cycle and whole-cell reconstruction

`PHENOTYPE UNRESOLVED — SPECIFIC MISSING EXPERIMENT OR FLUX IDENTIFIED`

This classification follows the full required hierarchy. State-resolved AE4
families were built and constrained; every conditional transporter survivor
was tested on the fixed chassis; Stage A failure triggered Stage B; the
missing-balance family was localized; the smallest source-supported cation
topology and a conservation-correct acid-base block were tested; and dynamics,
held-out validation, independent reproduction, and adversarial challenges
were completed. The conclusion is not that AE4 is biologically insufficient.
It is that the current evidence does not identify the whole-cell balance that
maps loss of AE4 into the measured secretion phenotype.

## 1. What is the smallest AE4 state-resolved mechanism compatible with the transporter-level evidence?

There are two levels of answer.

- **Smallest source-supported state graph:** SR2-112, one conserved carrier
  pool with explicit inward/outward states and separate reversible Na and K
  `1:1:2` branches sharing conformational states. It is the smallest graph
  that represents one mixed-cation transporter, direct transport of both
  cations, electroneutrality, reversal, cation-specific mutation limits, and
  finite-pool Na/K competition.
- **Smallest conditional realization that also expresses the reported
  cooperative dose summaries in the declared assay context:** SR2-123 at the
  slow common-transition gauge. It uses the unmeasured `1:2:3`
  stoichiometric alternative to generate near-quadratic cation dependence.
  SR5-112 is an observationally viable alternative: it preserves `1:1:2` net
  transport while adding an explicit nontransported catalytic-cation state.

Thus no unique microscopic mechanism is identified. Pure-cation dose curves
do not discriminate two transported cations from one transported plus one
catalytic cation, nor do they identify mixed-bath competition barriers. The
five numerical “passes” are conditional expressivity rows with an explicit
assay observation map, not formal likelihood winners. See
[`state_cycle_models.md`](state_cycle_models.md) and
[`transporter_fit.md`](transporter_fit.md).

## 2. Does cation-specific state structure explain anything that the Task-12 coarse-grained laws could not?

**Yes at transporter level, but not the held-out phenotype.** A shared carrier
pool gives exact cross-branch competition and permits a pure Na/K slip:

\[
J_{Cl}=J_{Na}+J_K,\qquad J_{Na}=-J_K\ne0
\quad\Longrightarrow\quad J_{anion}=0.
\]

That circulation cannot be represented by one scalar AE4 event rate and was
not present in the chosen Task-12 parametric families. It also makes one
branch depend on the other branch's rates through a common carrier
denominator.

The limitation is equally exact: carrier occupancies are solved quasi-steadily
inside each whole-cell RHS call. The whole-cell model therefore sees two
algebraic functions `J_Na(x,PKA)` and `J_K(x,PKA)`, not new dynamic carrier
states. A sufficiently flexible pair of scalar functions can reproduce the
same steady source, and the implemented state graph supplies no relaxation
memory or physical delay. Cation state resolution earned transporter-level
structure, not a delayed knockout secretion mechanism.

## 3. Can a source-supported state-resolved AE4 model reproduce the AE4-null secretion phenotype on the original best chassis?

**No.** Five conditional transporter/gauge rows entered Stage A. All five
failed before an admissible secretion prediction:

- `0/5` valid interior WT capacity calibrations;
- `0/5` WT physiology Gate-2 passes;
- joint WT-Cl attempts approached the derived `height=200 um` boundary with
  raw balance residuals near `3.05e-5` to `3.08e-5`, versus the required
  `1e-8`; and
- exact roots at the frozen attempted capacities had pH about `7.37845`,
  approximately `+6.69` reported WT SEM from `6.91 +/- 0.07`.

After deletion, every candidate has the same AE4-null vector field. Its common
bounded root is `Cl_i=48.782215 mM`, pH `7.378474`, and model-reference volume
`2.974178 pL`, with raw residual `5.77e-14`. It misses the released knockout
Cl and pH by `+7.676` and `+24.424` reported SEM. A remote positive root at
height around `10^14--10^15 um` is physiologically inadmissible but prevents a
global uniqueness claim. The valid statement is failure on the documented
physiological domain, not failure of AE4 biology.

No KO/WT secretion ratio was formed from a failed WT denominator.

## 4. If not, what exact balance/source direction is missing in the AE4-null state?

The released evidence does **not** define one exact vector. It defines an
exact five-dimensional family.

With knockout `Cl_i=36.50 mM`, pH `6.89`, and the retained electroneutrality
reduction,

\[
HCO_{3i}=h_i^*-36.50-X/H+Na_i+K_i,
\]

while `(Na_l,K_l,H,Na_i,K_i)` remain unmeasured. For each allowed state `y`,
the required correction is

\[
d(y)=-F_0(y).
\]

The derivative of the determined correction coordinates with respect to the
five free state coordinates has five nonzero singular values

`(0.40223543, 0.01887636, 0.00934814, 0.00163360, 0.000449923)`.

The required source therefore rotates with numerical rank five across states
that are equally compatible with released Cl/pH. No one- or two-module cell
signature is exact on both declared conditional slices, and the best pair
changes with the unmeasured coordinates. The exact missing object is the
state-dependent family of Na, K, Cl, acid-base, lumen/current, and water
corrections—not a unique protein vector.

In the expanded topology, the specifically unobserved cation-support
directions are

\[
\delta\sigma_P=(0_i;+3,-2,0,0;+1,-1,0),
\]

\[
\delta\sigma_K=(0_i;0,+1,0,0;+1,-1,0),
\]

for apical-versus-basal pump and K pathways. They have exact rank two in
lumen/current space and rank zero in intracellular chemistry.

## 5. Which known module or smallest module combination supplies that direction?

No known one- or two-module combination is identified as supplying the full
rank-five family.

The smallest independently supported omitted topology is **M1**, a conserved
redistribution of total Na/K-ATPase and total Ca-activated K conductance into
apical and basolateral copies. It supplies the two expanded lumen/current
directions above. The implementation uses local luminal K for the apical
pump, local membrane voltages for K currents, complete lumen Na/K terms, and
re-solves both membrane-current equations. It is exactly nested in the fixed
chassis at zero apical fractions.

M1 does not supply the invariant bicarbonate residual at either conditional
target slice. The smallest source-complete next candidate is therefore **M3:
M1 plus conserved CO2/HCO3 and buffer chemistry with shared,
genotype-invariant NHE1**. M3 has the required kinds of directions, but its
carbon totals, buffer pool, boundary fluxes, and cation localization fractions
are not measured. Released Cl/pH observe only rank one of the three declared
acid-base reaction directions, leaving flux nullity two before boundary terms
are added. M3 is a source-supported next chassis, not an identified repair.

## 6. Does adding that module, calibrated without AE4-null secretion, reproduce the secretion magnitude?

**No.** The cation split was tested at predeclared `0`, `0.2`, `0.4`, and
historical area-sensitivity fractions without reading knockout secretion.
Using frozen SR2-123 transporter parameters, M0, pump-only, K-only, and two
coupled split WT capacity attempts all approached height `200 um`, retained
raw closure error `3.05e-5--5.40e-5`, and had pH `7.520--7.534`. There were
zero valid capacities and zero WT Gate-2 passes.

In the knockout ionic sensitivity, the best displayed split (`f_P=0.4`,
`f_K=0.2`) moved the root only to `Cl_i=47.6512 mM`, pH `7.36742`, still
`+6.969` and `+23.871` reported SEM from the released targets. Because no M1
row provides a valid WT denominator and M3 cannot be independently calibrated,
the held-out secretion magnitude is `N/A`, not a failed or fitted ratio.

## 7. Does it reproduce the early-comparable then sustained-deficit time pattern?

**No admissible model can be scored.** There is no valid WT trajectory,
physical code-time-to-minute map, absolute flow calibration, measured PKA
activation time constant, or dynamic carrier occupancy. The current SR6 input
is an immediate code-coordinate step, so it cannot be relabelled as the
observed two-to-three-minute separation.

A candidate-independent knockout-only IVP was integrated with Radau and BDF
as a solver diagnostic. The solvers agree and conserve the declared balances,
but the trajectory has no WT denominator and no physical time axis. It is not
a held-out prediction. Endpoint, integral, first-2--3-minute, and sustained
late metrics remain explicitly `N/A` in
[`dynamic_trajectory_status.csv`](../../results/13_state_resolved_ae4/dynamic_trajectory_status.csv).

## 8. Which features are necessary, sufficient, and merely one parameterization?

| Status | Features or conclusions |
| --- | --- |
| **Necessary/exact within the declared model** | Carrier conservation; branch electroneutrality; local detailed balance; shared-pool competition/slip identity; genotype-zero annihilation of every AE4 contribution; candidate-independent fixed-chassis knockout field; five free knockout state coordinates after releasing Cl/pH; rank-two membrane-split signatures in expanded space and rank zero in cell chemistry. |
| **Necessary numerical gates** | An interior WT capacity, WT Cl/pH physiology, a source-calibrated chassis, positive/conserving dynamics, and a physical time/flow map before any secretion comparison. None is optional. |
| **Sufficient for the current unresolved classification** | Zero valid Stage-A and Round-5 WT rows; candidate-independent bounded null mismatch; rank-five required-source family; uncalibrated split/carbon directions; no licensed denominator or physical timing map; independent and adversarial reproduction of those facts. |
| **Not sufficient** | A conditional transporter dose fit; one good source cosine; a two-output local parameter fit; moving knockout Cl/pH in the desired direction; a successful KO-only integrator run; qualitative apical pump/K localization. |
| **Parameterization/sensitivity only** | The five headline transporter/gauge passes; fitted loaded-energy lumps and assay scale; imposed mutant branch limit; `f_P/f_K=0.2/0.4/0.398`; immediate PKA step; conditional Stage-B slices; the historical water/time/flow scale. |

No state-resolved mechanism, cation split, or acid-base block is presently
shown sufficient for the whole phenotype.

## 9. Are multiple successful mechanisms observationally equivalent?

There are **no successful whole-cell mechanisms** to call equivalent.

At transporter level, SR2-123 and SR5-112 are conditionally observationally
equivalent under the available pure-cation dose summaries, and pure-bath data
cannot distinguish SR1 from the corresponding SR2 branch. At chassis level,
many source or local-response decompositions are compatible with two released
ionic outputs, but none is a validated phenotype reconstruction. This is
nonidentifiability among failed or uncalibrated explanations, not multiple
successful models.

## 10. What single experiment would best discriminate them?

Run one predeclared, matched-tissue balance protocol in WT and systemic
`Ae4-/-` mouse **submandibular** acini (with a paired gland-flow arm) under the
exact `0.3 uM` carbachol plus `5 uM` isoproterenol stimulus for 10 minutes.
Measure, on the same time grid:

1. intracellular Na, K, Cl, pH, and three-dimensional cell volume;
2. luminal or effluent Na and K, so the free `Na_l,K_l` coordinates are not
   silently fixed;
3. apical and basolateral ouabain-sensitive pump current and Ca-activated K
   current, yielding `(P_a,P_b,J_Ka,J_Kb)`;
4. a cAMP/PKA activity trajectory and AE4-associated uptake; and
5. total inorganic carbon or CO2 flux plus an NHE/buffer-flux readout in a
   synchronized arm.

Use WT data only to freeze topology fractions, capacities, carbon/buffer
parameters, time constants, and physical unit maps. Then predict all knockout
ions and flow without retuning. The discrimination rule is predeclared:

- measured nonzero membrane-partitioned cation support that closes the
  extended balance and predicts the knockout ionic displacement under frozen
  carbon supports the M1 branch;
- normal or insufficient cation support with a resolved carbon/NHE residual
  supports the acid-base branch; and
- neither result moves the problem to a specifically measured remaining flux,
  rather than another secretion-only fit.

If carbon/NHE is not co-measured, the protocol is the first decisive
cation-topology branch test, not a complete mechanism identifier. The strict
secretion curve remains out of sample.

## 11. Which conclusions are facts, model-derived results, or new assumptions?

| Claim class | Content |
| --- | --- |
| **Direct experimental facts** | `35 +/- 4.7%` less total 10-minute saliva in AE4-null glands; first 2--3 minutes comparable and later sustained deficit; knockout resting `Cl_i=36.50 +/- 1.60 mM` and pH `6.89 +/- 0.02`; direct Na and K transport; electroneutral/voltage-independent behavior; effective Na/K EC50 and Hill summaries; reversal; PKA/H89/S173 association; T448/T756 cation nonidentity and double-mutant Na-loss/K-retention hierarchy; qualitative apical pump/K topology in parotid. |
| **Exact model deductions** | Cycle charge and affinities; local detailed balance; shared-pool/slip identities; SR3 `1:1:1:1` electroneutrality and conditional feasibility; candidate-independent knockout field; target-manifold dimension five; exact module-signature ranks; M1 nesting and current/lumen closure; carbon/buffer conservation. |
| **Numerical model results, scoped to declared domains** | Five conditional transporter rows; zero valid Stage-A and Round-5 WT capacities; common bounded null root and ionic errors; finite-grid and finite-optimizer wrong-sign chloride-source result; M1 knockout sensitivities; rank-five local correction family; two-solver knockout-only diagnostic; independent CTMC/root/rank reproduction. |
| **New assumptions/sensitivities** | Intracellular heterologous-assay tuple; observation-map Rmin nuisances and common scale; loaded-state energy lumps; discrete barrier gauges; catalytic-cation state; U/P layer and PKA edge placement; imposed mutant branch limit; QSS carrier occupancy; M1 fraction grid; historical seven-state baths, water, time, and flow scales; two conditional target slices. |
| **Unmeasured hypotheses retained** | Exact AE4 stoichiometry; bicarbonate versus carbonate; 2025 sequential state order; direct S173 phosphorylation; physiological PKA timing; submandibular apical pump/K fractions; conserved-carbon/NHE/buffer parameters; knockout Na/K/volume/lumen coordinates. |

The rounded “30%” is a later narrative restatement, not the primary endpoint;
the historical `24%` is a model output, not experimental evidence.

## 12. Does the result justify a new revisit paper, and what is the precise contribution without saying the old article was wrong?

**Yes, as a reconstruction and experimental-design paper—not yet as a solved
phenotype-mechanism paper.** Its precise contribution would be:

1. translate modern direct Na/K, mutational, carbonate, and PKA evidence into
   thermodynamically exact, state-resolved AE4 candidate graphs;
2. show which transporter features are identifiable and which reduce to
   steady scalar functions under current protocols;
3. demonstrate that all conditional state-resolved survivors fail the fixed
   historical-lineage chassis at the WT gate and share one knockout field;
4. prove that released knockout Cl/pH define a five-dimensional target family,
   preventing a unique missing-module claim;
5. implement and falsify the smallest source-supported cation-topology repair
   in a transparent finite domain; and
6. derive the matched-gland cation/lumen/carbon experiment that can close the
   balance before secretion is tested out of sample.

That framing neither says the 2018 article was wrong nor that AE4 is
insufficient. It says the historical comparison chassis was valuable but is
not source-complete enough, under current evidence and the tested domain, to
map molecular AE4 loss to the full 2015 secretion phenotype.

## Reproducibility and audit disposition

The independent route uses SVD and matrix-tree CTMC solvers, an independently
transcribed null field with a different nonlinear root parameterization, and
exact rational rank calculations. Maximum discrepancies are `5.71e-12` for
CTMC occupancy, `1.32e-12` for edge current, `1.67e-8` across null-root state
components, and `7.41e-13` for the saved rank-five correction Jacobian. The
adversarial audit also challenged joint versus sequential transporter fitting,
nine intracellular assay contexts, eight-decade barrier domains, remote roots,
and 25-start Round-5 WT fits. None supplied a validated phenotype model.

No Round-10 inference scan was performed because the prompt permits it only
after a physiologically validated model survives held-out validation and the
adversarial audit; the eligibility set is empty.

## Repository provenance note

The evidence freeze records local materialization commit `7846706645f9`.
The corresponding GitHub work was based on the requested branch head
`5a0eec995f244762ff138f16f68d17cac6956901`; the published Task-13 commit is
a direct child of that GitHub commit. This transport mapping changes no
evidence value, holdout assignment, equation, or result hash.
