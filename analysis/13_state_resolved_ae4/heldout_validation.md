# Round-7 held-out secretion validation

## Decision

**No candidate is eligible for the held-out secretion test.** The tested set
contains zero valid WT capacity calibrations, zero WT Gate-2 passes, and zero
independently calibrated Round-5 chassis extensions. Therefore endpoint,
cumulative, onset, and late KO/WT secretion metrics are reported as `N/A`, not
as ratios computed from a biologically invalid denominator.

The strict targets are recorded below to preserve their meaning, but they are
not joined to a candidate prediction. This validation gate does not state or
imply that AE4 is biologically insufficient. It establishes that the tested
mechanisms have no validated whole-cell denominator on which the knockout
phenotype can be predicted.

## 1. Exact held-out target semantics

| Evidence | Frozen observation | Correct validation meaning |
| --- | --- | --- |
| `E15-02` | Total saliva collected over 10 min was `35 +/- 4.7%` less in AE4-null glands (`n=6`/genotype, mean +/- SEM) | primary integrated endpoint |
| algebraic complement | KO/WT total `0.65 +/- 0.047`; one-reported-SEM display interval `[0.603,0.697]` | a derived ratio, not separately fitted by the paper or a reconstructed confidence interval |
| `E15-03`, onset | first 2--3 min comparable | qualitative shape constraint; not a free-delay target |
| `E15-03`, late | sustained deficit after about 3 min; later pointwise flow less than half WT | pointwise late constraint, not the 10-min integrated percentage |
| frozen digitization | minute means from the raster, with about `+/-0.4 uL/min` reading uncertainty per mean in addition to undigitized plotted SEM | envelope check only; pointwise fitting prohibited |
| later “approximately 30%” wording | rounded 2018 narrative restatement of the 2015 phenotype | not a second endpoint and not “KO equals 30% of WT” |
| historical 24% reduction | 2018 model output | not experimental validation evidence |

The primary `35 +/- 4.7%` integrated loss controls. The rounded 30% wording
must never replace it.

## 2. Eligibility set

Stage A produced five conditional transporter-assay survivors. All five WT
capacity attempts were invalid/bound-adjacent, and their exact frozen-capacity
roots failed WT pH by about `+6.69` reported SEM. There were therefore zero
WT Gate-2 passes.

Round 5 then tested the M0 control, pump-only, K-only, and two coupled pump/K
splits with the transporter parameters frozen. All five WT attempts again hit
height near `200 um`, retained raw closure errors from `3.05e-5` to `5.40e-5`,
and had pH `7.520--7.534`. No row provided a valid capacity or an evaluable WT
gate. The split fractions are sensitivity values rather than measured mouse
submandibular fractions. The acid-base alternative is structurally
nonidentifiable from the released Cl/pH pair and has no independent parameter
calibration.

Consequently the eligible set for Round 7 is empty. A successful integrator
run cannot override that static exclusion.

## 3. Required held-out metrics

| Required metric | Result | Reason |
| --- | --- | --- |
| 10-min total KO/WT secretion ratio | `N/A` | no valid WT denominator; no physical time/flow map |
| code-endpoint KO/WT ratio | `N/A` | no valid WT denominator; the experimental 0.65 is an integrated total, not an endpoint target |
| cumulative post-stimulus ratio | `N/A` | KO-only code integral exists, but no WT integral and no 10-min mapping are admissible |
| first 2--3 min comparability | `N/A` | code time cannot be mapped to minutes |
| sustained late deficit | `N/A` | same, and no valid trajectory pair |
| delayed separation mechanism | `N/A` | immediate PKA input, QSS carrier occupancy, and no independently constrained delay state |
| WT absolute/relative secretion | `N/A` | absolute flow scale uncertified; failed WT model cannot establish a relative trajectory |
| AE2-null secretion | `N/A` | upstream WT gate fails; primary evidence supplies no numerical equivalence margin |
| AE4-null resting ions | fail on fixed chassis | common root `Cl_i=48.7822 mM`, pH `7.37847`, versus `36.50 +/- 1.60 mM` and `6.89 +/- 0.02` |
| AE4-null volume | no experimental target | `2.974 pL` is a model diagnostic only |
| conservation/two-solver check | pass for KO-only code diagnostic | Radau/BDF reproduce the numerical IVP, but this is not held-out validation |

Blank ratio fields in
`results/13_state_resolved_ae4/dynamic_trajectory_status.csv` are intentional
machine-readable `N/A` values. They must not be imputed as zero, one, or a
failed numerical solve.

## 4. Diagnostic that was not promoted

The candidate-independent AE4-null IVP was integrated from its frozen common
root under the historical Ca step. Radau and BDF agree on endpoint code-native
flow to `8.3e-15` and on the post-step code integral to `3.6e-11`; the maximum
shared-grid `q_total` difference is `1.28e-9`. Conservation checks remain below
`4.3e-14`, and every state stays positive.

Those facts establish numerical reproducibility only. The diagnostic has:

- `wt_denominator_valid=false`;
- `physical_time_certified=false`;
- `heldout_prediction=false`;
- blank endpoint and integrated KO/WT ratios; and
- `strict_target_joined_to_prediction=false`.

It cannot be compared with the 10-min integrated target or the first-2--3-min
shape. Full setup and limitations are in `dynamic_reconstruction.md`.

## 5. Model selection disposition

There are no successful candidates to rank by mechanism count, parameter
count, data sets explained, robustness, held-out prediction, or conservation.
Ranking failed rows by proximity to a strict secretion target would be target
leakage. Ranking them by a KO-only flow would merely rank absent denominators.

The smallest independently supported next topology remains distributed
apical/basolateral Na/K pump plus K conductance, coupled if necessary to an
explicit conserved-carbon/NHE/buffer module. Round 5 shows that the presently
uncalibrated split alone does not repair WT closure. It is therefore a next
experimental/modeling target, not an accepted mechanism.

## 6. Specific missing experiment and formal outcome

The decisive missing information is the stimulus-dependent, membrane-resolved
cation/acid-base support and its physical time scale in matched mouse
submandibular tissue. Under `0.3 uM CCh + 5 uM IPR`, measure in WT preparations:

- apical and basolateral ouabain-sensitive pump current;
- apical and basolateral Ca-activated K current;
- time-resolved intracellular Na, K, Cl, pH, and cell volume;
- time-resolved luminal or collected-effluent Na and K;
- a calibrated cAMP/PKA activity trajectory and AE4-associated uptake;
- total inorganic carbon plus an NHE- or buffer-flux readout; and
- simultaneously scaled gland flow.

Use those WT data to freeze topology fractions, total capacities,
carbon/buffer parameters, recruitment/PKA time constants, and the physical
time/flow map. Only then run the AE4-null trajectory without changing a
parameter and compare its integrated ratio and delayed envelope with
`E15-02/E15-03`.

If total inorganic carbon/NHE/buffer flux cannot be co-measured, the cation
and luminal-Na/K experiment is still a decisive first branch test, but it does
not by itself identify the complete cation-versus-acid-base mechanism.

The required substantive conclusion is:

> **PHENOTYPE UNRESOLVED — SPECIFIC MISSING EXPERIMENT OR FLUX IDENTIFIED**
> The missing object is the independently calibrated, stimulus-dependent
> membrane-partitioned cation-current/carbon-support trajectory that supplies
> a valid WT denominator and physical time map, including luminal Na/K closure.

This conclusion continues through the mandated missing-balance and minimal
chassis stages; it does not stop at the failed state-resolved AE4 family.
