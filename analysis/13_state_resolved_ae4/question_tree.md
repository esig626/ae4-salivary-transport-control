# Task 13 exact question tree

**Frozen with the evidence table before fitting, 2026-08-27 UTC.** This is a
decision/decomposition tree, not a model design. A candidate may not jump from
transporter adequacy directly to the AE4-null secretion target.

## Holdout firewall

| Stage | Permitted evidence | Prohibited use |
|---|---|---|
| Transporter calibration | 2016 Na/K transport, electroneutrality, dose response; 2021 regulation; 2025 mutation hierarchy; assay contexts in `evidence_freeze.md` | Every 2015 AE4-null phenotype and every whole-gland secretion target |
| WT whole-cell construction | WT/rest and WT-stimulation evidence, source-supported topology, physical closure | AE4-null Cl, pH, volume, secretion magnitude, or secretion time course |
| Stage A reveal | Compare frozen predictions with all AE4-null evidence | Any retuning after seeing AE4-null values and calling the result out of sample |
| Stage B localization, only if all source-supported state-resolved families fail on the fixed chassis | AE4-null Cl and pH (and volume only if a direct value existed; none was found) as missing-balance localization constraints | Total/integrated AE4-null secretion, early/late secretion pattern, or any pointwise AE4-null flow value |
| Final validation | Frozen model prediction versus exact 10-min ratio and digitized/caption-level temporal envelope | Optimizing any parameter, bound, scale, or module choice to these secretion outcomes |

## Ordered gates

### Gate 1: transporter-level adequacy

**Question.** Does the candidate reproduce the transporter facts in their own
assay systems without using gland or knockout outcomes?

Required checks:

1. Na is transported with HCO3 opposite Cl, and K is also transported
   (E16-03, E16-05), rather than acting only as an allosteric ligand.
2. Net transport is consistent with the `<2 mV`, current, and -100/0 mV
   observations (E16-04).
3. Effective Na and K dependence is compatible with `EC50=49/62 mM` and Hill
   `2.0/1.8`, with missing fit uncertainty acknowledged (E16-07).
4. The candidate does not erase the 2025 Na/K functional nonidentity,
   especially T756A-T448I: no detectable Na-dependent activity but substantial
   K-dependent activity (E25-04).
5. If regulation is claimed, it is consistent with the approximately 25%
   forskolin effect, native IPR/H89 result, and S173 requirement, without
   inventing a measured phosphorylation time constant (E21-01-E21-06).
6. Reversal/direction and electroneutrality are demonstrated under each source
   bath; a fit to alkalinization alone is not sufficient.

**Fail branch.** Reject or revise the transporter family before whole-cell
integration. A failed family says nothing about whole-cell chassis adequacy.

**Pass branch.** Freeze the transporter parameter combinations and documented
nonidentifiabilities, then proceed to Gate 2.

### Gate 2: WT whole-cell adequacy

**Question.** With the transporter frozen, can a source-supported WT chassis
reach physical rest and reproduce independent WT stimulation behavior?

Required checks:

- positive concentrations and volume, electroneutrality, membrane-current
  closure, mass balance, and correct thermodynamic signs;
- WT resting Cl (`50.10 +/- 1.50 mM`) and pHi (`6.91 +/- 0.07`) in the exact
  Ae4-line context, plus independently justified Na, K, and volume ranges;
- the 2015 WT CCh+IPR protocol and flow scale without reading the paired null
  curve;
- topology is independently supported; abstract-level apical pump/K evidence
  cannot be converted into measured 30%/40% membrane fractions;
- all resting branches in the declared domain are retained and reported.

**Fail branch.** Localize the WT closure/residual first. Do not inspect or tune
against any AE4-null outcome. A transporter that passes Gate 1 cannot rescue an
invalid WT chassis by construction.

**Pass branch.** Freeze WT roots, branch-selection rule, protocol, scales, and
all non-AE4 parameters; proceed to Gate 3.

### Gate 3: AE2-null adequacy

**Question.** Does the frozen WT/chassis representation preserve the distinct
AE2-null phenotype without using AE4-null data?

Required checks:

- negligible change in 10-min secretion relative to the Ae2 line's own
  littermate control (E15-04), reported as a qualitative/no-detected-difference
  constraint because no equivalence margin is published;
- no detected changes in resting Cl, pHi, or combined-stimulus Cl uptake in
  that line (E15-11);
- no cross-line comparison of early WT kinetics between C57BL/6 Ae2 and
  BS/129Svj Ae4 backgrounds.

**Fail branch.** Localize AE2/topology/background assumptions. Do not compensate
by tuning toward the AE4-null secretion result.

**Pass branch.** Freeze the AE2 comparison and proceed to Gate 4.

### Gate 4: AE4-null ionic adequacy

**Question.** After the Stage A reveal, does deletion of AE4 predict the native
ionic phenotype without any AE4-null calibration?

Stage A checks include resting Cl (`36.50 +/- 1.60 mM`, 27 +/- 3.1% lower),
resting pHi (`6.89 +/- 0.02`, no detected change), combined CCh+IPR Cl-uptake
rate (`0.90 +/- 0.09 x 10^-3 s^-1`, 55.4 +/- 4.5% lower), IPR-specific versus
CCh-only behavior, and the assay-specific no-detected-change NKCC1 result.

**Pass branch.** Freeze the ionic prediction unchanged and proceed to Gate 5.

**Fail branch for one transporter family.** Record the ionic residual vector;
continue with every other source-supported state-resolved family.

**Fail branch for all source-supported families on the fixed chassis.** Enter
Stage B. Release AE4-null Cl/pH evidence only for missing-balance localization.
There is no direct volume datum to release. Reconstruct the smallest
independently supported chassis addition whose source direction addresses the
ionic residual, then repeat Gates 2-4. The secretion targets remain sealed.

### Gate 5: AE4-null secretion magnitude

**Question.** Does the fully frozen candidate predict the independent whole-
gland integral?

The only exact target is systemic Ae4-null total saliva over 10 min under 0.3
uM CCh + 5 uM IPR: **35 +/- 4.7% lower than littermate WT**, equivalently an
algebraic KO/WT complement of about `0.65 +/- 0.047` (E15-02). “30%” is only a
rounded later citation; 24% is a historical model output.

**Fail branch.** Preserve the magnitude residual and proceed to Gate 6 before
assigning mechanism. Do not optimize to the ratio. In Stage B, project the
remaining balance discrepancy onto independently supported modules without
opening the secretion target.

**Pass branch.** Magnitude alone is not success; proceed to Gate 6.

### Gate 6: time-resolved stimulation adequacy

**Question.** Does the same frozen prediction explain the temporal phenotype,
not only its integral?

Required shape:

- Ae4-null and WT flow comparable for the initial 2-3 min;
- divergence after about 3 min;
- sustained late deficit, with late null flow less than half WT;
- consistency with the uncertainty-labeled minute envelope in
  `evidence_freeze.md`, without pointwise fitting to digitized pixels.

**Magnitude pass / time fail.** The model has matched an integral but has not
explained delayed propagation; localize missing regulation, pool dynamics, or
slow homeostasis without retuning to the null curve.

**Magnitude fail / time pass.** The propagation mechanism is plausible but its
network gain is wrong; report this separately from a timing failure.

**Both pass.** Proceed to adversarial leakage, simplification, solver,
branch/continuation, and independent-reproduction audits. Only after those
audits may the mechanism be called successful.

## Interpretation matrix

| Last passed gate | Substantive conclusion allowed | Conclusion not allowed |
|---|---|---|
| 1 only | A transporter family is consistent with isolated/native transporter constraints. | It explains salivary secretion. |
| 2 only | The family can be embedded in a physically acceptable WT chassis. | Knockout propagation is correct. |
| 3 only | WT and the distinct AE2-null comparison are jointly adequate. | AE4-null ionic or flow adequacy. |
| 4 | AE4 deletion propagates into the observed ionic direction/magnitude without Stage A leakage (or after an explicitly labeled Stage B localization). | The secretion phenotype is explained. |
| 5 | The held-out 10-min magnitude is predicted. | Delayed divergence is explained. |
| 6 plus audits | Molecular transport, whole-cell propagation, magnitude, and timing are jointly supported within stated uncertainty and nonidentifiability. | Unique microscopic mechanism, unless alternatives have been falsified. |

The investigation must continue after a failed state-resolved family. It stops
only after a surviving candidate passes the ordered gates and audits, or after
the Stage B residual/localization and minimal chassis reconstruction establish
one of Task 13's required substantive conclusions.
