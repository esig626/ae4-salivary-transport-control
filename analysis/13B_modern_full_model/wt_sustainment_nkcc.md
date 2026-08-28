# WT sustainment repair: stimulus-dependent NKCC1

## Scope and firewall

This is a WT-only G5 reconstruction attempt on the frozen corrected-water
root `G2_BALANCED_APICAL_K_BIASED_G10:M00`. It investigates whether the
provisional model's scale-independent secretion decline is caused by leaving
NKCC1 capacity stimulus-independent. No AE4-null secretion magnitude, ratio,
or time course was accessed or used. The strict held-out ledger is not
imported by the implementation or runner.

The executable artifacts are:

- `src/modern_full_model/nkcc_stimulation.py`;
- `src/modern_full_model/run_nkcc_stimulation_profile.py`;
- `tests/test_modern_nkcc_stimulation.py`;
- `results/13B_modern_full_model/wt_nkcc_stimulation_profile.csv`;
- `results/13B_modern_full_model/wt_nkcc_stimulation_summary.json`; and
- `results/13B_modern_full_model/wt_nkcc_solver_crosscheck.json`.

## Evidence and observation-map limit

Peña-Münzenmayer et al. 2015 (P15, [DOI
10.1074/jbc.M114.612895](https://doi.org/10.1074/jbc.M114.612895)) provide the
controlling native mouse SMG observations:

| Observation | P15 result | Use here |
|---|---:|---|
| WT initial Cl uptake, 0.3 uM CCh | `2.18 +/- 0.20 x 10^-3 s^-1`, `n=7` | Centers the WT muscarinic-rate sensitivity |
| WT initial Cl uptake, 0.3 uM CCh + 5 uM IPR | `2.02 +/- 0.10 x 10^-3 s^-1`, `n=10` | Independent WT combined-arm scale |
| WT initial Cl uptake, 5 uM IPR | `0.40 +/- 0.07 x 10^-3 s^-1`, `n=4` | Shows that the strong uptake arm is muscarinic, not a reason to route beta input into NKCC1 |
| Isolated NKCC1 signal | 50 uM bumetanide reduced it `95.6 +/- 0.9%`; no detected genotype difference | Hard pathway-direction and genotype-invariance constraint |
| WT gland flow | graphically about `9--10 uL/min` at one-minute samples for 10 min | Scale-independent flow-shape gate |

The SPQ uptake experiment is a low/high-Cl depletion/repletion protocol and
reports a net fluorescence initial slope in a perturbed cell. The full model
does not yet contain the assay's depleted-state distribution or SPQ
fluorescence map. Therefore the reported uptake value is **not** fitted as an
NKCC1 capacity constant. For a transparent scale diagnostic only, the frozen
resting root's existing mechanistic value is

`2 J_NKCC1 / n_Cl = 1.1650653 x 10^-3 s^-1`.

An algebraic `1.75` multiplier would make that root-state diagnostic
`2.0389 x 10^-3 s^-1`, close to the P15 totals, but the biological contexts
are not exchangeable. Multipliers through `2.0` are consequently labelled
`P15_CENTERED_CONTEXT_MISMATCH_SENSITIVITY`, not an uncertainty interval.
Larger values are explicitly labelled `UNBOUNDED_FLOW_SHAPE_PROBE`.

## Smallest nested laws

Let

`u_Ca = clip((Ca - 0.058)/(0.55 - 0.058), 0, 1)`.

This maps the frozen resting calcium and historical-lineage stimulated
calcium sensitivity to zero and one without adding an unsupported Hill
coefficient or half-saturation.

| Law | Equation | New quantities | Exact nesting |
|---|---|---|---|
| `N0_BASELINE` | `m_NKCC=1` | none | Existing G4/G5 core |
| `N1_CCH_ALGEBRAIC` | `m_NKCC=1+(M-1)u_Ca` | fully activated multiplier `M` | `M=1` is N0 for every input; rest is N0 for every `M` |
| `N2_CCH_EFFECTIVE_DYNAMIC` | `dn/dt=(u_Ca-n)/tau_N`; `m_NKCC=1+(M-1)n` | one bounded effective state, `M`, and unmeasured `tau_N` | `M=1` is N0; `tau_N -> 0` approaches N1; rest `n=0` is exact |

The multiplier is applied only to the existing reversible NKCC1 capacity.
The ideal-solution affinity, reversal, `1 Na : 1 K : 2 Cl` conserved source,
and thermodynamic direction are unchanged. The extension cannot inject
charge, carbon, buffer, current, or water directly. It is composed with, and
does not replace, mandatory dynamic beta/cAMP/PKA regulation of AE4.

`N2` is a mechanistic sensitivity, not a measured WNK/SPAK or
phosphorylation state. P15 provides no native acute NKCC1 activation trace or
time constant under the matched gland protocol.

## Iterative reconstruction results

All rows use the same frozen WT root, transport parameters, 600-s physical
protocol, one-minute `9--10 uL/min` graphical-envelope gate, and production
Radau tolerances. No root or resting parameter was refitted.

| Attempt | Defect addressed | Result | Decision |
|---|---|---|---|
| `G5-N0` | Baseline diagnosis | Minute-flow max/min `1.733399`, versus allowed `1.111111`; endpoint Cl `20.0380 mM` | Reject as WT dynamic model; retain exact control |
| `G5-N1` | Add the smallest CCh/NKCC capacity link | `M=1.25--8` monotonically improved shape at first, but even `M=8` gave ratio `1.149869` and failed | Reject: a static capacity multiplier cannot close the sustainment defect |
| `G5-N2-source` | Test one effective activation state in the P15-centered `M<=2` scale | Best was `M=2`, `tau_N=120 s`: ratio `1.275461`; endpoint Na/K/Cl `25.7255/94.5849/31.9556 mM` | Reject: no source-centered member passes WT flow shape |
| `G5-N2-probe` | Determine whether the missing degree of freedom can repair shape numerically | Lowest profiled passing gain was `M=6`, `tau_N=60 s`: ratio `1.099181`, broad ion/volume/pH and conservation gates pass | Numerical repair only; reject for a physiological freeze because gain and time are unmeasured |

For the numerical `M=6`, `tau_N=60 s` probe:

| Quantity | 60 s | 600 s |
|---|---:|---:|
| NKCC1 capacity multiplier | `4.160603` | `5.999773` |
| Single-cell flow | `1.161091 x 10^-3 pL/s` | `1.056324 x 10^-3 pL/s` |
| Mechanistic `2 J_NKCC1/n_Cl` diagnostic | `4.966991 x 10^-3 s^-1` | `4.237590 x 10^-3 s^-1` |

Its 600-s cumulative flow is `0.6710135 pL/cell`. Endpoint cell state is
Na `35.1105 mM`, K `89.5703 mM`, Cl `46.6044 mM`, pH `6.88844`, and volume
`1.24604 pL`. These values remain positive and within the deliberately broad
predeclared dynamic gates, but positivity does not supply the missing native
NKCC1 gain evidence.

The selected-root profile contains 40 laws: 8 unbounded N2 probes pass the
graphical WT flow-shape gate, whereas zero of the 15 P15-centered controls and
sensitivities do. This preserves the distinction between numerical
sufficiency and physiological support.

## Numerical checks

Eight focused tests verify:

1. exact resting and `M=1` nesting;
2. unchanged NKCC1 `1:1:2` stoichiometry;
3. exact bulk-charge cancellation;
4. one and only one added N2 state;
5. finite parameter/domain guards;
6. explicit separation of the mechanistic NKCC1 rate from the SPQ target; and
7. Radau/BDF agreement.

For `N2 M=6, tau_N=60 s`, Radau versus BDF gives maximum relative state
difference `4.67 x 10^-7`, maximum relative flow difference
`1.30 x 10^-6`, and cumulative-flow relative difference `5.48 x 10^-8`.
The best P15-centered member has comparably small solver differences. All
declared conservation residuals pass their unit-specific tolerances.

## Conclusion and decision-critical measurement

On the selected frozen WT root, omission of stimulus-dependent NKCC1 is a
real model defect candidate but **not a source-constrained repair**. An
algebraic capacity change is insufficient. A delayed effective law can flatten
the flow curve only by using a large, unmeasured capacity amplification. It
must not be retained merely because it passes the WT graphical envelope.

The decision-critical experiment is a native mouse SMG-acinar,
bumetanide-sensitive NKCC1 flux time course during the matched `0.3 uM CCh +
5 uM IPR` protocol, with physical-second resolution and simultaneous
intracellular Na, K, Cl, and cell volume. That single experiment would decide
whether the approximately fourfold-at-60-s to sixfold-at-600-s activation
required by the numerical repair is physiologically available or whether a
different cation/homeostasis or epithelial-flow closure is required.

### Cross-root follow-up

The independent cross-root runner subsequently evaluated 144 primary rows
(16 nested NKCC laws on each of all 9 retained WT resting roots) at the fixed
historical calcium anchor. It found zero source-eligible candidates. Its best
nominal-gain row was `N2 M=1.75, tau_N=120 s` on
`G2_ASYMMETRIC_APICAL_K_HEAVY_G10:M00`, with flow max/min `1.304388` and a
failed shape gate. An unrestricted `M=8, tau_N=300 s` probe on that root
passed numerically with ratio `1.027331`, but was ineligible by its declared
out-of-scale status. See
`results/13B_modern_full_model/wt_g5_cross_root_repair_summary.json` and its
hashed profile. Thus the selected-root conclusion is not an accident of that
root: changing among the already-valid resting branches does not produce a
source-constrained NKCC-only repair over the tested panel.
