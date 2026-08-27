# Primary experimental evidence for Task 12

This ledger separates measurements from figure-derived quantities and from the
authors' mechanistic interpretations.  It is the experimental boundary for the
AE4 candidate family; the 2018 modelling paper is not used here to define AE4
biology.

## Evidence classes and modelling use

- **Direct measurement**: an ion-sensitive dye, electrophysiological recording,
  or gland secretion measurement reported by the primary study.  Reported error
  conventions and sample sizes are retained.
- **Reported calculation**: a percentage, rate, or fitted parameter calculated
  by the study authors from direct measurements.
- **Figure-derived approximation**: a ratio read from plotted group means when
  the article does not report that ratio numerically.  These values are useful
  only as coarse envelopes and are not exact calibration constants.
- **Interpretation**: a mechanistic conclusion or working model.  It is not
  promoted to a measurement.

The 2015 AE4-knockout secretion result is **held-out validation**.  It must not
be used to fit AE4 capacity or any non-AE4 chassis parameter.  WT resting
observables and the independent 2016/2021 transport assays may constrain a
candidate before knockout prediction.  The machine-readable companion is
`results/12_ae4_mechanism_reconstruction/phenotype_targets.csv`.

## Primary sources

1. Peña-Münzenmayer G, Catalán MA, Kondo Y, Jaramillo Y, Liu F, Shull GE,
   Melvin JE. 2015. “Ae4 (Slc4a9) Anion Exchanger Drives Cl−
   Uptake-dependent Fluid Secretion by Mouse Submandibular Gland Acinar
   Cells.” *J Biol Chem* 290:10677–10688.
   [DOI 10.1074/jbc.M114.612895](https://doi.org/10.1074/jbc.M114.612895),
   [PMCID PMC4409235](https://pmc.ncbi.nlm.nih.gov/articles/PMC4409235/).
2. Peña-Münzenmayer G, George AT, Shull GE, Melvin JE, Catalán MA. 2016.
   “Ae4 (Slc4a9) Is an Electroneutral Monovalent Cation-dependent
   Cl−/HCO3− Exchanger.” *J Gen Physiol* 147:423–436.
   [DOI 10.1085/jgp.201611571](https://doi.org/10.1085/jgp.201611571),
   [PMCID PMC4845690](https://pmc.ncbi.nlm.nih.gov/articles/PMC4845690/).
3. Peña-Münzenmayer G, Kondo Y, Salinas C, Sarmiento J, Brauchi S, Catalán
   MA. 2021. “Activation of the Ae4 (Slc4a9) Cation-driven Cl−/HCO3−
   Exchanger by the cAMP-dependent Protein Kinase in Salivary Gland Acinar
   Cells.” *Am J Physiol Gastrointest Liver Physiol* 321:G628–G638.
   [DOI 10.1152/ajpgi.00145.2021](https://doi.org/10.1152/ajpgi.00145.2021),
   [PMCID PMC8887885](https://pmc.ncbi.nlm.nih.gov/articles/PMC8887885/).

## 2015: organ phenotype and native acinar-cell constraints

The 2015 paper has the strongest phenotype evidence because it measures native
mouse glands and native acinar cells.  Values are mean ± SE unless noted in the
paper; `n` in Figure 1A/B is glands/mice, whereas Table 1 `n` is the number of
cell experiments used for the reported measurement.

### Secretion phenotype: Figure 1A/B and the corresponding Results subsection

- Ex vivo submandibular glands were perfused with physiological solution and
  stimulated for 10 min with **0.3 µM carbachol plus 5 µM isoproterenol**.
- Total saliva from AE4-null glands was **35 ± 4.7% lower** than littermate WT
  (`n=6` per genotype; Student's *t* test, `p<0.05`; Figure 1A inset).  The
  corresponding derived KO/WT target is **0.65 ± 0.047**.  The percentage is
  reported by the authors; the ratio is the linear complement `1 − 0.35`, so
  its absolute SEM remains `0.047`.  It is arithmetic, not an additional
  measurement or an independently estimated ratio error.
- The first 2–3 min were comparable, followed by a sustained AE4-null deficit
  after 3 min (Figure 1A time course).  Thus the phenotype is not just an
  instantaneous onset effect.
- Acinus-specific AE2-null and its own littermate control were essentially
  identical in flow rate and 10-min total saliva (`n=6` per group; Figure 1B).
  The paper does not report a precise AE2 KO/WT ratio, so “1.0” must not be
  treated as an exact measured number.
- AE4 and AE2 lines had different genetic backgrounds and different control
  time-course shapes.  The authors therefore compared each knockout only with
  its corresponding littermate control.

### Resting chloride, pH, and stimulated chloride: Table 1 and Figure 1C/D

| Observable | AE4 WT | AE4 KO | AE2 control | AE2 KO | Direct implication |
|---|---:|---:|---:|---:|---|
| Resting [Cl−]i, mM | 50.10 ± 1.50 (`n=7`) | 36.50 ± 1.60 (`n=6`, `p<0.001`) | 53.40 ± 1.80 (`n=6`) | 54.50 ± 1.80 (`n=6`) | AE4 loss lowers resting chloride; AE2 loss does not |
| Resting pHi | 6.91 ± 0.07 (`n=4`) | 6.89 ± 0.02 (`n=4`) | 6.87 ± 0.01 (`n=4`) | 6.95 ± 0.05 (`n=6`) | no detected resting-pH defect |
| Initial Cl− uptake, CCh+IPR, 10−3 s−1 | 2.02 ± 0.10 (`n=10`) | 0.90 ± 0.09 (`n=9`, `p<0.001`) | 2.11 ± 0.50 (`n=13`) | 2.39 ± 0.20 (`n=18`) | AE4 KO reduces combined-stimulus reuptake; AE2 KO does not |
| Initial Cl− uptake, CCh only, 10−3 s−1 | 2.18 ± 0.20 (`n=7`) | 2.30 ± 0.10 (`n=6`) | 2.20 ± 0.20 (`n=8`) | 2.30 ± 0.20 (`n=9`) | no AE4-specific defect under muscarinic stimulus alone |
| Initial Cl− uptake, IPR only, 10−3 s−1 | 0.40 ± 0.07 (`n=4`) | 0.20 ± 0.03 (`n=7`, `p<0.05`) | 0.70 ± 0.10 (`n=7`) | 0.58 ± 0.05 (`n=6`) | AE4-specific defect under beta-adrenergic stimulus |

The authors also report resting [Cl−]i as **27 ± 3.1% lower** and the combined
stimulus uptake rate as **55.4 ± 4.5% lower** in AE4-null cells.  Those are
reported calculations from the native-cell measurements, not independent
observations.  The SPQ initial chloride-exit signal was not detectably changed,
which argues against gross compensation through the apical TMEM16A response
(Figure 1C/D and Table 1), but does not measure channel abundance.

### Beta-adrenergic specificity and compensation tests

- In the Results text following Figures 3–5, the chloride-depletion/reuptake
  assay increased from **1.56 ± 0.19** (`n=10`) to **4.33 ± 0.91** (`n=13`)
  10−3 s−1 in control acini without versus with 5 µM isoproterenol
  (`p=0.016`).  AE4-null acini changed only from **0.95 ± 0.10** (`n=8`) to
  **1.08 ± 0.09** (`n=9`; `p=0.366`).  These direct native-cell data establish
  an AE4-associated beta-adrenergic response, not its molecular kinase step.
- Figure 2A/B isolates NKCC1-dependent uptake in HCO3−-free solution with
  30 µM ethoxyzolamide and 50 µM T16Ainh-A01.  NKCC1 functional activity was
  not detectably different in AE4- or AE2-null cells versus their controls;
  bumetanide reduced the response by **95.6 ± 0.9%**.  This rules out a large
  functional NKCC1 gain/loss under that assay, not every possible in vivo
  compensation.
- Stimulation-induced NHE-mediated alkalinization was not detectably altered
  in either knockout (Figure 2C and Results compensation paragraph).  Resting
  pH was also unchanged (Table 1).  The paper explicitly notes that these pH
  experiments do **not** directly measure carbonic-anhydrase activity.
- HCO3− secreted during the 10-min combined stimulus was **5.3 ± 1.4** versus
  **7.1 ± 0.8 µequiv/10 min** for AE4 WT and KO (`n=7` each, `p=0.28`).  AE2
  control and KO were **4.0 ± 0.4** (`n=5`) and **6.3 ± 1.8** (`n=6`), also
  `p=0.28`.  A model should not claim that reduced total bicarbonate secretion
  is the demonstrated cause of the AE4 flow phenotype.

### 2015 evidence boundary

The 2015 study directly establishes the organ phenotype, native chloride and
pH phenotypes, and beta-adrenergic dependence.  It did **not** establish AE4
ion stoichiometry, prove cation transport, identify PKA, or directly localize
AE4 protein in acinar membrane: antibody staining was not specific, and the
authors inferred a basolateral role from function while acknowledging that an
apical location could not be excluded.

## 2016: transported ions, reversibility, electroneutrality, and candidate
stoichiometries

The 2016 paper combines a native AE2-null acinar preparation (to isolate native
AE4) with heterologous mouse/human AE4 in CHO-K1 cells.  Results are mean ± SEM.
The CHO experiments establish transport capabilities and thermodynamic
constraints, but their absolute rates are not gland-scale capacities.

| Experimental fact | Exact source and conditions | `n` and uncertainty | Classification and constraint |
|---|---|---|---|
| Native AE4 is Na-dependent and HCO3-dependent | Figure 1A/B: BCECF in AE2-null acini; external Cl− 128.3→4 mM; Na replaced by NMDG; AE2/AE4 double KO control | Na condition `n=8`, Na-free `n=7`, double KO `n=8`; ≥5 acinar cells/experiment from ≥3 preparations; mean ± SEM; ANOVA/Bonferroni `p<0.05` | direct native-cell measurement; Na dependence alone does not distinguish transport from gating |
| Recombinant AE4 couples Cl− loss to HCO3− uptake and external Na | Figure 2A–D; mouse and human AE4 in CHO-K1; low external Cl, Na/NMDG and HCO3-free controls | figure-level mean ± SEM from ≥3 electroporations | direct dye measurements; establishes coupled dependencies |
| Na moves with HCO3 and opposite Cl | Figure 3: SBFI; low external Cl increases [Na]i in AE4-expressing cells and restoring external Cl reverses the response | AE4 `n=15`, nontransfected `n=13`; mean ± SEM | direct Na signal plus direction reversal; recovery could also include Na/K-ATPase, as the authors state |
| Transport is voltage independent within the tested range | Figure 4A–F: simultaneous SPQ/current clamp and BCECF/voltage clamp; ΔVm <2 mV; rates at −100 and 0 mV essentially identical | Cl-exit summary AE4 `n=7`, nontransfected `n=8`; voltage-clamp `n=5` per voltage and group; ≥3 electroporations | direct electrophysiology supports electroneutrality, not a unique stoichiometry |
| No AE4-specific macroscopic current | Figure 4 Results | at −100 mV: −1.68 ± 0.44 vs −1.27 ± 0.46 pA/pF (`n=5` each, `p=0.54`); at 0 mV: −0.28 ± 0.09 vs −0.21 ± 0.14 (`n=5` each, `p=0.69`) | direct negative result consistent with electroneutral transport |
| Cl-free “NBC-like” behavior is exchange | Figures 5–6 and Table 4: whole-cell control of intracellular Cl/HCO3; HCO3-dependent alkalinization without AE4-specific current | Table 4 gives mean ± SEM, `n=4–7` depending pipette solution | direct assay plus interpretation: HCO3− can substitute internally; not evidence for physiological net NBC cotransport |
| K and other monovalent cations support exchange | Figure 8A/B: Na, K, Li, Rb, Cs versus NMDG | Na `n=10`, K `n=8`, NMDG `n=6`; Li/Rb/Cs `n=6` each; mean ± SEM | direct cation-substitution measurements; AE4 is poorly selective among tested monovalent cations |
| K is transported, not only an external activator | Figure 8C/D: PBFI under Na-free conditions; inward Cl and outward K gradients | AE4 `n=6`, nontransfected `n=8`; mean ± SEM | direct K signal strongly supports K movement through the AE4-dependent process |
| Cation dose responses are cooperative | Figure 9 and Methods Eq. 4; external Na or K 5, 25, 50, 100, 125, 150 mM | ≥4 experiments/condition, ≥3 electroporations; plotted mean ± SEM; no fit-parameter errors reported | reported Hill fits constrain assay response shape, not transported ion count |

The Figure 9 fitted values are:

| Cation | EC50, mM | Hill `nH` | `Rmin` | `Rmax` |
|---|---:|---:|---:|---:|
| Na+ | 49 | 2.0 | 0.3 | 1.5 |
| K+ | 62 | 1.8 | 0.4 | 1.6 |

`Rmin` and `Rmax` are alkalinization rates in the paper's assay convention
(10−3 s−1).  No SEM or confidence interval is reported for EC50, Hill, Rmin, or
Rmax.  The similar Rmax values justify testing a near-equal Na/K capacity
partition, but do not prove that the Rmax ratio transfers unchanged to native
mixed-cation transport.

### Stoichiometry: what is and is not established

- Electroneutrality plus the measured transport directions permit, among other
  possibilities, the three integer candidates explicitly evaluated in Figure 7
  and Methods Eqs. 1–3: **Cl:Na:HCO3 = 1:1:2, 1:2:3, and 2:1:3**.
- These are author-considered thermodynamic candidates, **not measured
  stoichiometries**.  Figure 7 shows that, under the authors' Na-only reference
  concentrations, each requires nonphysiological conditions for chloride
  influx: Cl_i below about 7, 0.5, or 26 mM; Na_i above about 104, 143, or
  54 mM; or HCO3_i above about 49, 84, or 29 mM for 1:1:2, 1:2:3, or 2:1:3,
  respectively.
- Figure 10 then adopts **1:1:2 as a working model** after assuming AE4 does not
  distinguish Na from K and that total cation concentrations approximately
  balance across the membrane.  Under that assumption, Cl_i <73 mM or HCO3_i
  >16 mM favors chloride uptake.  Figure 10 does not experimentally select
  1:1:2 over 1:2:3 or 2:1:3.
- Hill values 2.0 and 1.8 do not determine cation stoichiometry.  The Discussion
  explicitly retains two explanations: more than one transported cation per
  cycle, or a transported cation plus additional cation-dependent allosteric
  regulation.
- Na- and K-sensitive fluorescence strongly supports transport of both cations;
  therefore a purely nontransported-cation gate is not sufficient as the sole
  mechanism.  A regulatory cation site in addition to transported cation remains
  viable.
- The transporter is reversible in the imposed-gradient assays.  Any candidate
  used as a physical mechanism must therefore have a thermodynamic reversal
  condition and cannot be a permanently inward chloride source.

## 2021: beta-adrenergic/cAMP/PKA regulation

The 2021 study supplies the strongest regulation evidence, but not a secretion
phenotype.  Figure 1 used 17 female mice aged 8–12 weeks: 10 acinus-specific
AE2-null mice and seven controls.  Both glands were processed, with at least
three animals per condition.  Native-cell Figure 1 captions label `n` as
experiments; CHO-cell figures label `n` as cells.  Error bars are **SD**, unlike
the SEM convention in the 2015/2016 papers.

| Result | Exact source, conditions, and `n` | Evidence strength |
|---|---|---|
| Beta-adrenergic stimulation increases total native exchange | Figure 1A/B: SPQ chloride depletion/reuptake in WT acini, 50 µM bumetanide + 50 µM T16Ainh-A01; control `n=11`, 5 µM IPR `n=12`; unpaired two-tailed t test, caption reports `p<0.001` | direct native measurement, but WT signal contains AE2 and AE4 |
| IPR activates native AE4 and H89 prevents it | Figure 1C in acinus-specific AE2-null acini: control `n=12`, 5 µM IPR `n=11`, IPR after 10 µM H89 `n=7`; Kruskal–Wallis, `p<0.0001` and `p<0.05` | strongest AE4-specific native regulation evidence; H89 is explicitly nonselective, so pharmacology alone is not kinase-site proof |
| Forskolin increases recombinant AE4 activity and H89 prevents it | Figure 2B/C: BCECF alkalinization in AE4-expressing CHO cells; control `n=40`, 10 µM forskolin `n=39`, forskolin + 10 µM H89 `n=41`; ≥3 electroporations; Kruskal–Wallis, `p<0.0001`/`p<0.01` | direct recombinant measurement linking adenylate cyclase/cAMP to macroscopic AE4 activity |
| Constitutively active PKAc increases both HCO3-linked and Cl flux | Figure 3: Ae4+PKAc versus Ae4+K73M dominant-negative PKAc; BCECF `n=18` vs `16`, SPQ `n=8` vs `9`; ≥3 electroporations; Mann–Whitney, `p<0.005` | two independent ion readouts support a PKA-dependent activity increase |
| S173 is required for the PKAc response in this construct | Figures 5–6: S173A and S273A reach membrane and remain functional; S173A+PKAc DN/PKAc `n=14/15`, S273A `n=12/12`; ≥3 electroporations; S173A no significant activation, S273A `p<0.01` | mutational evidence for S173 dependence; no direct phosphosite measurement, so “PKA phosphorylates S173” remains an interpretation |

### Fold envelope and its limitations

The paper reports group measurements and significance but does not publish a
single transferable sustained-capacity multiplier.  Approximate ratios read
from the plotted group means are **~1.6-fold** for WT total exchanger activity
(Figure 1B), **~3-fold** for native AE4-only activity (Figure 1C), and
**~1.25-fold** for recombinant forskolin stimulation (Figure 2C).  These are
figure-derived, rounded sensitivity values; their underlying measurements are
direct, but the ratios have no reported SEM or digitization uncertainty.

Consequently `(1.25, 1.6, 3.0)` is a defensible **response envelope**, not three
equally strong AE4 capacity calibrations:

- `1.6` is not AE4-specific because WT acini express AE2 and AE4;
- `3.0` is AE4-specific and native but is an initial exchanger-rate ratio, not a
  whole-gland sustained turnover multiplier;
- `1.25` is recombinant and may be host- or expression-dependent;
- none establishes activation timing relative to the 10-min combined
  carbachol/isoproterenol secretion protocol;
- the macroscopic assays cannot distinguish increased turnover, more active
  transporters at the plasma membrane, or a PKA-induced stoichiometric change.

The 2021 Discussion explicitly identifies those three unresolved mechanisms.
It also cites the 2015 observation that muscarinic stimulation alone does not
activate native AE4.  Calcium therefore cannot be silently substituted for the
beta/cAMP/PKA input in a candidate implementation.

## Cross-study constraints and nonclaims

The source-supported minimum candidate family must preserve all of the
following:

1. reversible chloride/cation-bicarbonate exchange;
2. electroneutrality;
3. transported Na and transported K branches or an equivalent mixed-cation law;
4. measured cooperative cation dependence as an occupancy/kinetic possibility,
   without equating Hill slope to stoichiometry;
5. all three explicitly considered electroneutral stoichiometries until an
   independent gate excludes them;
6. an independent beta/cAMP/PKA activation input in stimulated comparisons;
7. WT resting chloride and pH consistency before the held-out knockout test;
8. negligible AE2-knockout secretion effect and the native AE4-knockout
   chloride/pH phenotypes as independent validation.

The primary literature does **not** measure a unique AE4 stoichiometry, a native
Na:K flux partition, a gland-scale AE4 density, a sustained PKA capacity step,
or the stimulated intracellular Na/K/pH trajectories needed to distinguish all
candidate mechanisms.  Failure to reproduce the held-out flow phenotype after
enforcing these constraints is therefore evidence against AE4-law replacement
on the fixed chassis, not evidence that the primary measurements are wrong.
