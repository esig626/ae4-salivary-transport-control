# Task 13B primary-evidence and protocol ledger

**Status:** evidence roles frozen before Task 13B parameter optimization.  
**Scope:** primary experimental evidence, assay context, and historical-model
lineage for the modern full salivary acinar-cell reconstruction.  
**Firewall:** the AE4-null stimulated-saliva magnitude and time course occur
only in `results/13B_modern_full_model/heldout_targets.csv`. They are not
calibration or model-selection data. AE4-null resting chloride and pH are
conditional localization evidence, not WT calibration targets.

This ledger records what was measured, where it was measured, and what use is
licensed. It does not treat a molecular cartoon, an MD trajectory, a published
model output, or a failure to reject a null hypothesis as a direct biological
measurement.

## Evidence and use vocabulary

| Label | Meaning |
|---|---|
| `MEASUREMENT` | Direct assay readout or a quantity calibrated by the source's stated assay. |
| `REPORTED_FIT` | Parameter or summary obtained by the source from a regression, ratio, or calculation. |
| `INFERENCE` | Interpretation supported by measurements but not itself measured. |
| `MODEL_HYPOTHESIS` | Proposed stoichiometry, state order, molecular interaction, MD result, or mechanistic cartoon. |
| `HISTORICAL_ASSUMPTION` | Equation, target, parameter, or output of the 2018 mathematical lineage. |
| `NOT_MEASURED` | Quantity required by a candidate model but absent from the source. |
| `CALIBRATION` | May constrain transporter, regulatory, topology, WT resting, or WT stimulated parameters in the matching observation model. |
| `VALIDATION` | Must not be in the fitting objective; used to challenge a frozen module or model. |
| `CONDITIONAL_LOCALIZATION` | May be inspected only after a WT-valid model has failed the held-out prediction; may not be used to tune the secretion endpoint. |
| `STRICT_HELDOUT` | Never used for fitting, bounds, branch choice, time-map choice, regulatory kinetics, or model-generation selection. |

`n` is retained in the source's unit. Animals, glands, experimental
preparations, and cells are not interchangeable. Error conventions are also
retained: the 2015 and 2016 papers report mean +/- SEM, the 2021 paper reports
mean +/- SD, and the 2025 paper reports mean +/- SEM.

## Audited primary sources and access

| ID | Primary source | Direct full text used | Role and access limitation |
|---|---|---|---|
| `P15` | Peña-Münzenmayer et al., 2015, *JBC*, [DOI 10.1074/jbc.M114.612895](https://doi.org/10.1074/jbc.M114.612895), [PMC4409235](https://pmc.ncbi.nlm.nih.gov/articles/PMC4409235/) | Yes | Controlling native mouse SMG source for secretion, resting Cl/pH, AE2 comparison, NKCC1/NHE checks, and agonist-separated Cl handling. |
| `P16` | Peña-Münzenmayer et al., 2016, *JGP*, [DOI 10.1085/jgp.201611571](https://doi.org/10.1085/jgp.201611571), [publisher full text](https://rupress.org/jgp/article/147/5/423/43529/Ae4-Slc4a9-is-an-electroneutral-monovalent-cation), [PMC4845690](https://pmc.ncbi.nlm.nih.gov/articles/PMC4845690/) | Yes | Native and heterologous transporter mechanism. Pure-cation CHO assays are not native salivary concentrations. |
| `P21` | Peña-Münzenmayer et al., 2021, *AJP Gastrointestinal and Liver Physiology*, [DOI 10.1152/ajpgi.00145.2021](https://doi.org/10.1152/ajpgi.00145.2021), [PMC8887885](https://pmc.ncbi.nlm.nih.gov/articles/PMC8887885/) | Yes | Mandatory beta-adrenergic/cAMP/PKA regulation evidence. It supplies activity contrasts but no regulatory kinetic constants. |
| `P25` | Catalán et al., 2025, *AJP Cell Physiology*, [DOI 10.1152/ajpcell.00346.2024](https://doi.org/10.1152/ajpcell.00346.2024), [PMC12463136](https://pmc.ncbi.nlm.nih.gov/articles/PMC12463136/) | Yes | HEK-293 mutagenesis plus homology/MD hypotheses. No native acinar assay and no direct AE4 structure. |
| `A18` | Almássy et al., 2018, *Pflügers Archiv*, [DOI 10.1007/s00424-018-2109-0](https://doi.org/10.1007/s00424-018-2109-0), [publisher record](https://link.springer.com/article/10.1007/s00424-018-2109-0) | Abstract and publisher supplements | The version-of-record experimental Methods were subscription restricted. Only abstract-level experimental claims and explicit supplement-level model choices are used. Exact animal, bath, temperature, and `n` values are therefore `NOT_MEASURED` here. |
| `K15` | Kondo et al., 2015, *Journal of Dental Research*, [DOI 10.1177/0022034515570943](https://doi.org/10.1177/0022034515570943), [PMC4502782](https://pmc.ncbi.nlm.nih.gov/articles/PMC4502782/) | Yes | Mouse major-gland morphometry and physiology. The SMG acinar fraction is a geometry constraint, not a flow-fit parameter. |
| `K19` | Kondo et al., 2019, *AJP Cell Physiology*, [DOI 10.1152/ajpcell.00421.2018](https://doi.org/10.1152/ajpcell.00421.2018), [PMC6580159](https://pmc.ncbi.nlm.nih.gov/articles/PMC6580159/) | Yes | Same-laboratory mouse SMG/SLS ex vivo protocol record; confirms gland-specific ipsilateral collection. It is protocol context, not a replacement WT-flow target. |
| `C15V` | Catalán et al., 2015, *PNAS*, [DOI 10.1073/pnas.1415739112](https://doi.org/10.1073/pnas.1415739112), [PMC4343136](https://pmc.ncbi.nlm.nih.gov/articles/PMC4343136/) | Yes | Adult acinar-specific Tmem16A deletion, IPR-only secretion, cell-volume, inhibitor, and volume-regulated-anion-current evidence. |
| `PB94` | Poulsen and Bundgaard, 1994, *Pflügers Archiv*, [DOI 10.1007/BF00374318](https://doi.org/10.1007/BF00374318) | Abstract and quantitative record | Rat-parotid luminal and basolateral membrane area per cell volume. Used only for a declared parotid-to-mouse-SMG pump-routing sensitivity. |
| `G07/P10` | Gin et al., 2007, [DOI 10.1016/j.jtbi.2007.04.021](https://doi.org/10.1016/j.jtbi.2007.04.021); Palk et al., 2010, [DOI 10.1016/j.jtbi.2010.06.027](https://doi.org/10.1016/j.jtbi.2010.06.027) | Yes | Published parotid-model lineage for source-mapped whole-cell channel maxima and the Palk calcium gate. The cited Arreola and Thompson/Begenisich experiments retain their rat-parotid assay context. |
| `R10` | Romanenko et al., 2010, *JBC*, [DOI 10.1074/jbc.M109.068544](https://doi.org/10.1074/jbc.M109.068544), [PMC2857126](https://pmc.ncbi.nlm.nih.gov/articles/PMC2857126/) | Yes | Native mouse-SMG CaCC biophysics; its protocol- and voltage-dependent calcium sensitivity is diagnostic, not silently interchangeable with the Palk gate. |
| `V18` | Vera-Sigüenza et al., 2018, *Bulletin of Mathematical Biology*, [DOI 10.1007/s11538-017-0370-6](https://doi.org/10.1007/s11538-017-0370-6), [PMC5792321](https://pmc.ncbi.nlm.nih.gov/articles/PMC5792321/) | Yes | Historical equations, topology, calibration record, and model outputs only; not modern experimental ground truth. |

## Assay-context matrix

| Source | Native or heterologous context | Preparation and genotype | Bath, temperature, perturbation, and observation |
|---|---|---|---|
| `P15` | Native ex vivo whole gland and isolated native acini | Female and male mice, 2--4 months. Systemic `Ae4-/-` on BS/129Svj; acinus-specific `Ae2fl/fl-AQP5/ACID-Cre` on C57BL/6; each compared only with its own littermates. | SMG perfused via common carotid at 37 C with HCO3-containing high-Cl solution; 0.3 uM carbachol (CCh) + 5 uM isoproterenol (IPR), 10 min; flow read every minute. Acinar high-Cl bath: 4.3 mM KCl, 120 mM NaCl, 25 mM NaHCO3, 5 mM glucose, 10 mM HEPES, 1 mM CaCl2, 1 mM MgCl2, pH 7.4; low-Cl gluconate substitutions; 37 C; 95% O2/5% CO2 where HCO3 was present. SPQ and BCECF protocols used inhibitors stated below. |
| `P16` | Native SMG acini plus heterologous CHO-K1 | Native `Ae2-/-` and `Ae2-/-;Ae4-/-` mice, female and male, 2--4 months. CHO-K1 expressing mouse Ae4 variant 3 or human AE4 variant 2, assayed 18--20 h after electroporation. | Imaging at 37 C; simultaneous patch/imaging at 22 C. Core imaging solutions were pH 7.4 with 25 mM HCO3 and 128.3 or 4 mM Cl; Na protocols used 145 mM Na, K protocols 145 mM K; explicit NMDG/HCO3-free controls. Dose series: 5, 25, 50, 100, 125, 150 mM external Na or K. No secretagogue. |
| `P21` | Native SMG acini plus heterologous CHO-K1 | Seventeen female mice, 8--12 weeks: ten acinus-specific `Ae2-/-` and seven control animals. CHO-K1 expressing mouse Ae4 WT, S173A, or S273A with PKA constructs, 18--20 h after electroporation. | 37 C, 95% O2/5% CO2. High/low-Cl solutions matched the P15 composition. Native Cl-reuptake assay contained 50 uM bumetanide + 50 uM T16Ainh-A01; 5 uM IPR; 10 uM H89. CHO conditions used 10 uM forskolin, 10 uM H89, constitutively active PKAc, or K73M dominant-negative PKAc. The plotted 0--600 s traces show the Cl-depletion/reuptake assay, not beta-to-PKA activation onset. |
| `P25` | Heterologous only | Human AE4 isoform 2, WT and mutants, HEK-293 cells, 18--24 h post-transfection. | BCECF at room temperature; pH 7.4; 25 mM HCO3; high/low Cl 128.3/4 mM; Na- or K-dominant baths; NMDG control; 95% O2/5% CO2. Reported `n` is cells. Responding transfected cells were selected as ROIs, which limits population-level interpretation. At least three cells/experiment and three independent transfections; at least five transfections for the double-mutant figure. |
| `A18` | Native parotid acinar experiments plus a mathematical model | Parotid acinar cells; exact species and preparation details are not exposed in the accessible primary abstract. | Local Ca photolysis near the apical membrane with simultaneous electrophysiology/Ca imaging; pump immunolocalization. Exact bath, temperature, `n`, current magnitude, and localization quantification are unavailable from the accessible record. |
| `V18` | Published mathematical model | A generic salivary acinar cell assembled from mixed-gland and mixed-species literature targets. | A prescribed qualitative Ca input stimulates minutes 6--12; ten ODEs plus two algebraic QSS voltage constraints. These are model definitions, not a reconstruction of the 2015 CCh+IPR protocol. |

## `P15`: native gland and acinar evidence

All numerical uncertainties below are SEM. Table 1 `n` values are
experiment-level counts; the Methods state that cell experiments used separate
preparations from at least three mice.

| ID | Direct result | Value and source anchor | Evidence class | Task 13B role |
|---|---|---|---|---|
| `P15-01` | Ae4-line WT resting intracellular Cl | 50.10 +/- 1.50 mM, `n=7`; Table 1, calibrated SPQ | `MEASUREMENT` | WT resting `CALIBRATION` |
| `P15-02` | AE4-null resting intracellular Cl | 36.50 +/- 1.60 mM, `n=6`; 27 +/- 3.1% below matched WT, `p<0.001`; Table 1 | `MEASUREMENT` + `REPORTED_FIT` | `CONDITIONAL_LOCALIZATION`, never secretion tuning |
| `P15-03` | Ae4-line WT and null resting pHi | WT 6.91 +/- 0.07 (`n=4`); null 6.89 +/- 0.02 (`n=4`); Table 1; no detected difference | `MEASUREMENT` | WT `CALIBRATION`; null `CONDITIONAL_LOCALIZATION` |
| `P15-04` | WT combined-stimulus initial Cl uptake | 2.02 +/- 0.10 x10^-3 s^-1, `n=10`; 0.3 uM CCh + 5 uM IPR; Table 1 | `MEASUREMENT` | WT dynamic/transport `CALIBRATION` |
| `P15-05` | AE4-null combined-stimulus initial Cl uptake | 0.90 +/- 0.09 x10^-3 s^-1, `n=9`; reported 55.4 +/- 4.5% lower; Table 1 | `MEASUREMENT` + `REPORTED_FIT` | Genotype `VALIDATION`, not secretion calibration |
| `P15-06` | CCh-only initial Cl uptake | WT 2.18 +/- 0.20 (`n=7`); AE4-null 2.30 +/- 0.10 (`n=6`) x10^-3 s^-1; Table 1; no detected difference | `MEASUREMENT` | WT `CALIBRATION`; null `VALIDATION` |
| `P15-07` | IPR-only initial Cl uptake | WT 0.40 +/- 0.07 (`n=4`); AE4-null 0.20 +/- 0.03 (`n=7`) x10^-3 s^-1; `p<0.05`; Table 1 | `MEASUREMENT` | WT beta-arm `CALIBRATION`; null `VALIDATION` |
| `P15-08` | Initial CCh+IPR Cl exit | WT -0.21 +/- 0.01 (`n=10`); AE4-null -0.18 +/- 0.02 (`n=9`) normalized SPQ change; Table 1; no detected difference | `MEASUREMENT` | WT `CALIBRATION`; null negative-control `VALIDATION` |
| `P15-09` | AE2-null secretion | Kinetics and 10-min total described as essentially identical to its own matched controls; `n=6` glands/genotype; Fig. 1B | `MEASUREMENT` | Independent genotype `VALIDATION`; no equality margin was reported |
| `P15-10` | AE2-line resting Cl | Control 53.40 +/- 1.80 and null 54.50 +/- 1.80 mM; `n=6,6`; Table 1 | `MEASUREMENT` | AE2 `VALIDATION`; do not pool with the Ae4-line WT because strains differ |
| `P15-11` | AE2-line resting pHi | Control 6.87 +/- 0.01 (`n=4`) and null 6.95 +/- 0.05 (`n=6`); Table 1 | `MEASUREMENT` | AE2 `VALIDATION` |
| `P15-12` | AE2-line combined-stimulus Cl uptake | Control 2.11 +/- 0.50 and null 2.39 +/- 0.20 x10^-3 s^-1; `n=13,18`; Table 1 | `MEASUREMENT` | AE2 `VALIDATION` |
| `P15-13` | NKCC1 assay | Isolated uptake not detectably different among pooled controls (`n=16`), AE4-null (`n=8`), and AE2-null (`n=10`); 50 uM bumetanide reduced the isolated signal 95.6 +/- 0.9%; Fig. 2 | `MEASUREMENT` | Genotype-invariant NKCC1 constraint; not proof that all compensation is absent |
| `P15-14` | NHE-associated pH response | EIPA-sensitive stimulation-associated alkalinization showed no detected genotype effect; Fig. 2C-D | `MEASUREMENT` + `INFERENCE` | Genotype-invariant NHE constraint. The paper explicitly notes that this did not directly measure carbonic-anhydrase activity. |
| `P15-15` | Salivary HCO3 output | Ae4 WT 5.3 +/- 1.4 (`n=7`) vs null 7.1 +/- 0.8 (`n=7`) microequiv/10 min, `p=0.28`; Ae2 control 4.0 +/- 0.4 (`n=5`) vs null 6.3 +/- 1.8 (`n=6`), `p=0.28` | `MEASUREMENT` | Genotype `VALIDATION`; a nonsignificant difference is not exact equality |
| `P15-16` | AE4 membrane domain | Antibody specificity failed against null tissue; basolateral localization was inferred from transport direction and an apical contribution was not excluded | `INFERENCE` | Orientation uncertainty; no measured membrane fraction |
| `P15-17` | Isolated exchanger contributions | With bumetanide/T16Ainh-A01 isolation, activity was 38.2 +/- 6.6% slower in Ae4-null and 55.1 +/- 3.4% slower in Ae2-null than pooled controls (`n=7/8` Ae4 control/null; `n=8/10` Ae2 control/null); combined Ae2/Ae4 deletion reduced activity 87.7 +/- 3.1% (`n=8` control, `n=6` double null); Figs. 4--5 | `MEASUREMENT` + `REPORTED_FIT` | Genotype `VALIDATION`; residual exchange is not zero, and AE2 assay contribution must not be converted into a secretion contribution |
| `P15-18` | IPR regulation in isolated exchanger assay | Source control acini: 1.56 +/- 0.19 without IPR (`n=10`) versus 4.33 +/- 0.91 with 5 uM IPR (`n=13`) x10^-3 s^-1, `p=0.016`; Ae4-null: 0.95 +/- 0.10 (`n=8`) versus 1.08 +/- 0.09 (`n=9`), `p=0.366` | `MEASUREMENT` | WT/control beta-arm `CALIBRATION` and Ae4-null transporter `VALIDATION`; this is an endpoint assay contrast, not a cAMP/PKA time course |

The strict secretion result is deliberately absent from this table; it is
described only in the held-out file. The permitted WT-only flow constraint is
the Ae4-line WT trace in Fig. 1A: approximately 9--10 uL/min during the 10-min
CCh+IPR interval (`n=6`, graphical means with plotted SEM). The paper does not
tabulate minute-wise WT values.

## `P16`: AE4 Na/K transport mechanism

All uncertainties are SEM. Unless a native acinus is explicitly named, these
are CHO-K1 assays and require a context-specific observation model.

| ID | Result | Exact value or scope | Evidence class | Task 13B role |
|---|---|---|---|---|
| `P16-01` | Native Na dependence | In `Ae2-/-` SMG acini, low external Cl produced alkalinization with Na (`n=8`), little response with NMDG (`n=7`), and little response in `Ae2-/-;Ae4-/-` acini with Na (`n=8`); Fig. 1 | `MEASUREMENT` | Native transporter `CALIBRATION` |
| `P16-02` | Recombinant Cl/HCO3/Na coupling | Mouse Ae4 BCECF response with Na (`n=11`), not Na-free (`n=6`) or HCO3-free (`n=8`); SPQ Cl loss with Na (`n=10`), not HCO3-free (`n=7`) or Na-free (`n=6`); human AE4 remained Na dependent (`n=5`); Fig. 2 | `MEASUREMENT` | Transporter `CALIBRATION` |
| `P16-03` | Direct Na signal | SBFI showed intracellular Na increase with outward Cl gradient in Ae4 CHO (`n=15`) but not nontransfected cells (`n=13`); the signal returned after Cl restoration; Fig. 3 | `MEASUREMENT` + source-qualified `INFERENCE` | Direction/sign `CALIBRATION`; the return is not a clean AE4 reversal measurement because the authors explicitly note that Na/K-ATPase could produce it; no reversal potential was measured |
| `P16-04` | Near-electroneutral transport | Transport-associated current-clamp voltage change was <2 mV; Ae4 Cl-exit `n=7`, nontransfected `n=8`; Fig. 4 | `MEASUREMENT` + `INFERENCE` | Near-zero net-charge gate |
| `P16-05` | Voltage independence | Alkalinization rates at -100 and 0 mV were essentially identical (`n=5` each); Fig. 4E-F | `MEASUREMENT` | Electroneutrality gate |
| `P16-06` | Macroscopic current controls | At -100 mV: Ae4 -1.68 +/- 0.44 vs nontransfected -1.27 +/- 0.46 pA/pF (`n=5,5`, `p=0.54`). At 0 mV: -0.28 +/- 0.09 vs -0.21 +/- 0.14 pA/pF (`n=5,5`, `p=0.69`). | `MEASUREMENT` | Current-resolution constraint, not proof of a unique stoichiometry |
| `P16-07` | Direct K-associated signal | In Na-free conditions, PBFI showed sustained K loss under inward Cl/outward K gradients in Ae4 CHO (`n=6`) and little response in nontransfected cells (`n=8`); Fig. 8C-D. PBFI is only about 1.5-fold K-selective over Na. | `MEASUREMENT` | K direction `CALIBRATION`; not K stoichiometry |
| `P16-08` | Monovalent-cation permissiveness | Alkalinization with Na (`n=10`), K (`n=8`), Li/Rb/Cs (`n=6` each), not NMDG (`n=6`); Fig. 8 | `MEASUREMENT` | Qualitative cation-branch constraint |
| `P16-09` | Na dose fit | EC50 49 mM; Hill 2.0; Rmin 0.3 and Rmax 1.5 x10^-3 s^-1; Fig. 9 and text | `REPORTED_FIT` | Effective dose-response `CALIBRATION`; Hill value is not a microscopic site count |
| `P16-10` | K dose fit | EC50 62 mM; Hill 1.8; Rmin 0.4 and Rmax 1.6 x10^-3 s^-1; Fig. 9 and text | `REPORTED_FIT` | Effective dose-response `CALIBRATION` |
| `P16-11` | Dose-fit sampling limitation | 5, 25, 50, 100, 125, 150 mM cation; at least four experiments/concentration and three electroporations. No zero-cation point, parameter SE/CI, or covariance was reported. | `NOT_MEASURED` | Do not create precision or independent-fit errors |
| `P16-12` | Cl-free apparent NBC-like activity | Pipette O (15 mM Cl, no HCO3): 0.12 +/- 0.03 (`n=6`); P (4 mM Cl, no HCO3): 0.02 +/- 0.01 (`n=6`); Q (4 mM Cl, 15 mM HCO3): 0.08 +/- 0.01 (`n=7`) x10^-3 s^-1; P differs from O (`p=0.01`) and Q (`p=0.001`); Table 4 | `MEASUREMENT` + `INFERENCE` | Requires an internal exchangeable anion; OH or carbonate substitution was not excluded |
| `P16-13` | Candidate electroneutral stoichiometries | `Cl:Na:HCO3` = 1:1:2, 1:2:3, and 2:1:3 were evaluated in free-energy calculations | `MODEL_HYPOTHESIS` | Thermodynamic candidates only; no stoichiometry was measured |
| `P16-14` | Na-only thermodynamic thresholds | For the three candidates, Cl influx required respectively Cl_i <7, <0.5, <26 mM; Na_i >104, >143, >54 mM; or HCO3_i >49, >84, >29 mM under the source's fixed concentrations; Fig. 7 | `REPORTED_FIT` conditional on `MODEL_HYPOTHESIS` | Reversal/sign audit; not assay data |
| `P16-15` | Na/K-pooled thermodynamic proposal | Under assumed 1:1:2 and nearly equal Na/K handling, calculated Cl uptake was favorable for Cl_i <73 mM or HCO3_i >16 mM; Fig. 10 | `MODEL_HYPOTHESIS` + `REPORTED_FIT` | Must not be promoted to a shared-pool transport law without explicit state/thermodynamic testing |

The source does not measure the complete alternating-access sequence,
transporter density, microscopic rates, native stimulated Na/K trajectories,
or a unique HCO3-versus-CO3 substrate/stoichiometry.

## `P21`: mandatory beta/cAMP/PKA regulation

All uncertainties are SD. This source establishes the regulatory chain at the
level of perturbation-response logic; it does not identify a physical dynamic
subsystem.

| ID | Result | Exact scope | Evidence class | Task 13B role |
|---|---|---|---|---|
| `P21-01` | Native beta-adrenergic activation | WT SMG acini: control `n=11` experiments vs 5 uM IPR `n=12`; IPR increased isolated exchanger activity, `p<0.001`; Fig. 1B | `MEASUREMENT` | Native regulatory `CALIBRATION`; numeric means are figure-only and not tabulated |
| `P21-02` | Genetically isolated native AE4 response | `Ae2-/-` acini: control `n=12`, IPR `n=11`, IPR + 10 uM H89 `n=7`; IPR increased activity (`p<0.0001`) and H89 prevented the increase (`p<0.05` comparison); Fig. 1C | `MEASUREMENT` + `INFERENCE` | Mandatory beta/PKA hierarchy; H89 is explicitly nonselective |
| `P21-03` | Adenylate-cyclase perturbation | AE4 CHO: control `n=40`, 10 uM forskolin `n=39`, forskolin + 10 uM H89 `n=41`; the prose calls the forskolin increase approximately 25%, whereas a visual read of the plotted means is approximately 0.36 to 0.63 pH x10^-3 s^-1 (about 1.7-fold); Fig. 2 | `MEASUREMENT` + conflicting approximate `REPORTED_FIT` | Qualitative/sensitivity-only `CALIBRATION`. Preserve the prose-versus-figure discrepancy; neither approximation is licensed as an exact point target or universal AE4 multiplier. |
| `P21-04` | Constitutive PKAc perturbation | AE4+PKAc vs AE4+K73M dominant-negative PKAc: BCECF `n=18,16`; SPQ `n=8,9`; both readouts increased with PKAc (`p<0.005`); PKA constructs alone did not generate the flux; Fig. 3 | `MEASUREMENT` | PKA-dependent direction `CALIBRATION`; no numeric fold is tabulated |
| `P21-05` | Mutant basal function/localization | WT/S173A/S273A basal activity `n=42,37,46`; both mutants were functional and membrane localized but had reduced basal activity; Fig. 5 | `MEASUREMENT` | Qualitative transporter/regulation `VALIDATION` |
| `P21-06` | S173 dependence | S173A+DN/PKAc `n=14,15`; no detected PKAc activation; Fig. 6A-B | `MEASUREMENT` + `INFERENCE` | Regulatory model `VALIDATION`; S173 is necessary for this construct-level response |
| `P21-07` | S273 response retained | S273A+DN/PKAc `n=12,12`; PKAc activation retained (`p<0.01`); Fig. 6C-D | `MEASUREMENT` | Regulatory model `VALIDATION` |
| `P21-08` | Direct S173 phosphorylation | Not demonstrated; the authors explicitly say so | `NOT_MEASURED` | Do not state that PKA directly phosphorylates S173 |
| `P21-09` | Regulatory kinetics | No cAMP concentration trace, PKA activity trace, AE4 phosphorylation trace, activation/deactivation fit, washout/reversibility time constant, or reported physical delay. IPR/H89/forskolin pretreatment durations are not specified quantitatively. | `NOT_MEASURED` | The source cannot identify cAMP, PKA, phosphorylation, or dephosphorylation time constants |
| `P21-10` | Microscopic modulation site | Increased turnover, active surface number, altered stoichiometry, or another S173-dependent interaction remain alternatives in the Discussion | `MODEL_HYPOTHESIS` | Common-capacity versus transition-specific regulation is not resolved |

The time axes in Figs. 1--3 and 5--6 show the imposed high/low-Cl assay
response from which an initial slope was computed. They do **not** show AE4
activity beginning at the instant of beta stimulation and therefore cannot be
used to calibrate a beta-to-cAMP-to-PKA delay.

## `P25`: mutagenesis measurements versus structural hypotheses

All uncertainties are SEM. These HEK-293 results are retained primarily as
transporter-module validation because they were obtained in a different cell
line, expression construct, temperature, and observation model from the 2016
fits and native SMG physiology.

| ID | Result | Exact scope | Evidence class | Task 13B role |
|---|---|---|---|---|
| `P25-01` | WT AE4 activity/cation dependence | WT `n=17` vs nontransfected `n=16`; Na-dependent condition `n=6` vs nontransfected `n=5`; Fig. 3 | `MEASUREMENT` | Heterologous validation; readout is HCO3-linked alkalinization, not direct Na/K concentration |
| `P25-02` | Single-mutant hierarchy | Versus WT `n=17`: T448I (`n=17`) and T448G (`n=10`) about 50% lower; D709A (`n=19`), S446A (`n=8`), T756A (`n=9`) about 30% lower; I758G (`n=11`), I758R (`n=10`), D709N (`n=12`), S447A (`n=9`), T754A (`n=9`) no detected change; Fig. 4 | `MEASUREMENT` + approximate `REPORTED_FIT` | Mutation-rank `VALIDATION`; no unique microscopic rate mapping follows |
| `P25-03` | Double mutants in Na | S446A-T448I (`n=18`) and D709A-T448I (`n=19`) about 50% below WT; T756A-T448I (`n=18`) comparable to nontransfected; Fig. 5 | `MEASUREMENT` | Strong Na-context epistasis `VALIDATION` |
| `P25-04` | Membrane-expression control | Mander coefficient >0.8 for WT; no detected localization or immunoblot-abundance difference among tested mutants | `MEASUREMENT` | Supports functional rather than gross trafficking interpretation; not equal molecule counts at the surface |
| `P25-05` | Na/K nonidentity | WT, S446A, T448I, T448G similar under Na and K. D709A and T756A about 30% lower in K than Na in Fig. 6 caption (Discussion calls the extra decrement about 20%). T756A-T448I has no detectable Na response but substantial K response comparable to single mutants. NMDG abolishes responses. | `MEASUREMENT` + `INFERENCE` | Qualitative Na/K branch-asymmetry `VALIDATION`; preserve 20/30% wording uncertainty |
| `P25-06` | Voltage-independence control | Valinomycin/high-K depolarization did not change alkalinization rate; membrane-potential change was checked with DiO/DPA; Fig. S4 | `MEASUREMENT` | Electroneutrality validation |
| `P25-07` | Structural template and residues | Outward-open rat NDCBE cryo-EM structure PDB 7RTM was used as the template; Na/K-HCO3 pairs were manually placed. D709/T713 cation and G449/K879/T448/T756 bicarbonate contacts are simulation results. | `MODEL_HYPOTHESIS` | No direct AE4 structure or binding measurement |
| `P25-08` | MD sampling | Three 120-ns simulations/system, first 20 ns restrained; subsequent analysis used selected stable replicas/monomer chains, not every initial trajectory | `MODEL_HYPOTHESIS` | Qualitative hypothesis generation only |
| `P25-09` | MD H-bond occupancy | WT-Na D709/G449/K879/T448/T756 = 89.2/90.1/97.8/87.3/57.0%; double-mutant-Na = 66.5/54.8/58.6%; WT-K G449/K879/T448/T756 = 31.7/81.5/77.3/64.6%; double-mutant-K D709/G449/K879 = 42.0/20.0/86.6%; Table 2 | Simulation `REPORTED_FIT` | Not experimental calibration |
| `P25-10` | Substrate/stoichiometry ambiguity | 1:1:2 Cl:cation:HCO3 and 1:1:1 Cl:cation:CO3 remain possible. Proposed sequential 1:1:1:1 Cl:Na:HCO3:K and 1:1:1:2 Cl:Na:CO3:K cycles were not directly tested; no Cl-containing MD cycle or alternating-access trajectory was simulated. | `MODEL_HYPOTHESIS` | Do not hard-code the proposed sequential order as fact |
| `P25-11` | State ordering | The authors speculate that T448 and T756 may act at different stages | `MODEL_HYPOTHESIS` | Mutant nonadditivity is measured; temporal state order is not |

## `A18`: apical cation topology

The accessible primary abstract supports two qualitative measurements:

1. Local Ca photolysis close to the apical membrane of **parotid** acinar
   cells triggered a significant K current. This supports a functional apical
   Ca-activated K pathway, but it does not quantify its fraction of total
   conductance.
2. Na/K-ATPase labeling appeared distributed throughout the plasma membrane,
   including the apical pole. This supports a nonzero apical pump population,
   but it does not quantify functional pump capacity at either membrane.

The phrase that pump labeling is distributed throughout the membrane is a
localization/density statement; it does not mean that 50% of all pump molecules
or functional capacity is apical. `PB94` measured rat-parotid membrane area per
cell volume as `0.125 +/- 0.027 um2/um3` luminal and
`1.54 +/- 0.085 um2/um3` basolateral (`n=4` animals). Combining uniform
surface-density localization with those areas gives the derived transfer

\[
f_{P,a}=\frac{A_a}{A_a+A_b}=0.075075,
\]

with deliberately conservative low/high combinations `0.056878` and
`0.094586`. These are **derived parotid-to-mouse-SMG sensitivity values**, not
a direct measurement of functional pump turnover in the P15 preparation.

The official mathematical supplement says that the model **chose** 30:70
apical:basolateral Na/K-pump and 40:60 apical:basolateral K-channel
distributions. These are `HISTORICAL_ASSUMPTION` values motivated by the
observations, not measured fractions. They may define sensitivity cases or
priors, but may not be called source-measured calibration targets. Transfer
from parotid to the 2015 submandibular protocol is itself a declared modeling
decision.

The production pump sensitivity is therefore the derived
`0.056878/0.075075/0.094586` panel. The `0.20/0.30/0.40` apical-K routing panel
is a model/source bracket, not a measured confidence interval. Likewise the
AE4 mixed-bath Na routing values `0.05/0.10/0.20` are modeling assumptions,
not measured native fractions.

## Native whole-cell channel source map and calcium coordinate

The strict native-source reconstruction replaces the old common-conductance
scale with literal, separately sourced maxima:

| Quantity | Production interpretation | Evidence status |
|---|---:|---|
| apical CaCC maximum | `31.4 nS` | Palk/Gin parotid lineage, traced to the average whole-cell maximum reported from Arreola rat-parotid experiments |
| total Ca-activated K maximum | `14 nS` | Palk lineage, traced to Thompson/Begenisich whole-cell parotid experiments; partitioned between membranes without changing the total |
| Ca gate | `C^1.46 / (0.26^1.46 + C^1.46)` with `C` in `uM` | published Palk-lineage gate; voltage-independent production reference |

At `C=0.058, 0.10, 0.12, 0.15 uM`, that gate is respectively
`0.10062056`, `0.19860330`, `0.24437206`, and `0.30936970`; the corresponding
gated fractions of the `31.4 nS` maximum are `3.1595`, `6.2361`, `7.6733`,
and `9.7142 nS`. The calcium coordinate is therefore dimensionful. `0.058 uM`
is a source-linked rat-lineage resting transfer. `0.10`, `0.12`, and
`0.15 uM` are explicitly assumed spatially averaged WT sensitivity inputs;
they were not direct free-calcium measurements under the matched `0.3 uM`
CCh protocol and were never fitted to AE4-null secretion.

The native mouse-SMG CaCC calcium sensitivity in `R10` is retained only as a
diagnostic alternative because its approximately `0.126 uM` scale is
voltage/protocol dependent. It cannot replace the Palk gate in production
without an explicit voltage/observation map. The 2018 lineage's larger
calibrated conductances (`71.3 nS` CaCC and `30.4 nS` K) remain historical
model values. Multiplying either those values or provisional code literals by
one common scale is not a native conductance source map. No membrane-area
factor may be applied again to an already whole-cell conductance.

## Absolute gland geometry and single-cell flow gate

`K15` defines true acinar volume as gland weight times acinar fraction and
reports `66.91%` acinar tissue for mouse SMG. A deliberately generous
one-SMG upper sensitivity of `100 mg`, together with the declared
`1 mg approximately 1 uL` tissue-volume conversion and the model's `1.3 pL`
cell volume, gives

\[
V_{acini,max}=66.91\ \mathrm{uL},\qquad
N_{cell,max}=51,469,231.
\]

`K19` confirms that the corresponding ex vivo workflow separately collects
the ipsilateral gland; a paired-gland interpretation is retained only as a
sensitivity, not the production geometry. At the P15 WT graphical envelope
of `9--10 uL/min`, the one-SMG upper bound requires per-cell sustained flow
between `0.00291436` and `0.00323818 pL/s`, or specific flow between
`0.00224182` and `0.00249091 s^-1` at `1.3 pL`. A paired-gland sensitivity
halves those per-cell thresholds.

The superseded joint G5 ensemble implied `2.485e8--3.401e8` cells across its
retained rows, or approximately `323--442 uL` of modeled cells. Its two
selected profiles implied `294,706,403` and `260,704,116` cells. This is a
hard `ABSOLUTE_GLAND_GEOMETRY_SCALE_FAIL`, not a warning that can be waived by
normalization.

## `C15V`: beta/volume-sensitive anion exit

The adult acinar-specific Tmem16A study supplies a second, distinct beta-arm
constraint. Under `5 uM` IPR for `10 min`, secretion persisted after Tmem16A
loss; low extracellular Cl reduced the response by about `98%`, HCO3-free
conditions reduced it, and NPPB/DCPIB inhibited it. IPR caused a slow
`12.5 +/- 0.2%` acinar volume increase, whereas CCh caused a rapid
`10 +/- 1%` shrinkage. IPR- and swelling-activated macroscopic anion currents
had similar magnitude and DCPIB sensitivity. These are direct observations;
assigning one molecular VRAC identity or an exact apical conductance is not.

For a secreting acinar model this evidence requires a beta/volume-sensitive
anion-exit topology distinct from CaCC. PKA-dependent AE4 regulation is a
basolateral loading process and cannot substitute for that exit. However, the
published current magnitude is graph-only and no valid whole-cell
current-to-conductance map is available. In particular, the `15 +/- 2 pF`
capacitance in Romanenko et al. 2008, DOI `10.1152/ajpgi.90384.2008`, is a
salivary **duct-cell** value and is excluded from acinar scaling.

A nested `V_beta` term proportional to
`g_V * beta_input * max(V_i/V_rest - 1, 0)` is therefore diagnostic only. It
adds no arbitrary delay state. If an exact-rest initialization cannot produce
swelling before this term activates, the topology deadlocks and must be
rejected; an epsilon gate is not a physiological rescue. `V0` (no such exit)
is an exact numerical control but topology-incomplete for IPR-only secretion,
while a nonzero `V1` is production-ineligible until its absolute capacity is
measured or validly mapped.

## `V18`: historical mathematical lineage

The 2018 paper is a model record. Its historical state/equation inventory is:

- intracellular Na, K, Cl, HCO3, H, CO2 amounts;
- cell volume;
- luminal Na, K, Cl concentrations at fixed lumen volume;
- apical and basolateral membrane potentials replaced by two QSS current
  constraints for simulation;
- basolateral NKCC1, NHE1, AE2, AE4, Na/K pump, and Ca-activated K channel;
- apical Ca-activated Cl channel; paracellular Na/K/water; and membrane water
  fluxes;
- a prescribed qualitative Ca input rather than the matched 2015 CCh+IPR
  protocol.

The paper explicitly omitted apical K channels for simplicity and placed the
Na/K pump basolaterally. It calibrated conductances/transporter densities to a
mixed-source resting table. Historical targets included intracellular Cl 50.1
mM, K 120 mM, Na 20 mM (model result 25), HCO3 12 mM, pHi 6.8 (model result
6.91), cell volume 1.3 pL, apical voltage -50.2 mV, and basolateral voltage
-61.8 mV. These values retain their original mixed tissue/species/model
contexts; only P15's matched SMG Cl and pH values enter the modern calibration
ledger directly.

The historical model reported about a 24% AE4-knockout flow reduction and a
knockout timing mismatch. These are historical **model outputs**, not modern
validation data. The paper's narrative rounded restatement of the 2015
phenotype is not a second experimental endpoint; the controlling experimental
outcome remains sealed in `heldout_targets.csv`.

## Calibration/validation/holdout firewall

| Evidence group | Machine-readable location | Authorized use |
|---|---|---|
| WT resting Cl/pH, WT-only flow envelope, WT/native transporter activity, 2016 transporter assays, 2021 WT regulatory contrasts, qualitative apical topology | `calibration_targets.csv` | Matching-context calibration or hard/qualitative gates; no cross-system absolute-rate pooling without an observation model |
| AE2 deletion, conditional AE4-null ionic phenotypes, genotype-invariant NKCC/NHE checks, S173/S273 logic, 2025 mutation hierarchy | `validation_targets.csv` | Out-of-objective validation. Rows marked `CONDITIONAL_LOCALIZATION` require the staged release condition. |
| Strict AE4-null stimulated-saliva outcomes | `heldout_targets.csv` | `STRICT_HELDOUT`; evaluate only after WT state, WT dynamics, topology, regulatory dynamics, physical time map, solver settings, and parameter manifest are frozen |
| Reveal state and file hashes | `reveal_log.json` | Audit only; every future access/reveal must be appended with reason, actor, and frozen-model hash |

## Decision-critical missing measurements and explicit gaps

- No audited source directly measures native stimulated intracellular Na or K
  trajectories, AE4-null acinar volume, native AE4 surface density, or
  transporter turnover.
- `P16` does not identify a unique stoichiometry, complete transport sequence,
  microscopic rate constants, HCO3 versus CO3 substrate, or uncertainty of its
  Hill-fit parameters.
- `P21` requires dynamic beta/cAMP/PKA regulation in the biological model, but
  it does not identify a cAMP state, PKA activity state, direct S173
  phosphorylation, activation/deactivation times, washout kinetics, or whether
  PKA changes common capacity versus selected AE4 transitions. Its prose
  approximately 25% forskolin effect is also internally inconsistent with a
  visual reading of the plotted group means (approximately 1.7-fold); raw data
  are unavailable, so the contrast is qualitative/sensitivity-only.
- `P25` constrains mutation and Na/K asymmetry in HEK-293 cells. It does not
  provide a native structure, direct residue-ion binding, complete cycle, or
  physiological membrane orientation.
- `A18` supports apical topology in parotid acini, but the accessible primary
  record supplies no quantitative functional apical fraction and does not
  establish the same fraction in mouse SMG acini.
- `C15V` establishes a distinct beta/volume-sensitive anion-exit requirement,
  but not an absolute acinar conductance. This is a substantive full-model
  limitation rather than a license to tune a current from the AE4-null curve.
- The old G5 flow normalization violates the `K15` one-SMG acinar-volume
  bound and its common conductance scale violates the native channel source
  map. Those failures supersede its pre-reveal authorization.
- Exact WT minute-wise secretion values and their numerical SEM are not
  tabulated by `P15`; only the exact experimental protocol, graphical WT trace,
  and strict integrated genotype comparison are available.
- No negative finding is encoded as equality. "No detected difference" is a
  validation statement with the source's power and assay context, not a zero
  residual with arbitrary precision.
