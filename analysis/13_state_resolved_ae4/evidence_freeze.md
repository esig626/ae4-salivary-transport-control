# Task 13 experimental evidence freeze

**Status: FROZEN before Task 13 fitting.** Frozen 2026-08-27 UTC on branch
`codex/state-resolved-ae4-cycle`, starting commit `7846706645f9`. Downstream
work may add a separately documented amendment, but must not silently change a
value, classification, or holdout assignment in this table.

This document separates observations from their interpretation. It does not
specify an AE4 model.

## Classification and use keys

| Key | Meaning |
|---|---|
| **DM** | Direct assay-level measurement: the stated instrument readout or a calibrated quantity derived by the source's stated assay. This does not make the biological interpretation unique. |
| **RF** | A reported regression, fit, ratio, or calculation from measurements. |
| **DI** | Data-supported inference; plausible from the observations but not itself measured. |
| **MH** | Mechanistic, structural-model, MD, or physiological hypothesis. |
| **ND** | Not measured in the cited experiment. |
| **CAL-T** | May constrain/calibrate a transporter-level model. |
| **CAL-WT** | May constrain/calibrate the WT whole-cell protocol or state. |
| **VAL-A** | Stage A holdout: reveal only after the transporter and WT candidate are frozen. |
| **LOC-B** | Still withheld in Stage A; may be released only as a Stage B missing-balance localization constraint. |
| **STRICT** | Strict secretion holdout in both stages; never optimize to it. |
| **AUX-VAL** | Independent genotype/topology check, not a reason to tune an AE4-null secretion prediction. |
| **CONTEXT** | Qualitative structural/context constraint or limitation, not a fitted numerical target. |

`n` retains the source's unit. In particular, a source may report animals,
experiments, or cells; those are not interchangeable. Error conventions also
differ: 2015 and 2016 report mean +/- SEM; 2021 reports mean +/- SD; 2025
reports mean +/- SEM.

## Assay-context freeze

| Source | Genotype/species/system | Bath, temperature, and stimulation | Sampling and uncertainty |
|---|---|---|---|
| [Pena-Munzenmayer et al. 2015, JBC](https://doi.org/10.1074/jbc.M114.612895), Figs. 1-5 and Table 1 | Female and male mice, 2-4 months. Systemic `Ae4-/-` on BS/129Svj and acinus-specific `Ae2fl/fl-AQP5/ACID-Cre` on C57BL/6, each compared only with its own littermate controls; double-null acini were also used. Ex vivo perfused submandibular gland (SMG) and native isolated SMG acini. | Gland: common-carotid perfusion with HCO3-containing high-Cl solution at 37 C; 0.3 uM carbachol (CCh) + 5 uM isoproterenol (IPR) for 10 min. Acini: 37 C, 95% O2/5% CO2; high-Cl bath 4.3 KCl, 120 NaCl, 25 NaHCO3, 5 glucose, 10 HEPES, 1 CaCl2, 1 MgCl2, pH 7.4; low-Cl gluconate substitutions. SPQ and BCECF assays used protocol-specific inhibitors as stated below. | Flow recorded every minute; `n=6` glands/genotype. Cell experiments report experiment-level `n` in Table 1. Mean +/- SEM. The paper explicitly warns that Ae4-line and Ae2-line WT kinetics differ, likely because of strain/background, so cross-line WT pooling is prohibited. |
| [Pena et al. 2016, JGP](https://doi.org/10.1085/jgp.201611571), Figs. 1-10 and Tables 1-4 | Female and male mice, 2-4 months; native acini from acinus-specific `Ae2-/-` and `Ae2-/-;Ae4-/-` SMG. CHO-K1 heterologous mouse Ae4 variant 3 or human AE4 variant 2, 18-20 h after electroporation. | Imaging at 37 C; patch/imaging at 22 C. Standard high/low-Cl solutions contained 145 mM Na (or 145 mM K in K protocols), 25 mM HCO3, and 128.3/4 mM Cl at pH 7.4, with explicit NMDG and HCO3-free controls. Dose response used 5, 25, 50, 100, 125, and 150 mM Na or K. No secretagogue stimulation. | BCECF, SPQ, SBFI, PBFI, current clamp, and voltage clamp. At least three mice/electroporations per condition; figure captions report experiment or cell `n`. Mean +/- SEM. Fit-parameter confidence intervals were not reported. |
| [Pena et al. 2021, AJP-GI](https://doi.org/10.1152/ajpgi.00145.2021), Figs. 1-6 | Seventeen female mice, 8-12 weeks (`10` acinus-specific `Ae2-/-`, `7` controls); native isolated SMG acini. CHO-K1 expressing mouse Ae4, WT or S173A/S273A, with PKA constructs. | 37 C, 95% O2/5% CO2. Native high/low-Cl uptake in 50 uM bumetanide + 50 uM T16Ainh-A01; 5 uM IPR; 10 uM H89. CHO assays used 10 uM forskolin or constitutively active PKAc versus K73M dominant-negative PKAc. | Native captions report experiments and Methods require at least three animals/condition; CHO captions report cells from at least three electroporations. Mean +/- SD. No activation time constant was fitted or reported. |
| [Catalan et al. 2025, AJP Cell](https://doi.org/10.1152/ajpcell.00346.2024), Figs. 3-8 and Table 2 | Human AE4 isoform 2, WT and mutants, transiently expressed in HEK-293 cells for 18-24 h. No salivary acinar cells or knockout animals. | BCECF assay at room temperature. High/low-Cl baths: 145 mM Na (or K), 25 mM HCO3, and 128.3/4 mM Cl, pH 7.4, gassed 95% O2/5% CO2. NMDG substitution tested cation dependence. No receptor stimulation. | At least three cells/experiment and at least three independent transfections (at least five for Fig. 5); `n` in captions is cells. Mean +/- SEM. Regions of interest were responding transfected cells, which limits population-level interpretation. |
| [Almassy et al. 2018, Pflugers Archiv](https://doi.org/10.1007/s00424-018-2109-0) | Parotid acinar-cell single-cell electrophysiology/Ca imaging and Na/K-pump localization; the accessible abstract does not expose the species, exact bath, or sample sizes for these experiments. | Local Ca photolysis near the apical membrane; pump immunolocalization. Exact conditions could not be verified because the version-of-record main text is subscription restricted. | Publisher abstract and publisher supplements were accessible. The supplement is the mathematical model, not the missing experimental Methods. Consequently only abstract-level qualitative experimental claims are frozen here. |
| [Vera-Siguenza et al. 2018, Bull. Math. Biol.](https://doi.org/10.1007/s11538-017-0370-6) | Historical mathematical model and archived implementation; not new experimental data. | Historical model protocol. | **Historical model evidence only.** Its predictions are never calibration data for biological parameters. |

## Frozen evidence table

### 2015 native gland and acinar phenotypes

| ID | Quantity or constraint | Exact source anchor | Class | Numerical value, uncertainty, and `n` | Use |
|---|---|---|---|---|---|
| E15-01 | WT stimulated whole-gland flow scale and nearly sustained Ae4-line WT response | 2015 Fig. 1A | DM | About 9-10 uL/min through the 10-min CCh+IPR interval; `n=6`; plotted SEM. Raster-digitized values are listed below and are not source-tabulated. | **CAL-WT**, with the exact 2015 gland protocol and strain only |
| E15-02 | AE4-null total saliva deficit | 2015 Results, Fig. 1A inset | DM + RF | Total saliva collected over 10 min was **35 +/- 4.7% less** than littermate Ae4 WT; `n=6` glands/genotype, mean +/- SEM, Student t test `p<0.05`. Algebraic complement is KO/WT `0.65 +/- 0.047`, but that ratio/SEM was not separately fitted by the paper. | **STRICT** |
| E15-03 | AE4-null temporal phenotype | 2015 Results and Fig. 1A | DM | First 2-3 min comparable to WT; sustained deficit after 3 min; later pointwise flow less than half WT. Minute means are approximately digitized below; `n=6`/genotype, plotted SEM. | **STRICT** |
| E15-04 | AE2-null whole-gland secretion | 2015 Fig. 1B | DM | Kinetics and 10-min total described as essentially identical to the line's own littermate controls; `n=6`/genotype. No exact null/control ratio or equivalence margin is reported. | **AUX-VAL**; do not treat failure to reject as proof of exact equality |
| E15-05 | Resting intracellular Cl in Ae4-line WT | 2015 Table 1, calibrated SPQ | DM | `50.10 +/- 1.50 mM`, `n=7` | **CAL-WT** |
| E15-06 | Resting intracellular Cl in AE4-null acini | 2015 Table 1 | DM + RF | `36.50 +/- 1.60 mM`, `n=6`; reported as `27 +/- 3.1%` lower than WT, `p<0.001` | **LOC-B** (withheld in Stage A) |
| E15-07 | Resting pHi in Ae4-line WT and AE4 null | 2015 Table 1, BCECF | DM | WT `6.91 +/- 0.07` (`n=4`); null `6.89 +/- 0.02` (`n=4`); no detected difference | WT is **CAL-WT**; null is **LOC-B** |
| E15-08 | Combined CCh+IPR initial Cl uptake | 2015 Table 1, SPQ slope | DM + RF | WT `2.02 +/- 0.10`; AE4 null `0.90 +/- 0.09`, in `10^-3 s^-1`, `n=10,9`; reported `55.4 +/- 4.5%` lower, `p<0.001` | WT **CAL-WT**; null **LOC-B** |
| E15-09 | Agonist separation of AE4-null Cl uptake | 2015 Table 1 | DM | CCh only: WT `2.18 +/- 0.20`, null `2.30 +/- 0.10` (`n=7,6`; no detected difference). IPR only: WT `0.40 +/- 0.07`, null `0.20 +/- 0.03` (`n=4,7`; `p<0.05`). Units `10^-3 s^-1`. | WT **CAL-WT**; null **LOC-B**; supports beta/cAMP association, not a measured PKA reaction |
| E15-10 | Initial stimulation-induced Cl exit | 2015 Table 1 | DM | CCh+IPR normalized SPQ change: WT `-0.21 +/- 0.01`, null `-0.18 +/- 0.02` (`n=10,9`), no detected difference. | Null **LOC-B**; assay-level negative constraint on gross initial Cl exit, not proof that every apical-channel property is unchanged |
| E15-11 | Ae2-line resting Cl, pHi, and combined-stimulus Cl uptake | 2015 Table 1 | DM | Control/null: Cl `53.40 +/- 1.80` (`n=6`) / `54.50 +/- 1.80` (`n=6`); pHi `6.87 +/- 0.01` (`n=4`) / `6.95 +/- 0.05` (`n=6`); uptake `2.11 +/- 0.50` (`n=13`) / `2.39 +/- 0.20` (`n=18`) `10^-3 s^-1`; no detected differences. | **AUX-VAL**, within C57BL/6 Ae2 line only |
| E15-12 | NKCC1 functional compensation assay | 2015 Fig. 2 | DM | HCO3-free/EZA/T16Ainh-A01-isolated uptake was not detectably different among control (`n=16`), AE4 null (`n=8`), and AE2 null (`n=10`); 50 uM bumetanide reduced the signal `95.6 +/- 0.9%`. | **LOC-B/CONTEXT**. This is an assay-specific activity result, not a measurement of expression and not proof that all compensation is absent. |
| E15-13 | HCO3 secretion | 2015 Results | DM | Microequiv/10 min: Ae4 WT `5.3 +/- 1.4` (`n=7`) versus null `7.1 +/- 0.8` (`n=7`), `p=0.28`; Ae2 control `4.0 +/- 0.4` (`n=5`) versus null `6.3 +/- 1.8` (`n=6`), `p=0.28`. | **VAL-A only; not released for Stage B calibration.** This is an AE4-null whole-gland outcome, and a low-power null result is not evidence of exact equality. |
| E15-14 | NHE/acid-base compensation | 2015 pH protocol and Discussion | DM + DI | Stimulation-associated EIPA-sensitive pH response showed no detected genotype effect. That supports similar assay-level NHE-dependent alkalinization. The paper explicitly notes that pH experiments do **not** directly measure carbonic-anhydrase activity. | **LOC-B/CONTEXT**; do not encode “CA unchanged” as direct data |
| E15-15 | AE4 membrane domain in acini | 2015 Discussion | DI | AE4 antibodies failed specificity against null tissue. Basolateral localization was inferred from function; an apical contribution could not be excluded. | **CONTEXT**; no measured basolateral fraction |
| E15-16 | AE4-null cell volume and native stimulated Na/K trajectories | 2015 | ND | No direct AE4-null acinar-volume phenotype, and no native time-resolved intracellular Na or K trajectory, was reported. | No target exists; do not invent one for Stage B |

### 2016 transporter mechanism

| ID | Quantity or constraint | Exact source anchor | Class | Numerical value, uncertainty, and `n` | Use |
|---|---|---|---|---|---|
| E16-01 | Native AE4 is Na-dependent and HCO3/Cl coupled | 2016 Fig. 1 | DM + DI | In `Ae2-/-` native SMG acini, low external Cl caused alkalinization with Na (`n=8`), little response with NMDG (`n=7`), and little response in `Ae2-/-;Ae4-/-` acini despite Na (`n=8`); at least five acini/experiment and at least three preparations. | **CAL-T**; native-system anchor |
| E16-02 | Mouse and human AE4 Cl and HCO3 flux association | 2016 Fig. 2 | DM | Mouse Ae4 CHO: BCECF response with Na (`n=11`) but not Na-free (`n=6`) or HCO3-free (`n=8`); SPQ Cl loss with Na (`n=10`) but not HCO3-free (`n=7`) or Na-free (`n=6`). Human AE4 also Na dependent (`n=5`; mouse `n=8`; nontransfected `n=9`). | **CAL-T** |
| E16-03 | Direct Na signal associated with AE4 transport | 2016 Fig. 3 | DM | SBFI showed intracellular Na increase during outward-Cl-gradient transport in mouse-Ae4 CHO (`n=15`) but not nontransfected cells (`n=13`), and reversal after Cl restoration. | **CAL-T**. Direct ion-signal association; exact Na stoichiometry remains unidentified. |
| E16-04 | No detected electrogenicity/voltage dependence | 2016 Fig. 4 | DM + DI | Current-clamp transport-associated change `<2 mV`; Cl-exit `n=7` Ae4 versus `n=8` nontransfected. Alkalinization rates were essentially identical at `-100` and `0 mV` (`n=5` each). Current differences Ae4/nontransfected: `-1.68 +/- 0.44`/`-1.27 +/- 0.46 pA/pF` at -100 mV (`p=0.54`) and `-0.28 +/- 0.09`/`-0.21 +/- 0.14` at 0 mV (`p=0.69`). | **CAL-T**. The observations constrain net charge near zero at assay resolution; they do not measure unique stoichiometry. |
| E16-05 | Direct K signal associated with AE4 transport | 2016 Fig. 8C-D | DM | Under Na-free conditions, PBFI showed sustained K loss when inward Cl/outward K gradients were imposed in Ae4 CHO (`n=6`) but little signal in nontransfected cells (`n=8`). PBFI is only about 1.5-fold K selective over Na, hence the Na-free protocol. | **CAL-T**. Exact K stoichiometry remains unidentified. |
| E16-06 | Broad monovalent-cation permissiveness | 2016 Fig. 8A-B | DM | Robust alkalinization with Na (`n=10`), K (`n=8`), Li, Rb, and Cs (`n=6` each), but not NMDG (`n=6`). | **CAL-T/CONTEXT** |
| E16-07 | Na and K dose response | 2016 Fig. 9 and text | RF | Na: `EC50=49 mM`, Hill `nH=2.0`, `Rmin=0.3`, `Rmax=1.5`; K: `EC50=62 mM`, `nH=1.8`, `Rmin=0.4`, `Rmax=1.6`. `Rmin/Rmax` units are `10^-3 s^-1`; at least four experiments/concentration and at least three electroporations. No zero-cation point, fit SE/CI, or parameter covariance was reported. | **CAL-T** as effective affinities/cooperativity, not literal counts of microscopic sites |
| E16-08 | Alternative electroneutral stoichiometries | 2016 Fig. 7 and equations 1-3 | MH + RF | Paper evaluates `Cl:Na:HCO3` = `1:1:2`, `1:2:3`, and `2:1:3`. These are assumed candidate stoichiometries in thermodynamic calculations, not measured stoichiometry. | **CONTEXT**; all require independent thermodynamic testing |
| E16-09 | Na-only physiological-direction problem | 2016 Fig. 7 | RF + DI | Under the paper's fixed acinar concentrations, Na-only Cl influx requires, respectively, Cl_i `<7,<0.5,<26 mM`, Na_i `>104,>143,>54 mM`, or HCO3_i `>49,>84,>29 mM` for the three candidate stoichiometries. | **CAL-T/CONTEXT**; conditional on the assumed concentrations and stoichiometry |
| E16-10 | K/Na-pooled thermodynamic proposal | 2016 Fig. 10 | MH + RF | For assumed `1:1:2`, nearly equal Na/K handling, and `Na_i+K_i approximately Na_o+K_o`, calculated Cl uptake becomes favorable for Cl_i `<73 mM` or HCO3_i `>16 mM`. | **CONTEXT**, not a direct reversal measurement and not permission to pool Na+K in a state model |
| E16-11 | Physiological K recycling role | 2016 Discussion | MH | AE4 was proposed to import Cl while exporting K(Na) and HCO3 and possibly support Na/K-ATPase cation handling. | **CONTEXT**; physiological direction was not directly measured in stimulated native acini |
| E16-12 | Apparent Cl-free NBC-like activity requires an internal exchangeable anion | 2016 Figs. 5-6 and Table 4 | DM + DI | Whole-cell alkalinization rates with pipette O (15 mM Cl, no HCO3), P (4 mM Cl, no HCO3), and Q (4 mM Cl, 15 mM HCO3) were `0.12 +/- 0.03` (`n=6`), `0.02 +/- 0.01` (`n=6`), and `0.08 +/- 0.01` (`n=7`) `x10^-3 s^-1`; P differed from O (`p=0.01`) and Q (`p=0.001`). | **CAL-T/CONTEXT**. Supports electroneutral HCO3/Na-HCO3 exchange rather than net Na-HCO3 cotransport in this assay; OH or carbonate substitution was not excluded. |

### 2021 PKA regulation

| ID | Quantity or constraint | Exact source anchor | Class | Numerical value, uncertainty, and `n` | Use |
|---|---|---|---|---|---|
| E21-01 | Beta-adrenergic activation of total native exchanger activity | 2021 Fig. 1B | DM | WT native SMG acini: control `n=11` experiments versus 5 uM IPR `n=12`; significant increase (`p<0.001`), mean +/- SD. No numerical fold or fit uncertainty is printed in text. | **CAL-T/CONTEXT** |
| E21-02 | Beta/PKA-associated activation of native AE4 isolated genetically | 2021 Fig. 1C | DM + DI | `Ae2-/-` acini: control `n=12`, IPR `n=11`, IPR+10 uM H89 `n=7`; IPR increases activity and H89 prevents the increase (`p<0.0001` and `<0.05` comparisons), mean +/- SD. H89 is explicitly nonselective. | **CAL-T** for regulatory association; not proof of direct phosphorylation |
| E21-03 | Forskolin response in heterologous AE4 | 2021 Fig. 2 and Discussion | DM + RF | AE4 CHO: control `n=40`, forskolin `n=39`, forskolin+H89 `n=41`; source describes an **approximately 25%** alkalinization-rate increase, with H89 prevention (`p<0.0001`, `<0.01`), mean +/- SD. | **CAL-T** as an assay-specific effective change, not a universal turnover multiplier |
| E21-04 | PKAc-associated increase in both HCO3 and Cl readouts | 2021 Fig. 3 | DM | AE4+PKAc versus AE4+K73M dominant-negative PKAc: BCECF `n=18,16`; SPQ `n=8,9`; `p<0.005`, mean +/- SD. PKA constructs alone did not generate the flux. | **CAL-T** |
| E21-05 | S173-associated regulation | 2021 Figs. 5-6 | DM + DI | WT/S173A/S273A basal assays: `n=42,37,46`; both mutants functional and membrane localized but have altered basal activity. S173A+DN/PKAc `n=14,15`, no detected activation; S273A `n=12,12`, activation retained (`p<0.01`). | **CAL-T/CONTEXT**. S173 is necessary for the construct-level response under this assay. |
| E21-06 | Direct phosphorylation, microscopic regulatory transition, and activation kinetics | 2021 Discussion | ND + MH | The authors explicitly did not demonstrate direct phosphorylation of S173 and list increased turnover, active transporter number, or altered stoichiometry as alternatives. No PKA onset/decay time constant was measured. | **CONTEXT**; no defensible kinetic target exists |

### 2025 mutagenesis and structural hypotheses

| ID | Quantity or constraint | Exact source anchor | Class | Numerical value, uncertainty, and `n` | Use |
|---|---|---|---|---|---|
| E25-01 | WT human AE4 heterologous activity and cation dependence | 2025 Fig. 3 | DM | WT HEK `n=17` versus nontransfected `n=16`; Na-dependent condition `n=6` versus nontransfected `n=5`; mean +/- SEM. Direct readout is HCO3-linked alkalinization after Cl removal, not a direct Na or K concentration trace. | **CAL-T** |
| E25-02 | Single-mutant hierarchy | 2025 Fig. 4 | DM + RF | Versus WT `n=17`: T448I (`n=17`) and T448G (`n=10`) about 50% lower; D709A (`n=19`), S446A (`n=8`), and T756A (`n=9`) about 30% lower; I758G (`n=11`), I758R (`n=10`), D709N (`n=12`), S447A (`n=9`), T754A (`n=9`) no detected change. Mean +/- SEM; `p<0.001` where marked. | **CAL-T**, while recognizing expression-system and mutation-specific effects |
| E25-03 | T448/T756 double-mutant hierarchy in Na | 2025 Fig. 5 | DM | S446A-T448I (`n=18`) and D709A-T448I (`n=19`) remain about 50% below WT; T756A-T448I (`n=18`) falls to nontransfected activity. Membrane colocalization Mander coefficient `>0.8` and immunoblot abundance showed no detected variant differences. | **CAL-T** |
| E25-04 | Na/K functional nonidentity | 2025 Fig. 6 | DM + DI | WT, S446A, T448I, and T448G show comparable Na- and K-bath activity. D709A and T756A are about 30% lower in K than Na in the caption (Discussion calls the extra decrement about 20%). T756A-T448I has no detectable Na-dependent activity but substantial K-dependent activity comparable to the single mutants. NMDG abolishes WT response. | **CAL-T**. Freeze the qualitative rank reversal; the paper's 20/30% wording is approximate and no ratio CI is reported. |
| E25-05 | TM3-TM10 coordination site and residue contacts | 2025 homology model and 100-ns MD, Figs. 2,7,8 | MH | Model is based on outward-open rat NDCBE cryo-EM structure (PDB 7RTM), not an AE4 structure. Na/K-HCO3 pairs were manually placed. Selected stable replicas, not all initial trajectories, were analyzed. D709/T713 cation and G449/K879/T448/T756 bicarbonate contacts are simulation results. | **CONTEXT**, never direct structural calibration |
| E25-06 | Simulated HCO3 hydrogen-bond occupancies | 2025 Table 2 | MH + RF | Percent occupancies: WT-Na D709/G449/K879/T448/T756 `89.2/90.1/97.8/87.3/57.0`; double-mutant-Na `66.5/54.8/58.6`; WT-K G449/K879/T448/T756 `31.7/81.5/77.3/64.6`; double-mutant-K D709/G449/K879 `42.0/20.0/86.6`. | **CONTEXT**; simulation observables, not experiments |
| E25-07 | Bicarbonate versus carbonate and sequential Na/K cycle | 2025 Discussion | MH | `1:1:2 Cl:cation:HCO3` and `1:1:1 Cl:cation:CO3` both remain compatible. Proposed sequential variants include `1:1:1:1 Cl:Na:HCO3:K` and `1:1:1:2 Cl:Na:CO3:K`. No Cl-containing MD cycle or completed alternating-access trajectory was simulated. | **CONTEXT**; hypotheses requiring explicit thermodynamic tests |
| E25-08 | State ordering of T448 and T756 | 2025 Discussion | MH | Authors speculate that T448 and T756 may act at different stages. Mutant epistasis/nonidentity is measured; temporal state order is not. | **CONTEXT** |

### Independent cation-current architecture and historical model evidence

| ID | Quantity or constraint | Exact source anchor | Class | Numerical value, uncertainty, and `n` | Use |
|---|---|---|---|---|---|
| E18-01 | Functional apical K conductance | Almassy 2018 publisher abstract | DM + DI | Local Ca photolysis near the apical membrane triggered a significant K current in parotid acinar cells. The abstract supports functional apical K channels; it does not provide a conductance fraction, sample size, or bath. | **CONTEXT** |
| E18-02 | Na/K-pump membrane distribution | Almassy 2018 publisher abstract | DM + DI | Pump labeling appeared distributed around the plasma membrane, including the apical pole. The accessible source does not provide a quantitative apical fraction. | **CONTEXT** |
| E18-03 | Apical pump/K-channel fractions | Almassy 2018 mathematical supplement | MH | The supplement **chooses** 30:70 apical:basolateral Na/K pump and 40:60 apical:basolateral K-channel distributions. These are model allocations said to reflect observations, not measured fractions. | Historical/model sensitivity only; **not calibration data** |
| E18-04 | Historical AE4 knockout prediction | Vera-Siguenza 2018 Fig. 8 | Historical model output | Historical model reports around **24% reduction in salivary flow rate** for its AE4 knockout. | Never experimental calibration or validation target |

## Frozen digitization of the 2015 secretion time course

The paper supplies the exact integrated reduction but not a numeric table of
minute means. The values below were read from the 700 x 790-pixel PMC raster of
Fig. 1A by mapping symbol centers to the printed axes. They are suitable only as
a shape/envelope check. Conservative reading uncertainty is **+/-0.4 uL/min per
mean** (symbol overlap, line width, and pixel mapping), additional to the
plotted but not digitized SEM. They are not independent of E15-02 and must not
replace its exact integrated value.

| Minute of CCh+IPR | Ae4 WT, approx. uL/min | Ae4-null, approx. uL/min |
|---:|---:|---:|
| 1 | 9.5 | 9.9 |
| 2 | 10.1 | 9.8 |
| 3 | 9.6 | 8.2 |
| 4 | 9.8 | 7.9 |
| 5 | 9.4 | 6.6 |
| 6 | 9.9 | 5.3 |
| 7 | 10.1 | 4.8 |
| 8 | 9.9 | 4.0 |
| 9 | 9.9 | 4.4 |
| 10 | 9.3 | 4.1 |

The robust constraints are therefore: comparable onset, divergence after about
3 min, sustained late loss, and a 10-min integrated KO/WT ratio of about 0.65.
Pointwise fitting to the digitized means is prohibited; evaluate whether a
frozen prediction enters this envelope.

## Resolution of “30%” versus “35 +/- 4.7%”

These are **not two experimental endpoints**.

1. The primary 2015 Results sentence and Fig. 1A inset report **35 +/- 4.7%
   less total saliva accumulated over 10 min** in systemic Ae4-null glands than
   their littermate WT glands under 0.3 uM CCh + 5 uM IPR (`n=6` per genotype,
   mean +/- SEM). This is the frozen experimental number.
2. Vera-Siguenza et al. 2018 later described the same cited 2015 phenotype as
   an **“approximate 30% decrease in gland salivation.”** That is a rounded
   narrative restatement, not a separately measured quantity and not “KO flow
   equals 30% of control.” The primary value implies KO total near 65% of WT.
3. The same 2018 paper's **24% reduction** is its own model prediction in Fig.
   8, not experimental evidence.
4. The 2015 statement that late AE4-null flow is “less than half” WT is a
   pointwise late-time observation, not the integrated 10-min percentage.

Accordingly, validation must report both the predicted 10-min KO/WT integral
and the delayed time-course shape; it must not substitute a rounded 30% target.

## Explicit nonidentifiabilities and source-access limits

- No source above directly measures AE4 transport stoichiometry, a full
  alternating-access state sequence, microscopic binding/release rates,
  transporter surface density, native stimulated Na/K trajectories, or an
  AE4-null acinar-volume phenotype.
- The 2016 Na/K Hill fits identify effective assay-level dependence, not unique
  microscopic state counts; fit uncertainties/covariances were not reported.
- The 2021 data identify beta/cAMP/PKA association and an S173 requirement in a
  construct assay, but neither direct S173 phosphorylation nor an activation
  time constant.
- The 2025 mutation hierarchy identifies functional nonidentity of Na and K
  conditions. Residue coordination, state order, carbonate transport, and the
  proposed sequential Na/K cycles remain model/MD hypotheses.
- Full primary texts were available for 2015, 2016, 2021, and 2025. The Almassy
  2018 version-of-record main text was subscription restricted; only its
  publisher abstract and mathematical supplements were accessible. Therefore
  no unverified sample size, species, bath, or membrane fraction from that
  paper is treated as frozen experimental fact.
- Negative findings are recorded as “no detected difference,” not biological
  equality. System, temperature, sensor, background, and stimulation context
  travel with every number; cross-system absolute-rate pooling is prohibited.
