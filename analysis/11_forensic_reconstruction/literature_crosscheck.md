# Public-literature cross-check

## Scope and evidence boundary

This report records what can be established from public, primary sources independently of the historical files under `archive/legacy-2017/`. The public search was last repeated on **2026-08-27**. Historical-code observations are included only in the final comparison table and are explicitly labelled; they are not treated as published evidence.

The main public records are the [2018 AE4 article](https://doi.org/10.1007/s11538-017-0370-6) ([open full text, PMCID PMC5792321](https://pmc.ncbi.nlm.nih.gov/articles/PMC5792321/)), [Palk et al. 2010](https://doi.org/10.1016/j.jtbi.2010.06.027) ([open full text, PMCID PMC2954280](https://pmc.ncbi.nlm.nih.gov/articles/PMC2954280/)), and the [Palk corrigendum](https://doi.org/10.1016/j.jtbi.2012.10.027). Later model-lineage sources were checked at their publisher or full-text repository pages.

## Bottom line

1. Public evidence resolves two important conventions. Palk's volume balance is $d\omega_i/dt=q_b-q_a$, and the 2013 corrigendum supplies the NKCC coefficients in an M-based fourth-order convention. The equivalent mM coefficients are $a_2=2.0096\times10^{-5}$ and $a_4=1.3852\times10^{-6}$, as printed again in a later primary supplement.
2. Public evidence does **not** resolve the 2018 model's water recalibration, current/particle-flux normalization, acid/base closure, AE4 activity units, calcium input, or solver-time scaling. Later papers use related but changed models; they cannot be silently read back as corrections to the 2018 paper.
3. No correction, electronic supplement, or public source-code release attached to the 2018 article was located. This is a documented negative search, not proof that no private or unindexed implementation exists.
4. Later papers establish a clear scientific lineage from the 2018 secretion model to the 2019 single-cell, 2020 multicellular, 2021 in-vivo-calcium, and 2022 duct models. They also show that equations and parameterizations evolved. They support common ancestry, not identity of the archived MATLAB snapshot with the publication-generating code.

## Palk 2010 and its corrigendum

### What the 2010 article establishes

[Palk et al. 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2954280/) is the explicit source for much of the secretion-model scaffold inherited by the 2018 article.

- Equation 25 gives the cell-volume balance as
  \[
  \frac{d\omega_i}{dt}=q_b-q_a.
  \]
  Its model summary repeats the same sign as Eq. 48. This is both the conservation-consistent sign and the sign used by the historical seven-state function.
- The paper defines total luminal outflow as $q_{\mathrm{tot}}=q_a+q_{\mathrm{tight}}$; thus $q_a$ is cell-to-lumen water flow, whereas $q_b$ is interstitium-to-cell water flow.
- Palk states that the water permeabilities were fitted, and assumes the basolateral surface area to be 12.32 times the apical area, giving $L_{Pb}=12.32L_{Pa}$. Its Table 8 values are $L_{Pa}=1.68\times10^{-15}$, $L_{Pb}=2.07\times10^{-14}$, and $L_{Pt}=8.4\times10^{-17}\ \mathrm{L^2\,J^{-1}\,s^{-1}}$. Multiplication by $RT$ gives, to published rounding, the three 2018 permeability values reported in $\mathrm{L^2\,mol^{-1}\,s^{-1}}$.
- Appendix E constructs the two-state NKCC reduction from the Benjamin--Johnson transporter, while Appendix F similarly reduces the Smith--Crampin Na/K pump. The uncorrected Appendix E numerical form is dimensionally ambiguous when read without the later correction.

### What the 2013 corrigendum changes

The primary publisher record is [Palk et al., “Corrigendum to ‘A dynamic model of saliva secretion’,” *Journal of Theoretical Biology* 317:428](https://doi.org/10.1016/j.jtbi.2012.10.027). PubMed also flags an erratum against the [2010 record](https://pubmed.ncbi.nlm.nih.gov/20600135/). The publisher index exposes the corrected expression and coefficients, although the full ScienceDirect body returned HTTP 403 during this audit. The equivalent mM values were independently cross-checked against the 2019 publisher supplement. The corrected simplified NKCC expression is

\[
\nu_{\mathrm{NKCC}}=
\frac{\alpha_1-\alpha_2 K_iN_iL_i^2}
     {\alpha_3+\alpha_4 K_iN_iL_i^2},
\]

with $\alpha_1=157.55\ \mathrm{s^{-2}}$, $\alpha_2=2.0096\times10^7\ \mathrm{M^{-4}\,s^{-2}}$, $\alpha_3=1.0306\ \mathrm{s^{-1}}$, and $\alpha_4=1.3852\times10^6\ \mathrm{M^{-4}\,s^{-1}}$. For concentrations supplied in mM, the fourth-order coefficients must be multiplied by $10^{-12}$:

| Coefficient | Corrigendum, concentrations in M | Equivalent, concentrations in mM |
|---|---:|---:|
| $a_2$ | $2.0096\times10^7\ \mathrm{M^{-4}\,s^{-2}}$ | $2.0096\times10^{-5}\ \mathrm{mM^{-4}\,s^{-2}}$ |
| $a_4$ | $1.3852\times10^6\ \mathrm{M^{-4}\,s^{-1}}$ | $1.3852\times10^{-6}\ \mathrm{mM^{-4}\,s^{-1}}$ |

The [2018 article](https://pmc.ncbi.nlm.nih.gov/articles/PMC5792321/) prints the large numerical coefficients while labelling them with mM-based units. That combination is not the corrected convention. The historical code's large coefficients combined with an explicit $10^{-3}$ conversion of each concentration, and its later small-coefficient mM form, are algebraically equivalent to the corrigendum. This is a specific externally supported code convention; it does not authenticate the rest of the historical implementation.

## Upstream equation provenance

The 2018 article explicitly attributes its constitutive laws as follows. This table records source lineage, not a claim that every 2018 transcription is exact.

| 2018 component | Published source trail | What the trail supports |
|---|---|---|
| NKCC1 | [Benjamin and Johnson 1997](https://doi.org/10.1152/ajprenal.1997.273.3.F473), simplified by [Palk et al. 2010](https://doi.org/10.1016/j.jtbi.2010.06.027), then corrected in [Palk et al. 2013](https://doi.org/10.1016/j.jtbi.2012.10.027) | Two-state rational flux form and corrected M-to-mM coefficient conversion |
| Na/K ATPase | [Smith and Crampin 2004](https://doi.org/10.1016/j.pbiomolbio.2004.01.010), reduced by [Palk et al. 2010](https://doi.org/10.1016/j.jtbi.2010.06.027) | Reduced $K_e^2Na_i^3/(K_e^2+\alpha Na_i^3)$ dependence; it does not by itself determine the 2018 density normalization |
| AE2 and NHE1 | [Falkenberg and Jakobsson 2010](https://doi.org/10.1016/j.bpj.2009.11.045) | Saturating exchanger/antiporter source family cited by the 2018 appendices |
| Ca-activated K channel | [Takahata et al. 2003](https://doi.org/10.1152/ajpcell.00250.2002) | Experimental/current-model basis cited for the CaKC open-probability law |
| CO2/bicarbonate subsystem | [Sharp, Crampin, and Sneyd 2015](https://doi.org/10.1016/j.jtbi.2015.06.050) | Source cited for buffer and CO2 transport constants |
| AE4 stoichiometry | [Peña-Münzenmayer et al. 2016](https://doi.org/10.1085/jgp.201611571) | Experimental support that AE4 is a monovalent-cation-dependent Cl-/HCO3- exchanger; it does not select the paper's probabilistic Na/K partition rule or the historical code's sodium-only rule |

Three distinctions matter. First, a cited source for a functional form does not resolve the downstream paper's fitted activity or area/volume normalization. Second, the Palk corrigendum repairs NKCC units only. Third, the experimental AE4 paper establishes transport stoichiometry and cation dependence, but not the exact fourth-order mass-action implementation or how a mixed Na/K cycle should be split between the Na and K balances.

## The 2018 public record

The [Springer version of record](https://link.springer.com/article/10.1007/s11538-017-0370-6) is labelled “Original Article,” published online 5 December 2017 and assigned to volume 80 (2018). The [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC5792321/) supplies the complete equations, tables, and figures.

The public paper specifies, among other items:

- cell Na and K balances that allocate AE4 cation efflux using $Na_i/(Na_i+K_i)$ and $K_i/(Na_i+K_i)$, respectively;
- a squared intracellular-H gate in the NHE1 flux;
- resting $HCO_{3,i}=12.1$ mM, $HCO_{3,e}=42.9$ mM, and $CO_{2,i}=6.6$ mM;
- $K_{CaCC}=0.26\ \mu\mathrm{M}$, $K_{CaKC}=0.182\ \mu\mathrm{M}$, $G_{Ae2}=0.01807$ fmol/s, $G_{NHE1}=0.0305$ fmol/s, and $G_{Ae4}=0.66\ \mathrm{amol}/\mu\mathrm{m}^3$;
- an approximately 24% reduction in steady secretion for AE4 knockout and negligible change for AE2 knockout.

The public record does not supply executable code, a machine-readable parameter bundle, or a complete numerical protocol sufficient to resolve all ambiguities identified in the earlier audit.

### Search for a 2018 correction, supplement, or public code

The following were searched by DOI, exact title, and distinctive numerical constants: the Springer article and correction/Crossmark records; PMC and PubMed; publisher searches for “correction,” “erratum,” “corrigendum,” and “electronic supplementary material”; GitHub and general web-index queries for the title, DOI, and constants `20468119.8347071` and `18.9502607221547`.

No 2018 correction, attached electronic supplementary material, source-code repository, or downloadable model implementation was located. In contrast, the 2019 paper's publisher page has an explicit “Electronic supplementary material” section. The proper conclusion is **“not located in the indexed public sources searched as of 2026-08-27,”** not “does not exist.” The authors or a non-indexed private repository could still hold a prepublication implementation.

## Later published reuse and model evolution

| Publication | Direct public statement or equation evidence | Forensic implication |
|---|---|---|
| [Vera-Sigüenza et al. 2019, *A Mathematical Model of Fluid Transport in an Accurate Reconstruction of Parotid Acinar Cells*](https://doi.org/10.1007/s11538-018-0534-z) | The abstract says it combines a new Ca model with the existing 2018 secretion model. Its [publisher supplement](https://static-content.springer.com/esm/art:10.1007%2Fs11538-018-0534-z/MediaObjects/11538_2018_534_MOESM2_ESM.pdf) gives $d\omega_i/dt=J_b^w-J_a^w$, surface-area-scaled transporter fluxes, the corrected mM NKCC coefficients, and an AE4 mass-action form related to 2018. The main article's Appendix Eq. 21 instead prints a one-way saturating AE4 law. | Confirms explicit reuse and independently confirms the volume sign and NKCC conversion. Changed densities, geometry, water parameters, and internally differing AE4 forms mean it is not a corrigendum for the 2018 executable model. |
| [Vera-Sigüenza et al. 2020, *A Multicellular Model of Primary Saliva Secretion in the Parotid Gland*](https://doi.org/10.1007/s11538-020-00712-3) | Each cell uses the anatomically reconstructed secretion model. The appendices retain corrected mM NKCC coefficients and surface-area scaling, but use a different saturating bidirectional AE4 expression and later parameter values. | Establishes a continuing model family, with substantive constitutive evolution. |
| [Takano et al. 2021, *Highly localized intracellular Ca2+ signals promote optimal salivary gland fluid secretion*](https://doi.org/10.7554/eLife.66170) ([open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC8352588/)) | Methods state that parameters/equations remain those of the 2020 model except for enumerated Ca, channel-distribution, apical Na/K-pump, and PLC changes. They also state that anion-exchanger parameters came directly from the 2018 modelling work. | Strong lineage evidence, but the reference model is already the 2020 descendant, not the 2018 publication. The enumerated changes preclude exact identity. |
| [Su et al. 2022, *A Mathematical Model of Salivary Gland Duct Cells*](https://doi.org/10.1007/s11538-022-01041-3) | The paper says older acinar models are used with minor adjustments. It modifies the 2021 acinar model to add luminal HCO3-/H+, changes impermeant-solute concentration, adjusts parameters, and takes primary-saliva input from cell 4 of the 2020 acinus. | Downstream reuse is explicit, but it adds and refits components; it supplies no clarification of the disputed 2018 conventions. |

Later lineage therefore has an asymmetric evidential value: it strongly supports that the authors continued to reuse and revise the 2018 model family, and it independently validates the NKCC conversion and water-balance sign, but it weakens any attempt to infer a single immutable implementation from later equations.

## Published evidence versus historical implementation evidence

| Issue | External published evidence | Historical implementation evidence | What may be concluded |
|---|---|---|---|
| Cell-volume sign | Palk 2010 Eq. 25 and the 2019 supplement use inflow minus outflow | Seven-state code uses `qb-qa` | **Resolved publicly:** $d\omega_i/dt=q_b-q_a$ |
| NKCC fourth-order units | Palk corrigendum gives M-based coefficients; 2019/2020 give equivalent mM coefficients | `Original` multiplies each concentration by $10^{-3}$; outer code stores the small mM coefficients | **Resolved publicly:** the code convention is the corrected one |
| Water permeabilities | Palk and 2018 preserve a much larger basolateral than apical coefficient | Historical calibration multiplies the apical coefficient by 74.4 after forming the original ratios, leaving $b_2/b_1=0.0941$ | **Not a public correction:** this is an undocumented historical calibration and reverses the published ordering |
| AE4 cation allocation | 2018 explicitly partitions flux between Na and K; later papers use changed expressions | Historical dynamic code subtracts the full AE4 flux from Na and none from K | **Substantive model difference**, not resolved by later literature |
| NHE1 H dependence | 2018 and the 2019 public appendix use a squared H gate | Historical code uses first power | **Substantive model difference** |
| AE4 activity units | 2018 prose/table mix total, density, and per-volume language; 2019/2020 move to surface-flux conventions | Historical $g_4/H_{i0}$ numerically recovers 0.66 | Strong lineage clue, but **unit normalization remains incompletely published** |
| Acid/base bath and baseline | 2018 tables report 12.1 mM intracellular and 42.9 mM external bicarbonate | Historical baseline uses 10 mM, assigns external 40 mM but hard-codes 21 mM in exchanger terms | **Not resolved**; multiple historical/public values coexist |
| Calcium and time | 2018 plots a finite stimulation protocol; later models introduce spatial Ca and PLC systems | Historical reduced code uses a 0.05-to-0.55 step at code time 100 and hidden RHS factors of 10,000 or 1,000 | **Not resolved** by public sources |
| Electrical calibration | Public paper gives conductances and target voltages but no generating calibration script | Historical driver solves conductances from the target row | Historical provenance evidence only; no public confirmation of the normalization |

## Consequence for reconstruction

The literature cross-check validates applying the conservation volume sign and the corrected NKCC concentration convention in any independent reconstruction. It does not justify importing the historical water multiplier, sodium-only AE4 balance, first-power NHE1 gate, acid/base values, channel calibration, calcium step, or time multiplier as if they were published corrections.

The later papers should be used as lineage and consistency checks, not as retrospective replacements for 2018. Without a publication-generating source or a successful unchanged reproduction of the 2018 baseline and knockout phenotype, the public record does not raise the archived snapshot to `VERY LIKELY` identity and does not open the conditional clean-reconstruction gate.
