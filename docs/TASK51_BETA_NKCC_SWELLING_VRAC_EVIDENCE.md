# Task 51 source evidence: beta-NKCC -> swelling -> VRAC-like secretion

This document is binding source memory for the final mechanistic reconstruction attempt. It records the independent literature that motivates one predeclared causal model. It does **not** establish that the model succeeds in reproducing the AE4 phenotype.

## 1. Why Task 49 did not test the full literature mechanism

Task 49 tested a separate apical anion current with the predeclared gate

`G_V = g_V * beta * max(V_i / V_rest - 1, 0)`.

The inherited NKCC regulation was deliberately left unchanged. In AE4 knockout under IPR-only stimulation the inherited core therefore had no beta-responsive chemical source. Starting exactly at genotype rest, volume did not increase, the positive-swelling gate remained zero, and the VRAC-like branch could not self-start. The frozen result consequently predicted zero AE4-null IPR uptake, whereas Peña-Münzenmayer et al. 2015 report `0.20 +/- 0.03 x10^-3 s^-1`.

This rejects the **positive-swelling-only gate on the beta-blind parent core**. It does not reject the experimentally supported beta/cAMP secretory pathway.

## 2. Catalán et al. 2015: adult mouse SMG cAMP pathway

Primary source:

Marcelo A. Catalán et al., “A fluid secretion pathway unmasked by acinar-specific Tmem16A gene ablation in the adult mouse salivary gland,” *PNAS* 112(7), 2263-2268 (2015).

- DOI: `10.1073/pnas.1415739112`
- PMID: `25646474`
- PMCID: `PMC4343136`

Source-backed observations in adult mouse submandibular gland:

1. Muscarinic/Ca-dependent secretion requires TMEM16A in acinar cells.
2. Isoproterenol/beta-adrenergic cAMP-dependent secretion persists after acinar-specific TMEM16A deletion.
3. The IPR response is also not abolished by CFTR or ClC-2 loss.
4. DCPIB and NPPB markedly inhibit the IPR secretory response, supporting a volume-regulated-anion-conductance-like apical pathway.
5. IPR causes cell swelling rather than the shrinkage associated with muscarinic secretion. The paper reports approximately `12.5 +/- 0.2%` swelling; the Task 49 figure extraction retained the source discrepancy between this narrative value and the terminal graphical proxy rather than silently forcing them to agree.
6. The swelling and secretory response develop slowly on the scale of greater than roughly one minute; this is physiological context, not a hard timing target under `PHENOTYPE_TARGET_CONVENTION.md`.
7. IPR activates an outwardly rectifying chloride conductance whose reversal is close to the chloride equilibrium potential and whose current is partly DCPIB-sensitive.
8. Task 49 source metrology extracted an IPR-induced whole-cell current of approximately `208.6 +/- 35.2 pA` and a DCPIB-blocked component of approximately `107.8 +/- 17.2 pA` at the reported voltage condition. These are source measurements, not yet a directly identified model apical conductance; the blocked current need not be wholly apical and holding voltage is not by itself chloride driving force.
9. The authors interpret the beta-adrenergic pathway as involving a volume-regulated anion conductance and discuss beta-stimulated solute loading/cell swelling as the route that can activate it. The molecular identity of the apical beta-dependent channel was not established.

Important limitation: chronic TMEM16A deletion changes other protein abundances, including NKCC1/AQP5 in the source study. A model “TMEM16A off” experiment is therefore an idealised pathway-separation calculation, not an exact chronic knockout reproduction.

## 3. Independent salivary evidence for beta/cAMP activation of NKCC1

The beta-to-NKCC arm is independently supported and was not invented to repair Task 49.

### 3.1 Turner et al. 1992, rat parotid acini

“Beta-adrenergic upregulation of the Na+-K+-2Cl- cotransporter in rat parotid acinar cells.”

- DOI: `10.1172/JCI115695`
- PMID: `1313447`
- PMCID: `PMC442971`

The study reports that brief isoproterenol pretreatment increased the measured cotransporter-dependent recovery rate approximately **threefold**, with half-maximal isoproterenol concentration `K1/2 = 21.5 nM`. The response was beta1-adrenergic, prevented by protein-kinase inhibitors, and mimicked by cAMP-elevating interventions/cAMP analogues, supporting cAMP/PKA-dependent NKCC regulation.

Species/tissue qualification: rat parotid is not mouse SMG. The approximately threefold gain and 21.5 nM half-maximal concentration therefore cannot be copied as exact mouse-SMG constants. They establish sign, existence, approximate scale and a defensible cross-source prior/sensitivity range.

### 3.2 Tanimura et al. 1995, rat parotid acini

“Involvement of direct phosphorylation in the regulation of the rat parotid Na+-K+-2Cl- cotransporter.”

- DOI: `10.1074/jbc.270.42.25252`
- PMID: `7559664`

A brief isoproterenol exposure altered phosphorylation of the salivary NKCC protein. The half-maximal isoproterenol effect on phosphorylation was approximately `20 nM`, closely matching the earlier functional activation. The paper also reports an AlF4-induced cotransporter activation comparable with isoproterenol of approximately sixfold. This supports strong regulated capacity recruitment, but the “sixfold” value is not to be treated as a direct mouse-SMG numerical constant.

### 3.3 Kurihara et al. 1999 and 2002, rat parotid acini

1999:

“Characterization of a phosphorylation event resulting in upregulation of the salivary Na+-K+-2Cl- cotransporter.”

- DOI: `10.1152/ajpcell.1999.277.6.C1184`
- PMID: `10600770`

Beta-adrenergic stimulation increased NKCC activity, phosphorylation and the number of high-affinity bumetanide-binding sites, consistent with recruitment/activation of previously quiescent cotransporters rather than a different bumetanide sensitivity of already active transporters.

2002:

“Phosphorylation of the salivary Na+-K+-2Cl- cotransporter.”

- DOI: `10.1152/ajpcell.00352.2001`
- PMID: `11880270`

The study confirms that salivary NKCC1 activity is strongly upregulated by beta-adrenergic stimulation and that the upregulatory phosphorylation can be mimicked by cAMP in permeabilised acini and blocked by a PKA inhibitor. The authors conclude that PKA participates in the regulatory phosphorylation pathway, probably through additional cellular kinase/phosphatase machinery rather than simple direct phosphorylation in isolated membranes.

### 3.4 Rat submandibular gland evidence

“Activation of the Na+-K+(NH4+)-2Cl- cotransporter from rat submandibular glands in response to VIP.” PMID `9880083`.

In rat submandibular acinar cells, isoproterenol increased the measured cotransporter-mediated NH4+ influx approximately **2.5-fold**; forskolin mimicked the effect and inhibition of cAMP-dependent protein kinase blocked the response. This is especially useful as tissue-context evidence that beta/cAMP activation of NKCC is not restricted to parotid cells, although it remains rat rather than mouse.

## 4. Independent mouse-SMG evidence for beta/PKA activation of AE4

Gaspar Peña-Münzenmayer et al., “Activation of the Ae4 (Slc4a9) cation-driven Cl-/HCO3- exchanger by the cAMP-dependent protein kinase in salivary gland acinar cells,” *Am J Physiol Gastrointest Liver Physiol* 321:G628-G638 (2021).

- DOI: `10.1152/ajpgi.00145.2021`
- PMID: `34585968`
- PMCID: `PMC8887885`

In mouse SMG acinar cells, beta-adrenergic stimulation increases AE4-mediated exchanger activity. H89 prevents the activation; constitutively active PKA increases AE4 activity in heterologous expression; the S173A mutation abolishes PKA activation while S273A retains it. The present model already contains effective beta/PKA activation of AE4. Task 51 must not add this mechanism again.

## 5. Primary AE4 phenotype that remains held out from Task 51 calibration

Peña-Münzenmayer et al., JBC 2015, DOI `10.1074/jbc.M114.612895`, PMID `25745107`, PMCID `PMC4409235`:

- AE4 knockout reduces 10-minute CCh+IPR saliva by approximately `35 +/- 4.7%`.
- AE4-KO CCh-only chloride uptake is approximately preserved (`2.30 +/- 0.10 x10^-3 s^-1`).
- AE4-KO CCh+IPR uptake is strongly reduced (`0.90 +/- 0.09 x10^-3 s^-1`).
- AE4-KO IPR-only uptake remains positive (`0.20 +/- 0.03 x10^-3 s^-1`), approximately half the WT mean (`0.40 +/- 0.07`).

These AE4 phenotype values are **prediction/validation data for Task 51**, not calibration data for beta-NKCC or the VRAC-like conductance.

## 6. One literature-backed causal model permitted for Task 51

Task 51 is not a mechanism search. It tests exactly the following chain:

`beta/IPR -> increased NKCC1 activity -> additional Na/K/Cl solute loading -> cell swelling -> volume-sensitive apical anion conductance -> chloride/fluid secretion`

in parallel with the already represented

`beta/PKA -> AE4 activation`.

The intended causal interpretation is:

- In WT beta stimulation, both NKCC1 and AE4 can contribute to basolateral chloride loading while the beta-associated apical anion route exports chloride.
- In AE4 knockout, beta-NKCC can still start the IPR response and therefore remove Task 49’s exact self-start failure, but AE4 loading is absent. The same apical beta-secretory demand must therefore be sustained by the remaining loaders.
- If the remaining system cannot fully sustain that demand, intracellular chloride and/or secretory driving force fall and a substantial AE4-dependent secretion deficit can emerge **without an explicit AE4-expression multiplier on the apical conductance**.

This is a source-backed hypothesis, not an established result.

## 7. Parameter-identification rule

The Task 50 multiplier and the approximately 35% AE4 secretion deficit must not be used to choose Task 51 parameters.

Permitted independent information includes:

- direct beta-NKCC salivary studies above for sign and plausible activity scale;
- mouse-SMG IPR swelling magnitude from Catalán et al.;
- the IPR-induced/DCPIB-sensitive current-voltage information from Catalán et al. and its supplement, provided the conductance mapping is derived rather than guessed;
- source-backed mouse-SMG IPR secretion/channel-blockade data for validation of the beta branch, while preserving any acinus-to-gland scale limitation explicitly.

The rat beta-NKCC numbers are cross-species priors/sensitivity bounds, not exact mouse-SMG constants. Prefer mouse-SMG swelling/current data to identify the effective strength wherever possible.

If the independent data do not identify or tightly bound the beta-NKCC and VRAC parameters sufficiently for a held-out AE4 prediction, Task 51 must report that non-identifiability. It must not close the system with the AE4 35% target.

## 8. Benchmark relationship to Task 50

Task 50 is separately frozen in `docs/TASK50_EFFECTIVE_COUPLING_BENCHMARK.md` and on branch `archive/task-50-working-effective-coupling-benchmark` at commit `b16c30094b95f61a73d8f1977cd58e79c7bb50f6`.

Task 51 must **not** include the Task 50 multiplier in its equations. The benchmark exists to quantify the missing effective network influence that a mechanistic reconstruction should reproduce if the literature-backed mechanism is sufficient:

- AE4 5%: approximately `23.16%` cumulative deficit at 600 s;
- AE4 null: approximately `30.26%` cumulative deficit at 600 s;
- substantial deficit across the saved 60-600 s window;
- WT, REST, CCh-only null and AE2 parent nesting.

A Task 51 success would be stronger than Task 50 because the missing effective coupling would emerge from independently constrained physiology instead of a phenotype-calibrated AE4-dependent conductance multiplier.

A Task 51 failure would also be informative: it would show that even this fuller independently documented beta-secretory architecture is insufficient on the current conservation-explicit chassis, leaving Task 50 as the preserved proof-of-sufficiency benchmark for an additional unidentified network interaction.
