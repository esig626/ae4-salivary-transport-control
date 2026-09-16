# Task 49 source seed: cAMP-dependent salivary chloride pathway

Every item here must be verified from the primary source before numerical use. This is a navigation aid, not a substitute for source extraction.

## S1. Catalan et al., PNAS 2015

**Citation**
Marcelo A. Catalan et al., `A fluid secretion pathway unmasked by acinar-specific Tmem16A gene ablation in the adult mouse salivary gland`, *PNAS* 112(7):2263-2268 (2015).

- DOI: `10.1073/pnas.1415739112`
- PMID: `25646474`
- PMCID: `PMC4343136`

**Directly relevant source conclusions to verify**

- Acinar-specific TMEM16A loss abolishes muscarinic/Ca2+-dependent fluid secretion.
- IPR/beta-adrenergic secretion remains comparable to control after TMEM16A deletion.
- IPR secretion is not materially impaired by CFTR or ClC-2 loss.
- IPR secretion is strongly reduced by the VRAC blockers DCPIB and NPPB.
- IPR secretion is strongly chloride-dependent and bicarbonate-dependent.
- IPR produces slowly developing acinar-cell swelling; the source links the swelling time course to the delayed IPR secretory response.
- IPR activates a DCPIB-sensitive chloride conductance with VRAC-like properties.

**Quantitative items to extract before model calibration**

- WT/control IPR-only secretion time course and cumulative amount.
- TMEM16A-KO IPR-only secretion time course and cumulative amount.
- DCPIB and NPPB IPR-only secretion effects.
- Low-Cl and bicarbonate-free IPR secretion effects.
- IPR-induced cell-volume trajectory/fractional swelling and time scale.
- IPR-induced DCPIB-sensitive current/conductance data.
- exact agonist concentrations, temperature, bath composition and timing.

Where values are only figures, digitise with a saved calibration record and uncertainty. Do not type an eyeballed value directly into a parameter file.

## S2. Pena-Munzenmayer et al., JBC 2015

**Citation**
Gaspar Pena-Munzenmayer et al., `Ae4 (Slc4a9) Anion Exchanger Drives Cl- Uptake-dependent Fluid Secretion by Mouse Submandibular Gland Acinar Cells`, *J Biol Chem* 290:10677-10688 (2015).

- DOI: `10.1074/jbc.M114.612895`
- PMID: `25745107`
- PMCID: `PMC4409235`

Task 48 already verified the complete ledger. Reuse it rather than re-entering numbers.

**Task 49 held-out AE4 features**

- combined CCh+IPR total AE4-null saliva `35 +/- 4.7%` lower over 10 min;
- first roughly 2-3 min comparable, sustained deficit thereafter;
- AE4 KO initial Cl uptake under CCh: `2.30 +/- 0.10 x 10^-3 s^-1`;
- AE4 KO under CCh+IPR: `0.90 +/- 0.09 x 10^-3 s^-1`;
- AE4 KO under IPR only: `0.20 +/- 0.03 x 10^-3 s^-1`;
- corresponding WT, AE2, resting, NHE and bicarbonate constraints are already in Task 48.

Do not use the 35% secretion loss or AE4-null combined-stimulus uptake to set Task 49 VRAC or beta-NKCC parameters.

## S3. Pena-Munzenmayer et al., AJP GI 2021

**Citation**
`Activation of the Ae4 (Slc4a9) Cation-driven Cl-/HCO3- Exchanger by the cAMP-dependent Protein Kinase in Salivary Gland Acinar Cells`.

- DOI: `10.1152/ajpgi.00145.2021`
- PMCID: `PMC8887885`

**Relevant constraints**

- AE4 has beta/cAMP/PKA-dependent activation; the current model already contains an effective AE4 beta state.
- The paper explicitly distinguishes Ca-dependent TMEM16A secretion from cAMP-dependent secretion and cites the VRAC-like cAMP apical pathway.
- It also notes NKCC1 as a likely basolateral Cl loader during cAMP-dependent stimulation because beta-adrenergic activation of salivary NKCC1 has been reported.

Do not duplicate the AE4 PKA mechanism.

## S4. Beta/cAMP regulation of salivary NKCC1

A source to verify:

Kinji Kurihara et al., `Phosphorylation of the salivary Na+-K+-2Cl- cotransporter`, *Am J Physiol Cell Physiol* 282:C817-C823 (2002).

- DOI: `10.1152/ajpcell.00352.2001`
- PMID: `11880270`

The paper reports strong beta-adrenergic/cAMP-dependent upregulation/phosphorylation of NKCC1 in rat parotid acinar cells and evidence for PKA involvement. This establishes a plausible positive beta input to NKCC1 but is not automatically a quantitative mouse-SMG gain.

The evidence worker should inspect the original beta-NKCC studies cited by S1-S3 and determine whether a mouse-SMG numerical constraint exists. If not, freeze only a source-supported sign/envelope and state the cross-species limitation.

## S5. Resting bicarbonate transport issue deliberately reserved

Task 36's NBC is explicitly a **stimulus-recruited increment**, with zero activation at rest chosen for exact nesting. That design says it is not a claim that basal NBC transport is absent.

Independent literature exists for salivary bicarbonate transport, including mouse SMG studies, but Task 49 must not alter basal NBC because that would mix the chronic-rest reconstruction with the cAMP-branch test.

Keep this as a declared unresolved issue. If Task 49 fixes the stimulated phenotype but the Task 48 resting Cl/pH failure remains, that is a clean reason for a separate later reconstruction rather than permission to change NBC inside Task 49.

## Predeclared mechanistic interpretation

The source-backed hypothesis to test is:

1. CCh/Ca activates the existing TMEM16A branch.
2. Beta/cAMP activates AE4 and NKCC1-mediated loading.
3. Beta-driven solute accumulation causes slow swelling.
4. Swelling recruits a distinct VRAC-like apical Cl conductance.
5. Under CCh+IPR, the measured intracellular Cl reuptake slope reflects **net** loading minus both apical efflux routes, not NKCC/AE4 flux alone.
6. In AE4 KO, the beta/cAMP efflux branch remains but AE4 loading is absent, so stored chloride and remaining loading may sustain early flow but fail to sustain the later combined response.

This is a prospective hypothesis fixed before Task 49 numerical execution. It must be rejected if the independent source-calibrated implementation does not support it.
