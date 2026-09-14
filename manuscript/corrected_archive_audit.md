# Corrected 2018 archive audit

## Source archive

Archive reviewed: `Ae4_Article_Corrected(2).zip`

Primary source files used in this audit:

- `Saliva.tex`, the corrected source of the 2018 article;
- `mybib.bib`, the bibliography shipped with that source;
- `Letter of Reply/Reply Letter2.tex`, used to recover clarifications made during review.

The compiled paper, figures and model files in the archive were inspected as supporting material where useful. The corrected TeX and reviewer reply were treated as the primary record of what the original model actually claimed.

## Scientific points recovered or clarified

1. The 2018 model predicted an approximately 24% reduction in secretion after AE4 deletion. The present manuscript already stated this correctly, so no change to that number was required.
2. The original model reproduced the direction of the experimental timing but reached its lower AE4-null flow after about 2 min, compared with roughly 10 min experimentally.
3. The original model predicted only minimal NKCC1 activation after AE2 or AE4 deletion and explicitly noted agreement with the 2015 experimental NKCC1 assay.
4. Peña-Münzenmayer et al. 2015 isolated NKCC1-dependent chloride uptake under bicarbonate-free conditions and detected no genotype-dependent change in NKCC1 functional activity. The same paper found no genotype-dependent change in stimulation-induced NHE-dependent alkalinisation or resting intracellular pH.
5. The reconstructed Task 40 model instead produces strong integrated NKCC1 compensation during the full 600 s stimulated trajectory. These are not the same observable or protocol, so they are not a formal contradiction. They are nevertheless an experimental tension and the compensation must be described as a model prediction requiring direct testing.
6. The corrected 2018 source confirms that its acid-base architecture contained CO2 hydration, NHE1, AE2 and AE4 but no other proton or bicarbonate transport mechanism. This supports the present manuscript's interpretation of the missing alkalinity supply.
7. The original model held lumen volume constant, with total luminal water inflow equal to outflow. The reconstruction instead evolves a finite lumen with explicit outflow.
8. The old paper and reviewer response make the compartmental scope clearer than the first manuscript draft did: the model concerns primary acinar secretion. Ductal modification of final saliva is outside the model even though whole-gland flow is used as the experimental phenotype.
9. The reviewer response records resting osmolarities of 292.6 mM in the bath, 296.6 mM in the cell and 297.4 mM in the acinar lumen, and clarifies that transmembrane and paracellular water fluxes are water movement while luminal outflow carries dissolved solute. The reconstructed amount balances already make the latter distinction explicit.

## Manuscript changes made from the audit

### `main.tex`

- Added the two-stage acinar then ductal framing of salivary secretion.
- Made explicit that the model addresses primary acinar secretion rather than final duct-modified saliva.
- Restored the historical evidence for NKCC1 as the dominant chloride-loading pathway and the motivation for a parallel bicarbonate-dependent uptake route.
- Added the original model's prediction of only minimal NKCC1 compensation as an important contrast with the reconstruction.

### `discussion.tex`

- Added a direct comparison between the Task 40 NKCC1 compensation and the 2015 isolated NKCC1 and NHE assays.
- Stated explicitly that strong NKCC1 compensation is model-dependent and not experimentally established.
- Added a dedicated `Scope and limitations` subsection covering the single well-stirred acinar cell, fixed bath and finite lumen, absence of ductal transport, prescribed calcium input, simplified beta regulation, imported NHE1 kinetics, surrogate AE2 kinetics, structural rather than isoform-specific NBC, the fixed inherited kinetics of Task 42, the target-selected status of Task 41 and the finite 600 s validation horizon.

### `claims_and_evidence.md`

- Added historical constraints recovered from the corrected 2018 source.
- Added the 2015 NKCC1, NHE and resting-pH observations.
- Added explicit prohibitions against presenting modelled NKCC1 compensation as experimentally established or claiming the model reproduces final duct-modified saliva.

### `references.bib`

The current manuscript bibliography and `mybib.bib` were consolidated to 56 unique entries. Seven duplicate sources were collapsed onto the manuscript's existing keys:

- `silva1977mechanism` -> `Silva1977`
- `melvin2005regulation` -> `Melvin2005`
- `gin2007mathematical` -> `Gin2007`
- `benjamin1997quantitative` -> `Benjamin1997`
- `palk2010dynamic` -> `Palk2010`
- `pena2015ae4` -> `Pena2015`
- `pena2016Ae4` -> `Pena2016`

The `Pena2015` author list was corrected to Yasna Jaramillo and Frances Liu. The garbled author list for the 1973 Imai paper was corrected to Yusuke Imai, Hiroyasu Nishikawa, Kazuo Yoshizaki and Hiroshi Watari.

Legacy DOI fields that were incorrect were replaced rather than propagated. Corrections include Nauntofte 1992, Nguyen et al. 2004, Roussa et al. 2001, Wang et al. 2003 and Martinez and Cassity 1983. Verified missing DOI fields were also added for selected entries including Delporte and Steinfeld 2006, Mangos et al. 1973, Foskett and Melvin 1989, Bruce et al. 2002, Moldover et al. 1988, Dupont et al. 2016, Young 1968, Thaysen et al. 1954 and Benjamin and Johnson 1997.

## Files deliberately not changed

No model source, frozen trajectory, numerical checkpoint, figure source value or simulation result was modified during this audit. `model.tex` was inspected but not rewritten because the finite lumen, explicit luminal solute outflow, prescribed calcium stimulus, structural NBC interpretation and current closure are already stated there; the missing interpretive boundaries are now made explicit in the Discussion.

## Bottom line

The corrected archive strengthens the manuscript's central argument but also sharpens one important caveat. The 2018 model and the original experiment both supported little or no compensatory NKCC1 up-regulation under the conditions they examined, whereas the reconstruction relies on strong integrated NKCC1 compensation to erase most of the AE4-null secretion deficit. That discrepancy is now treated as a testable unresolved point rather than being hidden inside the newer model narrative.
