# Task 48 primary evidence advisory

Scientific starting SHA: `5c945c4a7857657b5f7c76d5bf9e3d3abc7d5721`. No model, parameter or branch changes.

All 28 Table 1 group means, SEs and sample sizes agree exactly with the seed. All primary 2015 numerical seed observations are confirmed against the full paper. `verified_constraints.json` preserves the values and adds source locators, units, verification status and caveats.

## Corrections and important qualifications

* Correct the 2016 JGP DOI to `10.1085/jgp.201611571`, PMCID `PMC4845690`, PMID `27114614`. The seed DOI is incorrect.
* The PKA paper is from 2021; its figure errors are SD, whereas all JBC 2015 results use SE.
* In JBC 2015, n denotes experiments from preparations using at least three mice per condition. It is not necessarily the independent animal count.
* Primary controls differ by strain: AE4 BS/129Svj and AE2 C57Bl/6. The paper explicitly reports different initial secretion kinetics between strains, and compares each knockout with its own control. Do not force both measured control curves to be identical.
* NHE and NKCC comparisons lack tabulated means. Keep them qualitative until a documented figure digitisation exists. No significance result authorises an arbitrary numerical band.
* SPQ uses F0/F. The reported exit is an amplitude change; uptake is a regression slope. BCECF uses the normalised 490/440 ratio F/F0, so the stimulated NHE assay is not a direct molar NHE flux.
* Reported bicarbonate amounts and uncertainties are confirmed, but the Methods phrase dividing concentration by saliva amount conflicts dimensionally with microequivalents. Preserve reported numbers and flag the inconsistency; do not silently invent a raw measurement correction.
* The primary source states output after about three minutes is below half WT, whilst the first two to three minutes are comparable. These are textual shape constraints until digitisation.
* Na free assays replace Na salts by NMDG/choline salts. They do not remove extracellular charge or osmoles.
* The 2016 Na/K Hill fits measure external cation dependence during reverse exchange in heterologous cells, not physiological realised routing or whole cycle stoichiometry. No uncertainties for fitted EC50/Hill coefficients are tabulated.
* The PKA activation amplitude is supported; no activation time constant is identified by the supplied means. H89 is nonselective and phosphorylation of S173 is proposed, while S173 dependence is observed.

## Source access and provenance

Full JBC HTML and JGP JATS were fetched through public primary repositories. Five unmodified JBC figures are available in this advisory directory for protocol inspection. No figure points have been digitised. No supplementary dataset link was identified in the fetched JBC HTML. PKA mechanisms were checked against the primary PubMed abstract and figure captions; its full PMC HTML returned a browser challenge.

| Source | Locator | Retrieval status |
|---|---|---|
| [JBC 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4409235/) | Table 1; Figures 1 to 5; Methods; Results | Full text and figures checked |
| [JGP 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4845690/) | Figures 2 to 4, 8, 9 | Full JATS checked |
| [PKA 2021](https://pubmed.ncbi.nlm.nih.gov/34585968/) | Abstract; Figures 1, 2, 3, 6 | All seed mechanism claims checked |
| [Evans 2000](https://pubmed.ncbi.nlm.nih.gov/10831596/) | Abstract | All seed claims checked |
| [Park 2001](https://pubmed.ncbi.nlm.nih.gov/11358967/) | Abstract | See JSON retrieval status |

## Table 1 verification

| Observable | AE4 WT mean ± SE (n) | AE4 KO | AE2 control | AE2 KO |
|---|---|---|---|---|
| rest_Cl_i | 50.1 ± 1.5 (7) | 36.5 ± 1.6 (6) | 53.4 ± 1.8 (6) | 54.5 ± 1.8 (6) |
| rest_pH_i | 6.91 ± 0.07 (4) | 6.89 ± 0.02 (4) | 6.87 ± 0.01 (4) | 6.95 ± 0.05 (6) |
| initial_Cl_exit_CCh_IPR | -0.21 ± 0.01 (10) | -0.18 ± 0.02 (9) | -0.16 ± 0.02 (13) | -0.15 ± 0.01 (18) |
| initial_Cl_exit_CCh | -0.18 ± 0.03 (7) | -0.15 ± 0.01 (6) | -0.14 ± 0.02 (8) | -0.17 ± 0.02 (9) |
| initial_Cl_uptake_CCh_IPR | 2.02 ± 0.1 (10) | 0.9 ± 0.09 (9) | 2.11 ± 0.5 (13) | 2.39 ± 0.2 (18) |
| initial_Cl_uptake_CCh | 2.18 ± 0.2 (7) | 2.3 ± 0.1 (6) | 2.2 ± 0.2 (8) | 2.3 ± 0.2 (9) |
| initial_Cl_uptake_IPR | 0.4 ± 0.07 (4) | 0.2 ± 0.03 (7) | 0.7 ± 0.1 (7) | 0.58 ± 0.05 (6) |

Units are mM for resting chloride, pH units for resting pH, dimensionless SPQ changes for exit, and 10⁻³ s⁻¹ for uptake.

## Source checksums

* `jbc.html` SHA256 `6a6278db90147a68fe86af03ad977a3a63899c21ee430bb12ed7b20f528db543`
* `jgp.xml` SHA256 `6b13788cac9a66edc331a183f7589a17cc4ac69f1c76fa3dfc27dbec8804b3e0`


Orchestrator integration: reviewed and accepted for the 48A evidence freeze. This document supplies evidence and restrictions, not a completed fit.
