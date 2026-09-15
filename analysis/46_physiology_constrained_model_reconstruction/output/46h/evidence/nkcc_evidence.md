# NKCC evidence review for Task 46H00

Read only review. No scientific calculations, fitting, model changes or new hypotheses were performed.

## Pinned sources

The Task 46 starting checkout is `/workspace/scratch/692dd1f74fa6/ae4-salivary-transport-control` at the requested `9553d04442aeef5b0a8e0ba5f9df07f0986d7074`. Its `analysis/46_physiology_constrained_model_reconstruction/source_pins.json` identifies Task 43 as `547f113d9123ab1976774639a69483faf40cef34` and Task 44 as `e5fa9840bf96be3147c6118daee4942468c78f8b`.

The exact pinned Task 44 checkout remains at `/workspace/scratch/32c3c5dd39a9/task44`; HEAD was verified read only. Primary evidence details are in `analysis/44_full_system_mathematics/literature_comparison.tex`, especially its first 60 lines. That report records a primary 2015 HTML audit with SHA256 `72dd59caaede1c477c643f818926390bf48e2966692b4dbd47b752b4f6fbe98d`, and source citation DOI `10.1074/jbc.M114.612895`, PMC4409235. This review read the preserved audit, not a newly fetched primary article.

The exact pinned Task 43 report is available at `/workspace/scratch/32c3c5dd39a9/task43/analysis/43_parameter_provenance/report.tex`; HEAD was verified as `547f113d9123ab1976774639a69483faf40cef34`. An older clean_task43_final copy was initially inspected and then all used numerical classification claims were confirmed in this exact pinned copy.

## What the isolated NKCC assay permits

The preserved primary source audit says the 2015 assay used bicarbonate free solution, 30 micromolar ethoxyzolamide and 50 micromolar T16Ainh A01. Bumetanide identified the NKCC1 dependent uptake component. Its observable was initial chloride uptake; no detectable genotype difference was reported for AE4 or AE2 deletion.

This is a protocol specific independent discriminator. It does not measure the full 600 s physiological stimulation integral, nor identify maximal carrier abundance or a numerical upper bound on stimulated compensation. No numerical equivalence tolerance is provided in the reviewed preserved audit. Do not convert lack of statistical difference into exact WT/KO flux equality, an arbitrary percentage band, or a normal stimulation cap.

Faithful translation requires the bicarbonate free, carbonic anhydrase inhibited, CaCC inhibited assay context. If that context cannot be represented by the compact CBM, retain the observation as a qualitative holdout with explicit protocol limitations. A hypothesis constraint about physiological NKCC recruitment must be labelled an assumption, not the assay result.

## Normal physiological stimulation and related observations

The gland stimulation used 0.3 micromolar carbachol plus 5 micromolar isoproterenol. The isolated uptake assay must remain separate from this combined stimulation. The same study reports lower resting chloride and slower chloride reuptake after AE4 deletion, whereas initial chloride exit and resting pH did not detectably change. Stimulation induced NHE dependent alkalinisation did not show a genotype difference. Neither unchanged pH nor alkalinisation identifies an equality of NHE net fluxes at different ionic states.

The gland total secretion phenotype is a separate downstream holdout and should not constrain H00 network capacities. Existing historical reports have already exposed it, so it must not be advertised as globally unseen evidence. It can remain withheld from the present structural screen's construction and predeclared selection procedures.

## Quantities that are model assumptions, not new measurements

Task 43 report, sections on NKCC scale and NBC capacity: active NKCC scale `0.017334746894096052 fmol/s` is derived from an inherited model resting flux; it is not a measured salivary carrier amount. The resting NKCC cycle flux `0.1281222023866353 fmol/s` is likewise a model reference, not a direct assay measurement. The Palk coefficients are effective for the adopted fixed bath, so they must not silently carry over to bicarbonate free or different bath conditions.

The 70/30 WT chloride loading partition used for NBC construction is explicitly an assumed partition. Task 39 `analysis/39_palk_nkcc1_full_validation/final_answer.md` identifies 65 to 75 percent NKCC loading as an experimental context band, not a fitting or trajectory failure gate. A physiological uncertainty distribution or hard experimental range has not been justified for this partition. If the CBM freezes a partition band, label its chosen status explicitly.

The inherited algebraic NKCC recruitment endpoint is reached at 0.10 micromolar calcium, whereas NBC is fully recruited at 0.25 micromolar. These are transferred implementation settings, not observations of equal genotype capacity. Shared WT/KO capacities represent a declared acute perturbation construction assumption, while established knockout animals can have different starting chloride stores.

## Existing baseline, only as diagnostic context

Task 46 `BASELINE_REPRODUCTION.md` already distinguishes relative integrated NKCC rise (23.1634495075 percent), missing AE4 chloride replaced by NKCC (88.5070351033 percent), and local WT equilibrium compensation gain (about 92.901074 percent). None is the isolated assay observable. Do not repeat Task 44 or Task 46G. The baseline used the same WT starting state for all genotypes, and thus does not represent established knockout resting physiology.

## Suggested H00 ledger classification

* Construction: charge and elemental stoichiometry, acute AE4 removal, explicit direction conventions and shared capacity assumption.
* Independent protocol specific holdout: isolated NKCC uptake observation; report only its supported qualitative comparison until a source justified numerical uncertainty is available.
* Separate phenotype holdouts: secretion magnitude/time course, resting chloride reduction and reuptake phenotype; do not use secretion to select the CBM construction.
* Context or assumptions: old NKCC abundance scale, stimulation multiplier, 70/30 partition and existing model compensation metrics.

Tasks 46D, 46E, 46G, 46H and 46I were subsequently read in full from the orchestrator supplied instructions directory. Under 46H, the qualitative substantial deficit may be used as a structural discriminator; the exact magnitude and dynamic details remain validation material. No new scientific family is proposed here. For 46H04 the distinct assay requires a separate diagnostic, not transfer of exact flux equality to normal stimulation.
