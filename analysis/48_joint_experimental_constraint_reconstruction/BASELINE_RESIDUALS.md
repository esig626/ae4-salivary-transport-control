# Inherited architecture residuals

The inherited architecture is not a defensible joint reconstruction under a common source faithful observation rule. This decision has a parameter independent component; it is not based on failure of a restricted optimiser.

## Structural gate

AE4 KO CCh and CCh plus IPR have identical core dynamics for every shared admissible parameter vector. Beta input acts only on the deleted AE4 contribution. The observed uptake means are 2.30 ± 0.10 and 0.90 ± 0.09 (10⁻³ s⁻¹). The minimum two row Gaussian mean objective over any common prediction is 108.28729282. This is a relaxed lower bound, not an optimised full model residual. Its contrast Jacobian is identically zero for all shared parameter directions. Existing AE4 activation gain and timing cannot change it.

This statement assumes the same source faithful optical normalisation and recovery extraction rule. Exact per experiment windows and calibrations were not published. Arbitrary cohort optical gains or selectively chosen windows would weaken inference, but have no evidential basis and are not added here. No calibrated global p value is claimed.

## Genotype resting states

| Genotype | Rest Cl, mM | Rest pH | Local stability |
| --- | ---: | ---: | --- |
| WT | 57.85495031 | 6.88312454 | stable |
| AE4_KO | 56.95413832 | 7.22459780 | stable |
| AE2_KO | 57.89058908 | 6.88133291 | stable |

Each root uses its actual deletion, exact cell and lumen charge constraints and one shared inherited vector. The saved WT state was a numerical seed only. Local stability and residual checks do not establish global uniqueness or chronic biological adaptation.

## Supported whole gland diagnostics

Under the source B+ projection, AE4 KO has a cumulative deficit of -0.50838720% (a small increase in secretion), compared with the measured 35 ± 4.7%. AE2 KO/control water ratio is 1.00002287. These are unfitted predictions from the genotype resting states, not Task 47 acute deletion trajectories.

| Window, s | AE4 KO/WT water ratio |
| --- | ---: |
| 0 to 120 | 1.01060164 |
| 120 to 180 | 1.00571043 |
| 180 to 600 | 1.00328709 |

## Full data accounting and limits

The residual table contains 48 source and derived rows, with 12 quantitatively mapped diagnostic rows. Unavailable predictions and residuals are null in JSON and blank in CSV. Raw gland bicarbonate amounts are retained but only their ratios are scored. No target is silently dropped or assigned a fabricated fluorescence flux conversion.

The independent group Gaussian interpretation uses SE as uncertainty of the reported mean. Correlations and individual animal clustering are unreported. If arbitrary covariance between the two contradictory means is allowed, their contrast SE is at most 0.19 and the discrepancy is still at least 7.3684 reported standard errors. This is not a finite sample p value.

A full numerical observation Jacobian and joint optimum are not defined because essential protocols and observation maps are missing. The exact zero sensitivity of the contradictory contrast is an early rejection gate. The 33 provenance eligible parameters were not treated as an identifiable fit subspace. No optimisation or parameter correction was performed.

NHE and EIPA pH diagnostics are retained, with their unverified timing and optical qualification. Isolated NKCC, isolated exchanger conditioning, sodium absent assays and exact SPQ values remain unavailable. The model therefore cannot claim to satisfy the isolated NKCC constraint or quantify compensation under it.

All numerical trajectories use the qualified positive B+ reservoir projection declared at 48A: exact transported ion concentrations, fixed bath pH, 17 mM spectator osmoles and explicit divalent/carbonate charge bookkeeping. The omitted HEPES/titrant chemistry and duct/acinar assay topology remain limitations. No inherited intracellular equation was changed.

Reproduce only in a separate review workspace using `run_baseline.py` then `summarise_baseline.py`. Existing per genotype and protocol records are reused after matching their input identities. No Tasks 46 or 47 calculation is rerun.
