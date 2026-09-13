# Task 30 handoff for Task 31

Published Task 30 head before Task 31 preparation: `4c74f74ff25103d11335cb38a26a9a3557aca7c6`.

Task 30 implemented Vera-Siguenza et al. (2018) Eq. 26 as an explicit selectable NHE1 evaluator. The old generic tanh law remained selectable but was not used for the repaired REST calculation.

Key completed findings:

- At the inherited saved WT states, the published Eq. 26 law with its nominal `G_NHE1 = 0.0305 fmol/s` produced substantially smaller NHE1 flux than the old generic law.
- R09 WT with the nominal published `G_NHE1` converged but at pH 6.766175, outside the declared WT pH interval.
- A permitted WT-only scalar activity calibration gave `G_NHE1 = 0.145419726007 fmol/s`, with R09 WT pH 6.872665, Na_i 11.565526 mM and volume 1.447604 pL.
- The first R09 exact-AE4-null solve was incomplete when the global 10,000 stationary residual-evaluation budget was exhausted. No accepted null state was obtained and no dynamic trajectory ran.
- R10 was not reached.

Task 31 must not treat the Task 30 scalar calibration as an established salivary NHE1 density. It is a handoff result showing that the 2018 equation shape alone plus a free overall activity did not complete the exact-null test within the previous budget.

Task 31 is instructed to move forward to a mechanistic NHE model with explicit ion-exchange and intracellular proton-modifier regulation, using Cha et al. (2009) as the preferred primary source.
