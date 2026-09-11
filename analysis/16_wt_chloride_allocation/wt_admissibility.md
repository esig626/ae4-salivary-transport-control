# WT-only admission checkpoint

Exactly ten baseline roots and twenty modified root/allocation paths were considered. All ten baseline roots and their thirty calcium trajectories are retained. No modified endpoint equilibrium was obtained inside the inherited continuation domains, so none is admitted to WT dynamics or genotype evaluation. This decision uses only WT evidence.

| Allocation condition | Roots attempted | Valid endpoint roots | WT-admissible root/calcium cases | New WT trajectories |
|---|---:|---:|---:|---:|
| Inherited baseline | 10 reused | 10 reused | 30 | 0 |
| Reference share 0.10 | 10 | 0 | 0 | 0 |
| Reference share 0.30 | 10 | 0 | 0 | 0 |

The baseline AE4 positive-loading share is 1.3553–1.8682%. At share 0.10, the required complete AE4 capacity multiplier is 5.3526–7.3785 and the NKCC1 multiplier is 0.912365–0.917134. At share 0.30 the multipliers are 16.0578–22.1356 and 0.709617–0.713327. No hidden cap is imposed on these scales.

Every modified path has converged, full-rank internal roots before failure. At its last connected internal state, Na is 59.8623–59.9989 mM and K is 61.8959–67.8579 mM, already far outside the WT gates of Na <=35 and K >=100 mM. The final failed local solves touch the inherited Na search ceiling of 60 mM. Their nonzero residuals are retained rather than treated as endpoint solutions. Last connected numerical residuals are 1.71e-12–1.65e-11 in the inherited scaled norm, and their conservation ratios are at most 0.000641 of the permitted limits. Thus the recorded preceding roots are numerically valid but physiologically inadmissible.

The accepted baseline resting AE4 cancellation fraction is 96.8019–99.0077%, with gross opposing-branch traffic 0.13713–0.38449 fmol/s. At the last saved internal numerical states, cancellation is 98.6096–99.5327% and gross traffic is 0.26661–0.68412 fmol/s. These are continuation diagnostics, not results at scientific allocation endpoints. The complete saved roots allow paired comparison without assigning those internal states new target shares. No realized equilibrium share, secretion total or genotype effect is reported for the unresolved 0.10 or 0.30 endpoints.

At an unchanged reference state, multiplying carrier amount increases both opposing AE4 branches and net chloride loading by the same factor, so the cancellation fraction cannot change. Once the state moves during continuation, the observed branch traffic and cancellation become worse in every path; the capacity intervention has not repaired the mixed-cation geometry. This is mechanistic diagnostic evidence, not proof that cancellation causes a secretion phenotype.

The result is consistent across all ten original roots, both routing families and both declared modified endpoints. Because the failures occur at rest, all three calcium cases are explicitly marked unavailable for each failed modified candidate. No calcium or root is selected using a genotype outcome.

The endpoint failures limit the conclusion. A continuation failure inside a frozen search domain is not a proof that no equilibrium exists anywhere, or that a connected branch could never return to a physiological region outside the permitted procedure. No domain expansion, disconnected-root hunt or share-axis refinement is performed. The primary Task 16 classification is therefore `TASK 16 NUMERICALLY INCONCLUSIVE`, rather than an unqualified impossibility theorem. Operationally, this declared study produces zero WT-admissible modified candidates, and licenses zero modified genotype simulations.
