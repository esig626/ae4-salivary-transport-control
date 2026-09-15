# Joint experimental reconstruction of AE4 dependent secretion

The inherited architecture does not provide a defensible joint reconstruction of the observed AE4 phenotype. A central failure is structural: in cells lacking AE4, the model cannot distinguish CCh stimulation from CCh plus IPR, whereas the measured chloride uptake differs substantially. The model already contains explicit AE4 activation by beta adrenergic input. Adding that activation again, or changing its strength, cannot repair a discrepancy that persists when AE4 is absent.

This is a negative result within the adopted equations and a common source faithful observation rule. No corrected architecture or completed joint numerical optimum is claimed. The evidence does not uniquely specify the one further equation repair permitted by Task 48, so no additional mechanism was introduced.

## Experimental reconstruction

The primary dataset is [Peña Münzenmayer et al. (2015)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4409235/). All 28 Table 1 group means, standard errors and sample sizes were checked against the full paper. The remaining numerical seed observations were also verified. The complete ledger retains the resting states, stimulus specific chloride measurements, isolated NKCC and exchanger assays, NHE response, AE2 controls, bicarbonate output and secretion timing.

The source reports means with SE. Its sample counts refer to experiments from preparations involving at least three animals per condition, not necessarily independent animal counts. AE4 and AE2 controls have different genetic backgrounds. Comparisons retain their own controls; the shared model does not introduce strain dependent fitted multipliers. In particular, it does not reproduce the different initial control secretion curves as separate strain predictions.

Two source details affect interpretation. SPQ observations use inverse fluorescence, F0/F, and uptake is a regression slope after initial chloride loss. Low chloride depletion and readdition belong to the isolated transporter protocols, not the Table 1 agonist experiment. The resting chloride values are already calibrated, but the paper does not supply the individual optical calibrations or all source regression windows. The reported bicarbonate amounts are retained despite a dimensionally inconsistent sentence in the Methods concerning their calculation.

The correct DOI for the 2016 biophysics paper is [10.1085/jgp.201611571](https://pmc.ncbi.nlm.nih.gov/articles/PMC4845690/), replacing the incorrect seed DOI. Its external cation curves do not measure the fixed physiological Na/K routing fraction used by the model. The [2021 PKA study](https://pubmed.ncbi.nlm.nih.gov/34585968/) supports the activation pathway and S173 dependence; it does not identify a native activation time constant. Its figure errors are SD, distinct from the primary 2015 SEs.

## Protocols, resting states and observations

All calculations use one shared inherited parameter vector. The genotype changes only the deleted transporter. CCh and IPR remain separate inputs. Acute bumetanide, EIPA and T16Ainh act on NKCC, NHE and the chloride channel respectively; they are not substituted for unrelated capacity changes. Their effects are applied through every public evaluation and integration entry.

The supported positive bath uses source Na 145, K 4.3, Cl 128.3 and bicarbonate 25 mM at pH 7.4 and 37 °C. Bicarbonate is converted to total inorganic carbon using the inherited carbonate equilibrium. Calcium, magnesium, glucose and HEPES contribute 17 mM external spectator particles. Divalent charge and the small carbonate/titrant accounting term are explicit. This is a qualified reservoir projection: unreported HEPES protonation and titrant quantities are not reconstructed. The original intracellular transport equations remain unchanged, including the limitations of their bath conditioned NKCC law.

WT, AE4 KO and AE2 KO resting states were solved separately on the exact cell and lumen charge constraints. The earlier WT state was only a numerical starting guess. Each solution is positive and locally asymptotically stable. Independent Jacobians in linear amount and volume coordinates agree with the calculation in logarithmic coordinates; the largest relative eigenvalue discrepancy is 1.23 × 10⁻⁷. The largest full resting RHS component is about 2.85 × 10⁻¹² in its native amount or volume units. The AE4 KO has a slow local decay time of about 106 minutes. These checks establish valid local equilibria of the adopted model, not global uniqueness, attraction from every initial state or physiological agreement.

Concentration observations include changing volume. In particular,

\[
\dot C_{\mathrm{Cl}}=(\dot n_{\mathrm{Cl}}-C_{\mathrm{Cl}}\dot V)/V.
\]

An isolated transporter flux is therefore not generally the measured concentration slope. Under the explicitly conditional quenching model F(C)=A/(1+KC), the inverse SPQ signal is

\[
Z(t)=\frac{1+KC(t)}{1+KC(0)}.
\]

The unknown factor does not cancel from genotype slope ratios when their starting chloride differs. Unsupported fluorescence conversions are left unavailable, rather than assigned an arbitrary calibration. The code requires an explicit calibration and recovery window before returning such a slope.

Gland secretion is compared through ratios of integrated outflow. Bicarbonate amount is the integral of outflow times luminal bicarbonate, not total inorganic carbon or alkalinity. This observation map assumes a common effective acinus to gland transfer and does not include ductal bicarbonate processing. Bicarbonate residuals use a Gaussian profile contrast with uncertainty in both reported group means; the raw amounts are not scored again.

## A parameter independent conflict

Write the twelve core variables as x and the AE4 activation state as a. Beta input enters the inherited core only through AE4 transport. With AE4 deleted,

\[
\dot x=F_{\mathrm{KO}}(x,u_{\mathrm{CCh}};\theta),
\]

which is independent of beta input and a. CCh alone and CCh plus IPR have the same CCh input and the same genotype resting state. Uniqueness of the initial value problem within the positive model domain therefore gives identical core trajectories for every admissible shared parameter vector. The two activation trajectories can differ, but neither affects the knockout core.

A common source faithful optical normalisation and recovery extraction rule must consequently predict the same expected uptake in both experiments. Table 1 instead reports:

| AE4 KO stimulus | Uptake, 10⁻³ s⁻¹ | SE | Experiments |
| --- | ---: | ---: | ---: |
| CCh | 2.30 | 0.10 | 6 |
| CCh plus IPR | 0.90 | 0.09 | 9 |

Even relaxing the model to allow any common predicted uptake m gives

\[
\min_m\left\{\frac{(m-2.30)^2}{0.10^2}+\frac{(m-0.90)^2}{0.09^2}\right\}
=\frac{1.40^2}{0.10^2+0.09^2}=108.2873.
\]

The minimising common value, 1.52652, is a relaxed statistical projection, not a fitted model prediction. The contrast differs by 10.406 reported standard errors under independent mean errors. Allowing any covariance consistent with the two reported SEs gives a contrast SE no larger than 0.19, leaving at least 7.368 reported standard errors. These are working comparisons of reported means and estimated SEs, not finite sample calibrated p values.

The contrast Jacobian is identically zero in all 33 provenance eligible shared parameter directions. No optimisation within those directions can remove this discrepancy. A second exact implication is zero IPR only uptake from an AE4 KO resting equilibrium, whereas the reported value is 0.20 ± 0.03. Again, a constant chemical state has zero optical slope for any fixed calibration.

The observation qualification is material. Arbitrary systematic calibration changes between cohorts, or different unreported windows selected from one identical nonlinear curve, could break measurement invariance. The paper does not establish such an explanation, and the task does not authorise fitted protocol specific optical factors. The conclusion is therefore conditional on a common source faithful observation rule and the adopted agonist input separation, not a theorem about every possible measurement process.

Runtime checks confirm the algebra: the knockout core RHS is bit identical under the two stimulus arms even when their activation states are set to zero and one. Separately integrated core trajectories agree within 4.16 × 10⁻⁹ in native state coordinates. No transporter capacity sweep was needed.

## Numerical joint diagnostics at the unchanged vector

The positive bath projection and genotype resting states give the following unfitted predictions. Residuals have sign model minus experiment.

| Quantity | Experiment, mean ± SE | Model | Standardised residual |
| --- | ---: | ---: | ---: |
| WT resting chloride, mM | 50.10 ± 1.50 | 57.85495 | 5.170 |
| AE4 KO resting chloride, mM | 36.50 ± 1.60 | 56.95414 | 12.784 |
| WT resting pH | 6.91 ± 0.07 | 6.88312 | −0.384 |
| AE4 KO resting pH | 6.89 ± 0.02 | 7.22460 | 16.730 |
| AE4 KO saliva reduction, percent | 35.0 ± 4.7 | −0.50839 | −7.555 |

Resting chloride falls by only 0.90 mM after deletion, accompanied by a pH rise of 0.341. This misses the measured low chloride state with preserved pH. The secretion prediction is a small increase, not the required sustained loss. It is distinct from the earlier 3.8575% result, which used an acute deletion from the saved WT state and a different bath convention. Task 47 was not rerun.

| Secretion summary | AE4 KO/WT |
| --- | ---: |
| Water, 0 to 120 s | 1.010602 |
| Water, 120 to 180 s | 1.005710 |
| Water, 180 to 600 s | 1.003287 |
| Total bicarbonate, 0 to 600 s | 1.014677 |

The source reports comparable early secretion followed by substantial later loss. These diagnostic windows show no later loss. They are summaries chosen to represent the source's qualitative timing, not digitised mean curves or invented equivalence bands.

The AE2 KO/control water ratio is 1.0000229, so the inherited point preserves near neutral AE2 secretion in this conditional comparison. Its bicarbonate ratio is 0.9999585. The observed bicarbonate ratios are imprecise: the corresponding scale profile residuals are −1.056 for AE4 and −1.247 for AE2. These findings do not establish that missing duct processing is negligible.

The model's first ten seconds of pH change are about 0.01826 in WT and 0.001824 in AE4 KO. Simultaneous EIPA still leaves a positive pH change in the tested projection. This is in tension with the qualitative NHE fingerprint, but source inhibitor preincubation, exact windows and BCECF calibration are not identified. The values are labelled conditional pH diagnostics, not measured NHE activity or a completed assay reconstruction.

## Why no further correction was accepted

The inherited AE4 regulatory state obeys da/dt=(beta−a)/30 s and supplies the capacity factor 1+0.25a. Its presence and use were verified in the executed model. The gain is an inherited effective setting, not an identified native gland parameter. However, changing either gain or timing affects no knockout core term. The first permitted correction therefore cannot address the decisive residual.

A further correction would need to explain beta associated behaviour that survives AE4 deletion, subject to the measurement qualification. The joint data do not uniquely identify the relevant physiological equation. They also expose independent assay limitations: the selected NKCC law has no extracellular substrate dependence during depletion and readdition; instantaneous carbonate equilibrium cannot represent carbonic anhydrase inhibition; sodium absent conditions lie outside the inherited constitutive domain; and isolated acinar boundary conditions are insufficiently specified. A collection of speculative fixes would violate the prescribed single repair rule.

No new mechanism, genotype multiplier, transporter cap, secretion penalty or protocol specific fitted gain was introduced. The full numerical likelihood and its Jacobian cannot be assembled from the currently available protocol and observation maps. It would be misleading to call a fit to the remaining subset a joint reconstruction. No such optimum, parameter uncertainty interval or global fit p value is reported.

The complete residual ledger contains 48 source or derived records, with 12 mapped diagnostic rows. Unavailable values are null in JSON and blank in CSV. Qualitative observations remain qualitative; a nonsignificant result is never converted into equality or a percentage band. Missing covariance, optical information and assay dynamics remain visible.

## Article interpretation and limits

The result does not establish that NKCC compensation in an isolated assay remains excessive: that assay cannot be faithfully simulated by the present equations. It establishes that the inherited architecture has the wrong conditional stimulus dependence after AE4 deletion, and that its unchanged parameter vector also fails the measured resting and sustained phenotype. Existing AE4 PKA regulation is insufficient to resolve these problems.

The [NKCC1 knockout study](https://pubmed.ncbi.nlm.nih.gov/10831596/) and [NHE1 knockout study](https://www.jbc.org/article/S0021-9258(20)89821-4/fulltext) were checked as secondary context. No numerical cross validation is claimed because there is no successful primary fit. Their parotid setting and chronic expression compensation cannot be silently equated with the present SMG shared deletion model.

The next scientific decision requires evidence that identifies the missing beta associated response and a reproducible observation map. This report does not start that work. It stops at Task 48 with a recoverable negative result, preserving the distinction between a conditional structural contradiction, an unfitted numerical failure and an unavailable assay.

## Reproducibility and publication

The scientific start was `5c945c4a7857657b5f7c76d5bf9e3d3abc7d5721` on the designated Task 48 branch. Five workers independently reviewed sources, protocols, observations, provenance and feasibility. Only the orchestrator integrated changes, accepted scientific conclusions, ran the three resting solves and twelve trajectories, and published the branch. The Task 43 and 44 outputs were read at their pinned commits; Tasks 46 and 47 were not repeated.

Every milestone contains CURRENT_STATUS.md and a machine readable receipt. Each remote commit was verified before dependent work. 48D records that no correction was accepted; 48E records blocked joint optimisation, not a successful fit. Final residuals are identical to the published baseline residuals. No retuning or new scientific run followed the final comparison.

The principal files are EXPERIMENTAL_CONSTRAINT_LEDGER.md, constraints.json, PROTOCOL_MAP.md, PARAMETER_ADMISSIBILITY.md, OBSERVATION_OPERATORS.md, BASELINE_RESIDUALS.md, CORRECTION_GATE.md, FINAL_INFERENCE_STATUS.md and RECONSTRUCTION_DECISION.md. Machine readable results, parameter provenance, saved states, trajectories, identifiability and verification records are under output/. Published output and source hashes must match before any cached work is reused. No frozen production source, previous numbered analysis, manuscript or main branch was modified.
