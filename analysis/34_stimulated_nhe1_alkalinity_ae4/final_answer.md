# Task 34 final answer

`NHE_STIM_TARGET_NOT_QUANTIFIABLE`

1. A quantitative WT initial alkalinisation target suitable for the model was not recoverable. Figure 2D supports a control slope of about `0.59 +/- 0.05 x 10^-3 s^-1`, with a conservative graphical interval of `[0.53, 0.65] x 10^-3 s^-1`, but this is the slope of normalised BCECF `F/F0`, not pH units per second.
2. No `G_NHE_stim` was identified. There were zero gain evaluations.
3. Yes. The accepted Task 31 R09 resting state and all model source files remain unchanged because the stop occurred before implementation or numerical work.
4. Calcium activated sustained secretion was not tested.
5. The sustained NKCC1, AE4, and AE2 positive chloride loading shares were not calculated.
6. AE4 at 5 percent and 0 percent was not simulated, so no secretion ratios exist.
7. Strengthening after 2 to 3 minutes was not assessed.
8. This task cannot decide whether stimulated NHE1 activation restores AE4 secretory leverage.
9. In addition to the mechanism question remaining untested, nothing was established about constitutive genotype specific resting equilibria because no expression trajectories were run.

The source methods specify normalised BCECF fluorescence for these stimulated traces and do not provide the optical ratio to pH conversion needed to compare Figure 2D with an instantaneous model pH slope. Substituting an assumed conversion would make the scalar gain appear identified by information absent from the permitted source.
