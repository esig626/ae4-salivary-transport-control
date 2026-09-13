# Task 34 source target assessment

## Primary source

Peña Münzenmayer et al. (2015), *Journal of Biological Chemistry* 290, 10677 to 10688. DOI `10.1074/jbc.M114.612895`; PMCID `PMC4409235`.

The primary article and its Figure 2 image were inspected directly. The downloaded NCBI Figure 2 image was 700 by 809 pixels and had SHA256 `9b4c95983ff33c806994ec3d6f81a597e44842e2ddf1783dfca5a8073f84c61a`.

Article: <https://pmc.ncbi.nlm.nih.gov/articles/PMC4409235/>

Figure image: <https://cdn.ncbi.nlm.nih.gov/pmc/blobs/1ab3/4409235/060bebbe6b0b/zbc0211514610002.jpg>

## What Figure 2C and 2D support

Figure 2C reports the response to 0.3 micromolar CCh plus 5 micromolar IPR. Figure 2D reports the control bar at approximately

`0.59 x 10^-3 s^-1`, with an error bar of approximately `0.05 x 10^-3 s^-1`.

The error bar endpoints read from the graph are approximately `0.54 x 10^-3 s^-1` and `0.64 x 10^-3 s^-1`. Allowing for pixel resolution and line thickness gives a conservative graphical interval of `[0.53, 0.65] x 10^-3 s^-1`.

These values were obtained by direct pixel digitisation of the control bar and its error bar against the labelled axis. The paper states that results are means plus or minus standard error, with `n = 13` for the pooled control bar.

## Why this is not a model pH slope target

The methods state that stimulated BCECF measurements were normalised fluorescence ratios, `F/F0`, except when intracellular pH itself was determined. Figure 2D therefore has units of inverse seconds because it is the slope of a dimensionless normalised optical signal. It is not an absolute slope in pH units per second.

For a BCECF ratio `R`, the measured quantity is locally

`d(R/R0)/dt = [R'(pH0)/R0] dpH/dt`.

The article does not provide `R'(pH0)/R0`, the stimulated trace calibration curve, or the raw 490/440 ratios needed to recover it. The separate high potassium and nigericin procedure used to determine resting pH does not supply this missing conversion for Figure 2C or 2D.

Treating `F/F0` as `pH/pH0`, or importing a generic BCECF calibration from another experiment, would add an unsupported calibration assumption. It would also violate the instruction to derive the numerical pH slope target from this primary figure alone.

## Decision

The plotted optical slope is quantitatively recoverable, but a defensible WT initial alkalinisation target in model pH units is not.

`NHE_STIM_TARGET_NOT_QUANTIFIABLE`

The required source stop was applied before any gain evaluation, model change, code path test, stationary solve, or dynamic integration.
