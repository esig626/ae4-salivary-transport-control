# Phenotype target and minimal-signalling convention

This file is binding for all future scientific reconstruction in this repository.

## Primary modelling objective

The primary AE4-loss target is the **magnitude and direction of the secretion phenotype**, not reproduction of the experimental clock.

The model should produce a substantial reduction in AE4-null secretion relative to WT under the stimulated condition, consistent with the reported approximately 35% ten-minute deficit and its experimental uncertainty/context.

The reduction should be a robust consequence of the model rather than an artefact of choosing one particular observation time. A candidate that only crosses the desired deficit at a finely selected instant is not acceptable. Conversely, a candidate must **not** be rejected merely because its onset, peak, or minute-by-minute time course does not reproduce the experimental trace exactly.

Experimental time is biological protocol information, not an exact mathematical trajectory constraint. The reported early-versus-late pattern may be used as qualitative context, but it is not a hard fitting, acceptance, or mechanism-selection criterion unless a future task explicitly supplies data precise enough to justify that use.

For future model comparison, prefer quantities such as:

- sustained or stationary WT versus AE4-null secretion ratio where the model admits such a regime;
- cumulative or time-averaged secretion ratio over a broad stimulation window;
- robustness of the secretion deficit across reasonable observation times;
- physiological state and conservation checks.

Do not optimise an onset time, half-time, minute landmark, or exact experimental trajectory shape.

## Minimal treatment of signalling

Do not build signalling networks merely because the experimental agonists are CCh, isoproterenol, beta-adrenergic, cAMP, PKA, or calcium linked.

If signalling is needed only to represent the **effective stimulated activity** of an existing transporter/channel, reduce it to the smallest justified parameter or algebraic multiplier. Prefer a static/effective protocol parameter over a signalling ODE, intermediate molecular states, delays, phosphorylation cycles, or cascades.

A dynamic signalling state is permitted only if the scientific result cannot be represented without memory/history dependence and there is independent evidence that identifies the added dynamics. It must not be introduced to reproduce an experimental time course.

In particular:

- existing source-supported AE4 PKA dependence may be represented by an effective stimulated AE4 activity factor;
- any beta/cAMP effect on another pathway should likewise be represented by the smallest source-supported effective multiplier if it is actually needed;
- do not add beta/cAMP/NKCC/VRAC signalling complexity merely to create the AE4-null secretion phenotype;
- no signalling parameter may be chosen from the approximately 35% AE4-null secretion target unless the task explicitly declares that target as calibration rather than validation, which is not the default project policy.

## What still matters experimentally

Experimental constraints on resting chloride, pH, transporter direction, genotype effects, inhibitor responses, and secretion magnitude remain scientifically important. They constrain whether a model is physiologically coherent.

The distinction is:

- **state/flux phenotype:** constrain the model;
- **exact experimental timing:** normally does not;
- **signalling detail:** collapse to an effective parameter unless independently required.

This convention supersedes any earlier task wording that treated the first 2-3 minutes, delayed divergence, or an exact signalling time constant as a hard success/failure gate for the final AE4-loss reconstruction.
