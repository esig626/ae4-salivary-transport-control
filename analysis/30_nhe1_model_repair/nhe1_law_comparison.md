# NHE1 law comparison

The repaired model uses Vera-Sigüenza et al. (2018), Eq. 26 exactly, with the squared proton saturation factors printed in that paper.

| Property | Historical comparator | Task 30 evaluator |
|---|---|---|
| Identifier | `legacy_tanh` | `vera_siguenza_2018_eq26` |
| Rate form | Capacity times `tanh` of an ideal activity affinity | Difference of forward and reverse sodium/proton saturation terms |
| Proton factor | Appears through a log activity ratio | Squared saturation on each proton term |
| Activity | Inherited background-specific placeholder | `0.0305 fmol/s` published default; only the permitted WT scalar calibration may change it |
| Constants | Shared log width | `K_H = 4.5e-4 mM`, `K_Na = 15 mM` |
| Status | Selectable only for historical reproduction | Explicitly selected in every Task 30 model |

Fixed state comparison at each inherited saved WT state, using the nominal published `G_NHE1 = 0.0305 fmol/s`:

| Background | Historical flux (fmol/s) | Eq. 26 flux (fmol/s) | Eq. 26 / historical |
|---|---:|---:|---:|
| R09 | 0.007084048618 | 0.00129316249306 | 0.182546 |
| R10 | 0.00700038325673 | 0.000447973396673 | 0.0639927 |

The selected R09 WT state used the permitted WT only scalar calibration, `G_NHE1 = 0.145419726007 fmol/s`. The nominal published value itself was tested first. R10 calibration was not reached before the numerical budget stop.

Positive flux is one sodium ion entering the cell and one proton leaving it. It therefore adds the same amount to intracellular sodium and total alkalinity, adds no TIC directly, and transports no net charge.
