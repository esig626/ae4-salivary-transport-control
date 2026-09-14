# Task 36 manual verification before Codex

No full ODE trajectory has been run yet. This note records the algebraic and local current-closure checks performed before handing anything to an external coding agent.

## Corrected transporter

The active surrogate is

`Na_o + 2 HCO3_o <-> Na_i + 2 HCO3_i`.

For positive inward cycle flux `J_B`, the cell conserved-coordinate source is

`(+J_B, 0, 0, +2 J_B, +2 J_B)`

for `(Na, K, Cl, TIC, TA)`.

Therefore the bulk charge source is exactly `-J_B` equivalents/s and the matching conventional basolateral current is `+F J_B` in the cell-to-bath direction.

## WT stimulated algebra

Using

- Task 31 R09 NKCC1 chloride loading `C = 0.2562444047732706 fmol/s`;
- WT stimulated NKCC1 share `0.70`;
- Task 31 NHE1 `0.007890580489555009 fmol/s`;
- source-fixed NHE1 stimulated factor `2.3`;

gives

`J_AE4 = (3/7) C = 0.10981903061711597 fmol/s`,

`J_NHE1,stim = 0.01814833512597652 fmol/s`,

and from

`J_NHE1 + 2 J_NBC - 2 J_AE4 = 0`

with negligible positive AE2,

`J_NBC = 0.1007448630541277 fmol/s`.

The non-compensated chloride-loading ratios then follow directly:

- AE4 100%: `1.000`;
- AE4 5%: `0.715`;
- AE4 0%: `0.700`.

These were not fitted to the knockout secretion phenotype.

## Independent local current-closure check

At the accepted Task 31 R09 chemical state, with `Ca = 0.25 uM`, the derived fixed NBC capacity

`G_NBC = 0.11570913197464398 fmol/s`

was inserted into the full two-membrane electrical equations with the NBC current included self-consistently.

The scalar basolateral voltage root gave

- `J_NBC = 0.10074486305412181 fmol/s`;
- NBC affinity `A_NBC = 2.6717127432789862`;
- `V_a = -27.482762203344512 mV`;
- `V_b = -80.34467497682778 mV`;
- `I_NBC = 9.720401571160862 pA`.

Current-closure residuals were

- apical: `-1.4217166778255163e-25 A`;
- basolateral: `3.877409121342317e-26 A`.

Combined source/current charge residuals were approximately

- cell: `9.41e-17 fmol/s`;
- lumen: `-5.13e-17 fmol/s`.

These are numerical round-off scale.

## REST nesting

The implementation delegates literally to the pre-existing membrane closure when the secretory activation coordinate is zero. At `Ca = 0.058 uM`, NBC flux and NBC current are exactly zero and the NHE1 multiplier is exactly one. Thus Task 31 REST is nested by construction rather than re-solved or refitted.

## What is not yet claimed

The focused pytest file has been rewritten for the corrected 1:2 electrogenic transporter, but the repository test suite has not been executed in this manual step. No WT, AE4 5%, or AE4 0% dynamic integration has been run yet.
