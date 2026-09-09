# Vbeta/VRAC pre-reveal diagnostic

## Evidence and scope

Catalan et al. (2015), *A fluid secretion pathway unmasked by acinar-specific
Tmem16A gene ablation in the adult mouse salivary gland*, PNAS,
doi `10.1073/pnas.1415739112`, supports an isoproterenol-associated apical
anion-exit pathway distinct from the muscarinic Ca-activated TMEM16A current.
The study reported slow IPR-associated acinar swelling and pharmacological
sensitivity consistent with a volume-sensitive anion pathway. These are
primary experimental constraints on topology and protocol context.

The source does not provide a certified whole-cell current-to-conductance map
for this model. Therefore an absolute Vbeta capacity is not licensed and is
not fitted. In particular, no AE4-null secretion observation enters this
diagnostic.

## Nested topology

- `V0_EXACT_CONTROL`: extra apical chloride conductance and current are
  identically zero.
- `V1_RECTIFIED_SWELLING_GATE`:
  `gV * beta_input * max(Vi / Vrest - 1, 0)` is the extra apical chloride
  conductance.

V1 adds no dynamic state. Its chloride current uses the same Nernst and
current-to-molar-flux conventions as the two-membrane model. Intracellular and
luminal chloride sources are equal and opposite, and their charge conversion
matches the conventional current exactly. The module is independent of the
`N_abs` and `AN_abs` source-class labels, so the same diagnostic composes with
both.

The rectified gate is a `NEW_MODELING_DECISION`, not a measured VRAC activation
law. Every positive `gV` is diagnostic-only because its absolute capacity is
unlicensed. V0 is an exact nesting control but omits the source-supported
IPR-associated pathway; V1 has a source-supported topology but is not
production eligible.

## Exact-rest anti-deadlock rule

V1 is evaluated on an IPR-only trajectory generated without V1. Starting from
the exact resting volume, some pre-V1 process must first make the cell swell.
V1 cannot be credited with producing the swelling required to activate
itself. An initially swollen state is not accepted as an escape from this
rule.

The earlier superseded reference root shrank from `1.3` to
`1.2737304841 pL` over the 600-s IPR-only protocol with resting Ca
`0.058 uM`, R1 beta/PKA/AE4 dynamics, and `tau=30 s`. Its maximum positive
swelling and V1 gate were exactly zero, so V1 deadlocked for every value of
`gV`.

The native strict/conditional-pump scan currently has 35 numerical roots and
zero WT-resting-gate passes. Radau and BDF both reproduced the IPR-only audit
for all 35 roots. Thirty roots had no positive gate. Five conditional
`pump=0.5`, `AE4NA05` roots showed only `0.143-0.183%` maximum transient
swelling, but all five fail the WT resting gates and cannot support a
production rescue. Solver gate classifications agreed for every root; the
maximum Radau/BDF difference in maximum relative swelling was
`1.94e-12`. The audit will be rerun after the `N_abs` and `AN_abs` root panels
are appended.

Machine-readable outputs are
`results/13B_modern_full_model/vbeta_native_exact_rest_diagnostic.csv` and
`results/13B_modern_full_model/vbeta_native_exact_rest_summary.json`.

## Decision

V1 is not a valid production repair in the present hierarchy. A WT-valid root
would still need to pass the anti-deadlock test, and a source-supported
absolute current map would then be required before V1 could affect the model.
The Catalan topology remains biologically important; rejecting this
unlicensed implementation is not evidence against that pathway.
