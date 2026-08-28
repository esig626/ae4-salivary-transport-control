# Task 13B whole-cell model architecture

## Status and scope

This document freezes the transparent whole-cell architecture implemented in
`src/modern_full_model/`. It is a Task 13B generation scaffold, not an accepted
WT calibration. No AE4-null secretion magnitude, ratio, or time course was
used in its equations, starting values, parameter values, topology, root
selection, or tests.

The architecture replaces the inherited seven-state comparison chassis with
conserved intracellular and luminal amount coordinates, finite volumes,
finite-buffer carbonate chemistry, two-membrane current closure, explicit
apical/basolateral cation topology, paracellular ions, osmotic water flux, and
a finite-lumen outflow. The dedicated AE4 and beta/cAMP/PKA modules connect
through narrow interfaces rather than being folded into a residual patch.

This work does **not** address reduction, GSPT, bifurcation analysis,
identifiability, or manuscript drafting.

## Claim and provenance policy

The implementation uses all five Task 13B provenance categories.

| Category | Use here |
|---|---|
| `PRIMARY_MEASUREMENT` | Qualitative transporter/topology requirements delegated to the evidence, AE4, and regulatory ledgers |
| `PUBLISHED_MODEL` | Historical-lineage pathway inventory and selected starting forms, never treated as physiological ground truth |
| `HISTORICAL_IMPLEMENTATION` | Comparison/nesting target only |
| `DERIVED_CONSTRAINT` | Stoichiometry, charge, current-to-flux conversion, Nernst reversals, carbon/alkalinity identities, conservation laws |
| `NEW_MODELING_DECISION` | Coordinate choice, effective reversible rate laws, finite-lumen closure, unmeasured capacities/fractions, executable starting values |

The exact parameter registry is produced by
`modern_full_model.parameters.parameter_records`. Default numerical values are
positive integration-test starting values. They are not a frozen calibration.
Parameter construction rejects nonpositive physical constants and half-scales,
negative capacities, conductances, permeabilities, hydraulic coefficients,
outflow rates, buffer/impermeant amounts, or invalid initial concentrations.
Zero remains legal for an optional pathway or nested-limit capacity.

## Full state vector

Amounts are in fmol and volumes in pL. Since `1 mM * 1 pL = 1 fmol`, each
concentration is exactly amount divided by volume in mM.

| Index | State | Meaning |
|---:|---|---|
| 0 | `na_i_fmol` | intracellular Na amount |
| 1 | `k_i_fmol` | intracellular K amount |
| 2 | `cl_i_fmol` | intracellular Cl amount |
| 3 | `tic_i_fmol` | intracellular total inorganic carbon, CO2 + HCO3 + CO3 |
| 4 | `alk_i_fmol` | intracellular total alkalinity equivalents |
| 5 | `volume_i_pL` | cell volume |
| 6 | `na_l_fmol` | luminal Na amount |
| 7 | `k_l_fmol` | luminal K amount |
| 8 | `cl_l_fmol` | luminal Cl amount |
| 9 | `tic_l_fmol` | luminal total inorganic carbon |
| 10 | `alk_l_fmol` | luminal total alkalinity equivalents |
| 11 | `volume_l_pL` | finite local lumen volume |
| 12 onward | regulatory-model-specific | cAMP/PKA/AE4 regulatory states from `CampPKAModel.state_names` |

The model therefore contains intracellular and luminal Na, K, Cl, carbon,
acid-base, and volume degrees of freedom. H, HCO3, CO2, CO3, buffer occupancy,
and pH are explicitly recoverable algebraic observables. Apical and
basolateral potentials are algebraic current-closure variables. This is a
deliberate fast-electrical/fast-acid-equilibrium closure, not an omission of
charge or carbon bookkeeping.

If independent WT dynamics show that carbonic-anhydrase relaxation or membrane
capacitance is not fast, one independent carbon species state or capacitive
voltage state must be added without changing the conserved amount rows.

## Compartment and voltage conventions

The bath is an open reservoir. The cell and lumen are finite compartments.

\[
V_a=\phi_i-\phi_l,\qquad V_b=\phi_i-\phi_e,\qquad
V_t=\phi_l-\phi_e=V_b-V_a.
\]

Positive conventional current points from the first named compartment to the
second. For valence \(z\),

\[
J_{1\to 2}=\frac{I_{1\to 2}}{zF}.
\]

Consequently a positive conventional chloride current corresponds to chloride
motion in the opposite physical direction. This conversion is centralized in
`membranes.current_to_fmol_s`.

With negligible capacitive charge on the modeled time scale, the two closure
equations are

\[
I_a-I_t=0,\qquad I_b+I_t=0,
\]

where \(I_t\) is positive lumen-to-bath paracellular current. The implementation
solves the resulting two-by-two linear system for \(V_a,V_b\) at every RHS
evaluation and reports both residuals.

## Whole-cell source assembly

Positive non-electrical transporter directions are:

| Flux | Intracellular amount signature `(Na,K,Cl,TIC,TA)` | Interpretation |
|---|---|---|
| `J_NKCC1` | `(1,1,2,0,0)` | bath-to-cell Na/K/2Cl cotransport |
| `J_NHE1` | `(1,0,0,0,1)` | Na influx and H extrusion; H removal raises TA |
| `J_AE2` | `(0,0,1,-1,-1)` | Cl influx and HCO3 extrusion |
| `J_AE4,Na` | `(-1,0,1,-2,-2)` | source-supported 1:1:2 Na branch used by retained SR2 family |
| `J_AE4,K` | `(0,-1,1,-2,-2)` | source-supported 1:1:2 K branch used by retained SR2 family |
| apical pump cycle | cell `(-3,+2,0,0,0)`, lumen `(+3,-2,0,0,0)` | 3 Na out/2 K in at the apical pole |
| basolateral pump cycle | cell `(-3,+2,0,0,0)` | 3 Na out/2 K in from open bath |
| apical K movement | cell K loss, matching lumen K gain | current converted with `I/F` |
| basolateral K movement | cell K loss to bath | current converted with `I/F` |
| apical Cl movement | matching cell/lumen Cl exchange | current converted with `I/(-F)` |
| paracellular ion movement | lumen loss to bath for positive `J_lumen->bath` | Na, K, Cl, and HCO3 paths |
| basolateral CO2 | cell TIC source only | neutral open-bath exchange |
| apical CO2 | equal/opposite cell and lumen TIC sources | neutral internal exchange |

Each NKCC1, NHE1, and AE2 law is reversible and vanishes at its ideal-solution
thermodynamic reversal. Their current constitutive `tanh(log activity ratio)`
forms are minimal bounded modeling decisions. They should be replaced or
parameterized only from independent transporter/WT evidence. Their source
signatures are exact and remain valid if the constitutive law changes.

The AE4 module is the separate QSS shared-carrier SR2 implementation in
`transporters.py`. It returns sources directly in the same conserved
coordinates and reports local-detailed-balance, occupancy, affinity, entropy,
and charge diagnostics. The whole-cell model passes the dynamic regulatory
capacity multiplier into AE4; it does not hide regulation inside a fitted
whole-cell coefficient.

## Charge identity

For the total-alkalinity definition used here, modeled bulk charge is exactly

\[
Q_i/F=n_{Na,i}+n_{K,i}-n_{Cl,i}-n_{TA,i}-n_{X,i},
\]

\[
Q_l/F=n_{Na,l}+n_{K,l}-n_{Cl,l}-n_{TA,l},
\]

where \(n_X\) is the fixed intracellular anion-equivalent amount. NKCC1,
NHE1, AE2, and both AE4 branches have zero source charge under this identity.
Pump/channel charge is reconciled to current through Faraday's constant.
Initial test coordinates are exactly electroneutral to roundoff; accepted WT
coordinates must independently satisfy the same gate.

Bath chloride, fixed-cell anion equivalents, and reference luminal chloride
are stored derived values. Model construction recomputes charge from the active
bath/reference tuple and rejects stale combinations. Thus changing bath
cations/TIC/pH or reference ions/TIC/pH/buffer/volume requires updating the
corresponding derived charge value in the same frozen parameter tuple.

## Water, dilution, and outflow

Ideal osmolarity includes Na, K, Cl, TIC, the finite buffer pool as one osmotic
particle per buffer molecule, and other impermeant osmoles. The parameter
`cell_impermeant_osmoles_fmol` means *other* impermeant particles and explicitly
excludes `cell_buffer_total_fmol`; the two are added exactly once. Protonation
changes buffer charge but not its particle count.
Water signs are

\[
q_b:e\to i,\quad q_a:i\to l,\quad q_t:e\to l,\quad q_{out}:l\to duct.
\]

Thus

\[
\dot V_i=q_b-q_a,\qquad
\dot V_l=q_a+q_t-q_{out}.
\]

All mobile luminal quantities leave as \(q_{out}c_l\). Because the states are
amounts, changing volume automatically changes concentration; no independently
coded dilution terms can double-count that effect. The current local outflow
law is a declared compliance closure and is not yet a secretion calibration.

## Dynamic stimulus and regulation hooks

A stimulus protocol returns two distinct physical inputs:

* intracellular calcium in uM, which gates Ca-activated K and Cl conductance;
* a dimensionless beta-adrenergic input, which drives the dedicated cAMP/PKA
  subsystem.

The regulatory subsystem contributes its own named state suffix, RHS, and AE4
capacity multiplier. Static `R0` remains only an unconfigured/control path.
The integrated model neither turns beta input into calcium automatically nor
uses an AE4-null timing feature to choose a regulatory time constant.

## Genotype operations

`AE4_NULL` sets the AE4 expression multiplier exactly to zero before AE4 is
evaluated. It does not alter NKCC1, NHE1, AE2, pump, conductance, water,
stimulus, regulatory kinetics, initial conditions, or solver settings.
`AE2_NULL` analogously zeros only AE2. Any genotype-specific adaptation would
require separate primary evidence and a new declared model generation.

## Stable executable API

The integration-facing interface is:

* `ModernFullModel.state_names`: complete named state tuple;
* `ModernFullModel.initial_state()`: positive reference state;
* `ModernFullModel.evaluate(t,y,genotype=...)`: RHS and all diagnostics;
* `ModernFullModel.rhs(t,y,genotype=...)`: solver RHS;
* `ModernFullModel.solve_resting_root(...)`: diagnostic bounded residual
  attempt only; it always reports `production_eligible=False` because it does
  not certify row rank, parameterize the charge manifold, track branches, or
  use a production multistart design;
* `ModernFullModel.solve_dynamics(...,method='Radau'|'BDF')`: physical-seconds
  stiff integration.

Calibration, validation, and run-analysis scripts remain separate by design.

## Exact nesting limits

1. `apical_pump_fraction = 0` removes apical pump capacity exactly.
2. `apical_k_fraction = 0` with zero apical K background removes the apical K
   path exactly.
3. AE4 expression zero removes every AE4 ionic and carbon source exactly.
4. AE2 expression zero removes only the AE2 source column.
5. Beta input zero is handled by the selected regulatory model's exact
   zero-input equilibrium; it does not delete basal AE4 capacity implicitly.
6. A stateless regulator gives the `R0` control, but cannot satisfy Task 13B's
   mandatory dynamic-regulation requirement for an accepted generation.

The amount-state architecture does not nest to the old seven-state system by
silently deleting conserved rows. A historical-like reduction would require
explicit limits: fixed lumen volume, algebraic lumen chloride, infinite/fast
acid-base constraints, and basolateral-only cation fractions. Such a limit is
for later verification, not a reason to retain the old chassis now.

## Implemented invariant tests

`tests/test_modern_full_model_whole_cell.py` currently certifies:

* amount/concentration units and full state naming;
* pH/species recovery from TIC/TA with a finite buffer;
* closed-reaction carbon, buffer-site, and charge invariants;
* NKCC1/NHE1/AE2 reversal and electroneutrality;
* rank and nesting of apical pump/K source directions;
* coupled apical/basolateral current closure and Faraday conversion;
* whole-model carbon, charge, buffer, and volume accounting;
* exact AE4 source deletion;
* positive short integrations and Radau/BDF agreement;
* rejection of out-of-domain physical parameters;
* explicit one-time osmotic counting of the finite buffer pool;
* rejection of stale derived electroneutrality values and acceptance of a
  coordinated nondefault charge tuple;
* absence of held-out secretion targets from the implementation.

These tests certify equations and interfaces, not WT physiology.

## Open integration and acceptance issues

1. Replace every executable starting value by a source-ledger value, a bounded
   WT calibration parameter, or an explicit removed parameter.
2. Establish cell area/volume and whole-gland scaling before reporting flow.
3. Calibrate only the regulatory subsystem from 2021 WT/regulatory data, then
   integrate its retained family.
4. Establish whether AE4 carrier QSS is independently fast relative to the WT
   regulatory and whole-cell time scales.
5. Calibrate NKCC1, NHE1, AE2, pump, conductance, buffer capacity, CO2 exchange,
   and water/outflow only against allowed WT/transporter evidence.
6. Demonstrate a bounded, non-boundary WT resting root; the reference state is
   not itself a root claim and the bundled diagnostic helper cannot satisfy
   this gate. The production root-geometry workflow must independently certify
   row rank, enforce the charge manifold, use non-collinear multistarts, and
   test branch uniqueness.
7. Demonstrate a physical-time WT stimulated trajectory before any normalized
   genotype comparison.
8. Freeze topology, parameters, inputs, time map, root rule, and solver
   tolerances before held-out AE4-null secretion is opened.

Until those gates pass, this architecture cannot justify later reduction or
any substantive Task 13B success classification.
