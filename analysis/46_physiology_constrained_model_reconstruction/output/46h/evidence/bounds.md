# WT physiology and bounds evidence review

Frozen base: `9553d04442aeef5b0a8e0ba5f9df07f0986d7074`.

The preserved Task 43 provenance audit provides no measured finite salivary capacity bound for the principal transport routes. The old model supplies a reference state and conditional effective parameters. It does not supply an independently measured feasible flux region. Task 43 is pinned to `547f113d9123ab1976774639a69483faf40cef34`; Task 44 is pinned to `e5fa9840bf96be3147c6118daee4942468c78f8b`.

## Numerical evidence and permissible use

| Quantity | Preserved value | Evidence status and CBM implication |
|---|---:|---|
| NKCC scale | 0.017334746894096052 fmol/s | Derived from modelled resting uptake divided by the adopted kinetic shape. Neither measured carrier number nor a sustained capacity. |
| Legacy NKCC capacity | 0.32 fmol/s | Inactive in the selected Palk law. Must not be imported as its physical cap. |
| NHE carrier amount | 0.00002339370005697548 fmol | Calibrated to resting pH 6.91 in one background; not a measured salivary abundance. |
| Stimulated NHE multiplier | 2.3 | Transferred effective setting; no justified salivary uncertainty interval. |
| NBC capacity | 0.11570913197464398 fmol/s | Derived from an assumed 70/30 chloride partition and electrical closure. Not an experimental capacity. |
| AE2 capacity | 0.005 fmol/s | Explicit small placeholder. |
| Na/K pump capacity | 0.08 fmol cycles/s | Effective model setting. |
| Apical pump fraction | 0.075075 | Localisation supported; quantitative fraction not measured. |
| Apical K fraction | 0.3 | Apical K current supported; quantitative fraction not measured. |
| AE4 carrier amount | 0.1501264265070657 fmol | Inherited effective model setting; thermodynamics fixes rate ratios, not abundance. |
| Basolateral/apical CO2 coefficients | 0.02 / 0.01 fmol/s/mM | Effective exchange settings. |
| CaCC conductance | 0.0000000314 S | Effective model setting, not a protocol matched maximal flux. |
| Total K conductance | 0.000000014 S | Effective model setting, not a protocol matched maximal flux. |

Exact entries, source paths and uncertainty classifications are in `worker_bounds.json`. The local Task 43 inventory was checked byte for byte against the pinned Git object.

Task 37 reported an NKCC positive loading context band of 65 to75 percent, but did not make it a hard physiological gate. The later Task 43 audit explicitly states that the 70/30 split used in NBC construction is assumed. Neither is an independently established sustained transporter capacity or sufficient grounds to impose a WT partition in this screen.

The inherited concentration gates are Na below40 mM, K50 to200 mM, Cl30 to80 mM, pH6.6 to7.3, and positive cell volume below3 pL. These are model acceptance criteria. A flux only CBM cannot verify concentration, pH or volume admissibility just because its amount balances close. Its result should be called sustained flux feasibility at stated direction assumptions, not full validation of WT physiology.

## Reference state is not a capacity measurement

The preserved WT equilibrium has intracellular Na18.0826 mM, K110.2372 mM, Cl53.2985 mM, pH7.05039, TIC8.26799 mM and volume1.38988 pL. Its NKCC cycle rate is0.1493911 fmol/s, AE4 rate0.07507245 fmol/s, NBC rate0.07168031 fmol/s, NHE rate0.00603620 fmol/s, pump rate0.06319047 fmol/s, signed AE2 chloride flux−0.00074806 fmol/s and CaCC export0.37310660 fmol/s. These are saved model observables. They may define a dimensional plotting scale or documented reference direction, but they are neither measurements nor upper capacity bounds. No trajectories or equilibrium solves were rerun.

## Conservative computational bounds

The orchestrator's proposed output coordinate J with WT reference J=1 is a suitable mathematical normalisation. Restricting J≤1 truncates the reporting region; it does not establish a physical output capacity. Report an attainable KO value1 as “the WT reference output is feasible”. The absolute maximum and possible output beyond1 remain unidentified.

Shared unknown flux caps M=100, with one predeclared M=1000 check, may be used as computational bounds. Every such entry must carry that classification. A finite jointly feasible rescue witness lying strictly below M demonstrates the same feasibility for all larger M. FVA upper endpoints that grow with M must remain labelled computational; invariance of every unconstrained endpoint should not be claimed. Exact identities and finite witness arguments can support the main structural conclusions even where absolute maxima are unbounded without the artificial caps.

For each predeclared single resource test, minimising the shared resource needed for WT J=1 supplies a conditional algebraic threshold. Applying that threshold to KO is a structural diagnostic. It becomes a defensible biological restriction only if independent evidence supports the required value. The normalisation or old WT reference flux must not be used to manufacture that evidence.

No LP, dynamics, source model, vendor core or canonical repository file was changed. No Git writes or CarbonScope access occurred.
