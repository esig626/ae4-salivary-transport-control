# Fixed structural and protocol inputs

This note makes explicit the fixed multipliers and stoichiometric constants alongside `parameter_provenance.csv`. It adds no parameter or correction after the prediction freeze.

| Input | Value or role | Classification |
| --- | --- | --- |
| NKCC expression | 1 in every case | Inherited protocol setting; no measured absolute abundance |
| NHE1 and AE2 expression | 1 in every case | Inherited protocol settings |
| AE4 expression | WT 1, null 0, deferred condition 0.05 | Defined interventions; not optimised values |
| NKCC transported Na:K:Cl | 1:1:2 inward | Literature derived transporter stoichiometry |
| Pump transported Na:K | 3 outward, 2 inward per cycle | Literature derived transporter stoichiometry |
| Pump Na and K saturation exponents | 3 and 2 | Inherited effective kinetic law; stoichiometry alone does not measure these Hill exponents |
| Equal AE4 cation export fractions | 0.5 Na, 0.5 K per cycle | Frozen Task 40 source assumption; not an independently measured salivary split |
| Intracellular concentrations | Amount divided by evolving volume | Derived state quantities, not fixed parameters |
| Voltage capacitance | No capacitive state | Inherited algebraic electrical approximation |
| Upper and lower voltage root bracket | Plus or minus 0.5 V | Numerical control only; no physiological capacity inference |
| Current conversion | `1 fmol/s = F*1e-15 A` for one outward positive charge | Dimensional conversion using the Faraday constant |
| Generic NKCC capacity | 0.32 fmol/s | Inactive legacy field under the selected Palk law |

The article conclusion remains conditional on the active effective kinetics. Neither the parameter ledger nor this complement supplies a measured uncertainty interval or authorises selecting a cap from the secretion phenotype.
