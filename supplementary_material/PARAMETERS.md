# Principal model parameter record

This document identifies the parameters used by `model/task52_model.py`. It is not a new calibration. The numerical implementation and the original JSON records retain their source precision. Values below are displayed compactly and must not replace the full precision inputs when reproducing a run.

## Where the full precision values are defined

| Quantity | Source inside this folder |
| --- | --- |
| Basic cell, bath, geometry, membrane and water parameters | `model/frozen_task37/src/modern_full_model/parameters.py` |
| Selected reference parameter construction | `model/frozen_task37/src/modern_full_model/task30_nhe1_repair.py` and `reference/native_wt_contract.json` within the same `frozen_task37` folder |
| NHE1 carrier amount and principal model assembly | `model/task52_model.py` and `model/frozen_task37/src/modern_full_model/task31_nhe1_repair.py` |
| NHE1 transition and binding constants | `model/frozen_task37/src/modern_full_model/nhe1_cha2009.py` |
| Ae4 carrier law and donor allocation | `model/frozen_task37/src/modern_full_model/ae4_routing_only.py` |
| NBC capacity and stimulated supply | `model/frozen_task37/src/modern_full_model/nbc_minimal.py` |
| NKCC1 stimulated recruitment | `model/frozen_task37/src/modern_full_model/nkcc_stimulation.py` |
| Ae4 activation variable | `model/frozen_task37/src/modern_full_model/camp_pka.py` |
| Complete projected onset states and simulation settings | `parameters/parameter_and_case_freeze.json` |

Some effective parameters are constructed from the reference inputs rather than listed as a single flat JSON vector. Run `scripts/export_parameters.py` to export the actual constructed parameter objects. The exporter also includes the central onset vectors and case definitions. The reference resting checkpoint is used by the constructor; it is not the initial state used for the measured onset simulations.

## Physical and acid base quantities

| Parameter | Displayed value | Units |
| --- | ---: | --- |
| Faraday constant | 96485.33 | C/mol |
| Gas constant | 8.31 | J/(mol K) |
| Absolute temperature | 310.15 | K |
| pK1 | 6.10 | dimensionless |
| pK2 | 10.30 | dimensionless |
| pKw | 14.00 | dimensionless |
| Buffer pK in cell and lumen | 7.00 | dimensionless |
| Bath sodium | 145.00 | mM |
| Bath potassium | 5.00 | mM |
| Bath chloride | 126.16 | mM |
| Bath total inorganic carbon | 25.00 | mM |
| Bath pH | 7.40 | dimensionless |
| Fixed cell charge | 78.06 | fmol equivalents |
| Other cell impermeant osmotic particles | 124.70 | fmol |
| Cell buffer sites | 30.00 | fmol |
| Luminal buffer sites | 1.00e-9 | fmol |

Concentrations in the transport equations are in mM, whereas carbonate dissociation constants are used with proton concentration in molar units when calculating the species fractions. Prescribed calcium is in micromolar units, not mM. The fixed cell charge is not an additional osmotic amount.

## Transport and water quantities

| Parameter | Displayed value | Units |
| --- | ---: | --- |
| NKCC1 cycle capacity G_N | 0.32 | fmol/s |
| Ae2 cycle capacity G_2 | 5.00e-3 | fmol/s |
| NKCC1 and Ae2 widths | 2.00 | dimensionless |
| NBC capacity G_B | 0.12 | fmol/s |
| NBC width | 2.00 | dimensionless |
| Basolateral CO2 permeability coefficient | 0.02 | fmol/(s mM) |
| Apical CO2 permeability coefficient | 0.01 | fmol/(s mM) |
| NHE1 carrier amount | 2.34e-5 | fmol |
| Ae4 carrier amount | 0.15 | fmol |
| Ae4 chloride attempt rate | 1.00 | 1/s |
| Ae4 sodium attempt rate | 0.10 | 1/s |
| Ae4 potassium attempt rate | 1.90 | 1/s |
| Ae4 activation time constant | 30.00 | s |
| Ae4 activation gain | 0.25 | dimensionless |
| Total potassium conductance | 14.00 | nS |
| Native apical chloride conductance | 31.40 | nS |
| Shared auxiliary chloride conductance | 0.00, 2.32, 4.49 | nS |
| Apical potassium conductance fraction | 0.30 | dimensionless |
| Apical pump fraction | 7.51e-2 | dimensionless |
| Total pump cycle capacity | 0.08 | fmol/s |
| Pump sodium half activation | 10.00 | mM |
| Pump potassium half activation | 1.50 | mM |
| Channel calcium half activation | 0.26 | micromolar |
| Channel calcium exponent | 1.46 | dimensionless |
| Paracellular sodium conductance | 0.30 | nS |
| Paracellular potassium conductance | 0.05 | nS |
| Paracellular chloride conductance | 0.25 | nS |
| Paracellular bicarbonate conductance | 0.05 | nS |
| Apical water coefficient | 4.32e-3 | pL/(s mOsm/L) |
| Basolateral water coefficient | 5.15e-2 | pL/(s mOsm/L) |
| Paracellular water coefficient | 2.60e-4 | pL/(s mOsm/L) |
| Dead luminal volume | 0.10 | pL |
| Luminal outflow rate constant | 0.20 | 1/s |

NHE1 transition constants k1+, k1-, k2+ and k2- are 10.50, 2.01e-1, 15.80 and 183.00 per millisecond. The numerical NHE1 law includes the factor 1000 converting turnover to per second. Its sodium binding constants are 195.00 and 16.20 mM outside and inside; proton binding constants are 1.62e-3 and 6.05e-4 mM; modifier constants are 4.80e-7 and 3.07e-5 mM outside and inside. The internal modifier exponent is 3.

Conductances must be converted from nS to S before multiplication by a membrane potential in volts. Positive transporter cycle rates and ionic fluxes use the sign conventions in the model source, not the sign of a negatively charged ion's conventional current.

## Initial states and stimulation

The central onset construction imposes intracellular chloride 50.10 mM and pH 6.91 in the control, and 36.50 mM and pH 6.89 in the knockout. Cell sodium concentration and volume are retained from the reference construction. Potassium, inorganic carbon and alkalinity follow from the charge, osmotic and chemical equations. The entire vectors, not only these two measured quantities, are provided in `parameters/parameter_and_case_freeze.json` under `projected_onsets.central`.

The principal protocol lasts 600.00 s. Calcium changes from 5.80e-2 to 0.25 micromolar and the beta input from 0 to 1. Ae4 expression is 1 in the control and 0 in the knockout. Full NKCC1 and NHE1 multipliers are 1.75 and 2.30. Both complete NKCC1 and Ae2 cycles are imposed from the control in the knockout, with multiplier 1 for all three central cases.

The retained source permits multipliers 1, 1.05 and 1.10. It does not implement the later continuous range through 1.49 and must not be silently changed to claim reproduction of that continuation.
