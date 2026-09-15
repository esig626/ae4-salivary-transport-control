# Carbon, alkalinity and charge review for 46H00

Reviewed only pinned sources, with Git blob checks; no dynamic rerun or scientific optimisation. Source pins and exact locations are in `worker_carbon.json`.

The parent steady cell balances, writing N for NKCC1, H for NHE1, E for AE2, B for inward 1Na:2HCO3 NBC, A for equal cation AE4, P for total pump, K for total K exit and C for apical Cl exit, are:

| Balance | Exact equation |
|---|---|
| Na | N + H + B − A/2 − 3P = 0 |
| K | N − A/2 + 2P − K = 0 |
| Cl | 2N + E + A − C = 0 |
| TIC | −E + 2B − 2A + Db + Da = 0 |
| TA | H − E + 2B − 2A = 0 |

Thus 2A = H + 2B − E, Db + Da = H, C = 6P − H, and K = 5P − H − B. These are exact stationary identities for this architecture, independent of inherited kinetic laws. They do not prove pH, volume, passive direction or saliva flow.

Db is signed CO2 bath to cell; Da is signed CO2 lumen to cell. CO2 carries carbon and no alkalinity. Forcing both inward is unjustified without concentration evidence. NHE proton extrusion adds TA without TIC. The finite buffer and closed acid/base reactions provide no sustained source of either conserved quantity. No free acid or TA sink exists in the parent architecture.

With jp positive lumen to bath, the five lumen sources before outflow are 3Pa − jpNa; −2Pa + Ka − jpK; C − jpCl; −jpHCO3 − Da; and −jpHCO3 for Na, K, Cl, TIC and TA respectively. Subtract the corresponding solute outflow from each. There is no transcellular apical HCO3 channel in this parent: lumen HCO3 is paracellular. Aggregating pump/K for cell equations is valid, but lumen equations still need the apical allocation Pa and Ka. Neutral outflow requires outNa + outK − outCl − outTA = 0. Lumen carbon bookkeeping additionally gives outTIC − outTA = −Da. TIC and TA outflows may be retained as conserved coordinates; unconstrained independent values must not be claimed to imply realistic carbonate composition.

TIC is CO2 + HCO3 + CO3. TA is HCO3 + 2CO3 + Bminus + OH − H. Cell bulk charge is Na + K − Cl − TA − fixed anion; the fixed anion has zero time derivative. Per cycle, NBC, pump and K exit supply −1 cell charge equivalent; Cl exit supplies +1. The other retained cell transporters and neutral CO2 have zero net charge. Current closure therefore agrees with C − P − K − B = 0. Current neutrality should be audited as a consequence or constraint of the entire transport system, not imposed on individual electrogenic columns.

Conditional loading/export directions are defensible as a declared screening sector. They are not global thermodynamic conclusions. The NBC affinity includes Vb/VT as well as ln(Nao HCO3o² / (Nai HCO3i²)); other transporter directions likewise depend on the actual state. A reference state sign cannot validate every point in an LP region. Equal AE4 routing is a retained model construction, not established microscopic stoichiometry.

The pinned CO2 permeabilities, AE2 capacity and pump/K localisation fractions are modelling values. NBC capacity is an algebraically constructed WT scale. The Task44 conditional NHE capacity inequality additionally uses inherited kinetics and pH >= 6.6; it is not an independent measurement. The [3,11] carbonate solver bracket is numerical, not physiological. Preserve these classifications when constructing bounds.

For the planned minimal lumen species bookkeeping, use nonnegative bicarbonate and CO2 outflows, set outTA = outHCO3 and outTIC = outHCO3 + outCO2, and explicitly neglect luminal carbonate, H/OH and the tiny buffer term. This avoids an independent carbon free alkalinity outlet. It induces Da = −outCO2, jpHCO3 = −outHCO3 and Db = H + outCO2. These signs follow the chosen chemistry/outflow sector and do not independently establish passive driving forces.
