# Baseline model source map

## Scope and source policy

The primary reconstruction source is the repository copy of the published article:

> Vera-Sigüenza E, Catalán MA, Peña-Münzenmayer G, Melvin JE, and Sneyd J (2018), “A Mathematical Model Supports a Key Role for Ae4 (Slc4a9) in Salivary Gland Secretion,” *Bulletin of Mathematical Biology* 80(2):255–282. DOI: 10.1007/s11538-017-0370-6. Repository copy: `literature/core/2018_AE4_model.pdf`.

Locations below use the published section, equation, table, figure, or appendix number. The historical files in `archive/legacy-2017/` are not scientific sources for this map and were not used to supply equations, values, or implementation details.

Provenance labels used here:

- **Direct**: the 2018 paper cites an experimental or standards source for the quantity.
- **Inherited model**: the 2018 paper adopts or simplifies a published model.
- **Paper-derived**: the 2018 paper derives or calibrates the item; it is published, but not independently measured.
- **Unresolved**: the published description is internally inconsistent or incomplete enough to block a faithful implementation decision.

## State equations and constraints

| Required component | Published location in the 2018 paper | Supporting published source named there | Status |
|---|---|---|---|
| Intracellular Cl− amount balance | §2.3, Eq. (1); collected system in §2.7 | Flux sources listed below | Published equation |
| Intracellular Na+ amount balance | §2.3, Eq. (2); §2.7 | Flux sources listed below | Published equation |
| Intracellular K+ amount balance | §2.3, Eq. (3); §2.7 | Flux sources listed below | Published equation |
| Intracellular HCO3− amount balance | §2.3, Eq. (4); §2.7 | Buffer, Ae2, and Ae4 sources below | Published equation |
| Intracellular H+ amount balance | §2.3, Eq. (5); §2.7 | Buffer and Nhe1 sources below | Published equation |
| Intracellular CO2 amount balance | §2.3, Eq. (6); §2.7 | Appendix 8; Sharp et al. (2015), DOI 10.1016/j.jtbi.2015.06.050 | Phase 10 conservation audit selects the sum of the two printed inward fluxes; retain the opposite printed aggregate as a sensitivity convention |
| Basolateral membrane potential | §2.4, Eq. (7); QSS form in §2.9 | Two-potential observation attributed to Young (1968), DOI 10.1007/BF00596392 | Published differential equation; QSS used for simulation |
| Apical membrane potential | §2.4, Eq. (8); QSS form in §2.9 | Young (1968) | Published differential equation; QSS used for simulation |
| Cell-volume balance | §2.5, Eq. (9); collected system in §2.7 | Palk et al. (2010), DOI 10.1016/j.jtbi.2010.06.027 | **Resolved in Phase 10:** conservation and the upstream equation require `dω_i/dt=q_b-q_a`; both printed signs give the same steady zero set |
| Apical, basolateral, and paracellular water fluxes | §2.5, Eqs. (10)–(13) | AQP5 context: Delporte and Steinfeld (2006), DOI 10.1016/j.bbamem.2006.01.022; permeabilities: Palk et al. (2010) | Formulas retained; `q_a` is cell-to-lumen and `q_b` bath-to-cell |
| Luminal Na+ balance | §2.6, Eq. (14); §2.7 | Tight-junction model in Appendix 10 | Published equation |
| Luminal K+ balance | §2.6, Eq. (15); §2.7 | Tight-junction model in Appendix 10 | Published equation |
| Luminal Cl− balance | §2.6, Eq. (16); §2.7 | CaCC model in Appendix 4 | Published equation; printed loss term drops a subscript |
| Intracellular electroneutrality and impermeant charge | Appendix 9, Eqs. (32)–(33) | Derived in the 2018 paper | Paper-derived constraint |
| Luminal and interstitial electroneutrality | Appendix 9, Eqs. (34)–(35) | Derived in the 2018 paper | Paper-derived constraints |
| QSS membrane constraints | §2.9 | DAE reference: Kunkel and Mehrmann (2006) | Published simulation choice; solver tolerances absent |

The QSS formulation is described as ten differential equations plus two algebraic membrane-potential constraints. See `states.md` and `equations.md` for the resulting inventory and the unresolved canonicalization gate.

## Flux and constitutive models

| Flux or constitutive relation | Published location in the 2018 paper | Upstream published source | Provenance/status |
|---|---|---|---|
| Nkcc1 cotransporter | §2.2; Appendix 1, Eq. (17), Table 3 | Palk et al. (2010), DOI 10.1016/j.jtbi.2010.06.027, and its 2013 [corrigendum](https://www.sciencedirect.com/science/article/pii/S0022519312005668), simplifying Benjamin and Johnson (1997), DOI 10.1152/ajprenal.1997.273.3.F473 | Upstream corrigendum resolves coefficient units as M-based; the 2018 density remains calibrated and lacks a whole-cell conversion |
| NaK-ATPase | §2.2; Appendix 2, Eq. (18), Table 4 | Palk et al. (2010), simplifying Smith and Crampin (2004), DOI 10.1016/j.pbiomolbio.2004.01.010 | Inherited model; density calibrated in the 2018 paper |
| Ca2+-activated K+ channel | §2.2; Appendix 3, Eqs. (19)–(21) | Takahata et al. (2003), DOI 10.1152/ajpcell.00250.2002; one constant taken from Palk et al. (2010) | Inherited/modified model; maximum conductance calibrated |
| Ca2+-activated Cl− channel | §2.2; Appendix 4, Eqs. (22)–(24) | Simplification based on Takahata et al. (2003); comparison to Arreola et al. (2002), DOI 10.1113/jphysiol.2002.021980; physiological context from Frizzell and Hanrahan (2012), DOI 10.1101/cshperspect.a009563 | Paper-selected simplification; maximum conductance calibrated |
| Ae2 exchanger | §2.2; Appendix 5, Eq. (25), Table 5 | Falkenberg and Jakobsson (2010), DOI 10.1016/j.bpj.2009.11.045 | Inherited model; activity calibrated |
| Nhe1 exchanger | §2.2; Appendix 6, Eq. (26), Table 6 | Falkenberg and Jakobsson (2010) | Inherited model; activity calibrated |
| Ae4 exchanger | §2.2; Appendix 7, Fig. 10, Eqs. (27)–(31), Table 7 | Transport evidence and stoichiometry: Peña-Münzenmayer et al. (2016), DOI 10.1085/jgp.201611571; Markov-state method: Dupont et al. (2016) | Paper-derived model and calibration; units and intermediate-state equations need audit |
| CO2 diffusion and HCO3−/H+ buffer | §2.3; Appendix 8, Table 8 | Sharp et al. (2015) | Inherited mass-action form and parameters; net-flux sign needs audit |
| Tight-junction Na+ and K+ fluxes | §2.4; Appendix 10, Eqs. (36)–(37), Nernst relations, and `V_t = V_a - V_b` | Derived in the 2018 paper | Paper-derived; `I`/`J` and `G`/`g` notation is inconsistent |
| Osmotic water transport | §2.5, Eqs. (10)–(13); Table 2 | Palk et al. (2010) for permeabilities | Inherited parameters with paper-specific three-compartment equations |
| Convective lumen outflow | §2.6, Eqs. (14)–(16) | Defined in the 2018 paper | Paper-derived |
| Prescribed intracellular Ca2+ drive | §2.8, Fig. 2 | Foskett and Melvin (1989), DOI 10.1126/science.2500708; Soltoff et al. (1989), DOI 10.1085/jgp.93.2.285; Bruce et al. (2002), DOI 10.1074/jbc.M106609200 | Qualitative fixed input only; no machine-readable function is published |

## Parameter locations

| Parameter family | Published location | Source status |
|---|---|---|
| Baseline concentrations, pH, membrane potentials, volume, and fitted osmotic terms | Table 1 | Mixed direct measurements and paper-derived values; preserve row-level labels in `parameters.md` and `observables.md` |
| Channel and tight-junction conductances; compartment volumes; water permeabilities; physical constants | Appendix 11, Table 2 | Conductances paper-derived; volumes/permeabilities from Palk et al. (2010); `F` and `R` from standards papers |
| Nkcc1 density and rate parameters | Appendix 1, Table 3 | Density paper-derived; corrected rates from the Palk et al. (2013) corrigendum are M-based and require explicit conversion for mM states |
| NaK density and rate parameters | Appendix 2, Table 4 | Density paper-derived; rates from Palk et al. (2010) |
| Ae2 activity and half-saturation constants | Appendix 5, Table 5 | Activity paper-derived; constants from Falkenberg and Jakobsson (2010) |
| Nhe1 activity and half-saturation constants | Appendix 6, Table 6 | Activity paper-derived; constants from Falkenberg and Jakobsson (2010) |
| Ae4 density and rate constants | Appendix 7, Table 7 | All marked as determined from the model |
| Buffer and CO2-transport rates | Appendix 8, Table 8 | Sharp et al. (2015) |

## Observable locations

| Observable | Published definition or target |
|---|---|
| Total primary saliva flow | `q_tot = q_a + q_t`, §2.5 Eq. (13); control and knockout time courses in Fig. 8 |
| Cell volume | Dynamic quantity in Eq. (9); resting target in Table 1 |
| Intracellular and luminal ion concentrations | Amount/volume states in Eqs. (1)–(6) and concentration states in Eqs. (14)–(16); resting targets in Table 1; time courses in Figs. 3–6 |
| Apical and basolateral membrane potentials | Eqs. (7)–(8), QSS constraints in §2.9; resting targets in Table 1; time courses in Fig. 7 |
| Transporter and channel fluxes | Constitutive equations in Appendices 1–8 and 10; knockout comparisons in Fig. 9 |
| Compartment osmolarities | Ionic sums and impermeant terms in §2.5; reported resting totals of 292.6 mM (interstitium), 296.6 mM (cell), and 297.4 mM (lumen) |
| pH | Resting `pH_i`, `pH_l`, and `pH_e` rows in Table 1; conversion and dynamical treatment are not fully specified |

## Independent-reconstruction gate

The published paper supplies the model topology, the balance equations, constitutive forms, and nominal tables without using the historical archive. It does **not** yet supply an unambiguous executable specification. Phase 10 resolves the volume sign and the upstream NKCC coefficient units, but flux/current and whole-cell scaling, incomplete numerical Ca2+ input, acid/base closure, model-derived calibration values, and missing numerical-method details remain. See `analysis/00_inventory/open_questions.md` and `analysis/10_identifiability_discrimination/reduction.md`.
