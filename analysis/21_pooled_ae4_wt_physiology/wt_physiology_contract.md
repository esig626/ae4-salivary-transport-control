# Task 21 WT physiology contract — declared before optimization

The model is Task 20's unchanged `POOLED_CATION_112_NO_SLIP`, with one
reversible 1:1:2 working cycle and donor-side Na/K partition. All ten inherited
roots supply independent numerical starts and capacity references. Their
conserved states, Na/K values and volumes are not exact targets.

| Quantity | WT reference | Task 21 screening envelope | Evidence role |
|---|---:|---:|---|
| Intracellular chloride | 50.10 ± 1.50 mM | 47.10–53.10 mM | Measured WT target; ±2 reported SEM |
| Intracellular pH | 6.91 ± 0.07 | 6.77–7.05 | Measured WT target; ±2 reported SEM |
| Intracellular Na | 15.5–25 mM in cited reference contexts | 2–60 mM | Broad physiology/provenance envelope |
| Intracellular K | 120 mM lineage reference; 140–150 mM physiological discussion | 60–200 mM | Broad physiology/provenance envelope |
| Cell volume | 1.3 pL published-model convention | 0.3–5.0 pL | Broad geometry/provenance envelope |

The chloride and pH targets follow the explicit Task 21 directive and the
matched WT calibration entries in `analysis/13B_modern_full_model/wt_calibration.md`.
Their bands are modeling screens, not population confidence intervals.

The provenance of the remaining references is:

- `model/parameters.md`, the Table 1 source mapping, attributes Na 20 mM to
  Grinstein and Foskett (1990); its 25 mM value is a model output. Neither is
  treated as a direct native mouse-SMG uncertainty interval.
- Peña-Münzenmayer et al. (2016), [JGP 147:423–436](https://doi.org/10.1085/jgp.201611571),
  uses Na 15.5 mM calculated/transferred from rat sublingual acini in its
  thermodynamic illustration and discusses intracellular K around 140–150 mM.
  These contextual numbers are not new measurements in the Task 21 preparation.
- `model/parameters.md` attributes the K 120 mM lineage reference to Pedersen
  and Petersen (1973). The source transfer and lack of a matched native-SMG
  uncertainty interval preclude treating 120 mM as an exact target.
- The same source map and `analysis/13B_modern_full_model/wt_calibration.md`
  identify 1.3 pL as a [Palk 2010](https://doi.org/10.1016/j.jtbi.2010.06.027)/2018
  modeling convention, not a direct native-SMG volume measurement.
- `analysis/13B_modern_full_model/evidence_ledger.md` labels historical lineage
  separately from primary biological evidence. The 2018 table is accessed here
  through that repository source mapping, not represented as a freshly verified
  original Na/K experiment.

No direct native-submandibular measurement interval for Na, K or single-cell
volume is established by these inputs. The selected endpoints are explicit,
broad modeling choices: Na spans substantially below the cited resting values
to severalfold above them; K spans half the lineage reference to a generous
margin above the physiological discussion; volume spans approximately
0.23–3.85 times the historical cell size. The volume envelope is wider than
the old numerical box and removes its narrow target. None is inferred from a
genotype result. These assumptions are frozen before optimization and will not
be revised in response to WT or genotype outcomes.

The full WT state may re-equilibrate. Fifteen existing capacities are adjustable:
AE4, NKCC1, AE2, NHE1, independent apical/basolateral pumps and K conductances,
CaCC, both neutral CO2 exchanges, and Na/K/Cl/HCO3 paracellular conductances.
Bath, chemistry, buffers, geometry, water/outflow parameters and equations,
regulation, source signs and stoichiometries remain inherited. There is no
chloride-loading share target or legacy AE4 flux match.

WT ranking is lexicographic: first normalized squared chloride/pH error, then
the largest absolute capacity log-fold departure, then summed squared log-fold
departure. No state-distance penalty or genotype quantity appears. Positive
pooled affinity and chloride loading, complete resting/current/conservation
closure, physical parameters, the declared physiology and all three production
WT dynamics are acceptance requirements. Duplicate states **and** calibrated
parameter vectors are identified at relative tolerance 1e-6; distinct immutable
backgrounds remain distinct. Every distinct admissible solution is retained.

This contract must be committed, pushed and remotely verified before inspecting
new optimization results. A second, complete WT checkpoint must include the
entire calibration ensemble and validated 600-s WT trajectories at calcium
0.10, 0.25 and 0.50 uM. It must be pushed and remotely verified before any Task
21 AE4-loss or AE2-loss calculation. No valid WT model means no genotype run.
