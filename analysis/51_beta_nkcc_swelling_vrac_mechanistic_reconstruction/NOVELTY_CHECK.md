# Task 51 novelty and no-repeat audit

## Controlling novelty claim

Task 51 is **not** a new VRAC search and **not** a new NKCC-core search.

Task 49 already implemented and conservation-tested the downstream law

`beta * positive swelling -> VRAC-like apical chloride current`

but deliberately kept beta-dependent NKCC regulation frozen. In AE4 knockout under IPR-only stimulation, the inherited core therefore had no beta-responsive chemical source, the exact genotype rest remained chemically static, positive swelling never began and the VRAC-like gate remained off. That exact self-start failure is closed and must not be rediscovered.

Task 51 adds exactly one source-backed upstream regulatory fact that changes the causal system:

`beta/cAMP -> increased NKCC1 activity`.

This permits IPR stimulation to alter Na/K/Cl loading and volume even when AE4 is absent, allowing the already tested swelling-sensitive apical branch to self-start without an AE4-expression multiplier.

## Independent source basis

Detailed citations and quantitative qualifications are in `docs/TASK51_BETA_NKCC_SWELLING_VRAC_EVIDENCE.md`.

Key independent sources:

- Catalán et al. 2015, adult mouse SMG, DOI `10.1073/pnas.1415739112`: TMEM16A-independent IPR secretion, IPR-induced swelling, DCPIB/NPPB-sensitive secretion and a VRAC-like chloride conductance.
- Turner et al. 1992, rat parotid, DOI `10.1172/JCI115695`: beta1/cAMP/PKA-associated NKCC activation, approximately 3-fold functional increase, K1/2 about 21.5 nM.
- Tanimura et al. 1995, DOI `10.1074/jbc.270.42.25252`: IPR-regulated NKCC phosphorylation, half-maximal effect about 20 nM and strong regulated activity recruitment.
- Kurihara et al. 1999/2002, DOI `10.1152/ajpcell.1999.277.6.C1184` and `10.1152/ajpcell.00352.2001`: beta/cAMP/PKA-associated salivary NKCC recruitment/phosphorylation.
- Rat submandibular acini PMID `9880083`: IPR increases measured NKCC-mediated NH4 influx about 2.5-fold and forskolin mimics the effect.
- Peña-Münzenmayer et al. 2021, DOI `10.1152/ajpgi.00145.2021`: beta/PKA activation of AE4 in mouse SMG; this arm is already represented and is not new in Task 51.

## Why this does not reopen excluded axes

Task 51 does not:

- replace or refit the Palk/Benjamin NKCC concentration-response core;
- alter NKCC 1Na:1K:2Cl stoichiometry;
- impose a hard NKCC capacity ceiling;
- search NKCC transporter families or regulatory mechanism families;
- invent another positive-swelling gate;
- add AE4-dependent TMEM16A recruitment;
- revisit calcium amplitude, NHE, NBC, pump, K handling, AE4 source stoichiometry or Na/K routing;
- fit timing;
- use the Task 50 multiplier in the mechanistic equations.

The only permitted new mechanism is the source-backed beta contribution to the existing NKCC activity regulator, upstream of the already tested swelling-sensitive current.

## Species qualification

Direct beta-NKCC functional gains are mainly rat parotid/submandibular measurements. They establish existence, sign and plausible magnitude, not an exact mouse-SMG constant. Mouse-SMG Catalán swelling/current data have priority for identifying the effective mouse strength wherever a defensible mapping exists.

## Benchmark separation

Task 50 is frozen separately at `archive/task-50-working-effective-coupling-benchmark` and specified by `docs/TASK50_EFFECTIVE_COUPLING_BENCHMARK.md`.

Task 50 is a **TARGET-CALIBRATED CONSTRUCTION** and provides the phenotype benchmark only. Task 51 must not inherit its multiplier or fit to its output.

## Held-out AE4 information

Before the immutable Task 51 prediction checkpoint, do not use to choose a parameter:

- the approximately `35 +/- 4.7%` AE4-null 10-minute secretion deficit;
- AE4-null CCh+IPR uptake `0.90 +/- 0.09 x10^-3 s^-1`;
- AE4-null IPR-only uptake `0.20 +/- 0.03 x10^-3 s^-1`;
- Task 50's `30.2612%` null and `23.1634%` 5%-AE4 deficits.

Those values are revealed only after the independently constrained mechanistic model is frozen and remotely published.

## Classification before execution

The literature chain is **SOURCE FACT + PREDECLARED HYPOTHESIS**.

Task 51 becomes a **NUMERICAL PREDICTION** of the held-out AE4 phenotype only if all new parameters are fixed/bounded from independent beta/NKCC/swelling/current information before phenotype reveal.
