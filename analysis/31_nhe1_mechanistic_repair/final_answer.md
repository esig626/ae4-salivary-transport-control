# Task 31: mechanistic NHE1 implementation and bounded REST outcome

Cha eight-state/Mod2 NHE1 is implemented with the exact published final Table S1 constants. The R09 WT-only scalar calibration converged. Neither exact-AE4-null attempt reached numerical REST closure within its predeclared limit; no dynamics were run.

## 1. Was the mechanistic model implemented correctly?

Yes, as the printed parameterized model: obligatory 1:1 Na inward/H outward exchange, shared forward/reverse cycle, intracellular proton modifier, finite rates, zero carbon and net charge source, and identical WT/null law. Independent rational equation reconstruction agrees to better than 1e-12 relative precision. The Task 30 Eq. 26 reference remains selectable.

The source has two unresolved precision inconsistencies. Printed k2−=183 ms^-1 differs by about 0.33% from the Eq. 6 value from its other rounded entries (183.60748); the resulting chemical-equilibrium zero shifts by 0.001439 pH units. The reconstructed source turnover is 2.52121 ms^-1, 0.8485% above the stated 2.5 upper restriction. All Table S1 entries are retained unchanged; exact thermodynamic consistency and exact reproduction of the source turnover restriction are not claimed.

## 2. What salivary density was obtained?

**2.339370005697548e-05 fmol of active carriers per cell**, equivalent to **14088.01546403 carriers**, or **0.0287569207 times** the Table S1 source-cell count N=489900. This is an effective carrier amount, not a measured membrane-area density. Five bounded scalar evaluations used only R09 WT pH, giving pH 6.9100002235 against target 6.91. The density was frozen before R10 WT and both null solves.

## 3–5. Did R09 WT, R09 exact null and R10 close?

| Background / genotype | pHi | Na_i (mM) | K_i (mM) | Cl_i (mM) | TIC / HCO3 (mM) | Volume (pL) | Max scaled residual | Outcome |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| R09 WT | 6.910000 | 11.6361 | 116.7632 | 60.3050 | 5.6717 / 4.9094 | 1.448426 | 5.104e-12 | Admissible; chloride context missed |
| R10 WT | 7.059426 | 7.8548 | 127.2206 | 66.8705 | 9.7817 / 8.8094 | 1.584337 | 1.699e-11 | Numerical root; pH context missed |
| R09 AE4_NULL | 7.062092 | 11.3693 | 119.2807 | 59.9588 | 11.6082 / 10.4606 | 1.563181 | 0.02215 | Unclosed candidate; local limit reached |
| R10 AE4_NULL | 7.139924 | 7.8143 | 128.3631 | 66.5973 | 13.5453 / 12.4051 | 1.670071 | 0.0131 | Unclosed candidate; local limit reached |

**R09 WT:** numerical, conservation/current, positivity, rank and capacity gates pass. pH, sodium and volume are admissible. Cl_i=60.305 mM misses the measured WT interval 48.6–51.6 mM; this was not a fitting target.

**R09 exact null:** no accepted REST state. Its final candidate has ordinary sodium and volume below 3 pL, but pH 7.0621 misses 6.87–6.91 and its scaled residual 0.02215 exceeds the 1e-7 gate. Candidate concentrations are not equilibrium predictions. All AE4 ion and charge fluxes are exactly zero.

**R10 with the frozen model:** WT closes numerically but pH 7.0594 misses 6.84–6.98; chloride also misses context. The exact null remains unclosed with scaled residual 0.01310 and pH 7.1399. No R10 or genotype recalibration was performed. Both null cases used their matched WT coordinates as their only starting point.

## 6. What happened dynamically?

No matched WT/null REST pair is admissible, so the requested conditional CCh+IPR/Ca=0.25 uM check was not started. Integrated secretion and null/WT secretion ratio are not available. No inference about the 30–35% secretion deficit is made.

## 7. What remains unresolved?

An admissible exact-null resting cell was not obtained under this bounded implementation test. Failure to converge within 2400 residual calls per null solve does not prove that a mathematical root is absent. Cross-cell-type transfer of Cha kinetics, the printed-precision discrepancies, and the WT chloride mismatch remain unresolved. No other transporter was changed and no additional search or model family was started.

## Verification and provenance

Compute: 5 scalar evaluations; 8 stationary calls; 8297 actual residual evaluations; 0 integrations; one worker and one BLAS thread; 2.309464 seconds of numerical execution. All global and per-solve limits were respected.

The previous 24 passing tests and one source-related skip remain in the original log without rerunning them. Seven new tests passed (three source tests and four continuation tests); the old source skip is now resolved. Completed calculations and accepted cases were not rerun. Only membranes.py (NHE dispatch) and parameters.py (NHE selector/amount validation) differ among the existing materialized files.

[Document S1 source](https://pmc.ncbi.nlm.nih.gov/articles/instance/2776256/bin/mmc1.pdf), Table S1, PDF page 5; 231386 bytes; SHA-256 `3fdff54677a5c3347045de9e3d1efd52b62cbc2455400a37db1595e9f2fc146c`. Exact values, units and provenance were recorded before production use in source_table_s1.json. The full source/equation note, calibration history, frozen parameters, per-case checkpoints, REST table, flux ledger, numerical log and verification JSON accompany this report.

Published through the connected GitHub contents integration to [codex/task-31-nhe1-mechanistic-repair](https://github.com/esig626/ae4-salivary-transport-control/tree/codex/task-31-nhe1-mechanistic-repair). All 28 source/test/result payload files were verified against local text at commit `0f4bc85c0bc09054490886defda45fc0c1827be6`, before adding this final report and verification record. All other prepared-tree files are unchanged. One checkpoint upload required retries after a timeout and an internal integration error; it succeeded without scientific recomputation. No merge to main.
