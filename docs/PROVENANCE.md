# Provenance record

## Historical material

The repository was initialised with a directory labelled `M. Golubitsky` containing an older AE4 analysis project. It includes an unpublished manuscript and associated build products, figures, MATLAB files, and a copy of the published AE4 model paper.

On 2026-08-27 the directory was moved without content modification to `archive/legacy-2017/` to separate historical collaborative material from the new work. The original state remains available in Git history.

A second repository pointer to the published 2018 AE4 paper is stored at `literature/core/2018_AE4_model.pdf` because that published paper is an appropriate scientific source for independent baseline reconstruction.

## Use policy

Historical unpublished material may be inspected to understand what was attempted and to prevent accidental duplication of old reasoning. It is not a source for new sole authored text or code.

The new implementation must be derived from published model equations and documented published evidence.

## Consultation log

### 2026-08-27 — Phase 00 inventory and source map

- Consulted `archive/README.md` to confirm the archive’s immutable historical status and its broad content categories.
- Consulted Git tree metadata (paths, file formats, and sizes) for `archive/legacy-2017/` to produce a complete file-level inventory. This metadata inspection did not open any archived file payload.
- Did **not** consult the contents of `Ae4_Basis.tex`, `Marty.tex`, `mybib.bib`, any historical MATLAB/MAT file, any archived figure, or any historical build output.
- Did **not** use the duplicate archived article. The scientific source was the separate published pointer `literature/core/2018_AE4_model.pdf` (cross-checked against the same open published article at PMCID PMC5792321).
- Purpose: classify historical holdings, prevent accidental reuse, and establish that the new `model/` specification depends only on published sources.
- Independently produced: state inventory, balance-equation map, flux/source map, parameter provenance table, observable definitions, and a list of published-record ambiguities. No Phase 01 code or manuscript prose was created.

### 2026-08-27 — Phase 10 identifiability and discrimination audit

- Consulted the open published article by Vera-Sigüenza et al., *Bulletin of
  Mathematical Biology* 80 (2018), DOI
  `10.1007/s11538-017-0370-6`, PMCID `PMC5792321`, as the primary equation and
  parameter source. The public PMC HTML was used because the repository PDF
  pointer was not locally retrievable in this environment.
- Consulted the Phase 00 `model/` specification, inventory, open questions,
  provenance record, and decision log as established inputs.
- Consulted the published upstream Palk et al. salivary model (2010), DOI
  `10.1016/j.jtbi.2010.06.027`, and its 2013 corrigendum to resolve the volume
  sign and NKCC1 coefficient units. Consulted Sharp et al. (2015), DOI
  `10.1016/j.jtbi.2015.06.050`, for the individually defined inward CO2
  fluxes, and Peña-Münzenmayer et al. (2016), DOI
  `10.1085/jgp.201611571`, for the experimental AE4 transport context.
- Did **not** open or inspect the payload of any file under `archive/`, and did
  not use unpublished historical prose, equations, code, data, or figures.
- Independently derived the steady balance identities, regular-rank results,
  two-condition criterion, and diagnostic code from the published equations.
  The implementation is deliberately narrower than a full simulator because
  the printed model does not specify a reproducible numerical baseline.
- The earlier affine diagnostic already present on the working branch was
  retained as an explicitly uncalibrated software fixture. It was not used as
  physiological evidence, and its activity domain was narrowed to normalized
  knockout-to-WT values.
- For the targeted novelty audit, consulted only public scholarly records and
  primary publications listed with DOI links in
  `analysis/10_identifiability_discrimination/novelty_audit.md`. No text or
  code was copied from those publications.

### 2026-08-27 — Phase 11 full forensic reconstruction

- Accounted for every one of the 67 files under `archive/legacy-2017/` and
  inspected every locally available payload as immutable provenance or
  implementation evidence; the next item records the one connector-omitted
  payload and its public-content substitute.
  The complete per-file role, dependency, chronology, confidence, and integrity
  record is `analysis/11_forensic_reconstruction/archive_audit.md`; its
  machine-readable integrity companion is
  `results/11_forensic_reconstruction/archive_hashes.csv`.
- Consulted the duplicate archived published article as a remote-Git identity
  and inspected its scientific content through the canonical open article at
  PMCID `PMC5792321`. The connector omitted the 1,108,802-byte PDF payload, so
  the authoritative Git blob and public full text—not the zero-byte transport
  placeholder—were used as evidence.
- Consulted all seven MATLAB/MAT scientific files:
  `Ae4_Dynamics_Project/Original/Parameters.m`,
  `Original/Saliva_Ae4.m`, `Par.m`, `Par.mat`, `Salivary.m`,
  `Salivary2.m`, and `Salivary_ex.m`. Their roles are, respectively, raw
  calibration parameters, baseline calibration/driver, serialized-parameter
  script, binary parameter handoff, seven-state reduced engine, five-state
  fixed-lumen fork, and exploratory execution/plot driver.
- Loaded `Par.mat` non-destructively with SciPy. It contains only one `1 x 1`
  `par` struct with 45 scalar-double fields; every value exactly matches the 46
  ordered assignments to 45 distinct fields in `Par.m`. It contains no stored
  state trajectory, solver history, activity sweep, or knockout result.
- Consulted `Ae4_Basis.tex`, its compiled 19-page `Ae4_Basis.pdf`, `Marty.tex`,
  and `mybib.bib` as unpublished manuscript/provenance evidence. Their role is
  a Golubitsky/Sneyd/Vera-Sigüenza mathematical reduction and continuation
  analysis of the secretion-model family, not a published implementation or a
  source of new manuscript prose.
- Visually inspected all 20 `Figure1.eps` through `Figure20.eps`, all 20
  corresponding `FigureN-eps-converted-to.pdf` files, and all six standalone
  `Figure3.pdf` through `Figure8.pdf` files. The EPS/PDF pairs are 2017
  continuation and parameter-sensitivity assets for `Ae4_Basis`; the
  standalone PDFs are alternate/later `Marty.tex` panels. None is treated as a
  duplicate of a published 2018 result merely because of visual similarity.
- Parsed `Ae4_Basis.aux`, `Ae4_Basis.log`, and `Ae4_Basis.synctex.gz` for
  cross-references, exact figure dependencies, build chronology, and the
  historical `/Documents/Work/Marty/` source path. Inspected `mchem.sty`,
  `spbasic.bst`, `spmpsci.bst`, `spphys.bst`, `svglov3.clo`, and `svjour3.cls`
  to distinguish stock build support from scientific content. The connector-
  normalized `.log` is identified by its authoritative remote Git blob.
- Re-evaluated archived arithmetic with the labelled forensic tools in
  `analysis/11_forensic_reconstruction/tools/`. Neither MATLAB nor Octave was
  available, and the archived Original driver has a stale dependency
  signature, so no native-execution claim is made. Exact calibration checks are
  separated from translated SciPy-BDF trajectory evidence.
- Consulted, as public evidence independent of the archive, Vera-Sigüenza et
  al. (2018), DOI `10.1007/s11538-017-0370-6`; Palk et al. (2010), DOI
  `10.1016/j.jtbi.2010.06.027`, and its corrigendum, DOI
  `10.1016/j.jtbi.2012.10.027`; Benjamin and Johnson (1997), DOI
  `10.1152/ajprenal.1997.273.3.F473`; Smith and Crampin (2004), DOI
  `10.1016/j.pbiomolbio.2004.01.010`; Falkenberg and Jakobsson (2010), DOI
  `10.1016/j.bpj.2009.11.045`; Sharp et al. (2015), DOI
  `10.1016/j.jtbi.2015.06.050`; and Pena-Munzenmayer et al. (2016), DOI
  `10.1085/jgp.201611571`.
- Cross-checked later primary model-lineage records from 2019, 2020, 2021, and
  2022, as listed with direct links in
  `analysis/11_forensic_reconstruction/literature_crosscheck.md`. They confirm
  reuse and evolution of the model family and independently corroborate the
  volume sign and corrected NKCC concentration convention; they are not used as
  retrospective replacements for the 2018 equations.
- Searched the public article/correction/supplement records and indexed code
  records by DOI, title, and distinctive constants. No corrected 2018
  supplement or public source-code release was located as of 2026-08-27. This
  is recorded as a negative search, not proof that private or unindexed code
  does not exist.
- Independently produced the file ledger, dependency/lineage analysis,
  equation/parameter concordance, translated execution harness, blocker
  reassessment, and identity decision. No unpublished prose or historical
  source code was ported into a manuscript or a new scientific implementation.
  Nothing under `archive/` was modified.

### 2026-08-27 — Phase 12 AE4 mechanism reconstruction

- Read all completed Phase 00, Phase 10, and Phase 11 reports, the canonical
  `model/` specification, and the project provenance, decision, and results
  ledgers before defining the Task-12 comparison. The complete Phase-11
  67-file historical inventory and integrity evidence were treated as the
  controlling archive map.
- Reconsulted the immutable historical scientific implementation evidence
  needed to define the smallest internally closed comparison chassis:
  `Original/Parameters.m`, `Original/Saliva_Ae4.m`, `Par.m`, `Par.mat`,
  `Salivary.m`, `Salivary2.m`, and `Salivary_ex.m`. Historical expressions were
  independently rederived and implemented outside `archive/`; unpublished
  prose or source was not copied into a manuscript or attributed to one
  author.
- Verified the archive before and during the investigation with the Phase-11
  file ledger and tests. The aggregate local SHA-256 checkpoint over every
  materialized `archive/` file remained
  `e92e5a586d3f46d6823d19571558285923a9c1539b94afcf83ab1f48cb70715f`.
- Consulted the primary knockout study, Peña-Münzenmayer et al. (2015),
  [DOI 10.1074/jbc.M114.612895](https://doi.org/10.1074/jbc.M114.612895), for
  the saliva, resting chloride, pH, exchanger-activity, timing, AE2, NKCC1,
  channel, and NHE observations.
- Consulted the primary transport study, Peña-Münzenmayer et al. (2016),
  [DOI 10.1085/jgp.201611571](https://doi.org/10.1085/jgp.201611571), for
  direct Na and K movement, electroneutrality, reversibility, cation
  dose-response parameters, and the three explicitly considered
  stoichiometries. These measurements take precedence over the 2018 model law.
- Consulted the primary regulation study, Peña-Münzenmayer et al. (2021),
  [DOI 10.1152/ajpgi.00145.2021](https://doi.org/10.1152/ajpgi.00145.2021),
  for beta-adrenergic/cAMP/PKA activation, H89 inhibition, and the S173
  perturbation. Figure-read activation folds are explicitly distinguished
  from values stated in text and from whole-cell secretion data.
- Consulted Vera-Sigüenza et al. (2018), DOI
  `10.1007/s11538-017-0370-6`, only as a lower-priority published model record,
  not as ground truth for AE4 biology. Its pooled Na/K law and resting-model
  values were reconstructed as candidates or independent sensitivity targets
  and were allowed to fail.
- Independently produced the swappable AE4 implementation, C1–C8 candidate
  registry, exact stoichiometric/thermodynamic derivations, knockout-blind
  root and capacity screens, PKA diagnostics, adversarial tests, and final
  localization. Every equation and parameter in the Task-12 reports is
  labelled as primary experimental evidence, published model record,
  historical implementation evidence, or a new modelling decision.
- The approximately 35% knockout secretion decrement was held out from all
  calibration, capacity, root, and variant-selection functions. It was read
  only after candidates were frozen for validation. No file under `archive/`
  was modified.
- The moving-root continuation used a finite `1e-8`--`1e4` capacity domain.
  Its later +/-10% cell-volume acceptance band was a new, post-hoc sensitivity
  convention selected and screened without knockout outputs; it is not a
  reported experimental uncertainty or a global root-uniqueness claim.

### 2026-08-28 — Task 13B native-source supersession audit

- Read `AGENTS.md` and the complete Task-13B multi-agent prompt before the
  audit.  The AE4-null secretion holdout remained mechanically unavailable to
  source mapping, WT calibration, model generation, and the decision recorded
  here; no held-out plaintext was opened or used.
- Re-audited Gin et al. (2007), DOI `10.1016/j.jtbi.2007.04.021`, and Palk et
  al. (2010), DOI `10.1016/j.jtbi.2010.06.027`, for the published
  rat-parotid-lineage whole-cell channel map: `31.4 nS` maximum CaCC, `14 nS`
  total Ca-activated K, and the dimensionful Palk calcium gate with
  `K_Ca=0.26 uM`, exponent `1.46`.  These are cross-gland transfers with their
  assay lineage preserved, not new matched mouse-SMG measurements.
- Consulted Romanenko et al. (2010), DOI `10.1074/jbc.M109.068544`, for native
  mouse-SMG CaCC biophysics.  Its voltage/protocol-dependent approximately
  `0.126 uM` sensitivity is retained as a diagnostic and is not silently
  substituted for the Palk production gate.  The `15 +/- 2 pF` capacitance in
  Romanenko et al. (2008), DOI `10.1152/ajpgi.90384.2008`, is from salivary
  duct cells and was excluded from acinar whole-cell scaling.
- Re-audited Almassy et al. (2018), DOI `10.1007/s00424-018-2109-0`, and
  Poulsen and Bundgaard (1994), DOI `10.1007/BF00374318`.  The former licenses
  distributed pump topology, not a 50% capacity assignment; the latter's
  rat-parotid luminal/basolateral area densities yield an explicitly assumed
  uniform-density transfer of `0.075075`, with
  `0.056878/0.094586` conservative sensitivities.  The historical 30:70 pump
  and 40:60 K choices remain model assumptions rather than measurements.
- Consulted Kondo et al. (2015), DOI `10.1177/0022034515570943`, for the
  `66.91%` mouse-SMG acinar fraction and Kondo et al. (2019), DOI
  `10.1152/ajpcell.00421.2018`, for ipsilateral collection context.  A
  deliberately generous one-SMG construction gives an upper bound of
  `51,469,231` modeled `1.30 pL` acinar cells; this independently rejects the
  old G5 absolute scale.
- Consulted Catalan et al. (2015), DOI `10.1073/pnas.1415739112`, for
  beta-associated, volume-sensitive apical anion-exit evidence after adult
  acinar Tmem16A deletion.  The topology is required, but a quantitative
  whole-cell `gV` is not source licensed.  A WT-only diagnostic opening law
  was therefore allowed to fail rather than fitted as a rescue.
- Completed the WT-only native hydraulic reroot diagnostics without consulting
  a genotype secretion target.  Canonical H1 contains 45 panels/135 attempts/
  35 numerical roots; HW and HWQ each contain 225 panels/675 attempts/175
  numerical roots.  All three modes have zero WT-passing roots, for a combined
  495 panels, 1,485 attempts, and 385 numerical roots.
- Canonicalized mode identity explicitly: H1 is the production spelling;
  scale-one HW and HWQ are exact H1 alias controls rather than duplicated
  candidates.  Nonunit HW scales the three hydraulic coefficients, and
  nonunit HWQ also scales local outflow; every such row remains
  diagnostic-only.  At each tested H their displayed intracellular ranges
  coincide, a numerical observation that is not promoted to an analytic
  theorem.
- Preserved the old G5 gate, manifest, and hash table byte-for-byte and added a
  separate `g5_supersession.json` controlling record.  No file under
  `archive/` was modified, and no Task-13B final classification or reveal was
  made by this audit.

### 2026-08-28 — Task 13B absolute-source tier and engineering gate

- Completed the predeclared WT-only three-start production grid for
  `N_ABS_NKCC` and `AN_ABS_AE4_NHE` without consulting the protected genotype
  secretion targets.  The tier contains 630 panels and 1,890 attempts: 1,847
  report optimizer success, all 1,890 are admissible, and 981 meet the
  numerical-root criterion.  Within-panel clustering leaves 354 one-root
  panels and 276 rootless panels.
- Recorded the full numerical landscape rather than pruning at the first pass.
  N root counts at scales `0.5/2/4/5/6/7/8` are
  `40/35/33/31/30/30/30`; AN counts are `45/25/20/10/10/10/5`.  AN has no WT
  pass.  N has ten preliminary passes, all at scale 4: each of the five
  topologies at pump-capacity scale 1/`AE4NA05` and at pump-capacity scale
  2/`AE4NA20`.
- Preserved the WT-gate margins and numerical diagnostics for those ten rows.
  Cl spans `51.0172--53.0228 mM` with minimum margin `0.07723 mM`; pH spans
  `6.83860--7.04827` with minimum margin `0.001728`; Na and K span
  `14.3306--18.9814 mM` and `105.935--118.116 mM`.  Cell volume is `1.3 pL`,
  apical/basolateral voltages span `-36.1641` to `-34.5019 mV` and `-73.5186`
  to `-66.5896 mV`, and OTHER spans `110.024--126.008 fmol`.  Maximum scaled
  residuals are `1.01e-12--2.41e-11`; all normalized Jacobians have rank
  11/nullity 0 and no boundary hit.
- The combined H1/HW/HWQ/N/AN artifact set now contains 1,125 panels, 3,375
  attempts, and 739 numerical roots.  The controlling files are
  `native_source_panel.csv`, `native_source_root_attempts.csv`,
  `native_source_wt_roots.csv`, and `native_source_wt_summary.json`.  Their
  current state is a preliminary screen, not a frozen native model.
- Paused promotion before seven-start confirmation for an engineering audit of
  the hierarchy runner.  The resulting safeguards make panel reruns replace
  attempts and root clusters atomically, preserve both quartile starts in the
  three-/seven-start designs, reject duplicate IDs and fail-open eligibility,
  cover every licensed scale-one alias, require every upstream panel ID, and
  reject duplicate, failed/ineligible, or unconfirmed 33-start intake.  Tests
  in `tests/test_modern_native_source_panel.py` encode these execution
  invariants.  No equation, source scalar, WT interval, parameter bound, or
  scientific calibration was changed by this audit.
- All ten passes remain pre-seven-start and pre-33-start.  No root is eligible
  for native dynamics, independent final reproduction, a new pre-reveal
  manifest, or held-out reveal.  No protected target/validation content was
  accessed for this documentation update, and no file under `archive/` was
  modified.

## Future entries

Whenever archived material materially informs a research decision, append the date, file consulted, purpose, and what was independently rederived.
