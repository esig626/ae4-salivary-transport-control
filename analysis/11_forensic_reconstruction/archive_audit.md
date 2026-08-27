# Archive audit

## Scope, handling, and inspection method

This is a content-level audit of all 67 files below `archive/legacy-2017/`. The archive was treated as immutable. Text and source files were read directly; the MATLAB files were inspected line by line; `Par.mat` was opened non-destructively with SciPy; PDF metadata and text were extracted with Poppler; every page of `Ae4_Basis.pdf` and every distinct figure rendering was rendered to PNG and visually checked; EPS headers and embedded labels were inspected; and the LaTeX auxiliary, log, SyncTeX, bibliography, class, style, and bibliography-style files were parsed for build dependencies and chronology. Nothing below `archive/` was changed.

Full integrity values are repeated in the ledger and in `results/11_forensic_reconstruction/archive_hashes.csv`. There are two transport caveats:

- The archived published-paper PDF is a 1,108,802-byte remote blob (`04ba5733c02e2cbb0f4ed642b9a8944c0e4ba267`) identical in the remote tree to `literature/core/2018_AE4_model.pdf`. The connector omitted its payload, leaving a zero-byte local placeholder. The placeholder SHA-256 is **not** archive evidence. The article was inspected via the public PMC copy (PMCID `PMC5792321`).
- The authoritative `Ae4_Basis.log` is the 51,957-byte remote blob `aebe7db9e214d3c3c5a3c40d840dfcd13f57cb7f`. Connector text normalization produced a 51,961-byte local copy with SHA-256 `a7b30496edf7aa35229fee4763802e988d8794b2f91bd7992377ca99452bc36e`; that local digest is recorded only as transport evidence.

All other 65 files were available byte-exact locally and have authoritative SHA-256 values.

## Archive-level findings

1. The archive contains three scientifically distinct layers, not one coherent release:
   - the published 2018 article;
   - a 2017-2018, three-state mathematical reduction/continuation manuscript coauthored by Martin Golubitsky, with its figures and LaTeX build products; and
   - a separate full/reduced MATLAB dynamics project whose stored `Par.mat` was created on 2018-10-10.
2. `Ae4_Basis.tex` explicitly says it analyses the Sigüenza et al. secretion model, reproduces the full model balances, then freezes acid/base, lumen, volume, and voltage quantities to derive a three-state Na/K/Cl system. It is therefore an **analysis of** the 2018 model, not its published executable implementation. `Marty.tex` is a later/alternate fragment of the same reduction.
3. The MATLAB project is unquestionably in the same scientific model family: it contains the same compartments, seven dynamic quantities after algebraic elimination, NKCC1/NaK/Ae2/Ae4/NHE1/buffer fluxes, quasi-steady voltages, water flow, calcium switching, and knockout activity parameters. That structural relationship is high-confidence. Identity with the code that generated the 2018 publication is not established by the archive ledger alone.
4. The MATLAB set is not a self-contained historical release. `Original/Saliva_Ae4.m` loads `Par.mat`, calls `Salivary` with a signature inconsistent with the archived outer `Salivary.m`, and writes an absent `Conductances.mat`. `Salivary_ex.m` later calls `Salivary2` with too few arguments. Those are concrete signs of signature drift and missing dependencies.
5. `Par.mat` supplies a hard lower bound on chronology: it is a MATLAB v5 file created on 2018-10-10 and stores one `par` structure with 45 scalar fields. Thus the outer parameter bundle, as archived, postdates the online/print 2018 article; it may preserve, reorganise, or revise earlier values, but cannot by itself prove prepublication provenance.
6. The historical analysis figures do not duplicate the published article's figures. They plot continuation determinants, normalised states/flow, and one-at-a-time parameter sensitivities. Their direct link to 2018 is analytical: they use its reduced equations, calibrated values, and AE2/AE4 knockout interpretation.

## Dependency and chronology evidence

| Evidence | What it establishes |
|---|---|
| EPS `%%Creator`, `%%Title`, and `%%CreationDate` | All 20 EPS files were exported by MATLAB R2014b on 2017-06-09 from a `Desktop/New_Paper` working directory. |
| `Ae4_Basis.log`, PDF metadata, and `.aux` | A 19-page manuscript build completed on 2017-09-19 (PDF metadata: 2017-09-18 PDT), using all 20 EPS-derived PDFs, `mybib.bib`, `spbasic.bst`, `svjour3.cls`, and `svglov3.clo`. |
| `Ae4_Basis.synctex.gz` | The working source path was `/Users/esig526/Documents/Work/Marty/Ae4_Basis.tex`; this independently links the build to the Golubitsky/Marty project. |
| `Ae4_Basis.tex` author line | Names Elías Vera-Sigüenza, Martin Golubitsky, and James Sneyd; direct Golubitsky provenance. |
| `mybib.bib` | Contains 73 records and a still-year-2017 record for the AE4 paper, consistent with an online-first/pre-print-build interval. |
| Standalone `Figure3.pdf`-`Figure7.pdf` metadata | Re-rendered analysis panels created on 2018-01-08 and referenced by `Marty.tex`; these postdate the September manuscript build. |
| `Figure8.pdf` metadata/content | An additional direct-PDF MATLAB plot from June 2017, later used by `Marty.tex`; it is not the same content as `Figure8.eps`. |
| `Par.mat` header | MATLAB v5, platform `MACI64`, created 2018-10-10 18:18:55; definite post-publication parameter bundle. |

## Figure roles and duplicate/derived relationships

The role codes used in the ledger are content-based. Each EPS/PDF role was verified from rendered axes, formula labels, legends, and its `\includegraphics` position.

| Code | Scientific content |
|---|---|
| `R1` | `Det(Jx)` and `Det(Jy)` versus AE2 activity `G2`. |
| `R2` | `Det(Jx)` and `Det(Jy)` versus AE4 activity `G4`. |
| `R3` | Normalised fluid flow `Q/Q0` versus normalised `G2` and `G4`. |
| `R4` | Normalised states `x,y,z` versus `G4/G4_0`. |
| `R5` | Normalised states `x,y,z` versus `G2/G2_0`. |
| `R6` | `Q/Q0` versus activity for the alternate Hill-form model. |
| `S7` | Flow/activity response under `K_CaKC` variation (90%, 100%, 103%). |
| `S8`-`S9` | AE4 kinetic sensitivity to `k1` and `k2`. |
| `S10`-`S11` | AE2 sensitivity to `KCl` and `KB`. |
| `S12`-`S13` | NHE1 sensitivity to `KNa` and `KH`. |
| `S14`-`S17` | NKCC1 sensitivity to `a1`, `a2`, `a3`, and `a4`. |
| `S18`-`S19` | NaK sensitivity to `r` and `alpha1`. |
| `S20` | Flow/activity response under `K_CaCC` variation (90%, 100%, 103%). |
| `M3` | Later black-line replot of `R3`, used as panel A in `Marty.tex`. |
| `M4` | Later black-line replot of `R1`, used as panel E in `Marty.tex`. |
| `M5` | Later black-line replot of `R2`, used as panel F in `Marty.tex`. |
| `M6` | Later black-line replot of `R4`, used as panel D in `Marty.tex`. |
| `M7` | Later black-line replot of `R5`, used as panel C in `Marty.tex`. |
| `M8` | Normalised `dQ/dG4` versus AE4 coefficient `lambda11` at 50%, 100%, 150%, used as panel B in `Marty.tex`; not the content of `Figure8.eps`. |

For each number 1-20, `FigureN-eps-converted-to.pdf` is a direct `epstopdf` derivative of `FigureN.eps`; the build log records the source, output, conversion command, sizes, and inclusion lines. The six standalone PDFs are not byte/render conversions of those EPS files. Five (`M3`-`M7`) are later styled replots of the same numerical relationships as `R3`, `R1`, `R2`, `R4`, and `R5`; `M8` is a distinct lambda11 analysis.

## File-level ledger (67/67)

Relationship codes: `P` = the published 2018 source itself; `I` = full-model implementation evidence/candidate; `R` = explicit reduction or analysis of the published model; `B` = build/provenance only. Golubitsky codes: `D` = direct named/path/manuscript link; `A` = asset included by that manuscript; `N` = no direct link in the file. Confidence applies to the stated role and relationship, not to identity as the publication-generating code.

| # | Path | Type; bytes | Integrity evidence | Inspected content and scientific role | Dependencies | 2018 / Golubitsky | Confidence and basis |
|---:|---|---|---|---|---|---|---|
| 1 | `A Mathematical Model Supports a Key Role for Ae4 (Slc4a9) in Salivary Gland Secretion.pdf` | PDF; 1,108,802 | Git blob `04ba5733c02e2cbb0f4ed642b9a8944c0e4ba267` | Published article: full acinar-cell AE4 specification, parameter tables, baseline, calcium protocol, WT/AE2-KO/AE4-KO outputs. | None within archive; duplicate of canonical literature blob. | `P` / `N` | High: title, byte identity, and public PMCID. |
| 2 | `Ae4_Basis.aux` | LaTeX AUX; 8,530 | SHA-256 `bc11eefafbaded3c6b7e6f5c6bc3c758ece3449f49b68caa122a71f8339737f2` | Cross-reference/citation ledger: four figures, six tables, nine resolved references, repeated-label warnings visible through labels. | Generated from `Ae4_Basis.tex`; bibliography `mybib.bib`. | `B` / `D` | High: explicit generated labels and citations. |
| 3 | `Ae4_Basis.log` | pdfTeX log; 51,957 authoritative | Git blob `aebe7db9e214d3c3c5a3c40d840dfcd13f57cb7f` | Complete 2017-09-19 build trace; proves all 20 EPS conversions/inclusions, 19-page output, TeX Live 2014, and multiply-defined labels. | TeX/class/style/bib and all 20 EPS/PDF pairs. | `B` / `D` | High: tool-emitted dependency and timestamp record. |
| 4 | `Ae4_Basis.pdf` | PDF; 431,602 | SHA-256 `4a18dcb03f1af0cf820adfaad62e1f8713f3f060b287b4cf8d9f267bf4ad2d0a` | Rendered 19-page compiled unpublished manuscript: full model summary, three-state reduction, continuation, sensitivities, transporter appendix. | `Ae4_Basis.tex`, figures, bibliography, LaTeX support. | `R` / `D` | High: visually and textually inspected; named authors. |
| 5 | `Ae4_Basis.synctex.gz` | gzip SyncTeX; 152,846 | SHA-256 `42286857f58a813d74072b00bcf4ed9dfbf3f4819864f5d5c9de0f4c41fe08f2` | Source/PDF coordinate map; reveals `/Documents/Work/Marty/` source path and full TeX dependency list. | Generated by `Ae4_Basis.tex` build. | `B` / `D` | High: decompressed mapping inspected. |
| 6 | `Ae4_Basis.tex` | LaTeX source; 46,907 | SHA-256 `8b3eb8b7840978fcdb611dbc9dc7bc8f556d5bebac3363a3ded5f42d2a43e4ff` | Main Golubitsky manuscript source. Repeats full balances, then fixes acid/base/lumen/volume/voltages and derives three Na/K/Cl ODEs, 21 lambdas, continuation determinants, and sensitivity panels. | `mybib.bib`, `spbasic.bst`, `svjour3.cls`, `svglov3.clo`, figures 1-20. | `R` / `D` | High: explicit authorship, equations, includes, and references to Sigüenza 2017. |
| 7 | `Ae4_Dynamics_Project/Original/Parameters.m` | MATLAB script; 1,822 | SHA-256 `fa19bfabd63f0798bbbf72e269a4d469cf00b2b41434bc7d32f614163db20eef` | Full-model baseline/constants and calibration setup: geometry, channel gates, buffer, NKCC1, baths, states, voltages, exchanger kinetics, and deliberately rescaled water permeabilities (`b1*74.4`). | Run by `Original/Saliva_Ae4.m`. | `I` / `N` | High role; medium chronology: content is explicit, no embedded date. |
| 8 | `Ae4_Dynamics_Project/Original/Saliva_Ae4.m` | MATLAB script; 4,621 | SHA-256 `7d72a055ab240f658287643520558a592a618845945c01bbdc4cb454f84fdcf3` | Baseline-calibration driver. Solves conductances/activities from prescribed balance, forms 12 residuals, saves missing `Conductances.mat`, then attempts a 7-state `ode15s` run and plots states. | `Parameters.m`, `Par.mat`, a `Salivary` function; archived signature is incompatible. | `I` / `N` | High role; medium implementation identity because dependency drift prevents literal provenance claim. |
| 9 | `Ae4_Dynamics_Project/Par.m` | MATLAB script; 1,086 | SHA-256 `06ae901e1c045e234ec08d808bf909539b5d6d8cf009ebd98797791b7d313eff` | Serialises 45 scalar fields (one `kn` reassignment) for buffers, conductances, transporter activities, water coefficients, bath constants, and geometry into `Par.mat`. | Produces `Par.mat`; consumed by `Salivary_ex.m`. | `I` / `N` | High: assignments match the MAT structure field-for-field. |
| 10 | `Ae4_Dynamics_Project/Par.mat` | MATLAB v5 data; 773 | SHA-256 `fea3bdbe49dd794d43204bb234bd2b31e6215d784d920f46d3b4300f7068e1ca` | Hidden numerical state: one 1x1 `par` struct with 45 numeric scalar fields; materially stores calibrated conductances and activities. Header date 2018-10-10. | Produced by `Par.m`; loaded by both drivers. | `I` (post-publication bundle) / `N` | High: non-destructive field/type/value dump plus MAT header. |
| 11 | `Ae4_Dynamics_Project/Salivary.m` | MATLAB function; 2,317 | SHA-256 `a6524879c8d948d8be41cdc081c2423b1bda8510943ddd3ce76fab41dc87ed8e` | Seven-state full/reduced hybrid: lumen Na/K, height, cell Na/K/Cl/HCO3; eliminates H, CO2, lumen Cl and voltages; switches Ca after scaled `t=100`; includes AE2/AE4 parameters and multiplies derivatives by 10,000. | `par` struct; called correctly by `Salivary_ex.m`, incorrectly by `Original/Saliva_Ae4.m`. | `I` / `N` | High model-family role; medium publication identity. |
| 12 | `Ae4_Dynamics_Project/Salivary2.m` | MATLAB function; 2,160 | SHA-256 `94e8ddc12e26d4cac0437f3e3d946014c91a1dcffa8fac6072aed56ac010b923` | Five-state reduction fixing lumen Na/K/Cl; same transport core, QSS voltages, volume equation, and acid/base eliminations; derivatives multiplied by 1,000. Formal `g4` is unused. | `par`; intended second call in `Salivary_ex.m`, but called with one argument missing. | `I` (reduction) / `N` | High role; high evidence of incomplete/signature-drifted archive. |
| 13 | `Ae4_Dynamics_Project/Salivary_ex.m` | MATLAB script; 4,007 | SHA-256 `5f3ecaa2473a3bb409a1ff3dac980dcf356897f42923c50d092ae847f1c6f760` | Execution/plot driver. Runs 7-state model with AE2 set to zero, computes pH/CO2/flows, plots states and fluxes, then attempts the broken five-state comparison. | `Par.mat`, `Salivary.m`, `Salivary2.m`. | `I` / `N` | High: explicit calls and knockout assignment; no saved figure dependency. |
| 14 | `Figure1-eps-converted-to.pdf` | PDF figure; 7,551 | SHA-256 `9c588ff3a4f6815ecb78c124d41ed4874b286fe9908da25ae2b11e112c1dc86b` | Direct PDF rendering of `R1`. | Derived from `Figure1.eps`; included in `Ae4_Basis.pdf`. | `R` / `A` | High: log mapping, render, axes. |
| 15 | `Figure1.eps` | MATLAB EPS; 80,204 | SHA-256 `680408895252ecf617ff370793e35f0177dafb726db5b5ef18f8790f2be9a4a9` | `R1`; determinant response to AE2 activity. Created 2017-06-09 11:45:42. | Included by `Ae4_Basis.tex`; source of converted PDF. | `R` / `A` | High: MATLAB header, rendered content, TeX include. |
| 16 | `Figure10-eps-converted-to.pdf` | PDF figure; 13,467 | SHA-256 `b0bd96b71111c6213885610c0d24838b9eebbb4ae37ae9ae255054b66854f14f` | Direct PDF rendering of `S10`. | Derived from `Figure10.eps`; included in manuscript build. | `R` / `A` | High. |
| 17 | `Figure10.eps` | MATLAB EPS; 313,380 | SHA-256 `361814a57ca55e5c6b15e1009956bfda42862c08679c52740d4c8bf58a7416b8` | `S10`; AE2 `KCl` sensitivity, with exchanger formula printed above axes. Created 2017-06-09 13:51:22. | `Ae4_Basis.tex`; converted PDF. | `R` / `A` | High. |
| 18 | `Figure11-eps-converted-to.pdf` | PDF figure; 13,200 | SHA-256 `ccc0c2424ff6e7f5a59abe61ecc5da711d99349a4af71ac675f3e852f371760f` | Direct PDF rendering of `S11`. | Derived from `Figure11.eps`; included in manuscript build. | `R` / `A` | High. |
| 19 | `Figure11.eps` | MATLAB EPS; 313,247 | SHA-256 `970862f79341820d6ca482e0d9d01f17984039ed39fd8daaf145e853abf597c7` | `S11`; AE2 `KB` sensitivity. Created 2017-06-09 13:55:46. | `Ae4_Basis.tex`; converted PDF. | `R` / `A` | High. |
| 20 | `Figure12-eps-converted-to.pdf` | PDF figure; 13,290 | SHA-256 `bfcf2fb1ac10257b4cade6fe1a50034125e25b535bca1070d19d0806ce86dcd2` | Direct PDF rendering of `S12`. | Derived from `Figure12.eps`; included in manuscript build. | `R` / `A` | High. |
| 21 | `Figure12.eps` | MATLAB EPS; 311,765 | SHA-256 `4bee7dd6dc8fd088fed4888c81714bdaa601c37f00ebd0b15a7c2e2aeab65541` | `S12`; NHE1 `KNa` sensitivity. Created 2017-06-09 14:03:52. | `Ae4_Basis.tex`; converted PDF. | `R` / `A` | High. |
| 22 | `Figure13-eps-converted-to.pdf` | PDF figure; 13,529 | SHA-256 `3f92b63f84d3bbf5cecf7c94b580c8272148e9a25b65f80bb6afefcec36cf1a7` | Direct PDF rendering of `S13`. | Derived from `Figure13.eps`; included in manuscript build. | `R` / `A` | High. |
| 23 | `Figure13.eps` | MATLAB EPS; 312,569 | SHA-256 `417f5253744f14bfe77598bb482136a34f7ab4723654f7f7c8ef388bf3ac15e3` | `S13`; NHE1 `KH` sensitivity. Created 2017-06-09 14:08:39. | `Ae4_Basis.tex`; converted PDF. | `R` / `A` | High. |
| 24 | `Figure14-eps-converted-to.pdf` | PDF figure; 13,955 | SHA-256 `906f7adbeb8de30221d7f68e32ba77e2e8afe60b06cf4b5f95c9e70adbdf2575` | Direct PDF rendering of `S14`. | Derived from `Figure14.eps`; included in manuscript build. | `R` / `A` | High. |
| 25 | `Figure14.eps` | MATLAB EPS; 309,114 | SHA-256 `3e0d277efade28f9798c48fc84e6dfd087fdf88454f8fcbc349e71a7d8af5a0e` | `S14`; NKCC1 `a1` sensitivity. Created 2017-06-09 14:54:02. | `Ae4_Basis.tex`; converted PDF. | `R` / `A` | High. |
| 26 | `Figure15-eps-converted-to.pdf` | PDF figure; 14,287 | SHA-256 `f70a55515cc245019ba992fa92537f8850c3d4db564cf0cf80e99f405bb049d4` | Direct PDF rendering of `S15`. | Derived from `Figure15.eps`; included in manuscript build. | `R` / `A` | High. |
| 27 | `Figure15.eps` | MATLAB EPS; 309,768 | SHA-256 `0b8295d99e88745f0453ef81c50671f0de2de6522822a05a6aed7e4bfded430e` | `S15`; NKCC1 `a2` sensitivity. Created 2017-06-09 14:52:10. | `Ae4_Basis.tex`; converted PDF. | `R` / `A` | High. |
| 28 | `Figure16-eps-converted-to.pdf` | PDF figure; 14,113 | SHA-256 `48b344ade25a487e650d0866e4a0db87a546ede1ed5957d4ca26b9f5f66031d1` | Direct PDF rendering of `S16`. | Derived from `Figure16.eps`; included in manuscript build. | `R` / `A` | High. |
| 29 | `Figure16.eps` | MATLAB EPS; 309,713 | SHA-256 `83f9b7a53af01d89d8ae6c66164b4fdc891aa40985e82b2fc58867a626c2f203` | `S16`; NKCC1 `a3` sensitivity. Created 2017-06-09 14:47:07. | `Ae4_Basis.tex`; converted PDF. | `R` / `A` | High. |
| 30 | `Figure17-eps-converted-to.pdf` | PDF figure; 14,102 | SHA-256 `3949f2971f938278c32f9d06487d212bbc9e1cf04c97ed87b9f9924134b27b04` | Direct PDF rendering of `S17`. | Derived from `Figure17.eps`; included in manuscript build. | `R` / `A` | High. |
| 31 | `Figure17.eps` | MATLAB EPS; 309,366 | SHA-256 `afe170c88c282c19bccab94eb24ab5e3756edec69b20eb98fac2f36c3082e46a` | `S17`; NKCC1 `a4` sensitivity. Created 2017-06-09 14:45:32. | `Ae4_Basis.tex`; converted PDF. | `R` / `A` | High. |
| 32 | `Figure18-eps-converted-to.pdf` | PDF figure; 13,447 | SHA-256 `5bdb496032f1377e3d1992097a503b1dd2e653e95c3a9efb72f360ead4f22b63` | Direct PDF rendering of `S18`. | Derived from `Figure18.eps`; included in manuscript build. | `R` / `A` | High. |
| 33 | `Figure18.eps` | MATLAB EPS; 307,923 | SHA-256 `3c52a2353baa43a0f88c08a58c9b81dc10ef3988e2038b99639d5cc33c5486ce` | `S18`; NaK turnover `r` sensitivity. Created 2017-06-09 14:43:46. | `Ae4_Basis.tex`; converted PDF. | `R` / `A` | High. |
| 34 | `Figure19-eps-converted-to.pdf` | PDF figure; 12,722 | SHA-256 `de83fc558394eed7e0406f5c69c7d1b36730379d5ec94408d7b4a29f23858859` | Direct PDF rendering of `S19`. | Derived from `Figure19.eps`; included in manuscript build. | `R` / `A` | High. |
| 35 | `Figure19.eps` | MATLAB EPS; 245,747 | SHA-256 `4994dd488a43d5ef33f3a2b3e6ec889d78b485150f2dc1e419af28b087d0cdd5` | `S19`; NaK `alpha1` sensitivity. Created 2017-06-09 14:41:48. | `Ae4_Basis.tex`; converted PDF. | `R` / `A` | High. |
| 36 | `Figure2-eps-converted-to.pdf` | PDF figure; 12,807 | SHA-256 `73196bb9c90f82f788c5b7591a273b33a4df5a12ee20d3d1df25f4de34783999` | Direct PDF rendering of `R2`. | Derived from `Figure2.eps`; included in `Ae4_Basis.pdf`. | `R` / `A` | High. |
| 37 | `Figure2.eps` | MATLAB EPS; 92,595 | SHA-256 `4d911b5ff2d930fd8296c169033fb3931bb36d78ac27efbf4bf874fc2f0b74e5` | `R2`; determinant response to AE4 activity. Created 2017-06-09 11:36:51. | `Ae4_Basis.tex`; converted PDF. | `R` / `A` | High. |
| 38 | `Figure20-eps-converted-to.pdf` | PDF figure; 13,792 | SHA-256 `ba084a6f9801f83e8378b0030e27ebaccc846705f585090bb35361948185b4b9` | Direct PDF rendering of `S20`. | Derived from `Figure20.eps`; included in manuscript build. | `R` / `A` | High. |
| 39 | `Figure20.eps` | MATLAB EPS; 306,067 | SHA-256 `4b9bf37c63a5ba505aac9e25e9c1c6155ecdf2a1194aa8cb11580c1d4afaa423` | `S20`; CaCC gate constant sensitivity. Created 2017-06-09 15:06:15. | `Ae4_Basis.tex`; converted PDF. | `R` / `A` | High. |
| 40 | `Figure3-eps-converted-to.pdf` | PDF figure; 9,793 | SHA-256 `6be8438ec3e423e0ee762c2355c2b3110f2ed66c7d7e7e9629265041909dab08` | Direct PDF rendering of `R3`. | Derived from `Figure3.eps`; included in `Ae4_Basis.pdf`. | `R` / `A` | High. |
| 41 | `Figure3.eps` | MATLAB EPS; 23,903 | SHA-256 `d2e057e7ee5e53bb2db19656d80efaf3af09d80f6341494c24b7ad04eea1e9b6` | `R3`; normalised knockout-flow comparison. Created 2017-06-09 11:32:37. | `Ae4_Basis.tex`; converted PDF; later replotted as `Figure3.pdf`. | `R` / `A` | High. |
| 42 | `Figure3.pdf` | PDF figure; 125,976 | SHA-256 `d3420e8f980bf7f630e0cbe64a708e6b1a7e9df6131e3852b4e5eb86336b91ca` | `M3`; later styled flow/activity panel. Created 2018-01-08 12:54:13 PST. | Included by `Marty.tex`; semantically derived from `R3`. | `R` / `A` | High. |
| 43 | `Figure4-eps-converted-to.pdf` | PDF figure; 14,322 | SHA-256 `4dd22d6e808e809520af8d577210579ab0afa4a33ab12679c212d6007fde34e8` | Direct PDF rendering of `R4`. | Derived from `Figure4.eps`; included in `Ae4_Basis.pdf`. | `R` / `A` | High. |
| 44 | `Figure4.eps` | MATLAB EPS; 37,069 | SHA-256 `58d74947a3ef76dc7a4bfc8a0155e58fcc270601a2fe39bf582817f87a9a483e` | `R4`; normalised states along AE4 activity. Created 2017-06-09 11:34:24. | `Ae4_Basis.tex`; converted PDF; later replotted as `Figure6.pdf`. | `R` / `A` | High. |
| 45 | `Figure4.pdf` | PDF figure; 109,536 | SHA-256 `fa01cfa15992ba7a10913502fe2b823223a52992970decfcd3fd90ad409fbb8c` | `M4`; later styled AE2 determinant panel. Created 2018-01-08 13:46:06 PST. | Included by `Marty.tex`; semantically derived from `R1`. | `R` / `A` | High. |
| 46 | `Figure5-eps-converted-to.pdf` | PDF figure; 7,988 | SHA-256 `4012ef2d8234a8c9317761b056c670fe09df3754a041e8dd6bb18dc571cf24c5` | Direct PDF rendering of `R5`. | Derived from `Figure5.eps`; included in `Ae4_Basis.pdf`. | `R` / `A` | High. |
| 47 | `Figure5.eps` | MATLAB EPS; 21,075 | SHA-256 `8376f9d97c32da261c369af5ec536d0dd1bef70ca9e280723a14f09cc2f4e06a` | `R5`; normalised states along AE2 activity. Created 2017-06-09 11:35:10. | `Ae4_Basis.tex`; converted PDF; later replotted as `Figure7.pdf`. | `R` / `A` | High. |
| 48 | `Figure5.pdf` | PDF figure; 140,010 | SHA-256 `e0fa7129e547a02752f6f99b69c1a6493a8c170f26f1b193283febcec16f0e74` | `M5`; later styled AE4 determinant panel. Created 2018-01-08 14:50:30 PST. | Included by `Marty.tex`; semantically derived from `R2`. | `R` / `A` | High. |
| 49 | `Figure6-eps-converted-to.pdf` | PDF figure; 6,466 | SHA-256 `45308498c3ee13f6bfc8722e645826eb76063a2267ba1834b472fa2f15f8fa47` | Direct PDF rendering of `R6`. | Derived from `Figure6.eps`; included in `Ae4_Basis.pdf`. | `R` / `A` | High. |
| 50 | `Figure6.eps` | MATLAB EPS; 15,348 | SHA-256 `c22cde381803e7c604bcd8e4d43240bc971f4532a4283c1003499797e19c80bc` | `R6`; Hill-form robustness comparison. Created 2017-06-09 11:30:47. | `Ae4_Basis.tex`; converted PDF. | `R` / `A` | High. |
| 51 | `Figure6.pdf` | PDF figure; 141,131 | SHA-256 `4ac675e5b525a0c1776aa8f0c193f50ff9eced5d21da0495768aa8eebf25b45b` | `M6`; later styled states-versus-AE4 panel. Created 2018-01-08 15:58:23 PST. | Included by `Marty.tex`; semantically derived from `R4`. | `R` / `A` | High. |
| 52 | `Figure7-eps-converted-to.pdf` | PDF figure; 13,694 | SHA-256 `d1c8d4c50c329f06b53dbeec1523304480da3e19af590ffc1ec3d01cf9c5f379` | Direct PDF rendering of `S7`. | Derived from `Figure7.eps`; included in manuscript build. | `R` / `A` | High. |
| 53 | `Figure7.eps` | MATLAB EPS; 305,943 | SHA-256 `c719f3cdf1545bce359fe4a7d4551dbc1c0d2e7ea319acb053bccc9e48bc9459` | `S7`; CaKC gate-constant sensitivity. Created 2017-06-09 15:08:02. | `Ae4_Basis.tex`; converted PDF. | `R` / `A` | High. |
| 54 | `Figure7.pdf` | PDF figure; 118,896 | SHA-256 `e78be5dc79f5dc638f67d19f40247aa3ebdbd7f99171fc6fcc882912b7b26b90` | `M7`; later styled states-versus-AE2 panel. Created 2018-01-08 16:14:28 PST. | Included by `Marty.tex`; semantically derived from `R5`. | `R` / `A` | High. |
| 55 | `Figure8-eps-converted-to.pdf` | PDF figure; 13,383 | SHA-256 `0178b35f1f3c2903af9121de679f1fe70990d702b84a0396718c5221cbcbb058` | Direct PDF rendering of `S8`. | Derived from `Figure8.eps`; included in manuscript build. | `R` / `A` | High. |
| 56 | `Figure8.eps` | MATLAB EPS; 307,409 | SHA-256 `a7e4a9893481be58ebc96c65ecbec7312599a848760ce01309741eb88becbcc4` | `S8`; AE4 forward coefficient `k1` sensitivity. Created 2017-06-09 13:37:47. | `Ae4_Basis.tex`; converted PDF. | `R` / `A` | High. |
| 57 | `Figure8.pdf` | PDF figure; 120,273 | SHA-256 `e15c95152c67cc6af2345b4867e6f0aaa362fa6f5ceb77f3c0efc8dbb499b647` | `M8`; `dQ/dG4` sensitivity to `lambda11`; distinct from `Figure8.eps`. PDF metadata is June 2017. | Included by `Marty.tex`. | `R` / `A` | High: visual/text extraction and include order. |
| 58 | `Figure9-eps-converted-to.pdf` | PDF figure; 13,644 | SHA-256 `5c3b4931c7405911c4e75c11a9f8eab73bcca8348ccddecead79427b4979b72f` | Direct PDF rendering of `S9`. | Derived from `Figure9.eps`; included in manuscript build. | `R` / `A` | High. |
| 59 | `Figure9.eps` | MATLAB EPS; 307,108 | SHA-256 `8a3fe940b87e9fa41f9f7a457d1608f5af246c57443d101f3ef6a9a8b9d050ec` | `S9`; AE4 reverse coefficient `k2` sensitivity. Created 2017-06-09 13:40:26. | `Ae4_Basis.tex`; converted PDF. | `R` / `A` | High. |
| 60 | `Marty.tex` | LaTeX fragment; 24,580 | SHA-256 `2d9929f8c92f2f4b5fe5d505f7f79ce44c5021b31566ba3c95abccd3e930a2da` | Alternate/later mathematical-analysis section. Uses the same three-state reduction and lambda table, explicitly frames it as simplifying the published model, and includes standalone `M3`-`M8` panels. | Requires a missing parent TeX preamble plus `Figure3.pdf`-`Figure8.pdf`; cites `mybib` keys. | `R` / `D` | High scientific relationship; medium exact document chronology because fragment has no build record. |
| 61 | `mchem.sty` | LaTeX style; 130,083 | SHA-256 `4eee0a8e2d9cbd638f2383c53985b1c77f08c87364f3f12c530619f47ae28bf7` | Stock `mhchem` v4.04 chemistry-typesetting package content saved under the nonmatching name `mchem.sty`; no model content. The build log used system `mhchem.sty`, not this file. | None in archived build as named. | `B` / `N` | High: package header and log dependency path. |
| 62 | `mybib.bib` | BibTeX; 26,762 | SHA-256 `04d7a19bd65ed02c25989a54e592924e96eea703eaed98e8b7682ca9a5e27a9d` | 73-record scientific bibliography. Includes Palk 2010/2012, experimental sources, transporter sources, and the AE4 article under key `siguenza2017Ae4`; provenance, not executable data. | Used by `Ae4_Basis.tex`; keys also cited by `Marty.tex`. | `R/B` / `D` | High: all entries and cited-key coverage inspected. |
| 63 | `spbasic.bst` | BibTeX style; 33,250 | SHA-256 `29b1212310e3c58b705725e97a8b56d7b904f8583106424b88d0c9a9ca0469db` | Stock Springer author-year bibliography style (2004); no scientific content. | Selected by `Ae4_Basis.tex`. | `B` / `A` | High: header and `\bibliographystyle` match. |
| 64 | `spmpsci.bst` | BibTeX style; 30,143 | SHA-256 `013368de6cc985cad2d30f76e6b2fa275822ccfa1e74751bc9a576812ec2542e` | Stock Springer numerical style for mathematics/physical sciences; unused alternative. | No archived source selects it. | `B` / `N` | High: style header and repository-wide reference search. |
| 65 | `spphys.bst` | BibTeX style; 28,600 | SHA-256 `763ca46a4846468618846bcc8a8a431762b714709b635d7306a86788ce2b03bc` | Stock Springer physics numerical style; unused alternative. | No archived source selects it. | `B` / `N` | High. |
| 66 | `svglov3.clo` | LaTeX class option; 3,809 | SHA-256 `a0b9f3522e6ef83f3dd8d5901c74745d515355ccd61f0a148668283aabdb5467` | Springer SVJour3 global journal layout option v3.2 (2009-12-18); no model content. | Auto-loaded by `svjour3.cls`; confirmed by log and SyncTeX. | `B` / `A` | High. |
| 67 | `svjour3.cls` | LaTeX class; 47,679 | SHA-256 `7334bcfda97ba34d06e28dc546a7bd7801711a1ba7732f03481dafbe0d04f892` | Springer journal document class v3.2 (2007-05-08); no model content. | Declared by `Ae4_Basis.tex`, loads `svglov3.clo`. | `B` / `A` | High. |

## Completeness and limits

- Ledger rows: **67**, matching the repository tree exactly.
- Machine check: `python -c "import csv; r=list(csv.DictReader(open('results/11_forensic_reconstruction/archive_hashes.csv'))); assert len(r)==67==len({x['path'] for x in r})"` passes. The numbered Markdown ledger was independently parsed as 67 unique paths numbered 1-67, with no path difference from the CSV.
- Content-bearing MATLAB/TeX/BibTeX/MAT files: inspected directly.
- Figures: all 20 EPS sources, all 20 converted PDFs, and all 6 standalone PDFs were accounted for; all 26 PDF renderings were visually inspected, and EPS/PDF derivation was verified from the build log rather than assumed from names.
- Build artefacts: all 10 were inspected for dependency, chronology, and provenance content.
- No archive file was edited. The final archive-integrity check must compare the authoritative Git blob evidence for the two connector caveats and SHA-256 values for the other 65 files; hashing the two local transport placeholders would answer a different question.

This audit establishes roles and provenance relationships. Numerical identity of the MATLAB variants with the publication-generating implementation requires the separate equation concordance and reproduction analyses; it is not inferred here from shared names or directory placement.
