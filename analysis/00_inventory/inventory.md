# Phase 00 historical inventory

## Method and handling

This inventory was produced from the Git tree paths, file sizes/extensions, the repository’s `archive/README.md`, and the commit that moved the historical directory intact. No file content under `archive/legacy-2017/` was opened. Classification is therefore by repository context, filename, and format; ambiguous scientific content is not inferred.

The historical tree contains 67 files:

| Classification | Count |
|---|---:|
| Published source | 1 |
| Unpublished manuscript material | 3 |
| Historical code | 7 |
| Figure | 46 |
| Build artefact | 10 |
| **Total** | **67** |

`Ae4_Dynamics_Project/` and `Ae4_Dynamics_Project/Original/` are directory containers for the historical code. They contain no additional unlisted files.

## File-level inventory

| Path relative to `archive/legacy-2017/` | Classification | Basis |
|---|---|---|
| `A Mathematical Model Supports a Key Role for Ae4 (Slc4a9) in Salivary Gland Secretion.pdf` | Published source | Published article copy; canonical working pointer is `literature/core/2018_AE4_model.pdf` |
| `Ae4_Basis.aux` | Build artefact | LaTeX auxiliary file |
| `Ae4_Basis.log` | Build artefact | LaTeX build log |
| `Ae4_Basis.pdf` | Build artefact | Compiled output of historical unpublished source |
| `Ae4_Basis.synctex.gz` | Build artefact | TeX synchronization output |
| `Ae4_Basis.tex` | Unpublished manuscript material | Historical TeX source; content not consulted |
| `Ae4_Dynamics_Project/Original/Parameters.m` | Historical code | MATLAB source; content not consulted |
| `Ae4_Dynamics_Project/Original/Saliva_Ae4.m` | Historical code | MATLAB source; content not consulted |
| `Ae4_Dynamics_Project/Par.m` | Historical code | MATLAB source; content not consulted |
| `Ae4_Dynamics_Project/Par.mat` | Historical code | MATLAB data/workspace file; content not consulted |
| `Ae4_Dynamics_Project/Salivary.m` | Historical code | MATLAB source; content not consulted |
| `Ae4_Dynamics_Project/Salivary2.m` | Historical code | MATLAB source; content not consulted |
| `Ae4_Dynamics_Project/Salivary_ex.m` | Historical code | MATLAB source; content not consulted |
| `Figure1-eps-converted-to.pdf` | Figure | Converted figure output |
| `Figure1.eps` | Figure | Encapsulated PostScript figure |
| `Figure10-eps-converted-to.pdf` | Figure | Converted figure output |
| `Figure10.eps` | Figure | Encapsulated PostScript figure |
| `Figure11-eps-converted-to.pdf` | Figure | Converted figure output |
| `Figure11.eps` | Figure | Encapsulated PostScript figure |
| `Figure12-eps-converted-to.pdf` | Figure | Converted figure output |
| `Figure12.eps` | Figure | Encapsulated PostScript figure |
| `Figure13-eps-converted-to.pdf` | Figure | Converted figure output |
| `Figure13.eps` | Figure | Encapsulated PostScript figure |
| `Figure14-eps-converted-to.pdf` | Figure | Converted figure output |
| `Figure14.eps` | Figure | Encapsulated PostScript figure |
| `Figure15-eps-converted-to.pdf` | Figure | Converted figure output |
| `Figure15.eps` | Figure | Encapsulated PostScript figure |
| `Figure16-eps-converted-to.pdf` | Figure | Converted figure output |
| `Figure16.eps` | Figure | Encapsulated PostScript figure |
| `Figure17-eps-converted-to.pdf` | Figure | Converted figure output |
| `Figure17.eps` | Figure | Encapsulated PostScript figure |
| `Figure18-eps-converted-to.pdf` | Figure | Converted figure output |
| `Figure18.eps` | Figure | Encapsulated PostScript figure |
| `Figure19-eps-converted-to.pdf` | Figure | Converted figure output |
| `Figure19.eps` | Figure | Encapsulated PostScript figure |
| `Figure2-eps-converted-to.pdf` | Figure | Converted figure output |
| `Figure2.eps` | Figure | Encapsulated PostScript figure |
| `Figure20-eps-converted-to.pdf` | Figure | Converted figure output |
| `Figure20.eps` | Figure | Encapsulated PostScript figure |
| `Figure3-eps-converted-to.pdf` | Figure | Converted figure output |
| `Figure3.eps` | Figure | Encapsulated PostScript figure |
| `Figure3.pdf` | Figure | PDF figure output |
| `Figure4-eps-converted-to.pdf` | Figure | Converted figure output |
| `Figure4.eps` | Figure | Encapsulated PostScript figure |
| `Figure4.pdf` | Figure | PDF figure output |
| `Figure5-eps-converted-to.pdf` | Figure | Converted figure output |
| `Figure5.eps` | Figure | Encapsulated PostScript figure |
| `Figure5.pdf` | Figure | PDF figure output |
| `Figure6-eps-converted-to.pdf` | Figure | Converted figure output |
| `Figure6.eps` | Figure | Encapsulated PostScript figure |
| `Figure6.pdf` | Figure | PDF figure output |
| `Figure7-eps-converted-to.pdf` | Figure | Converted figure output |
| `Figure7.eps` | Figure | Encapsulated PostScript figure |
| `Figure7.pdf` | Figure | PDF figure output |
| `Figure8-eps-converted-to.pdf` | Figure | Converted figure output |
| `Figure8.eps` | Figure | Encapsulated PostScript figure |
| `Figure8.pdf` | Figure | PDF figure output |
| `Figure9-eps-converted-to.pdf` | Figure | Converted figure output |
| `Figure9.eps` | Figure | Encapsulated PostScript figure |
| `Marty.tex` | Unpublished manuscript material | Historical TeX source; content not consulted |
| `mchem.sty` | Build artefact | TeX style dependency |
| `mybib.bib` | Unpublished manuscript material | Bibliography database associated with historical manuscript; content not consulted |
| `spbasic.bst` | Build artefact | BibTeX style dependency |
| `spmpsci.bst` | Build artefact | BibTeX style dependency |
| `spphys.bst` | Build artefact | BibTeX style dependency |
| `svglov3.clo` | Build artefact | TeX class option dependency |
| `svjour3.cls` | Build artefact | TeX document class dependency |

## Use decision

- The duplicate published article is scientifically usable, but baseline reconstruction uses the separate canonical pointer `literature/core/2018_AE4_model.pdf`.
- Unpublished manuscript material, historical code, and historical figures are immutable evidence only. They supplied no prose, equations, values, figures, or code to the new model files.
- Build artefacts are retained for provenance and have no role in reconstruction.

The new source map is therefore independent of the unpublished archive. Its remaining blockers arise from the published record and are listed in `open_questions.md`.
