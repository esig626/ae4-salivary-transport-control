# Figure and manuscript provenance

## Scope and evidence classes

This report treats the unpublished manuscript, figures, and build products in
`archive/legacy-2017/` as immutable historical evidence. It does not treat them
as published scientific evidence. Claims below are tagged conceptually as:

- **Published**: the final article, Vera-Sigüenza et al., *Bulletin of
  Mathematical Biology* 80, 255–282, DOI
  [10.1007/s11538-017-0370-6](https://doi.org/10.1007/s11538-017-0370-6),
  online 5 December 2017 and in the February 2018 issue.
- **Historical**: content or metadata preserved under `archive/legacy-2017/`.
- **New inference**: a comparison or interpretation made in this audit.

In this connector-materialized working tree, the two paths for the published
article appear as zero-byte placeholders:

- `archive/legacy-2017/A Mathematical Model Supports a Key Role for Ae4
  (Slc4a9) in Salivary Gland Secretion.pdf`
- `literature/core/2018_AE4_model.pdf`

This is a checkout/transport artifact, not a property of the authoritative
repository objects: remote-object verification in this audit found the two
paths to be identical 1,108,802-byte PDF blobs. Because the local placeholders
cannot be rendered, the published comparison used the canonical [PubMed
Central record
PMC5792321](https://pmc.ncbi.nlm.nih.gov/articles/PMC5792321/), including all
ten full-resolution figure images and their captions. The historical
`Ae4_Basis.pdf`, all twenty EPS/PDF plot pairs, and the six standalone PDFs
were rendered to temporary PNG files outside `archive/`; PDF text and EPS DSC
metadata were also extracted. No archive file was modified.

The machine-readable file-level mapping is
`results/11_forensic_reconstruction/figure_mapping.csv`.

## Bottom line

1. **`Ae4_Basis.tex` is explicitly a derivative mathematical analysis of the
   model that became the published article, not the article manuscript or its
   executable model.** It fixes calcium, acid/base variables, luminal
   concentrations, volumes, and membrane potentials, reduces the model to
   three intracellular variables, and studies steady states and parameter
   derivatives in the Ae2 and Ae4 activity parameters.
2. **The twenty EPS plots and their direct conversion PDFs are a coherent
   June 2017 analysis set.** Their embedded MATLAB R2014b timestamps are 9 June
   2017, and the archived `Ae4_Basis.pdf` was built in September 2017. They
   therefore predate the final article's 5 December 2017 online publication.
3. **None of the historical mathematical-analysis plots duplicates or closely
   reproduces a published panel.** The published figures are a cell schematic,
   a prescribed calcium input, full-model time courses, knockout time courses,
   transporter-flux time courses, and an Ae4 Markov scheme. The historical
   plots use exchanger activity or another parameter on the horizontal axis
   and show reduced steady states, determinants, or normalized derivatives.
4. **There is one strong numerical phenotype correspondence, but it is not
   executable provenance.** Historical `Figure3.eps` and standalone
   `Figure3.pdf` show normalized flow falling to about 0.75–0.76 when the
   reduced-model Ae4 activity is taken to zero, while the Ae2 branch stays near
   one. Published Figure 8 reports about a 24% Ae4-knockout flow reduction and
   negligible Ae2 effect. The plots are not visual duplicates or time courses,
   and no archived MATLAB source generates the continuation plots. Separate
   execution of the available historical models does not recover that endpoint
   phenotype.
5. **`Marty.tex` is a later, incomplete revision/extract of the same reduced
   analysis.** It shares the equations, all 21 reduced parameters, activity
   ranges, and stimulated reference state with `Ae4_Basis.tex`, but embeds six
   later monochrome/redrawn PDFs. Its filename and the `Marty` directory name
   recorded in SyncTeX are project-path evidence only; they are not a basis for
   inferring sole authorship.
6. **The exact analysis-producing source is missing.** `Ae4_Basis.tex` names
   MATCONT and `Marty.tex` names Mathematica, but there is no archived MATCONT
   continuation system, saved continuation data, or Mathematica notebook.
   Across the six archived `.m` sources, every figure call is time-domain;
   there is no `print`, `saveas`, `exportgraphics`, or MATCONT call.

## Chronology

Dates below are embedded creator/build metadata, not filesystem modification
times. They establish ordering, but not authorship. Where PDF and log times use
different time zones, the metadata is reported rather than forced into a
single presumed local time.

| Date | Historical/published item | Evidence and implication |
|---|---|---|
| 8 June 2017, 20:17:43 PDT | Embedded creation metadata in standalone `Figure8.pdf` | MATLAB R2014b/Ghostscript metadata and a title referring to `Figure8.eps`. The PDF was later modified on 26 March 2018, so the creation field most safely dates its underlying graphic lineage, not necessarily the final edited file. |
| 9 June 2017, 11:30–15:08 | `Figure1.eps`–`Figure20.eps` | Every EPS identifies MATLAB R2014b on Mac OS X and `/Users/esig526/Desktop/New_Paper/FigureN.eps`; creation times span one work session. The converted PDFs are the corresponding Ghostscript 9.10/`epstopdf` products. |
| 18 September 2017, 22:40:26 PDT | `Ae4_Basis.pdf` | PDF metadata: pdfTeX 1.40.15, 19 pages. |
| 19 September 2017, 17:40 | `Ae4_Basis.log` | TeX Live 2014 build log; successful 19-page build using the twenty converted EPS PDFs and an `Ae4_Basis.bbl` that is no longer archived. |
| 5 December 2017 | Published article | Final edited article published online; DOI metadata and PMC record. The June plots and September analysis PDF precede this date. |
| 8 January 2018 | Standalone `Figure3.pdf`–`Figure7.pdf` | Adobe Illustrator CC 2014 creation metadata. These are later monochrome/redrawn versions of five June plots. |
| February 2018 | Published article issue date | Volume 80, issue 2. |
| 26 March 2018 | Standalone `Figure3.pdf`–`Figure8.pdf` | All six have this modification date. These edited assets postdate online publication and are used by the later `Marty.tex` fragment. |

This ordering supports a specific but limited conclusion: the colored
continuation/sensitivity figure family and the complete `Ae4_Basis` analysis
precede the published article, while the polished standalone family postdates
it. Ordering alone does not show that any figure was generated by the code used
for the published full-model simulations.

## Manuscript roles and lineage

### `Ae4_Basis.tex` and `Ae4_Basis.pdf`

`Ae4_Basis.tex` is a complete 783-line Springer `svjour3` document titled
"Mathematical analysis of a secretion model." It names Elías Vera-Sigüenza,
Martin Golubitsky, and James Sneyd as authors. Its introduction states that it
analyses the earlier Sigüenza et al. secretion model and uses numerical
continuation to explain its Ae2/Ae4 knockout predictions.

The mathematical role is unambiguous from the source:

- Lines 158–159 fix intracellular calcium, bicarbonate, carbon dioxide,
  hydrogen, all three luminal ions, intracellular and luminal volumes, and
  both membrane potentials.
- Lines 162–199 retain only intracellular Na, K, and Cl, renamed \(x,y,z\),
  and write a three-state dimensionless system.
- Lines 203–226 define the reduced flow
  \(Q=\lambda_{20}-\lambda_{21}(x+y)\), 21 condensed parameters, reference
  state \(x_0=12.29,\ y_0=89.19,\ z_0=37.13\), and activities
  \(G_2=0.347,\ G_4=6.74\times10^{-5}\).
- Line 282 states that the \(G_2\in[0,0.347]\) and
  \(G_4\in[0,6.74\times10^{-5}]\) branches were solved with MATCONT.
- Figures 1–4 then show continuation determinants, normalized reduced states,
  normalized flow, and local parameter sweeps. They are not outputs of a
  transient full-cell simulation.

The analysis is structurally descended from the archived full/reduced MATLAB
equations, but its numerical system is not stored in those `.m` files. A
repository search finds none of the distinctive reduced literals `0.347`,
`6.74e-5`, `4589.648`, `12.29`, `89.19`, or `37.13` in archived MATLAB. Some
condensed values retain visible ancestry—for example λ3 = 1.03 approximately
matches `a3 = 1.0306`, λ4 = 1,385,200 matches `a4 = 1.3852e6`, λ7 = 0.641
matches `alpha1`, λ10 = 15 matches `KNa`, and λ19 = 5.6 matches `KCl`—but
this is a derived parameterization rather than a file-to-file executable
match.

The historical text must not be treated as a corrected specification. It
contains internal or paper-facing inconsistencies, including:

- `Ae4_Basis.tex` gives \(d\omega_i/dt=q_b-q_a\), whereas published Eq. 9
  gives \(q_a-q_b\).
- The reduced Ae4 term retains an intracellular Na contribution but not the K
  contribution present in published Eq. 31.
- Its Ae2 appendix lists \(K_B=10^4\) mM, as the published appendix does, while
  both historical code and `Figure11.eps` use a baseline near \(10^{-4}\) mM,
  an eight-order discrepancy.
- Figure and table labels are repeatedly reused, three main captions are
  empty, and the source contains draft red text and several algebraic/notation
  inconsistencies.
- The plotted \(G_2\to0\) state branch visibly places normalized \(z\) near
  0.967, although the later `Marty.tex` caption says all three variables vary
  by less than 1%. Similar prose/caption percentages for the \(G_4\) state
  branch are not exactly the plotted values.

These features are consistent with an unpublished analysis draft and argue
against using it as an authoritative replacement for the published model.

### `Marty.tex`

`Marty.tex` is a 284-line fragment: it has no document class, preamble,
`\begin{document}`, or `\end{document}`, and cannot build standalone. It
repeats the same three-state reduced equations, exact λ1–λ21 table,
activity values, reference state, flow expression, determinant analysis, and
knockout interpretation as `Ae4_Basis.tex`. It explicitly calls the fixed-cell-
volume assumption "incorrect" but useful for reducing the analysis to three
differential equations. It also says the Jacobian was calculated with
Wolfram Mathematica, but no notebook is present.

The fragment's six-panel figure uses the later standalone PDFs and supplies a
complete caption. Structural identity and the later asset metadata make it a
later revision or extraction of the same collaborative analysis with high
confidence. The filename `Marty.tex` does not establish who wrote or owned all
of its content.

### `mybib.bib`

`mybib.bib` is a 73-entry project bibliography, not a minimal frozen build
database. The `siguenza2017Ae4` entry is a prepublication-style record with
title, authors, journal, year, and publisher but no DOI, volume, or pages. The
file also contains `almassy2018new`, showing that the bibliography may have
been changed after the September 2017 build.

The archived `.aux` contains resolved `\bibcite` records for
`martinez1966micropuncture` and `holappa2001identification`, but neither entry
exists in the current `mybib.bib`. Because `Ae4_Basis.bbl` is also absent, the
exact bibliography state that produced the archived PDF is not recoverable
from the archive as stored.

## TeX build records and dependencies

### Historical build

`Ae4_Basis.log` proves a successful build:

- pdfTeX 1.40.15, TeX Live 2014;
- local `svjour3.cls` and `svglov3.clo`;
- system `mhchem.sty` v3.17, not the archived `mchem.sty`;
- `spbasic.bst` and a generated `Ae4_Basis.bbl`;
- all twenty `FigureN-eps-converted-to.pdf` files embedded on pages 11–14;
- output `Ae4_Basis.pdf`, 19 pages, 431,602 bytes.

The archive supplies `svjour3.cls`, `svglov3.clo`, `spbasic.bst`, the EPS
sources, and the already converted PDFs. `spmpsci.bst` and `spphys.bst` are
not referenced. `mchem.sty` is orphaned/misnamed relative to
`\usepackage{mhchem}` and was not used by the historical build.

`Ae4_Basis.synctex.gz` adds dependency and path evidence. It records the main
source and local class files under
`/Users/esig526/Documents/Work/Marty/`, plus TeX Live 2014 system packages.
The directory name is a build-location fact, not author identity evidence.

The historical log reports a successful artifact despite warnings: `cite`
used with `natbib`, unsupported caption class handling, no explicit mhchem
version, many multiply-defined labels, font substitutions, and over/underfull
boxes. The `.aux` confirms that Figure 1 has a complete caption on page 11;
Figures 2–4 have empty main captions on pages 12–14.

### Current reproducibility check

An unchanged copy of the source and dependencies was tested only in a
temporary directory outside `archive/`:

1. With the current TeX installation, the build stops because `mhchem.sty` is
   unavailable.
2. The archived `mchem.sty` internally provides mhchem v4.04 but has the wrong
   filename for the source and depends on obsolete `l3regex` interfaces that
   are incompatible with the current TeX Live 2023 packages.
3. The manuscript contains no `\ce` call or other mhchem command. Supplying a
   temporary no-op `mhchem.sty` compatibility stub allows `latexmk` to produce
   a 19-page PDF. This demonstrates that the document and figures are
   otherwise buildable; it is not an exact historical rebuild.
4. Current BibTeX emits missing-entry warnings for
   `martinez1966micropuncture` and `holappa2001identification`. The archived
   `.aux` and PDF resolved both, confirming that an exact historical
   bibliography dependency is missing.

Thus the archived PDF itself is the authoritative rendering of the September
2017 document. A layout-equivalent current build is possible with a package
shim, but the exact original dependency closure is incomplete.

## Exact TeX-to-figure mapping

### `Ae4_Basis.tex`

| Manuscript figure | Panel | Included file | Visual content |
|---|---:|---|---|
| 1 | A | `Figure1.eps` | Det(\(J_{x_2}\)) and Det(\(J_{y_2}\)) versus \(G_2\) |
| 1 | B | `Figure2.eps` | Det(\(J_{x_4}\)) and Det(\(J_{y_4}\)) versus \(G_4\) |
| 1 | C | `Figure5.eps` | normalized \(x,y,z\) versus \(G_2/G_{2,0}\) |
| 1 | D | `Figure4.eps` | normalized \(x,y,z\) versus \(G_4/G_{4,0}\) |
| 1 | E | `Figure3.eps` | normalized \(Q\) versus normalized exchanger activity |
| 1 | F | `Figure6.eps` | analogous normalized \(Q\) curves for a Hill-type alternative model |
| 2 | A | `Figure7.eps` | sensitivity to Ca-activated K-channel \(K_{CaKC}\) |
| 2 | B | `Figure20.eps` | sensitivity to Ca-activated Cl-channel \(K_{CaCC}\) |
| 2 | C | `Figure14.eps` | sensitivity to Nkcc1 \(a_1\) |
| 2 | D | `Figure15.eps` | sensitivity to Nkcc1 \(a_2\) |
| 2 | E | `Figure16.eps` | sensitivity to Nkcc1 \(a_3\) |
| 2 | F | `Figure17.eps` | sensitivity to Nkcc1 \(a_4\) |
| 3 | A | `Figure18.eps` | sensitivity to NaK rate \(r\) |
| 3 | B | `Figure19.eps` | sensitivity to NaK half-saturation α1 |
| 3 | C | `Figure8.eps` | sensitivity to Ae4 forward coefficient \(k_1\) |
| 3 | D | `Figure9.eps` | sensitivity to Ae4 reverse coefficient \(k_2\) |
| 3 | E | `Figure10.eps` | sensitivity to Ae2 \(K_{Cl}\) |
| 3 | F | `Figure11.eps` | sensitivity to Ae2 \(K_B\) |
| 4 | A | `Figure12.eps` | sensitivity to Nhe1 \(K_{Na}\) |
| 4 | B | `Figure13.eps` | sensitivity to Nhe1 \(K_H\) |

The historical pdfTeX run actually consumes the direct conversion PDFs, but
the TeX source names the EPS originals and `epstopdf` resolves each to its
`FigureN-eps-converted-to.pdf` sibling. Visual rendering and text extraction
show that every converted PDF is the same logical plot as its EPS source.

Only Figure 1 has a substantive manuscript caption. Figures 2–4 have blank
main captions, so the parameter names in the plot equations and axes are the
only archived descriptions of those panels. The source reuses labels such as
`fig:det_G2`, `fig:Q`, and `fig:results` across multiple figures; the `.aux`
and log record the resulting multiply-defined-label warnings.

### `Marty.tex`

| Panel | Included file | Relation to June 2017 plot |
|---:|---|---|
| A | `Figure3.pdf` | monochrome edited/redrawn counterpart of `Figure3.eps`; \(G_4\to0\) branch ends near 0.75 and \(G_2\to0\) branch near 1 |
| B | `Figure8.pdf` | three-point 50/100/150% sweep of reduced λ11; same conceptual Ae4-forward-rate sensitivity as `Figure8.eps`, but reparameterized and not an exact conversion |
| C | `Figure7.pdf` | monochrome edited/redrawn counterpart of `Figure5.eps` (\(x,y,z\) versus \(G_2\)) |
| D | `Figure6.pdf` | monochrome edited/redrawn counterpart of `Figure4.eps` (\(x,y,z\) versus \(G_4\)) |
| E | `Figure4.pdf` | monochrome edited/redrawn counterpart of `Figure1.eps` (determinants versus \(G_2\)) |
| F | `Figure5.pdf` | monochrome edited/redrawn counterpart of `Figure2.eps` (determinants versus \(G_4\)) |

The standalone files are not merely the June conversion PDFs renamed. Their
fonts, line styling, legends, crop/page sizes, and embedded producer metadata
differ. `Figure3.pdf`–`Figure7.pdf` were edited with Illustrator in January
and March 2018. Standalone `Figure8.pdf` retains MATLAB/Ghostscript creation
metadata, was modified in March 2018, changes the horizontal variable from
\(k_1\) to λ11, and displays only three \(dQ/dG_4\) samples.

## Visual and numeric inventory of all plot families

### Continuation, states, flow, and determinants (`Figure1`–`Figure6`)

| File | Horizontal range | Visible result |
|---|---|---|
| `Figure1.eps` | \(G_2=0\) to about 0.30 | Det(\(J_{x_2}\)) is of order \(10^{28}\) and rises slightly; Det(\(J_{y_2}\)) lies visually near zero on that shared scale. |
| `Figure2.eps` | \(G_4=0\) to about \(6.7\times10^{-5}\) | Both determinants decline; the two displayed branches are approximately \(2.2\to1.35\times10^5\) and \(1.18\to0.98\times10^5\). |
| `Figure3.eps` | \(G_i/G_{i,0}=0\) to 1 | The branch obtained while holding \(G_2=0.34\) fixed (therefore varying \(G_4\)) rises from about 0.76 to 1. The branch holding \(G_4=6.74\times10^{-5}\) fixed (therefore varying \(G_2\)) stays about 0.98–1. The legend names the fixed companion parameter, so reading the legend as the varied parameter would reverse the biological interpretation. |
| `Figure4.eps` | \(G_4/G_{4,0}=0\) to 1 | At \(G_4=0\), normalized \(x\approx1.15\), \(y\approx1.03\), and \(z\approx0.93\); all meet at 1 at baseline activity. |
| `Figure5.eps` | \(G_2/G_{2,0}=0\) to 1 | Small changes: at \(G_2=0\), \(x\approx1.006\), \(y\approx1.003\), and \(z\approx0.967\). |
| `Figure6.eps` | \(G_i/G_{i,0}=0\) to 1 | Hill-model alternative with essentially the same qualitative split: a branch near 0.76 at zero activity and a nearly flat branch near one. |

The standalone equivalents preserve those curve shapes but simplify the
styling. In particular, standalone `Figure3.pdf` removes the ambiguous
fixed-parameter legend and labels the near-flat dashed branch \(G_2\) and the
flow-losing solid branch \(G_4\).

### Sensitivity family (`Figure7`–`Figure20`)

All fourteen plots display normalized responses labelled as \(dQ/dG_4\) and
\(dQ/dG_2\) against one model parameter. Across the family the \(G_4\)
derivative is normalized near one and the \(G_2\) derivative lies near zero.
Most use 50%, 100%, and 150% markers; the two channel-calcium plots show points
at horizontal values 90, 100, and 103 and annotate them as 90%, 100%, and 103%
of nominal. The plot/code unit convention is not sufficiently reliable to add
a physical unit. The visual sweep inventory is:

| File | Parameter | Approximate 50–150% range or displayed range | Visible behavior |
|---|---|---|---|
| `Figure7.eps` | \(K_{CaKC}\) | displayed 90–103; markers labelled 90/100/103% | modest positive change in normalized \(dQ/dG_4\); \(dQ/dG_2\approx0\) |
| `Figure20.eps` | \(K_{CaCC}\) | displayed 90–103; markers labelled 90/100/103% | nearly flat to slightly positive; \(dQ/dG_2\approx0\) |
| `Figure14.eps` | Nkcc1 \(a_1\) | about 78.8–236.3 s⁻¹ | small negative slope |
| `Figure15.eps` | Nkcc1 \(a_2\) | about \(1.00–3.01\times10^7\) mM⁻⁴s⁻¹ | nearly flat |
| `Figure16.eps` | Nkcc1 \(a_3\) | about 0.515–1.546 s⁻¹ | small positive slope |
| `Figure17.eps` | Nkcc1 \(a_4\) | about \(0.693–2.078\times10^6\) mM⁻⁴s⁻¹ | small positive slope |
| `Figure18.eps` | NaK \(r\) | about \(0.653–1.958\times10^6\) mM⁻³s⁻¹ | modest negative slope |
| `Figure19.eps` | NaK α1 | about 0.321–0.962 mM⁻¹ | nearly flat |
| `Figure8.eps` | Ae4 \(k_1\) | about 0.0096–0.0288 | dominant response, approximately 0.5→1.4 in normalized \(dQ/dG_4\) |
| `Figure9.eps` | Ae4 \(k_2\) | about \(0.67–2.01\times10^{-5}\) | near one with a small negative slope |
| `Figure10.eps` | Ae2 \(K_{Cl}\) | 2.8–8.4 mM | small variation in \(dQ/dG_4\); \(dQ/dG_2\approx0\) |
| `Figure11.eps` | Ae2 \(K_B\) | \(0.5–1.5\times10^{-4}\) mM | small variation; uses the code-scale \(10^{-4}\), not the manuscript/published \(10^4\) |
| `Figure12.eps` | Nhe1 \(K_{Na}\) | 7.5–22.5 mM | nearly flat |
| `Figure13.eps` | Nhe1 \(K_H\) | \(2.25–6.75\times10^{-4}\) mM | nearly flat |

The family supports the unpublished manuscript's intended local conclusion—
the reduced output derivative is much more responsive to the Ae4 forward
coefficient than to the other tested parameters over these chosen slices.
It is not a global robustness analysis, does not show a sampled physiological
uncertainty distribution, and cannot establish identifiability.

## Parameter-sweep provenance

The sensitivity axes link the plots strongly to the parameter conventions of
the published appendices and `Original/Parameters.m`, even though the exact
plotting program is absent:

| Quantity | Figure baseline | Published appendix | `Original/Parameters.m` | Outer `Par.m` | Inference |
|---|---:|---:|---:|---:|---|
| Nkcc1 \(a_1\) | 157.55 | 157.5 | 157.55 | 157.55 | same raw scale |
| Nkcc1 \(a_2\) | \(2.0096\times10^7\) | \(2.0096\times10^7\) | \(2.0096\times10^7\) | \(2.0096\times10^{-5}\) | figure follows published/Original scale, not outer rescaling |
| Nkcc1 \(a_3\) | 1.0306 | 1.0306 | 1.0306 | 1.0306 | exact |
| Nkcc1 \(a_4\) | \(1.3852\times10^6\) | \(1.3852\times10^6\) | \(1.3852\times10^6\) | \(1.3852\times10^{-6}\) | figure follows published/Original scale, not outer rescaling |
| NaK \(r\) | \(1.305\times10^6\) | \(1.305\times10^6\) | \(1.305\times10^6\) | same literal family | same raw scale |
| NaK α1 | 0.641 | 0.641 | 0.641 | 0.641 | exact |
| Ae4 \(k_1\) | \(1.92\times10^{-2}\) | \(1.92\times10^{-2}\) | \(1.92\times10^{-2}\) | \(1.92\times10^{-2}\) | exact |
| Ae4 \(k_2\) | about \(1.3\times10^{-5}\) | \(1.3\times10^{-5}\) | \(1.34184\times10^{-5}\) | \(1.59306\times10^{-3}\) | figure follows published/Original scale |
| Ae2 \(K_B\) | \(10^{-4}\) | \(10^4\) | \(10^{-4}\) | \(10^{-4}\) | figure exposes the code/manuscript sign-of-exponent contradiction |

This is strong evidence that the June plots analyse a reduction parameterized
from the same model family and many of the same raw constants as the final
article. It is not evidence that the archived time-domain MATLAB files emitted
the plots. In particular, the reduced activity values \(G_2=0.347\) and
\(G_4=6.74\times10^{-5}\) do not equal outer `Par.m` values
`g2 = 0.3991` and `g4 = 18.9503`, and the complete λ table appears only in
the TeX fragments.

The source states that MATCONT generated the activity branches. The parameter
sensitivity plots likely evaluate the Cramer-rule derivative expressions
along those branches or their reference steady state, but no saved numerical
table or plotting script is present. That exact computational path therefore
remains **unresolved**.

## Comparison with every published 2018 figure

| Published figure | Published content | Historical relation |
|---:|---|---|
| 1 | Salivary acinar-cell schematic with transporters, compartments, currents, and water paths | No historical analysis-panel match. `Ae4_Basis.tex` describes a related network in equations, but none of `Figure1`–`Figure20` is a schematic. |
| 2 | Prescribed Ca²⁺ input over 0–20 min; 58 nM rest, stimulation between minutes 6 and 12 | No match. Calcium is fixed/eliminated in the three-state analysis, and no historical figure has time or calcium on an axis. |
| 3 | Six WT intracellular time courses | No direct match. Historical \(x,y,z\) plots are normalized reduced steady-state branches versus exchanger activity. |
| 4 | WT luminal Na/K/Cl, cell volume, and apical/basolateral voltage time courses | No match; these variables are fixed or eliminated in the reduced analysis. |
| 5 | WT flow-rate time course with about 13-fold stimulation | No match. Historical flow plots are normalized steady summaries against exchanger activity and do not show the rest/stimulation protocol. |
| 6 | Six WT/Ae2-KO/Ae4-KO intracellular time-course comparisons | Directional qualitative correspondence only: historical `Figure4.eps` shows increased normalized Na and K and decreased Cl as \(G_4\to0\), while `Figure5.eps` shows much smaller \(G_2\) effects. It is not a reproduced panel or trajectory. |
| 7 | Knockout luminal concentrations and membrane-potential time courses | No match; those variables are absent from the reduced plots. |
| 8 | WT/Ae2-KO/Ae4-KO flow time courses; about 24% Ae4 reduction | Strong endpoint-summary correspondence to historical `Figure3.eps`/`Figure3.pdf` (about 24–25% loss on the \(G_4\to0\) branch and negligible \(G_2\) loss), but different axes, curves, protocol, and graphic. Not a duplicate. |
| 9 | Ae4, Ae2, Nkcc1, and Nhe1 flux time courses under WT/knockout conditions | No match. Historical sensitivity panels print some transporter formulas but plot normalized derivatives, not transporter flux trajectories. |
| 10 | Four-state bidirectional Ae4 Markov scheme | No match. The analysis uses a condensed mass-action Ae4 term and has no state-scheme graphic. |

The graphical vocabularies are plainly different. Published Figures 2–9 use
time, generally 0–20 or 0–30 minutes, and label an agonist interval; the
historical analysis plots use \(G_2\), \(G_4\), normalized exchanger activity,
or a parameter value. This excludes direct panel identity independently of
curve-color or phenotype similarity.

## Can the archived code explain the 24–25% plot?

No archived `.m` file contains the continuation system or emits any of the
archived continuation/sensitivity graphics:

- `Original/Saliva_Ae4.m` and `Salivary_ex.m` use only time-domain
  `figure(1)`/`figure(2)` plotting; `Salivary.m` and `Salivary2.m` are
  right-hand-side routines and emit no figures.
- There is no MATCONT call, continuation callback, print/export command, or
  stored continuation result.
- The three-state λ parameterization exists only in `Ae4_Basis.tex` and
  `Marty.tex`.

The independent historical-execution audit further separates numerical
correspondence from provenance. Under the archived protocols, a translated
SciPy-BDF harness for the seven-state `Salivary.m` gives endpoint flow ratios
AE2-KO/WT = 0.999052 and AE4-KO/WT = 0.998556. At the first sampled point after
the calcium step (t = 101 in that harness), the ratios are 0.997382 and
1.101877. The fixed-lumen five-state `Salivary2.m` gives endpoint ratios
1.000000 and 1.059019, with first post-step ratios 1.002023 and 0.903297.
`Original/Saliva_Ae4.m` is a WT time-course driver whose final historical call
does not match the archived function signature, and `Salivary_ex.m` hard-codes
`g2 = 0` despite a "Full Model" comment. None produces the 0.75–0.76 Ae4
endpoint in the historical figure.

Therefore the 24–25% continuation result most likely comes from the separate
three-state MATCONT analysis described in TeX, whose executable source/data are
absent. It cannot be used to prove that any runnable archived MATLAB variant is
the exact published simulation implementation.

## Provenance conclusions and confidence

| Question | Conclusion | Confidence and basis |
|---|---|---|
| Does the unpublished manuscript explicitly reduce the published model family? | Yes. It says so and enumerates fixed variables before deriving the three-state system. | High; direct historical text and equation structure. |
| Did the colored analysis plots precede the published article? | Yes, according to embedded June 2017 figure metadata and the September 2017 TeX build. | High for ordering; metadata does not identify the operator. |
| Do any archive figures duplicate published panels? | No. | High; all 46 historical figure files and all ten published figures were rendered/inspected, with incompatible axes and content. |
| Do archive plots encode the published Ae4/Ae2 phenotype summary? | Yes, the reduced flow and state branches encode a similar endpoint/directional summary. | High for what is visible; only moderate as model-implementation provenance. |
| Were the June plots generated by the archived seven-/five-state MATLAB files? | No evidence supports that, and code/parameter/output evidence argues against it. | High that no generator is archived; the missing MATCONT source prevents identifying the actual generator. |
| Is `Marty.tex` later than `Ae4_Basis.tex`? | Very likely. | High from exact structural reuse plus January/March 2018 edited figure metadata. |
| Does a filename or visual resemblance establish individual authorship or executable identity? | No. | Evidentiary constraint; no identity inference is made. |

The most defensible characterization is: **the June/September historical
materials directly precede and mathematically analyse the same model project,
and they preserve a reduced steady-state explanation of the published
knockout claim; they neither reproduce the published figure panels nor supply
the missing full-model simulation pipeline.**

## Inspection coverage

The following families were inspected, not merely inventoried:

- `Ae4_Basis.tex`, its 19-page rendered `Ae4_Basis.pdf`, `.aux`, `.log`, and
  decompressed `.synctex.gz` dependency records;
- `Marty.tex` and all of its `\includegraphics` targets;
- all 73 records in `mybib.bib` at the dependency/provenance level;
- `svjour3.cls`, `svglov3.clo`, `spbasic.bst`, `spmpsci.bst`, `spphys.bst`,
  and `mchem.sty` for used/orphan build roles;
- all twenty EPS files, their twenty direct conversion PDFs, and all six
  standalone PDFs, by render plus axis/text/metadata extraction;
- all ten published figures and captions from the canonical PMC article;
- source-wide searches for the reduced parameter literals, figure export
  commands, and MATCONT/Mathematica dependencies;
- a temporary current-environment TeX build test and the separate historical
  execution results summarized above.
