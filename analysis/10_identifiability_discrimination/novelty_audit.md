# Targeted novelty audit

## Search date, question, and scope

Search completed **2026-08-27**. This is a focused claim-level audit, not a
systematic review. It asks whether the strongest results currently available in
this repository are new enough to support the proposed quick paper:

1. generic nonidentifiability of two transporter activities from one scalar
   steady-state output;
2. local-rank and minimal-measurement results for epithelial or pump--leak
   models;
3. structural identifiability of ion-channel or ion-transporter models;
4. recovery of identifiability by perturbations or multiple experimental
   conditions; and
5. stoichiometric-signature or balance-law rank criteria.

The searches used combinations of `pump-leak`, `epithelial transport`,
`steady state`, `inverse`, `identifiability`, `minimal output set`, `flow
measurement`, `ion transport`, `multiple experiments`, `constant input`,
`perturbation design`, `stoichiometric matrix`, `reaction network`, and
`structural sensitivity`. Primary journal pages, DOI records, PubMed/PMC, and
author or institutional copies were preferred. Citation trails from the
closest papers were checked. The search included papers available online by
the search date, including two directly relevant 2026 pump--leak papers.

The audit deliberately distinguishes:

- transporter **activities** `(G2,G4)` from the resulting cycle **fluxes**
  `(J2,J4)`;
- a rank calculation at one numerical baseline from a generic structural
  identifiability result;
- steady-state snapshots from dynamic time-series identifiability; and
- a stoichiometric forcing signature from the projected parameter-to-observation
  map.

No paper located in this focused search gives the exact combination “formal
steady-state identifiability of AE2 and AE4 activities from physiologically
plausible outputs of the 2018 salivary acinar model.” Absence of that exact
application is not, by itself, mathematical novelty.

## What the current repository result actually establishes

The current machine-readable outputs are intentionally limited. They certify
the rank-two flux-coordinate map

\[
A=J_2+J_4,\qquad B=J_2+2J_4,
\]

with exact inverse

\[
J_4=B-A,\qquad J_2=2A-B,
\]

and the rank-two AE2/AE4 source signatures in the printed intracellular
balances, provided both transporter turnover factors are nonzero. Each scalar
aggregate alone has rank one. The outputs do **not** provide a reconstructed
physiological steady-state map, a secretion Jacobian, or a certified rank for
`Q` plus an intracellular ion. The full-model baseline remains blocked by
inconsistencies in the published executable specification.

This distinction is decisive. If `J_k=G_k phi_k(u)`, recovery of `(J2,J4)` does
not recover `(G2,G4)` unless the state-dependent factors `phi_k(u)` are known
and nonzero. Likewise, rank two of the unobserved forcing matrix does not imply
rank two after the steady-state dynamics and measurement projection. At a
regular steady state,

\[
D_\theta H=h_\theta-h_u F_u^{-1}F_\theta,
\]

so a claim about an observation panel must concern this projected Jacobian,
not the stoichiometric columns alone.

## Closest primary precedents

| Primary precedent | Closest prior result | Implication for this project |
| --- | --- | --- |
| Vera-Sigüenza et al. (2018), [DOI 10.1007/s11538-017-0370-6](https://doi.org/10.1007/s11538-017-0370-6), [open full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC5792321/) | The source salivary model already reports that AE2 knockout leaves flow nearly unchanged, AE4 knockout lowers flow, and loss of either exchanger increases the activity of the other. | Do not claim discovery of AE2/AE4 compensation, or merely that secretion can mask a transporter change. A formal equivalence geometry or certified minimal physiological panel could still be new. |
| Mori (2012), [DOI 10.1007/s00285-011-0483-8](https://doi.org/10.1007/s00285-011-0483-8) | Proves existence, uniqueness, stability, and global behavior for a general class of pump--leak models. | A “general pump--leak theorem” must be clearly inverse/observational and must add more than forward steady-state existence or stability. |
| Ouellet et al. (2025; online 2024), [DOI 10.1007/s00285-024-02163-z](https://doi.org/10.1007/s00285-024-02163-z) | Gives a formalism and conditions for existence and uniqueness of steady states in a broad pump--leak--cotransport class. | The nearby general mathematical terrain already includes cotransporters. A new class theorem must engage observation projection or activity recovery, not just transporter stoichiometry and steady-state uniqueness. |
| Catacuzzeno, Cavaliere, and Michelucci (2026), [DOI 10.1016/j.bpj.2026.01.021](https://doi.org/10.1016/j.bpj.2026.01.021), [PubMed](https://pubmed.ncbi.nlm.nih.gov/41536064/) | Introduces an explicitly inverse, steady-state, experimentally driven pump--leak construction. Intracellular Na, K, Cl, membrane potential, and volume constrain parameters; linear programming explores feasible ranges and progressively restricts them toward unique solutions. | This is the closest direct precedent. “First steady-state parameter inference for a pump--leak model” is false. A defensible contribution would need a formal AE2/AE4 equivalence/minimal-panel theorem that is not supplied by the linear-programming construction. |
| Tarrant, Kay, and Aminzare (2026), [DOI 10.1007/s11538-026-01683-7](https://doi.org/10.1007/s11538-026-01683-7), [PubMed](https://pubmed.ncbi.nlm.nih.gov/42446753/) | Develops a two-compartment epithelial pump--leak model, derives analytical steady states, and studies stability, global sensitivity, and robustness. | An epithelial geometry and analytical steady-state formulas are no longer a novel setting. This paper does not appear to present inverse identifiability, leaving a narrow opening for a genuine observation theorem. |
| Anguelova, Karlsson, and Jirstrand (2012), [DOI 10.1016/j.mbs.2012.04.005](https://doi.org/10.1016/j.mbs.2012.04.005), [PubMed](https://pubmed.ncbi.nlm.nih.gov/22609467/) | Establishes methodology for finding minimal output sets for identifiability. | Selecting a smallest full-rank output panel is an established problem, not a new mathematical concept. Novelty must come from a nontrivial model-specific classification or theorem. |
| Joubert, Stigter, and Molenaar (2018), [DOI 10.1371/journal.pone.0207334](https://doi.org/10.1371/journal.pone.0207334) | Gives an iterative structural-identifiability algorithm that finds potentially multiple minimal output sets, including restrictions to measurable candidates. | A numerical enumeration of `Q` plus ion panels is routine unless it exposes and proves a new structural obstruction or complete model-specific classification. |
| Haus, Drengstig, and Thorsen (2023), [DOI 10.1371/journal.pcbi.1011398](https://doi.org/10.1371/journal.pcbi.1011398) | Studies 3,648 concentration/flow measurement combinations across 128 biomolecular controller motifs. For those motifs, one measurement is insufficient and measurements related to different species are necessary; the paper explicitly discusses transmembrane ion-flow measurements and epithelial Na/K homeostasis. | Broad claims about the first identifiability study involving ion flows, flow outputs, or epithelial homeostasis are untenable. Its data are dynamic time series, so a salivary steady-state-only theorem could remain distinct. |
| Csercsik, Hangos, and Szederkényi (2012), [DOI 10.1016/j.neucom.2011.09.006](https://doi.org/10.1016/j.neucom.2011.09.006); Walch and Eisenberg (2016), [DOI 10.1016/j.neucom.2016.03.027](https://doi.org/10.1016/j.neucom.2016.03.027) | Derive structural nonidentifiability and identifiable combinations for Hodgkin--Huxley-type ion-channel models and use those results to inform voltage-clamp design. | “First structural identifiability analysis of an ion-transport model” is false. These are ion-channel dynamics rather than a coupled epithelial transporter network, so the exact biological scope differs. |
| Gross and Blüthgen (2020), [DOI 10.1093/bioinformatics/btaa404](https://doi.org/10.1093/bioinformatics/btaa404) | Assumes relaxation to a stable steady state, applies the implicit function theorem to obtain the steady-state response derivative `-J^{-1}S`, uses rank to quantify solution-space dimension, and designs perturbations using graph/matroid structure. | The IFT steady-state sensitivity formula, rank restoration by perturbation, and minimal perturbation design are established. A paper here needs an AE2/AE4-specific or pump--leak-specific theorem beyond reusing this template. |
| Villaverde et al. (2019), [DOI 10.1109/LCSYS.2018.2868608](https://doi.org/10.1109/LCSYS.2018.2868608), [author copy](https://wrap.warwick.ac.uk/id/eprint/110005/) | Characterizes input-dependent structural identifiability and shows that a time-varying input can sometimes be replaced by multiple experiments with constant inputs. | “One condition unidentifiable, two conditions identifiable” is not a new general principle. Novelty would require a sharp minimal perturbation result tied to this transport network and experimentally available controls. |
| Ligon et al. (2018), [DOI 10.1093/bioinformatics/btx735](https://doi.org/10.1093/bioinformatics/btx735) | Implements multi-experiment structural-identifiability analysis for nonlinear SBML models. | Merely stacking observation Jacobians from several conditions is standard methodology. |
| Ovchinnikov et al. (2022), [DOI 10.1137/21M1389845](https://doi.org/10.1137/21M1389845), [preprint](https://arxiv.org/abs/2011.10868) | Gives an algorithm for the exact number of experiments needed for multi-experiment local identifiability and an almost-exact global bound for ODE models. | Any broad claim about determining a minimal number of experiments is precluded. A useful result here must exploit the AE2/AE4 or pump--leak structure. |
| Frøysa, Skaug, and Alendal (2020), [DOI 10.1016/j.mbs.2019.108291](https://doi.org/10.1016/j.mbs.2019.108291) | Treats pointwise and global identifiability, minimal measurements, modified stoichiometric matrices, and D-optimal input design for steady-state linear metabolic networks. | This is a close precedent for a stoichiometric steady-state rank/design result. Independent transporter signature columns alone are insufficient novelty. |
| Palombo et al. (2023), [DOI 10.3390/sym15020368](https://doi.org/10.3390/sym15020368) | Analyzes kinetic-parameter recovery from multiple nonlinear steady-state measurements, including Jacobian nonsingularity/local identifiability and information-count conditions. | Multiple nonlinear steady states plus a nonsingular Jacobian is established territory. A new result must be sharper than a condition count or numerical determinant. |
| Brehm and Fiedler (2018), [DOI 10.1002/mma.4668](https://doi.org/10.1002/mma.4668), [preprint](https://arxiv.org/abs/1606.00279) | Derives qualitative steady-state responses to rate perturbations from reaction-network stoichiometry under generic nondegeneracy conditions. | Stoichiometry-based perturbation signatures are established. A new criterion must address what survives conservation, kinetics, and partial observation in a pump--leak network. |
| Craciun and Pantea (2008), [DOI 10.1007/s10910-007-9307-x](https://doi.org/10.1007/s10910-007-9307-x), [author copy](https://people.math.wisc.edu/~craciun/PAPERS_NEW/Craciun_Pantea_IdentifChem_2008.pdf) | Treats identifiability of chemical reaction networks as a structural network problem. | A broad “first stoichiometric/network identifiability criterion” claim is unsafe. Any theorem here must state its pump--leak assumptions and observation operator precisely. |

## Claim-level novelty assessment

### Defensible now

- **Search statement:** “This focused audit located no prior formal
  characterization of steady-state AE2/AE4 activity identifiability in the
  published salivary acinar model.” This is a report of the search, not proof of
  priority.
- **Narrow algebraic statement:** the printed AE2/AE4 balances have independent
  full-state forcing signatures, and the two flux aggregates `(A,B)` invert the
  two exchanger cycle fluxes. This can be presented as a checked,
  model-specific lemma, not as a major identifiability theorem.
- **Limitation statement:** secretion alone cannot generically separate two
  independent activities at a regular point, and the current published record
  is insufficient to certify a physiological two-output panel. The latter is a
  useful model-audit conclusion even though the dimension argument is standard.

### Potentially defensible only after additional proof

- “First model-specific characterization of the local AE2/AE4 equivalence
  manifold for physiologically plausible steady-state measurements in the
  salivary acinar model,” if the actual observation map is reconstructed and
  the manifold is proved or certified over a stated regular domain.
- A complete minimal-panel classification for a fixed, biologically measurable
  candidate set, if generic full rank or structural rank loss is proved rather
  than checked at one baseline.
- A concrete perturbation theorem identifying an experimentally controllable
  second condition that restores rank, including the smallest required
  perturbation set and nuisance assumptions.
- A general pump--leak criterion only if it concerns the projected matrix
  `h_u F_u^{-1} F_theta` (and any direct term), supplies necessary or sufficient
  conditions under stated nondegeneracy assumptions, and is not merely the
  independence of raw stoichiometric columns.

### Not defensible

- First inverse or steady-state parameter-inference method for pump--leak
  models.
- First identifiability analysis of ion transport, transporter networks, flow
  measurements, or epithelial homeostasis.
- Novelty of the generic fact that one regular scalar output leaves a
  one-dimensional level set in a two-parameter problem.
- Novelty of the implicit-function sensitivity formula, local Jacobian-rank
  test, minimal-output-set problem, or stacking multiple experimental
  conditions.
- Structural identifiability of `(G2,G4)` from rank two of the raw activity
  signatures or from the invertible `(A,B)` flux map.
- A certified `Q`-plus-ion minimal panel from the current outputs; its
  physiological observation Jacobian has not been reconstructed.
- Global injectivity from a local rank calculation or a sampled parameter box.
- Mechanism discrimination from unequal knockout phenotypes alone. The 2018
  model already reports those phenotypes and compensation, and discrimination
  requires nonoverlapping predicted observation sets over the stated parameter
  classes.
- A new two-condition theorem whose only content is that two scalar gradients
  are noncollinear, or a new stoichiometric theorem whose only content is that
  two columns are linearly independent.

## Mathematical wording constraints

At a regular point of a scalar map `H:R^2 -> R`, where `grad H != 0`, the
regular-level-set theorem gives a locally one-dimensional observational
equivalence curve. A rank-zero point does not by itself imply a two-dimensional
fiber; for example, `H(x,y)=x^2+y^2` has an isolated zero fiber. Fiber dimension
claims therefore require regularity or a locally constant-rank hypothesis.

For two outputs, a nonzero determinant of `D_theta H` is sufficient for a local
inverse. For more than two outputs, full column rank two gives local
identifiability. A numerical rank at one physiological point is only a
pointwise local result. Calling it **structural** requires showing that an
appropriate minor is generically nonzero (for example, not identically zero on
the admissible model class). Raw determinants and singular values also depend
on units; panel comparisons should use dimensionless/log sensitivities or a
noise-normalized information matrix.

## Viability implication: `STOP`

The novelty audit supports **`STOP` for the proposed paper in its current
state**.

1. The only completed exact result is a two-by-two linear inversion in
   unobserved flux coordinates, plus independent source signatures in the full
   state balances. It does not establish activity identifiability from
   physiological measurements.
2. The scalar-output dimension result is standard differential topology, and
   minimal-output/rank and multi-condition designs already have substantial
   prior literatures.
3. The 2026 Catacuzzeno et al. paper is a direct steady-state inverse
   pump--leak precedent, while the 2026 Tarrant et al. paper occupies nearby
   analytical epithelial pump--leak territory.
4. The 2018 salivary paper already contains the qualitative AE2/AE4
   compensation and knockout-flow observations that might otherwise motivate
   the result biologically.
5. No physiological `H(G2,G4)`, full-rank biologically plausible measurement
   panel, perturbation-induced rank restoration, global equivalence theorem, or
   separated mechanism image sets have been validated. In particular, no
   physiological minimal-panel result is frozen in `docs/RESULTS_LEDGER.md`.

This classification does not assert that a stronger result is impossible. It
means that drafting a paper from the current result would require presenting
standard mathematics or an unmeasured flux inversion as the headline. Reopen
the publication decision only if at least one of the following is completed:

- an independently validated physiological steady-state map and a proved
  AE2/AE4 equivalence or minimal-panel theorem;
- a concrete, experimentally plausible perturbation that provably changes
  rank one to rank two and remains discriminating with stated nuisance
  parameters; or
- a genuinely general observation-level pump--leak theorem that goes beyond
  raw stoichiometric column rank and is shown not to be covered by the
  precedents above.
