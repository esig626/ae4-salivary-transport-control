# Targeted novelty audit

Search date: 2026-08-27.

## Exact candidate claims searched

1. A regular steady-state observation Jacobian is obtained from the implicit-function derivative, and full column rank gives local identifiability.
2. One scalar steady-state observation generically leaves a one-dimensional equivalence set for two parameters.
3. Additional observables or experimental conditions can restore rank; two scalar conditions suffice exactly when their stacked sensitivity rows are independent.
4. Competing mechanism classes are discriminable only when experimental design separates their predicted observation sets.
5. The projected transporter-stoichiometry factorization might constitute a new pump--leak-specific theorem.

Queries combined terms such as `steady state`, `implicit function theorem`, `observation Jacobian`, `structural identifiability`, `multiple experiments`, `model discrimination`, `pump-leak`, `epithelial`, `AE2`, and `AE4`. Searches were scoped to primary papers and official journal pages. No search result was treated as proof of absence.

## Closest prior art

| Source | Relevant established result | Effect on novelty |
| --- | --- | --- |
| Raue et al. (2014), [doi:10.1093/bioinformatics/btu006](https://doi.org/10.1093/bioinformatics/btu006) | Structural/practical identifiability methods; identifiability analysis guides model simplification and additional experiments | Rank/experiment recommendations are established systems-biology practice |
| Kiss et al. (2023), [doi:10.1007/s11538-023-01121-y](https://doi.org/10.1007/s11538-023-01121-y) | Uses an observability-identifiability Jacobian and the implicit function theorem for local unique parameter determination; discusses characteristic steady/asymptotic observations | The IFT/full-rank core is explicitly prior art in the same journal family |
| Villaverde et al. (2022), [doi:10.1093/bib/bbab387](https://doi.org/10.1093/bib/bbab387) | Recommends resolving nonidentifiable manifolds and adding observables/conditions; gives an example where multiple experiments restore structural identifiability | “One experiment insufficient, multiple experiments sufficient” is not new in general |
| Myung and Pitt (2009), [doi:10.1037/a0016104](https://doi.org/10.1037/a0016104) | Formal optimal experimental design for model discrimination | Observation-set separation/design is an established research area |
| Vanlier et al. (2014), [doi:10.1186/1752-0509-8-20](https://doi.org/10.1186/1752-0509-8-20) | Optimal experiment design for selecting among biochemical-network models | A salivary mechanism-discrimination calculation would be an application unless it revealed special geometry |
| Vera-Sigüenza et al. (2018), [doi:10.1007/s11538-017-0370-6](https://doi.org/10.1007/s11538-017-0370-6) | Published AE2/AE4 salivary model and knockout predictions | Supplies the application, not an identifiability theorem |
| Tarrant, Kay, and Aminzare (2026), [arXiv:2601.08975](https://arxiv.org/abs/2601.08975) | Analytical equilibria, stability, and robustness for a generalized epithelial cell/lumen pump--leak model | Shows active current work on general pump--leak structure; does not make this audit's IFT rank factorization novel |

The salivary-specific searches found the 2018 model, experimental AE4/AE2 papers, and other salivary transport models, but no paper centered on identifying `G2,G4` from steady-state measurement panels. That application-level gap is not enough: the candidate mathematical results are standard, and the physiological instance has not passed reproduction.

## Novelty verdict

The factorization

`J_h = -h_u F_u^{-1}[s2 phi2, s4 phi4]`

is a useful explanatory organization of the salivary inverse problem, but it is direct IFT/linear algebra. The scalar equivalence curve, full-rank panel criterion, and stacked-condition determinant are standard consequences of established theory. No new global, pump--leak-specific, or testing result was found.

Verdict: **no mathematical novelty sufficient for a quick paper in the present state**. A future paper would need a reproducible physiological map plus a special global/perturbational result, not merely this rank calculation.
