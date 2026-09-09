# Mandatory escalation scan

## Decision summary

The scan found a useful exact stoichiometric decomposition but no validated route from transporter activities to physiological steady-state observations. The strongest certified result is the globally invertible internal-flux map

\[
(J_2,J_4)\longmapsto(A,B)=(J_2+J_4,J_2+2J_4).
\]

This is elementary balance algebra, not by itself a publishable model-specific identifiability result. Every proposed escalation that depends on the nonlinear physiological steady-state map is blocked by the absence of a dimensionally coherent, baseline-consistent reconstruction from the published record.

The evidence-based classification is therefore:

> **STOP** for the requested model-specific quick paper under the currently available published specification.

This is not evidence that Ae2 and Ae4 are biologically indistinguishable. It means that the requested mathematical claims cannot be established from the present executable evidence without introducing unvalidated modelling choices.

## Evidence gate

The scan uses only the 2018 published model, the Phase 00 source map and blocker record, and new independent algebra and code. It does not use unpublished archived code, derivations, text, or figures.

The exact certified pieces are:

- the Ae2 and Ae4 stoichiometric forcing columns;
- the rank-two aggregate map `(J2,J4) -> (A,B)` and its inverse;
- the forward-flux wedge `A >= 0`, `A <= B <= 2A`;
- the generic rank bound for a scalar observation of two parameters; and
- the local and two-condition Jacobian criteria.

The numerical audit is a diagnostic evaluation of the printed equations, not a reproduced full model.

## Published-model blockers affecting every nonlinear escalation

| Blocker | Direct consequence |
| --- | --- |
| Eq. (9) has opposite dynamic signs in two locations | At steady state the conflict disappears because both require `q_a = q_b`; it cannot explain the remaining baseline mismatch. |
| Reported osmolarities and permeabilities give `q_b/q_a` about 59.6 | The printed resting row is not a steady cell-volume state under Eqs. (10)--(11). |
| Printed resting electrical data leave about 67.8 pA apical QSS residual | The reported voltages and ionic row cannot be used as a consistent algebraic initialization. |
| AE4 fourth-order kinetics combine a printed mM-rate law with inconsistent density/amount prose and an unknown whole-cell conversion | Reusing the printed constants with M-valued states changes the turnover by `10^12`; direct mM use is literal, but `F_theta` still lacks a commensurate whole-cell scale. |
| The upstream Palk corrigendum makes NKCC coefficients source-correctable, but no equivalent NaK correction and no fitted-density conversion are published | NKCC's sign is repaired by an explicit conversion; the complete residual still requires documented conventions for the other terms. |
| CO2 aggregate sign is the negative of the sum of the two printed membrane influxes | Acid-base steady-state direction is ambiguous. |
| Buffer reverse-rate and CO2 permeability dimensions do not close the amount balances as printed | Flux magnitudes cannot be combined in one residual without an added convention. |
| Lumen acid-base closure and consistent initial conditions are absent | The reported ten-state QSS system is not closed as an executable specification. |
| The stimulated calcium drive and stimulated steady state are not machine-readable | The secretion regime relevant to the Ae4 phenotype cannot be reconstructed from a numeric table. |

The direct audit also shows that straightforward unit conversions produce flux scales separated by many orders of magnitude. Fitting hidden scale factors to force the printed baseline would make that baseline true by construction and would not constitute independent reproduction.

## A. Nonlinear observation geometry and global injectivity

### Question examined

Does a physiological map such as

\[
(G_2,G_4)\mapsto Q
\quad\text{or}\quad
(G_2,G_4)\mapsto(Q,[X]_i)
\]

have folds, self-intersections, exact equivalence sets, or a globally injective image on an admissible parameter region?

### What is supported

For the internal aggregate map,

\[
(J_2,J_4)\mapsto(A,B),
\]

global injectivity is exact because the determinant is one. Each scalar aggregate has exact straight equivalence lines. For nonnegative forward cycle rates, the full image is the wedge `A >= 0`, `A <= B <= 2A`.

### What is not supported

The aggregate map is not the nonlinear activity-to-physiology map. No validated `u*(G2,G4)` is available, so there is no defensible computation of physiological folds, self-intersections, near-overlaps, image boundaries, singular sets, or global injectivity. Numerical sampling of a fitted surrogate would answer a different model and cannot repair this gap.

No global physiological identifiability theorem is claimed.

## B. Perturbation-induced distinguishability

### Sharp criterion

For a scalar secretion measurement under controlled condition `xi`, define

\[
g(\xi)=D_{(\log G_2,\log G_4)}Q(\xi).
\]

Baseline and perturbed secretion measurements have local rank two exactly when

\[
\det[g(\xi_0);g(\xi_1)]\neq0.
\]

For a small perturbation, the first-order criterion is

\[
\det[g(\xi_0);\partial_\xi g(\xi_0)]\neq0.
\]

### Candidate perturbations considered

- extracellular chloride or bicarbonate;
- extracellular sodium or potassium;
- agonist/calcium level; and
- selective partial inhibition of one exchanger.

These quantities occur in the published equations, so they are structurally meaningful candidates. However, the paper does not supply validated perturbation ranges, a numerical calcium function, or a coherent baseline Jacobian. The required rows `g(xi)` therefore cannot be computed or compared.

The scan supports the criterion, not the claim that one extra condition succeeds.

## C. Experimental design and uncertainty

A defensible local design would compare, in noise-scaled coordinates,

- a `2 x 2` determinant;
- the smallest singular value;
- or the determinant of a Fisher information matrix formed from stacked conditions.

The scale of each physiological observation must be set by explicit measurement noise or a stated nondimensionalization. Raw singular values of secretion and millimolar concentrations are not comparable.

No optimization was performed because the prerequisites are absent. In particular:

- the AE4 unit ambiguity alone changes a formal turnover by `10^12`;
- the baseline is not an algebraic steady state under the printed data;
- the admissible activity and perturbation domains are unpublished; and
- local derivatives have not been verified against re-solved physiological equilibria.

An optimizer applied after choosing hidden scale factors would optimize those choices, not the published model. No optimal intervention or minimal perturbation count is claimed.

## D. General pump--leak structure

The scan found one modest general statement. If two transporter activities enter a balance system through independent stoichiometric columns and the full steady-state Jacobian is nonsingular, then local output identification depends on whether the selected observation rows preserve two independent columns of

\[
-F_u^{-1}F_\theta.
\]

For the two exchangers here, the cycle signatures are independent because AE4 has an additional monovalent-cation signature and a second bicarbonate equivalent. The internal aggregate pair `(A,B)` consequently separates the cycle fluxes exactly.

This does not imply that adding one arbitrary ion concentration to secretion must identify transporter activities. Network feedback can rotate, attenuate, or align the activity columns after multiplication by `F_u^{-1}` and projection through `h_u`. A useful general physiological theorem would require assumptions on the leak Jacobian and observation operator that have not been established.

The present algebra is too elementary and too remote from measured outputs to justify `PUSH HARDER` as a general pump--leak theorem.

## E. Mechanism variants and hypothesis testing

The pure forward Ae2 and Ae4 rays in aggregate space are

\[
B=A\quad\text{and}\quad B=2A,
\]

respectively. They intersect only at zero in the two-dimensional `(A,B)` coordinate system. This is exact separation of idealized internal cycle types.

It is not a physiological mechanism-discrimination result because `(A,B)` has not been established as a measured panel. Under a scalar projection, unknown cycle amplitudes can overlap. Under physiological outputs, the two class images cannot be computed until the full model and parameter classes are defined.

Accordingly, the scan found no justified positive minimum separation, least-favourable pair, Gaussian error calculation, or nontrivial composite testing problem. Adding Renyi or finite-blocklength theory would not solve the missing observation map and is not warranted.

## F. Inference-relevant singularities

No fold, branch collision, bifurcation, or rank-changing physiological singularity has been demonstrated. The aggregate surrogate is affine and therefore has none. A continuation package would merely continue an unvalidated model choice, so AUTO or MATCONT is not justified.

If a corrected residual becomes available, the first diagnostic should be a branch-continuous steady-state solve with monitoring of

\[
\det F_u,\quad \sigma_{\min}(F_u),\quad
\text{and}\quad \sigma_{\min}(D_\theta H),
\]

not a general bifurcation catalogue.

## Systematic escalation status

| Direction | Status | Reason |
| --- | --- | --- |
| Nonlinear physiological geometry | Not established | No validated steady-state map. |
| Global injectivity | Certified only for internal `(A,B)` coordinates | The physiological map is unavailable. |
| Exact equivalence | Certified for either scalar aggregate | No corresponding exact result for `Q`. |
| Extra perturbations | Criterion established; success untested | Perturbation Jacobians cannot be computed. |
| Optimal design | Not attempted | Noise scales, domain, and sensitivities are missing. |
| General pump--leak theorem | No escalation | Only routine stoichiometric rank algebra is supported. |
| Mechanism variants | Separated only in idealized aggregate coordinates | Physiological class images are unavailable. |
| Hypothesis testing | Not justified | No validated observation-set separation. |
| Inference-relevant singularity | None found | No coherent nonlinear branch is available. |
| Uncertainty robustness | Fails at specification level | Unit and closure choices dominate any numerical uncertainty analysis. |

## Decision and restart conditions

The current result is insufficient for a short mathematical biology paper about physiological identification of Ae2 and Ae4 activities. It supplies a correct generic scalar-observation limitation and an exact internal-flux decomposition, but not a model-specific minimal panel, perturbation result, nonlinear observation geometry, or discrimination analysis.

The project should be restarted only after at least the following are obtained from a published correction, supplement, or author-confirmed executable specification:

1. one coherent unit ledger for every amount flux, current, water flux, and density;
2. a complete resting or stimulated algebraically consistent steady state, including acid-base closure and membrane potentials; and
3. a numerical calcium/input condition and an independently reproducible baseline/knockout check.

After that gate, the next three technical steps are:

1. solve the reduced steady branch and verify all residuals and conservation identities;
2. verify implicit derivatives against finite differences of independently re-solved equilibria; and
3. evaluate single-panel and two-condition rank determinants over a source-supported parameter and perturbation domain.

Until those conditions are met, the defensible classification remains `STOP`, and the generated aggregate figure must remain explicitly labeled illustrative.
