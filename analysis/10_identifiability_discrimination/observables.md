# Observables, flux coordinates, and local identification criteria

## Scope and status

This note separates three objects that must not be conflated:

1. exact aggregate coordinates of the Ae2 and Ae4 cycle fluxes;
2. physiological steady-state outputs such as secretion rate and intracellular ion concentrations; and
3. an illustrative stoichiometric map used to display the exact aggregate geometry.

The first object is certified directly from the published balance stoichiometry. The second requires a validated steady-state map for the full cell model and is not currently numerically available. The third is useful for explaining the algebra, but it is not a reconstructed physiological model and does not establish a measurement recommendation.

Throughout, let

\[
F(u;\theta,\xi)=0,\qquad \theta=(G_2,G_4),
\]

where `u` is a steady state, `G_2` and `G_4` are Ae2 and Ae4 activities, and `xi` denotes a possible controlled condition. A physiological observation panel is

\[
H(\theta,\xi)=h(u^*(\theta,\xi),\theta,\xi).
\]

No archived code, derivation, prose, or figure was used for this note.

## Certified transporter forcing signatures

In the intracellular balance ordering

\[
(\mathrm{Cl}_i,\mathrm{Na}_i,\mathrm{K}_i,\mathrm{HCO}_{3,i},
  \mathrm{H}_i,\mathrm{CO}_{2,i}),
\]

one forward Ae2 cycle and one forward Ae4 cycle have the published signatures

\[
s_2=(1,0,0,-1,0,0)^T,
\]

and

\[
s_4=(1,-\alpha,-\beta,-2,0,0)^T,
\qquad
\alpha=\frac{[\mathrm{Na}]_i}{[\mathrm{Na}]_i+[\mathrm{K}]_i},
\quad
\beta=1-\alpha.
\]

At the printed resting concentrations, `alpha = 25/145` and `beta = 120/145`. The two columns are linearly independent whenever the intracellular monovalent-cation sum is positive. In particular, their chloride/bicarbonate minor is

\[
\det\begin{pmatrix}1&1\\-1&-2\end{pmatrix}=-1.
\]

If activities are parameterized logarithmically and the turnover laws are valid, their formal columns in `F_theta` are

\[
F_{\log G_2}=J_2s_2,
\qquad
F_{\log G_4}=J_4s_4.
\]

The two formal columns are independent only where both cycle fluxes are nonzero. Log-activity coordinates are not defined at a knockout; activity multipliers must be used there instead. These are statements about how the transporter cycles enter the balance equations, not yet statements about the response of a measured output: feedback through all other transporters is contained in `F_u^{-1}`.

## Certified aggregate flux coordinates

Define

\[
A=J_2+J_4,
\qquad
B=J_2+2J_4.
\]

`A` is the combined Ae-mediated chloride-cycle contribution and `B` is the corresponding bicarbonate-cycle demand. Algebraically,

\[
\binom{A}{B}
=
\begin{pmatrix}1&1\\1&2\end{pmatrix}
\binom{J_2}{J_4}.
\]

The determinant is one, so the map is globally invertible:

\[
J_2=2A-B,
\qquad
J_4=B-A.
\]

For nonnegative forward cycle rates, the image is exactly the wedge

\[
A\geq 0,
\qquad
A\leq B\leq 2A.
\]

Thus the pair `(A,B)` discriminates the two cycle rates exactly. Either scalar aggregate alone cannot do so: its level sets are straight one-dimensional equivalence sets. These conclusions are certified algebraic statements.

However, `A` and `B` are internal flux combinations, not established independent experimental measurements in the 2018 model. The exact inverse therefore does not by itself solve the physiological minimal-measurement problem.

## Local Jacobian criterion for physiological panels

Where `F_u` is nonsingular, the implicit function theorem gives

\[
D_\theta u^*=-F_u^{-1}F_\theta
\]

and

\[
D_\theta H
=h_u(-F_u^{-1}F_\theta)+h_\theta.
\]

The correct local criteria are:

- a panel locally identifies `(G_2,G_4)` only if `rank(D_theta H) = 2`;
- one scalar steady-state output has rank at most one and cannot generically locally identify two independent activities;
- at a regular rank-one point, the constant-rank theorem gives a local one-dimensional observational-equivalence set;
- at rank zero, both parameter directions are invisible to first order; and
- for a panel of two or more outputs, at least one scaled `2 x 2` minor must be nonzero.

Rank must be assessed after fixing dimensionless parameter coordinates and defensible observation scales or noise units. A raw determinant comparing quantities with unrelated units is not a biological measure of information.

## Why `Q` plus an ion concentration is not certified

Candidate panels include

\[
(Q,[\mathrm{Na}]_i),\quad
(Q,[\mathrm{K}]_i),\quad
(Q,[\mathrm{Cl}]_i),\quad
(Q,[\mathrm{HCO}_3]_i),\quad
(Q,\omega_i).
\]

None can currently be certified as full rank for the published AE4 cell model. Doing so requires a validated baseline `u*`, a nonsingular and dimensionally coherent `F_u`, and a verified secretion observable. The printed record presently fails those prerequisites:

- Eqs. (10)--(11), evaluated at the paper's reported osmolarities and permeabilities, give `q_b/q_a` about 59.6 although either printed volume-balance sign requires `q_a = q_b` at steady state.
- The printed Table 1 voltages, lumen concentrations, resting calcium value, and conductances leave an apical QSS residual of about 67.8 pA.
- The Ae4 fourth-order rate has an unresolved M-versus-mM numerical convention, changing its turnover by a factor of `10^12`.
- The CO2 aggregate sign conflicts with the two individually printed membrane influxes; the acid-base rate units and lumen closure are incomplete.
- No machine-readable stimulated steady state or calcium input is published.

Consequently, substituting the transporter signatures directly for `D_theta H` would omit the factor `-F_u^{-1}` and would be mathematically incorrect. The present result set records `Q` alone only as a generic rank-at-most-one case and records every `Q`-plus-ion panel as not certified.

## Two-condition and perturbation criteria

Let a scalar secretion measurement at condition `xi` have log-activity gradient

\[
g(\xi)=D_{(\log G_2,\log G_4)}Q(\xi).
\]

Two scalar measurements, at `xi_0` and `xi_1`, locally identify the activities precisely when

\[
\det
\begin{pmatrix}
g(\xi_0)\\
g(\xi_1)
\end{pmatrix}
\neq 0.
\]

For a small one-parameter perturbation, the sharp first-order rotation test is

\[
\det\begin{pmatrix}
g(\xi_0)\\
\partial_\xi g(\xi_0)
\end{pmatrix}\neq 0.
\]

For vector observation panels, replace the two rows by the corresponding stacked observation Jacobians and require column rank two. These criteria show exactly what would need to be calculated for extracellular-ion changes, agonist level, or partial inhibition. They do not assert that any published perturbation passes the test.

## Derivative verification requirement

If a coherent residual is obtained later, an implicit derivative must be checked against independently re-solved equilibria:

1. compute `F_u` and `F_theta` in dimensionless coordinates;
2. evaluate `-solve(F_u, F_theta)`;
3. re-solve the steady state at centered parameter perturbations, using branch-continuous warm starts;
4. compare the resulting finite-difference Jacobian over a decreasing step sequence; and
5. report the condition number of `F_u` and both absolute and relative derivative errors.

The current tests verify the software identity on an independent nonlinear toy equilibrium. They do not validate an implicit derivative of the physiological AE4 model.

## Scope of the generated geometry

`results/10_identifiability_discrimination/stoichiometric_geometry.csv` and `figures/identifiability_stoichiometric_geometry.png` plot the exact linear map

\[
(J_2,J_4)\longmapsto(A,B)
\]

over normalized nonnegative cycle rates. The plot correctly displays the forward-flux wedge, the scalar equivalence lines, and the separation obtained by retaining both aggregate coordinates.

It is an illustrative stoichiometric figure. Its axes are normalized internal cycle rates and aggregate flux coordinates, not `(G_2,G_4)` and physiological observations. Because the map is affine, the figure cannot establish or exclude folds, self-intersections, near-overlaps, or rank changes in the nonlinear physiological observation map.

## Current minimal-measurement conclusion

The following statements are supported:

- one scalar steady-state observation cannot generically identify two activities;
- the exact pair of internal aggregates `(A,B)` identifies `(J_2,J_4)`;
- Ae2 and Ae4 have independent stoichiometric forcing signatures; and
- a second experimental condition is useful exactly when its sensitivity row is not collinear with the first.

The following statements are not supported:

- that `Q` and any particular ion concentration form a minimal identifying physiological panel;
- that a one- or two-condition secretion experiment identifies `(G_2,G_4)` in the published cell model;
- that the selected physiological observation map is globally injective; or
- that the illustrated aggregate geometry is a reconstruction of the 2018 model.
