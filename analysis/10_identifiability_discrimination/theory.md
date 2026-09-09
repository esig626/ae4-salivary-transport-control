# Identifiability theory for the steady-state AE2/AE4 problem

## Scope and status

This note records the analytical conclusions that follow from the published
steady-state balances without requiring a numerically reconstructed baseline.
The unknown activity parameters are

\[
\theta=(G_2,G_4),
\]

where \(G_2\) and \(G_4\) multiply the AE2 and AE4 constitutive driving terms.
The results below distinguish four logically different questions:

1. regular local identifiability of a chosen observation map;
2. exact recovery of transporter fluxes from balance summaries;
3. uniqueness of the activities from a complete steady state;
4. global or singular behavior of a partial observation map.

Only the first three admit analytical statements from the currently available
published equations. No model-specific partial panel such as
\((Q,[\mathrm{Na}]_i)\), \((Q,[\mathrm{K}]_i)\), or
\((Q,[\mathrm{Cl}]_i)\) is numerically certified in this repository. In
particular, the structural rank results below must not be reported as a
numerical rank result for one of those panels.

The proofs use no historical unpublished derivation or code. They use the
published intracellular and luminal balances recorded in `model/equations.md`.

## Assumptions and notation

Fix an experimental condition \(\xi\), such as the bath composition and
prescribed calcium level. Let

\[
F(u;\theta,\xi)=0
\]

be a square steady-state system with state \(u\in U\subset\mathbb R^n\). The
state may include algebraic variables, provided the selected square system has
an invertible Jacobian \(F_u\) at the point under study. Unless direct parameter
dependence is written explicitly, an observation is

\[
H(\theta;\xi)=h(u^*(\theta;\xi),\xi).
\]

The activity parameterization is assumed to be genuinely multiplicative:

\[
J_2=G_2\varphi_2(u;\xi),\qquad
J_4=G_4\varphi_4(u;\xi).
\]

Here \(J_2=J_{\rm Ae2}\) and \(J_4=J_{\rm Ae4}\) are signed transporter cycle
fluxes in the common amount-per-time convention required by the intracellular
balances. This common convention is an assumption of the symbolic theory, not
a resolution of the published AE4 unit inconsistency. Numerical use of the
formulas requires that inconsistency to be resolved first.

For the AE4 cation partition, define

\[
\alpha(u)=\frac{[\mathrm{Na}]_i}
{[\mathrm{Na}]_i+[\mathrm{K}]_i},\qquad
\beta(u)=\frac{[\mathrm{K}]_i}
{[\mathrm{Na}]_i+[\mathrm{K}]_i},
\]

so that \(\alpha+\beta=1\) whenever the denominator is nonzero.

## Regular local identifiability

### Proposition 1. Implicit sensitivity and regular fibers

Suppose \(F\) and \(h\) are continuously differentiable,

\[
F(u_0;\theta_0,\xi)=0,
\qquad \det F_u(u_0;\theta_0,\xi)\ne0.
\]

Then a unique local steady branch \(u^*(\theta;\xi)\) exists and

\[
D_\theta u^*=-F_u^{-1}F_\theta.
\]

For an observation that may also depend directly on \(\theta\),

\[
D_\theta H
=h_\theta-h_uF_u^{-1}F_\theta.
\]

If \(D_\theta H\) has constant rank \(r\) in a neighborhood of \(\theta_0\),
the local fiber

\[
\{\theta:H(\theta;\xi)=H(\theta_0;\xi)\}
\]

is a smooth submanifold of dimension \(2-r\).

**Proof.** The derivative formulas are the implicit-function theorem and the
chain rule. The dimension statement is the constant-rank theorem. \(\square\)

### Consequences and critical-point caveat

At a regular point of a scalar observation, \(r=1\), so its local equivalence
set is a one-dimensional curve. Thus a single scalar steady-state observable
cannot *regularly* identify two independent activities.

For a vector observation with two or more components, rank two is sufficient
for local injectivity. For a two-component panel \(H=(H_1,H_2)\), the regular
criterion is

\[
\Delta_H
=\frac{\partial H_1}{\partial G_2}
  \frac{\partial H_2}{\partial G_4}
-\frac{\partial H_1}{\partial G_4}
  \frac{\partial H_2}{\partial G_2}
\ne0.
\]

Rank failure at a critical point does not by itself prove pointwise
nonidentifiability. For example, the scalar map
\(H(x,y)=x^2+y^2\) has an isolated zero fiber even though its derivative
vanishes at the origin. Such critical identification is nongeneric and not
robust under perturbation. Accordingly, all rank claims in this note are
claims of regular, first-order local identifiability unless explicitly stated
otherwise.

## Activity-linear balance rank

### Proposition 2. Full-state local rank and exact common-state uniqueness

Suppose, at a fixed condition and with all nuisance parameters fixed,

\[
F(u;\theta)=f_0(u)+A(u)\theta.
\]

Then:

1. if \(F_u\) is invertible at a steady state,
   \[
   \operatorname{rank}D_\theta u^*
   =\operatorname{rank}A(u^*);
   \]
2. if \(A(u)\) has full column rank, the same complete state \(u\) cannot be a
   steady state for two distinct activity vectors;
3. if one activity vector \(\theta_0\) makes a fixed state \(u\) steady, all
   activity vectors making that same state steady are exactly
   \[
   \theta_0+\ker A(u),
   \]
   intersected with the admissible parameter domain.

**Proof.** At a regular steady state,

\[
D_\theta u^*=-F_u^{-1}A(u^*),
\]

and multiplication by an invertible matrix preserves rank. If the same \(u\)
solves the equations at \(\theta\) and \(\widetilde\theta\), subtraction gives

\[
A(u)(\theta-\widetilde\theta)=0.
\]

The second and third conclusions follow immediately. \(\square\)

The second conclusion is an exact, global inverse statement for observation of
the complete state. It does not assert existence or uniqueness of a steady
state for each activity vector, and it does not establish global injectivity of
any partial observation map.

### Corollary 2.1. AE2 and AE4 have independent balance signatures

In the intracellular Cl, Na, K, and HCO3 balance rows, the two activity columns
of \(F_\theta\) are

\[
a_2=\varphi_2
\begin{pmatrix}
1\\0\\0\\-1
\end{pmatrix},
\qquad
a_4=\varphi_4
\begin{pmatrix}
1\\-\alpha\\-\beta\\-2
\end{pmatrix}.
\]

The Cl/HCO3 minor is

\[
\begin{pmatrix}
\varphi_2&\varphi_4\\
-\varphi_2&-2\varphi_4
\end{pmatrix},
\qquad
\det=-\varphi_2\varphi_4.
\]

Therefore, wherever

\[
\varphi_2\varphi_4\ne0,
\]

the activity forcing has rank two. If \(F_u\) is also invertible, the complete
steady-state map has regular rank two. Moreover, a fixed complete state cannot
be shared by two distinct pairs \((G_2,G_4)\).

This conclusion is independent of the disputed cell-volume sign and of the
AE4 Na/K partition, because the nonzero minor uses only the Cl and HCO3
stoichiometries. It is also invariant under any nonzero rescaling used to
repair the AE4 activity units. It fails on a transporter reversal surface
\(\varphi_j=0\), where the corresponding activity disappears from the
steady-state equations at that state.

Because \(D_\theta u^*\) is an \(n\times2\) matrix of rank two, at least one pair
of state coordinates has a nonzero \(2\times2\) sensitivity minor at each
regular active point. The proof is existential: it does not identify the pair,
and the successful pair can change with the state.

## Exact steady-state flux decomposition

Use the main-text signed particle-flux convention

\[
C=\frac{I_{\rm CaCC}}{Fz_{\rm Cl}},\qquad
N=J_{\rm Nkcc1},\qquad B=J_{\rm Buffer}.
\]

For outward chloride, \(I_{\rm CaCC}<0\) and \(z_{\rm Cl}=-1\), so \(C>0\).
This \(C\) must not be confused with the Appendix expression labelled
\(J_{\rm CaCC}=I_{\rm CaCC}/F\), which has the opposite sign for outward
chloride.

At steady state, the intracellular Cl and HCO3 balances give

\[
C=2N+J_2+J_4,
\qquad
B=J_2+2J_4.
\]

The H balance gives

\[
B=J_{\rm Nhe1}.
\]

The printed CO2 balance also gives \(B=J_{\rm CO2}\), but the published
direction of \(J_{\rm CO2}\) is unresolved. None of the decomposition below
requires that disputed equality.

Define the exchanger-mediated chloride aggregate

\[
A_{\rm tot}=C-2N=J_2+J_4.
\]

The two exchanger fluxes are then exactly

\[
J_4=B-A_{\rm tot}=B-C+2N,
\]

\[
J_2=2A_{\rm tot}-B=2C-4N-B.
\]

If the missing luminal subscript in the published Cl loss term is corrected,
the steady luminal Cl balance gives

\[
C=Q[\mathrm{Cl}]_l.
\]

Under that correction,

\[
J_4=B-Q[\mathrm{Cl}]_l+2N,
\qquad
J_2=2Q[\mathrm{Cl}]_l-4N-B.
\]

Thus \(Q\) and luminal Cl determine the total apical chloride flux \(C\), but
they do not by themselves decompose AE2 and AE4 because \(N\) and \(B\) remain
nuisance fluxes.

If \(\varphi_2\varphi_4\ne0\), the activities reconstruct as

\[
G_2=\frac{2A_{\rm tot}-B}{\varphi_2(u)},
\qquad
G_4=\frac{B-A_{\rm tot}}{\varphi_4(u)}.
\]

These are exact symbolic identities. Their numerical use requires a complete
state sufficient to evaluate the driving terms, all background fluxes in a
common unit, and a resolved AE4 scaling convention.

### Forward-flux wedge

If \(J_2\ge0\) and \(J_4\ge0\), then

\[
A_{\rm tot}\le B\le2A_{\rm tot}.
\]

For \(A_{\rm tot}>0\), the exchanger fractions are

\[
\frac{J_4}{J_2+J_4}=\frac{B}{A_{\rm tot}}-1,
\qquad
\frac{J_2}{J_2+J_4}=2-\frac{B}{A_{\rm tot}}.
\]

The boundary \(B/A_{\rm tot}=1\) is pure AE2 flux, and
\(B/A_{\rm tot}=2\) is pure AE4 flux. This wedge statement is not valid when
one of the net exchanger fluxes reverses.

### Independent cation-balance check

Let

\[
P=J_{\rm NaK},\qquad
K_{\rm out}=\frac{I_{\rm CaKC}}{Fz_K}.
\]

At steady state, the Na and K balances are

\[
0=N+B-\alpha J_4-3P,
\]

\[
0=N+2P-K_{\rm out}-\beta J_4.
\]

Adding them, using \(\alpha+\beta=1\), and substituting the Cl/HCO3 identities
gives

\[
C=P+K_{\rm out}.
\]

This is a useful exact regression check for any eventual implementation.

## A general stoichiometric rank criterion

### Theorem 3. Target fluxes modulo linear nuisance fluxes

Suppose observable or computable steady balance summaries satisfy

\[
r=S v+T w,
\]

where \(v\in\mathbb R^p\) are target transporter fluxes and \(w\) are
unrestricted nuisance fluxes. Treat \(S\) and \(T\) as fixed at the state and
condition under consideration. Then \(v\) is uniquely determined modulo the
nuisance fluxes if and only if

\[
\operatorname{rank}[T\ S]=\operatorname{rank}T+p.
\]

Equivalently, the induced map from target fluxes into the quotient balance
space \(\mathbb R^q/\operatorname{im}T\) is injective.

**Proof.** Two target values \(v\) and \(\widetilde v\) can produce the same
summary after changing nuisance fluxes precisely when

\[
S(v-\widetilde v)\in\operatorname{im}T.
\]

Thus uniqueness holds precisely when no nonzero vector in
\(\operatorname{im}S\) arising from a target-flux difference belongs to
\(\operatorname{im}T\). This is equivalent to the stated rank increment. \(\square\)

With no nuisance fluxes, the criterion reduces to
\(\operatorname{rank}S=p\). If

\[
v=\operatorname{diag}(\varphi_1,\ldots,\varphi_p)\theta
\]

and every \(\varphi_j\ne0\), flux identifiability transfers to activity
identifiability. If the quotient rank is \(r<p\) and the variables lie in the
interior of an unrestricted local domain, the linear equivalence family has
dimension \(p-r\).

This theorem is conditional on fixed linearized stoichiometric matrices and
unrestricted nuisance differences. Bounded flux domains, nonlinear nuisance
constraints, and state feedback can reduce an equivalence set. Conversely,
the theorem does not by itself establish identifiability from an arbitrary
partial state observation \(h(u)\), because such an observation need not make
the balance summary \(r\) or its nuisance fluxes known.

### Two stoichiometric signatures

For two forward transporters that both carry one unit of a reporter species
but carry \(k_1\ne k_2\) units of a second species, let

\[
A=v_1+v_2,
\qquad
D=k_1v_1+k_2v_2.
\]

Then

\[
v_1=\frac{k_2A-D}{k_2-k_1},
\qquad
v_2=\frac{D-k_1A}{k_2-k_1}.
\]

For nonnegative fluxes, \(D/A\) lies between \(k_1\) and \(k_2\). The AE2/AE4
decomposition is the case \((k_1,k_2)=(1,2)\).

## Minimal measurement criteria

Three levels of measurement must be kept separate.

1. **Number of scalar outputs.** At least two scalar outputs are necessary for
   regular local identification of two activities. Two outputs are not
   sufficient without a rank-two Jacobian.
2. **Balance-summary panel.** The pair \((A_{\rm tot},B)\) is exactly minimal
   for decomposing \((J_2,J_4)\), because its stoichiometric matrix
   \[
   \begin{pmatrix}1&1\\1&2\end{pmatrix}
   \]
   has determinant one. Either summary alone leaves a one-dimensional flux
   equivalence family.
3. **Physiological observation panel.** For a candidate panel such as
   \(h=(Q,[\mathrm{Cl}]_i)\), the required matrix is
   \[
   S_h=-h_uF_u^{-1}
   \begin{pmatrix}a_2&a_4\end{pmatrix},
   \]
   together with any direct \(h_\theta\) term. Stoichiometric rank does not
   determine the rank of this projected matrix.

The full-state rank result proves that some pair of state coordinates has a
nonzero minor at every regular active point. It neither names that pair nor
proves that it belongs to the biologically plausible candidate list. A
candidate-pair determinant must be computed only after a consistent baseline
and \(F_u\) have been reconstructed.

For numerical conditioning, derivatives should be taken with respect to
\((\log G_2,\log G_4)\) when the activities are positive, and observation rows
should be scaled by measurement noise or an explicit physical scale. Such
invertible scalings do not change rank but do change singular values and
condition numbers.

## Perturbation-induced distinguishability

Let \(H(\theta;\xi)\) be a scalar observation at condition \(\xi\), and define

\[
g(\xi)=D_\theta H(\theta_0;\xi).
\]

Measurements at conditions \(\xi_0\) and \(\xi_1\) have stacked sensitivity

\[
S_{\rm stack}=
\begin{pmatrix}
g(\xi_0)\\g(\xi_1)
\end{pmatrix}.
\]

They give regular first-order local identifiability if and only if

\[
\det S_{\rm stack}\ne0.
\]

Thus one extra condition breaks a regular baseline equivalence curve exactly
when its sensitivity vector is not parallel to the baseline sensitivity.
If \(\xi_1=\xi_0+\varepsilon d\), then

\[
\det[g(\xi_0);g(\xi_1)]
=\varepsilon\det[g_0;D_\xi g_0d]+O(\varepsilon^2).
\]

Consequently,

\[
\det[g_0;D_\xi g_0d]\ne0
\]

is a sufficient nondegenerate criterion for an arbitrarily small perturbation
in direction \(d\) to create first-order distinguishability. Its failure does
not rule out separation at second or higher order.

If the two scalar observations have independent Gaussian errors with variances
\(\sigma_0^2\) and \(\sigma_1^2\), the Fisher information satisfies

\[
\det\mathcal I
=\frac{\det(S_{\rm stack})^2}{\sigma_0^2\sigma_1^2}.
\]

Therefore, with one condition fixed, the natural two-condition D-optimal
design criterion is the noise-weighted area spanned by the two sensitivity
vectors. This analytical reduction does not show that any allowed salivary
perturbation produces a nonzero area; that is a model-specific calculation not
currently available.

## Global, discrimination, and singularity claims

### What is global

The complete-state common-state uniqueness in Proposition 2 is exact and
global in the activity variables at a fixed state and condition. The flux
decomposition from \((A_{\rm tot},B)\) is also exact and global on its stated
domain.

### What is not global

No chosen physiological partial observation map is proved globally injective.
Even a nonzero Jacobian determinant throughout a two-dimensional domain is not
by itself a complete global-univalence proof; an appropriate boundary,
properness, monotonicity, or global-univalence argument would also be needed.
Numerical sampling cannot supply that proof.

### Mechanism rays

After known background fluxes are removed, a pure AE2 mechanism gives the
balance-summary ray

\[
(A_{\rm tot},B)=v(1,1),
\]

whereas pure AE4 gives

\[
(A_{\rm tot},B)=v(1,2).
\]

The nonzero rays are distinct, but their closures meet at the origin. Hence a
positive minimum separation requires a nonzero lower flux bound or a justified
normalization. If the admissible flux can approach zero, the composite testing
distance is zero.

AE4 variants that differ only in Na/K partition have signatures

\[
(1,-\alpha,-\beta,-2).
\]

Their Cl/HCO3 projections are identical. No anion-only balance panel can
discriminate such variants; a cation-sensitive output, cation balance, or
cation perturbation is necessary.

### Inference-relevant singular sets

Three different singularities must not be conflated:

1. \(\varphi_2=0\) or \(\varphi_4=0\): an activity multiplier disappears from
   the balance equations at that state;
2. \(\det F_u=0\): the local steady branch fails the implicit-function
   hypothesis, as at a possible fold;
3. \(\det D_\theta H=0\) for a two-output panel: the chosen projection loses
   regular first-order rank even though the complete-state map may remain
   rank two.

Near the first type, practical information about the corresponding activity
vanishes. Near the second, sensitivities can become large or ill-conditioned,
which is not the same as reliable identifiability. The third type is
panel-specific. Continuation would be justified only after a consistent model
shows that one of these sets intersects the physiological domain.

## Certification boundary imposed by the published baseline

The analytical identities above are certified under their stated sign, unit,
and activity-linearity assumptions. The current published-equation audit does
not provide a consistent numerical steady state from which \(F_u\), a
physiological observation Jacobian, or perturbation sensitivities can be
computed. In particular, the audit recorded:

- a factor of approximately \(59.6\) between the two water fluxes that must be
  equal at steady cell volume;
- an apical quasi-steady current residual of approximately \(68\) pA;
- incompatible AE4 concentration and activity units, including a
  fourth-power M-versus-mM ambiguity of \(10^{12}\).

The machine-readable printed-equation diagnostics are in
`results/10_identifiability_discrimination/baseline_printed_equation_audit.json`.
They are an audit of the printed equations, not a reproduced model.

Accordingly, the following claims are **not certified**:

- full local rank of \(Q\) plus any named intracellular ion;
- a minimal biologically plausible state-observation panel;
- the image geometry of a physiological \((G_2,G_4)\) domain;
- perturbation-induced distinguishability for any bath, calcium, agonist, or
  inhibition condition;
- global injectivity or exact global equivalence of a partial panel;
- positive separation of full model-class image sets;
- an inference-relevant fold or bifurcation.

Any later numerical rank calculation must first exhibit a consistent baseline,
state the chosen resolutions of the published sign and unit conflicts, verify
the steady residuals, and check implicit derivatives against numerical
derivatives. Until then, the exact stoichiometric results are the strongest
defensible model-specific mathematics.
