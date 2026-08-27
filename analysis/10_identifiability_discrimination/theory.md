# Identifiability and equivalence theory

## Setting

Let `F: R^n x R^2 -> R^n` and `h: R^n -> R^m` be continuously differentiable. Write `theta=(G2,G4)`. Suppose

`F(u0;theta0)=0` and `det F_u(u0;theta0) != 0`.

The implicit function theorem gives a unique local steady branch `u*(theta)` and

`D_theta u* = -F_u^{-1} F_theta`.

For `y(theta)=h(u*(theta))`, the observation Jacobian is

`J_h(theta0) = -h_u F_u^{-1} F_theta`.

If an observable depends directly on `theta`, its direct derivative must be added. None is assumed here.

## Proposition 1: one regular scalar observable leaves a curve

For `m=1`, `rank J_h <= 1`, so a scalar cannot supply two independent local parameter directions. More precisely, if `grad_theta y(theta0) != 0`, the regular-level-set theorem implies that

`{theta : y(theta)=y(theta0)}`

is locally a one-dimensional `C^1` curve. Its tangent vectors `v` satisfy

`J_h v = 0`.

If the gradient vanishes, the scalar has no first-order identifying information; higher-order behavior may create an isolated nongeneric point, but it does not support a generic two-parameter identification claim. Independently of derivatives, a continuous scalar map cannot be one-to-one on an open two-dimensional parameter region. Therefore secretion rate alone is insufficient for generic identification of independent AE2 and AE4 activities.

## Proposition 2: local identification by a measurement panel

If `rank J_h(theta0)=2`, two components of `h` have a nonsingular `2 x 2` parameter Jacobian. The inverse function theorem applied to those components gives a locally unique inverse for `theta`. Hence rank two is sufficient for local identifiability.

If `rank J_h=r<2` is constant near `theta0`, the constant-rank theorem gives a local observational-equivalence manifold of dimension `2-r`.

This is a local result. It does not imply global injectivity, and a small nonzero singular value is a practical-conditioning warning rather than a failure of structural rank.

## Proposition 3: projected transporter-signature criterion

For the independently recovered reduction

`F(u;theta)=f0(u)+s2 G2 phi2(u)+s4 G4 phi4(u)`,

the parameter derivative at a fixed regular state has columns

`F_theta = [s2 phi2, s4 phi4]`.

Define the panel projection `P=-h_u F_u^{-1}`. If both activity factors are nonzero, then

`rank J_h = rank [P s2, P s4]`.

Thus a panel identifies `(G2,G4)` locally exactly when the two transporter signatures remain linearly independent after propagation through the steady-state restoring operator and projection into measured outputs.

Consequences:

1. Distinct stoichiometric columns are not by themselves enough; `F_u^{-1}` and the measurement projection matter.
2. Two named concentrations are not guaranteed to identify the activities until their two projected columns are checked.
3. Structural deficiency occurs when a transporter is inactive (`phi_j=0`), when the regular branch fails, or when the panel collapses the projected columns onto the same line.

## Proposition 4: consistent residual-row signs do not affect sensitivity

Let `D` be any constant invertible diagonal matrix, including a matrix of equation-sign choices, and define `F_tilde=D F`. Then

`-F_tilde_u^{-1} F_tilde_theta = -(D F_u)^{-1}D F_theta = -F_u^{-1}F_theta`.

Therefore a complete, consistent equation-row sign change cannot affect identifiability. A physical redefinition of flux direction or of an observable is different and must be audited separately.

## Two conditions and perturbation-induced rank

For experimental conditions `c=1,...,k` with condition-specific regular states and scalar outputs, stack the rows

`J_stack = [J_1; ...; J_k]`.

Two activities are locally identifiable from the combined experiments if and only if `rank J_stack=2`. For two scalar conditions this is equivalent to

`det([J_1; J_2]) != 0`.

This gives the exact local form of “one condition insufficient, two conditions sufficient”: the perturbation must rotate the sensitivity row rather than merely rescale it. The statement is general experimental-design linear algebra; without a closed salivary model it cannot identify which calcium, bath, or channel perturbation produces the required rotation.

## Exact results for the affine diagnostic only

For `F_d=A(x-1)-S(theta-1)`, every observation panel is affine:

`y(theta)=y0+H A^{-1}S(theta-1)`.

Therefore:

- a scalar panel with nonzero gradient has exact straight-line equivalence sets;
- a two-output panel is globally injective on every domain if its `2 x 2` Jacobian is nonsingular; and
- the observation image is a parallelogram for a rectangular activity domain.

These global statements are validation fixtures for the code. They are not global claims about the nonlinear published salivary model.

## What has and has not been proved

Proved for any regular steady-state model: the implicit derivative, scalar regular equivalence curve, rank-two local criterion, projected-signature criterion, consistent-row-sign invariance, and stacked-condition rank criterion.

Not proved: regularity of the 2018 model at a stimulated steady state, numerical rank of any real measurement panel, global physiological injectivity, separation of AE2/AE4 mechanism classes, or an inference result under measurement noise.
