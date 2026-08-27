# Structural explanation of the Task-13 result

## Controlling conclusion

Task 13 separates three statements that are easy to conflate.

1. Explicit shared AE4 states add real transporter structure: finite-pool
   Na/K competition, multicycle coupling, and a possible Na/K slip current.
2. In the implemented whole-cell model those carrier occupancies are solved
   quasi-steadily, so they enter the seven-state chassis only through
   algebraic functions `J_Na(x,PKA)` and `J_K(x,PKA)`. They add no carrier
   memory and cannot by themselves generate a delayed whole-cell response.
3. Deleting AE4 annihilates both functions exactly. Every implemented AE4
   family then has the same fixed-chassis knockout vector field, whose bounded
   root misses the released knockout Cl/pH measurements.

Thus the investigation does not establish that AE4 biology is dispensable or
that a molecular cycle is irrelevant. It establishes that state resolution,
as embedded quasi-steadily in this unchanged seven-state chassis, cannot by
itself repair the knockout balance. Stage B then proves that the two released
ionic measurements do not uniquely identify the missing chassis module.

The numerical chain is recorded in
[`transporter_fit_summary.csv`](../../results/13_state_resolved_ae4/transporter_fit_summary.csv),
[`stage_a_summary.json`](../../results/13_state_resolved_ae4/stage_a_summary.json),
[`stage_a_roots.csv`](../../results/13_state_resolved_ae4/stage_a_roots.csv),
[`stage_a_null_roots.csv`](../../results/13_state_resolved_ae4/stage_a_null_roots.csv),
and
[`module_projection.csv`](../../results/13_state_resolved_ae4/module_projection.csv).

## Claim classes

| Claim | Class and scope |
| --- | --- |
| Carrier conservation, `J_Cl=J_Na+J_K`, the branch source vector, local detailed balance, and the Na/K slip identity | **Necessary exact consequences** of the declared passive shared-state graph. |
| The implemented carrier contributes no new whole-cell state or delay | **Necessary exact consequence** of solving CTMC occupancy QSS at every RHS evaluation. It is not a statement about a future model with dynamic carrier occupancies. |
| All AE4 candidates have the same knockout field | **Necessary exact consequence** of the genotype-zero short circuit while every non-AE4 equation and parameter is fixed. |
| Zero valid WT capacity calibrations and zero WT Gate-2 passes open Stage B | **Sufficient for the frozen question-tree transition**, but not sufficient to identify the missing module or reject AE4 biology. |
| Every transporter-passing SR5/SR6-core point in the saved 17-point gauge grid has the wrong chloride-source sign | **Finite-domain parameterization result**, not a global theorem. The wider adversarial optimizer search is also finite-domain. |
| Released knockout Cl/pH define a five-dimensional target-state manifold | **Exact identifiability result** for the retained seven-state/electroneutrality reduction. |
| The associated required-source family has local rank five | **Numerical local result**, reproduced at the mean and four one-SEM corners. It rules out a unique missing vector at this evidence level; it is not a global rank theorem. |
| Apical/basal pump and K redistribution add rank two in expanded lumen/current space and rank zero in intracellular chemical projection | **Exact structural result** for the declared per-event signatures. It proves that the old projection could not test the split, not that the split is sufficient to repair the phenotype. |
| Conserved carbon plus NHE1/buffer flux is the needed repair | **Not established.** Cl/pH directly resolve only rank one of three acid-base reaction directions. |

## 1. What shared carrier states change

For the minimal shared graph, Na and K cycles use one conserved carrier pool.
In the reduced two-conformation notation, the Cl edge has rates `a,b`, the Na
return edge has `c_Na,d_Na`, and the K return edge has `c_K,d_K`. With

\[
C=c_{Na}+c_K,\qquad D=d_{Na}+d_K,\qquad
\Sigma=a+b+C+D,
\]

the exact QSS currents are

\[
J_s=\frac{c_s(a+D)-d_s(b+C)}{\Sigma},\qquad
J_{Cl}=\frac{aC-bD}{\Sigma}=J_{Na}+J_K. \tag{1}
\]

Equation (1) contains the structural difference from Task 12. Both branch
currents depend on the same denominator and on the other branch's rates.
Changing K can therefore change Na turnover without inserting an arbitrary
cross-gate. One branch may even run against its isolated affinity while the
coupled network still satisfies

\[
J_{Na}\mathcal A_{Na}+J_K\mathcal A_K\ge 0. \tag{2}
\]

The intracellular source for the `1:1:2` member is exactly

\[
S_{AE4}=(-J_{Na},-J_K,J_{Na}+J_K,-2(J_{Na}+J_K))
\tag{3}
\]

in `(Na,K,Cl,HCO3)` order. If `J_Na=-J_K`, the anion entries vanish but the
cation entries do not: the network carries a pure Na/K exchange slip. A
single scalar `J_AE4` cannot express that circulation.

These facts do not mean that the current experiments identify shared states.
In a pure-Na or pure-K bath the unused branch disappears, and in the
rare-occupancy limit the shared denominator is effectively constant. The
model then reduces to the branch-resolved Task-12 scalar law. More generally,
at any one whole-cell state, two scalar functions chosen to return the same
pair `(J_Na,J_K)` produce the same source (3). Shared states become
observationally distinct only across mixed-cation competition/reversal
protocols or through measured carrier-state transients. The available
pure-cation dose summaries constrain effective loaded-state energies and
attempt-rate combinations, not microscopic binding order or individual edge
rates. See [`state_cycle_models.md`](state_cycle_models.md) and
[`thermodynamic_derivation.md`](thermodynamic_derivation.md).

## 2. Why explicit states did not add delayed propagation

The implemented CTMC stationary distribution is recalculated at every
whole-cell RHS call. Its occupancies are diagnostics, not ODE states. After
elimination, every candidate is exactly a static map

\[
x\longmapsto (J_{Na}(x,u_{PKA}),J_K(x,u_{PKA}))
\longmapsto S_{AE4}(x,u_{PKA}). \tag{4}
\]

Consequently the whole-cell state remains seven-dimensional. Explicit state
competition can change the concentration dependence and the instantaneous
Jacobian, but it supplies no relaxation eigenmode of its own. SR6 likewise
uses explicit U/P layers in the stationary graph, yet without a calibrated
phosphorylation-state ODE or activation time constant it cannot support a
physical two-to-three-minute timing claim. This exact reduction is audited as
`T08` in
[`adversarial_checks.csv`](../../results/13_state_resolved_ae4/adversarial_checks.csv)
and is independently inspectable from the frozen Q matrices in
[`ctmc_snapshot_manifest.json`](../../results/13_state_resolved_ae4/ctmc_snapshot_manifest.json).

This does not prove that carrier memory is physiologically absent. It says
only that the present implementation did not test it. A dynamic carrier model
would also have to account for transient bound-ion storage and obtain its time
scale from an independent transporter protocol before it could be credited
with gland-level delay.

## 3. Transporter compatibility versus whole-cell direction

Thirty-three transporter/gauge rows were retained, of which five conditionally
express the permitted dose-response, electroneutrality, mutation-limit, and,
where claimed, PKA-fold constraints. They are conditional observation-map
fits rather than identified mechanisms. Their fitted quantities include two
loaded-state energy combinations, an Na/K attempt-rate ratio, a common assay
scale, two fixed `Rmin` nuisances, a discrete barrier gauge, and for SR6 a
regulated-layer factor. The mutant hierarchy is imposed as an expressivity
limit rather than predicted out of sample.

The salivary insertion exposes a separate issue. In the saved KO-blind
17-point log-spaced SR5/SR6-core barrier grid from `10^-4` to `10^4`, all 11
points that pass the transporter summary have negative unit-capacity chloride
source (`-0.0257267` to `-0.00443125`), whereas the fixed WT row demands
positive chloride loading. At accepted broad-gauge points, the source ratios
approach

\[
Na/Cl\simeq-1.0205,\qquad K/Cl\simeq+0.0205,\qquad HCO_3/Cl=-2. \tag{5}
\]

The sequential fit at gauge `0.1` is optimizer-order sensitive, but a joint
fit can match the effective dose summaries and still gives chloride source
`-0.0111649851` per capacity. A wider finite-domain differential-evolution
challenge found no gate-feasible positive source; its best value was
`-1.75670e-5` at an energy bound. These are strong finite-domain exclusions,
not an exhaustive proof over all microscopic parameterizations. The complete
grid is
[`barrier_gauge_source_scan.csv`](../../results/13_state_resolved_ae4/barrier_gauge_source_scan.csv),
and the scope limitation is recorded as checks `T05-T07` in
[`adversarial_checks.csv`](../../results/13_state_resolved_ae4/adversarial_checks.csv).

The Stage-A conclusion does not depend on promoting that finite grid to an
invariant. All five conditional rows failed the fixed WT chassis independently:
every capacity attempt was invalid, and the exact bounded roots obtained after
freezing those attempts missed WT pH by about `6.69` reported SEM. Hence there
is no admissible frozen WT denominator for a secretion prediction. Details are
in [`transporter_fit.md`](transporter_fit.md) and
[`stage_a_results.md`](stage_a_results.md).

## 4. Why the knockout field is candidate independent

Let `F_0(x)` denote all fixed non-AE4 balances and let candidate `j` contribute
`S_j J_j(x)`. The implementation has the exact genotype-zero short circuit

\[
F_j^{KO}(x)=F_0(x)+0\,S_jJ_j(x)=F_0(x). \tag{6}
\]

Equation (6) is independent of state graph, Na/K barrier gauge, PKA placement,
and fitted AE4 capacity. It is stronger than observing numerically similar
roots: the vector fields are identical. Therefore changing only the AE4 cycle
cannot change any knockout root of this fixed non-AE4 chassis.

Within the documented positive domain with cell height between `2` and
`200 um`, the common 20-start root is

\[
Cl_i=48.78221503\ {\mathrm{mM}},\qquad pH_i=7.37847391,
\]

with raw residual `5.77e-14`. It misses the released targets by `+7.676` and
`+24.424` reported SEM. A remote exact branch exists near height
`1.25e14 um`; it is physiologically inadmissible but prevents any global
uniqueness claim. Thus the exact statement is candidate independence of the
field, plus failure of every root found in the declared physiological domain,
not uniqueness over all positive coordinates.

Candidate independence is sufficient to show why another QSS AE4 graph alone
cannot repair the fixed-chassis knockout state. It is not sufficient to say
which surrounding transport or acid-base process is missing.

## 5. Why Stage B yields a family rather than one residual vector

The retained state is

\[
x=(Na_l,K_l,H,Na_i,K_i,Cl_i,HCO_{3i}).
\]

Only knockout `Cl_i` and pH are released. With fixed impermeant charge `X` and
intracellular proton eliminated by electroneutrality,

\[
h_i=HCO_{3i}+Cl_i+X/H-Na_i-K_i. \tag{7}
\]

At `Cl_i=36.50 mM` and `pH_i=6.89`, Eq. (7) fixes

\[
HCO_{3i}=h_i^*-36.50-X/H+Na_i+K_i, \tag{8}
\]

but leaves `(Na_l,K_l,H,Na_i,K_i)` free. Positivity restricts this exact
five-dimensional target manifold without reducing its dimension. No knockout
volume, Na, K, or lumen composition measurement exists to select one member.

For each admissible member `x`, the required missing source is
`d(x)=-F_0(x)` in the expanded amount/current/water coordinates. Its
finite-difference derivative with respect to the five free coordinates has
five nonzero singular values

\[
(0.40223543,0.01887636,0.00934814,0.00163360,0.000449923), \tag{9}
\]

so the local required-source family has rank five. The rank remains five at
all four one-SEM Cl/pH corners. The source therefore rotates across states
that are observationally indistinguishable under the released evidence; the
ambiguity is not one unknown magnitude multiplying a fixed vector.

This explains why one- and two-module rankings change between the
`model_anchored` and `wt_reference_cations` slices and why no one/two cell
signature is exact on both. It also explains why the local `(Cl,pH)` response
matrix can have rank two while remaining nonidentifying: its 11 parameter
columns have a nine-dimensional nullspace, so many parameter pairs can fit two
outputs algebraically, often only through changes outside local validity.
Neither a best cosine nor a two-output solve is a reconstructed mechanism.
The full derivation is in [`missing_direction.md`](missing_direction.md).

## 6. What the membrane split proves—and does not prove

An apical and a basolateral Na/K-pump cycle have the same intracellular
signature `(-3,+2,0,0,0,0)`. Apical and basolateral K efflux likewise share
the intracellular signature `(0,-1,0,0,0,0)`. Their differences appear only
after lumen and membrane-current rows are retained:

\[
\delta\sigma_P=(0_i;+3,-2,0,0;+1,-1,0), \tag{10}
\]

\[
\delta\sigma_K=(0_i;0,+1,0,0;+1,-1,0). \tag{11}
\]

Equations (10)-(11) have exact rank two in expanded space and rank zero in the
six intracellular chemical rows. This is why Task 12's intracellular source
projection could not test membrane redistribution. It is also why the
released intracellular Cl/pH pair cannot identify the apical pump fraction
or apical K fraction.

The rank result does not make the split inert. After re-solving lumen Na/K,
the two voltage closures, tight-junction currents, channel driving forces, and
water flux, rank-zero direct cell signatures can have indirect equilibrium
effects. The nested Round-5 implementation conserves total pump and K capacity,
uses local luminal K for apical pump turnover, inserts explicit luminal terms,
and reduces algebraically to the fixed chassis when both apical fractions are
zero. Those are necessary implementation checks. They are not evidence that
particular fractions are calibrated in mouse submandibular acini or that M1 is
sufficient to repair Gates 2-4. The exact signatures and ranks are exported in
[`module_projection.csv`](../../results/13_state_resolved_ae4/module_projection.csv),
and the source boundary is documented in
[`chassis_source_audit.md`](chassis_source_audit.md).

## 7. Acid-base identifiability

The minimal explicit reaction block contains three independent event
directions: NHE1 (`Na` entry/`H` exit), carbon hydration
(`CO2 -> HCO3 + H`), and buffer dissociation (`BH -> B- + H`). Their full
signature matrix has exact rank three. Projecting them onto the two released
observations `(Cl_i,H_i)` gives a zero Cl row and one H row, hence rank one and
reaction-flux nullity two. The membrane split has direct acid-base rank zero.

Therefore the same resting pH can arise from distinct NHE1, hydration, buffer,
and boundary-flux combinations. Adding unknown total carbon, buffer pool,
luminal carbon, or CO2 exchange only enlarges the ambiguity. A conserved
carbon/buffer module is structurally preferable to the historical algebraic
CO2 formula, but Cl/pH alone cannot calibrate its rates or assign a
genotype-specific change. The smallest discriminating data are total
inorganic carbon and buffer-pool constraints or a matched NHE/carbon flux
measurement, alongside the currently missing native Na/K/volume coordinates.

## 8. Reconstruction boundary and next discriminator

The current substantive conclusion is unresolved by proof, not by a stopped
optimizer:

- state-resolved transporter families were carried through the fixed-chassis
  failure rather than stopping after one family;
- the knockout-field identity shows why more AE4-only QSS parameterization
  cannot move the fixed-chassis knockout field;
- exact target-manifold geometry and rank-five source variation show why the
  released Cl/pH pair cannot select one correction;
- expanded signatures establish the apical pump/K split as a genuinely new,
  source-supported topology, while also proving that its two fractions are
  invisible to the released intracellular projection; and
- acid-base rank leaves at least two reaction-flux directions unresolved even
  before unknown pools and boundary fluxes are added.

Accordingly, no further unconstrained inference scan over pump/K fractions or
acid-base gains can establish a mechanism. Such a scan would choose among
unmeasured coordinates using the same two released outputs and would convert
structural nonidentifiability into a numerical optimum. The admissible next
test is a predeclared sensitivity or a fit to independent WT/native
measurements, followed by Gates 2-4 with all secretion outcomes still sealed.

The decisive missing measurement is the matched-gland, stimulus-dependent
tuple `(P_a,P_b,J_Ka,J_Kb)` together with intracellular Na, K, and cell volume,
luminal or collected-effluent Na and K, and retained Cl and pH. Ideally the
same protocol also measures total inorganic carbon plus an NHE- or buffer-flux
readout. If that extended cation-current support predicts the ionic
displacement under a frozen carbon model, it tests M1. If the tuple is normal
yet the Cl/pH displacement remains unexplained, the residual moves
specifically to conserved carbon/NHE/buffer flux. Without the carbon arm this
is the first decisive cation-topology branch test, not a complete mechanism
identifier. Until one of those branches is independently constrained, there
is no licensed WT dynamics or held-out secretion validation: no candidate has
passed WT Gate 2, and the historical physical time/flow scales remain
uncertified.
