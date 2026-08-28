# Dynamic beta/cAMP/PKA regulation of AE4

## Result

Task 13B now has four executable, nested regulatory families (`R0`--`R3`)
with a common API and physical seconds as the time coordinate. `R1` is the
smallest dynamic family. `R2` makes an effective cAMP coordinate explicit.
`R3` makes upstream effective PKA activity and a reversible AE4 regulated
fraction explicit. The 2021 evidence establishes the pathway direction and
construct-level response logic, but does not identify any acute kinetic
constant or uniquely select `R1`, `R2`, or `R3`.

All retained families regulate the modern SR2 AE4 module through one common
active-capacity multiplier. Transition-specific `R4` regulation is rejected:
the source does not identify which microscopic carrier edge would be altered.

Implementation: `src/modern_full_model/camp_pka.py`  
Focused tests: `tests/test_modern_camp_pka.py`

## What the 2021 primary evidence does and does not establish

The source statements used here are frozen in `evidence_ledger.md` and in the
pre-optimization target tables. They are not reinterpreted as kinetic data.

| Observation | Evidential use in this module | What is not inferred |
|---|---|---|
| Isoproterenol increases native exchanger/AE4 readout | beta input must increase AE4 capacity | receptor occupancy, an acute onset rate, or a universal fold |
| Forskolin increases heterologous AE4 readout | an adenylate-cyclase/cAMP-linked arm lies upstream | a measured cAMP concentration or cAMP time constant |
| H89 prevents the detected increase | qualitative PKA-dependent hierarchy | H89 specificity, a quantitative inhibition constant, or direct action on AE4 |
| Constitutively active PKAc increases both bicarbonate- and chloride-linked AE4 readouts | PKA activity is sufficient to increase the construct readout in that assay | a physiological PKA waveform or phosphorylation rate |
| S173A has no detected PKAc increment | S173 dependence of the construct-level response | direct phosphorylation of S173 |
| S273A retains the PKAc increment | S273 is not required for the detected increment | equality of WT and S273A kinetics or basal activity |
| Prose describes about a 25% forskolin increase, while the plotted means visually imply about a 1.7-fold increase | two declared gain sensitivities, not a point fit | pooling either number across native and heterologous contexts |
| The assay figures show depletion/reuptake protocols and endpoint initial slopes | direction and hierarchy gates | beta-to-cAMP, PKA, or AE4-regulatory onset timing |

There is no reported acute cAMP concentration trace, PKA activity trace,
AE4 regulated-fraction trace, activation/deactivation fit, quantitative
washout time constant, or direct demonstration that PKA phosphorylates S173.
Every kinetic default below is therefore a `NEW_MODELING_DECISION` used only
for sensitivity analysis.

## Common interface and capacity map

Each family object exposes

- `state_names`;
- `initial_state(beta_input)`, which returns the analytic constant-input
  equilibrium;
- `rhs(time_s, state, beta_input)`;
- `evaluate(time_s, state, beta_input)`.

The evaluation contains `family`, normalized `activity`, and the dimensionless
`capacity_multiplier` consumed by the full model. Let

- \(u\in[0,1]\) be normalized effective beta input;
- \(a\in[0,1]\) be a family's AE4 regulatory activity;
- \(M_b\) be basal AE4 capacity;
- \(\Delta M\) be the fully activated increment;
- \(s_c\) be a declared coupling sensitivity;
- \(s_m\) be the construct-level response switch.

The only retained coupling is

\[
M_{AE4}=M_b+\Delta M\,s_c\,s_m\,a.
\]

In `transporters.py`, this multiplier scales active carrier amount before both
Na and K branch fluxes are evaluated. It does not change branch stoichiometry,
affinity, reversal, detailed balance, or one selected transition rate.

The native WT construct and S273A use \(s_m=1\); S173A uses \(s_m=0\) for the
PKA-dependent increment. This switch deliberately says “response depends on
S173,” not “S173 is directly phosphorylated.” Approximate basal mutant
differences are not converted into kinetic parameters.

The default \(M_b=1\), \(\Delta M=0.25\), and \(s_c=1\) records the 2021 prose
reading as one sensitivity. A second registry entry uses
\(\Delta M=0.70\), corresponding to the approximate visual figure reading.
Both are tagged `SENSITIVITY_ONLY_NOT_POINT_TARGET`; neither is a universal
capacity measurement. Setting \(s_c=0\) is a qualitative PKA-dependence
control, not a quantitative or specific model of H89.

## Nested families

### R0 — instantaneous static control

`R0` has no state:

\[
a(t)=u(t).
\]

It is a static source-supported fold comparison and a negative control for
the mandatory dynamic investigation. It is not an adequate final substitute
for beta/cAMP/PKA dynamics.

### R1 — one effective beta/PKA/AE4 state

With effective activation \(x\),

\[
\frac{dx}{dt}=\frac{u-x}{\tau_x},\qquad a=x.
\]

The default \(\tau_x=30\) s is an unmeasured sensitivity. This state aggregates
adenylate-cyclase/cAMP, PKA, and the AE4 response. For constant \(u\),

\[
x^*=u.
\]

At \(x=0\), \(dx/dt\ge0\); at \(x=1\), \(dx/dt\le0\). Thus `[0,1]` is forward
invariant for every admissible input.

### R2 — effective cAMP plus downstream activation

Let \(c\) be a normalized effective cAMP coordinate and \(x\) the downstream
AE4 activation coordinate:

\[
\frac{dc}{dt}=\frac{u-c}{\tau_c},\qquad
\frac{dx}{dt}=\frac{c-x}{\tau_x},\qquad a=x.
\]

The defaults \(\tau_c=10\) s and \(\tau_x=30\) s are unmeasured sensitivities.
The analytic constant-input equilibrium is

\[
(c^*,x^*)=(u,u).
\]

Both state faces point inward, so the unit square is forward invariant. Under
constant input, the manifold \(c=u\) is invariant. On that manifold the
second R2 equation is exactly R1, not an approximate numerical resemblance.

The coordinate \(c\) is not a cAMP concentration and must not be compared to
a biochemical cAMP assay without an observation map.

### R3 — effective PKA plus reversible AE4 regulated fraction

Let \(p\) be effective PKA activity and \(q\) an AE4 regulated fraction:

\[
\frac{dp}{dt}=\frac{u-p}{\tau_p},
\]

\[
\frac{dq}{dt}=k_f p(1-q)-k_r(1-p)q,\qquad a=q.
\]

The coordinate \(q\) is phosphorylation/dephosphorylation-*like*, but the
model and code use the neutral names `forward_regulation_rate_s` and
`reverse_regulation_rate_s`. This prevents an unmeasured direct S173
phosphorylation step from becoming an apparent fact.

For constant input,

\[
p^*=u,\qquad
q^*=\frac{k_f u}{k_f u+k_r(1-u)}.
\]

At \(q=0\), its derivative is \(k_fp\ge0\); at \(q=1\), its derivative is
\(-k_r(1-p)\le0\). Together with the inward PKA equation, the unit square is
forward invariant.

The defaults are \(\tau_p=10\) s and \(k_f=k_r=1/30\) s\(^{-1}\). These
choices are not estimates. They make

\[
\frac{dq}{dt}=\frac{p-q}{30\ \mathrm{s}},
\]

so default R3 is exactly dynamically isomorphic to default R2 under
\(p\leftrightarrow c\) and \(q\leftrightarrow x\), for any beta waveform.
This exact collapse is useful: endpoint evidence cannot distinguish two
differently named but equation-equivalent pathway interpretations.

## Parameter and dimensional ledger

| Family | Quantity | Unit | Default | Status |
|---|---|---:|---:|---|
| all | beta input \(u\) | dimensionless fraction | protocol-dependent | normalized model input, not receptor occupancy measurement |
| all | activity and regulatory states | dimensionless fraction | analytic basal equilibrium | derived bounded coordinates |
| all | basal capacity \(M_b\) | dimensionless multiplier | 1.0 | normalization decision |
| all | fully activated increment \(\Delta M\) | dimensionless multiplier | 0.25; 0.70 sensitivity | conflicting approximate 2021 readings; not fitted |
| all | coupling scale \(s_c\) | dimensionless | 1.0 | modeling switch; zero only for qualitative dependence check |
| R1 | `tau_activation_s` | s | 30 | unmeasured kinetic sensitivity |
| R2 | `tau_camp_s` | s | 10 | unmeasured kinetic sensitivity |
| R2 | `tau_activation_s` | s | 30 | unmeasured kinetic sensitivity |
| R3 | `tau_pka_s` | s | 10 | unmeasured kinetic sensitivity |
| R3 | `forward_regulation_rate_s` | s\(^{-1}\) | 1/30 | unmeasured kinetic sensitivity |
| R3 | `reverse_regulation_rate_s` | s\(^{-1}\) | 1/30 | unmeasured kinetic sensitivity |

No hidden conversion from historical code time is used. `rhs` returns
fraction per second and the full stiff integrators receive seconds directly.
The analytic vector fields preserve the unit interval exactly. Radau/BDF can
nevertheless query an internal, unaccepted trial stage just outside that
interval. State evaluation therefore clamps only a declared numerical trial
band of `1e-4`; larger excursions fail. This band is not used for beta input,
parameters, or analytic equilibrium validation, which retain the exact
`[0,1]` domain. A short integrated full-model Radau step is a regression test
for this separation.

## Family discrimination and R4 decision

The available 2021 comparisons constrain a monotone beta/cAMP/PKA-dependent
capacity increase and the S173/S273 response hierarchy. They provide endpoint
or pretreated assay contrasts, not a resolved response to an acute beta step.
Consequently:

1. `R0` is the exact static equilibrium map but cannot represent timing.
2. `R1` is the parsimonious dynamic family.
3. `R2` adds an unobserved effective cAMP state; its two time constants cannot
   be separated by the existing endpoint evidence.
4. `R3` adds a mechanistic interpretation of the regulated fraction; its
   forward and reverse rates are likewise unconstrained, and its symmetric
   form is exactly equivalent to `R2`.
5. `R4`, modulation of selected AE4 transition rates, is
   `REJECTED_NOT_RETAINED`. No transporter experiment in the audited 2021
   source identifies a transition-specific regulatory edge, so retaining R4
   would add flexibility without an independent discriminator.

An acute native time course measuring beta stimulation together with either
cAMP, PKA activity, or an AE4 regulatory-state readout would be needed to
separate upstream and downstream kinetic scales. Until such a measurement is
available, family choice and time constants must remain declared WT-only
sensitivity decisions and may not be selected from an unrelated phenotype
trajectory.

## Verification

The focused suite establishes:

- exact analytic equilibrium for R0--R3;
- inward vector fields and bounded numerical step responses;
- exact R0/R1 equilibrium nesting;
- exact R2-to-R1 collapse on the constant-input cAMP manifold;
- exact R2/R3 dynamic isomorphism under the symmetric default rates;
- the S173-dependent/S273-retained qualitative response logic;
- preservation of both 2021 gain readings as sensitivities;
- direct compatibility with `ModernFullModel` without an adapter;
- an explicit, non-executable R4 rejection.

Focused command:

```bash
PYTHONPATH=src python -m unittest tests.test_modern_camp_pka -v
```

Result: 20 tests passed.
