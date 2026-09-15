"""Typeset the completed, independently verified Task 41 inverse analysis."""
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent
D=HERE/'inverse_detail'
read=lambda name:json.loads((D/name).read_text())
check=read('independent_verification.json')
assert check['status']=='PASS'
cross=read('crossings.json');selected=read('selected_reproduction.json')
long=read('long_time_summary.json');sens=read('threshold_sensitivities.json')
a=cross['adaptive'];s=cross['sampled']
names={'NBC_capacity':'NBC capacity','AE4_carrier':'AE4 carrier amount','NKCC_scale':'NKCC effective scale',
       'NHE_carrier':'NHE1 carrier amount','pump_capacity':'Pump capacity','CaCC_conductance':'CaCC conductance',
       'buffer_pool':'Cell buffer amount','water_apical':'Apical water permeability'}
rows=[]
for name,title in [('wt','WT'),('ae4_5pct',r'5\% AE4'),('ae4_null','AE4 null')]:
    r=selected[name]
    deficit='0' if name=='wt' else f"{100*r['sampled_deficit']:.6f}"
    rows.append(f"{title} & {r['cumulative_sampled_pL']:.9f} & {r['cumulative_adaptive_pL']:.9f} & {deficit} \\\\")
sr=[]
for r in sens:
    if r['observable']=='adaptive':
        sr.append(f"{names[r['parameter']]} & {r['crossing_log_elasticity']:+.7f} & {r['residual_recruitment_log_elasticity']:+.7f} \\\\")
text=r'''\section{The prescribed Task 41 inverse family}
\label{sec:task41inverse}

The Task 41 construction is target selected. Its added hypothesis changes
stimulated CaCC conductance before the full electrical and NBC closure:
\begin{equation}
 g_{\mathrm{CaCC}}(e,a;\rho)
 =g_{\mathrm{parent}}\{1-a\rho(1-e)\},\qquad
 \rho=1-b,\quad 0\leq\rho<1.
 \label{eq:inversefamily}
\end{equation}
Here $e$ is AE4 expression and $a$ is the protocol normalised recruitment
fraction used by NBC:
\[
 a=\min\!\left\{1,\max\!\left\{0,
 \frac{[\mathrm{Ca}]_i-0.058\,\mu\mathrm M}{(0.25-0.058)\,\mu\mathrm M}
 \right\}\right\}.
\]
The parent CaCC Hill gate remains a distinct retained factor.
The full conserved system, carbonate closure, electrical closure, both
volumes and AE4 regulation remain active. WT nests the parent exactly for
every state and every $\rho$. At full activation, the null retains fraction
$b=1-\rho$ of its parent CaCC conductance.

Every 600 s trajectory uses the identical frozen WT initial conserved state
and original stimulus protocol. No genotype resting solve, new mechanism or
parameter fit is introduced. All inherited positivity, concentration, pH,
volume, source and conservation gates are checked at accepted solver endpoints
and at the original observation times. The original WT stimulation requirement
is also retained; it is not imposed retrospectively on mutants.

For $X\in\{\mathrm{A},\mathrm{S}\}$ define
\begin{equation}
 D_X(\rho,p)=1-\frac{I_{0,X}(\rho,p)}{I_{1,X}(p)}.
 \label{eq:inverseobservable}
\end{equation}
$I_{e,\mathrm A}$ is the adaptive time integral of outflow and
$I_{e,\mathrm S}$ is the original trapezoidal observable on integer seconds,
including $t=0$ and $10^{-6}$ s. These are distinct numerical observations.
The prescribed comparison threshold is $D_X=0.303$: it is one reported
standard error below the experimental mean of 35\%, whose reported standard
error is 4.7 percentage points. It is not a confidence limit, a measured
biological range or an admissibility gate.

\subsection{Reproduction and numerical crossing}
The selected Task 41 parameter is
$b=0.10511872843288446$, or $\rho=0.8948812715671155$.
The three selected trajectories reproduce the frozen sampled cumulative
outputs to within '''+f"{max(abs(v['reproduction_sampled_error_pL']) for v in selected.values()):.2g}"+r''' pL.
\begin{center}
\begin{tabular}{lrrr}
\hline
Case & Sampled output (pL) & Adaptive output (pL) & Sampled deficit (\%) \\
\hline
'''+ '\n'.join(rows)+r'''
\hline
\end{tabular}
\end{center}
In particular, the frozen selected null deficit is slightly below 30.3\%.
It must not be described as exactly attaining this new comparison threshold.

An ordered scan of 17 recruitment fractions, followed by bracketed numerical
root refinement in the prescribed family, gives:
\begin{center}
\begin{tabular}{lrrr}
\hline
Observable & Crossing $\rho$ & AE4 dependent recruitment (\%) & Residual $b$ \\
\hline
Adaptive & '''+f"{a['rho']:.10f} & {100*a['rho']:.7f} & {a['b']:.10f}"+r''' \\
Original sampled & '''+f"{s['rho']:.10f} & {100*s['rho']:.7f} & {s['b']:.10f}"+r''' \\
\hline
\end{tabular}
\end{center}
Both crossings pass every inherited 600 s gate. At the adaptive crossing,
the original sampled deficit is '''+f"{100*a['check']['sampled_deficit']:.7f}"+r'''\%;
at the sampled crossing, the adaptive deficit is
'''+f"{100*s['check']['adaptive_deficit']:.7f}"+r'''\%.
Refinement at relative solver tolerances $10^{-8}$ and $2\times10^{-9}$,
with respective maximum steps 2 s and 1 s, changes the two crossing fractions
by '''+f"{abs(a['solver_tolerance_rho_difference']):.2g}"+r''' and
'''+f"{abs(s['solver_tolerance_rho_difference']):.2g}"+r''', respectively.
Independent integration of the 13 production states, followed by adaptive
quadrature of dense flow, also reproduces the selected cases and crossing.

The deficit increases at the sampled points, but monotonicity between them
and exclusion of all earlier unsampled crossings have not been proved.
These are the first numerical crossings found in this prescribed family,
not global lower bounds, universal recruitment requirements or independent
validation of the coupling hypothesis.

\subsection{The separate extension beyond 600 s}
Continuing constant full stimulation to 3600 s gives the following null
diagnostics. Each run passes all inherited gates through 600 s. Intracellular
pH is the only later gate failure; all other monitored gates continue to pass.
\begin{center}
\begin{tabular}{lrrr}
\hline
Case & First pH $=7.3$ time (s) & pH at 3600 s & Limiting local root pH \\
\hline
'''
for name,title in [('selected','Selected Task 41'),('crossing','Adaptive crossing')]:
    r=long[name];eq=r.get('diagnostic_equilibrium',{}).get('observables',{})
    value=f"{eq['ph_i']:.7f}" if 'ph_i' in eq else 'Unestablished'
    text+=f"{title} & {r['precise_ph_gate_crossing_s']:.5f} & {r['maximum']['ph_i']:.7f} & {value} \\\\\n"
text+=r'''\hline
\end{tabular}
\end{center}
The limiting roots found from these late states are mathematical continuations
outside the inherited pH domain. They do not establish globally unique
equilibria or physiological long term validity. These separate diagnostics do
not change any frozen 600 s Task 41 result.

\subsection{Local uncertainty of the crossing}
For a locally regular numerical crossing with $\partial_\rho D_X\ne0$,
\begin{equation}
 \frac{\mathrm d\rho_X}{\mathrm d\log p}
 =-\frac{\partial_{\log p}D_X}{\partial_\rho D_X},\qquad
 \frac{\mathrm d\log b_X}{\mathrm d\log p}
 =-\frac{1}{1-\rho_X}\frac{\mathrm d\rho_X}{\mathrm d\log p}.
 \label{eq:thresholdsensitivity}
\end{equation}
Eight effective or architecture derived quantities were selected using the
Task 43 provenance classification. A perturbed WT denominator is used with
the same perturbed parameter in the null. The conserved initial state is
always held at the inherited WT value. Thus changing the buffer amount also
changes the algebraic initial pH; it does not trigger a new resting solve.

The following elasticities refer to the adaptive crossing. The companion
sampled results are retained in the machine readable crosswalk.
\begin{center}
\begin{tabular}{lrr}
\hline
Parameter & $\mathrm d\log\rho/\mathrm d\log p$ & $\mathrm d\log b/\mathrm d\log p$ \\
\hline
'''+ '\n'.join(sr)+r'''
\hline
\end{tabular}
\end{center}
The derivative calculation is checked against independently solved nearby
crossings at logarithmic perturbation sizes $10^{-3}$ and $5\times10^{-4}$.
Every perturbed WT and null trajectory passes its inherited gates. The maximum
relative discrepancy between implicit crossing derivatives and nearby root
differences is '''+f"{max(r['maximum_validation_relative_error'] for r in sens):.2g}"+r'''.
The larger relative sensitivity of residual recruitment $b$ is relevant when
interpreting a conductance assay. None of these eight quantities has an
established experimental uncertainty interval in the audited configuration.
Accordingly, these calculations establish local sensitivity only; they do
not prove robustness over a finite range or identify a finite parameter
change that overturns the conclusion.

\subsection{A discriminating experiment and scope}
At the numerical comparison crossing the null must retain approximately
10.49\% of stimulated CaCC conductance in this family. A direct test would
compare stimulated CaCC conductance and channel recruitment in WT and AE4
loss cells while matching calcium activation, voltage, intracellular chloride
and pH, so that a changed driving force is distinguished from a changed
conductance. Surface channel abundance and channel open probability would
help distinguish recruitment from gating. Concurrent measurements of pH,
chloride, cell volume and secretion over the first ten minutes and after
approximately 25 minutes would test the predicted coupled response and its
later pH failure. A secretion measurement alone does not identify this
specific dependency.

'''
(HERE/'inverse_results.tex').write_text(text)
print('Inverse report fragment generated from verified results')
