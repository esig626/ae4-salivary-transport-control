# Task 31: Cha NHE1 source, equations and reconstruction

Primary source: [Cha et al., Biophysical Journal 97 (2009), 2674–2683](https://doi.org/10.1016/j.bpj.2009.08.053), [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC2776256/).
The article identifies the supplement as “Document S1. One figure and three tables.”
The adopted family is eight-state ion exchange, Mod2, n_H=3 and m_H=1.

## Retrieved supplement and exact transcription

- File: `mmc1.pdf`, seven pages, 231386 bytes, retrieved through the browser on 2026-09-13.
- Actual source URL: https://pmc.ncbi.nlm.nih.gov/articles/instance/2776256/bin/mmc1.pdf
- Original supplement endpoint: http://www.biophysj.org/biophysj/supplemental/S0006-3495(09)01444-1
- SHA-256: `3fdff54677a5c3347045de9e3d1efd52b62cbc2455400a37db1595e9f2fc146c`
- Table S1, PDF page 5: extracted text was cross-checked against the rendered page.
- Exact values, units, URL and hash were recorded in `results/31_nhe1_mechanistic_repair/source_table_s1.json` before using the production constants.

| Table S1 symbol | Final value | Units |
|---|---:|---|
| k1+ | 10.5 | ms^-1 |
| k1− | 0.201 | ms^-1 |
| k2+ | 15.8 | ms^-1 |
| k2− | 183 | ms^-1 |
| K_H_i | 6.05e-4 | mM |
| K_H_o | 1.62e-3 | mM |
| K_Na_i | 16.2 | mM |
| K_Na_o | 195 | mM |
| n_H | 3 | dimensionless |
| m_H | 1 | dimensionless |
| K_i | 3.07e-5 | mM |
| K_o | 4.8e-7 | mM |
| N | 489900 | molecules per source cell |

N is recorded as the source normalization only. No cardiac carrier count,
geometry, passive flux, pump, calcium mechanism or other background mechanism
is imported. Transfer of kinetic shape from the source preparations to salivary
acinar cells remains an unvalidated cross-cell-type assumption.

## Equations and amount conversion

Let u_o=Na_o/K_Na_o, v_o=H_o/K_H_o, u_i=Na_i/K_Na_i, v_i=H_i/K_H_i.
Define D_o=(1+u_o)(1+v_o), D_i=(1+u_i)(1+v_i), and

```
a = k1_plus * u_o / D_o
b = k2_plus * v_i / D_i
c = k1_minus * u_i / D_i
d = k2_minus * v_o / D_o
j_exchange = (a*b-c*d)/(a+b+c+d)               [Eq. 3, ms^-1]
Mod2 = 1 / (1 + (1+H_o/K_o)*(K_i/H_i)^3)       [Eq. 5]
J_amount = carrier_amount_fmol * 1000 * Mod2 * j_exchange
```

The final line is Eq. 1 in the repository's amount coordinates. Concentrations
are mM, so H=10^(3-pH). Carrier amount is N/Avogadro times 1e15; one fmol
corresponds to 6.02214076e8 carriers. No division by cell volume is required
when evolving amounts. The one salivary parameter is an effective active carrier
amount, combining expression and an overall activity normalization.

Each forward cycle adds one Na and removes one H. The existing amount ledger
adds (+J,0,0,0,+J) to cell (Na,K,Cl,TIC,alkalinity): zero carbon and net charge.
The same cycle supports reversal. No genotype term or new acid/carbon source
is present. The Task 30 Vera-Siguenza Eq. 26 evaluator remains selectable.

Binding and modifier probabilities are evaluated in logarithms. For finite
fixed rates, |j_exchange| <= max(min(k1+,k2+),min(k1−,k2−)); the modifier lies
between zero and one. The production conservative turnover bound is 10500 s^-1.

## Printed precision and source reconstruction

Eq. 6 requires

```
k2_minus = k1_plus*k2_plus*K_H_o*K_Na_i / (k1_minus*K_H_i*K_Na_o)
```

With the other printed entries, Eq. 6 gives **183.6074807146718 ms^-1**,
0.3319567% above the printed **183**. Production retains every printed Table
S1 constant exactly. It evaluates Eq. 3's actual product difference, rather
than substituting ideal chemical affinity and silently changing the equation.
The resulting cycle affinity is ln(Na_o H_i / (Na_i H_o)) + ln(183.60748/183).
This shifts its zero by 0.00143928 pH units and leaves a small nonzero flux at
exact chemical equilibrium. Consequently exact thermodynamic consistency is
not claimed for the rounded printed set. Synthetic fixtures continue to use
Eq. 6 exactly; they are explicitly separate from production parameters.

An independent rational-arithmetic evaluation reconstructs the printed Eq. 3
and Mod2 at the main-text turnover conditions and the Fig. 4 C/D forward/reverse
conditions to relative agreement better than 1e-12. Source Na=0 limits use
1e-12 mM. No cardiac-cell figure dynamics or unreported fitted coordinates
were reconstructed.

At pHi=6.0, pHo=7.4, Na_i approximately zero and Na_o=140 mM, the reconstructed
regulated turnover is **2.521212543779516 ms^-1**. The main text states a fitted
restriction of 2.3–2.5 ms^-1: the printed parameter reconstruction exceeds its
upper endpoint by 0.8485%. This restriction is therefore **not reproduced
exactly**. Equation/transcription tests pass; the source-band check is recorded
as false. Neither discrepancy was corrected by changing kinetic parameters.

## Single salivary density and bounded execution

WT R09 alone sets the active carrier amount using a safeguarded one-dimensional
secant in log(amount), target pH 6.91, root tolerance 1e-5. The numerical starting
amount was 1e-4 fmol and the predeclared bounds were [1e-8,0.1] fmol; these are
solver choices, not source-density measurements. Only five density evaluations
were needed. Chloride, null measurements and secretion were not used.

The frozen amount is **2.339370005697548e-5 fmol**, or **14088.015464032636
carriers per salivary cell**, **0.028756920726745534 times** Table S1's source-cell
N. This is an amount/count ratio, not a membrane-area density ratio or a claim
of measured salivary expression. The same amount and all ten printed kinetic
constants were used in R10 WT and both exact AE4-null solves.

REST uses the inherited ten-coordinate charge manifold, coordinate bounds,
residual scales, conservation/current/charge/positivity/rank gates and unchanged
production transporter laws. Each solve permits at most 2400 actual residual
calls including finite-difference and final-gate evaluations. There are no
alternate starts or null-density adjustments. Chloride context is reported,
not used as a gate or calibration objective. The null pH diagnostic is
6.87–6.91; WT is 6.84–6.98. TIC/HCO3 are reported individually and checked
against the inherited broad 100 mM carbon ceiling.

Eight stationary calls used 8297 residual evaluations. Both null attempts
reached their local limits without numerical closure. The final candidates
are not resting states and do not establish nonexistence of a root. No matched
pair was admissible, so no stimulated integration ran. The bounded outcome
is reported without another search or transporter repair.
