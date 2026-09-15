# Milestone 01 — unchanged Task 44 failure baseline reproduced

All 1,794 recorded numerical comparisons with the pinned Task 40/44 results
match bit for bit. No equation, parameter, initial state, stimulus or source
vector was changed. The original production calculation uses 13 dynamic
coordinates, including cell/lumen carbon, alkalinity, volumes and regulation;
pH and membrane potentials are solved algebraically. The stationary problem
retains 11 independent coordinates on the two exact charge constraints.

## Reproduced secretion and compensation

| Quantity | WT | 5% AE4 | AE4 null |
|---|---:|---:|---:|
| Cumulative secretion, 0–600 s (pL) | 0.992524544081835 | 0.958604247509151 | 0.954238045797850 |
| Cumulative deficit (%) | 0 | 3.417577608 | 3.857486297 |
| Stationary secretion (pL/s) | 0.00158406711086 | 0.00157144589556 | 0.00156907909882 |
| Stationary Na_i (mM) | 18.08263257 | 17.17284377 | 17.08443833 |
| Stationary Cl_i (mM) | 53.29853853 | 52.28506820 | 52.12870969 |
| Stationary pH_i | 7.05038628 | 7.24518552 | 7.27231514 |
| Stationary cell volume (pL) | 1.38987955 | 1.51684999 | 1.53780201 |
| Slowest local decay (s) | 345.81857188 | 442.08904428 | 459.91603104 |

The stationary null secretion deficit is **0.9461727943%**, a different
observable from the 600 s cumulative deficit.

Over **60–600 s**, NKCC1 supplies 36.3370075332 fmol additional chloride after
AE4 deletion, against 41.0555019619 fmol missing AE4 chloride. Thus the
replacement fraction is **88.5070351033%**, while NKCC's relative integrated
increase is **23.1634495075%**. NKCC cycles carry two chloride ions. These are
signed flux integrals on the full one-second grid, not integrals of the compact
60-second display table.

The local WT equilibrium compensation gain, computed with the unchanged
independent Task 44 IFT verifier, is **92.9010741779%**. The original Task 44
analysis calculation reports 92.9010742738%; the small difference is between
two numerical differentiation implementations, not a model change. The
independent implementation is matched to its own pinned reference.

## What was rerun and verified

- Original 13-state production Radau stepping for WT, 5% and null, all from the
  identical frozen WT rest. Each has 934 checked states: accepted endpoints
  plus the one-second grid and the two onset records.
- Original solver controls: rtol 1e-7, amount/regulation atol 1e-10, volume
  atol 1e-12, maximum step 2 s. Resting observation at t=0 and stimulated
  integration beginning at 1e-6 s are preserved.
- All full inherited physiological, conservation, NBC source and equal-routing
  source-signature checks pass. WT activation checks pass.
- A separate constant-stimulus WT integration supplies the equilibrium seed.
  The original C0 continuation uses eight expressions. All eight roots are
  physiological and locally stable on the charge manifold. Four independently
  solved nearby roots verify the WT IFT derivative at two steps.
- All state/output references and individually matched eigenvalues are
  compared. The regulatory boundary derivative is inward, preserving the
  correct 30 s regulatory mode.
- Independent saved-output review reevaluates retained states and checks
  volume extrema explicitly. Its additional counts are recorded separately in
  `output/checkpoint_01_independent_verification.json`: 67,067 numerical
  comparisons, 3,886 logical checks, 1,806 fresh state evaluations and zero
  failures. Seven independently summed integral differences are roundoff,
  at most 2.85e-14 fmol; the other 67,060 numbers are bit-identical.

The successful replay made four integrations (three production plus one
constant-input seed) and twelve equilibrium solves (eight continuation plus
four nearby verification). A preceding implementation attempt completed one
additional WT integration but failed when serialising a NumPy boolean to JSON.
Its two CSV outputs and failure record are preserved under
`output/checkpoint_01_baseline/` and `output/checkpoint_01_attempt01_failure.json`.
Only JSON scalar conversion was repaired. That attempt is neither a scientific
failure nor a successful complete replay. In total, five integrations were
performed in this milestone. No optimisation or parameter search was run.

## Scope and unresolved limitations

These calculations reproduce the known failure baseline; they do not accept
it as a mechanistic reconstruction. Thermodynamic opposition of C0's
substituted ensemble source vector at the WT equilibrium remains an inherited
qualification. Physiological gates do not remove it. All genotypes share WT
rest; this is not a model of established knockout animals at their own resting
states. The isolated NKCC experiment has a different protocol and has not been
simulated or used for fitting here. Local stable eigenvalues do not establish
global stability, uniqueness or positivity. No justified parameter uncertainty
interval is inferred.

The primary replay comparison record contains 1,794 matches, all exact, rather
than the historical Task 44 full-replay count of 331,157. The raw trajectory,
state, summary, root, IFT and comparison files are in
`output/checkpoint_01_replay/`. Both pinned source checkouts remain unchanged.

Next: milestone 02, causal diagnosis of the excessive compensation, without
changing equations or searching new models.
