# Task52 architecture and numerical freeze — 52A

Operational start: `eacfe603faa3ce0b2aa26814b348a1ae44417cdc` on
`analysis/task-52-chloride-reservoir-final-test`.
Prepared canonical merge: `d2b14c2821dd7e38703620389097cc8f1052977a`.
Scientific chassis: accepted Task37 snapshot
`2ba0e4e0c2a9e940aa96ed86c3c72456afd0e3aa`.

The complete CODEX_START_HERE binding stack, including R51H and R51I,
has been read. R51H supersedes the R51G prospective beta-NKCC/swelling
recommendation. Task52 tests only the measured chloride reservoir, ordinary
WT AE4 supply, matched WT-derived NKCC/AE2 supply, and shared beta demand.

## Exact scientific inheritance

After verified 52A, materialise a byte-exact Task37 dependency snapshot under
this task's `frozen_task37/` directory. Its isolated import path must exclude
later production modules. This avoids importing optional Palk code through
current generic selectors. Original source and frozen results remain intact.
Use the Task31 R09 factory and accepted rest used by Task37, its mechanistic
Cha NHE1, recruited NBC, ordinary beta/PKA AE4 regulation, original
donor-weighted AE4 cation routing, pre-Palk NKCC, and all other inherited laws.

**Provenance qualification:** Task37 AE4 cation allocation is donor-weighted,
not the later Task40 equal allocation. The Task52 prohibition on changing
routing requires retaining it. R51G's equal-routing stationary pump/NHE
identity is not imposed on this chassis. The required chloride-reservoir
identity is unchanged and exact regardless of cation routing.

The paired KO receives the WT's instantaneous signed NKCC and AE2 cycle
fluxes through their complete electroneutral stoichiometric vectors. This is
a MODEL IDEALISATION, not experimental proof of exact native flux equality.
WT AE4 expression is one; KO AE4 expression is zero. All other local transport,
water, voltage and lumen responses evolve independently within the pair.
The same auxiliary law `I_aux=beta*G_aux*(V_a-E_Cl)` applies to both cells.

Source Cl/pH values are imposed only by the prescribed central and alternate
charge/osmotic projections. No pre-stimulation relaxation or stationary solve
may remove them. These are onset states, not predicted chronic equilibria.

## Frozen numerics and reporting

- Radau, rtol `1e-7`; amount/regulatory atol `1e-10`, volume atol `1e-12 pL`;
  maximum step `2 s`. Coordinate tolerances are tiled for the paired system.
- Inherited exact-zero REST/right-limit stimulation convention; start the
  continuous integration at `1e-6 s` from the unchanged projected onset state.
- One scientific process, one thread per BLAS pool; no parallel trajectories.
- Preserve every accepted solver endpoint and integer-second dense sample.
  Apply unchanged Task37 positivity, broad physiological and conservation
  gates to both genotypes. Stop a case at its first detected failed sample;
  preserve partial data, do not substitute an invalid 600 s result.
- The inherited BDF fallback is permitted once only for an explicit Radau
  solver-status failure, never for a physical/conservation gate or exception.
  No changes of inputs or tolerances accompany it. Every attempt is retained.
- Headline cumulative flow/export: inherited one-second trapezoidal quadrature
  with separate 0 and 1e-6 onset samples. Report 0–600 s, mean 60–600 s,
  endpoint flow, and cumulative observations at 60,120,...,600 s. No timing fit.
- Mandatory reservoir audit: independently evaluate signed J,N,A,E at fixed
  8-point Gauss–Legendre nodes on every accepted dense-output step, and use
  a 4-point rule as a numerical cross-check. Neither integral is inferred
  from chloride amount changes. Require absolute identity error <=`1e-5 fmol`
  and 8-versus-4 quadrature discrepancy <=`1e-6 fmol` for every budget term.
  The negligible unsimulated [0,1e-6] interval is explicitly zero-duration
  state evolution, with separately reported onset quadrature contribution.
  Its maximum flux-times-duration bound must fit the same `1e-5 fmol` gate.
  Report the coarser trapezoidal balance error separately; it is not a
  pointwise conservation tolerance.
- Projection gates: unchanged physiological ranges; cell charge <=`1e-9 fmol`,
  bath-isotonic residual <=`1e-9 mOsm`, imposed pH error <=`1e-10`.
- No model evaluation at 52A; software/fixed-state work only at 52B;
  immutable inputs at 52C; central cases at 52D; fixed controls/sensitivities
  at 52E; experimental/Task50 comparison and final classification at 52F.

The central-projection failure rule in the prompt is a hard stop before
production. A later trajectory failure is preserved as a failed case and
does not license rescue or additional cases. Remaining predeclared cases
may be executed at their authorised checkpoint, with invalid endpoints null.
The final result must not manufacture a successful classification from
failed or unavailable trajectories.

The eight cases are exactly D52-07 and `output/case_matrix_52A.json`.
The frozen Task50 archive is benchmark-only. No AE4 target, optimisation,
additional mechanism, new regulatory state, or post-output retuning is allowed.
