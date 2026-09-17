# Task52 independent production-runner static review

Scope: `run_frozen_cases.py` and its already verified `diagnostics.py`, after
verified 52B and before any Task52 production trajectory. The reviewer read
the new R52B ledger entry and the binding 52A numerical freeze. No scientific
module was imported or evaluated; no model test, projection or trajectory ran
in this review. This is a software inspection, not a numerical result.

Reviewed diagnostics SHA-256:
`11b5ff66e35052672072be76a64a52a0d9bc265e78f2da3eb36482b9a7b35f30`.
Final reviewed runner SHA-256:
`d5988680b99638fa91ed0e6e82a70f6b1d6bca581ee102438060d6e1025909ce`.

## Scientific dispatch and frozen controls

The stiff solver receives `PairedModel.rhs` for the complete 26-coordinate
state. WT and KO therefore advance together and the WT-derived signed supply
is recomputed at every paired evaluation. Onset vectors come directly from
the immutable central or alternate projection record. Only the three 52D
and five 52E cases in the frozen matrix are selected. BLAS limits are set
before numerical imports and reinforced through `threadpool_limits(1)`.

The branch, prior publication receipt, local HEAD, freeze hash and immutable
artifact hashes are checked before production. 52D requires the verified 52C
receipt; 52E requires the verified 52D receipt. These are consumers of the
separate independent publication-verification process, not a replacement for
that remote verification. The 52C freeze and receipt are not yet present at
this review boundary and must include the final runner hash. Stage and case
start records reject accidental repeated execution.

## Physical failure and persistence

Each accepted step's integer samples and endpoint are checked in increasing
time. A failed physical/conservation sample stops further solver stepping.
It preserves its state, failure details and the last passing sample time.
The prescribed BDF fallback is used only after an explicit Radau solver-status
failure; physical failures, exceptions and reservoir failures cannot trigger it.
Incomplete or invalid cases cannot populate the 600-second headline fields.

The first draft could overwrite a physical failure if subsequent quadrature
raised an exception, and could leave one quadrature rule advanced while the
other was incomplete. Before production, the root corrected this by preserving
the first failure, recording secondary exceptions separately, and committing
the two quadrature increments atomically. Quadrature coverage is explicit;
incomplete coverage forces the budget gate to fail and leaves the complete
`gauss8_integrals_0_T` unavailable. Partial integrals have their own endpoint.

## Independent chloride-reservoir audit

The four- and eight-node Gauss-Legendre rules independently evaluate actual
paired fluxes on the solver's dense polynomial for each accepted interval.
Neither rule infers export from chloride storage. Signed total apical J,
NKCC cycles N, ordinary AE4 chloride A and signed AE2 E enter the specified
WT-minus-KO identity; the two-N stoichiometric factor is retained.

The root added an explicit eight-minus-four discrepancy for `2*D_N` as well
as N cycles, applying the unchanged `1e-6 fmol` gate to the actual chloride
mass term. Individual-cell and differential discrepancies are both checked.
The independent budget error must remain within `1e-5 fmol`.

The unsimulated [0,1e-6] interval holds the onset state fixed. Its right-limit
flux contribution is separately recorded and its conservative differential
flux-times-duration bound enters the same reservoir error gate. The finite
onset convention is not silently represented as an exactly integrated state.

## Reporting windows and source labels

Headline secretion and chloride export use only integer-second samples plus
the separate zero/epsilon onset samples, with the frozen trapezoidal rule.
They report 0–600 seconds, 60–600-second mean flow, endpoint readouts and all
predeclared cumulative observations. A partial final point is used only in
explicitly partial output. The coarser trapezoidal budget error is reported
separately from the stricter Gauss audit and pointwise conservation residuals.

Total chloride current is converted with anion valence minus one. Native
TMEM16A and auxiliary components remain separate in trajectories and integrated
flux records. Reservoir delta uses the two chloride **amount** coordinates;
the concentration derivative separately includes dilution by changing volume.

## Review disposition

The inspected logic preserves the prescribed model, solver and acceptance
rules. The exception/coverage and two-N bookkeeping corrections above occurred
before production and change no scientific parameter, equation or tolerance.
The root also implemented the requested preservation correction: every raw
accepted solver endpoint is retained separately from observation rows, with
a validation mask. An endpoint beyond an earlier failed integer sample remains
stored but unvalidated, with no further diagnostic evaluation or integration.
It does not enter headline or valid-window results. No blocking static-review
finding remains in the final runner identified above. Runtime numerical and
conservation outcomes remain pending the remotely verified 52C boundary.
