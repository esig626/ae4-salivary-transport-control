# Dynamic reconstruction handoff

## Status

This note closes the current CBM phase as a mechanistic handoff. It does **not** claim completion of Tasks 46H05–46H07. The published branch contains H00–H04 plus preserved H04 parallel work. No reserved AE4 secretion phenotype is used to derive the identity below.

## Meaningful structural result

Using the frozen canonical CBM notation

- `J`: apical chloride export,
- `N`: NKCC1 cycle flux,
- `H`: NHE1 flux,
- `A`: AE4 flux,
- `P = P_a + P_b`: total Na/K-ATPase cycle flux,
- `K = K_a + K_b`: total K-channel efflux,

intracellular steady-state balances give

`J = 2N + H + 2B - A`

and

`P = (H + J)/6`.

Substitution into the potassium balance gives the exact identity

`K = J/3 + N + H/3 - A/2`.

For AE4 knockout, `A = 0`, so

`K_KO = J/3 + N + H/3 >= J/3`.

This provides the mechanistic lead for the next dynamic reconstruction. At a fixed chloride secretory output, AE4 activity reduces the K-recycling requirement by `A/2`. Removing AE4 removes that relief. Any increase in NKCC1 used to replace missing AE4 chloride adds a one-for-one additional K-recycling burden; NHE1 activity adds a further `H/3` burden.

This does **not** prove that K recycling is the physiological limiting mechanism. The CBM contains no membrane potential, concentration-dependent electrochemical driving force, channel gating, or kinetic pump saturation. It only identifies the coupling that a dynamic model must test.

## What H00–H04 established

1. A single conservation-explicit CBM was constructed and audited for conserved Na, K, Cl, total inorganic carbon, total alkalinity and charge bookkeeping.
2. WT reference chloride export is feasible.
3. AE4 deletion alone still permits the full WT reference chloride export under unchanged shared computational capacities. Therefore the conserved network by itself does not force a secretion deficit.
4. The independent NKCC assay cannot defensibly be converted into an exact WT/KO stimulated NKCC equality or an invented percentage tolerance. It measured initial chloride uptake under bicarbonate-free, carbonic-anhydrase-inhibited and CaCC-inhibited conditions, whereas the CBM is a sustained zero-storage model.
5. Consequently the full-rescue CBM witness survives H04. The missing physiology must be sought in a dynamic constraint rather than manufactured as an unsupported NKCC cap.

## Article narrative

The intended narrative is therefore:

1. The conservation-explicit dynamic model reproduces much of WT physiology but permits excessive compensation after AE4 loss.
2. A phenotype-blind structural CBM shows why: AE4 deletion does not make normal chloride export stoichiometrically infeasible. The remaining network can rearrange fluxes and preserve the reference output.
3. Independent NKCC uptake data do not justify imposing an arbitrary stimulated NKCC cap.
4. The CBM instead exposes a specific systems-level cost of compensation: the exact K-recycling identity above.
5. The next question is dynamic and physiological: whether K-channel conductance, membrane voltage, intracellular K, Na/K-ATPase kinetics and NKCC1 jointly impose a finite recycling ceiling that prevents unlimited AE4 replacement.

This is a mechanistic hypothesis to test, not a parameter-fitting instruction. The AE4-null secretion magnitude must not be used to tune that ceiling.
