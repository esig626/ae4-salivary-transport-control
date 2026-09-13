# Task 36 design: minimal Na/HCO3 alkalinity pathway

## Why this exists

Task 33 established that increasing routed AE4 capacity does not materially increase AE4 chloride loading. The steady alkalinity balance explains why:

`J_AE2 = J_NHE1 - 2 J_AE4`

and therefore the net HCO3-dependent chloride loading is

`J_AE2 + J_AE4 = J_NHE1 - J_AE4`.

A large AE4 chloride-loading flux therefore requires a comparably large independent source of intracellular alkalinity. Reverse AE2 cannot supply that alkalinity without simultaneously exporting chloride.

For a target NKCC1 fraction `0.70`, using the accepted Task 31 R09 resting chloride-loading scale

`J_NKCC1,Cl = 0.2562444047732706 fmol/s`,

the required positive non-NKCC chloride-loading flux is

`(3/7) J_NKCC1,Cl = 0.10981903061711597 fmol/s`.

Using the source-fixed stimulated-NHE hypothesis `J_NHE1 = 2.3 * 0.007890580489555009 = 0.01814833512597652 fmol/s` and allowing AE2 to contribute at most its existing positive capacity `0.005 fmol/s`, AE4 must still provide at least

`J_AE4 = 0.10981903061711597 - 0.005 = 0.10481903061711596 fmol/s`.

Steady alkalinity balance with an added chloride-neutral alkalinity source `J_NBC` is

`J_NHE1 + J_NBC - J_AE2 - 2 J_AE4 = 0`,

so the minimum required new alkalinity influx on that reference scale is

`J_NBC = 0.1964897261082554 fmol/s`.

This is a feasibility requirement, not a measured NBC flux.

## Minimal transporter

The new module `src/modern_full_model/nbc_minimal.py` implements one deliberately simple structural surrogate:

- basolateral;
- reversible;
- electroneutral `1 Na : 1 HCO3` cotransport;
- positive direction bath to cell;
- no direct chloride source;
- no membrane-current term;
- one fixed reference capacity, no optimisation.

The flux law is

`J_NBC = G_NBC tanh(A_NBC / w)`

with

`A_NBC = log((Na_o HCO3_o)/(Na_i HCO3_i))`.

At the Task 31 R09 chemical state this affinity is `4.100247118100945`; with `w=2`, the required `0.1964897261082554 fmol/s` flux corresponds to the fixed reference capacity

`G_NBC = 0.20311053520552247 fmol/s`.

That capacity is labelled explicitly as a **derived feasibility scale**, not a measurement and not an inferred transporter abundance.

## What this model is not

It is not yet an assignment to NBCe1, NBCn1, or another SLC4 protein. Evidence for Na/HCO3 transport exists in salivary acinar preparations, especially parotid, but the specific identity and abundance relevant to the mouse submandibular AE4 phenotype are not established strongly enough here to pretend otherwise.

It is also not a reinterpretation of AE4 as an independent Na/HCO3 cotransporter. Peña-Münzenmayer et al. 2016 support an electroneutral cation-dependent Cl/HCO3 exchanger. In the physiologically proposed chloride-loading direction, AE4 imports chloride while exporting cation plus bicarbonate. The reverse direction can import Na/HCO3 while exporting chloride, but that cannot explain the lower intracellular chloride in AE4-null secretory cells if it were the dominant physiological direction.

The purpose of the surrogate is therefore narrower: test whether the missing degree of freedom is a chloride-neutral alkalinity-loading pathway. If adding this pathway allows AE4 to carry substantial productive chloride flux without forcing reverse AE2 cancellation, the acid-base obstruction has been isolated.

## No search

Do not sweep this capacity. Do not combine it with alternative stoichiometries or additional transporters in the same task. First test this one structural intervention. If it fails, report that failure before considering a different mechanism.
