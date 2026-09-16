# Task52: accepted Task37 provenance and dependency boundary

**52A file/Git/static-source audit: PASS.** No model import, model evaluation, trajectory, stationary solve, fit or parameter search occurred. This audit applies the Task52 prompt and D52-03/04/08 under R36–R40, R51H, R51I and R52A.

The accepted pre-Palk WT publication is `2ba0e4e0c2a9e940aa96ed86c3c72456afd0e3aa`, tree `daf5e2cd3042d0dbd7734bcfd2fca4dc711f634e`, prepared from `d925785512a33b241ec24843660ff3ef88458549`. That historical commit was fetched into the local Git object store without changing the Task52 branch. All 16 transitive scientific modules match their exact published blobs in Task37 `source_verification.json`. Input and saved-result bytes are verified in `output/task37_provenance_audit.json`.

The exact historical import closure is:

`__init__`, `acid_base`, `ae4_routing_only`, `camp_pka`, `membranes`, `model`, `nbc_minimal`, `nhe1_cha2009`, `nkcc_stimulation`, `parameters`, `states`, `task30_nhe1_repair`, `task31_nhe1_repair`, `transporters`, `validation`, `water`.

The historical closure contains no Palk NKCC, Task41 recruitment, Task50 coupling or later equal-cation-routing dependency. Static AST inspection finds no dynamic-import, `exec` or `eval` calls. Current `membranes.py`, `model.py`, `nbc_minimal.py` and `nkcc_stimulation.py` differ from Task37; current modules are therefore unsuitable for a literal no-Palk dependency claim. After remotely verified 52A, the plan is to copy the exact historical modules under Task52 `frozen_task37/src/modern_full_model` and load them as an isolated package. The historical scientific source remains byte-exact.

The loader additionally needs exact historical `reference/R09_wt_rest.json` and `reference/native_wt_contract.json` at `frozen_task37/reference`. Its existing relative path then resolves correctly. The entrypoint is:

```python
background = task31_nhe1_repair.load_background("R09")
parent = task31_nhe1_repair.build_model(
    background, carrier_amount_fmol=2.339370005697548e-05, stimulated=True
)
model = MinimalNbcModel(parent)
```

The accepted onset core comes from `results/31_nhe1_mechanistic_repair/rest_checkpoints/185f5cee3bcec7a023cf82e4092872114f2292732ad9c4bab68ecf4e062f76e3.json`, labelled `R09_WT_density_5`, followed by regulatory coordinate zero. The earlier reference JSON and `model.initial_state()` do **not** supply the accepted core. Static canonical-JSON reconstruction independently matches the active parameter hash without constructing a model.

| Inherited object | SHA-256 |
| --- | --- |
| Accepted 12-state core | `d77c1907e3e7124061721d606148a733f858cce2ef591c4e9f16e5cb7fc2469d` |
| Accepted core plus regulator | `879a650761e9044691ab6b2a9b1bfc3be8f951ece176925117596cc2baa31abe` |
| Active whole-cell parameters | `ac05ebb3d2b2dde30a94eb8d928b363ba128c8c5e2a099374373be5446146ca0` |
| AE4 parameters | `6c634251c09a6a06d697aa7a3b4ec20d4ee69db039537df1506c22d38a561e5a` |
| NBC parameters | `4fd085686c98ae3e1aa227f38670d28f1aafb8059841b68557b077e43642eaad` |

Task37 uses reversible bounded NKCC/AE2 `capacity*tanh(log(activity ratio)/2)`, exact NKCC `1 Na:1 K:2 Cl` and AE2 `1 Cl:1 HCO3`. NKCC capacity is `0.32 fmol/s`; the inherited calcium multiplier reaches `1.75`. Cha NHE1 retains the WT-pH-calibrated carrier amount above and source-fixed stimulated multiplier `2.3`. Stimulus-recruited NBC retains `1 Na:2 HCO3`, capacity `0.11570913197464398 fmol/s`, width `2`; its capacity is a derived WT feasibility scale, not a measured transporter abundance. AE4 keeps the inherited effective beta gain `1+0.25*a`, with its existing 30-second activation state initially zero. No new signalling state is added.

**Routing qualification:** Task37 retains donor-weighted AE4 cation allocation `r=Na/(Na+K)`; donor is intracellular for nonnegative net chloride loading and extracellular for negative loading. Its source is `(-r*A, -(1-r)*A, A, -2*A, -2*A)`. This is not the later fixed 50:50 routing. The historical `AE4NA05` root label is metadata, not 5% AE4 expression. WT expression is exactly one. Consequently the later R51G/Task40 `J=6P-H` stationary identity must not be transferred unchanged to Task37. The mandatory chloride-reservoir identity is independent of cation allocation and remains applicable.

The production protocol is CCh `0.3 uM`, IPR `5 uM`, calcium `0.25 uM` and beta occupancy one, for 600 seconds. It is REST at exactly zero and stimulated at positive time; the inherited integration begins at `1e-6 s` with the unchanged initial state. Radau uses `rtol=1e-7`, amount/regulatory `atol=1e-10`, volume `atol=1e-12 pL`, `max_step=2 s`. The historical numerical fallback is BDF at the same tolerances after explicit solver failure only, never after physiological failure. Task52 freezes its own execution/failure policy before production.

Inherited production physical limits are finite states/RHS/diagnostics; all twelve conserved core coordinates positive; `Na<40 mM`, `50<=K<=200 mM`, `30<=Cl<=80 mM`, `6.6<=pH<=7.3`, `0<HCO3<100 mM`, `TIC>0`, and `0<Vcell<3 pL`. Accepted endpoints and integer-second dense points were checked. Eleven conservation tolerances remain: `1e-10 fmol/s` for charge/carbon/buffer/homeostasis/AE4 residuals, `1e-12 pL/s` for water, `1e-20 A` for each current closure, and `1e-9 mM` for cell/lumen alkalinity speciation. NBC diagnostic source fields equal `-J_NBC` charge and `2*J_NBC` carbon; they are not zero residuals.

Task37 headline quadrature used the one-second trapezoidal grid plus zero and `1e-6`, with windows 0–600 and 60–600 seconds. Only its 60-second readout table was retained; the complete one-second reference trajectory was not saved. Task52's independent reservoir quadrature is separately predeclared and must not be described as an inherited Task37 gate.

The accepted Task37 rest has Na `11.636125748680639 mM`, volume `1.4484260832245268 pL`, Cl `60.30496692587399 mM` and pH `6.910000223494689`. Task52 must still apply the prescribed measured-state projections; these are onset states, not claimed equilibria of the incomplete chronic model. Task37 was accepted with its NKCC chloride partition outside the contextual band, and its prior calibration/provenance limitations remain. No archived solve/fit utility is authorised by its presence in the copied package.
