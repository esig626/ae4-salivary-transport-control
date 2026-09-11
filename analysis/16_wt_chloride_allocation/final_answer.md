# Task 16 completed bounded study

Primary classification: **TASK 16 NUMERICALLY INCONCLUSIVE**.

All twenty modified continuation paths stop at the inherited Na search boundary before obtaining a valid endpoint equilibrium. Their preceding full-rank numerical states already violate the WT Na/K physiological gates. No modified candidate is WT-admissible, so no modified genotype simulation is licensed. This is a complete accounting of the declared experiment, with an unresolved endpoint-equilibrium limitation rather than a proof of global infeasibility.

| Condition | Reference AE4 share | Complete AE4 capacity multiplier | Complete NKCC1 capacity multiplier | WT-admissible root/calcium cases |
|---|---:|---:|---:|---:|
| Inherited baseline | 1.3553–1.8682% | 1 | 1 | 30/30 |
| AE4 share 0.10 | 10% target | 5.3526–7.3785 | 0.912365–0.917134 | 0/30 |
| AE4 share 0.30 | 30% target | 16.0578–22.1356 | 0.709617–0.713327 | 0/30 |

These targets preserve the reference AE4+NKCC1 positive chloride-loading sum. AE2's negative counterflux is reported separately. They are sensitivity coordinates; no measured physiological allocation is asserted. The realized resting shares at the unresolved modified endpoints are unavailable.

Capacity scaling alone leaves cancellation unchanged at a fixed state. Along every saved modified continuation path, gross opposing Na/K traffic and cancellation increase. Last connected internal states have Na 59.86–60.00 mM and K 61.90–67.86 mM, compared with inherited WT gates of 10–35 and 100–160 mM. Their cancellation is 98.61–99.53%, versus baseline resting cancellation of 96.80–99.01%. These internal states are solver diagnostics, not additional scientific allocation conditions.

After the WT-only checkpoint was pushed, all thirty baseline genotype comparisons were reused. At 5% AE4, total 0–600 s secretion **increases by 8.23–21.66%**. AE2 loss remains comparatively neutral, with **+0.0318–0.1601%** secretion changes. Those baseline findings cover every root, both routing families and calcium 0.10/0.25/0.50 uM. The modified effects are explicitly unavailable; there is no supported sign-change result or threshold estimate.

The declared capacity-only intervention produces no admissible repair. The mixed-cation mechanism merits further investigation, but these bounded continuation failures do not establish its causal role or prove that changing it is necessary. Task 16 does not change that mechanism or begin a further task.

The allocation contract was published at `b1f2697bb20cfe7456bb4ea1d45060bee6669638`. The WT-only checkpoint was pushed at `f5df88448b1a669c7c7b3986d83ac17360b0338d` before genotype evaluation. All ninety root/share/calcium combinations remain in the final tables. Thirty WT and thirty 5% AE4 baseline trajectories, plus validated AE2 ratios, are reused; no new ODE trajectories are generated. The prescribed internal continuation calculations and all unsuccessful attempts are preserved.

Validation: **56 tests passed** (16 Task 16 checks and 40 inherited checks), including coupled source scaling, conservation, thermodynamic invariance, fixed parameters, the phenotype firewall, checkpoint integrity and complete ensemble accounting. All **187 frozen inputs remain unchanged**. Inherited scientific model/source files and `archive/` remain unchanged. Work is confined to `codex/task-16-wt-chloride-allocation`; no merge to `main` is performed.
