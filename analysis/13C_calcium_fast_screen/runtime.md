The primary screen contains exactly twenty trajectories. Its serial pilot is counted once and reused; the remaining nineteen run through a process pool with nine workers and one BLAS thread per worker. Each request runs the existing 600 s loose Radau protocol. A single 0.50 µM root then receives production Radau and BDF confirmation, for twenty two new integrations overall.

| Run | Wall time (s) | nfev | njev | nlu | Samples |
| --- | --- | --- | --- | --- | --- |
| Serial 0.25 pilot | 0.266816 | 1013 | 2 | 252 | 122 |
| production_radau decisive 0.50 | 0.591432 | 2409 | 3 | 604 | 122 |
| production_bdf decisive 0.50 | 0.269730 | 828 | 1 | 70 | 122 |


The pilot integrator itself took 0.244071 s. The full pilot took 0.266816 s, far below the predeclared 120 s profiling threshold. No numerical optimisation was needed. Counter capture wraps the existing solve call and observes its return values without changing its arguments or equations.

The measured remaining screen phase took 1.169619 s, including cached pilot loading and output work. Adding the earlier pilot gives 1.436435 s for primary computation. The confirmation phase took 0.999082 s, giving 2.435517 s overall. These are numerical phase timings; repository retrieval, implementation, interpreter startup, documentation and tests are excluded. Per trajectory parallel wall times include CPU contention and must not be added to claim elapsed process pool time.

The original nine point primary design required ninety production Radau trajectories. It was not run. The exact primary trajectory reduction is 90 to 20, or 77.8%; including the two confirmations, 22 trajectories are 75.6% fewer than that old primary stage alone. A rough same machine estimate uses the observed production Radau cost 0.591432 s and ten waves of nine workers: about 5.914 s for the original primary stage, excluding process and output overhead. On that estimate, primary computation saves about 4.478 s (75.7%), or about 3.479 s (58.8%) when all fast screen confirmations are counted. This is an estimate based on one production root and calcium level, not a measured counterfactual; costs can vary with calcium, root and concurrent load. No unperformed control or ensemble runs have been counted as measured savings.

Recommended invocation, allowing the worker count to follow available CPU affinity and the ten root cap:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
python -m src.modern_full_model.run_calcium_fast_screen
```

This run used `--workers 9` with separate `--phase pilot`, `--phase screen` and `--phase confirm` calls. The default `--phase all` follows the same order. Saved completed requests are reused with hashes checked. Generate the reports from saved files with `python -m src.modern_full_model.report_calcium_fast_screen`.
