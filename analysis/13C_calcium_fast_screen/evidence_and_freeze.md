Task 13C uses the fast screen prompt on branch `codex/task-13c-calcium-input-sweep`, read after `AGENTS.md`. The user's explicit fast screen instruction controls the two new amplitudes, preferred regulator and selective confirmation. It supersedes the broader grid, extra protocol and ensemble provisions in the older prompt.

Branch intake was `289a313fc0c7ba103d490681c43ab3987ca1b62e`. The controlling final Task 13B manifest has SHA256 `6d283961421a6718a54fd7c86e6575bd5c032a93e649923777b02295cee97f29`. Its own hash ledger, saved dynamics profile, group gates and solver comparisons agree. All ten embedded root objects, twelve component core states, whole cell parameter payloads and AE4 parameter payloads pass their individual saved hashes.

Inherited provenance discrepancies were found before simulation. The current native roots CSV is truncated at 1,048,606 bytes and contains none of the ten final roots; its hash and the current native summary hash disagree with the final manifest. The old native source hash ledger also describes a different revision. The complete final manifest contains the exact retained roots and is therefore the intake source. No missing row was regenerated and no root was solved again.

The current historical dynamic runner also differs from its saved source hash. The new screen neither imports nor invokes it. Every other equation source file named in the final manifest matches its saved hash, including model, chemistry, membranes, transporters, water, regulation, numerical validation, native model construction and the NKCC wrapper. The new runner calls these unchanged implementations and verifies its reconstructed parameter objects against the frozen payloads.

| Inherited discrepant path | Expected SHA256 | Actual SHA256 |
| --- | --- | --- |
| `src/modern_full_model/run_native_dynamic_contract.py` | `522fa945e54755ef5bf42519b9c64130516942b7e18e127e9b7b4e9a3fccf187` | `f4380a27f859b6d76368630259d5b3988a3a693ffd8a58520c42457d46966c4f` |
| `results/13B_modern_full_model/native_source_wt_roots.csv` | `c2901a349736c7ba8376ede8499925062282e76614c48c3277251beb16332a99` | `9c3f34e69c54b073f84f0f446b7b026a8238350bd70ee748def13273c11bdec3` |
| `results/13B_modern_full_model/native_source_wt_summary.json` | `d4cf74b59f4c6b861e0c60ac563060083c3f318796195ccc5ee69cd87119b522` | `d6ffd436de28a6dbfdb72f75f0add40d195ca0e0a1840ec742be9e350da283a7` |


Only the protocol's stimulated calcium changes. Resting calcium stays at 0.058 µM. The Ca gate stays C^1.46 / (0.26^1.46 + C^1.46). The NKCC normalisation reference remains 0.10 µM and its gain remains 1.75; both new inputs saturate that same existing NKCC arm. Whole cell transporter capacities, NKCC source scale 4, AE4 routing, membrane conductances, water, geometry, bath and preferred R1 regulation remain unchanged.

The physical protocol remains WT CCH_IPR for 600 s, 0.3 µM CCh and 5 µM IPR, with the frozen beta step. The grid is basal zero, the right limit at 0.000001 s, then 5 s increments through 600 s. No regulatory ensemble, nearby states, alternative agonist arm, null genotype or additional calcium value is simulated.

The confirmation rule was written to frozen_manifest.json before the pilot. A pass, a minimum minute flow at least 90% of 9/3088.15386, or inconclusive numerical status triggers production Radau across all ten roots at that amplitude. A clear failure at 0.25 adds no confirmation; a clear failure at 0.50 adds only its best flowing root. Best means greatest minimum minute flow, then mean flow, then ascending root ID. BDF is restricted to the decisive production representative, with a second representative only for heterogeneous production classifications. The actual decision and selected requests are saved in confirmation_decisions.json.

The calcium domain is motivated by the parotid context recorded in AGENTS.md and the original prompt: [Vera Sigüenza et al. (2019)](https://doi.org/10.1007/s11538-018-0534-z) and [Vera Sigüenza et al. (2020)](https://doi.org/10.1007/s11538-020-00712-3). These are not matched SMG calcium measurements and their per cell flow values are not calibration targets for this screen.

The new runner never opens a target ledger, an opaque target hash path, the genotype evaluation driver or archive content. The held out AE4 phenotype remains sealed. The existing manuscript and archive are unchanged.

The single independent audit checked root and parameter hashes, raw saved initial states and time grids, minute interpolation, trapezoidal totals, interval algebra, exact trajectory inventory and solver agreement. It confirmed the absolute deficit and the five versus five flow shape split. Review corrected the non calcium hash metadata to exclude the experimental calcium input and required both numerical gates in the BDF comparison. These changes required no new ODE integrations.

Focused execution tests replay genuine saved solver outputs through the complete stateless worker in serial and process modes, and recompute all twenty saved trajectory summaries in both modes. This tests execution and reduction identity, not a repeated independent ODE integration. No extra trajectory is charged to this verification. Test results and file hashes are saved separately.
