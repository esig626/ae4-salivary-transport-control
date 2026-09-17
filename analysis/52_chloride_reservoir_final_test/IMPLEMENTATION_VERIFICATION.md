# Task52 implementation verification — 52B

**PASS: 16 fixed-state checks; no production trajectory or stationary solve.**

The candidate uses the exact 19-file Task37 source/input snapshot in an isolated
package. Runtime origins and byte hashes pass; no Palk, Task41 recruitment,
Task50 coupling or later equal-routing module is imported. Historical source
bytes are unchanged. Static independent review found no blocking defect.

The deterministic charge/osmotic projectors retain the measured genotype
chloride and pH at onset. All four projected states pass the inherited
physiological, speciation and conservation checks.

| Projection | Genotype | Cl (mM) | pH | K (mM) | TIC (mM) | Volume (pL) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| central | WT | 50.10000000 | 6.91000000 | 116.02297004 | 16.59734405 | 1.4484260832 |
| central | KO | 36.50000000 | 6.89000000 | 114.83072231 | 31.38959178 | 1.4484260832 |
| alternate | WT | 50.10000000 | 6.91000000 | 114.13698474 | 5.67174623 | 1.2932920128 |
| alternate | KO | 36.50000000 | 6.89000000 | 110.43778403 | 5.67174623 | 1.1298857751 |

The central projection keeps Na=11.636125748680639 mM and
V=1.4484260832245268 pL for both genotypes. Bath osmolarity is
301.1615928036609 mOsm. Largest analytical charge residual is
1.4211e-14 fmol and largest osmotic residual is 5.6844e-14 mOsm.
The changed K/TIC/volume coordinates are explicit MODEL IDEALISATIONS;
they are not measurements or stationary transport predictions. In particular,
the central projection introduces a substantial TIC difference. No relaxation,
new chronic-state solve or secretion fit was applied.

Three original-state WT comparisons give bit-identical parent RHS, water and
AE4 sources when the auxiliary conductance is zero. Exact-null AE4 sources
vanish. The central and x1.05/x1.10 KO supplies equal the prescribed WT cycle
values exactly; independent source-vector checks pass and local NBC/NHE1,
water and electrical responses remain unchanged at each fixed state.
Unlicensed conductance/supply arguments are rejected.

All three declared auxiliary scales use the same current law and parameter
across genotypes. Their total/native/auxiliary split closes. Beta-zero dispatch
and CCh-only paired RHS nest exactly. IPR-only auxiliary current enters the
zero-NBC-recruitment closure. Two explicitly synthetic fixed-state fixtures
check opposite current directions and cell/lumen chloride cancellation; they
are not physiological onset states or additional production cases.

All inherited eleven conservation tolerances are retained. The NBC physical
charge/carbon source fields remain -B/+2B. Independent conserved cell-source
assembly and osmole accounting pass. Full evidence, software versions and
counter scope are in output/verification_attempt_01.json; projected vectors,
recovered pH and onset diagnostics are in output/projected_onsets.json.

There was one successful verification attempt: 36 top-level local evaluations,
9 paired evaluations and 19 separately counted historical constructor core
evaluations. Nested wrapper recomputations are not labelled as trajectories.
No scientific gate/tolerance, transporter parameter or projection was adjusted
after these outputs. Production remains prohibited until verified52C.
