# Frozen baseline potassium diagnostics

Exactly one WT and one AE4 null trajectory used the identical frozen Task 40 parameters and WT rest. All inherited gates passed. There was no CBM solve, resting solve, optimisation or scientific parameter change. The WT metadata recovery is documented in `output/EXECUTION_NOTES.md`; its stored extrema are for the saved grid. The null has 333 accepted steps and 934 checked states.

## Integrated transport, 60 to 600 s

| Quantity (fmol) | WT | AE4 null |
| --- | ---: | ---: |
| K channel efflux | 135.495501 | 159.975608 |
| Apical K efflux | 35.343994 | 33.626488 |
| Basolateral K efflux | 100.151507 | 126.349120 |
| NKCC K entry | 78.436089 | 96.604592 |
| Pump cycles | 32.721489 | 33.329779 |
| NBC cycles | 39.110273 | 5.393927 |
| NHE1 cycles | 4.056565 | 1.762180 |
| Na net entry before pump | 101.075175 | 103.760699 |
| Na pump clearance | 98.164467 | 99.989336 |
| K amount storage | -12.144185 | 3.288542 |

Total K channel efflux increases by 18.0671%, whereas pump cycles increase by 1.8590%. Increased NKCC recruitment and removal of the AE4 K exit add 38.696255 fmol to the relative K burden. Extra pump influx adds 1.216579 fmol. Channels carry 24.480107 fmol more and the difference in K storage is 15.432727 fmol. The latter is a difference between the two trajectories, not the absolute null storage.

The Na burden does not grow in proportion to NKCC alone. NBC influx falls by 33.716346 fmol and NHE1 by 2.294385 fmol, partly offsetting the extra NKCC Na entry and removal of AE4 Na export. Net Na entry before the pump changes by only 2.685524 fmol. These are the complete Na and K balances, including changing storage.

## Endpoint mechanism

| Quantity at 600 s | WT | AE4 null |
| --- | ---: | ---: |
| Intracellular K (mM) | 110.925054 | 113.531289 |
| Intracellular Na (mM) | 17.554236 | 17.461719 |
| Intracellular Cl (mM) | 53.979919 | 51.784216 |
| Basolateral voltage (mV) | -79.203783 | -78.602441 |
| Basolateral K reversal (mV) | -82.837058 | -83.457750 |
| Basolateral K drive (mV) | 3.633276 | 4.855309 |
| Apical K drive (mV) | 2.996502 | 2.939599 |
| Total K channel efflux (fmol/s) | 0.242586 | 0.301668 |
| Pump / nominal capacity | 0.779394 | 0.777453 |
| Pump / external K conditioned limit | 0.843978 | 0.841879 |
| Dynamic storage correction to CBM (fmol/s) | 0.004296 | -0.003799 |

The basal voltage changes by 0.601342 mV and the K reversal changes by -0.620692 mV. Together these increase the K driving voltage by 1.222033 mV. Both genotypes have exactly the same stimulated active K conductance: 6.799637 nS, split into 2.039891 apical and 4.759746 basal nS. No channel recruitment change is needed for the extra current in this model.

K feedback on NKCC exists and opposes compensation. At 600 s its contribution to log X is 0.023224; Na contributes -0.005284 and twice log Cl contributes -0.083053. The total -0.065114 is negative, so the chloride fall dominates K accumulation and N increases through the unchanged decreasing concentration law. This is an exact endpoint concentration product decomposition, not an independent intervention proving the origin of each state change.

The channel active law ratio is one by construction. The ratio to the fully open current is 0.485688; it is fixed by the unchanged calcium input and is not free additional reserve. The maximum K conditioned pump utilisation is 0.843978 in the WT saved samples and 0.843250 over all monitored null states. The pump does not approach its kinetic endpoint because the extra net sodium burden is small.

Both ideal bath and selected NKCC affinities remain positive at every sampled state. The reversal products nevertheless differ: 7837380.573248 versus 11539641.936599 mM^4. The fixed bath law is not globally thermodynamically validated. No observed flux has the wrong sign under the actual bath in these runs.

Exact voltage closure derivatives were checked for K conductance and pump capacity at 60, 300 and 600 s in both genotypes, using two finite difference steps. These 24 checks are local state conditioned calculations, not alternate trajectories or a parameter sweep. They verify electrical coupling but do not identify physiological parameter values. Pump voltage dependence remains absent from the inherited law.

## Decision before phenotype comparison

The existing finite recycling block accommodates the extra burden without pump saturation or a sampled physiology failure. This rejects an already binding K recycling ceiling as the explanation in the frozen architecture. It does not show that its unmeasured effective conductance is correct in real acini, nor rule out a differently constrained biological recycling system.

No uniquely supported correction to K conductance, gating, distribution or pump capacity has been identified. A missing detailed pump voltage law or mixed BK/IK gating requires kinetic and abundance information, not an invented scalar reduction. The fixed bath NKCC discrepancy is recorded explicitly; all realised signs agree, and a bath dependent kinetic replacement is not uniquely specified by its reversal alone. No numerical correction is made. The raw WT and null prediction outputs are ready to be frozen before the phenotype comparison; the 5% case has not been evaluated.
