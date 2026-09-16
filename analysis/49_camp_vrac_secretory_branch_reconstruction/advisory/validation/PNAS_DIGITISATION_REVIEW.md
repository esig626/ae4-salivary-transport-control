# Independent PNAS digitisation review

Advisory review of `/workspace/scratch/957a06f68b6a/advisory/pnas`, 2026-09-16. No model evaluation, parameter inference, canonical edit, commit or push was performed.

**Recommendation:** the corrected version 2 package is ready for Task 49C evidence provenance, subject to the explicit observation-map limitations. It does not identify a calibrated apical conductance.

## Findings and correction

The first extraction clipped the Fig. 5C marker envelope because its narrow vertical search corridor intersected the visible upper trace. For example, the original 200 s interval retained pixel rows 129–144 although the same-strip trace extends to row 123; at 300 s it retained 95–106 although the trace extends to 89. The PNAS worker preserved the original derived files under `digitisation_v1/`, expanded the extraction regions, and recorded all changed envelopes in `digitisation_correction_v2.json`.

The corrected extraction was independently recomputed from the original JPEG. Every saved version 2 pixel bound and converted midpoint/lower/upper value agrees with the image/axis calculation. No retained dark support touches an extraction-region boundary. The revised overlay covers the upper trace. Early overlapping symbols remain explicitly flagged, and the raster envelopes are not labelled biological SEM.

The misleading ledger key `current_clamp_mV` was corrected to `holding_potential_mV`; the unchanged value is -60 mV. The recording is a voltage-clamp experiment.

## Checks passed

- All 35 entries in the regenerated manifest match their SHA256 hashes.
- Every digitised non-volume mean and upper/lower SEM reproduces exactly from its saved pixels and affine axis mapping. Those CSV values are byte-identical to version 1.
- Original Fig. 4 and Fig. 5 JPEGs remain byte-identical to the preserved version 1 sources.
- Visual inspection found no control/knockout swap, axis-unit confusion, or confusion between Fig. 2 capacitance-normalised CaCC currents and Fig. 5 absolute current changes.
- Fig. 5E grey bars are correctly recorded as the magnitudes blocked by DCPIB, not currents remaining after blockade. The IPR-induced and blocked entries are respectively 208.59375 ± 35.15625 pA and 107.8125 ± 17.1875 pA from the saved picks; source SEM and the separate ±2.34375 pA point-pick interval are not conflated.

## Limits that remain necessary

The many printed digits are affine-conversion output, not measurement precision. Source SEM, pixel uncertainty and declared axis uncertainty remain separate. Unresolved SEM is null. The apparent zero NPPB trace and near-zero low-Cl bar must not become exact-zero biological constraints.

The Fig. 5C endpoint/textual swelling discrepancy is preserved; do not resolve it by retuning. Separate current and volume experiments do not supply a paired current-per-swelling coefficient. The missing verified SI details remain an access limitation, not evidence that those details are absent from the paper. Partial DCPIB blockade and possible DCPIB-insensitive current preclude treating the total IPR current as known VRAC current. No whole-cell-to-apical-model conductance map is established by this metrology review.
