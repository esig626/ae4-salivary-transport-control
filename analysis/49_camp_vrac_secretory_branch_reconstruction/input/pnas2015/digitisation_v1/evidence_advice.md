# Advisory evidence extraction — Catalán PNAS 2015

Source: [Catalán et al., DOI 10.1073/pnas.1415739112](https://pmc.ncbi.nlm.nih.gov/articles/PMC4343136/). Primary HTML and PMC-hosted figures were retrieved. Supplement access returned a download challenge rather than a PDF; publisher and Europe PMC routes were unsuccessful. No supplementary protocol quantity is invented.

## Primary observations

IPR was 5 μM; gland collection lasted 10 min. Figure 4 used six control and six Tmem16A-null glands (total-volume comparison p=0.092). Figure 5 tested DCPIB 50 μM (n=8) and NPPB 100 μM (n=6), with control and Clcn2-null groups n=8. The source reports IPR swelling 12.5±0.2% and CCh shrinkage 10±1%. Figure 5E gives absolute current changes at −60 mV in pA, with IPR n=4 and hypotonic n=3. Grey bars are DCPIB-blocked current magnitudes, not currents remaining after blockade. Current is outward rectifying; DCPIB only partially blocks it. Other current contributions remain possible. Low extracellular chloride reduced secretion by 98%. CFTR and ClC-2 loss did not significantly affect IPR secretion. VRAC identity is physiological/pharmacological inference, not demonstrated molecular LRRC8 identity.

## Reproducible extraction

Run `digitise_pnas.py` with Pillow/NumPy. It reads original JPEGs, saves all axis calibration and raw pixel selections, and emits CSVs plus visual overlays. Point uncertainty is ±1.5 original pixels; axis-anchor uncertainty is separately declared ±1 pixel. Published SEM caps are digitised separately. Hidden or inseparable SEM is blank, never assigned zero. Figure 5C returns a marker/trace envelope, not a biological SEM.

- `digitised_points.csv`: Figure 4A/B and Figure 5A/B/E.
- `digitised_volume_trace.csv`: Figure 5C IPR trace envelope at fixed 25-s locations.
- `digitisation_metadata.json`: complete calibration, source checksums, raw selections and limitations.
- `figure4_picks.png`, `figure5_picks.png`: red mean/trace selections; blue SEM-cap selections.

Source text and Figure 5C are not numerically identical: the final raster trace lies near 1.16, while the article narrative states 12.5% swelling. This is preserved as an unresolved source inconsistency. Do not choose one silently or use the difference to infer a transport parameter.

## Advice on the one-scale map

For the predeclared branch

`I_apical = g * beta * max(V/Vrest - 1, 0) * (Va - ECl)`,

a single paired observation could identify g only if the blocked apical current, beta normalisation, swelling at that same instant, driving force and correspondence between the sampled cell and model cell were fixed independently. Current and swelling here come from separate experiments. The primary main text does not give the complete matched state or identify all blocked whole-cell current as apical current. No verified SI bath/pipette chloride, temperature, hypotonicity, capacitance/geometry or numerical reversal-potential value has been extracted in this run. A missing verified quantity is an extraction limitation, not proof that the original supplement lacks it.

The digitised blocked IPR current is approximately 108 pA. Dividing this by −60 mV would be invalid because −60 mV is clamp voltage, not chloride driving force. Dividing the total IPR current (~209 pA) by the swelling fraction would also mix current identities and unmatched observations. The figure contains a nonzero independent signal but does not by itself certify an absolute model apical conductance.

Whole-gland flow is in μL/min, not pL/s per model acinar cell. No gland-to-model-cell scale is provided by these figure axes. Controls in Figure 4 and Figure 5 have materially different total outputs, so do not combine cohorts as if they shared a fixed physical scale. A nonsignificant Tmem16A genotype comparison does not establish exact equality.

Any inference from a zero model swelling gate is algebraic: if the inherited trajectory remains at or below its own resting volume, the new current is exactly zero for every finite g. Therefore source current/flow cannot identify g on that protocol, and selecting an arbitrary positive g would not be a source calibration. This is advice for orchestrator verification, not an accepted scientific conclusion.

## Restrictions observed

No AE4 held-out output or uptake was used. No canonical repository files, model trajectories, fit, commit or push were changed or executed. The script performs image metrology only.
