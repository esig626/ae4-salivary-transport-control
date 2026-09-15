"""Generate the article narrative from fixed, verified outputs only."""
from pathlib import Path
import json

HERE=Path(__file__).resolve().parent

def main():
    summaries={case:json.loads((HERE/'output'/case/'summary.json').read_text()) for case in ['wt','ae4_null','ae4_5pct']}
    validation=json.loads((HERE/'output/final_validation.json').read_text())
    attribution=json.loads((HERE/'output/baseline_interpretation.json').read_text())
    verification=json.loads((HERE/'output/final_verification.json').read_text())
    w,k=summaries['wt']['endpoint'],summaries['ae4_null']['endpoint']
    metrics={r['case']:r for r in validation['results']}
    source_base='https://github.com/esig626/ae4-salivary-transport-control'
    prediction=validation['prediction_commit']
    predictions='\n'.join(f"| {'WT' if case=='wt' else ('AE4 null' if case=='ae4_null' else '5% AE4')} | {metrics[case]['cumulative_outflow_pL']:.9f} | {metrics[case]['cumulative_deficit_percent']:.4f}% | {metrics[case]['nkcc_compensation_percent']:.4f}% |" for case in ['wt','ae4_null','ae4_5pct'])
    physiology='\n'.join(f"| {label} | {w[key]*scale:.6f} | {k[key]*scale:.6f} | {summaries['ae4_5pct']['endpoint'][key]*scale:.6f} |" for label,key,scale in [('Na (mM)','na_i_mM',1),('K (mM)','k_i_mM',1),('Cl (mM)','cl_i_mM',1),('pH','ph_i',1),('Cell volume (pL)','volume_i_pL',1),('Apical voltage (mV)','v_apical_V',1000),('Basolateral voltage (mV)','v_basolateral_V',1000)])
    text=rf'''# Dynamic potassium recycling after AE4 loss

**Conclusion.** The inherited dynamic system does not impose a potassium recycling ceiling that prevents NKCC1 compensation. It carries the extra K burden through increased basolateral efflux and altered ion storage, while the pump experiences little additional net sodium demand. No independently specified correction to the recycling block was identified. This is a negative result within the frozen model; the available evidence does not establish that its effective channel and pump capacities are correct in real acini.

## 1. Problem

The conservation explicit reconstruction preserves primary secretion after AE4 deletion because NKCC1 replaces most of the missing chloride loading. The scientific question is whether that replacement has a hidden physiological cost: NKCC1 imports K with Na and Cl, and the cell must store or export that K while maintaining voltage and pump function. A substantial potassium burden could therefore exist even when chloride transport remains stoichiometrically feasible.

Task 47 tests this one explanation. It retains the exact Task 40 WT parameterisation and initial state, uses the completed Task 46 H00 to H04 results, and performs no CBM rerun, fitting, cap selection, mechanism search or parameter sweep. The precomparison decision was to make no scientific parameter or equation change.

## 2. Structural result

The completed CBM shows that AE4 deletion alone does not exclude full reference chloride export. Its feasible flux rearrangement survives the independent NKCC assay evidence. The isolated initial uptake assay cannot be converted into an exact stimulated WT/null NKCC equality or an invented tolerance. A cap chosen to remove the full rescue solution would add a new assumption rather than follow from those data.

The CBM nevertheless provides an exact mechanistic lead. With J the apical chloride export, N the NKCC cycle flux, H the NHE1 flux, A the AE4 flux and K the sum of both K channel effluxes, its zero storage balances give

$$K=J/3+N+H/3-A/2.$$

Deleting AE4 removes the negative A/2 term; recruiting an additional NKCC cycle adds one K to the recycling requirement. This establishes a transport cost. It does not supply a measured upper limit on K recycling.

## 3. Mechanistic lead in the dynamic model

The dynamic amount balances retain finite storage. Their exact counterpart is

$$K=J/3+N+H/3-A/2+S,$$
$$S=\frac{{-2\dot n_{{Na}}+\dot n_{{TA}}+\dot n_{{Cl}}}}{{3}}-\dot n_K
  =-\frac{{\dot n_{{Na}}+2\dot n_K}}{{3}},$$

where the last equality uses charge conservation. K concentration also depends on changing cell volume. The steady identity was never imposed on the ODE; both forms of the dynamic identity were checked against its independently assembled fluxes.

The active K currents already use finite calcium gated conductances and Nernst reversal potentials. The pump already has a finite capacity and Na and extracellular K saturation. NKCC already responds to intracellular K through `X = [Na_i][K_i][Cl_i]^2`, with a strictly decreasing flux law in X. There is no missing potassium conservation term, pump stoichiometric coefficient or passive channel reversal in this block.

Voltage is algebraic, obtained from current balance as amounts and volumes evolve. Including NBC gives the further exact instantaneous relation `K + P + B = J`, where P is total pump turnover and B is the inward NBC cycle flux. NBC must therefore remain in the causal account of the electrical response.

## 4. Dynamic test

One WT and one null trajectory were run for the inherited 600 s stimulus, with the same saved WT rest and six identical parameter hashes. The analysis records Na and K amounts and concentrations, volume dilution, apical and basal K currents and reversal potentials, pump utilisation, NKCC driving terms, chloride export and water output. It retains all inherited physiological gates. Exact local electrical derivatives were checked at three saved times with two numerical difference steps; no altered parameter trajectory was run.

![Dynamic potassium recycling diagnostics](figures/dynamic_recycling.png)

*Figure 1. The frozen model increases NKCC recruitment and K efflux after AE4 loss, with modest K concentration and driving voltage changes and no pump saturation. The kinetic limit is internal to the adopted pump law, not an experimentally measured capacity. Panel E starts at the stimulated right limit. The saved data retain the exact resting onset sample.*

During 60 to 600 s, total K channel efflux rises from 135.495501 to 159.975608 fmol, an increase of **{attribution['k_channel_increase_percent']:.4f}%**. Most of the increase is basal: basal efflux rises from 100.151507 to 126.349120 fmol, while apical efflux falls slightly. The fixed conductance fractions therefore do not fix realised efflux fractions.

At 600 s, basal K driving voltage increases from {1000*w['k_drive_basolateral_V']:.6f} to {1000*k['k_drive_basolateral_V']:.6f} mV. A {attribution['endpoint_basal_voltage_change_mV']:.6f} mV depolarisation and a {abs(attribution['endpoint_basal_reversal_change_mV']):.6f} mV more negative K reversal jointly supply the extra **{attribution['endpoint_extra_basal_drive_mV']:.6f} mV**. The stimulated K conductance is unchanged at {(w['g_k_apical_S']+w['g_k_basolateral_S'])*1e9:.6f} nS in both cases. The model accommodates the extra current by changing its driving force, without additional channel recruitment.

An Ohmic current divided by its instantaneous active conductance times driving voltage is identically one. Dividing instead by the fully open conductance gives the calcium gate, {w['k_calcium_gate']:.6f}. Neither ratio independently measures spare physiological capacity. A global K flux ceiling would require additional bounds on voltage, concentration and channel behaviour; the numerical voltage bracket cannot supply those bounds.

## 5. Result and complete balance of the compensation

The proposed sequence does not reach a limiting step in the retained model. K feedback opposes compensation, but does not dominate the chloride response. Endpoint K is higher in the null by {k['k_i_mM']-w['k_i_mM']:.6f} mM, while Cl is lower by {w['cl_i_mM']-k['cl_i_mM']:.6f} mM. The endpoint changes in log X are {attribution['endpoint_log_product_attribution']['k']:.6f} from K, {attribution['endpoint_log_product_attribution']['na']:.6f} from Na and {attribution['endpoint_log_product_attribution']['cl']:.6f} from twice log Cl. Their sum is negative. The declining chloride concentration therefore outweighs K feedback in the NKCC law. This decomposition is an algebraic attribution of the observed concentration product, not an independent intervention on each ion.

![Na and K balance attribution](figures/ion_balance_attribution.png)

*Figure 2. Differences in integrated null and WT balances over 60 to 600 s. Each panel sums to zero. The storage term is moved to the removal side for bookkeeping; it denotes retention rather than export.*

The complete K account is particularly informative. Extra NKCC imports 18.168504 fmol K, loss of AE4 removes 20.527751 fmol K export, and the slightly increased pump adds 1.216579 fmol K influx. The additional 39.912834 fmol burden is balanced by 24.480107 fmol extra channel efflux and 15.432727 fmol more K storage relative to WT. Actual null storage in this window is only 3.288542 fmol; much of the relative storage difference arises because WT loses K. Treating the zero storage CBM identity as a transient constraint would miss this distinction.

The Na account explains why the pump is not exhausted. Extra NKCC Na entry and removal of AE4 Na export add the same 38.696255 fmol burden. NBC Na entry falls by 33.716346 fmol and NHE1 entry by 2.294385 fmol. Net sodium entry before the pump therefore rises by only 2.685524 fmol. Pump Na clearance increases by 1.824869 fmol and Na storage by 0.860655 fmol. Total pump turnover rises by just **{attribution['pump_increase_percent']:.4f}%**, despite **{attribution['nkcc_increase_percent']:.4f}%** greater NKCC turnover.

At 600 s, pump use is {100*w['pump_nominal_utilisation']:.4f}% of nominal capacity in WT and {100*k['pump_nominal_utilisation']:.4f}% in the null. Relative to the more relevant limit conditioned on the instantaneous extracellular K values, the fractions are {100*w['pump_k_conditioned_utilisation']:.4f}% and {100*k['pump_k_conditioned_utilisation']:.4f}%. The null never exceeds 84.3251% over the monitored trajectory. The pump thus retains kinetic headroom in this parameterisation because it receives little extra net Na burden.

These results reject an already binding K recycling ceiling as the explanation for failed rescue in this frozen architecture. They do **not** demonstrate that its reserve is biologically correct. Task 43 identifies the effective conductance, localisation fractions, gating and pump capacity as uncertain model settings. Apical BK and IK currents are supported experimentally in parotid acini, but regional channel density does not determine the present whole cell conductance partition or a matched sustained capacity: [Almassy et al. (2012)](https://doi.org/10.1085/jgp.201110718).

The audit also prevents a stronger claim of complete physiological validation. The pump has no direct voltage or ATP dependence, and the K law lumps BK and IK gating. Consequently the full proposed voltage to pump feedback is represented only indirectly through ionic state changes. The Palk coefficients also imply a fixed bath NKCC reversal product that differs from the active bath. Both selected and ideal bath affinities remain positive in every saved state, but the law is not globally validated against bath thermodynamics. Correct reversal alone does not identify its kinetics. None of these limitations supplies a unique, independently constrained numerical repair. Reducing a conductance or pump capacity to obtain the desired phenotype would therefore be an unsupported construction, and was not done.

## 6. Phenotype validation after freezing

The prediction checkpoint [{prediction[:7]}]({source_base}/commit/{prediction}) was published and its remote SHA verified before the deferred 5% trajectory and numerical phenotype comparison. All {verification['frozen_files_unchanged']} frozen files remain unchanged. No retuning followed comparison.

| Condition | Cumulative water output, 0 to 600 s (pL) | Cumulative deficit | NKCC increase, 60 to 600 s |
| --- | ---: | ---: | ---: |
{predictions}

The direction is correct, but the null deficit is only {metrics['ae4_null']['cumulative_deficit_percent']:.4f}%, leaving a descriptive gap of {validation['null_magnitude_gap_percentage_points']:.4f} percentage points from the central approximately 35% whole gland benchmark reported by [Peña-Münzenmayer et al. (2015)](https://doi.org/10.1074/jbc.M114.612895). The 5% condition is a computational expression intervention; no separate experimental magnitude is claimed for it. These comparisons do not equate an acute deletion from a shared WT state with an established knockout animal or introduce a fitted whole gland observation model.

The null deficit is {metrics['ae4_null']['early_0_60_deficit_percent']:.4f}% in the early 0 to 60 s window and {metrics['ae4_null']['late_60_600_deficit_percent']:.4f}% in the 60 to 600 s window. Thus the later window has the greater loss, but neither its magnitude nor the cumulative loss approaches the benchmark. NKCC still replaces {100*attribution['replaced_ae4_chloride_fraction']:.4f}% of the missing integrated AE4 chloride load. Compensation has not been reduced by a new recycling mechanism.

| Intracellular or electrical readout at 600 s | WT | AE4 null | 5% AE4 |
| --- | ---: | ---: | ---: |
{physiology}

All monitored Na, K, Cl, pH and volume values satisfy the inherited bounds; tracked amounts remain positive and current, charge, carbon and water checks pass. WT remains acceptable under its original activation checks. Voltages remain finite and negative and channel fluxes follow their driving forces, but the inherited protocol does not supply a separately measured voltage acceptance interval. Passing those model gates is not proof of every omitted kinetic dependency or of behaviour beyond 600 s.

The article implication is therefore specific. Conservation alone permits rescue, and the existing finite recycling dynamics do not remove it. The additional K burden is real, but it is handled through basal efflux, storage and redistribution of sodium entry. The present evidence neither identifies a physiological recycling correction nor establishes which different coupling limits secretion in the gland. Measurements of K current versus voltage and calcium, together with pump turnover, Na, K and volume during the same stimulation protocol, would distinguish excessive model reserve from another missing mechanism. Task 47 stops at this result without opening another mechanism search.

## Reproducibility and provenance

The equation audit, complete parameter provenance, baseline diagnostics, freeze manifest, remote receipt and machine readable validation are retained in this directory. Three production integrations were performed, with no resting solve, CBM solve, optimisation, scientific numerical retry or altered parameter trajectory. The independent archive verifier checks 1,806 saved states without importing the production model; 24 local electrical derivative checks have maximum scaled error {verification['max_scaled_local_derivative_error']:.3g}. All 1,882 files in the starting repository tree are preserved by the publication checks.

One WT summary JSON serialization failure was repaired from its already saved CSVs without repeating the integration. The WT summary therefore reports saved grid extrema and does not invent lost solver counters; the original accepted endpoint gate checks completed before saving. This limitation is documented in `output/EXECUTION_NOTES.md`. Figures are generated only from saved predictions and were visually inspected. The production source, manuscript and previous numbered analyses were not changed.
'''
    (HERE/'DYNAMIC_POTASSIUM_RECYCLING_REPORT.md').write_text(text)
    print('Article report written from verified fixed results.')

if __name__=='__main__':main()
