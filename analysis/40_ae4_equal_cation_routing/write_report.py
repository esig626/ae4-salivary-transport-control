"""Render the final report from saved results only; no scientific execution."""
from summarize_results import ROOT, OUT, HERE, CASES, read_json, read_csv


def table(headers, rows):
    return '\n'.join(['| '+' | '.join(headers)+' |',
                      '| '+' | '.join('---' for _ in headers)+' |']+
                     ['| '+' | '.join(map(str,row))+' |' for row in rows])


def main():
    summary=read_json(OUT/'comparison_summary.json')
    rest=read_json(OUT/'wt_rest.json');budget=summary['budget']
    verification={c:read_json(OUT/(c+'_verification.json')) for c in CASES}
    series={c:[{k:float(v) for k,v in row.items()} for row in read_csv(OUT/(c+'_timeseries.csv'))] for c in CASES}
    endpoint={c:series[c][-1] for c in CASES}
    p=summary['perturbations'];five,null=p['ae4_5pct'],p['ae4_null']
    parts=[f"Task 40 — fixed 50:50 AE4 cation routing\n\n**{summary['status']}**\n",
    "The single WT resting-state solve and all three production trajectories passed. Equal routing reduced AE4-loss NKCC1 compensation, but did not preserve higher intracellular Na concentration at 600 s or produce a large secretion deficit. This is a completed structural test with a partial compensation reduction; the proposed sustained-Na outcome is unsupported in this frozen architecture.",
    table(['Measure','AE4 = 0.05','AE4 = 0.00'],[
        ['600 s Na_i minus WT (mM)',f"{five['endpoint_differences']['na_i_mM']:+.6f}",f"{null['endpoint_differences']['na_i_mM']:+.6f}"],
        ['60–600 s NKCC1 compensation',f"+{five['nkcc1_compensation_percent']:.5f}%",f"+{null['nkcc1_compensation_percent']:.5f}%"],
        ['Task 39 compensation',f"+{five['task39_nkcc1_compensation_percent']:.5f}%",f"+{null['task39_nkcc1_compensation_percent']:.5f}%"],
        ['Compensation change (percentage points)',f"{five['compensation_change_percentage_points']:+.5f}",f"{null['compensation_change_percentage_points']:+.5f}"],
        ['600 s Cl_i minus WT (mM)',f"{five['endpoint_differences']['cl_i_mM']:+.6f}",f"{null['endpoint_differences']['cl_i_mM']:+.6f}"],
        ['Cumulative secretion / Task 40 WT',f"{five['cumulative_ratio']:.10f}",f"{null['cumulative_ratio']:.10f}"],
        ['Cumulative secretion deficit',f"{five['cumulative_deficit_percent']:.6f}%",f"{null['cumulative_deficit_percent']:.6f}%"],
        ['Task 39 cumulative deficit',f"{five['task39_cumulative_deficit_percent']:.6f}%",f"{null['task39_cumulative_deficit_percent']:.6f}%"]]),
    "The binding inputs were read in order: [AGENTS.md](../../AGENTS.md), [source contract](source_contract.md), then [Task 40 prompt](../../prompts/40_ae4_equal_cation_routing.md). The prepared head was `a6a58ebabcd2aedb1702b52dfdd86359c60cec32`, on `codex/task-40-ae4-equal-cation-routing`. The scientific parent was merged Task 39, `afd101448763439f369f2682c467f069ea18a442`. The preparation changed only these three instruction files relative to that parent. The downloaded scientific inputs were verified against the prepared Git tree; the manifest contains 150 unchanged files. The old donor-weighted evaluator remains selectable.",
    "The new named evaluator calls the same inherited QSS evaluator for J4 and then returns the fixed source vector. There is no routing-fraction parameter or donor-dependent choice. For each evaluated state, expression and regulation, the J4 law is identical to Task 39. J4 can differ between evolved trajectories because their states differ; it was never clamped to its reference value.",
    "`(S_Na, S_K, S_Cl, S_TIC, S_TA) = (-J4/2, -J4/2, J4, -2J4, -2J4)`. Thus `S_Na + S_K = -J4`, and `S_Na + S_K - S_Cl - S_TA = -J4/2 - J4/2 - J4 + 2J4 = 0`. Negative J4 reverses every source automatically. The focused tests verify actual forward and reversed environments, zero cycle, exact inherited J4 and anion sources, expression scaling, and unchanged nonrouting diagnostics. This is the stipulated ensemble source allocation, not a separately derived microscopic equal-routing carrier.",
    "Independent 50-digit decimal arithmetic and full-model evaluation reproduced the preflight:",
    table(['Reference quantity','Reproduced value'],[
        ['Na_i / (Na_i + K_i)','0.09062450161183984'],
        ['Inherited J4 (fmol/s)','0.004937067441546469'],
        ['Expected delta S_Na (fmol/s)','−0.0020211144444590447'],
        ['Expected delta S_K (fmol/s)','+0.0020211144444590447'],
        ['Actual full-RHS Na/K shifts (fmol/s)','−0.0020211144444590434 / +0.0020211144444590434'],
        ['Direct Cl/TIC/TA shifts','Exactly zero'],
        ['Charge perturbation (fmol/s)','−4.34e−19 (round-off)']]),
    "All five focused tests passed after one documented literal correction to a floating-point scaling assertion. The original error was 2.0817e−17 fmol/s against a 2e−17 allowance; a scale-aware 32-machine-epsilon comparison replaced that assertion. Exact identity to inherited J4 remained required. No evaluator, parameter, production tolerance or scientific gate was changed. Both logs and [the correction record](source_correction.md) are retained. [Source verification](../../results/40_ae4_equal_cation_routing/source_verification.json) records the algebra and parameter hashes.",
    "Exact Task 31 REST nesting was not expected and was not imposed. The WT rest used only the accepted Task 31 R09 state. A single square ten-equation `root/hybr` solve used the inherited charge manifold, residual scales, coordinate bounds and acceptance tolerances. A logit coordinate transform enforced those numerical bounds. No least-squares fit, calibration, multistart, continuation, other root family, or genotype rest was used. There were 86 root residual evaluations and 12 final audit RHS evaluations, with no retry. The sole initial residual used the literal accepted state; its coordinate reconstruction differs by less than 1e−12.",
    "The accepted core state and its shift are below. Amounts are fmol and volumes pL. The sole regulatory coordinate is zero at rest. The frozen full-state SHA-256 is `"+rest['state_sha256']+"`.",
    table(['Core state','Task 31','Task 40','Change'],[
        [row['state'],f"{float(row['task31']):.12g}",f"{float(row['task40']):.12g}",f"{float(row['task40_minus_task31']):+.8g}"]
        for row in read_csv(OUT/'wt_rest_core_comparison.csv')]),
    table(['Rest observable','Task 31','Task 40','Change'],[
        [key,f"{rest['observables'][key]-rest['observed_shift_from_task31'][key]:.10g}",
         f"{rest['observables'][key]:.10g}",f"{rest['observed_shift_from_task31'][key]:+.7g}"]
        for key in ('na_i_mM','k_i_mM','cl_i_mM','ph_i','hco3_i_mM','tic_i_mM','volume_i_pL',
                    'q_out_pL_s','v_apical_V','v_basolateral_V','v_transepithelial_V','osm_lumen_mOsm')]),
    "The complete rest flux ledger, water fluxes, voltages, acid/base speciation, and current/charge/carbon residuals are in [wt_rest.json](../../results/40_ae4_equal_cation_routing/wt_rest.json); every ledger shift is in [wt_rest_flux_comparison.csv](../../results/40_ae4_equal_cation_routing/wt_rest_flux_comparison.csv). Rest J4 is 0.004933930166208064 fmol/s; each cation source is −0.002466965083104032 fmol/s. Rest NBC cycle and current are zero. The small change from the old J4 is a state consequence, not a rate-law change.",
    f"Rest passed finite positive tracked amounts/volumes, Na_i < 40 mM, K_i 50–200 mM, Cl_i 30–80 mM, pH_i 6.6–7.3, positive TIC/HCO3 with HCO3 < 100 mM, volume < 3 pL, and all inherited production conservation tolerances. The maximum scaled stationary residual was {rest['max_abs_scaled_independent_rhs']:.3e} against 1e−7; maximum amount RHS {rest['max_abs_amount_rhs_fmol_s']:.3e} fmol/s; maximum volume RHS {rest['max_abs_volume_rhs_pL_s']:.3e} pL/s. The independent Jacobian had rank 10, no coordinate boundary was hit, the omitted charge rows and regulatory RHS passed, and the resting RHS was time invariant.",
    "The rest was frozen before the three trajectories. WT, 5%, and null have bitwise identical initial vectors and identical scientific parameter hashes. Only `genotype.ae4_expression` differs. All used the inherited 600 s standard stimulus, Radau rtol 1e−7, amount atol 1e−10 fmol, volume atol 1e−12 pL, regulatory atol 1e−10, and maximum step 2 s. The inherited onset convention keeps the exact t=0 resting-stimulus readout and begins integration at 1e−6 s from the same vector. Every accepted endpoint and integer-second dense sample was checked: 934 states per case. All gates passed, including WT secretion activation; none was rescued.",
    "All requested observables at 0, 60, …, 600 s are in [WT](../../results/40_ae4_equal_cation_routing/wt_timeseries.csv), [5%](../../results/40_ae4_equal_cation_routing/ae4_5pct_timeseries.csv), and [null](../../results/40_ae4_equal_cation_routing/ae4_null_timeseries.csv). They include q_out, cumulative secretion, intracellular composition/volume, all requested transporter fluxes, separate AE4 Na/K sources and export branches, NBC current, voltages, lumen osmolarity, and Palk X. Integrals use the inherited one-second trapezoidal grid, with onset samples; they were not recomputed from the coarser report rows.",
    table(['600 s observable','WT','5%','Null'],[
        [key,*[f"{endpoint[c][key]:.9g}" for c in CASES]] for key in
        ('na_i_mM','k_i_mM','cl_i_mM','ph_i','hco3_i_mM','tic_i_mM','volume_i_pL',
         'nkcc1_x_i_mM4','q_out_pL_s','cumulative_outflow_0_t_pL','v_apical_V','v_basolateral_V',
         'v_transepithelial_V','osm_lumen_mOsm')]),
    "Na_i is higher after AE4 loss at the reported 60–480 s samples for 5%, and 60–540 s samples for null. The largest tabulated differences occur at 180 s: +0.9186 and +1.0317 mM. The 5% difference changes sign between 480 and 540 s; null between 540 and 600 s. At 600 s, both are below WT. These intervals use the saved 60 s readouts; no finer crossing time is claimed.",
    "Positive NKCC1 chloride loading integrated over 60–600 s rises from 156.8722 fmol in WT to 189.5319 / 193.2092 fmol. Compensation is smaller than Task 39 by 8.8350 / 9.9379 percentage points, a roughly 30% relative reduction in the compensation percentage. It neither disappears nor reverses. The NKCC1 increase replaces 88.86% / 88.51% of lost AE4 chloride loading. Endpoint Palk X is lower by 5.648% / 6.304%; this is consistent with increased loading under the unchanged decreasing Palk concentration-response law, despite the transient Na elevation. These comparisons use each task's own frozen WT; Task 39 trajectories were not rerun.",
    "Cl_i remains lower at every positive-time reported checkpoint. At 600 s the deficits are 1.8943 / 2.1957 mM (3.509% / 4.068%). They are smaller than the 180 s differences of 3.1704 / 3.5016 mM, so compensation partly closes the concentration gap without removing it within 600 s. This is a modest persistent concentration deficit, not evidence of a large depletion of cellular chloride content. Larger cell volumes matter: final Cl amounts are "+
        f"{endpoint['wt']['cl_i_mM']*endpoint['wt']['volume_i_pL']:.4f} / {five['endpoint_cell_cl_amount_fmol']:.4f} / {null['endpoint_cell_cl_amount_fmol']:.4f} fmol (WT/5%/null), while Na amounts are {endpoint['wt']['na_i_mM']*endpoint['wt']['volume_i_pL']:.4f} / {five['endpoint_cell_na_amount_fmol']:.4f} / {null['endpoint_cell_na_amount_fmol']:.4f} fmol. Concentration and amount therefore have different comparisons.",
    "Integrated fluxes over 60–600 s are shown below. Percentages for positive AE2 loading are undefined relative to its zero WT value; signed AE2 changes are described by direction. NBC and NHE1 stimulation settings remain unchanged despite their flux responses.",
    table(['Integrated quantity','Unit','WT','5%','Null','5% change','Null change'],[
        [row['quantity'],row['unit'],f"{float(row['wt']):.8g}",f"{float(row['ae4_5pct']):.8g}",f"{float(row['ae4_null']):.8g}",
         ('sign reversal' if row['quantity']=='signed_ae2_cl_flux' else f"{float(row['ae4_5pct_relative_change_percent']):+.3f}%" if row['ae4_5pct_relative_change_percent'] else 'undefined'),
         ('sign reversal' if row['quantity']=='signed_ae2_cl_flux' else f"{float(row['ae4_null_relative_change_percent']):+.3f}%" if row['ae4_null_relative_change_percent'] else 'undefined')]
        for row in read_csv(OUT/'integrated_flux_comparison.csv')]),
    "AE2 switches from signed export in WT to a small positive loading route after AE4 loss. Integrated positive loading is 0 / 0.43510 / 0.51148 fmol, only 0 / 0.22397% / 0.26403% of the positive chloride-loading pool. It remains negligible in that quantitative sense, rather than exactly absent. WT's NKCC1 positive-loading share is 79.257%, outside the inherited 65–75% contextual band; that band is explicitly validation context and was not an acceptance gate.",
    "Loss of AE4 base export is accompanied by higher TIC, bicarbonate and pH. Integrated NBC cycles, bicarbonate influx and current fall by 77.683% / 86.208%; NHE1 influx falls by 53.096% / 56.560%. Integrated pump cycles rise modestly, by 1.588% / 1.859%, while their endpoint rates are slightly below WT (0.062119 / 0.062196 versus 0.062352 fmol/s). Thus a transient rise in Na concentration does not imply a sustained rise: reduced Na export through AE4 coexists with markedly reduced NBC/NHE1 Na entry, pump feedback, and volume change. This interpretation follows the recorded coupled balances and is not an additional intervention.",
    "CaCC chloride export falls by 3.687% / 4.162%, and paracellular chloride return falls by 3.335% / 3.750%. At 600 s, NBC cycle flux is 0.071918 / 0.009343 / 0.001679 fmol/s and current is 6.939 / 0.901 / 0.162 pA (WT/5%/null); NHE1 is 0.006245 / 0.001836 / 0.001567 fmol/s. CaCC export is 0.376856 / 0.367142 / 0.365543 fmol/s, and paracellular return is 0.140254 / 0.136569 / 0.136006 fmol/s. The complete endpoint fluxes are in the trajectory CSVs.",
    "Cumulative secretion deficits are `100 × (1 − Q_case(0,t)/Q_WT(0,t))`; absolute loss is `Q_WT − Q_case`. At t=0, both cumulatives are zero and their ratio/percentage deficit is undefined.",
    table(['Time (s)','WT cumulative (pL)','5% deficit','Null deficit','5% lost (pL)','Null lost (pL)'],[
        [int(float(row['time_s'])),f"{float(row['wt_cumulative_pL']):.8f}",
         f"{float(row['ae4_5pct_cumulative_deficit_percent']):.5f}%" if row['ae4_5pct_cumulative_deficit_percent'] else 'undefined',
         f"{float(row['ae4_null_cumulative_deficit_percent']):.5f}%" if row['ae4_null_cumulative_deficit_percent'] else 'undefined',
         f"{float(row['ae4_5pct_cumulative_lost_pL']):.8f}",f"{float(row['ae4_null_cumulative_lost_pL']):.8f}"]
        for row in read_csv(OUT/'flow_and_cumulative_comparison.csv')]),
    "The relative cumulative deficits grow early, reach their largest reported values at 360 s (3.6871% / 4.1007%), and then decline modestly to 3.4176% / 3.8575%. They persist through 600 s rather than collapsing toward the small Task 39 endpoints. Absolute cumulative lost volume continues to increase across the reported checkpoints. Final flow rates remain 2.4073% / 2.8283% below WT, so the late percentage decline is not recovery of the lost absolute secretion. No infinite-time or steady-state conclusion is claimed.",
    "Final cumulative secretion is 0.9925245441 pL for WT, 0.9586042475 pL at 5%, and 0.9542380458 pL in null. The ratios are 0.9658242239 and 0.9614251370, with losses of 0.0339202966 and 0.0382864983 pL. Task 39's deficits were only 0.48216% / 0.62732%, but its WT cumulative was 0.9586989000 pL. The larger Task 40 deficit partly reflects its changed WT reference; the reported relative effects do not imply the same change in absolute perturbation secretion between tasks.",
    f"The held-out approximately 35% experimental null deficit is much larger than the predicted {null['cumulative_deficit_percent']:.4f}% cumulative deficit, a descriptive gap of {summary['held_out_validation']['gap_percentage_points']:.4f} percentage points. The endpoint flow deficit is only 2.8283%. The experimental number was never a fit target, gate, root selector, or trigger for another mechanism; exact matching of experimental protocol and observation model is not asserted.",
    "Residual maxima and unchanged tolerances for rest and all trajectories are in [conservation_summary.csv](../../results/40_ae4_equal_cation_routing/conservation_summary.csv). Across the trajectories, the largest cell charge-rate residual is 2.572e−13 fmol/s against 1e−10; carbon residual 5.551e−17 fmol/s against 1e−10; water residual 4.337e−19 pL/s against 1e−12; basolateral current residual 2.478e−23 A against 1e−20; and speciation residual 6.128e−12 mM against 1e−9. AE4 charge residual is exactly zero. These are sampled production checks, not a continuous-time mathematical proof.",
    f"Execution used one scientific agent/worker, one BLAS thread, one WT stationary solve, exactly three production integrations, zero numerical retries, zero calibration/optimisation/sweeps, and no genotype-specific rest. Recorded numerical execution was {budget['numerical_execution_seconds']:.3f} seconds, below 900 seconds. The one source-test assertion correction is recorded separately from numerical retries. No split other than 50:50 was tested; no scientific parameter or second mechanism was changed. Scientific execution is locked complete. Compact UTF-8 source, tests and outputs are published only on the designated branch through the connected GitHub integration; main is not modified or merged.",
    "Reproduction order on a fresh checkout without Task 40 outputs: `python analysis/40_ae4_equal_cation_routing/run_preflight.py`, then `run_rest.py`, followed by `run_case.py wt`, `run_case.py ae4_5pct`, and `run_case.py ae4_null` in the same directory. A saved completion ledger prevents reruns. The report-only scripts `summarize_results.py` and `write_report.py` read saved outputs and do not import or execute the model."
    ]
    (HERE/'final_answer.md').write_text('\n\n'.join(parts)+'\n')
    print('Wrote final_answer.md from completed results.')


if __name__=='__main__':main()
