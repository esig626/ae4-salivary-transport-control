"""Read only scientific interpretation of the two completed baseline cases."""
from common import *
import gzip
import math

def read_rows(path):
    path=Path(path)
    opener=path.open if path.exists() else lambda: gzip.open(str(path)+'.gz','rt')
    with opener() as handle:
        return [{k:float(v) for k,v in row.items()} for row in csv.DictReader(handle)]

def main():
    out=HERE/'output'
    wt=json.loads((out/'wt/summary.json').read_text());ko=json.loads((out/'ae4_null/summary.json').read_text())
    wi,ki=wt['integrals']['60_600'],ko['integrals']['60_600']
    w,k=wt['endpoint'],ko['endpoint']
    delta=lambda field:ki[field]-wi[field]
    dp=delta('pump_total_cycles_fmol_s')
    attribution=dict(
        extra_nkcc_k_fmol=delta('k_nkcc_influx_fmol_s'),
        removed_ae4_k_export_fmol=wi['k_ae4_export_fmol_s']-ki['k_ae4_export_fmol_s'],
        extra_pump_k_influx_fmol=2*dp,
        extra_k_channel_export_fmol=delta('k_total_efflux_fmol_s'),
        extra_k_storage_integrated_rhs_fmol=delta('k_storage_fmol_s'))
    assert abs(attribution['extra_nkcc_k_fmol']+attribution['removed_ae4_k_export_fmol']+attribution['extra_pump_k_influx_fmol']-attribution['extra_k_channel_export_fmol']-attribution['extra_k_storage_integrated_rhs_fmol'])<1e-8
    log_parts={ion:math.log(k[ion+'_i_mM']/w[ion+'_i_mM'])*(2 if ion=='cl' else 1) for ion in ['na','k','cl']}
    affinity_product={'selected_reversal_mM4':157.5/2.0096e-5,
        'ideal_bath_reversal_mM4':145*5*126.16159280366088**2}
    data=dict(
        window_s=[60,600],k_channel_increase_percent=100*delta('k_total_efflux_fmol_s')/wi['k_total_efflux_fmol_s'],
        pump_increase_percent=100*dp/wi['pump_total_cycles_fmol_s'],
        nkcc_increase_percent=100*delta('k_nkcc_influx_fmol_s')/wi['k_nkcc_influx_fmol_s'],
        replaced_ae4_chloride_fraction=2*delta('k_nkcc_influx_fmol_s')/wi['ae4_cl_inward_fmol_s'],
        potassium_attribution=attribution,endpoint_log_product_attribution=log_parts,
        endpoint_log_product_change=sum(log_parts.values()),
        endpoint_extra_basal_drive_mV=1000*(k['k_drive_basolateral_V']-w['k_drive_basolateral_V']),
        endpoint_basal_voltage_change_mV=1000*(k['v_basolateral_V']-w['v_basolateral_V']),
        endpoint_basal_reversal_change_mV=1000*(k['e_k_basolateral_V']-w['e_k_basolateral_V']),
        source_reversal_products=affinity_product,
        no_parameter_change=True,no_null_magnitude_comparison=True)
    # Integrated storage is checked against actual endpoint amounts, independently
    # of the differential residual. Small residuals are quadrature error.
    storage_checks={}
    for case in ['wt','ae4_null']:
        rows=read_rows(out/case/'trajectory.csv');start=next(r for r in rows if r['time_s']==60)
        end=rows[-1];summary=wt if case=='wt' else ko
        s=summary['integrals']['60_600']
        errors={ion:s[ion+'_storage_fmol_s']-(end[ion+'_amount_fmol']-start[ion+'_amount_fmol']) for ion in ['na','k','cl']}
        assert max(abs(v) for v in errors.values())<5e-4
        storage_checks[case]=errors
    data['storage_quadrature_error_fmol']=storage_checks
    dump(out/'baseline_interpretation.json',data)
    table=[]
    for label,field in [('K channel efflux','k_total_efflux_fmol_s'),('Apical K efflux','k_apical_efflux_fmol_s'),('Basolateral K efflux','k_basolateral_efflux_fmol_s'),('NKCC K entry','k_nkcc_influx_fmol_s'),('Pump cycles','pump_total_cycles_fmol_s'),('NBC cycles','nbc_cycle_inward_fmol_s'),('NHE1 cycles','nhe1_inward_fmol_s'),('Na net entry before pump','na_net_transport_entry_fmol_s'),('Na pump clearance','na_pump_clearance_fmol_s'),('K amount storage','k_storage_fmol_s')]:
        table.append(f'| {label} | {wi[field]:.6f} | {ki[field]:.6f} |')
    text='''# Frozen baseline potassium diagnostics

Exactly one WT and one AE4 null trajectory used the identical frozen Task 40 parameters and WT rest. All inherited gates passed. There was no CBM solve, resting solve, optimisation or scientific parameter change. The WT metadata recovery is documented in `output/EXECUTION_NOTES.md`; its stored extrema are for the saved grid. The null has 333 accepted steps and 934 checked states.

## Integrated transport, 60 to 600 s

| Quantity (fmol) | WT | AE4 null |
| --- | ---: | ---: |
'''+ '\n'.join(table)+f'''

Total K channel efflux increases by {data['k_channel_increase_percent']:.4f}%, whereas pump cycles increase by {data['pump_increase_percent']:.4f}%. Increased NKCC recruitment and removal of the AE4 K exit add {attribution['extra_nkcc_k_fmol']+attribution['removed_ae4_k_export_fmol']:.6f} fmol to the relative K burden. Extra pump influx adds {attribution['extra_pump_k_influx_fmol']:.6f} fmol. Channels carry {attribution['extra_k_channel_export_fmol']:.6f} fmol more and the difference in K storage is {attribution['extra_k_storage_integrated_rhs_fmol']:.6f} fmol. The latter is a difference between the two trajectories, not the absolute null storage.

The Na burden does not grow in proportion to NKCC alone. NBC influx falls by {wi['nbc_cycle_inward_fmol_s']-ki['nbc_cycle_inward_fmol_s']:.6f} fmol and NHE1 by {wi['nhe1_inward_fmol_s']-ki['nhe1_inward_fmol_s']:.6f} fmol, partly offsetting the extra NKCC Na entry and removal of AE4 Na export. Net Na entry before the pump changes by only {delta('na_net_transport_entry_fmol_s'):.6f} fmol. These are the complete Na and K balances, including changing storage.

## Endpoint mechanism

| Quantity at 600 s | WT | AE4 null |
| --- | ---: | ---: |
| Intracellular K (mM) | {w['k_i_mM']:.6f} | {k['k_i_mM']:.6f} |
| Intracellular Na (mM) | {w['na_i_mM']:.6f} | {k['na_i_mM']:.6f} |
| Intracellular Cl (mM) | {w['cl_i_mM']:.6f} | {k['cl_i_mM']:.6f} |
| Basolateral voltage (mV) | {1000*w['v_basolateral_V']:.6f} | {1000*k['v_basolateral_V']:.6f} |
| Basolateral K reversal (mV) | {1000*w['e_k_basolateral_V']:.6f} | {1000*k['e_k_basolateral_V']:.6f} |
| Basolateral K drive (mV) | {1000*w['k_drive_basolateral_V']:.6f} | {1000*k['k_drive_basolateral_V']:.6f} |
| Apical K drive (mV) | {1000*w['k_drive_apical_V']:.6f} | {1000*k['k_drive_apical_V']:.6f} |
| Total K channel efflux (fmol/s) | {w['k_total_efflux_fmol_s']:.6f} | {k['k_total_efflux_fmol_s']:.6f} |
| Pump / nominal capacity | {w['pump_nominal_utilisation']:.6f} | {k['pump_nominal_utilisation']:.6f} |
| Pump / external K conditioned limit | {w['pump_k_conditioned_utilisation']:.6f} | {k['pump_k_conditioned_utilisation']:.6f} |
| Dynamic storage correction to CBM (fmol/s) | {w['cbm_storage_correction_fmol_s']:.6f} | {k['cbm_storage_correction_fmol_s']:.6f} |

The basal voltage changes by {data['endpoint_basal_voltage_change_mV']:.6f} mV and the K reversal changes by {data['endpoint_basal_reversal_change_mV']:.6f} mV. Together these increase the K driving voltage by {data['endpoint_extra_basal_drive_mV']:.6f} mV. Both genotypes have exactly the same stimulated active K conductance: {(w['g_k_apical_S']+w['g_k_basolateral_S'])*1e9:.6f} nS, split into {w['g_k_apical_S']*1e9:.6f} apical and {w['g_k_basolateral_S']*1e9:.6f} basal nS. No channel recruitment change is needed for the extra current in this model.

K feedback on NKCC exists and opposes compensation. At 600 s its contribution to log X is {log_parts['k']:.6f}; Na contributes {log_parts['na']:.6f} and twice log Cl contributes {log_parts['cl']:.6f}. The total {sum(log_parts.values()):.6f} is negative, so the chloride fall dominates K accumulation and N increases through the unchanged decreasing concentration law. This is an exact endpoint concentration product decomposition, not an independent intervention proving the origin of each state change.

The channel active law ratio is one by construction. The ratio to the fully open current is {w['k_calcium_gate']:.6f}; it is fixed by the unchanged calcium input and is not free additional reserve. The maximum K conditioned pump utilisation is {wt['maxima']['pump_k_conditioned_utilisation']:.6f} in the WT saved samples and {ko['maxima']['pump_k_conditioned_utilisation']:.6f} over all monitored null states. The pump does not approach its kinetic endpoint because the extra net sodium burden is small.

Both ideal bath and selected NKCC affinities remain positive at every sampled state. The reversal products nevertheless differ: {affinity_product['selected_reversal_mM4']:.6f} versus {affinity_product['ideal_bath_reversal_mM4']:.6f} mM^4. The fixed bath law is not globally thermodynamically validated. No observed flux has the wrong sign under the actual bath in these runs.

Exact voltage closure derivatives were checked for K conductance and pump capacity at 60, 300 and 600 s in both genotypes, using two finite difference steps. These 24 checks are local state conditioned calculations, not alternate trajectories or a parameter sweep. They verify electrical coupling but do not identify physiological parameter values. Pump voltage dependence remains absent from the inherited law.

## Decision before phenotype comparison

The existing finite recycling block accommodates the extra burden without pump saturation or a sampled physiology failure. This rejects an already binding K recycling ceiling as the explanation in the frozen architecture. It does not show that its unmeasured effective conductance is correct in real acini, nor rule out a differently constrained biological recycling system.

No uniquely supported correction to K conductance, gating, distribution or pump capacity has been identified. A missing detailed pump voltage law or mixed BK/IK gating requires kinetic and abundance information, not an invented scalar reduction. The fixed bath NKCC discrepancy is recorded explicitly; all realised signs agree, and a bath dependent kinetic replacement is not uniquely specified by its reversal alone. No numerical correction is made. The raw WT and null prediction outputs are ready to be frozen before the phenotype comparison; the 5% case has not been evaluated.
'''
    (HERE/'BASELINE_DIAGNOSTICS.md').write_text(text)
    print(json.dumps({k:data[k] for k in ['k_channel_increase_percent','pump_increase_percent','endpoint_extra_basal_drive_mV']}))

if __name__=='__main__':main()
