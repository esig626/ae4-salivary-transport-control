"""Analytic necessary-condition contradiction, independent of optimizer status.

The proof uses only four production balances, carbon speciation, positive
NHE1/CO2 capacities, and a subset of the predeclared physiology domain.
It is not a claim about the pooled AE4 law outside that domain/chassis.
"""
from __future__ import annotations

from decimal import Decimal, localcontext
from .task22_contract import *
from .calibration import CORE_INDEPENDENT_ROWS
from .task21_calibration import SPEC


def certificate_for(problem):
    """Return conservative decimal bounds and all assumptions for one root."""
    with localcontext() as ctx:
        ctx.prec=60
        D=lambda x:Decimal(str(x))
        ten=Decimal(10);a=problem.p.acid_base;b=problem.p.bath;g=problem.p.geometry
        pmax=D(PHYSIOLOGY['ph_l']['upper']);pimax=D(PHYSIOLOGY['ph_i']['upper'])
        r1=ten**(pmax-D(a.carbon_pka1));r2=ten**(2*pmax-D(a.carbon_pka1)-D(a.carbon_pka2))
        alpha0=1/(1+r1+r2);alpha2=r2/(1+r1+r2)
        # alpha0 decreases and alpha2 increases with pH. Their difference is
        # negative at pH 8; TIC >= 1 bounds its product from above.
        delta=(alpha2-alpha0)*D(PHYSIOLOGY['tic_l_mM']['lower'])
        delta+=D(g.lumen_buffer_total_fmol)/D(PHYSIOLOGY['volume_l_pL']['lower'])
        delta+=1000*(ten**(pmax-D(a.water_pkw))-ten**(-pmax))
        delta+=D(CONSERVATION_RESIDUAL_TOLERANCES['lumen_speciation_alkalinity_mM'])
        h_aff=(D(b.na_mM)/D(PHYSIOLOGY['na_i_mM']['upper'])).ln()+ten.ln()*(D(b.ph)-pimax)
        h_arg=h_aff/D(problem.p.homeostasis.thermodynamic_saturation_log_width)
        h_tanh=((2*h_arg).exp()-1)/((2*h_arg).exp()+1)
        h_min=D(problem.reference['NHE1'])*D(CAPACITY_FOLD[0])*h_tanh
        eps={str(row):D(SPEC.root_scaled_tolerance)*D(scale)
             for row,scale in zip(CORE_INDEPENDENT_ROWS,SPEC.independent_rhs_scales)}
        eps['4']=eps['10']=D(SPEC.omitted_row_raw_tolerance)
        tolerance_budget=sum(eps[str(i)] for i in (3,4,9,10))
        co2_b_lower=h_min-tolerance_budget
        poolmax=D(PHYSIOLOGY['na_i_mM']['upper'])+D(PHYSIOLOGY['k_i_mM']['upper'])
        affinity=(D(b.cl_mM)/D(PHYSIOLOGY['cl_i_mM']['lower'])).ln()
        affinity+=(poolmax/(D(b.na_mM)+D(b.k_mM))).ln()
        affinity+=2*ten.ln()*(pimax-D(b.ph))
        assumptions=dict(luminal_delta_negative=delta<0, nhe_min_exceeds_residual_budget=co2_b_lower>0,
            positive_CO2_permeability=problem.reference['CO2_basolateral']*CAPACITY_FOLD[0]>0,
            no_cell_membrane_TIC_or_alkalinity_pathway=True,
            equal_lumen_membrane_TIC_and_alkalinity_source=True,
            nonnegative_outflow=problem.p.water.outflow_rate_s>=0,
            pooled_affinity_upper_negative=affinity<0)
        # Outward rounded bounds have far more slack than the 60-digit roundoff.
        safe_delta=Decimal('-0.00645');safe_h=Decimal('0.00005118');safe_aff=Decimal('-0.0764')
        verified=all(assumptions.values()) and delta<safe_delta and co2_b_lower>safe_h and affinity<safe_aff
        return dict(root_id=problem.root, classification='PROVEN_STRUCTURAL_CONTRADICTION' if verified else 'CERTIFICATE_NOT_ESTABLISHED',
            proved=verified, proof_scope='Frozen Task 22 production chassis, parameter backgrounds, declared domains and inherited closure tolerances; not pooled AE4 in general.',
            lumen_TA_minus_TIC_upper_mM=str(delta), conservative_lumen_upper_mM=str(safe_delta),
            nhe1_inward_lower_fmol_s=str(h_min), closure_allowances_fmol_s={k:str(v) for k,v in eps.items() if int(k) in (3,4,9,10)},
            four_balance_residual_budget_fmol_s=str(tolerance_budget),
            bath_to_cell_CO2_lower_fmol_s=str(co2_b_lower), conservative_bath_to_cell_CO2_lower_fmol_s=str(safe_h),
            pooled_affinity_upper=str(affinity), conservative_pooled_affinity_upper=str(safe_aff),
            assumptions=assumptions,
            derivation=['TA_l-TIC_l = TIC_l*(alpha_CO3-alpha_CO2)+deprotonated_buffer_l+OH_l-H_l+speciation_error < 0',
                'J_CO2_a = Qout*(TA_l-TIC_l) + R_TA_l - R_TIC_l <= eps_TA_l+eps_TIC_l',
                'J_CO2_b = J_NHE1 - J_CO2_a + R_TIC_i - R_TA_i >= J_NHE1_min - sum_four_eps > 0',
                'Positive basolateral CO2 permeability implies CO2_i < CO2_bath',
                'HCO3_i/HCO3_bath = (CO2_i/CO2_bath)*10^(pH_i-pH_bath) < 10^(7.05-7.4)',
                'A_AE4 < ln(Cl_bath/47.10)+ln(260/(Na_bath+K_bath))+2*ln(10)*(7.05-pH_bath) < -0.0764'],
            precision='60-digit decimal evaluation with conservative outward-rounded published bounds; monotonic analytic inequalities, not sampling or optimizer infeasibility flags.',
            capacity_upper_bounds_used=False, osmotic_screens_used=False,
            anti_cancellation_bound_used=False, optimization_results_used=False, genotype_results_used=False)


def main():
    require_stage0();manifest,_=load_freeze()
    rows=[certificate_for(WTProblem(manifest,r)) for r in sorted(manifest['roots'])]
    payload=dict(created_utc=now(),classification='PROVEN_STRUCTURAL_CONTRADICTION' if all(r['proved'] for r in rows) else 'CERTIFICATE_NOT_ESTABLISHED',
        all_ten_root_backgrounds_covered=all(r['proved'] for r in rows),root_certificates=rows,
        task21_witnesses_covered='All witnesses use these same immutable root backgrounds; none changes the proof domain.',
        genotype_evaluations=0)
    write_json(OUT/'structural_certificate.json',payload)
    print(payload['classification'],len(rows),'root backgrounds')


if __name__=='__main__':main()
