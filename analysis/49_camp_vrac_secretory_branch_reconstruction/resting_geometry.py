"""49G: one complete fixed source matrix, one SVD, no fitted parameters.

The four central states are conditional algebraic lifts of measured Cl/pH.
Centred differences below are local measurement/nuisance derivatives only;
they never rebuild the source matrix, optimise or generate trajectories.
"""
from __future__ import annotations
import csv
import json
import math
import sys
import numpy as np
from checkpoint import ROOT, HERE, OUT, dump, sha, verified

PARENT = ROOT/'analysis/48_joint_experimental_constraint_reconstruction'
sys.path.insert(0, str(PARENT))
from protocol_layer import build_model, GENOTYPES
from modern_full_model.acid_base import _carbon_fractions, speciate, total_alkalinity_mM
from modern_full_model.nhe1_cha2009 import cha_nhe1_flux_fmol_s, published_cha_kinetics
from modern_full_model.nbc_minimal import nbc_affinity_log

ROWS = ['Na_i','K_i','Cl_i','TIC_i','TA_i','V_i',
        'Na_l','K_l','Cl_l','TIC_l','TA_l','V_l',
        'charge_i','charge_l','current_a_minus_para','current_b_plus_para']
COLS = ['NKCC1','NHE1','AE2','AE4','NBC','pump_a','pump_b','K_a','K_b',
        'TMEM16A_Cl','VRAC_Cl','CO2_b','CO2_a','para_Na','para_K','para_Cl',
        'para_HCO3','water_b','water_a','water_para','outflow']
LATENTS = ['Na_i_mM','Na_l_mM','K_l_mM','Cl_l_mM','TIC_l_mM','V_l_pL']


class Analysis:
    def __init__(self):
        self.dependency = verified('49E')
        self.model = build_model('REST')
        self.p = self.model.parameters
        assert self.p.homeostasis.nhe1_model == 'cha_2009_eight_state_mod2'
        assert self.p.membranes.g_apical_background_S == 0
        assert self.p.membranes.g_basolateral_background_S == 0
        self.evaluations = 0
        self.matrix_assemblies = 0
        self.svd_calls = 0
        p = self.p
        self.bath = speciate(total_carbon_mM=p.bath.tic_mM,
            total_alkalinity_mM_value=total_alkalinity_mM(ph=p.bath.ph,
                total_carbon_mM=p.bath.tic_mM, buffer_total_mM=p.bath.buffer_total_mM,
                buffer_pka=p.acid_base.lumen_buffer_pka, parameters=p.acid_base),
            buffer_total_mM=p.bath.buffer_total_mM,
            buffer_pka=p.acid_base.lumen_buffer_pka, parameters=p.acid_base)
        self.osm_b = p.bath.na_mM+p.bath.k_mM+p.bath.cl_mM+p.bath.tic_mM+p.bath.untracked_osmolyte_mM
        self.pb = p.homeostasis.co2_basolateral_permeability_fmol_s_mM
        self.pa = p.homeostasis.co2_apical_permeability_fmol_s_mM

    def lift(self, cl, ph, latent):
        """Enforce measured Cl/pH, both charges, C-H=0 and cell water=0."""
        p = self.p
        na, nal, kl, cll, ticl, vl = latent
        tal = nal+kl-cll
        lab = speciate(total_carbon_mM=ticl,total_alkalinity_mM_value=tal,
            buffer_total_mM=p.geometry.lumen_buffer_total_fmol/vl,
            buffer_pka=p.acid_base.lumen_buffer_pka,parameters=p.acid_base)
        H = cha_nhe1_flux_fmol_s(na_i_mM=na,h_i_mM=10**(3-ph),
            na_o_mM=p.bath.na_mM,h_o_mM=self.bath.h_mM,
            kinetics=published_cha_kinetics(),
            carrier_amount_fmol=p.homeostasis.nhe1_cha_carrier_amount_fmol,
            activity_scale=1.)
        cbar = (self.pb*self.bath.co2_mM+self.pa*lab.co2_mM)/(self.pb+self.pa)
        co2 = cbar-H/(self.pb+self.pa)
        a0,a1,a2 = _carbon_fractions(ph,p.acid_base)
        tic = co2/a0
        a = a1+2*a2
        b = 1/(1+10**(p.acid_base.cell_buffer_pka-ph))
        w = 10**(3+ph-p.acid_base.water_pkw)-10**(3-ph)
        osm_l = nal+kl+cll+ticl
        lb,la = p.water.basolateral_hydraulic_pL_s_mOsm,p.water.apical_hydraulic_pL_s_mOsm
        osm_i = (lb*self.osm_b+la*osm_l)/(lb+la)
        B,X,M = (p.geometry.cell_buffer_total_fmol,
                 p.geometry.fixed_cell_anion_equivalents_fmol,
                 p.geometry.cell_impermeant_osmoles_fmol)
        v = (X+B*b+M+B)/(osm_i-2*cl-(a+1)*tic-w)
        k = cl+a*tic+w+(B*b+X)/v-na
        ta = a*tic+w+B*b/v
        if min(co2,tic,v,k,na,cl,vl,ticl,tal) <= 0:
            raise ValueError('Nonpositive conditional lift; no clipping or new fit permitted')
        y = np.array([na*v,k*v,cl*v,tic*v,ta*v,v,
                      nal*vl,kl*vl,cll*vl,ticl*vl,tal*vl,vl,0.])
        return y,dict(H=H,Cbar_mM=cbar,CO2_i_mM=co2,CO2_l_mM=lab.co2_mM,
                      osm_i_target_mOsm=osm_i)

    def evaluate(self, cl, ph, latent, genotype):
        y, lift = self.lift(cl,ph,latent)
        ev = self.model.evaluate(0,y,genotype=GENOTYPES[genotype])
        self.evaluations += 1
        d = ev.diagnostics
        f = ev.rhs[:12]
        conv = 1e15/self.p.constants.faraday_C_mol
        rates = np.r_[f, f[0]+f[1]-f[2]-f[4], f[6]+f[7]-f[8]-f[10],
                      d.membranes.current_residuals_A['apical']*conv,
                      d.membranes.current_residuals_A['basolateral']*conv]
        C = sum(d.co2_fluxes_fmol_s.values())
        assert abs(C-d.homeostasis.nhe1_inward_fmol_s) < 1e-10
        assert abs(f[5]) < 1e-10
        assert max(map(abs,d.state_charge_fmol.values())) < 1e-10
        assert abs(d.observables.cell_acid_base.ph-ph) < 1e-10
        assert abs(d.observables.cell_concentrations_mM['cl']-cl) < 1e-10
        assert abs(ev.rhs[12]) < 1e-12
        return y,lift,ev,rates

    def source_coefficients(self, ev):
        d=ev.diagnostics; h=d.homeostasis; m=d.membranes; w=d.water
        cv=1e15/self.p.constants.faraday_C_mol; i=m.currents_A
        return np.array([h.nkcc1_inward_fmol_s,h.nhe1_inward_fmol_s,h.ae2_inward_fmol_s,
            d.ae4.cl_cell_fmol_s,0.,m.pump_apical_fmol_s,m.pump_basolateral_fmol_s,
            i['k_apical']*cv,i['k_basolateral']*cv,-i['cl_apical']*cv,0.,
            d.co2_fluxes_fmol_s['bath_to_cell'],d.co2_fluxes_fmol_s['lumen_to_cell'],
            i['para_na']*cv,i['para_k']*cv,-i['para_cl']*cv,-i['para_hco3']*cv,
            w.bath_to_cell_pL_s,w.cell_to_lumen_pL_s,w.bath_to_lumen_pL_s,w.lumen_outflow_pL_s])

    def assemble_once(self, central):
        assert self.matrix_assemblies == 0
        self.matrix_assemblies += 1
        S=np.zeros((4*16,4*21))
        # Every signature is predeclared; no column selection, capacity fit or search.
        v=np.eye(12)
        signatures=[v[0]+v[1]+2*v[2],v[0]+v[4],v[2]-v[3]-v[4],
            -.5*v[0]-.5*v[1]+v[2]-2*v[3]-2*v[4],v[0]+2*v[3]+2*v[4],
            -3*v[0]+2*v[1]+3*v[6]-2*v[7],-3*v[0]+2*v[1],
            -v[1]+v[7],-v[1],-v[2]+v[8],-v[2]+v[8],v[3],v[3]-v[9],
            -v[6],-v[7],-v[8],-v[9]-v[10],v[5],-v[5]+v[11],v[11]]
        ia=np.zeros(21); ib=np.zeros(21); ip=np.zeros(21)
        ia[[5,7]]=1; ia[[9,10]]=-1
        ib[[4,6,8]]=1
        ip[[13,14]]=1; ip[[15,16]]=-1
        for k,c in enumerate(central):
            y=c['y']; g=c['genotype']
            block=np.zeros((16,21))
            block[:12,:20]=np.array(signatures).T
            block[6:11,20]=-y[6:11]/y[11]; block[11,20]=-1
            if g=='AE4_KO': block[:,3]=0
            if g=='AE2_KO': block[:,2]=0
            block[12]=block[0]+block[1]-block[2]-block[4]
            block[13]=block[6]+block[7]-block[8]-block[10]
            block[14]=ia-ip; block[15]=ib+ip
            S[k*16:(k+1)*16,k*21:(k+1)*21]=block
        return S


def main():
    if (OUT/'resting_geometry.json').exists():
        raise RuntimeError('Result already exists; reuse published output rather than recomputing')
    a=Analysis(); ledger=json.loads((PARENT/'constraints.json').read_text())
    cohorts=[('AE4_WT','WT','ae4_table1','WT'),('AE4_KO','AE4_KO','ae4_table1','AE4_KO'),
             ('AE2_control','WT','ae2_table1','control'),('AE2_KO','AE2_KO','ae2_table1','AE2_KO')]
    central=[]; references={}
    for cohort,g,table,key in cohorts:
        rcl=next(r for r in ledger[table] if r['observable']=='rest_Cl_i')[key]
        rph=next(r for r in ledger[table] if r['observable']=='rest_pH_i')[key]
        path=PARENT/'output/rests'/f'{g}.json'
        cached=json.loads(path.read_text()); ref=np.array(cached['state'])
        latent=np.r_[ref[0]/ref[5],ref[6:10]/ref[11],ref[11]]
        y,lift,ev,rates=a.evaluate(rcl[0],rph[0],latent,g)
        central.append(dict(cohort=cohort,genotype=g,cl=rcl,ph=rph,latent=latent,
                            y=y,lift=lift,ev=ev,rates=rates,j=a.source_coefficients(ev)))
        references[cohort]=dict(path=str(path.relative_to(ROOT)),sha256=sha(path),
                               saved_observables=cached['observables'])
    # This is the sole complete matrix assembly and sole SVD in this analysis.
    S=a.assemble_once(central)
    native=np.concatenate([c['j'] for c in central])
    F=np.concatenate([c['rates'] for c in central])
    decomposition_error=float(np.max(np.abs(S@native-F)))
    assert decomposition_error < 1e-10
    row_scale=np.tile(np.array([1.,1,1,1,1,a.osm_b,1,1,1,1,1,a.osm_b,1,1,1,1]),4)
    SW=row_scale[:,None]*S
    norms=np.sqrt(np.sum(SW*SW,axis=0))
    column_scale=np.where(norms>0,norms,1.)
    A=SW/column_scale
    U,s,Vh=np.linalg.svd(A,full_matrices=True); a.svd_calls+=1
    tol=max(A.shape)*np.finfo(float).eps*s[0]
    rank=int(np.count_nonzero(s>tol))
    target=-F; bw=row_scale*target
    projection=U[:,:rank]@(U[:,:rank].T@bw)
    orthogonal=bw-projection
    # SVD representative is algebra, not a constrained fit or a proposed law.
    z=Vh[:rank].T@((U[:,:rank].T@bw)/s[:rank])
    delta=z/column_scale
    nullspace=Vh[rank:].T/column_scale[:,None]
    assert np.max(np.abs(S@delta-target)) < 1e-10
    assert np.max(np.abs(S@nullspace)) < 1e-10
    # Fixed-state native addition cone: directions at actual affinities,
    # with exact-zero gates/deletions unavailable. This is a derived orientation
    # of the SAME matrix, not another basis or another SVD.
    orientation=np.sign(native)
    C=S*orientation
    separator=np.zeros(64); separator[16+4]=1.
    cone_products=separator@C
    separator_target=float(separator@target)
    assert np.min(cone_products) >= -1e-14
    assert separator_target < 0
    # Explicit analytic coefficient null relations, audited using fixed S.
    identities=[]
    for k,c in enumerate(central):
        if c['genotype']=='WT':
            v=np.zeros(84); v[k*21+3]=1;v[k*21+2]=-2;v[k*21]=.5
            identities.append(dict(cohort=c['cohort'],name='AE4 - 2 AE2 + 0.5 NKCC1',
                max_source_error=float(np.max(np.abs(S@v)))))
        v=np.zeros(84);v[k*21+4]=1;v[k*21+1]=-2;v[k*21+11]=-2
        v[k*21+6]=-1/3;v[k*21+8]=-2/3
        identities.append(dict(cohort=c['cohort'],name='NBC - 2 NHE1 - 2 CO2_b - pump_b/3 - 2 K_b/3',
            max_source_error=float(np.max(np.abs(S@v)))))
    measured_jacobian=np.zeros((64,8)); latent_jacobian=np.zeros((64,24))
    summaries=[]
    for k,c in enumerate(central):
        cl,ph=c['cl'][0],c['ph'][0];latent=c['latent'];g=c['genotype']
        # Fixed numerical differentiation scales, not sampled candidate values.
        for j,h in enumerate([1e-3,1e-5]):
            plus=[cl,ph];minus=[cl,ph];plus[j]+=h;minus[j]-=h
            fp=a.evaluate(*plus,latent,g)[3];fm=a.evaluate(*minus,latent,g)[3]
            measured_jacobian[k*16:(k+1)*16,k*2+j]=(fp-fm)/(2*h)
        for j,value in enumerate(latent):
            h=max(abs(value)*1e-5,1e-8)
            plus=latent.copy();minus=latent.copy();plus[j]+=h;minus[j]-=h
            fp=a.evaluate(cl,ph,plus,g)[3];fm=a.evaluate(cl,ph,minus,g)[3]
            latent_jacobian[k*16:(k+1)*16,k*6+j]=(fp-fm)/(2*h)
        d=c['ev'].diagnostics;ci=d.observables.cell_concentrations_mM
        H=d.homeostasis.nhe1_inward_fmol_s;E=d.homeostasis.ae2_inward_fmol_s
        ccrit=(cl*a.bath.hco3_mM/a.p.bath.cl_mM)*10**(a.p.acid_base.carbon_pka1-ph)
        bound_co2=(a.pb+a.pa)*(c['lift']['Cbar_mM']-ccrit)
        nbca=nbc_affinity_log(na_i_mM=ci['na'],hco3_i_mM=ci['hco3'],
            na_o_mM=a.p.bath.na_mM,hco3_o_mM=a.bath.hco3_mM,
            v_basolateral_V=d.membranes.v_basolateral_V,
            thermal_voltage_V=a.p.constants.thermal_voltage_V)
        summaries.append(dict(cohort=c['cohort'],genotype=g,
            observed_Cl_mean_SE_n=c['cl'],observed_pH_mean_SE_n=c['ph'],
            retained_nuisance=dict(zip(LATENTS,latent.tolist())),
            lifted_state=c['y'].tolist(),cell_concentrations_mM=dict(ci),
            cell_pH=d.observables.cell_acid_base.ph,
            lumen_concentrations_mM=dict(d.observables.lumen_concentrations_mM),
            cell_volume_pL=float(c['y'][5]),lift=c['lift'],
            residual=dict(zip(ROWS,c['rates'].tolist())),
            missing_source=dict(zip(ROWS,(-c['rates']).tolist())),
            native_source_coefficients=dict(zip(COLS,c['j'].tolist())),
            cell_alkalinity_balance=dict(NHE1=H,AE2=E,AE4=d.ae4.cl_cell_fmol_s,
                NBC=0.,residual=H-E-2*d.ae4.cl_cell_fmol_s),
            ae2_affinity_log=d.homeostasis.affinities['AE2'],
            NBC_affinity_audit_only=nbca,
            NBC_audit_limit='Thermodynamic affinity at unchanged native voltage; gate remains zero. No basal NBC model was run.',
            KO_scalar_compatibility=dict(
                applicable=g=='AE4_KO',Ccrit_mM=ccrit,Cbar_mM=c['lift']['Cbar_mM'],
                H_fmol_s=H,AE2_fmol_s=E,frozen_AE2_capacity_fmol_s=a.p.homeostasis.ae2_capacity_fmol_s,
                carbon_sign_upper_bound_fmol_s=bound_co2,
                necessary_positive_H_upper_bound_fmol_s=min(a.p.homeostasis.ae2_capacity_fmol_s,bound_co2),
                scalar_residual_H_minus_AE2=H-E),
            voltage_V=dict(apical=d.membranes.v_apical_V,basolateral=d.membranes.v_basolateral_V)))
    se=np.array([value for c in central for value in [c['cl'][1],c['ph'][1]]])
    covariance=(measured_jacobian*se)@(measured_jacobian*se).T
    sigma=np.sqrt(np.diag(covariance))
    covariance_bound=np.sum(np.abs(measured_jacobian)*se,axis=1)
    for k,summary in enumerate(summaries):
        summary['residual_linearised_SE_independence_working_assumption']=dict(zip(ROWS,sigma[k*16:(k+1)*16].tolist()))
        summary['residual_linearised_SE_upper_bound_unknown_Cl_pH_correlation']=dict(zip(ROWS,covariance_bound[k*16:(k+1)*16].tolist()))
    row_names=[c['cohort']+':'+n for c in central for n in ROWS]
    column_names=[c['cohort']+':'+n for c in central for n in COLS]
    dump(OUT/'fixed_source_matrix.json',dict(rows=row_names,columns=column_names,matrix=S.tolist(),
        raw_row_units=['pL/s' if n.endswith(':V_i') or n.endswith(':V_l') else
            'fmol-charge/s' if ':charge_' in n or ':current_' in n else 'fmol/s' for n in row_names],
        raw_coefficient_units='fmol/s for solute cycles; pL/s for water/outflow',
        audit_row_units='fmol-charge/s; multiply by F*1e-15 for current in amperes',
        row_multipliers=row_scale.tolist(),column_divisors=column_scale.tolist(),
        metric='Water rows times bath osmolarity; all other rows unity; then normalise each nonzero column.',
        coefficient_meaning='State-specific physical fluxes, NOT independently fitted genotype multipliers.',
        shared_law_constraint='All evaluations use the same frozen parameters and genotype masks; no tying matrix is fitted. The block matrix is a relaxed source-space feasibility test. Cone failure in one cohort survives any additional shared-law/current constraints.',
        native_flux_coefficients=native.tolist(),native_addition_cone_orientation=orientation.tolist(),
        zero_availability='Native resting NBC and VRAC gates are zero; genotype deletions are zero columns.',
        source_decomposition_max_error=decomposition_error))
    dump(OUT/'source_svd_projection.json',dict(rank=rank,shape=list(S.shape),
        nullity=84-rank,tolerance=tol,singular_values=s.tolist(),
        target_missing_source=target.tolist(),scaled_span_projection=projection.tolist(),
        scaled_orthogonal_residual=orthogonal.tolist(),
        scaled_orthogonal_residual_norm=float(np.linalg.norm(orthogonal)),
        scaled_target_norm=float(np.linalg.norm(bw)),
        min_norm_coefficient_representative=delta.tolist(),
        min_norm_definition='Minimises Euclidean norm of z=column_divisors*delta in the column-normalised matrix, not unweighted raw delta.',
        coefficient_nullspace=nullspace.tolist(),
        exact_identity_checks=identities,
        warning='Representative coefficients are not an accepted pathway law. F=Sj makes span membership automatic; arbitrary delta=-j or all capacities zero is not a physiological solution.',
        cone_certificate=dict(functional=separator.tolist(),
            functional_times_oriented_columns=cone_products.tolist(),
            functional_times_target=separator_target,
            lower_bound_scaled_cone_distance=abs(separator_target),
            conclusion='Required source lies outside the native fixed-state addition cone.',
            stronger_fixed_state_statement='Nonnegative reweighting of total NHE/AE2 capacities cannot balance KO TA unless both nonzero fluxes are suppressed to zero; this is not a valid inherited-law reconstruction.',
            limit='Conditional states/voltages and native affinity signs. Not a global impossibility theorem on the unmeasured-state manifold.'),
        constrained_least_squares_calls=0))
    dump(OUT/'resting_uncertainty.json',dict(residual_rows=row_names,
        units='Jacobian entries: raw residual-row unit divided by measurement/nuisance-coordinate unit. Covariance entries: product of residual-row units. Raw units are listed in fixed_source_matrix.json.',
        roundoff_warning='Derivatives and propagated errors in exactly imposed water, charge and current identities are numerical noise, not biological uncertainty.',
        measured_columns=[c['cohort']+':'+n for c in central for n in ['Cl_i_mM','pH_i']],
        measurement_jacobian=measured_jacobian.tolist(),reported_SE=se.tolist(),
        residual_covariance_working_independent_errors=covariance.tolist(),
        maximum_linearised_SE_unknown_within_cohort_correlation=covariance_bound.tolist(),
        nuisance_columns=[c['cohort']+':'+n for c in central for n in LATENTS],
        nuisance_jacobian=latent_jacobian.tolist(),
        nuisance_errors='Unreported; no invented prior or uncertainty assigned. The two eliminated cell degrees are fixed by C-H=0 and cell-water=0 for this conditional lift; six retained degrees remain.',
        derivative_steps=dict(Cl_mM=.001,pH=.00001,nuisance='max(abs(reference)*1e-5,1e-8)'),
        interpretation='Local sensitivities only, not parameter sweeps, confidence intervals, physiological bounds or proof of nuisance-robust cone separation. Missing-source derivatives are the negatives of these residual derivatives.'))
    result=dict(checkpoint='49G',verified_49E_sha=a.dependency['sha'],
        status='SPAN_CONTAINS_RESIDUAL; NATIVE_FIXED_STATE_CONE_FAILS; UNIQUE_SHARED_CORRECTION_NOT_IDENTIFIED',
        input_ledger_sha256=sha(PARENT/'constraints.json'),cached_references=references,
        cohorts=summaries,rank=rank,nullity=84-rank,
        distinct_fixed_source_bases=1,
        accepted_execution_matrix_assemblies=a.matrix_assemblies,accepted_execution_svd_calls=a.svd_calls,
        execution_provenance=dict(attempts=2,total_matrix_assemblies=2,total_svd_calls=2,
            failed_attempt_RHS_evaluations=20,
            failure='Initial identical calculation reached output export and raised KeyError for diagnostic name ae2; corrected to source key AE2. No outputs survived that process. Identical basis/inputs were replayed once; no alternative basis, refit or model selection.'),
        full_RHS_evaluations=a.evaluations,production_trajectories=0,
        fitted_parameters=0,constrained_least_squares_calls=0,model_edits=0,
        exact_balance='TA_dot=H+2B-E-2A; TIC_dot=C+2B-E-2A; rest requires C=H. Native AE4-KO requires H=E=C.',
        exact_KO_scalar_equation='H=G2*tanh(log((Cbar-H/D)/Ccrit)/w); 0<H<min(G2,D*(Cbar-Ccrit)).',
        correction_gate='No uniquely identified shared correction; 49H and49I inapplicable.',
        scope='Cl/pH each leave eight charge-compatible hidden state degrees; the declared lift uses two stationary identities and retains six cached nuisance coordinates. Conditional source vectors are not directly measured fluxes.',
        source_hash_verification='output/immutable_source_check_49G.json; separate byte-hash and branch-scope audit before publication')
    dump(OUT/'resting_geometry.json',result)
    with (OUT/'resting_residuals.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['cohort','coordinate','residual','required_source','working_SE','unknown_correlation_SE_upper_bound'])
        for k,c in enumerate(central):
            for j,row in enumerate(ROWS):
                i=k*16+j;writer.writerow([c['cohort'],row,F[i],-F[i],sigma[i],covariance_bound[i]])
    print(json.dumps(dict(rank=rank,nullity=84-rank,orthogonal_norm=float(np.linalg.norm(orthogonal)),
        decomposition_error=decomposition_error,KO_TA_required=separator_target,
        central=[dict(cohort=s['cohort'],TA_residual=s['residual']['TA_i'],H=s['cell_alkalinity_balance']['NHE1'],
                      AE2=s['cell_alkalinity_balance']['AE2'],NBC_affinity=s['NBC_affinity_audit_only']) for s in summaries],
        RHS_evaluations=a.evaluations)))

if __name__=='__main__': main()
