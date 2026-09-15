"""Two frozen production cases before comparison; 5% only after remote freeze.

Run from the repository root. Outputs cannot be overwritten. No optimiser,
stationary solver, phenotype target or alternate production parameter occurs.
"""
from common import *
import argparse
import math
import time
import scipy
from scipy.integrate import Radau, cumulative_trapezoid
from modern_full_model.validation import PRODUCTION_RADAU
from modern_full_model.nkcc1_palk2010 import A1, A2_MM, A3, A4_MM

DIAGNOSE = inherited_diagnose()

def diagnostics(model, t, y, genotype):
    ev = model.evaluate(float(t), y, genotype=genotype)
    row, failures, ratios = DIAGNOSE(model, t, y, ev)
    if failures:
        raise RuntimeError(f'Inherited gate failure at {t}: {failures}')
    d, mp = ev.diagnostics, model.parameters.membranes
    m, h, reg = d.membranes, d.homeostasis, d.regulatory
    ci, li = d.observables.cell_concentrations_mM, d.observables.lumen_concentrations_mM
    constants, bath = model.parameters.constants, model.parameters.bath
    f = constants.faraday_C_mol * 1e-15
    vt = constants.thermal_voltage_V
    ca = model.stimulus(float(t)).calcium_uM
    gate = ca**mp.calcium_hill/(ca**mp.calcium_hill + mp.calcium_half_uM**mp.calcium_hill)
    g_a = mp.g_k_total_S*mp.apical_k_fraction*gate + mp.g_apical_background_S
    g_b = mp.g_k_total_S*(1-mp.apical_k_fraction)*gate + mp.g_basolateral_background_S
    ek_a, ek_b = vt*math.log(li['k']/ci['k']), vt*math.log(bath.k_mM/ci['k'])
    drive_a, drive_b = m.v_apical_V-ek_a, m.v_basolateral_V-ek_b
    ka, kb = m.currents_A['k_apical']/f, m.currents_A['k_basolateral']/f
    pa, pb = m.pump_apical_fmol_s, m.pump_basolateral_fmol_s
    p, n, a, b, hh, j = pa+pb, h.nkcc1_inward_fmol_s, d.ae4.cl_cell_fmol_s, reg['minimal_nbc_cycle_inward_fmol_s'], h.nhe1_inward_fmol_s, -m.currents_A['cl_apical']/f
    g_na = ci['na']**3/(ci['na']**3+mp.nak_na_half_mM**3)
    g_ka = li['k']**2/(li['k']**2+mp.nak_k_half_mM**2)
    g_kb = bath.k_mM**2/(bath.k_mM**2+mp.nak_k_half_mM**2)
    pa_limit = mp.nak_capacity_fmol_s*mp.apical_pump_fraction*g_ka
    pb_limit = mp.nak_capacity_fmol_s*(1-mp.apical_pump_fraction)*g_kb
    ceff = pa_limit+pb_limit
    x = ci['na']*ci['k']*ci['cl']**2
    scale = model.nkcc1_kinetics.alpha_eff_fmol_s*reg['nkcc1_capacity_multiplier']
    ndlogx = -scale*x*(A2_MM*A3+A1*A4_MM)/(A3+A4_MM*x)**2
    cbm = j/3+n+hh/3-a/2
    storage = (-2*ev.rhs[0]+ev.rhs[4]+ev.rhs[2])/3-ev.rhs[1]
    independent_na = n+hh+b-a/2-3*p
    independent_k = n-a/2+2*p-ka-kb
    # Independent unit/sign/current and conserved-row checks.
    residuals = dict(
        na_rhs=ev.rhs[0]-independent_na, k_rhs=ev.rhs[1]-independent_k,
        cl_rhs=ev.rhs[2]-(2*n+h.ae2_inward_fmol_s+a-j),
        ta_rhs=ev.rhs[4]-(hh+2*b-h.ae2_inward_fmol_s-2*a),
        dynamic_cbm=ka+kb-cbm-storage,
        storage_charge=storage+(ev.rhs[0]+2*ev.rhs[1])/3,
        current_sum=ka+kb+p+b-j,
        channel_a=ka-g_a*drive_a/f, channel_b=kb-g_b*drive_b/f,
        pump_a=pa-pa_limit*g_na, pump_b=pb-pb_limit*g_na)
    assert max(abs(v) for v in residuals.values()) < 1e-10
    assert m.currents_A['k_apical']*drive_a >= -1e-24
    assert m.currents_A['k_basolateral']*drive_b >= -1e-24
    # No assumed physiological voltage interval or ATP free energy is introduced.
    row.update(dict(
        na_amount_fmol=float(y[0]), k_amount_fmol=float(y[1]), cl_amount_fmol=float(y[2]),
        na_storage_fmol_s=float(ev.rhs[0]), k_storage_fmol_s=float(ev.rhs[1]),
        cl_storage_fmol_s=float(ev.rhs[2]), ta_storage_fmol_s=float(ev.rhs[4]),
        volume_rate_pL_s=float(ev.rhs[5]),
        k_concentration_rate_mM_s=float((ev.rhs[1]-ci['k']*ev.rhs[5])/y[5]),
        k_dilution_rate_mM_s=float(-ci['k']*ev.rhs[5]/y[5]),
        k_l_mM=float(li['k']),
        k_nkcc_influx_fmol_s=float(n), k_pump_influx_fmol_s=float(2*p),
        k_ae4_export_fmol_s=float(a/2), k_apical_efflux_fmol_s=float(ka),
        k_basolateral_efflux_fmol_s=float(kb), k_total_efflux_fmol_s=float(ka+kb),
        k_calcium_gate=float(gate), g_k_apical_S=float(g_a), g_k_basolateral_S=float(g_b),
        e_k_apical_V=float(ek_a), e_k_basolateral_V=float(ek_b),
        k_drive_apical_V=float(drive_a), k_drive_basolateral_V=float(drive_b),
        k_current_apical_A=float(m.currents_A['k_apical']),
        k_current_basolateral_A=float(m.currents_A['k_basolateral']),
        k_active_law_ratio_apical=float(m.currents_A['k_apical']/(g_a*drive_a)),
        k_active_law_ratio_basolateral=float(m.currents_A['k_basolateral']/(g_b*drive_b)),
        k_full_open_ratio_apical=float(m.currents_A['k_apical']/((mp.g_k_total_S*mp.apical_k_fraction+mp.g_apical_background_S)*drive_a)),
        k_full_open_ratio_basolateral=float(m.currents_A['k_basolateral']/((mp.g_k_total_S*(1-mp.apical_k_fraction)+mp.g_basolateral_background_S)*drive_b)),
        pump_apical_cycles_fmol_s=float(pa), pump_basolateral_cycles_fmol_s=float(pb),
        pump_apical_available_fmol_s=float(pa_limit), pump_basolateral_available_fmol_s=float(pb_limit),
        pump_available_fmol_s=float(ceff), pump_nominal_utilisation=float(p/mp.nak_capacity_fmol_s),
        pump_k_conditioned_utilisation=float(p/ceff),
        pump_apical_nominal_utilisation=float(pa/(mp.nak_capacity_fmol_s*mp.apical_pump_fraction)),
        pump_basolateral_nominal_utilisation=float(pb/(mp.nak_capacity_fmol_s*(1-mp.apical_pump_fraction))),
        pump_k_conditioned_reserve_fmol_s=float(ceff-p), pump_na_elasticity=float(3*(1-g_na)),
        na_pump_clearance_fmol_s=float(3*p), na_nkcc_entry_fmol_s=float(n),
        na_nbc_entry_fmol_s=float(b), na_nhe_entry_fmol_s=float(hh),
        na_ae4_export_fmol_s=float(a/2), na_net_transport_entry_fmol_s=float(n+b+hh-a/2),
        nkcc_selected_affinity=float(math.log(A1/(A2_MM*x))),
        nkcc_ideal_bath_affinity=float(math.log(bath.na_mM*bath.k_mM*bath.cl_mM**2/x)),
        nkcc_kinetic_positive_bound_fmol_s=float(scale*A1/A3),
        nkcc_dflux_dlog_k_fmol_s=float(ndlogx),
        nkcc_dflux_dlog_na_fmol_s=float(ndlogx), nkcc_dflux_dlog_cl_fmol_s=float(2*ndlogx),
        cbm_zero_storage_k_fmol_s=float(cbm), cbm_storage_correction_fmol_s=float(storage),
        independent_balance_max_abs_fmol_s=float(max(abs(v) for v in residuals.values())),
        pump_apical_electrochemical_work_kJ_mol=float((constants.gas_constant_J_mol_K*constants.temperature_K*(3*math.log(li['na']/ci['na'])+2*math.log(ci['k']/li['k']))-constants.faraday_C_mol*m.v_apical_V)/1000),
        pump_basolateral_electrochemical_work_kJ_mol=float((constants.gas_constant_J_mol_K*constants.temperature_K*(3*math.log(bath.na_mM/ci['na'])+2*math.log(ci['k']/bath.k_mM))-constants.faraday_C_mol*m.v_basolateral_V)/1000)))
    return row, ratios

def local_checks(model, t, y, genotype):
    ev=model.evaluate(t,y,genotype=genotype); d=ev.diagnostics; mp=model.parameters.membranes
    row,_=diagnostics(model,t,y,genotype)
    f=model.parameters.constants.faraday_C_mol*1e-15
    vt=model.parameters.constants.thermal_voltage_V
    ga,gb=row['g_k_apical_S'],row['g_k_basolateral_S']
    gcl=mp.g_cl_apical_S*row['k_calcium_gate']
    gp=mp.g_para_na_S+mp.g_para_k_S+mp.g_para_cl_S+mp.g_para_hco3_S
    r=d.regulatory; nb=model.nbc_parameters
    gnb=f*nb.capacity_fmol_s*r['minimal_nbc_activation_fraction']/(nb.log_width*vt)/math.cosh(r['minimal_nbc_affinity_log']/nb.log_width)**2
    mat=np.array([[ga+gcl+gp,-gp],[-gp,gb+gnb+gp]])
    results=[]
    for name in ['g_k_total_S','nak_capacity_fmol_s']:
        if name=='g_k_total_S':
            direct=np.array([mp.g_k_total_S*mp.apical_k_fraction*row['k_calcium_gate']*row['k_drive_apical_V'],mp.g_k_total_S*(1-mp.apical_k_fraction)*row['k_calcium_gate']*row['k_drive_basolateral_V']])
        else:
            direct=f*np.array([row['pump_apical_cycles_fmol_s'],row['pump_basolateral_cycles_fmol_s']])
        dv=-np.linalg.solve(mat,direct)
        expected=np.r_[dv,-gcl*dv[0]/f,gnb*dv[1]/f]
        for step in [1e-4,5e-5]:
            values=[]
            for direction in [-1,1]:
                mm=clone_parameters(model,**{name:getattr(mp,name)*math.exp(direction*step)})
                ee=mm.evaluate(t,y,genotype=genotype).diagnostics
                values.append(np.array([ee.membranes.v_apical_V,ee.membranes.v_basolateral_V,-ee.membranes.currents_A['cl_apical']/f,ee.regulatory['minimal_nbc_cycle_inward_fmol_s']]))
            fd=(values[1]-values[0])/(2*step)
            scaled=np.max(abs(fd-expected)/(1e-9+abs(expected)))
            assert scaled < 2e-5, (name,t,scaled)
            results.append(dict(time_s=t,parameter=name,log_step=step,
                dVa_dlog_parameter_V=float(expected[0]),dVb_dlog_parameter_V=float(expected[1]),
                dCl_export_dlog_parameter_fmol_s=float(expected[2]),dNBC_dlog_parameter_fmol_s=float(expected[3]),
                max_scaled_derivative_error=float(scaled)))
    return results

def run_case(case):
    expression={'wt':1.0,'ae4_null':0.0,'ae4_5pct':0.05}[case]
    out=HERE/'output'/case
    assert not out.exists(), f'Immutable output already exists: {out}'
    if case=='ae4_5pct':
        receipt=json.loads((HERE/'REMOTE_PREDICTION_RECEIPT.json').read_text())
        assert receipt['remote_verified'] and len(receipt['prediction_commit'])==40
        frozen=json.loads((HERE/'PREDICTION_MANIFEST.json').read_text())
        for p,sha in frozen['sha256'].items():
            assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==sha,p
    model,y0=load_model()
    source_before=source_hashes()
    genotype=replace(WT,ae4_expression=expression)
    started=time.monotonic(); minima={}; maxima={}; ratios_max={}
    rows=[]; states=[]; checks=0; accepted=0
    def inspect(t,y,save=False):
        nonlocal checks
        row,ratios=diagnostics(model,t,y,genotype);checks+=1
        for k,v in row.items():
            minima[k]=min(minima.get(k,v),v);maxima[k]=max(maxima.get(k,v),v)
        for k,v in ratios.items():ratios_max[k]=max(ratios_max.get(k,0),float(v))
        if save:
            rows.append(row);states.append(dict(zip(['time_s',*model.state_names],[float(t),*map(float,y)])))
    inspect(0.,y0,True);inspect(1e-6,y0,True)
    spec=PRODUCTION_RADAU
    solver=Radau(lambda t,y:model.rhs(t,y,genotype=genotype),1e-6,y0.copy(),600.,
                 rtol=spec.rtol,atol=spec.atol_vector(model.state_names),max_step=spec.max_step_s)
    next_time=1.
    while solver.status=='running':
        solver.step()
        if solver.status=='failed':raise RuntimeError(f'{case}: Radau failed')
        accepted+=1;dense=solver.dense_output()
        while next_time<=solver.t:
            inspect(next_time,dense(next_time),True);next_time+=1.
        if solver.t!=next_time-1:inspect(float(solver.t),solver.y)
    assert solver.t==600.
    ts=np.array([r['time_s'] for r in rows]);assert np.array_equal(ts,np.r_[0.,1e-6,np.arange(1.,601.)])
    cumulative=cumulative_trapezoid([r['q_out_pL_s'] for r in rows],ts,initial=0.)
    for r,q in zip(rows,cumulative):r['cumulative_outflow_pL']=float(q)
    fields=[k for k in rows[0] if k.endswith('_fmol_s') or k in ('q_out_pL_s','nbc_current_A')]
    integrals={}
    for start,end in [(0,60),(60,600),(0,600)]:
        mask=(ts>=start)&(ts<=end)
        integrals[f'{start}_{end}']={k:float(np.trapezoid([r[k] for r,keep in zip(rows,mask) if keep],ts[mask])) for k in fields}
    rest=model.evaluate(0,y0).diagnostics.water.lumen_outflow_pL_s
    wt_gates=dict(mean=integrals['60_600']['q_out_pL_s']/540>1.1*rest,
                  cumulative=cumulative[-1]>1.1*600*rest,endpoint=rows[-1]['q_out_pL_s']>=rest)
    assert case!='wt' or all(wt_gates.values())
    local=[]
    if case!='ae4_5pct':
        for t in [60.,300.,600.]:
            s=states[list(ts).index(t)]
            local.extend(local_checks(model,t,np.array([s[k] for k in model.state_names]),genotype))
    assert source_before==source_hashes()
    summary=dict(case=case,expression=expression,parameter_hashes=hashes(model),
        initial_state_sha256=sha256_object(y0.tolist()),initial_state=y0.tolist(),final_state=solver.y.tolist(),
        accepted_steps=accepted,checked_states=checks,nfev=solver.nfev,
        elapsed_seconds=time.monotonic()-started,python=sys.version.split()[0],numpy=np.__version__,scipy=scipy.__version__,
        solver=dict(method='Radau',rtol=spec.rtol,atol=spec.atol_vector(model.state_names).tolist(),max_step_s=spec.max_step_s),
        all_inherited_gates_passed=True,maximum_conservation_ratios=ratios_max,minima=minima,maxima=maxima,
        integrals=integrals,endpoint=rows[-1],wt_activation=wt_gates,
        local_derivative_cases=len(local),production_source_unchanged=True)
    csv_write(out/'trajectory.csv',rows);csv_write(out/'states.csv',states)
    if local:csv_write(out/'local_derivatives.csv',local)
    dump(out/'summary.json',summary)
    print(json.dumps(dict(case=case,status='PASS',checked_states=checks,accepted_steps=accepted,
                         cumulative_pL=float(cumulative[-1]),max_balance_residual=maxima['independent_balance_max_abs_fmol_s'])))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('case',choices=['wt','ae4_null','ae4_5pct'])
    run_case(parser.parse_args().case)
