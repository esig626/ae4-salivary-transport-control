"""Full model analysis, with no changes to historical numerical outputs."""
from core import *
from time import perf_counter
from scipy.linalg import expm
from scipy.optimize import brentq
from modern_full_model.ae4_catalan2025 import forward_affinity
from modern_full_model.nhe1_cha2009 import published_cha_kinetics
from modern_full_model.nkcc1_palk2010 import A1,A2_MM,A3,A4_MM
OUT=Path(__file__).resolve().parent/'output';OUT.mkdir(exist_ok=True)

def main():
    start=perf_counter();rest,stim,y0,task,audit=reference();m=clone(stim);sc=Scaling(m,y0)
    write_json(OUT/'source_audit.json',dict(base_main=BASE_SHA,files=audit,parameter_hashes=task.parameter_hashes(stim)))
    write_json(OUT/'active_parameters.json',{n:asdict(getattr(stim,n)) for n in ['parameters','ae4_parameters','nbc_parameters','nkcc1_kinetics','regulatory_model','stimulus']})
    groups=sc.payload();mp=m.parameters.membranes;wp=m.parameters.water
    groups.update(theta_ae4=30/sc.t0,nbc_capacity_ratio=m.nbc_parameters.capacity_fmol_s/sc.J0,
      pump_capacity_ratio=mp.nak_capacity_fmol_s/sc.J0,ae2_capacity_ratio=m.parameters.homeostasis.ae2_capacity_fmol_s/sc.J0,
      outflow_rate_dimensionless=wp.outflow_rate_s*sc.t0,water_apical=wp.apical_hydraulic_pL_s_mOsm*sc.C0**2/sc.J0,
      water_basolateral=wp.basolateral_hydraulic_pL_s_mOsm*sc.C0**2/sc.J0,
      water_paracellular=wp.paracellular_hydraulic_pL_s_mOsm*sc.C0**2/sc.J0,
      nkcc_reversal_product_mM4=A1/A2_MM,nkcc_forward_bound=m.nkcc1_kinetics.alpha_eff_fmol_s*1.75*A1/A3,
      nkcc_reverse_bound=m.nkcc1_kinetics.alpha_eff_fmol_s*1.75*A2_MM/A4_MM)
    vt=m.parameters.constants.thermal_voltage_V;far=m.parameters.constants.faraday_C_mol
    groups['thermal_voltage_V']=vt
    for k in ['g_k_total_S','g_cl_apical_S','g_para_na_S','g_para_k_S','g_para_cl_S','g_para_hco3_S']:
        groups[k+'_dimensionless']=getattr(mp,k)*vt*1e15/(far*sc.J0)
    groups['NHE_global_bound']=published_cha_kinetics().turnover_bound_per_s*m.parameters.homeostasis.nhe1_cha_carrier_amount_fmol*2.3
    groups['AE4_without_NBC_global_bound']=(groups['NHE_global_bound']+m.parameters.homeostasis.ae2_capacity_fmol_s)/2
    kk=published_cha_kinetics();ab=m.evaluate(1.,y0).diagnostics.observables.bath_acid_base
    no=m.parameters.bath.na_mM;ho=ab.h_mM;hi=1e3*10**(-6.6)
    aa=kk.k1_plus_per_ms*no/(no+kk.k_na_o_mM)*kk.k_h_o_mM/(ho+kk.k_h_o_mM)
    bb=kk.k2_plus_per_ms*hi/(hi+kk.k_h_i_mM)
    dd=kk.k2_minus_per_ms*ho/(ho+kk.k_h_o_mM)*kk.k_na_o_mM/(no+kk.k_na_o_mM)
    mod=1/(1+(1+ho/kk.modifier_k_o_mM)*(kk.modifier_k_i_mM/hi)**3)
    hmax=m.parameters.homeostasis.nhe1_cha_carrier_amount_fmol*2.3*1000*mod*aa*bb/(aa+bb+dd)
    groups.update(physiological_NHE_upper_fmol_s=hmax,
      physiological_AE4_without_NBC_upper_fmol_s=(hmax+.005)/2,
      NBC_lower_for_AE4_demand_0075_fmol_s=.075-(hmax+.005)/2,
      NKCC_eta2=A2_MM*sc.C0**4/A1,NKCC_eta4=A4_MM*sc.C0**4/A3,
      gamma_AE4_reference=m.ae4_parameters.carrier_amount_fmol*m.ae4_parameters.common_cl_attempt_rate_s/sc.J0)
    write_json(OUT/'dimensionless_groups.json',groups)
    wt=integrate(m,y0);ds=integrate(m,y0,scaled=sc)
    ts=np.r_[0.,1e-6,np.arange(1.,601.)];yy=wt.sol(ts);xx=ds.sol(ts/sc.t0)
    obs=[observe(m,t,yy[:-1,i]) for i,t in enumerate(ts)]
    orig=list(csv.DictReader((ROOT/'results/40_ae4_equal_cation_routing/wt_timeseries.csv').open()))
    ot=np.array([float(r['time_s']) for r in orig]);oq=np.array([float(r['q_out_pL_s']) for r in orig])
    nq=np.array([observe(m,t,wt.sol(t)[:-1])['q_pL_s'] for t in ot]);qsample=np.trapezoid([r['q_pL_s'] for r in obs],ts)
    ver=dict(dimensionless_max_scaled_error=float(np.max(abs(yy[:-1]-xx[:-1]*sc.D[:,None])/sc.D[:,None])),
      canonical_flow_max_error=float(abs(oq-nq).max()),original_sampled_cumulative_pL=float(orig[-1]['cumulative_outflow_0_t_pL']),
      reproduced_sampled_cumulative_pL=float(qsample),adaptive_integral_pL=float(wt.y[-1,-1]))
    assert ver['dimensionless_max_scaled_error']<1e-7 and ver['canonical_flow_max_error']<1e-8
    assert abs(ver['reproduced_sampled_cumulative_pL']-ver['original_sampled_cumulative_pL'])<1e-8
    write_json(OUT/'rescaling_verification.json',ver);write_csv(OUT/'wt_trajectory.csv',obs)
    proj=[];rng=np.random.default_rng(44);err=0.
    for name,nu in CLASSES.items():
        mm=with_mechanism_class(m,name);v=nu.conserved_coefficients;kappa=v[2]-2*v[0]+v[4]
        for _ in range(12):
            z=sc.independent(y0);z[:10]*=np.exp(rng.normal(0,.003,10));z[-1]=.5;y=sc.expand(z)
            ev=mm.evaluate(1,y);r=observe(mm,1,y)
            expected=6*r['P']-r['H']+kappa*r['J4']+2*ev.rhs[0]-ev.rhs[4]-ev.rhs[2]
            err=max(err,abs(expected-r['CaCC']))
        d=mm.evaluate(1,y0).diagnostics;ci=dict(d.observables.cell_concentrations_mM);ci['co3']=d.observables.cell_acid_base.co3_mM
        b=m.parameters.bath;ab=d.observables.bath_acid_base
        co=dict(na=b.na_mM,k=b.k_mM,cl=b.cl_mM,hco3=ab.hco3_mM,co3=ab.co3_mM)
        af=forward_affinity(nu,ci,co,voltage_V=d.membranes.v_basolateral_V,thermal_voltage_V=vt)['affinity_over_rt']
        proj.append(dict(mechanism=name,nu_Na=v[0],nu_K=v[1],nu_Cl=v[2],nu_T=v[3],nu_A=v[4],charge=float(nu.exact_charge),kappa=kappa,affinity_reference=af))
    assert err<1e-11
    write_csv(OUT/'source_projections.csv',proj);write_json(OUT/'algebra_verification.json',dict(states_per_class=12,classes=7,max_identity_error_fmol_s=err))
    results=[];failures=[];seed=wt.y[:-1,-1];roots={}
    for e in [1.,.9,.75,.5,.25,.1,.05,0.]:
        rr=stationary(m,sc,seed,Genotype('AE4',ae4_expression=e));rr.update(mechanism='C0',ae4_expression=e)
        results.append(rr);roots[e]=rr;seed=np.array(rr['state'])
    for name in ['C1','C2','C3a','C3b','C4a','C4b']:
        mm=with_mechanism_class(m,name);seed=np.array(roots[1.]['state'])
        for e in [1.,.05,0.]:
            try:
                if e<1.:seed=np.array(roots[e]['state'])
                rr=stationary(mm,sc,seed,Genotype(name,ae4_expression=e));rr.update(mechanism=name,ae4_expression=e)
                results.append(rr);seed=np.array(rr['state'])
            except Exception as ex:failures.append(dict(mechanism=name,ae4_expression=e,error=str(ex),status='root_not_established'))
    write_json(OUT/'equilibria.json',dict(roots=results,failures=failures))
    write_csv(OUT/'equilibrium_summary.csv',[dict(mechanism=r['mechanism'],ae4_expression=r['ae4_expression'],
        physiological=r['physiological'],failed_gates=';'.join(r['failed_gates']),slowest_decay_s=r['slowest_decay_s'],
        max_real_eigenvalue_per_s=r['max_real_eigenvalue_per_s'],residual_max=r['residual_max'],
        jacobian_relative_step_error=r['jacobian_relative_step_error'],condition_number=r['condition_number'],**r['observables']) for r in results])
    print('EQUILIBRIA',len(results),'failures',len(failures),flush=True)
    def quantities(mm,z,g):
        o=observe(mm,1,sc.expand(z),g);return np.array([o[k] for k in ['q_pL_s','N','J4','P','H','ph_i']])
    sens=[];terms=[]
    for e in [1.,.5,.05]:
        rr=roots[e];z=np.array(rr['z']);j=np.array(rr['jacobian_dimensionless']);h=1e-4
        gp=Genotype('plus',ae4_expression=e+h);gm=Genotype('minus',ae4_expression=e-h)
        fmu=(sc.reduced_rhs(m,z,gp)-sc.reduced_rhs(m,z,gm))/(2*h);dz=np.linalg.solve(j,-fmu)
        qx=jacobian(lambda a:quantities(m,a,Genotype('base',ae4_expression=e)),z)
        direct=(quantities(m,z,gp)-quantities(m,z,gm))/(2*h);total=direct+qx@dz
        rp=stationary(m,sc,np.array(rr['state']),gp);rm=stationary(m,sc,np.array(rr['state']),gm)
        fd=(quantities(m,np.array(rp['z']),gp)-quantities(m,np.array(rm['z']),gm))/(2*h)
        row=dict(ae4_expression=e,chloride_compensation_gain=-2*total[1]/total[2],
          q_elasticity=e*total[0]/rr['observables']['q_pL_s'],validation_relative_error=np.linalg.norm(total-fd)/np.linalg.norm(fd))
        for k,a,b in zip(['Q','N','J4','P','H','pH'],total,fd):row['d'+k+'_de']=float(a);row['FD_d'+k+'_de']=float(b)
        sens.append(row)
        for name,value in zip([m.state_names[k] for k in KEEP],qx[0]*dz):terms.append(dict(ae4_expression=e,coordinate=name,dQ_de=float(value)))
    write_csv(OUT/'expression_sensitivities.csv',sens);write_csv(OUT/'coordinate_attribution.csv',terms)
    builders={
      'NBC_capacity':lambda h:clone(stim,nbc_parameters=replace(m.nbc_parameters,capacity_fmol_s=m.nbc_parameters.capacity_fmol_s*np.exp(h))),
      'AE4_carrier':lambda h:clone(stim,ae4_parameters=replace(m.ae4_parameters,carrier_amount_fmol=m.ae4_parameters.carrier_amount_fmol*np.exp(h))),
      'NKCC_scale':lambda h:clone(stim,nkcc1_kinetics=replace(m.nkcc1_kinetics,alpha_eff_fmol_s=m.nkcc1_kinetics.alpha_eff_fmol_s*np.exp(h))),
      'NHE_carrier':lambda h:clone(stim,parameters=replace(m.parameters,homeostasis=replace(m.parameters.homeostasis,nhe1_cha_carrier_amount_fmol=m.parameters.homeostasis.nhe1_cha_carrier_amount_fmol*np.exp(h)))),
      'pump_capacity':lambda h:clone(stim,parameters=replace(m.parameters,membranes=replace(mp,nak_capacity_fmol_s=mp.nak_capacity_fmol_s*np.exp(h)))),
      'CaCC_conductance':lambda h:clone(stim,parameters=replace(m.parameters,membranes=replace(mp,g_cl_apical_S=mp.g_cl_apical_S*np.exp(h)))),
      'buffer_pool':lambda h:clone(stim,parameters=replace(m.parameters,geometry=replace(m.parameters.geometry,cell_buffer_total_fmol=m.parameters.geometry.cell_buffer_total_fmol*np.exp(h)))),
      'water_apical':lambda h:clone(stim,parameters=replace(m.parameters,water=replace(wp,apical_hydraulic_pL_s_mOsm=wp.apical_hydraulic_pL_s_mOsm*np.exp(h))))}
    pars=[];rr=roots[1.];z=np.array(rr['z']);j=np.array(rr['jacobian_dimensionless']);h=1e-4
    for name,build in builders.items():
        pm,mm=build(h),build(-h);fmu=(sc.reduced_rhs(pm,z,WT)-sc.reduced_rhs(mm,z,WT))/(2*h)
        dz=np.linalg.solve(j,-fmu);qx=jacobian(lambda a:quantities(m,a,WT),z)
        total=(quantities(pm,z,WT)-quantities(mm,z,WT))/(2*h)+qx@dz
        rp=stationary(pm,sc,np.array(rr['state']));rm=stationary(mm,sc,np.array(rr['state']))
        fd=(quantities(pm,np.array(rp['z']),WT)-quantities(mm,np.array(rm['z']),WT))/(2*h)
        pars.append(dict(parameter=name,Q_log_elasticity=total[0]/rr['observables']['q_pL_s'],pH_per_log_parameter=total[5],FD_relative_error=np.linalg.norm(total-fd)/max(np.linalg.norm(fd),1e-12)))
    write_csv(OUT/'parameter_sensitivity.csv',pars)
    def frozen_integrals(case):
        return {r['quantity']:float(r['value']) for r in csv.DictReader((ROOT/('results/40_ae4_equal_cation_routing/'+case+'_integrated_fluxes.csv')).open())}
    wfi,nfi=frozen_integrals('wt'),frozen_integrals('ae4_null')
    extra=2*(nfi['nkcc1_cycles']-wfi['nkcc1_cycles'])
    lost=wfi['ae4_signed_cl_flux']-nfi['ae4_signed_cl_flux']
    write_json(OUT/'compensation_comparison.json',dict(signed_integrated_compensation=extra/lost,
      relative_NKCC_increase=extra/(2*wfi['nkcc1_cycles']),extra_NKCC_chloride_fmol=extra,
      lost_AE4_chloride_fmol=lost,source='Frozen one-second integrated fluxes, not the compact 60-second display table.',
      steady_null_deficit=1-roots[0.]['observables']['q_pL_s']/roots[1.]['observables']['q_pL_s']))
    rr=roots[1.];z=np.array(rr['z']);j=np.array(rr['jacobian_dimensionless'])/sc.t0;h=1e-4
    fp=(sc.reduced_rhs(m,z,Genotype('plus',ae4_expression=1+h))-sc.reduced_rhs(m,z,Genotype('minus',ae4_expression=1-h)))/(2*h*sc.t0)
    qx=jacobian(lambda a:quantities(m,a,WT),z)[0];final=-np.linalg.solve(j,fp)
    modes=[]
    for t in [0,30,60,120,300,600,1200,2400]:
        E=expm(j*t);modes.append(dict(time_s=t,linear_dQ_de=float(qx@((np.eye(len(z))-E)@final)),
          asymptotic_dQ_de=float(qx@final),scaled_state_amplification=float(np.linalg.norm(E,2))))
    write_csv(OUT/'local_step_response.csv',modes)
    samples=[];cache={};wtQ=float(wt.y[-1,-1]);g0=Genotype('null',ae4_expression=0.)
    def trial(rho,tight=False):
        key=(float(rho),tight)
        if key in cache:return cache[key]
        mm=cacc(m,float(rho));sol=integrate(mm,y0,g0,rtol=2e-9 if tight else 1e-8,max_step=1. if tight else 2.)
        failures={};ratiosmax=0.;phmin=100.;phmax=0.;clmin=1e9
        for t in np.unique(np.r_[ts,sol.t]):
            y=sol.sol(t)[:-1];ev=mm.evaluate(t,y,genotype=g0);r,ff,ratios=task.parent.diagnose(mm,t,y,ev)
            for f in ff:failures.setdefault(f['gate'],float(t))
            ratiosmax=max(ratiosmax,max(ratios.values()));phmin=min(phmin,r['ph_i']);phmax=max(phmax,r['ph_i']);clmin=min(clmin,r['cl_i_mM'])
        q=[observe(mm,t,sol.sol(t)[:-1],g0)['q_pL_s'] for t in ts]
        rec=dict(rho=float(rho),b=1-float(rho),adaptive_deficit=1-float(sol.y[-1,-1])/wtQ,
          sampled_deficit=1-float(np.trapezoid(q,ts))/qsample,cumulative_adaptive_pL=float(sol.y[-1,-1]),
          physiological=not failures,failed_gates=failures,ph_min=phmin,ph_max=phmax,cl_min=clmin,
          max_conservation_ratio=ratiosmax,nfev=int(sol.nfev),tight=tight)
        cache[key]=rec;samples.append(rec);write_json(OUT/'inverse_samples.json',samples)
        print('INVERSE',rho,rec['adaptive_deficit'],rec['physiological'],flush=True);return rec
    for rho in [0.,.25,.5,.75,.85,.89,1-.10511872843288446,.90,.925,.95]:trial(rho)
    threshold=brentq(lambda r:trial(r)['adaptive_deficit']-.303,.89,.90,xtol=2e-7)
    check=trial(threshold,True);nearby=[trial(threshold+d) for d in [-1e-4,1e-4]]
    write_json(OUT/'inverse_threshold.json',dict(target_lower_deficit=.303,target_upper_deficit=.397,
      numerical_first_crossing_rho=threshold,tight_check=check,nearby=nearby,
      selected_task41_reproduction=trial(1-.10511872843288446),
      monotone_at_sampled_points=bool(np.all(np.diff([r['adaptive_deficit'] for r in sorted([v for k,v in cache.items() if not k[1]],key=lambda r:r['rho'])])>0)),
      scope='First numerical crossing in sampled connected family, not a proven global or universal lower bound.'))
    legacy_threshold=brentq(lambda r:trial(r)['sampled_deficit']-.303,.89,.90,xtol=2e-7)
    write_json(OUT/'inverse_legacy_threshold.json',dict(numerical_rho=legacy_threshold,check=trial(legacy_threshold,True),quadrature='Original one-second grid including 0 and 1e-6 s'))
    for rho in [0.,.5,.85,threshold]:
        try:res=stationary(cacc(m,rho),sc,np.array(roots[0.]['state']),g0)
        except Exception as ex:res=dict(error=str(ex),status='not_established')
        write_json(OUT/('inverse_equilibrium_'+str(round(rho,7))+'.json'),res)
    mm=cacc(m,threshold);long=integrate(mm,y0,g0,end=3600.,max_step=5.)
    lg=[]
    for tt in np.arange(0.,3601.,5.):lg.append(observe(mm,tt,long.sol(tt)[:-1],g0))
    write_csv(OUT/'inverse_long_time.csv',lg)
    crossing=None
    for a,b in zip(lg,lg[1:]):
        if a['ph_i']<7.3<=b['ph_i']:
            crossing=brentq(lambda tt:observe(mm,tt,long.sol(tt)[:-1],g0)['ph_i']-7.3,a['time_s'],b['time_s']);break
    write_json(OUT/'inverse_long_time_summary.json',dict(ph_gate_crossing_s=crossing,end_ph=lg[-1]['ph_i'],end_s=3600.,rho=threshold,protocol='constant stimulus extension, not Task 41 validation'))
    nb=[];seed=np.array(roots[1.]['state'])
    for fraction in [1.,.75,.5,.25,0.]:
        try:
            nm=clone(stim,nbc_parameters=replace(m.nbc_parameters,capacity_fmol_s=m.nbc_parameters.capacity_fmol_s*fraction))
            rr=stationary(nm,sc,seed);seed=np.array(rr['state']);nb.append(dict(NBC_fraction=fraction,physiological=rr['physiological'],failed_gates=';'.join(rr['failed_gates']),**rr['observables']))
        except Exception as ex:nb.append(dict(NBC_fraction=fraction,error=str(ex)))
    write_csv(OUT/'nbc_continuation.csv',nb)
    assert len(results)==26 and not failures
    assert all(r['jacobian_relative_step_error']<1e-4 for r in results)
    assert all(r['validation_relative_error']<1e-4 for r in sens)
    assert all(r['FD_relative_error']<1e-4 for r in pars)
    write_json(OUT/'execution_summary.json',dict(status='complete',elapsed_s=perf_counter()-start,equilibrium_rows=len(results),
      root_failures=len(failures),expression_sensitivities=len(sens),parameter_sensitivities=len(pars),inverse_integrations=len(cache),
      production_sources_modified=False,frozen_results_modified=False,numpy=np.__version__,scipy=__import__('scipy').__version__))
    print('COMPLETE',perf_counter()-start,flush=True)
if __name__=='__main__':main()
