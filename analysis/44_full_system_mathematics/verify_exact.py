"""Independent algebra, closure and dimensional RHS checks for Task 44.

Only the active parameter objects and named kinetic constants are shared with
production. The chemistry, electrical solve, flux laws and amount equations
below are independently assembled for checking the equations in exact_results.tex.
"""
from core import *
from fractions import Fraction
from scipy.optimize import brentq
from modern_full_model.nhe1_cha2009 import published_cha_kinetics
from modern_full_model.nkcc1_palk2010 import A1, A2_MM, A3, A4_MM

OUT = Path(__file__).resolve().parent / 'output'


def chemistry(T, A, V, buffer, pka, ab):
    def pieces(ph):
        h = 10.**(-ph)
        k1, k2 = 10.**(-ab.carbon_pka1), 10.**(-ab.carbon_pka2)
        z = h*h+k1*h+k1*k2
        a = np.array([h*h, k1*h, k1*k2])/z
        b = 1/(1+10.**(pka-ph))
        hm, ohm = 1000*h, 1000*10.**(ph-ab.water_pkw)
        return a, b, hm, ohm
    def residual(ph):
        a, b, hm, ohm = pieces(ph)
        return T*(a[1]+2*a[2])+buffer*b+V*(ohm-hm)-A
    ph = brentq(residual, ab.ph_lower, ab.ph_upper, xtol=1e-13)
    a, b, hm, ohm = pieces(ph)
    mean = a[1]+2*a[2]
    derivative = np.log(10)*(T*(a[1]+4*a[2]-mean**2)+buffer*b*(1-b)+V*(hm+ohm))
    return dict(ph=ph, h=hm, hco3=T/V*a[1], co2=T/V*a[0], co3=T/V*a[2],
                buffer_fraction=b, mean_charge=mean, water_charge=ohm-hm,
                Gp=derivative, ph_amount_derivatives=np.array([-mean, 1., -(ohm-hm)])/derivative)


def independent(model, y, nu, genotype=WT):
    p, nb, a4 = model.parameters, model.nbc_parameters, model.ae4_parameters
    mp, hp, wp, ab = p.membranes, p.homeostasis, p.water, p.acid_base
    ct = p.constants
    ni, nl, vi, vl = y[:5], y[6:11], y[5], y[11]
    ci, cl = ni/vi, nl/vl
    ai = chemistry(ni[3], ni[4], vi, p.geometry.cell_buffer_total_fmol, ab.cell_buffer_pka, ab)
    al = chemistry(nl[3], nl[4], vl, p.geometry.lumen_buffer_total_fmol, ab.lumen_buffer_pka, ab)
    bath = p.bath
    hh = 10.**(-bath.ph); ka1, ka2 = 10.**(-ab.carbon_pka1), 10.**(-ab.carbon_pka2)
    zz = hh*hh+ka1*hh+ka1*ka2
    be = dict(na=bath.na_mM, k=bath.k_mM, cl=bath.cl_mM,
              hco3=bath.tic_mM*ka1*hh/zz, co2=bath.tic_mM*hh*hh/zz, h=1000*hh)
    stimulus = model.stimulus(1.)
    ca = stimulus.calcium_uM
    u = np.clip((ca-nb.resting_calcium_uM)/(nb.fully_recruited_calcium_uM-nb.resting_calcium_uM), 0, 1)
    rg = model.regulatory_model
    kn = rg.nkcc1_regulatory_model
    un = np.clip((ca-kn.resting_calcium_uM)/(kn.stimulated_calcium_uM-kn.resting_calcium_uM), 0, 1)
    mn = 1+(kn.fully_activated_multiplier-1)*un
    mh = 1+(nb.nhe1_stimulated_multiplier-1)*u
    product = ci[0]*ci[1]*ci[2]**2
    N = genotype.nkcc1_expression*model.nkcc1_kinetics.alpha_eff_fmol_s*mn*(A1-A2_MM*product)/(A3+A4_MM*product)
    k = published_cha_kinetics()
    occ = lambda c,K: c/(c+K)
    aa = k.k1_plus_per_ms*occ(be['na'],k.k_na_o_mM)*(1-occ(be['h'],k.k_h_o_mM))
    bb = k.k2_plus_per_ms*occ(ai['h'],k.k_h_i_mM)*(1-occ(ci[0],k.k_na_i_mM))
    cc = k.k1_minus_per_ms*occ(ci[0],k.k_na_i_mM)*(1-occ(ai['h'],k.k_h_i_mM))
    dd = k.k2_minus_per_ms*occ(be['h'],k.k_h_o_mM)*(1-occ(be['na'],k.k_na_o_mM))
    modifier = 1/(1+(1+be['h']/k.modifier_k_o_mM)*(k.modifier_k_i_mM/ai['h'])**3)
    H = hp.nhe1_cha_carrier_amount_fmol*genotype.nhe1_expression*mh*1000*modifier*(aa*bb-cc*dd)/(aa+bb+cc+dd)
    E = genotype.ae2_expression*hp.ae2_capacity_fmol_s*np.tanh(np.log(be['cl']*ai['hco3']/(ci[2]*be['hco3']))/hp.thermodynamic_saturation_log_width)
    lc = np.log(be['cl']/ci[2])
    ln = np.log(ci[0]/be['na'])+2*np.log(ai['hco3']/be['hco3'])
    lk = np.log(ci[1]/be['k'])+2*np.log(ai['hco3']/be['hco3'])
    cf, cr = a4.common_cl_attempt_rate_s*np.exp(np.array([lc, -lc])/2)
    nf, nr = a4.na_loaded_attempt_rate_s*np.exp(np.array([ln, -ln])/2)
    kf, kr = a4.k_loaded_attempt_rate_s*np.exp(np.array([lk, -lk])/2)
    sigma = cf+cr+nf+nr+kf+kr
    pi = (cf+nr+kr)/sigma
    po = (cr+nf+kf)/sigma
    gain = rg.ae4_regulatory_model.gain
    r = y[12]
    capacity_gain = gain.basal_capacity_multiplier+gain.fully_activated_increment*gain.coupling_scale*r
    J4 = genotype.ae4_expression*a4.carrier_amount_fmol*capacity_gain*((nf+kf)*pi-(nr+kr)*po)
    pump = lambda fraction, ko: mp.nak_capacity_fmol_s*fraction*ci[0]**3/(ci[0]**3+mp.nak_na_half_mM**3)*ko**2/(ko**2+mp.nak_k_half_mM**2)
    Pa, Pb = pump(mp.apical_pump_fraction,cl[1]), pump(1-mp.apical_pump_fraction,be['k'])
    gate = ca**mp.calcium_hill/(ca**mp.calcium_hill+mp.calcium_half_uM**mp.calcium_hill)
    gka = mp.g_k_total_S*gate*mp.apical_k_fraction+mp.g_apical_background_S
    gkb = mp.g_k_total_S*gate*(1-mp.apical_k_fraction)+mp.g_basolateral_background_S
    gca = mp.g_cl_apical_S*gate
    gs = np.array([mp.g_para_na_S,mp.g_para_k_S,mp.g_para_cl_S,mp.g_para_hco3_S])
    valences = np.array([1,1,-1,-1])
    vt, conversion = ct.thermal_voltage_V, 1e15/ct.faraday_C_mol
    eka, ekb, eca = vt*np.log(cl[1]/ci[1]), vt*np.log(be['k']/ci[1]), -vt*np.log(cl[2]/ci[2])
    epara = vt/valences*np.log(np.array([be['na'],be['k'],be['cl'],be['hco3']])/np.r_[cl[:3],al['hco3']])
    ga, gp = gka+gca, sum(gs)
    caa, cpa = -gka*eka-gca*eca+Pa/conversion, -np.dot(gs,epara)
    va = lambda vb: (gp*vb-caa+cpa)/(ga+gp)
    ell = np.log(be['na']*be['hco3']**2/(ci[0]*ai['hco3']**2))
    nbc = lambda vb: nb.capacity_fmol_s*u*np.tanh((ell+vb/vt)/nb.log_width)
    balance = lambda vb: gkb*(vb-ekb)+Pb/conversion+nbc(vb)/conversion+gp*(vb-va(vb))+cpa
    vb = brentq(balance,-.5,.5,xtol=1e-14,rtol=1e-13)
    v_a = va(vb)
    B = nbc(vb)
    ka = conversion*gka*(v_a-eka)
    kb = conversion*gkb*(vb-ekb)
    cla = -conversion*gca*(v_a-eca)
    jp = conversion*gs*(vb-v_a-epara)/valences
    cb = hp.co2_basolateral_permeability_fmol_s_mM*(be['co2']-ai['co2'])
    c_a = hp.co2_apical_permeability_fmol_s_mM*(al['co2']-ai['co2'])
    oi = sum(ci[:4])+(p.geometry.cell_impermeant_osmoles_fmol+p.geometry.cell_buffer_total_fmol)/vi
    ol = sum(cl[:4]); oe = be['na']+be['k']+be['cl']+bath.tic_mM+bath.untracked_osmolyte_mM
    qa = wp.apical_hydraulic_pL_s_mOsm*(ol-oi)
    qb = wp.basolateral_hydraulic_pL_s_mOsm*(oi-oe)
    qt = wp.paracellular_hydraulic_pL_s_mOsm*(ol-oe)
    qo = wp.outflow_rate_s*max(vl-wp.lumen_dead_volume_pL,0)
    rhs = np.zeros(13)
    rhs[:5] = np.array([N+H+B-3*(Pa+Pb),N+2*(Pa+Pb)-ka-kb,2*N+E-cla,-E+2*B+cb+c_a,H-E+2*B])+J4*np.array(nu)
    rhs[5] = qb-qa
    rhs[6:11] = np.array([3*Pa-jp[0],-2*Pa+ka-jp[1],cla-jp[2],-jp[3]-c_a,-jp[3]])-qo*cl
    rhs[11] = qa+qt-qo
    rhs[12] = (stimulus.beta_input-r)/rg.ae4_regulatory_model.tau_activation_s
    slope = gkb+gp*ga/(gp+ga)+nb.capacity_fmol_s*u/(conversion*nb.log_width*vt)/np.cosh((ell+vb/vt)/nb.log_width)**2
    return rhs, dict(ai=ai,al=al,N=N,H=H,E=E,B=B,J4=J4,Pa=Pa,Pb=Pb,CaCC=cla,
                      voltage_slope=slope, voltage_slope_fd=(balance(vb+1e-6)-balance(vb-1e-6))/2e-6,
                      passive_slope=gkb+gp*ga/(gp+ga),Va=v_a,Vb=vb,
                      voltage_bracket=[balance(-.5),balance(.5)],carrier_time_s=1/sigma)


def main():
    _, stim, y0, _, audit = reference()
    m = clone(stim); sc = Scaling(m,y0); rng = np.random.default_rng(45044)
    maxima = {name:0. for name in ['rhs_scaled','carbonate_ph','carbonate_derivative_relative','electrical_voltage','electrical_slope_relative','structural_identity','carbon_alkalinity_identity','charge_rate']}
    min_slope = np.inf; phs=[]; rows=[]; count=0
    for name, mechanism in CLASSES.items():
        mm = with_mechanism_class(m,name)
        nu = mechanism.conserved_coefficients
        for sample in range(12):
            z = sc.independent(y0); z[:10] *= np.exp(rng.normal(0,.008,10)); z[-1]=rng.uniform(.1,.9)
            y = sc.expand(z); gg = Genotype('check',ae4_expression=float(rng.uniform(0,1)))
            rr,d = independent(mm,y,nu,gg)
            canonical = mm.evaluate(1.,y,genotype=gg)
            maxima['rhs_scaled'] = max(maxima['rhs_scaled'],float(np.max(np.abs((rr-canonical.rhs)*sc.t0/sc.D))))
            maxima['carbonate_ph'] = max(maxima['carbonate_ph'],abs(d['ai']['ph']-canonical.diagnostics.observables.cell_acid_base.ph),abs(d['al']['ph']-canonical.diagnostics.observables.lumen_acid_base.ph))
            maxima['electrical_voltage'] = max(maxima['electrical_voltage'],abs(d['Va']-canonical.diagnostics.membranes.v_apical_V),abs(d['Vb']-canonical.diagnostics.membranes.v_basolateral_V))
            maxima['electrical_slope_relative'] = max(maxima['electrical_slope_relative'],abs(d['voltage_slope_fd']/d['voltage_slope']-1))
            min_slope=min(min_slope,d['voltage_slope']); phs.append(d['ai']['ph'])
            assert d['voltage_bracket'][0]<0<d['voltage_bracket'][1]
            kappa = nu[2]-2*nu[0]+nu[4]
            identity = 6*(d['Pa']+d['Pb'])-d['H']+kappa*d['J4']+2*rr[0]-rr[4]-rr[2]
            maxima['structural_identity']=max(maxima['structural_identity'],abs(identity-d['CaCC']))
            maxima['carbon_alkalinity_identity']=max(maxima['carbon_alkalinity_identity'],abs(rr[4]-(d['H']-d['E']+2*d['B']+nu[4]*d['J4'])))
            maxima['charge_rate']=max(maxima['charge_rate'],abs(rr[0]+rr[1]-rr[2]-rr[4]),abs(rr[6]+rr[7]-rr[8]-rr[10]))
            for offset,buffer,pka in [(0,m.parameters.geometry.cell_buffer_total_fmol,m.parameters.acid_base.cell_buffer_pka),(6,m.parameters.geometry.lumen_buffer_total_fmol,m.parameters.acid_base.lumen_buffer_pka)]:
                point=np.array([y[offset+3],y[offset+4],y[offset+5]])
                closure=chemistry(*point,buffer,pka,m.parameters.acid_base)
                # Volume derivative is tiny because common amount dilution cancels;
                # use larger volume steps and absolute error scaled by 1/Gp.
                for h in [1e-5,5e-6]:
                    deriv=[]
                    for j in range(3):
                        delta=np.zeros(3);delta[j]=h*point[j]
                        a=chemistry(*(point+delta),buffer,pka,m.parameters.acid_base)['ph']
                        b=chemistry(*(point-delta),buffer,pka,m.parameters.acid_base)['ph']
                        deriv.append((a-b)/(2*delta[j]))
                    error=max(abs(np.array(deriv)-closure['ph_amount_derivatives']))*closure['Gp']
                    maxima['carbonate_derivative_relative']=max(maxima['carbonate_derivative_relative'],float(error))
            count += 1
        nf=tuple(map(Fraction,nu)); charge=nf[0]+nf[1]-nf[2]-nf[4]
        kap=nf[2]-2*nf[0]+nf[4]
        assert charge==0
        rows.append(dict(mechanism=name,nu=[str(v) for v in nf],chloride=str(nf[2]),carbon=str(nf[3]),alkalinity=str(nf[4]),charge=str(charge),stationary_projection=str(kap)))
    assert maxima['rhs_scaled']<1e-9
    assert maxima['carbonate_ph']<1e-10 and maxima['carbonate_derivative_relative']<1e-6
    assert maxima['electrical_voltage']<1e-10 and maxima['electrical_slope_relative']<1e-7
    assert maxima['structural_identity']<1e-12 and maxima['charge_rate']<1e-11
    _,d0=independent(m,y0,CLASSES['C0'].conserved_coefficients)
    zz=sc.independent(y0);zz[-1]=1.
    reg_errors=[]
    for step in [1e-5,5e-6]:
        jr=jacobian(lambda z:sc.reduced_rhs(m,z,WT),zz,step=step)/sc.t0
        reg_errors.append(abs(jr[-1,-1]+1/30))
    assert max(reg_errors)<1e-10
    geom=m.parameters.geometry; ab=m.parameters.acid_base; ct=m.parameters.constants
    groups=dict(cell_buffer_scaled=geom.cell_buffer_total_fmol/sc.n0,lumen_buffer_scaled=geom.lumen_buffer_total_fmol/sc.nl0,
      fixed_cell_charge_scaled=geom.fixed_cell_anion_equivalents_fmol/sc.n0,
      cell_other_impermeant_scaled=geom.cell_impermeant_osmoles_fmol/sc.n0,
      Ka1_over_C0=1000*10.**(-ab.carbon_pka1)/sc.C0,Ka2_over_C0=1000*10.**(-ab.carbon_pka2)/sc.C0,
      Kw_over_C0_squared=1e6*10.**(-ab.water_pkw)/sc.C0**2,
      cell_buffer_K_over_C0=1000*10.**(-ab.cell_buffer_pka)/sc.C0,
      lumen_buffer_K_over_C0=1000*10.**(-ab.lumen_buffer_pka)/sc.C0,
      co2_basolateral=m.parameters.homeostasis.co2_basolateral_permeability_fmol_s_mM*sc.C0/sc.J0,
      co2_apical=m.parameters.homeostasis.co2_apical_permeability_fmol_s_mM*sc.C0/sc.J0,
      AE4_to_NKCC_cycle_initial=d0['J4']/d0['N'],AE4_to_NKCC_chloride_initial=d0['J4']/(2*d0['N']),
      carrier_QSS_time_ratio_initial=d0['carrier_time_s']/sc.t0,
      electrical_relaxation_group=None,
      electrical_relaxation_reason='No capacitance or voltage differential equation in the frozen model; conductance groups measure currents and do not establish an electrical relaxation time.')
    # Exact projections all five jointly have full rank: Cl, T, A, q, kappa.
    # Their common kernel is {0}, since Cl=T=A=0, q=Na+K=0,
    # kappa=-2Na=0 imply Na=K=0.
    result=dict(status='passed',random_seed=45044,states=count,states_per_class=12,
      maxima=maxima,minimum_electrical_slope_S=float(min_slope),regulatory_derivative_error_per_s=reg_errors,
      intracellular_ph_range=[float(min(phs)),float(max(phs))],source_vectors=rows,
      joint_projection_rank=5,joint_projection_kernel='zero',dimensionless_additional_groups=groups,
      regulatory_jacobian_caveat='At r=1 the executable trial clip is not differentiable. A centred difference halves its column. Use the inward derivative or the stated smooth analytical extension dr/dt=(1-r)/30 for full spectrum and modal calculations.',
      proof_scope='Exact identities and conditional closure proofs; random checks establish implementation agreement, not global physiology or global stability.')
    write_json(OUT/'exact_verification.json',result)
    print('EXACT_CHECKS passed:',count,'states; scaled RHS error',format(maxima['rhs_scaled'],'.3g'))


if __name__=='__main__':
    main()
