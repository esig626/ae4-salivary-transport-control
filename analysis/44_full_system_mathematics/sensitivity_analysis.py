"""Task 45: local uncertainty crosswalk on the unchanged full Task 40 system.

Every numerical perturbation is an analysis copy. No physiological uncertainty
interval is inferred from a finite difference step or from root admissibility.
"""
from core import *
from modern_full_model.nhe1_cha2009 import published_cha_kinetics

OUT = Path(__file__).resolve().parent / 'output'
STEPS = (1e-3, 5e-4)
OBS = ('q_pL_s', 'ph_i', 'N', 'J4', 'P', 'H', 'B', 'E')


def parameter_builders(stim, m):
    """Selected uncertain effective parameters, agreed with the Task 43 audit."""
    p, mp, wp = m.parameters, m.parameters.membranes, m.parameters.water
    reg = m.regulatory_model
    ae = reg.ae4_regulatory_model
    rows = {}

    def add(name, path, value, build, experiment):
        rows[name] = dict(inventory_path=path, value=float(value), build=build,
                          experiment=experiment)

    add('NBC_capacity', 'nbc.capacity_fmol_s', m.nbc_parameters.capacity_fmol_s,
        lambda h: clone(stim, nbc_parameters=replace(m.nbc_parameters, capacity_fmol_s=m.nbc_parameters.capacity_fmol_s*np.exp(h))),
        'Measure basolateral sodium coupled alkalinity influx at controlled intracellular pH, sodium and bicarbonate.')
    add('AE4_carrier', 'ae4.carrier_amount_fmol', m.ae4_parameters.carrier_amount_fmol,
        lambda h: clone(stim, ae4_parameters=replace(m.ae4_parameters, carrier_amount_fmol=m.ae4_parameters.carrier_amount_fmol*np.exp(h))),
        'Measure functional AE4 chloride turnover and carrier abundance under the same stimulation protocol.')
    add('NKCC_scale', 'nkcc1.alpha_eff_fmol_s', m.nkcc1_kinetics.alpha_eff_fmol_s,
        lambda h: clone(stim, nkcc1_kinetics=replace(m.nkcc1_kinetics, alpha_eff_fmol_s=m.nkcc1_kinetics.alpha_eff_fmol_s*np.exp(h))),
        'Measure NKCC chloride uptake together with intracellular sodium, potassium and chloride in WT and AE4 null cells.')
    add('NHE_carrier', 'parameters.homeostasis.nhe1_cha_carrier_amount_fmol', p.homeostasis.nhe1_cha_carrier_amount_fmol,
        lambda h: clone(stim, parameters=replace(p, homeostasis=replace(p.homeostasis, nhe1_cha_carrier_amount_fmol=p.homeostasis.nhe1_cha_carrier_amount_fmol*np.exp(h)))),
        'Measure acid extrusion versus intracellular pH in the salivary preparation, separating basal and stimulated capacity.')
    add('pump_capacity', 'parameters.membranes.nak_capacity_fmol_s', mp.nak_capacity_fmol_s,
        lambda h: clone(stim, parameters=replace(p, membranes=replace(mp, nak_capacity_fmol_s=mp.nak_capacity_fmol_s*np.exp(h)))),
        'Measure total sodium potassium pump turnover and its response to intracellular sodium.')
    add('CaCC_conductance', 'parameters.membranes.g_cl_apical_S', mp.g_cl_apical_S,
        lambda h: clone(stim, parameters=replace(p, membranes=replace(mp, g_cl_apical_S=mp.g_cl_apical_S*np.exp(h)))),
        'Measure apical chloride conductance with calcium, voltage and intracellular chloride controlled.')
    add('buffer_pool', 'parameters.geometry.cell_buffer_total_fmol', p.geometry.cell_buffer_total_fmol,
        lambda h: clone(stim, parameters=replace(p, geometry=replace(p.geometry, cell_buffer_total_fmol=p.geometry.cell_buffer_total_fmol*np.exp(h)))),
        'Measure intrinsic buffer capacity with simultaneous cell volume, carbonate and intracellular pH.')
    add('water_apical', 'parameters.water.apical_hydraulic_pL_s_mOsm', wp.apical_hydraulic_pL_s_mOsm,
        lambda h: clone(stim, parameters=replace(p, water=replace(wp, apical_hydraulic_pL_s_mOsm=wp.apical_hydraulic_pL_s_mOsm*np.exp(h)))),
        'Measure apical water permeability and cell and lumen volume responses to controlled osmotic steps.')
    add('AE4_tau', 'regulation.ae4_regulatory_model.tau_activation_s', ae.tau_activation_s,
        lambda h: clone(stim, regulation=replace(reg, ae4_regulatory_model=replace(ae, tau_activation_s=ae.tau_activation_s*np.exp(h)))),
        'Measure the time course of AE4 activation separately from whole cell storage and secretion.')
    add('AE4_gain', 'regulation.ae4_regulatory_model.gain.fully_activated_increment', ae.gain.fully_activated_increment,
        lambda h: clone(stim, regulation=replace(reg, ae4_regulatory_model=replace(ae, gain=replace(ae.gain, fully_activated_increment=ae.gain.fully_activated_increment*np.exp(h))))),
        'Resolve the discrepant AE4 activation gain readings under matched experimental conditions.')
    add('apical_pump_fraction', 'parameters.membranes.apical_pump_fraction', mp.apical_pump_fraction,
        lambda h: clone(stim, parameters=replace(p, membranes=replace(mp, apical_pump_fraction=mp.apical_pump_fraction*np.exp(h)))),
        'Quantify functional pump turnover at each membrane rather than inferring a fraction from localisation alone.')
    add('NHE_stimulated_multiplier', 'nbc.nhe1_stimulated_multiplier', m.nbc_parameters.nhe1_stimulated_multiplier,
        lambda h: clone(stim, nbc_parameters=replace(m.nbc_parameters, nhe1_stimulated_multiplier=m.nbc_parameters.nhe1_stimulated_multiplier*np.exp(h))),
        'Measure the stimulation dependent change in NHE acid extrusion in salivary cells.')
    add('CO2_basolateral_permeability', 'parameters.homeostasis.co2_basolateral_permeability_fmol_s_mM', p.homeostasis.co2_basolateral_permeability_fmol_s_mM,
        lambda h: clone(stim, parameters=replace(p, homeostasis=replace(p.homeostasis, co2_basolateral_permeability_fmol_s_mM=p.homeostasis.co2_basolateral_permeability_fmol_s_mM*np.exp(h)))),
        'Measure intracellular carbonate and pH relaxation after a controlled basolateral carbon dioxide step.')
    add('AE2_capacity', 'parameters.homeostasis.ae2_capacity_fmol_s', p.homeostasis.ae2_capacity_fmol_s,
        lambda h: clone(stim, parameters=replace(p, homeostasis=replace(p.homeostasis, ae2_capacity_fmol_s=p.homeostasis.ae2_capacity_fmol_s*np.exp(h)))),
        'Measure AE2 chloride bicarbonate exchange capacity and direction over the observed pH range.')
    add('NBC_log_width', 'nbc.log_width', m.nbc_parameters.log_width,
        lambda h: clone(stim, nbc_parameters=replace(m.nbc_parameters, log_width=m.nbc_parameters.log_width*np.exp(h))),
        'Measure the shape of the sodium bicarbonate flux versus chemical driving force relation.')
    add('NKCC_stimulated_multiplier', 'regulation.nkcc1_regulatory_model.fully_activated_multiplier', reg.nkcc1_regulatory_model.fully_activated_multiplier,
        lambda h: clone(stim, regulation=replace(reg, nkcc1_regulatory_model=replace(reg.nkcc1_regulatory_model, fully_activated_multiplier=reg.nkcc1_regulatory_model.fully_activated_multiplier*np.exp(h)))),
        'Measure the NKCC stimulation gain in the same calcium and secretagogue protocol used for secretion.')
    return rows


def quantities(m, sc, z, g):
    o = observe(m, 1., sc.expand(z), g)
    return np.array([o[k] for k in OBS])


def expression_gain(m, sc, rr, validate=False):
    """IFT chloride compensation, with independent nearby expression roots."""
    z = np.array(rr['z'])
    j = np.array(rr['jacobian_dimensionless'])
    qx = jacobian(lambda a: quantities(m, sc, a, WT), z, 5e-6)
    rows = []
    for h in (1e-4, 5e-5):
        gp, gm = Genotype('plus', ae4_expression=1+h), Genotype('minus', ae4_expression=1-h)
        fp = (sc.reduced_rhs(m, z, gp)-sc.reduced_rhs(m, z, gm))/(2*h)
        dz = np.linalg.solve(j, -fp)
        direct = (quantities(m, sc, z, gp)-quantities(m, sc, z, gm))/(2*h)
        total = direct+qx@dz
        rec = dict(expression_step=h, gain=float(-2*total[2]/total[3]),
                   dN_de=float(total[2]), dJ4_de=float(total[3]))
        if validate:
            rp = stationary(m, sc, np.array(rr['state']), gp)
            rm = stationary(m, sc, np.array(rr['state']), gm)
            fd = (quantities(m, sc, np.array(rp['z']), gp)-quantities(m, sc, np.array(rm['z']), gm))/(2*h)
            rec.update(FD_gain=float(-2*fd[2]/fd[3]),
                       FD_relative_error=float(np.linalg.norm(total-fd)/max(np.linalg.norm(fd), 1e-12)))
        rows.append(rec)
    return rows[-1]['gain'], rows


def nhe_physiological_upper(m):
    """Upper bound for pH_i >= 6.6, positive sodium and the fixed bath.

    In Cha Eq. 3, c >= 0, b <= bmax, and (ab-cd)/(a+b+c+d)
    <= a*bmax/(a+bmax+d). Mod2 is at most its value at pH_i=6.6.
    This bounds influx even when exchange is negative.
    """
    k = published_cha_kinetics()
    no = m.parameters.bath.na_mM
    ho = 1e3*10**(-m.parameters.bath.ph)
    hi = 1e3*10**(-6.6)
    a = k.k1_plus_per_ms*no/(no+k.k_na_o_mM)*k.k_h_o_mM/(ho+k.k_h_o_mM)
    b = k.k2_plus_per_ms*hi/(hi+k.k_h_i_mM)
    d = k.k2_minus_per_ms*ho/(ho+k.k_h_o_mM)*k.k_na_o_mM/(no+k.k_na_o_mM)
    mod = 1/(1+(1+ho/k.modifier_k_o_mM)*(k.modifier_k_i_mM/hi)**3)
    return m.parameters.homeostasis.nhe1_cha_carrier_amount_fmol*m.nbc_parameters.nhe1_stimulated_multiplier*1000*mod*a*b/(a+b+d)


def resolved_decay(m, sc, rr, g):
    """Fourth order chemical differentiation reduces noise in tiny elasticities.

    The full Jacobian is block upper triangular: the regulatory equation depends
    only on its own state. Its exact eigenvalue is retained, never frozen away.
    """
    z = np.array(rr['z'])
    estimates = []
    for step in (2e-4, 1e-4):
        cols = []
        for k in range(10):
            d = np.zeros(11)
            d[k] = step*max(abs(z[k]), .1)
            f = lambda a: sc.reduced_rhs(m, a, g)[:10]
            cols.append((f(z-2*d)-8*f(z-d)+8*f(z+d)-f(z+2*d))/(12*d[k]))
        chemical = np.linalg.eigvals(np.stack(cols, axis=-1)/sc.t0)
        regulatory = -1/m.regulatory_model.ae4_regulatory_model.tau_activation_s
        dominant = max(max(chemical.real), regulatory)
        assert dominant < 0
        estimates.append(float(-1/dominant))
    return estimates[-1], dict(fourth_order_steps=[2e-4, 1e-4], decay_estimates_s=estimates,
      relative_step_disagreement=abs(estimates[0]-estimates[1])/estimates[1],
      relative_difference_from_second_order=abs(estimates[-1]-rr['slowest_decay_s'])/estimates[-1])


def combined_metrics(m, sc, wt, null, validate_gain=False):
    o, n = wt['observables'], null['observables']
    gain, checks = expression_gain(m, sc, wt, validate_gain)
    hn = nhe_physiological_upper(m)
    upper = (hn+m.parameters.homeostasis.ae2_capacity_fmol_s)/2
    tw, tw_check = resolved_decay(m, sc, wt, WT)
    tn, tn_check = resolved_decay(m, sc, null, Genotype('null', ae4_expression=0.))
    checks.append(dict(eigenvalue_validation={'WT':tw_check, 'null':tn_check}))
    assert o['J4'] > 0
    vals = dict(WT_Q=o['q_pL_s'], WT_pH=o['ph_i'],
                null_stationary_deficit=1-n['q_pL_s']/o['q_pL_s'],
                NBC_alkalinity_share=o['B']/o['J4'],
                noNBC_observed_support_share=(o['H']-o['E'])/(2*o['J4']),
                AE4_alkalinity_demand=2*o['J4'],
                NHE_physiological_upper=hn,
                noNBC_AE4_capacity_upper=upper,
                noNBC_capacity_to_current_demand=upper/o['J4'],
                noNBC_capacity_gap=o['J4']-upper,
                chloride_compensation_gain=gain,
                WT_slowest_decay=tw,
                null_slowest_decay=tn)
    return {k: float(v) for k, v in vals.items()}, checks


def main():
    _, stim, y0, _, _ = reference()
    m = clone(stim)
    sc = Scaling(m, y0)
    roots = json.loads((OUT/'equilibria.json').read_text())['roots']
    seed = {r['ae4_expression']: r for r in roots if r['mechanism']=='C0'}
    wt = stationary(m, sc, np.array(seed[1.]['state']))
    gn = Genotype('null', ae4_expression=0.)
    null = stationary(m, sc, np.array(seed[0.]['state']), gn)
    baseline, gain_checks = combined_metrics(m, sc, wt, null, True)
    selections = parameter_builders(stim, m)
    detail, flat = [], []
    for name, selection in selections.items():
        build = selection['build']
        by_step, root_checks = [], []
        for h in STEPS:
            models = [build(h), build(-h)]
            perturbed = []
            for sign, mm in zip((1, -1), models):
                rw = stationary(mm, sc, np.array(wt['state']))
                rn = stationary(mm, sc, np.array(null['state']), gn)
                metrics, cg = combined_metrics(mm, sc, rw, rn, True)
                perturbed.append((rw, rn, metrics))
                for genotype, rr in [('WT', rw), ('null', rn)]:
                    root_checks.append(dict(log_step=sign*h, genotype=genotype,
                      residual_max=rr['residual_max'], physiological=rr['physiological'],
                      failed_gates=rr['failed_gates'], jacobian_step_error=rr['jacobian_relative_step_error'],
                      condition_number=rr['condition_number']))
                root_checks.append(dict(log_step=sign*h, genotype='WT_expression_gain', checks=cg))
            dw = []
            dn = []
            validation = []
            # Fixed-state buffer terms nearly cancel against the state response.
            # Use a smaller derivative step than the independent root displacement.
            derivative_h = h/100
            derivative_models = [build(derivative_h), build(-derivative_h)]
            for genotype, rr, g, index in [('WT', wt, WT, 0), ('null', null, gn, 1)]:
                z = np.array(rr['z'])
                j = np.array(rr['jacobian_dimensionless'])
                fp = (sc.reduced_rhs(derivative_models[0], z, g)-sc.reduced_rhs(derivative_models[1], z, g))/(2*derivative_h)
                dz = np.linalg.solve(j, -fp)
                qx = jacobian(lambda a: quantities(m, sc, a, g), z, 5e-6)
                direct = (quantities(derivative_models[0], sc, z, g)-quantities(derivative_models[1], sc, z, g))/(2*derivative_h)
                total = direct+qx@dz
                fd = (quantities(models[0], sc, np.array(perturbed[0][index]['z']), g)-quantities(models[1], sc, np.array(perturbed[1][index]['z']), g))/(2*h)
                absolute = float(np.max(abs(total-fd)))
                relative = float(np.linalg.norm(total-fd)/max(np.linalg.norm(fd), 1e-8))
                validation.append(dict(genotype=genotype, IFT_log_parameter_step=derivative_h, observables=list(OBS),
                  IFT_derivatives=total.tolist(), independently_solved_FD_derivatives=fd.tolist(),
                  absolute_max_error=absolute, relative_error=relative))
                if genotype=='WT': dw=total
                else: dn=total
            deriv = {k: (perturbed[0][2][k]-perturbed[1][2][k])/(2*h) for k in baseline}
            # Explicit IFT results for secretion, pH and the derived knockout effect.
            deriv['WT_Q_IFT'] = float(dw[0])
            deriv['WT_pH_IFT'] = float(dw[1])
            deriv['null_stationary_deficit_IFT'] = float(-dn[0]/baseline['WT_Q']+null['observables']['q_pL_s']*dw[0]/baseline['WT_Q']**2)
            by_step.append(dict(log_parameter_step=h, derivatives=deriv, IFT_validation=validation))
        small, large = by_step[-1]['derivatives'], by_step[0]['derivatives']
        row = dict(parameter=name, inventory_path=selection['inventory_path'], active_value=selection['value'],
                   interval_status='No justified physiological uncertainty interval; local analysis only.')
        for k, b in baseline.items():
            row[k+'_per_log_parameter']=float(small[k])
            row[k+'_log_elasticity']=float(small[k]/b) if abs(b)>1e-12 else None
            row[k+'_elasticity_step_disagreement']=float(abs(small[k]-large[k])/abs(b)) if abs(b)>1e-12 else None
        row.update(WT_pH_per_log_parameter=float(small['WT_pH_IFT']),
                   WT_H_concentration_log_elasticity=float(-np.log(10)*small['WT_pH_IFT']),
                   null_deficit_percentage_points_per_log_parameter=float(100*small['null_stationary_deficit_IFT']),
                   step_absolute_max_error=float(max(abs(small[k]-large[k]) for k in baseline)),
                   max_IFT_relative_error=float(max(v['relative_error'] for st in by_step for v in st['IFT_validation'])),
                   all_nearby_roots_physiological=all(r.get('physiological', True) for r in root_checks))
        detail.append(dict(**{k:v for k,v in selection.items() if k!='build'}, parameter=name,
                           results=row, by_step=by_step, nearby_root_checks=root_checks))
        flat.append(row)
        print('SENSITIVITY', name, 'Q', f"{row['WT_Q_log_elasticity']:.5g}",
              'deficit', f"{row['null_stationary_deficit_log_elasticity']:.5g}", flush=True)
    result = dict(scope='Local sensitivities on constant stimulus equilibria of the full retained system. No fitted parameters or inferred uncertainty intervals.',
      independent_log_parameter_steps=list(STEPS), baseline=baseline,
      baseline_expression_gain_validation=gain_checks, parameters=detail,
      limitations=['Finite difference step sizes are numerical checks, not admissible uncertainty ranges.',
        'Stationary total alkalinity support equals one by balance, so robustness is reported for its components and a conditional no NBC capacity comparison.',
        'No finite parameter change is established to overturn a conclusion by these local derivatives.',
        'Positive rates have a logical positive domain; the apical pump fraction has logical domain [0,1]. Neither is an experimentally justified uncertainty interval.',
        'Buffer sensitivity holds fixed the independent charge and impermeant osmole pools. It is an explicit partial sensitivity, not a refitted resting calibration.',
        'AE4 activation time has zero local stationary effect while its regulatory eigenvalue remains separate from the slower chemical mode.',
        'Task 41 threshold sensitivity is calculated separately in inverse_analysis.py.'])
    write_json(OUT/'sensitivity_extended.json', result)
    write_csv(OUT/'sensitivity_extended.csv', flat)
    write_fragment(result)
    verify_saved_results()


def verify_saved_results():
    """Audit saved results and check the capacity bound against production flux."""
    from modern_full_model.nhe1_cha2009 import cha_nhe1_flux_fmol_s
    r = json.loads((OUT/'sensitivity_extended.json').read_text())
    _, stim, _, _, _ = reference()
    m = clone(stim)
    rows = r['parameters']
    assert len(rows)==16
    max_ift = max(p['results']['max_IFT_relative_error'] for p in rows)
    max_step = max(v for p in rows for k,v in p['results'].items() if k.endswith('_elasticity_step_disagreement'))
    max_gain = 0.
    max_eigen = 0.
    max_residual = 0.
    nearby = 0
    for p in rows:
        assert p['results']['all_nearby_roots_physiological']
        for check in p['nearby_root_checks']:
            if 'residual_max' in check:
                max_residual = max(max_residual, check['residual_max'])
                nearby += 1
            for c in check.get('checks', []):
                if 'FD_relative_error' in c:
                    max_gain = max(max_gain, c['FD_relative_error'])
                for ev in c.get('eigenvalue_validation', {}).values():
                    max_eigen = max(max_eigen, ev['relative_step_disagreement'], ev['relative_difference_from_second_order'])
    assert max_ift < 1e-4
    assert max_step < 1e-4
    assert max_gain < 1e-4
    assert max_eigen < 1e-5
    assert max_residual < 1e-7
    rng = np.random.default_rng(45043)
    maximum_bound_ratio = 0.
    bound = nhe_physiological_upper(m)
    for _ in range(128):
        ph = rng.uniform(6.6, 7.3)
        sodium = 10**rng.uniform(-6, np.log10(40.))
        actual = cha_nhe1_flux_fmol_s(na_i_mM=sodium, h_i_mM=1e3*10**(-ph),
          na_o_mM=m.parameters.bath.na_mM, h_o_mM=1e3*10**(-m.parameters.bath.ph),
          kinetics=published_cha_kinetics(), carrier_amount_fmol=m.parameters.homeostasis.nhe1_cha_carrier_amount_fmol,
          activity_scale=m.nbc_parameters.nhe1_stimulated_multiplier)
        maximum_bound_ratio = max(maximum_bound_ratio, actual/bound)
        assert actual <= bound*(1+1e-12)
    b = r['baseline']
    assert abs(b['NBC_alkalinity_share']+b['noNBC_observed_support_share']-1) < 1e-9
    assert b['noNBC_capacity_to_current_demand'] < 1
    write_json(OUT/'sensitivity_verification.json', dict(status='passed', parameter_count=len(rows),
      independent_nearby_WT_null_roots=nearby, log_parameter_steps=list(STEPS),
      max_IFT_relative_error=max_ift, max_elasticity_step_disagreement=max_step,
      max_expression_gain_FD_relative_error=max_gain,
      max_eigenvalue_method_relative_difference=max_eigen,
      max_equilibrium_residual=max_residual, NHE_bound_random_checks=128,
      largest_sampled_NHE_to_bound_ratio=maximum_bound_ratio,
      bounds_scope='Random checks supplement the algebraic inequality; they are not its proof.',
      uncertainty_interval_inferred=False))


def write_fragment(result):
    labels = {'NBC_capacity':'NBC capacity', 'AE4_carrier':'AE4 amount', 'NKCC_scale':'NKCC scale',
      'NHE_carrier':'NHE amount', 'pump_capacity':'Pump capacity', 'CaCC_conductance':'CaCC conductance',
      'buffer_pool':'Cell buffer amount', 'water_apical':'Apical water permeability', 'AE4_tau':'AE4 activation time',
      'AE4_gain':'AE4 activation gain', 'apical_pump_fraction':'Apical pump fraction',
      'NHE_stimulated_multiplier':'NHE stimulation gain', 'CO2_basolateral_permeability':'Basolateral CO$_2$ permeability',
      'AE2_capacity':'AE2 capacity', 'NBC_log_width':'NBC saturation width', 'NKCC_stimulated_multiplier':'NKCC stimulation gain'}
    lines = [r'\subsection{Sensitivity crosswalk with Task 43}',
      r'The sixteen selected parameters are uncertain effective settings or inherited calibrated capacities. Each is perturbed independently by $p\exp(\pm h)$, with $h=10^{-3}$ and $5\times10^{-4}$; both genotypes are solved afresh at each perturbed point. These steps are numerical checks and do not define uncertainty intervals. No defensible physiological interval is available from the provenance audit, so no global robustness claim is made.',
      r'\begin{center}\small', r'\begin{tabular}{lrrrrr}\hline',
      r'Parameter & $E_Q$ & $\partial\mathrm{pH}/\partial\log p$ & $E_D$ & $E_{C_N}$ & $E_{\tau_{\rm WT}}$\\\hline']
    for rec in result['parameters']:
        r = rec['results']
        vals = [r['WT_Q_log_elasticity'], r['WT_pH_per_log_parameter'], r['null_stationary_deficit_log_elasticity'], r['chloride_compensation_gain_log_elasticity'], r['WT_slowest_decay_log_elasticity']]
        lines.append(labels[rec['parameter']]+' & '+' & '.join(f'{v:.4g}' for v in vals)+r'\\')
    lines += [r'\hline\end{tabular}\end{center}',
      r'Here $E_Y=\partial\log Y/\partial\log p$ and $D=1-Q_{\rm null}/Q_{\rm WT}$ is the stationary deficit. Its small baseline value makes $E_D$ potentially large; the machine readable table also supplies percentage point derivatives. pH is already logarithmic, so its unnormalised derivative is shown; the hydrogen concentration elasticity is $-\log(10)\partial\mathrm{pH}/\partial\log p$. All derivatives are partial sensitivities at the inherited parameter point, without recalibration.',
      r'Pump capacity has the largest local effect on WT secretion and on the stationary null deficit. The latter is also sensitive to NBC capacity and its assumed saturation width, and to effective NKCC capacity. These are priorities for experimental constraint. The calculations do not show that a finite parameter change overturns a conclusion, and cannot establish that the small null phenotype persists throughout an unknown uncertainty range. The exact balance identities survive these capacity changes within the specified architecture; their numerical consequences need not.',
      r'The equilibrium total alkalinity support ratio is identically one and cannot demonstrate robustness. The saved results instead give the NBC share $J_B/J_4$, the observed non NBC share $(J_H-J_2)/(2J_4)$, and a conditional capacity ratio $(H_{\max}+J_{2,\max})/(2J_4)$. The bound assumes the fixed bath, positive sodium and $\mathrm{pH}_i\geq6.6$. It permits AE2 to supply alkalinity at its maximal reverse capacity. Its local sensitivity is useful, but it does not prove the available capacity over an unknown parameter interval.',
      r'The AE4 activation time has no local effect on constant stimulus equilibria. Its eigenvalue is $-1/\tau_4$ and is separated from the slower chemical mode at this point; a zero local slow mode sensitivity does not establish that activation timing is irrelevant to finite time secretion. The buffer calculation holds the fixed charge and impermeant osmole pools fixed while changing only buffer amount. NKCC capacity and its stimulated multiplier produce the same stationary perturbation because only their product enters at constant full stimulation.',
      r'Independent root perturbations validate the IFT derivatives at both steps, and expression perturbations validate the chloride compensation gain at both steps for every nearby parameter point. Fixed state parameter derivatives use the smaller steps $10^{-5}$ and $5\times10^{-6}$ to resolve cancellation in the buffer response. Decay sensitivities additionally use fourth order differentiation of the chemical block at two steps, together with the exact regulatory eigenvalue; the block triangular structure retains all eleven dynamical eigenvalues. The supplied files retain root residuals, physiological gate results, Jacobian step checks and all derivative comparisons. The separate inverse family calculation supplies the Task 41 crossing sensitivity.']
    (Path(__file__).resolve().parent/'sensitivity_results.tex').write_text('\n'.join(lines)+'\n')


if __name__ == '__main__':
    main()
