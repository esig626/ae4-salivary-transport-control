"""Conservative source and equation annotations, without numerical sensitivity runs."""
from __future__ import annotations
import ast
import hashlib
import json
from pathlib import Path

# These are functions in the selected production path, not legacy alternatives.
ROLE = {
 'constants': ('membranes.py','nernst_voltage_V','Thermal voltage RT/F, Nernst potentials, and conversion between current and ionic amount flux.'),
 'acid_base': ('acid_base.py','speciate','Algebraic TA = carbonate alkalinity + finite buffer base + OH - H; pH is solved from TIC, TA and volume.'),
 'bath': ('model.py','_decode_observables','Fixed reservoir composition determines bath pH species, osmolarity, transporter affinities and reversal potentials.'),
 'geometry': ('model.py','_decode_observables','Finite buffer and impermeant pools enter algebraic pH, cell osmolarity, or the neutral charge constraint.'),
 'homeostasis': ('membranes.py','evaluate_homeostasis','NHE1 supplies Na and TA; AE2 imports Cl and exports carbon/TA; CO2 transfers carbon only.'),
 'membranes': ('nbc_minimal.py','evaluate_membrane_closure_with_nbc','Pump and channel currents determine two algebraic voltages and conserved cell/lumen source vectors.'),
 'water': ('water.py','evaluate_water_fluxes','Q_b=L_b(O_i-O_b), Q_a=L_a(O_l-O_i), Q_p=L_p(O_l-O_b), Q_out=k_out max(V_l-V_dead,0).'),
 'initial': ('model.py','initial_state','Constructor reference concentrations and volumes define a seed and charge consistency check, not the frozen production initial condition.'),
 'ae4': ('transporters.py','evaluate_ae4_qss','Shared carrier QSS determines scalar J4; equal routing imposes (-J4/2,-J4/2,J4,-2J4,-2J4) in cell Na,K,Cl,TIC,TA.'),
 'nbc': ('nbc_minimal.py','evaluate_minimal_nbc','JB=G_B u(Ca) tanh(affinity/log_width), with electrogenic 1 Na and 2 bicarbonate inward per cycle.'),
 'nkcc': ('nkcc1_palk2010.py','palk_cycle_flux_fmol_s','JN=alpha*m_N*(A1-A2[Na]i[K]i[Cl]i^2)/(A3+A4[Na]i[K]i[Cl]i^2); inward 1 Na, 1 K, 2 Cl per cycle.'),
 'cha': ('nhe1_cha2009.py','cha_nhe1_rate','JH=N_H*m_H*Mod2*(a*b-c*d)/(a+b+c+d), with rate conversion from ms to s and printed fourth reverse rate retained.'),
 'regulation': ('camp_pka.py','capacity_multiplier','da/dt=(beta-a)/tau; AE4 multiplier=basal+increment*coupling*construct_scale*a.'),
 'nkcc_activation': ('nkcc_stimulation.py','normalized_calcium_arm','m_N=1+(m_max-1)*clip((Ca-Ca_rest)/(Ca_sat-Ca_rest),0,1).'),
 'protocol': ('validation.py','__call__','Declared input arm and time window select calcium and beta activation steps; no agonist dose response map is present.'),
 'Task41': ('ae4_cacc_recruitment.py','recruitment_factor','In the prescribed inverse family, residual stimulated CaCC recruitment is b and the AE4 dependent fraction is rho=1-b.'),
 'genotype': ('model.py','Genotype','Externally passed expression multipliers scale the corresponding transport flux; zero AE4 expression is an exact deletion of its source.'),
 'nbc_design': ('nbc_minimal.py','DERIVED_REFERENCE_NBC_CAPACITY_FMOL_S','Upstream construction: J4_req=((1-s)/s)*L_N; JB_req=J4_req-H_stim/2; capacity includes the self consistent voltage.'),
 'frozen_initial': ('model.py','evaluate','Full conserved state and regulatory initial value at t=0 for the finite duration protocol; all amounts and both volumes evolve subsequently.'),
}

EXPERIMENTS = {
 'nbc.capacity_fmol_s': 'Measure Na dependent bicarbonate loading or alkalinity recovery during stimulation with calibrated pH, TIC and volume, separating NHE1 and AE2 contributions; establish molecular identity independently.',
 'homeostasis.nhe1_cha_carrier_amount_fmol': 'Measure resting and stimulated Na/H exchange turnover during controlled acid loads, with pH and cell volume; distinguish carrier abundance from kinetic activity.',
 'nbc.nhe1_stimulated_multiplier': 'Measure the ratio of stimulated to resting NHE1 flux in the same salivary preparation and calcium protocol.',
 'nkcc.alpha_eff_fmol_s': 'Measure acute intact cell NKCC chloride uptake during AE4 loss with simultaneous intracellular Na, K and Cl; compare separately with the isolated 2015 uptake assay.',
 'membranes.nak_capacity_fmol_s': 'Measure pump current or turnover versus intracellular Na during the secretion protocol and its apical distribution.',
 'membranes.apical_pump_fraction': 'Quantify apical versus basolateral functional pump turnover, rather than inferring the fraction from localisation alone.',
 'membranes.g_cl_apical_S': 'Measure CaCC current versus calcium and voltage in matched WT and AE4 deficient cells, including onset and sustained recruitment.',
 'geometry.cell_buffer_total_fmol': 'Measure intracellular titration capacity together with TIC and cell volume across the experimental pH interval.',
 'geometry.cell_impermeant_osmoles_fmol': 'Measure resting and stimulated cell volume and osmolarity to separate fixed buffer sites from other impermeant osmotic particles.',
 'ae4.carrier_amount_fmol': 'Measure absolute AE4 abundance and exchange turnover in the same preparation; expression alone does not fix per carrier turnover.',
 'regulation.tau_ae4_s': 'Acquire acute beta input, cAMP/PKA activity and AE4 transport time series to distinguish regulation from slow ionic and volume modes.',
 'regulation.fully_activated_increment': 'Measure the AE4 transport increment under matched stimulation and substrate gradients; resolve the discrepant approximate 2021 gain readings.',
 'Task41.b': 'Measure the fraction of stimulated CaCC recruitment lost with acute AE4 depletion while controlling calcium, pH, chloride gradient, voltage and expression; test the prescribed family rather than use the secretion target again.',
}

def source_location(root, filename, symbol):
    path=root/'src/modern_full_model'/filename
    lines=path.read_text().splitlines()
    for i,line in enumerate(lines,1):
        if ('def '+symbol+'(' in line or 'class '+symbol+':' in line or line.lstrip().startswith(symbol+' =')):
            return f'{path.relative_to(root)}:{i}'
    # No guessed line: a file level source remains exact when a symbol is absent.
    return str(path.relative_to(root))


def enrich_inventory(rows,root):
    crosswalk=root/'analysis/43_parameter_provenance/output/task44_sensitivity_crosswalk.json'
    crosswalk_data=json.loads(crosswalk.read_text()) if crosswalk.exists() else {}
    sensitivity={r['inventory_parameter']:r for r in crosswalk_data.get('parameters',[])}
    selected=set(x[0] for x in ROLE.values())|{'ae4_equal_cation_routing.py','task31_nhe1_repair.py','parameters.py'}
    source_lines={str((root/'src/modern_full_model'/f).relative_to(root)):(root/'src/modern_full_model'/f).read_text().splitlines() for f in selected}
    for r in rows:
        name=r['name'];group=name.split('.')[0];leaf=name.split('.')[-1]
        filename,symbol,role=ROLE[group]
        r['equation_role']=role
        r['equation_source']=source_location(root,filename,symbol)
        tokens={leaf}
        if group=='cha' and leaf.endswith('_modifier_hill'):tokens={leaf.upper()}
        if group=='nkcc' and leaf in {'A2','A4'}:tokens={leaf+'_MM'}
        if name=='regulation.tau_ae4_s':tokens={'tau_activation_s'}
        locations=[]
        for path,lines in source_lines.items():
            for i,line in enumerate(lines,1):
                if any(token in line for token in tokens):locations.append(f'{path}:{i}')
        r['code_locations']=sorted(set(locations))
        r['selected_law_usage']='Used by selected production law or its initialisation/closure.' if r['active'] else 'Not a dynamic parameter of the parent production law.'
        if name in {'protocol.cch_dose_uM','protocol.ipr_dose_uM'}:
            r['active']=False
            r['selected_law_usage']='Protocol dose metadata only. __call__ returns prescribed calcium and beta values independently of these dose fields.'
            r['justification']='Nominal experiment dose labels retained for provenance. The selected model has no dose response map from these fields to calcium or beta activation.'
        if name=='protocol.beta_occupancy_on':r['units']='dimensionless fraction'
        if name=='protocol.ae4_expression_cases':
            r['equation_role']='Task40 calls the same model with e4=1,0.05,0; eNKCC=eNHE=eAE2=1 and all other parameters fixed.'
            r['equation_source']='analysis/40_ae4_equal_cation_routing/run_case.py:78'
            r['code_locations']=['analysis/40_ae4_equal_cation_routing/run_case.py:21','analysis/40_ae4_equal_cation_routing/run_case.py:78']
        if name.startswith('frozen_initial.'):
            r['selected_law_usage']='Used as an initial condition. This is an output of the frozen WT root solve, not a kinetic coefficient.'
            r['code_locations']=['results/40_ae4_equal_cation_routing/wt_rest.json:state_vector']
        if name.startswith('nbc_design.'):
            r['selected_law_usage']='Used upstream in constructing the frozen NBC capacity; not independently read by the production RHS.'
        if name.startswith('initial.'):
            r['selected_law_usage']='Constructor reference and consistency checks only; not the Task 40 initial state.'
        if name=='geometry.fixed_cell_anion_equivalents_fmol':
            r['selected_law_usage']='Reference charge consistency and neutral charge manifold. The charge amount is fixed; it is not an additional transport source.'
        if name=='nbc.nhe1_stimulated_multiplier':
            r['equation_role']='JH_stim=[1+(multiplier-1)*u(Ca)]*JH_Cha; shares the NBC recruitment endpoints.'
            r['equation_source']=source_location(root,'nbc_minimal.py','nhe1_stimulation_multiplier')
        if name.startswith('homeostasis.co2_'):
            r['equation_role']='J_CO2,b=P_b([CO2]bath-[CO2]cell); J_CO2,a=P_a([CO2]lumen-[CO2]cell); modifies carbon amounts, not alkalinity.'
            r['equation_source']=source_location(root,'nbc_minimal.py','evaluate')
        r['literature_source']='No independently verified primary measurement attached to this numerical setting; model lineage is the source.'
        if group=='cha':r['literature_source']='Cha et al. 2009, Eqs. 3 and 5, Table S1; DOI 10.1016/j.bpj.2009.08.053; repository results/31_nhe1_mechanistic_repair/source_table_s1.json.'
        elif group=='nkcc':r['literature_source']='Palk et al. 2010/Benjamin and Johnson 1997 law, reproduced as Vera Siguenza et al. 2018 Eq. 17; coefficients in src/modern_full_model/nkcc1_palk2010.py.'
        elif group=='water' and 'hydraulic' in name:r['literature_source']='Palk 2010/2018 inherited whole cell water coefficients; dimensional conversion is recorded in parameters.py.'
        elif group=='regulation':r['literature_source']='Pena Munzenmayer et al. 2021 gives pathway and endpoint activation evidence; acute kinetic times and exact common capacity gain remain effective settings.'
        elif name in {'protocol.cch_dose_uM','protocol.ipr_dose_uM','protocol.duration_s'}:r['literature_source']='2015 protocol lineage recorded by WTProtocol.provenance in validation.py; implementation transfer retains dose labels and duration.'
        elif r['classification']=='physical_constant':r['literature_source']='Standard physical constant, rounded to the stored precision; not a physiological estimate.'
        r['what_constrains_it']=r['justification']
        r['what_does_not_justify_it']='An accepted physiological point does not identify this setting or supply a numerical uncertainty interval. No knockout fit or new measurement is performed in this audit.'
        if r['classification']=='published_kinetic_coefficient':r['what_does_not_justify_it']='Publication of the kinetic law does not measure salivary carrier abundance or establish an uncertainty interval for transfer to this preparation.'
        if name.startswith('nbc_design.') or name=='nbc.capacity_fmol_s':r['what_does_not_justify_it']='The 70/30 partition is assumed, not measured. Solving the coupled current balance does not establish the physiological truth or unique identification of that partition.'
        if name=='Task41.b':r['what_does_not_justify_it']='Matching the target secretion comparison is inverse selection, not independent validation or experimental estimation of CaCC recruitment.'
        cls=r['classification']
        r['uncertainty_status']='genuinely_unknown_for_physiological_transfer'
        if cls=='physical_constant':r['uncertainty_status']='fixed_standard_constant_at_stored_precision'
        elif cls in {'inactive_legacy','constructor_seed','law_selection','numerical_bracket','numerical_regularisation'}:r['uncertainty_status']='not_a_measured_uncertainty_parameter'
        elif group=='frozen_initial':r['uncertainty_status']='conditional_root_value_no_experimental_interval'
        elif group=='nbc_design':r['uncertainty_status']='conditional_design_value_no_experimental_interval'
        elif cls in {'WT_constraint_derived','reference_constraint_derived','WT_architecture_derived'}:r['uncertainty_status']='conditional_derivation_no_experimental_interval'
        r['admissible_interval']=None
        r['interval_source']=None
        value=r['value'];r['logical_domain']='Selected discrete setting or stated numerical reference.'
        if isinstance(value,(int,float)) and not isinstance(value,bool):
            r['logical_domain']='Finite stored value; any physical or constructor sign constraint is necessary, not sufficient for physiological admissibility.'
            if value>0:r['logical_domain']='Positive for the retained law; no finite upper physiological bound established.'
            elif value==0:r['logical_domain']='Exact zero records a declared omitted contribution or protocol normalisation; it is not a confidence bound.'
        if ('fraction' in leaf or name in {'protocol.beta_occupancy_on','Task41.b'} or leaf.endswith('SHARE')):
            r['logical_domain']='[0,1] by fractional definition only; not a measured uncertainty interval.'
        if name.startswith('acid_base.') and 'ph_' not in name:r['logical_domain']='Finite p scale dissociation constant; environmental applicability must be measured.'
        if name in {'acid_base.ph_lower','acid_base.ph_upper'}:r['logical_domain']='Ordered numerical inversion endpoints [3,11]; not the inherited physiological pH gate [6.6,7.3].'
        if name=='nbc.nhe1_stimulated_multiplier':r['logical_domain']='At least 1 by the selected model constructor; this is a structural choice, not an experimental lower bound.'
        if name=='Task41.b':r['logical_domain']='0 < b <= 1, equivalently 0 <= rho < 1, by the prescribed family constructor. b=0 is excluded; this is not a measured uncertainty interval.'
        if group=='genotype' and leaf!='name':r['logical_domain']='Nonnegative expression scale; values 1, 0.05 and 0 are specified interventions, not a fitted interval.'
        r['manuscript_dependence']='Numerical state, current and secretion claims are conditional on this active setting; exact stoichiometric identities are independent of its numerical value.' if r['active'] else 'No independent parent RHS sensitivity. Used only for the explicitly stated metadata, initial reference, design lineage or inverse family role.'
        if group=='regulation':r['manuscript_dependence']='Regulatory transient claims depend on the selected family/timescale and gain; constant full stimulus equilibrium is independent of tau within R1.'
        if name=='Task41.b':r['manuscript_dependence']='Directly controls the target selected Task 41 CaCC recruitment mechanism and its 600 s inverse crossing; outside the parent Task 40 model.'
        r['constraint_experiment']=EXPERIMENTS.get(name,
            'Measure this quantity or its identifiable combination in the same preparation and stimulus protocol, with independent WT data; determine physiological ranges before making global robustness claims.')
        if cls=='physical_constant':
            r['constraint_experiment']='No physiological experiment required; retain standard value and unit convention.'
            r['what_does_not_justify_it']='No physiological fitting or salivary uncertainty interval applies to a standard physical constant.'
        elif cls in {'numerical_bracket','numerical_regularisation'}:
            r['constraint_experiment']='Verify numerical convergence and successful closure on the stated domain; this is not an abundance measurement.'
        elif cls in {'inactive_legacy','constructor_seed'}:
            r['constraint_experiment']='No new production parameter experiment is implied for this inactive field; measure the actual corresponding state or active pathway instead.'
        r['sensitivity_crosswalk']='See output/task44_sensitivity_crosswalk.json when available; a missing numerical derivative is not evidence of insensitivity.'
        if name in sensitivity:
            z=sensitivity[name]['results']
            r['sensitivity_crosswalk']={'source':'output/task44_sensitivity_crosswalk.json',
                 'source_commit':crosswalk_data['source_commit'],
                 'status':'local_numerical','WT_secretion_log_elasticity':z['WT_Q_log_elasticity'],
                 'WT_pH_per_log_parameter':z['WT_pH_per_log_parameter'],
                 'null_stationary_deficit_log_elasticity':z['null_stationary_deficit_log_elasticity'],
                 'chloride_compensation_gain_log_elasticity':z['chloride_compensation_gain_log_elasticity'],
                 'WT_slowest_decay_log_elasticity':z['WT_slowest_decay_log_elasticity'],
                 'noNBC_capacity_ratio_log_elasticity':z['noNBC_capacity_to_current_demand_log_elasticity'],
                 'finite_change_overturning_a_claim_established':False}
            r['manuscript_dependence']+=' Local numerical dependence is quantified in the linked crosswalk; it establishes neither a measured uncertainty range nor global robustness.'


def write_source_manifest(rows,root,out):
    paths={str(p.relative_to(root)) for p in (root/'src/modern_full_model').glob('*.py')}
    paths.update({'reference/native_wt_contract.json','results/40_ae4_equal_cation_routing/frozen_inputs.json',
                  'results/40_ae4_equal_cation_routing/wt_rest.json',
                  'analysis/40_ae4_equal_cation_routing/validation_common.py',
                  'analysis/39_palk_nkcc1_full_validation/validation_common.py'})
    manifest=[{'path':p,'sha256':hashlib.sha256((root/p).read_bytes()).hexdigest()} for p in sorted(paths)]
    (out/'source_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
