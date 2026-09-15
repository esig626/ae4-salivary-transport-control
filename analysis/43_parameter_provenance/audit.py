"""Inventory actual production values; do not fit or alter model parameters."""
from __future__ import annotations
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import ast,csv,importlib.util,json,sys,hashlib
from pathlib import Path
from dataclasses import asdict,fields,is_dataclass
from collections import Counter
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent/'output';OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT/'src'))
from modern_full_model.parameters import parameter_records,FullModelParameters
from modern_full_model.nhe1_cha2009 import (published_cha_kinetics,
    INTRACELLULAR_MODIFIER_HILL,EXTRACELLULAR_MODIFIER_HILL)
from modern_full_model.nkcc1_palk2010 import A1,A2_MM,A3,A4_MM
from modern_full_model.validation import sha256_object
spec=importlib.util.spec_from_file_location('task40_readonly',ROOT/'analysis/40_ae4_equal_cation_routing/validation_common.py')
task=importlib.util.module_from_spec(spec);spec.loader.exec_module(task)
_,_,_,rest,model,_=task.models_and_reference()
freeze=json.loads((ROOT/'results/40_ae4_equal_cation_routing/frozen_inputs.json').read_text())
assert task.parameter_hashes(model)==freeze['stimulus_parameter_hashes']
rows=[];defaults={r.name:r.value for r in parameter_records(FullModelParameters())}
legacy={'homeostasis.nkcc1_capacity_fmol_s','homeostasis.nhe1_capacity_fmol_s',
        'homeostasis.nhe1_published_g_fmol_s','homeostasis.nhe1_published_k_h_mM','homeostasis.nhe1_published_k_na_mM'}
for rec in parameter_records(model.parameters):
    r=asdict(rec);r['group']='whole_cell';r['active']=True;r['classification']='effective_assumption'
    r['source_path']='src/modern_full_model/parameters.py';r['justification']=rec.note
    r['uncertainty']='No justified numerical uncertainty interval established by this audit.'
    r['default_value']=defaults[rec.name];r['differs_from_default']=defaults[rec.name]!=rec.value
    n=rec.name
    if n in legacy:
        r.update(active=False,classification='inactive_legacy',justification='Not used by selected Palk NKCC1 or Cha NHE1 flux. Retained for older models and payload compatibility.')
    elif n.startswith('initial.'):
        r.update(active=False,classification='constructor_seed',justification='Constructor reference tuple, not the actual production initial state. Task 40 uses the independently frozen WT state; this tuple also defines charge consistency checks.')
    elif n in ['acid_base.ph_lower','acid_base.ph_upper']:
        r.update(classification='numerical_bracket',justification='Numerical pH inversion bracket, not an experimental physiological limit.')
    elif n=='geometry.lumen_buffer_total_fmol':
        r.update(classification='numerical_regularisation',justification='Positive approximation to no fixed luminal buffer. Its size is numerical, not measured.')
    elif n in ['constants.gas_constant_J_mol_K','constants.faraday_C_mol']:
        r.update(classification='physical_constant',uncertainty='Not a fitted physiological parameter.')
    elif n=='homeostasis.nhe1_model':
        r.update(classification='law_selection',justification='Production explicitly selects Cha 2009 eight-state Mod2; the default field note refers to an older selection.')
    elif n=='homeostasis.nhe1_cha_carrier_amount_fmol':
        r.update(classification='WT_constraint_derived',source_path='src/modern_full_model/task31_nhe1_repair.py',
          justification='One salivary carrier amount calibrated against a resting pH target 6.91 on the selected R09 background. It is not the cardiac carrier count or a direct salivary abundance measurement.')
    elif n=='geometry.fixed_cell_anion_equivalents_fmol':
        r.update(classification='reference_constraint_derived',justification='Derived from the reference initial cations, chloride and total alkalinity. Conditional on that reference tuple, not independently measured.')
    elif n in ['bath.cl_mM']:
        r.update(classification='reference_constraint_derived',justification='Bath electroneutrality using the adopted bath sodium, potassium and carbonate speciation.')
    elif n=='geometry.cell_impermeant_osmoles_fmol':
        r.update(classification='inherited_effective_setting',source_path='reference/native_wt_contract.json',justification='Active R09 background contains 124.69936937124379 fmol, not the dataclass default. An effective osmotic closure setting, not a directly measured molecular pool.')
    elif n.startswith('water.') and 'hydraulic' in n:
        r.update(classification='published_model_conversion',justification=rec.note+' Inherited model coefficient, not a new direct measurement.')
    elif rec.provenance in ['PUBLISHED_MODEL','Provenance.PUBLISHED_MODEL']:
        r.update(classification='inherited_model_setting',justification=rec.note+' Source metadata identifies a model lineage; independent experimental measurement is not established here.')
    if rec.value==0 and n not in legacy and not n.startswith('initial.'):
        r['justification']+=' An exact zero is an explicit structural omission, not a small positive fitted value.'
    rows.append(r)

def add(name,value,units,cls,source,why,active=True):
    rows.append(dict(name=name,value=value,units=units,provenance='audited_active_object',note='',group=name.split('.')[0],
      active=active,classification=cls,source_path=source,justification=why,
      uncertainty='No measured uncertainty interval supplied; vary only with stated conditional assumptions.',
      default_value='',differs_from_default=''))
for n,v in asdict(model.ae4_parameters).items():
    add('ae4.'+n,v,'fmol' if 'amount' in n else 's^-1' if 'rate' in n else 'configuration','inherited_effective_setting',
        'reference/native_wt_contract.json; src/modern_full_model/transporters.py',
        'Shared-carrier model setting. Thermodynamic construction fixes forward/reverse ratios, not this absolute abundance or attempt rate.',active=v is not None)
for n,v in asdict(model.nbc_parameters).items():
    cls='WT_architecture_derived' if n=='capacity_fmol_s' else 'effective_assumption'
    why='Inherited NBC recruitment parameter; not an isoform-specific measurement.'
    if n=='capacity_fmol_s':why='Derived from an assumed 70/30 positive chloride loading architecture at the Task 31 reference state, a 2.3 NHE1 multiplier and coupled electrical closure. Derived conditional on those assumptions, not measured or identified by knockout data.'
    if n=='nhe1_stimulated_multiplier':why='Inherited source-fixed 2.3 multiplier. Fixing it before a test prevents retuning but does not make its transfer to this tissue/protocol uncertainty-free.'
    add('nbc.'+n,v,'fmol/s' if n=='capacity_fmol_s' else 'uM' if 'calcium' in n else 'dimensionless',cls,'src/modern_full_model/nbc_minimal.py; analysis/36_minimal_nahco3_alkalinity/design.md',why)
add('nkcc.alpha_eff_fmol_s',model.nkcc1_kinetics.alpha_eff_fmol_s,'fmol/s','WT_constraint_derived','analysis/39_palk_nkcc1_full_validation/validation_common.py',
    'J_N,rest divided by the published shape at the accepted resting concentrations. The resting flux is itself an inherited model value, not a direct new measurement.')
for n,v,u in [('A1',A1,'dimensionless'),('A2',A2_MM,'mM^-4'),('A3',A3,'dimensionless'),('A4',A4_MM,'mM^-4')]:
    add('nkcc.'+n,v,u,'published_kinetic_coefficient','src/modern_full_model/nkcc1_palk2010.py','Palk/Benjamin effective two-state law in the fixed bath, with fourth-order concentration factors converted to mM. Do not vary bath composition without revisiting these fixed-bath coefficients.')
for n,v in asdict(published_cha_kinetics()).items():
    add('cha.'+n,v,'ms^-1' if 'per_ms' in n else 'mM','published_kinetic_coefficient','src/modern_full_model/nhe1_cha2009.py; results/31_nhe1_mechanistic_repair/source_table_s1.json',
      'Printed Cha 2009 Table S1 coefficient, transferred as a kinetic law. The printed fourth reverse rate is retained despite its documented rounding discrepancy from exact microscopic reversibility.')
reg=model.regulatory_model
add('regulation.tau_ae4_s',reg.ae4_regulatory_model.tau_activation_s,'s','effective_assumption','src/modern_full_model/camp_pka.py','A 30 s effective response time, not a resolved measurement of the intermediate signalling kinetics.')
for n,v in asdict(reg.ae4_regulatory_model.gain).items():
    add('regulation.'+n,v,'dimensionless','effective_assumption','src/modern_full_model/camp_pka.py','Effective normalisation or gain; pathway activation evidence does not identify each coefficient independently.')
for n,v in asdict(reg.nkcc1_regulatory_model).items():
    add('nkcc_activation.'+n,v,'uM' if 'calcium' in n else 'dimensionless','inherited_effective_setting','src/modern_full_model/nkcc_stimulation.py','Inherited algebraic capacity modulation. Its 0.10 uM saturation endpoint is distinct from the NBC 0.25 uM endpoint; at the production input both are fully recruited.')
for n,v in asdict(model.stimulus).items():
    add('protocol.'+n,v,'s' if n.endswith('_s') else 'uM' if n.endswith('_uM') else 'configuration','protocol_setting','src/modern_full_model/validation.py','Declared experimental-input representation, including a step waveform; this is not evidence for measured intracellular calcium kinetics.')
add('Task41.b',.10511872843288446,'dimensionless','target_selected','src/modern_full_model/task41_selected.py','Residual CaCC recruitment in the deliberately inverse Task 41 construction. Not part of the parent Task 40 parameter set.',active=False)
for n,v in [('intracellular_modifier_hill',INTRACELLULAR_MODIFIER_HILL),
            ('extracellular_modifier_hill',EXTRACELLULAR_MODIFIER_HILL)]:
    add('cha.'+n,v,'dimensionless','published_kinetic_coefficient',
        'src/modern_full_model/nhe1_cha2009.py',
        'Cha 2009 Eq. 5 Mod2 exponent; a fixed published exponent, not an identified salivary abundance.')
add('nkcc.law',model.nkcc1_kinetics.law,'model identifier','law_selection',
    'src/modern_full_model/nkcc1_palk2010.py','Explicit Palk/Benjamin law selected by the Task 39 builder and inherited by Task 40.')
add('ae4.law',model.ae4_evaluator.__name__,'model identifier','law_selection',
    'src/modern_full_model/ae4_equal_cation_routing.py','Equal cation routing fixes source shares at one half each while retaining the inherited total scalar cycle flux.')
add('regulation.family',reg.ae4_regulatory_model.family,'model identifier','law_selection',
    'src/modern_full_model/camp_pka.py','One retained effective AE4 activation state; no measured separation of cAMP and PKA kinetics.')
add('regulation.construct',reg.ae4_regulatory_model.construct.value,'construct identifier','protocol_setting',
    'src/modern_full_model/camp_pka.py','WT construct selects unit response scale in the gain map.')
add('nkcc_activation.family',reg.nkcc1_regulatory_model.family,'model identifier','law_selection',
    'src/modern_full_model/nkcc_stimulation.py','Selected instantaneous calcium recruitment family; no added NKCC dynamic state.')
from modern_full_model import nbc_minimal as nbc
for n in ['TARGET_NKCC1_POSITIVE_CL_SHARE','TASK31_R09_NKCC1_CL_FMOL_S','TASK31_R09_NHE1_FMOL_S',
          'DERIVED_AE4_CL_FMOL_S','DERIVED_STIMULATED_NHE1_FMOL_S','DERIVED_REQUIRED_NBC_CYCLE_FMOL_S']:
    add('nbc_design.'+n,getattr(nbc,n),'dimensionless' if n.endswith('SHARE') else 'fmol/s',
        'effective_assumption' if n.endswith('SHARE') else 'WT_architecture_derived',
        'src/modern_full_model/nbc_minimal.py; analysis/36_minimal_nahco3_alkalinity/design.md',
        'Design input or conditional derived quantity upstream of the frozen NBC capacity; not independently evaluated in the production RHS. The 70/30 share is assumed, and the reference fluxes are modelled.',active=False)
frozen_rest=json.loads((ROOT/'results/40_ae4_equal_cation_routing/wt_rest.json').read_text())
assert sha256_object(frozen_rest['state_vector'])==freeze['frozen_rest_state_sha256']
for n,v in zip(model.state_names,frozen_rest['state_vector']):
    add('frozen_initial.'+n,v,'pL' if 'volume' in n else 'dimensionless' if 'fraction' in n else 'fmol',
        'WT_constraint_derived','results/40_ae4_equal_cation_routing/wt_rest.json',
        'Solved, accepted WT resting state used as the shared initial condition for production genotype trajectories; not the constructor seed and not an independent measurement.')
from modern_full_model.model import WT
for n,v in asdict(WT).items():
    add('genotype.'+n,v,'genotype label' if n=='name' else 'dimensionless expression scale',
        'protocol_setting','src/modern_full_model/model.py:Genotype',
        'Externally passed WT expression scale, equal to one for all retained transporters; not a measured abundance. Task 40 alters only the AE4 expression scale.',active=n!='name')
case_tree=ast.parse((ROOT/'analysis/40_ae4_equal_cation_routing/run_case.py').read_text())
cases=next(ast.literal_eval(n.value) for n in case_tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='CASES' for t in n.targets))
add('protocol.ae4_expression_cases',list(cases.values()),'dimensionless expression scales','protocol_setting',
    'analysis/40_ae4_equal_cation_routing/run_case.py:21',
    'Predeclared WT, 5 percent and null interventions; only the AE4 expression field differs. These are protocol conditions, not estimated abundance uncertainty.')
from provenance_metadata import enrich_inventory,write_source_manifest
enrich_inventory(rows,ROOT)
write_source_manifest(rows,ROOT,OUT)
from modern_full_model.validation import PRODUCTION_RADAU, CONSERVATION_RESIDUAL_TOLERANCES
from modern_full_model.camp_pka import INTEGRATION_TRIAL_TOLERANCE
(OUT/'numerical_controls.json').write_text(json.dumps({
    'classification':'numerical controls and unit conventions, not physiological uncertainty',
    'production_solver':asdict(PRODUCTION_RADAU),
    'conservation_tolerances':dict(CONSERVATION_RESIDUAL_TOLERANCES),
    'regulatory_trial_tolerance':INTEGRATION_TRIAL_TOLERANCE,
    'ph_inversion':{'lower':model.parameters.acid_base.ph_lower,'upper':model.parameters.acid_base.ph_upper,
                    'xtol':1e-12,'rtol':1e-13,'source':'src/modern_full_model/acid_base.py:speciate'},
    'voltage_inversion':{'lower_V':-.5,'upper_V':.5,'xtol_V':1e-14,'rtol':1e-13,
                         'source':'src/modern_full_model/nbc_minimal.py:evaluate_membrane_closure_with_nbc'},
    'unit_conventions':{'mM_times_pL_equals_fmol':True,'fmol_per_mol':1e15,'ms_per_s':1000},
    'stimulus_discontinuity':{'initial_observable_time_s':0.0,'integration_start_s':1e-6,
                             'shared_initial_state_preserved':True,
                             'source':'analysis/40_ae4_equal_cation_routing/run_case.py:152'},
    'note':'Tolerances are audited from retained source or production specifications. They do not bound model error or parameter uncertainty.'
},indent=2,sort_keys=True)+'\n')
with (OUT/'parameter_inventory.csv').open('w',newline='') as f:
    keys=list(dict.fromkeys(k for r in rows for k in r));w=csv.DictWriter(f,keys,lineterminator='\n');w.writeheader()
    w.writerows({k:json.dumps(v,sort_keys=True) if isinstance(v,(list,dict)) else v for k,v in r.items()} for r in rows)
(OUT/'parameter_inventory.json').write_text(json.dumps(rows,indent=2,default=str,allow_nan=False)+'\n')
summary=dict(base_main='c4d4f207f702d8eee09ab94d2d90ae90552dd641',records=len(rows),active_records=sum(r['active'] for r in rows),
  classifications=dict(Counter(r['classification'] for r in rows)),active_parameters_hashes=task.parameter_hashes(model),
  calibrated_uncertainty_intervals=0,parameter_fits_performed=0,model_or_frozen_outputs_modified=False,
  frozen_initial_state_sha256=freeze['frozen_rest_state_sha256'],
  coverage='All fields of the actual nested production objects, fixed Cha exponents, selected laws, NBC design assumptions, and frozen initial state. Arithmetic unit conversions and numerical solver tolerances are documented separately, not counted as physiological parameters.',
  corrections=['Cha Mod2 exponents added','Agonist doses are protocol metadata, not inputs used by the selected RHS',
               'Law identifiers, WT construct, NBC construction quantities and actual frozen initial state added'],
  uncertainty_interpretation='Logical domains and physiology gates are not measured uncertainty intervals; all effective parameter intervals remain unestablished.')
(OUT/'audit_summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
print(json.dumps(summary,indent=2))
