"""Rebuild and verify provenance outputs without fitting or new dynamics."""
from __future__ import annotations
from pathlib import Path
from dataclasses import asdict,replace
import csv,hashlib,importlib.util,json,subprocess,sys
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OUT=HERE/'output'
tracked_outputs=['parameter_inventory.json','parameter_inventory.csv','audit_summary.json','source_manifest.json','numerical_controls.json']
before={n:(OUT/n).read_bytes() for n in tracked_outputs}
for r in json.loads(before['source_manifest.json']):
    assert hashlib.sha256((ROOT/r['path']).read_bytes()).hexdigest()==r['sha256'],r['path']
subprocess.run([sys.executable,str(HERE/'audit.py')],check=True,stdout=subprocess.DEVNULL)
assert all((OUT/n).read_bytes()==before[n] for n in tracked_outputs),'Saved audit outputs differ from clean recomputation.'
rows=json.loads(before['parameter_inventory.json']);by={r['name']:r for r in rows}
assert len(by)==len(rows)==137
required=['name','value','units','active','classification','source_path','equation_role','equation_source','code_locations',
          'selected_law_usage','literature_source','what_constrains_it','what_does_not_justify_it','uncertainty_status',
          'admissible_interval','logical_domain','manuscript_dependence','constraint_experiment']
for r in rows:
    assert all(k in r for k in required),r['name']
    assert all(r[k] for k in ['equation_role','selected_law_usage','what_constrains_it','what_does_not_justify_it']),r['name']
    assert r['admissible_interval'] is None,'No physiological uncertainty interval is sourced by this audit.'
    for location in [r['equation_source']]+r['code_locations']:
        assert (ROOT/location.split(':')[0]).is_file(),location
assert len(list(csv.DictReader((OUT/'parameter_inventory.csv').open())))==len(rows)
sys.path.insert(0,str(HERE))
# The import executes the same source reconstruction and hash checks used by audit.py.
import audit
model=audit.model
from modern_full_model.parameters import parameter_records
for r in parameter_records(model.parameters):assert by[r.name]['value']==r.value,r.name
for prefix,obj in [('ae4',model.ae4_parameters),('nbc',model.nbc_parameters),('nkcc',model.nkcc1_kinetics),('nkcc_activation',model.regulatory_model.nkcc1_regulatory_model),('protocol',model.stimulus)]:
    for n,v in asdict(obj).items():
        assert prefix+'.'+n in by,prefix+'.'+n
        assert by[prefix+'.'+n]['value']==v,prefix+'.'+n
assert by['cha.intracellular_modifier_hill']['value']==3
assert by['cha.extracellular_modifier_hill']['value']==1
assert by['regulation.construct']['value']=='WT'
from modern_full_model.model import WT
for n,v in asdict(WT).items():assert by['genotype.'+n]['value']==v
assert by['protocol.ae4_expression_cases']['value']==[1.,.05,0.]
assert by['protocol.beta_occupancy_on']['units']=='dimensionless fraction'
assert by['Task41.b']['logical_domain'].startswith('0 < b <= 1')
assert sum(r['classification']=='inactive_legacy' for r in rows)==5
assert sum(r['classification']=='constructor_seed' for r in rows)==12
assert len([r for r in rows if r['name'].startswith('frozen_initial.')])==model.layout.size==13
for n in ['cch_dose_uM','ipr_dose_uM']:
    changed=replace(model.stimulus,**{n:2*getattr(model.stimulus,n)})
    for t in [0,1e-6,30,600,601]:assert changed(t)==model.stimulus(t)
    assert by['protocol.'+n]['active'] is False
assert by['Task41.b']['classification']=='target_selected' and not by['Task41.b']['active']
assert by['nbc_design.TARGET_NKCC1_POSITIVE_CL_SHARE']['classification']=='effective_assumption'
assert by['nbc.capacity_fmol_s']['classification']=='WT_architecture_derived'
# These protected sources must be byte-identical to the recovered Task 43 checkpoint.
protected=subprocess.run(['git','diff','--name-only','693f44c0b05e92faec737aa7793c3a99ec76445c','--','src','manuscript','model','results','archive','reference'],cwd=ROOT,text=True,capture_output=True,check=True).stdout.splitlines()
assert not protected,protected
all_changes=subprocess.run(['git','diff','--name-only','693f44c0b05e92faec737aa7793c3a99ec76445c'],cwd=ROOT,text=True,capture_output=True,check=True).stdout.splitlines()
assert all(p.startswith('analysis/43_parameter_provenance/') or p=='.github/workflows/task43-analysis.yml' for p in all_changes),all_changes
crosswalk=OUT/'task44_sensitivity_crosswalk.json'
if crosswalk.exists():
    cross=json.loads(crosswalk.read_text())
    assert len(cross['parameters'])==16
    assert sum(len(r['inverse_family_sensitivities']) for r in cross['parameters'])==16
    for r in cross['parameters']:
        assert r['inventory_parameter'] in by,r['inventory_parameter']
        assert r['active_value']==by[r['inventory_parameter']]['value']
    assert cross['stationary_verification']['status']=='passed'
    assert cross['stationary_verification']['independent_nearby_WT_null_roots']==128
    assert cross['inverse_verification']['status']=='PASS'
    assert cross['inverse_verification']['sensitivity_rows']==16
    assert cross['source_commit']=='f783785a863440df459e9c5530beba4c7eb59110'
    subprocess.run([sys.executable,str(HERE/'import_crosswalk.py'),'--check','--source-commit',cross['source_commit']],check=True)
    assert cross['no_NBC_equilibrium']['physiological']
    assert cross['no_NBC_equilibrium']['secretion_decrease_fraction']<.031
    for row in cross['parameters']:
        assert by[row['inventory_parameter']]['sensitivity_crosswalk']['source_commit']==cross['source_commit']
tables=['audit_counts.tex','inventory_table.tex','sensitivity_table.tex']
saved_tables={n:(OUT/n).read_bytes() for n in tables}
subprocess.run([sys.executable,str(HERE/'prepare_report.py')],check=True,stdout=subprocess.DEVNULL)
assert all((OUT/n).read_bytes()==saved_tables[n] for n in tables),'Report tables differ from their machine-readable sources.'
checks={
 'passed':True,'records':len(rows),'active_records':sum(r['active'] for r in rows),
 'production_parameter_hashes_verified':True,'frozen_initial_state_hash_verified':True,
 'all_actual_nested_dataclass_fields_covered':True,'cha_fixed_exponents_verified':True,
 'metadata_only_agonist_doses_verified':True,'machine_readable_schema_verified':True,
 'saved_outputs_match_recomputation':True,'all_record_source_paths_resolve':True,
 'production_source_manifest_entries':len(json.loads(before['source_manifest.json'])),
 'protected_source_changes':protected,'all_branch_changes_within_task43':True,'sensitivity_crosswalk_present':crosswalk.exists(),
 'sensitivity_source_commit':cross['source_commit'],
 'sensitivity_source_file_hashes_reverified':len(cross['imported_source_files']),
 'sensitivity_crosswalk_matches_committed_source':True,'report_tables_match_machine_readable_sources':True,
 'parameter_fits':0,'new_sensitivity_simulations_on_task43':0,
 'scope':'Inventory reconstruction, selected input law checks, original parameter and state hashes, source files, uncertainty semantics and CSV/JSON consistency.'}
(OUT/'verification.json').write_text(json.dumps(checks,indent=2,sort_keys=True)+'\n')
print(json.dumps(checks,indent=2))
