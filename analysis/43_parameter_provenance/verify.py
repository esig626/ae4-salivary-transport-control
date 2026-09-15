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
assert len(by)==len(rows)==131
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
crosswalk=OUT/'task44_sensitivity_crosswalk.json'
if crosswalk.exists():
    cross=json.loads(crosswalk.read_text())
    assert len(cross['parameters'])>0
    for r in cross['parameters']:assert r['inventory_parameter'] in by,r['inventory_parameter']
checks={
 'passed':True,'records':len(rows),'active_records':sum(r['active'] for r in rows),
 'production_parameter_hashes_verified':True,'frozen_initial_state_hash_verified':True,
 'all_actual_nested_dataclass_fields_covered':True,'cha_fixed_exponents_verified':True,
 'metadata_only_agonist_doses_verified':True,'machine_readable_schema_verified':True,
 'saved_outputs_match_recomputation':True,'all_record_source_paths_resolve':True,
 'production_source_manifest_entries':len(json.loads(before['source_manifest.json'])),
 'protected_source_changes':protected,'sensitivity_crosswalk_present':crosswalk.exists(),
 'parameter_fits':0,'new_sensitivity_simulations_on_task43':0,
 'scope':'Inventory reconstruction, selected input law checks, original parameter and state hashes, source files, uncertainty semantics and CSV/JSON consistency.'}
(OUT/'verification.json').write_text(json.dumps(checks,indent=2,sort_keys=True)+'\n')
print(json.dumps(checks,indent=2))
