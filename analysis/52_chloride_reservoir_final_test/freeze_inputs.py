"""File/import-only immutable52C freeze; no scientific model construction."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
import ast
import hashlib
import json
import platform
import subprocess
import numpy
import scipy
from task52_model import dependency_audit,COUNTS
from diagnostics import LIMITS

TASK=Path(__file__).resolve().parent
ROOT=TASK.parents[1]
OUT=TASK/'output'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return json.loads(p.read_text())
def main():
    target=OUT/'parameter_and_case_freeze.json'
    assert not target.exists(),'Freeze is immutable'
    pub=read(OUT/'publication_52B.json')
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    assert pub['pass'] and pub['commit_sha']==head
    verification=read(OUT/'verification_attempt_01.json')
    assert verification['status']=='PASS'
    for path,expected in verification['source_sha256'].items(): assert sha(TASK/path)==expected
    files=list((TASK/'frozen_task37').rglob('*.py'))+list((TASK/'frozen_task37').rglob('*.json'))
    files += list(TASK.glob('*.py'))
    files += [OUT/n for n in ('projected_onsets.json','case_matrix_52A.json','task37_provenance_audit.json','verification_attempt_01.json','binding_inputs_at_start.json')]
    files += [TASK/n for n in ('ARCHITECTURE_AND_NUMERICS_FREEZE.md','PROJECTION_DERIVATION.md','TASK37_PROVENANCE.md')]
    original=read(OUT/'binding_inputs_at_start.json')
    for path,expected in original['sha256'].items():
        if path=='docs/MANDATORY_RESEARCH_LEDGER.md': continue
        assert sha(ROOT/path)==expected
        files.append(ROOT/path)
    # Parse every candidate source. Runtime scientific imports must all resolve
    # within the hash-verified historical closure; neither operation evaluates it.
    for p in files:
        if p.suffix=='.py': ast.parse(p.read_text(),filename=str(p))
    audit=dependency_audit();assert audit['pass'] and all(v==0 for v in COUNTS.values())
    audit['scientific_evaluations']=dict(COUNTS)
    (OUT/'dependency_audit_52C.json').write_text(json.dumps(audit,indent=2,sort_keys=True)+'\n')
    files.append(OUT/'dependency_audit_52C.json')
    ledger=ROOT/'docs/MANDATORY_RESEARCH_LEDGER.md'
    data={'checkpoint':'52C','parent_verified_commit':head,'production_trajectories_before_freeze':0,
      'retuning_permitted':False,'case_matrix':read(OUT/'case_matrix_52A.json'),
      'projected_onsets':read(OUT/'projected_onsets.json'),
      'immutable_artifact_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(set(files))},
      'append_only_ledger_prefix':{'path':str(ledger.relative_to(ROOT)),'length_bytes':ledger.stat().st_size,'sha256':sha(ledger)},
      'physical_limits':LIMITS,'positive_core_coordinates':list(range(12)),
      'extra_diagnostic_tolerances':{'independent_cell_assembly_fmol_s':1e-10,'tracked_cell_osmoles_fmol_s':1e-10,'aux_current_split_A':1e-20},
      'projection_tolerances':{'charge_fmol':1e-9,'bath_isotonicity_mOsm':1e-9,'pH':1e-10},
      'auxiliary_source_interpretation':'Shared beta demand: 0S native control;2.32e-9S source-central equivalent;4.49e-9S predeclared upper equivalent. No genotype/channel coupling or target fit.',
      'chloride_flux_convention':'J outward apical total Cl;N inward NKCC cycles;A signed AE4 cell Cl source;E inward AE2 cycles;delta is WT minus KO cell Cl amounts.',
      'driving_force_convention':'E_Cl minus V_a; positive means outward chloride flux. I_aux=beta*G_aux*(V_a-E_Cl).',
      'diagnostic_record':'Every accepted endpoint and integer second: paired13states and RHS, cell/lumen concentrations/speciation, flow, voltages, Cl drive/current components, all transport/cycle fluxes, independent conserved-source assembly, delta(t), cumulative8-point flux quadrature.',
      'failed_case_policy':'First sampled physical/conservation failure stops case. Preserve partial trajectory/budget at recorded T;600s endpoints null. Only explicitRadau status failure permits one identical-inputBDF attempt. No reruns, rescue or added cases.',
      'software':{'python':platform.python_version(),'numpy':numpy.__version__,'scipy':scipy.__version__},
      'dependency_audit':audit}
    target.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'freeze_sha256':sha(target),'immutable_artifacts':len(data['immutable_artifact_sha256']),'scientific_evaluations':dict(COUNTS)}))
if __name__=='__main__': main()
