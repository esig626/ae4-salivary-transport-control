"""Disclosed52D recorder-only recovery. Frozen scientific files remain intact."""
import run_frozen_cases as runner
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
from threadpoolctl import threadpool_limits,threadpool_info
import numpy as np

original_diagnose=runner.diagnose
def fixed_schema_diagnose(*args,**kwargs):
    row,faults,ratios=original_diagnose(*args,**kwargs)
    # Rest omits these current aliases. The explicit NBC current is zero there.
    # Active inherited aliases retain their original values exactly.
    row.setdefault('current_nbc_basolateral_A',row['I_nbc_A'])
    row.setdefault('current_basolateral_old_A',row['current_basolateral_total_A']-row['I_nbc_A'])
    return dict(sorted(row.items())),faults,ratios

def verified_software_boundary():
    out=runner.OUT;root=runner.ROOT
    pub=json.loads((out/'publication_52D_recording_boundary.json').read_text())
    assert pub['pass'] and pub['commit_sha']==runner.git('rev-parse','HEAD')
    receipt=json.loads((out/'software_compatibility/receipt.json').read_text())
    for path,expected in receipt['artifacts'].items(): assert runner.sha(root/path)==expected
    freeze=json.loads((out/'parameter_and_case_freeze.json').read_text())
    for path,expected in freeze['immutable_artifact_sha256'].items(): assert runner.sha(root/path)==expected
    assert runner.git('branch','--show-current')=='analysis/task-52-chloride-reservoir-final-test'
    prior=json.loads((out/'publication_52C.json').read_text())
    assert prior['pass'] and prior['commit_sha']==receipt['parent_commit']
    subprocess.run(['git','merge-base','--is-ancestor',prior['commit_sha'],'HEAD'],cwd=root,check=True)
    prefix=freeze['append_only_ledger_prefix'];data=(root/prefix['path']).read_bytes()
    assert hashlib.sha256(data[:prefix['length_bytes']]).hexdigest()==prefix['sha256']
    assert runner.dependency_audit()['pass']
    return freeze

def main():
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['52D','52E']);args=parser.parse_args()
    runner.diagnose=fixed_schema_diagnose
    if args.stage=='52E':
        # The original CLI verifies the independently published52D and all
        # immutable52C inputs. Its case/attempt guards remain in force.
        runner.main();return
    freeze=verified_software_boundary()
    out=runner.OUT;path=out/'execution_52D_recovery.json'
    assert not path.exists(),'Only one explicitly disclosed software recovery is permitted'
    original=json.loads((out/'execution_52D.json').read_text())
    assert not original['complete'] and original['cases']==[]
    assert not list((out/'cases/case_01').glob('*summary.json'))
    record={'stage':'52D','parent_commit':runner.git('rev-parse','HEAD'),
      'freeze_sha256':runner.sha(out/'parameter_and_case_freeze.json'),
      'compatibility_source_sha256':runner.sha(Path(__file__)),
      'original_failed_attempt':'output/software_compatibility/failed_attempt.json',
      'original_attempt_evaluation_counts':'unavailable after process exit; not zero',
      'replayed_case_ids':['case_01'],'distinct_scientific_case_ids':['case_01','case_02','case_03'],
      'cases':[],'complete':False}
    runner.dump(path,record)
    with threadpool_limits(limits=1):
        record['threadpools']=threadpool_info();assert all(p['num_threads']==1 for p in record['threadpools'])
        for case in freeze['case_matrix']['cases'][:3]:
            directory=out/'cases'/case['id']
            if case['id']=='case_01': directory=directory/'recording_recovery'
            directory.mkdir(parents=True,exist_ok=False)
            runner.dump(directory/'started.json',{'case':case,'freeze_sha256':record['freeze_sha256'],
              'recording_adapter_sha256':record['compatibility_source_sha256']})
            states=freeze['projected_onsets'][case['projection']]
            y0=np.r_[states['WT']['state'],states['KO']['state']]
            model=runner.PairedModel(protocol=case['protocol'],G_aux_S=case['G_aux_S'],ko_supply_multiplier=case['ko_supply_multiplier'])
            attempts=['Radau']
            result=runner.run_attempt(case,model,y0,'Radau',directory,freeze['case_matrix']['solver'])
            if result['status']=='solver_status_failure':
                result=runner.run_attempt(case,model,y0,'BDF',directory,freeze['case_matrix']['solver']);attempts.append('BDF')
            summary=directory/f'{attempts[-1]}_summary.json'
            entry={'id':case['id'],'status':result['status'],'attempts':attempts,'selected_attempt':attempts[-1],
              'summary_path':str(summary.relative_to(runner.ROOT)),'summary_sha256':runner.sha(summary)}
            record['cases'].append(entry);runner.dump(path,record)
            print(json.dumps({**entry,'T_s':result['recorded_end_s'],'mass_budget_pass':result['reservoir_budget']['pass']}),flush=True)
        record['complete']=True;record['evaluation_counts']=dict(runner.COUNTS)
        record['dependency_audit']=runner.dependency_audit();assert record['dependency_audit']['pass']
        runner.dump(path,record)
if __name__=='__main__': main()
