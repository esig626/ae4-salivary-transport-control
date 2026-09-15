"""Recompute the published Task 44 numerical checkpoint from an empty output area.

Run only in the dedicated clean checkout supplied by the root agent. The
published outputs are moved to an immutable scratch comparison snapshot first;
all numerical products are then regenerated, with no inverse trial cache reads.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from time import perf_counter
from threading import Lock

PARAMETERS=('NBC_capacity','AE4_carrier','NKCC_scale','NHE_carrier',
            'pump_capacity','CaCC_conductance','buffer_pool','water_apical')
REL='analysis/44_full_system_mathematics'
# Explicitly metadata only. Parameter/state/source hashes remain compared.
IGNORE_FIELDS={
 ('output/execution_summary.json',('elapsed_s',)):'execution elapsed time only',
 ('output/completion_verification.json',('analysed_commit',)):'precommit provenance versus clean numerical checkpoint revision; ancestry independently verified'
}
HISTORICAL_FILES={'output/recovery_verification.json'}
SELF_FILES={'output/clean_recomputation_verification.json'}
POLICY={
 'default':{'rtol':1e-8,'atol':1e-14},
 'electrical_current_ampere':{'rtol':1e-8,'atol':1e-19},
 'conservation_tolerance_ratio_roundoff':{'rtol':1e-8,'atol':1e-5},
 'small_equilibrium_or_charge_residual':{'rtol':1e-8,'atol':1e-10},
 'rationale':'The same runtime should reproduce deterministic solves nearly exactly. Tiny differences of conservation cancellation residuals are allowed far below their physiological gates. No failed gate, boolean, class, root count or experimental observable is excluded.'}


def dump(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n')


def git(checkout,*args):
    return subprocess.check_output(['git','-C',str(checkout),*args],text=True).strip()


def numeric_files(folder):
    result={}
    for base in ['output','inverse_detail']:
        for p in (folder/base).glob('*'):
            if p.is_file() and p.suffix in ['.json','.csv']:
                relative=str(p.relative_to(folder))
                if relative not in SELF_FILES:result[relative]=p
    return result


def serialise_csv(path):
    with path.open() as source:
        rows=list(csv.DictReader(source))
    def convert(value):
        if value=='':return value
        if value in ('True','False'):return value=='True'
        try:return float(value)
        except ValueError:return value
    return [{k:convert(v) for k,v in row.items()} for row in rows]


def tolerance(path):
    labels=[str(p).lower() for p in path]
    joined='/'.join(labels)
    if 'conservation_ratio' in joined or 'conservation_tolerance_ratio' in joined:
        return 'conservation_tolerance_ratio_roundoff'
    last=next((s for s in reversed(labels) if not s.isdigit()),'')
    if last.endswith('_a') or 'current_residual_a' in joined:
        return 'electrical_current_ampere'
    if last in {'residual_max','full_state_rhs_max','independent_rhs_max','charge_coordinate_recovery_max',
                'charge_max_fmol','max_equilibrium_residual','max_amount_rhs','max_volume_rhs',
                'charge_rate','structural_identity','carbon_alkalinity_identity'}:
        return 'small_equilibrium_or_charge_residual'
    return 'default'


class Comparison:
    def __init__(self):
        self.failures=[];self.skipped=[];self.numeric=0;self.exact=0;self.tolerated=0
        self.largest=[];self.tolerance_uses={};self.file_results=[]

    def fail(self,file,path,reason,expected,actual):
        self.failures.append(dict(file=file,path=list(path),reason=reason,expected=expected,actual=actual))

    def walk(self,file,path,left,right):
        if (file,path) in IGNORE_FIELDS:
            self.skipped.append(dict(file=file,path=list(path),reason=IGNORE_FIELDS[(file,path)]));return
        if isinstance(left,bool) or isinstance(right,bool):
            if type(left)!=type(right) or left!=right:self.fail(file,path,'boolean differs',left,right)
        elif isinstance(left,(int,float)) and isinstance(right,(int,float)):
            self.numeric+=1
            if not (math.isfinite(left) and math.isfinite(right)):
                self.fail(file,path,'nonfinite number',str(left),str(right));return
            if left==right:self.exact+=1;return
            name=tolerance(path);rule=POLICY[name];difference=abs(left-right)
            allowed=rule['atol']+rule['rtol']*max(abs(left),abs(right))
            self.tolerance_uses[name]=self.tolerance_uses.get(name,0)+1
            entry=dict(file=file,path=list(path),expected=left,actual=right,absolute_difference=difference,
                       allowed_difference=allowed,tolerance_class=name,fraction_of_tolerance=difference/allowed)
            self.largest.append(entry)
            self.largest=sorted(self.largest,key=lambda r:r['fraction_of_tolerance'],reverse=True)[:25]
            if difference>allowed:self.fail(file,path,'numeric tolerance exceeded',left,right)
            else:self.tolerated+=1
        elif isinstance(left,dict) and isinstance(right,dict):
            if set(left)!=set(right):self.fail(file,path,'dictionary keys differ',sorted(left),sorted(right))
            for key in sorted(set(left)&set(right)):self.walk(file,path+(key,),left[key],right[key])
        elif isinstance(left,list) and isinstance(right,list):
            if len(left)!=len(right):self.fail(file,path,'list lengths differ',len(left),len(right))
            for n,(a,b) in enumerate(zip(left,right)):self.walk(file,path+(n,),a,b)
        elif type(left)!=type(right) or left!=right:
            self.fail(file,path,'string/null/type differs',left,right)

    def compare(self,expected,actual):
        left=numeric_files(expected);right=numeric_files(actual)
        for name in sorted(set(left)-set(right)):self.fail(name,(),'missing regenerated output',True,False)
        for name in sorted(set(right)-set(left)):self.fail(name,(),'unexpected regenerated output',False,True)
        for name in sorted(set(left)&set(right)):
            failures=len(self.failures);numbers=self.numeric
            loader=serialise_csv if name.endswith('.csv') else lambda p:json.loads(p.read_text())
            self.walk(name,(),loader(left[name]),loader(right[name]))
            self.file_results.append(dict(path=name,numeric_values=self.numeric-numbers,
                status='passed' if len(self.failures)==failures else 'failed',
                historical_source_record_only=name in HISTORICAL_FILES,
                expected_sha256=hashlib.sha256(left[name].read_bytes()).hexdigest(),
                recomputed_sha256=hashlib.sha256(right[name].read_bytes()).hexdigest()))
        return dict(status='passed' if not self.failures else 'failed',
            authoritative_files=len(left),files_compared=len(self.file_results),numeric_values_compared=self.numeric,
            bit_identical_numeric_values=self.exact,numeric_values_with_tolerated_difference=self.tolerated,
            tolerance_policy=POLICY,tolerance_classes_used=self.tolerance_uses,
            exact_metadata_exclusions=self.skipped,
            failures=self.failures,largest_numeric_differences=self.largest,files=self.file_results,
            trial_cache_scope='Inverse trial caches are regenerated from zero and are not authoritative report files. Every top level inverse summary/CSV and every output JSON/CSV is compared. No cache value is read by a solver.')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--checkout',required=True);parser.add_argument('--sha',required=True)
    parser.add_argument('--snapshot',required=True);parser.add_argument('--logs',required=True)
    args=parser.parse_args()
    checkout=Path(args.checkout).resolve();snapshot=Path(args.snapshot).resolve();logs=Path(args.logs).resolve()
    assert git(checkout,'rev-parse','HEAD')==args.sha
    assert not git(checkout,'status','--porcelain'), 'Recomputation must start from a clean checkout.'
    assert not snapshot.exists(), 'Do not overwrite a comparison snapshot.'
    snapshot.mkdir(parents=True);logs.mkdir(parents=True,exist_ok=True)
    analysis=checkout/REL
    source_manifest={str(p.relative_to(checkout)):hashlib.sha256(p.read_bytes()).hexdigest()
                     for folder in [checkout/'src',analysis] for p in sorted(folder.rglob('*.py'))
                     if '__pycache__' not in str(p)}
    protected_tree=git(checkout,'ls-tree','-r','HEAD','src','results','manuscript','archive')
    initial=dict(checkout=str(checkout),source_sha=args.sha,clean_git_status=True,
        source_python_sha256=source_manifest,source_acquisition='Fresh checkout of the remotely published numerical checkpoint; exact HEAD and clean status checked before removing generated outputs.',
        initial_inverse_trial_cache_files=sum(1 for _ in (analysis/'inverse_detail'/'trials').glob('*.json')),
        runtime_python=sys.version)
    dump(logs/'initial_source.json',initial)
    for dirname in ['output','inverse_detail']:
        shutil.move(str(analysis/dirname),str(snapshot/dirname))
        (analysis/dirname).mkdir()
    for path in HISTORICAL_FILES:
        if (snapshot/path).exists():shutil.copyfile(snapshot/path,analysis/path)
    assert not any((analysis/'inverse_detail').rglob('*.json'))
    dump(logs/'progress.json',dict(status='started',completed=[]))
    records=[];progress_lock=Lock()
    env=dict(os.environ)
    for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS']:
        env[key]='1'
    def run(label,args2):
        command=[sys.executable,*args2];start=perf_counter()
        with (logs/(label+'.log')).open('w') as log:
            proc=subprocess.run(command,cwd=checkout,env=env,stdout=log,stderr=subprocess.STDOUT)
        rec=dict(label=label,command=command,returncode=proc.returncode,elapsed_s=perf_counter()-start,
                 log_sha256=hashlib.sha256((logs/(label+'.log')).read_bytes()).hexdigest())
        with progress_lock:
            records.append(rec)
            dump(logs/'progress.json',dict(status='running',completed=records))
        print('CLEAN_STAGE',label,'PASS' if proc.returncode==0 else 'FAIL',round(rec['elapsed_s'],1),'s',flush=True)
        if proc.returncode:raise RuntimeError(label+' failed; inspect '+str(logs/(label+'.log')))
        return rec
    run('baseline',[REL+'/run_analysis.py','--skip-inverse'])
    run('exact',[REL+'/verify_exact.py'])
    run('equilibrium',[REL+'/equilibrium_review.py'])
    run('sensitivity',[REL+'/sensitivity_analysis.py'])
    run('inverse_baseline',[REL+'/inverse_analysis.py','--phase','baseline','--recompute'])
    run('inverse_long',[REL+'/inverse_analysis.py','--phase','long','--recompute'])
    with ThreadPoolExecutor(max_workers=3) as pool:
        pending=[pool.submit(run,'inverse_'+name,[REL+'/inverse_analysis.py','--phase','sensitivity','--parameter',name,'--recompute']) for name in PARAMETERS]
        for done in as_completed(pending):done.result()
    run('inverse_combine',[REL+'/inverse_analysis.py','--phase','combine','--recompute'])
    run('inverse_independent',[REL+'/inverse_verify.py'])
    run('tests',[REL+'/verify.py'])
    run('completion',[REL+'/verify_completion.py'])
    assert git(checkout,'ls-tree','-r','HEAD','src','results','manuscript','archive')==protected_tree
    changed=git(checkout,'diff','--name-only').splitlines()
    assert all(p.startswith(REL+'/') for p in changed),changed
    for name,digest in source_manifest.items():
        assert hashlib.sha256((checkout/name).read_bytes()).hexdigest()==digest,name
    previous_commit=json.loads((snapshot/'output/completion_verification.json').read_text())['analysed_commit']
    ancestor=subprocess.run(['git','-C',str(checkout),'merge-base','--is-ancestor',previous_commit,args.sha]).returncode==0
    assert ancestor, 'Published precommit verification provenance must be an ancestor of the clean checkpoint.'
    assert json.loads((analysis/'output/completion_verification.json').read_text())['analysed_commit']==args.sha
    result=Comparison().compare(snapshot,analysis)
    result['verification_revision_provenance']=dict(published_precommit_check_sha=previous_commit,clean_recomputed_check_sha=args.sha,published_check_is_ancestor_of_clean_checkpoint=ancestor)
    result.update(independent_clean_source=initial,all_commands=records,
        regenerated_inverse_trial_cache_files=sum(1 for _ in (analysis/'inverse_detail'/'trials').glob('*.json')),
        numerical_cache_reuse=False,source_files_unchanged=True,protected_paths_unchanged=True,
        comparison_helper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    dump(logs/'comparison.json',result)
    dump(logs/'progress.json',dict(status=result['status'],completed=records))
    print('CLEAN_COMPARISON',result['status'],'files',result['files_compared'],'numeric',result['numeric_values_compared'],'differences',len(result['failures']),flush=True)
    if result['failures']:sys.exit(1)


if __name__=='__main__':main()
