"""Orchestrator-only receipt preparation and connected-publication synchronisation."""
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OUT=HERE/'output'
BRANCH='analysis/task-51-beta-nkcc-swelling-vrac-mechanistic-reconstruction'
START='ad0e7344c8ba79552cbf54714b42e7131bcc1142'

def git(*args):
    return subprocess.check_output(['git',*args],cwd=ROOT,text=True).strip()

def dump(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(data,indent=2,sort_keys=True,allow_nan=False)+'\n')

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def verified(name):
    records=json.loads((OUT/'publication_verification.json').read_text())
    row=next(r for r in records if r['checkpoint']==name)
    if not row['verified']:raise RuntimeError('Unverified dependency')
    return row

def integrity():
    inputs=json.loads((HERE/'input/frozen_inputs.json').read_text())['sha256']
    wrong=[p for p,h in inputs.items() if sha(ROOT/p)!=h]
    if wrong:raise RuntimeError(f'Frozen scientific input changed: {wrong}')
    return len(inputs)

def prepare(name,status,next_step):
    if git('branch','--show-current')!=BRANCH:raise RuntimeError('Wrong branch')
    receipt=OUT/'checkpoints'/f'{name}.json'
    if receipt.exists():raise RuntimeError('Receipt already exists')
    parent=git('rev-parse','HEAD')
    records=json.loads((OUT/'publication_verification.json').read_text())
    if name=='51A':assert parent==START
    else:assert records[-1]['sha']==parent and records[-1]['verified']
    count=integrity()
    (HERE/'CURRENT_STATUS.md').write_text(
        f'# Task 51 current status\n\nCheckpoint: **{name}**\n\n{status}\n\n'
        f'Branch: `{BRANCH}`. Requested start: `{START}`.\n\n'
        f'Previous verified remote head: `{parent}`.\n\nNext: {next_step}\n\n'
        'One source-backed beta-NKCC/swelling/VRAC mechanism only. Task 50 is a separate frozen benchmark. '
        'No phenotype calibration, new signalling states or post-reveal retuning. '
        'Only the orchestrator executes scientific calculations and publishes.\n\n'
        'Every dependent stage requires the preceding remote publication receipt. '
        'A failed identification gate is not a completed numerical prediction. Stop after verified 51F; Task 52 remains unstarted.\n')
    paths=[p for p in HERE.rglob('*') if p.is_file() and '__pycache__' not in p.parts and 'checkpoints' not in p.parts]
    paths += [ROOT/'docs/MANDATORY_RESEARCH_LEDGER.md']
    dump(receipt,dict(checkpoint=name,status=status,next_step=next_step,parent_sha=parent,
        requested_start_sha=START,branch=BRANCH,created_utc=datetime.now(timezone.utc).isoformat(),
        verified_dependencies=records,frozen_input_hashes_verified=count,
        artifact_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},
        publication='Independent remote ref and commit/tree/receipt are read after publication. Verification enters the next boundary; final SHA is externally recorded.'))

def bundle(name,destination):
    receipt=json.loads((OUT/'checkpoints'/f'{name}.json').read_text());parent=receipt['parent_sha']
    assert git('rev-parse','HEAD')==parent
    changes=git('status','--porcelain').splitlines()
    allowed=(str(HERE.relative_to(ROOT))+'/', 'docs/MANDATORY_RESEARCH_LEDGER.md')
    if any(not any(a in s for a in allowed) for s in changes):raise RuntimeError('Unexpected change')
    git('add',str(HERE.relative_to(ROOT)),'docs/MANDATORY_RESEARCH_LEDGER.md')
    message=f'Task 51 checkpoint {name}: {receipt["status"][:140]}'
    git('commit','-m',message)
    entries=[]
    for p in git('diff','--name-only',parent,'HEAD').splitlines():
        f=ROOT/p
        if not f.is_file():raise RuntimeError('Deletion forbidden')
        entries.append(dict(path=p,mode='100644',type='blob',content=f.read_text(),expected_sha=git('hash-object',p)))
    dump(Path(destination),dict(checkpoint=name,parent=parent,base_tree=git('rev-parse',parent+'^{tree}'),
        expected_tree=git('rev-parse','HEAD^{tree}'),message=message,entries=entries))
    print(json.dumps(dict(checkpoint=name,files=len(entries),local_tree=git('rev-parse','HEAD^{tree}'))))

def sync(name,metadata):
    m=json.loads(Path(metadata).read_text());assert git('rev-parse','HEAD^{tree}')==m['tree']['sha']
    assert m['parents'][0]['sha']==git('rev-parse','HEAD^')
    lines=['tree '+m['tree']['sha']]+['parent '+x['sha'] for x in m['parents']]
    for kind in ['author','committer']:
        a=m[kind];stamp=int(datetime.fromisoformat(a['date'].replace('Z','+00:00')).timestamp())
        lines.append(f'{kind} {a["name"]} <{a["email"]}> {stamp} +0000')
    raw=('\n'.join(lines)+'\n\n'+m['message']).encode();found=False
    for b in [raw,raw+b'\n']:
        if hashlib.sha1(b'commit '+str(len(b)).encode()+b'\0'+b).hexdigest()==m['sha']:
            subprocess.run(['git','hash-object','-t','commit','-w','--stdin'],cwd=ROOT,input=b,check=True,stdout=subprocess.DEVNULL);found=True;break
    if not found:raise RuntimeError('Remote raw commit reconstruction mismatch')
    git('reset','--soft',m['sha']);git('update-ref','refs/remotes/origin/'+BRANCH,m['sha'])
    record=dict(checkpoint=name,sha=m['sha'],remote_sha=m['sha'],tree=m['tree']['sha'],verified=True,
        verified_utc=datetime.now(timezone.utc).isoformat(),transport='Connected Git Data API; independent ref and commit/tree/receipt verification')
    records=json.loads((OUT/'publication_verification.json').read_text());records.append(record)
    # Final verification is outside the branch, avoiding a seventh post-51F commit.
    dest=ROOT.parent/'publication_51F_verified.json' if name=='51F' else OUT/'publication_verification.json'
    dump(dest,records)
    print(json.dumps(record))

if __name__=='__main__':
    {'prepare':prepare,'bundle':bundle,'sync':sync}[sys.argv[1]](*sys.argv[2:])
