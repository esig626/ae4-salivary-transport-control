"""Local half of connected GitHub publication: bundle, then verify/synchronise.

The orchestrator sends the bundle through GitHub create_tree/create_commit/
update_ref(force=False), reads the remote ref independently and supplies the
verified commit metadata here. It does not obtain or expose credentials.
"""
import base64
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
import sys
from checkpoint import HERE, ROOT, OUT, BRANCH, git, dump

def bundle(name, destination):
    receipt=json.loads((OUT/'checkpoints'/f'{name}.json').read_text())
    parent=receipt['parent_sha']
    if git('rev-parse','HEAD')!=parent:raise RuntimeError('Unexpected local parent')
    changes=git('status','--porcelain').splitlines()
    if any('analysis/49_camp_vrac_secretory_branch_reconstruction/' not in s for s in changes):
        raise RuntimeError('Out-of-task change')
    git('add',str(HERE.relative_to(ROOT)))
    message=f'Task 49 checkpoint {name}: {receipt["status"]}'
    git('commit','-m',message)
    elements=[]
    for path in git('diff','--name-only',parent,'HEAD').splitlines():
        p=ROOT/path
        if not p.exists():raise RuntimeError('Task publication must not delete files')
        e=dict(path=path,mode='100644',type='blob',expected_sha=git('hash-object',path))
        try:e['content']=p.read_bytes().decode('utf-8')
        except UnicodeDecodeError:e['content_b64']=base64.b64encode(p.read_bytes()).decode()
        elements.append(e)
    dump(Path(destination),dict(checkpoint=name,parent=parent,
        base_tree=git('rev-parse',parent+'^{tree}'),expected_tree=git('rev-parse','HEAD^{tree}'),
        message=message,entries=elements))

def sync(name, metadata):
    m=json.loads(Path(metadata).read_text())
    if git('rev-parse','HEAD^{tree}')!=m['tree']['sha']:
        raise RuntimeError('Connected remote tree differs from local tree')
    # GitHub dates omit the original timezone. Reuse the previous API commit's
    # timezone only if the resulting raw object passes its exact Git SHA.
    parent=git('rev-parse','HEAD^')
    old=git('cat-file','-p',parent)
    tz=next(s.split()[-1] for s in old.splitlines() if s.startswith('author '))
    lines=['tree '+m['tree']['sha']]+['parent '+p['sha'] for p in m['parents']]
    for kind in ['author','committer']:
        a=m[kind];t=int(datetime.fromisoformat(a['date'].replace('Z','+00:00')).timestamp())
        lines.append(f"{kind} {a['name']} <{a['email']}> {t} {tz}")
    raw=('\n'.join(lines)+'\n\n'+m['message']).encode()
    found=False
    for b in [raw,raw+b'\n']:
        h=subprocess.check_output(['git','hash-object','-t','commit','--stdin'],cwd=ROOT,input=b).decode().strip()
        if h==m['sha']:
            subprocess.run(['git','hash-object','-t','commit','-w','--stdin'],cwd=ROOT,input=b,check=True,stdout=subprocess.DEVNULL)
            found=True;break
    if not found:
        subprocess.run(['git','fetch','--depth=1','origin',BRANCH],cwd=ROOT,check=True)
        if git('rev-parse',m['sha']+'^{tree}')!=m['tree']['sha']:raise RuntimeError('Object fetch mismatch')
    git('reset','--soft',m['sha'])
    git('update-ref','refs/remotes/origin/'+BRANCH,m['sha'])
    path=OUT/'publication_verification.json'
    records=json.loads(path.read_text())
    if any(r['checkpoint']==name for r in records):raise RuntimeError('Checkpoint already recorded')
    records.append(dict(checkpoint=name,sha=m['sha'],remote_sha=m['sha'],verified=True,
        verified_utc=datetime.now(timezone.utc).isoformat(),
        transport='Connected GitHub API fast-forward; independent GET verified; exact local tree match'))
    dump(path,records)
    print(json.dumps(records[-1]),flush=True)

if __name__=='__main__':
    if sys.argv[1]=='bundle':bundle(*sys.argv[2:])
    elif sys.argv[1]=='sync':sync(*sys.argv[2:])
    else:raise ValueError('bundle or sync required')
