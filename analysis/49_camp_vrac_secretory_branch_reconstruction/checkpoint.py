"""Orchestrator-only checkpoint preparation, publication and remote verification.

The receipt records its parent and content hashes: a commit cannot contain its
own hash. The independently observed remote SHA is appended after publication
and included in the next commit. Final verification is a separate closeout.
"""
from __future__ import annotations
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BRANCH = 'analysis/task-49-camp-vrac-secretory-branch-reconstruction'
START = '2bf53b773c14f3c49f539957ea3f1964bc8d227f'
OUT = HERE / 'output'

def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT, text=True).strip()

def dump(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False)+'\n')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def verified(checkpoint):
    records = json.loads((OUT/'publication_verification.json').read_text())
    row = next(r for r in records if r['checkpoint'] == checkpoint)
    if not row['verified']:
        raise RuntimeError('Dependent work cannot precede remote verification')
    return row

def prepare(name, status, next_step):
    if git('branch','--show-current') != BRANCH:
        raise RuntimeError('Wrong branch')
    receipt = OUT/'checkpoints'/f'{name}.json'
    if receipt.exists():
        raise RuntimeError('Published checkpoints are immutable; resume instead')
    parent = git('rev-parse','HEAD')
    pub = OUT/'publication_verification.json'
    records = json.loads(pub.read_text()) if pub.exists() else []
    if name != '49A' and (not records or records[-1]['sha'] != parent):
        raise RuntimeError('Parent must be the last verified checkpoint')
    files = {str(p.relative_to(ROOT)):sha(p) for p in HERE.rglob('*')
             if p.is_file() and '__pycache__' not in p.parts
             and p.name != 'CURRENT_STATUS.md' and 'checkpoints' not in p.parts}
    source_files = [*sorted((ROOT/'src/modern_full_model').glob('*.py')),
                    ROOT/'analysis/48_joint_experimental_constraint_reconstruction/protocol_layer.py',
                    ROOT/'analysis/47_dynamic_potassium_recycling_reconstruction/input/active_parameters.json']
    dump(receipt,dict(checkpoint=name, status=status, next_step=next_step,
        created_utc=datetime.now(timezone.utc).isoformat(), orchestrator_only=True,
        branch=BRANCH, requested_start_sha=START, parent_sha=parent,
        verified_dependencies=records, artifact_sha256=files,
        frozen_source_sha256={str(p.relative_to(ROOT)):sha(p) for p in source_files},
        publication='Remote SHA is independently read after commit/push; see next receipt or final closeout.'))
    (HERE/'CURRENT_STATUS.md').write_text(
        f'# Task 49 current status\n\nCheckpoint: **{name}**\n\n{status}\n\n'
        f'Branch: `{BRANCH}`. Requested start: `{START}`.\n\n'
        f'Previous verified remote head: `{parent}`.\n\nNext: {next_step}\n\n'
        'The 49a addendum controls: VRAC only in Stage I, no beta-NKCC change; '
        'fixed full source basis in Stage II, no pathway subsets or searches. '
        'No AE4 phenotype calibration or post-reveal retuning. '
        'Only the orchestrator writes canonical files, runs inference/trajectories and publishes.\n\n'
        'Restart: inspect output/publication_verification.json and checkpoint receipts; '
        'verify remote head and reuse all published results. '
        'A blocked checkpoint is a recorded gate failure, not a successful fit or simulation.\n')

def publish(name):
    receipt = json.loads((OUT/'checkpoints'/f'{name}.json').read_text())
    remote = git('ls-remote','origin','refs/heads/'+BRANCH).split()[0]
    parent = git('rev-parse','HEAD')
    if remote != parent or parent != receipt['parent_sha']:
        raise RuntimeError('Concurrent remote change; do not overwrite')
    changed = git('status','--porcelain').splitlines()
    if any('analysis/49_camp_vrac_secretory_branch_reconstruction/' not in x for x in changed):
        raise RuntimeError('Unexpected out-of-task modification')
    git('add',str(HERE.relative_to(ROOT)))
    git('commit','-m',f'Task 49 checkpoint {name}: {receipt["status"][:120]}')
    head = git('rev-parse','HEAD')
    subprocess.run(['git','push','origin','HEAD:refs/heads/'+BRANCH],cwd=ROOT,check=True)
    remote = git('ls-remote','origin','refs/heads/'+BRANCH).split()[0]
    if remote != head:
        raise RuntimeError('Remote SHA verification failed')
    pub=OUT/'publication_verification.json'
    records=json.loads(pub.read_text()) if pub.exists() else []
    records.append(dict(checkpoint=name,sha=head,remote_sha=remote,verified=True,
        verified_utc=datetime.now(timezone.utc).isoformat()))
    dump(pub,records)
    print(json.dumps(records[-1]),flush=True)

if __name__ == '__main__':
    if sys.argv[1] == 'prepare': prepare(*sys.argv[2:])
    elif sys.argv[1] == 'publish': publish(sys.argv[2])
    else: raise ValueError('prepare or publish required')
