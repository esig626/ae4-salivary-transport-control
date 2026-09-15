"""Read-only guard for Task 46 publication and immutable source boundaries."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]
AREA = 'analysis/46_physiology_constrained_model_reconstruction/'
BRANCH = 'analysis/task-46-physiology-constrained-model-reconstruction'

def git(*args):
    return subprocess.check_output(['git', '-C', str(ROOT), *args], text=True).strip()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--published', action='store_true')
    args = parser.parse_args()
    pins = json.loads((ROOT / AREA / 'source_pins.json').read_text())
    assert git('branch', '--show-current') == BRANCH
    for source in pins['sources'].values():
        assert git('cat-file', '-t', source) == 'commit'
    base = pins['base_sha']
    changed = set(git('diff', '--name-only', base).splitlines())
    untracked = set(git('ls-files', '--others', '--exclude-standard').splitlines())
    assert all(p.startswith(AREA) for p in changed | untracked), sorted(changed | untracked)
    for name, blob in pins['specification_blobs'].items():
        data = (ROOT / name).read_bytes()
        actual = hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()
        assert actual == blob, name
    head = git('rev-parse', 'HEAD')
    if args.published:
        assert not git('status', '--porcelain'), 'Working tree must be clean at publication.'
        remote = git('ls-remote', '--heads', 'origin', 'refs/heads/' + BRANCH).split()[0]
        assert remote == head, (remote, head)
    print(json.dumps({'status': 'passed', 'head': head, 'published': args.published,
                      'changed_paths': len(changed), 'untracked_paths': len(untracked),
                      'protected_base_paths_unchanged': True}, sort_keys=True))

if __name__ == '__main__':
    main()
