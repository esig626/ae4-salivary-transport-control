"""File-only Task52 checkpoint receipts; never imports a scientific model."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess

ROOT = Path(__file__).resolve().parents[2]
TASK = Path(__file__).resolve().parent
BRANCH = 'analysis/task-52-chloride-reservoir-final-test'

def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT, text=True).strip()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('checkpoint', choices=['52A','52B','52C','52D','52E','52F'])
    parser.add_argument('--description', required=True)
    args = parser.parse_args()
    assert git('branch', '--show-current') == BRANCH
    receipt_path = TASK / 'output/checkpoints' / (args.checkpoint+'.json')
    assert not receipt_path.exists(), 'Never overwrite a published receipt'
    preceding = ['52A','52B','52C','52D','52E','52F'].index(args.checkpoint)
    if preceding:
        prior = ['52A','52B','52C','52D','52E'][preceding-1]
        verified = json.loads((TASK/'output'/('publication_'+prior+'.json')).read_text())
        assert verified['pass'] and verified['commit_sha'] == git('rev-parse','HEAD')
    files = sorted(q for q in TASK.rglob('*') if q.is_file()
                   and '__pycache__' not in q.parts and '.pytest_cache' not in q.parts
                   and q.suffix not in {'.pyc'} and q != receipt_path)
    files += [ROOT/'docs/MANDATORY_RESEARCH_LEDGER.md']
    receipt = dict(checkpoint=args.checkpoint, branch=BRANCH,
                   parent_commit=git('rev-parse','HEAD'), description=args.description,
                   commit_reference='The commit containing this immutable receipt',
                   artifacts={str(q.relative_to(ROOT)):hashlib.sha256(q.read_bytes()).hexdigest()
                              for q in files})
    receipt_path.write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'checkpoint':args.checkpoint,'artifact_count':len(files),
                      'receipt_sha256':hashlib.sha256(receipt_path.read_bytes()).hexdigest()}))

if __name__ == '__main__':
    main()
