"""Operational branch adaptation only; execute the byte-unchanged 52C runner."""
from pathlib import Path
import hashlib
import inspect
import subprocess
import sys

TASK = Path(__file__).resolve().parents[1]
ROOT = TASK.parents[1]
BRANCH = 'analysis/task-52D-recovery-from-52C'
BASE = '36f7c3d0bb2353bf55ac0ed9e22582720097c829'
RUNNER_SHA256 = 'd5988680b99638fa91ed0e6e82a70f6b1d6bca581ee102438060d6e1025909ce'

def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT, text=True).strip()

assert sys.argv[1:] == [], 'This launcher authorises only 52D'
assert git('branch', '--show-current') == BRANCH
assert git('rev-parse', 'HEAD') == BASE
assert hashlib.sha256((TASK/'run_frozen_cases.py').read_bytes()).hexdigest() == RUNNER_SHA256
sys.path.insert(0, str(TASK))
import run_frozen_cases as runner

# User explicitly authorises the recovery branch in place of the original.
# Replace one metadata constant in the input-verification function in memory.
# All scientific functions and every other integrity gate remain unchanged.
source = inspect.getsource(runner.verify_inputs)
old = "'analysis/task-52-chloride-reservoir-final-test'"
assert source.count(old) == 1
adapted = source.replace(old, repr(BRANCH))
assert adapted.replace(repr(BRANCH), old) == source
exec(compile(adapted, str(__file__) + ':branch_gate', 'exec'), runner.__dict__)
sys.argv = [str(TASK/'run_frozen_cases.py'), '52D']
runner.main()
