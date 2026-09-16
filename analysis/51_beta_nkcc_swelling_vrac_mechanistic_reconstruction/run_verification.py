"""Run once after verified51A; save all test outcomes and core-evaluation count."""
import contextlib
import io
import json
import unittest
import sys
from unittest.mock import patch
from checkpoint import HERE,OUT,verified,integrity,dump,sha

def main():
    verified('51A');integrity()
    attempt=sys.argv[1] if len(sys.argv)>1 else '01'
    resultfile=OUT/f'implementation_verification_attempt_{attempt}.json'
    if resultfile.exists():raise RuntimeError('Do not silently rerun verification')
    from test_task51 import Task51Tests
    from modern_full_model.model import ModernFullModel
    original=ModernFullModel.evaluate
    count=0
    def counted(self,*args,**kwargs):
        nonlocal count
        count+=1
        return original(self,*args,**kwargs)
    output=io.StringIO()
    with patch.object(ModernFullModel,'evaluate',counted):
        suite=(unittest.TestSuite([Task51Tests(sys.argv[2])]) if len(sys.argv)>2
               else unittest.defaultTestLoader.loadTestsFromTestCase(Task51Tests))
        result=unittest.TextTestRunner(stream=output,verbosity=2).run(suite)
    (OUT/f'implementation_tests_attempt_{attempt}.txt').write_text(output.getvalue())
    record=dict(passed=result.wasSuccessful(),tests=result.testsRun,failures=len(result.failures),
        errors=len(result.errors),core_evaluations=count,production_trajectories=0,parameter_fits=0,
        source_sha256={p.name:sha(p) for p in [HERE/'task51_model.py',HERE/'test_task51.py']})
    dump(resultfile,record);print(output.getvalue());print(json.dumps(record))
    if not result.wasSuccessful():raise SystemExit(1)

if __name__=='__main__':main()
