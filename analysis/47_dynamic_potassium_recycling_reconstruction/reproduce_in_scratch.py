"""Optional explicit trajectory replay without changing any frozen output.

Not executed during Task 47. Run with wt, ae4_null or ae4_5pct.
"""
import argparse
import csv
import gzip
import json
from pathlib import Path
import shutil
import tempfile
import numpy as np
import run_diagnostics as run

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('case',choices=['wt','ae4_null','ae4_5pct'])
    case=parser.parse_args().case
    original=run.HERE
    with tempfile.TemporaryDirectory(prefix='ae4_task47_reproduction_') as directory:
        output=Path(directory)
        for name in ['PREDICTION_MANIFEST.json','REMOTE_PREDICTION_RECEIPT.json']:
            shutil.copyfile(original/name,output/name)
        # load_model and source_hashes belong to common.py and retain their
        # canonical input paths. Only this runner's output directory moves.
        run.HERE=output
        run.run_case(case)
        with gzip.open(original/'output'/case/'states.csv.gz','rt') as handle:
            expected=list(csv.DictReader(handle))
        with (output/'output'/case/'states.csv').open() as handle:
            actual=list(csv.DictReader(handle))
        fields=list(expected[0]);a=np.array([[float(r[k]) for k in fields] for r in actual]);b=np.array([[float(r[k]) for k in fields] for r in expected])
        assert a.shape==b.shape
        assert np.allclose(a,b,rtol=2e-6,atol=1e-9)
        print(json.dumps(dict(case=case,archived_prediction_reproduced=True,max_absolute_state_difference=float(np.max(abs(a-b))),canonical_outputs_changed=False)))
    run.HERE=original

if __name__=='__main__':main()
