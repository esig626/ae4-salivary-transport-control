"""Verify published scientific records without repeating any model solve."""
from pathlib import Path
import hashlib
import json
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
START="5c945c4a7857657b5f7c76d5bf9e3d3abc7d5721"
BRANCH="analysis/task-48-joint-experimental-constraint-reconstruction"


def git(*args):
    return subprocess.check_output(["git",*args],cwd=ROOT)


def main():
    assert git("branch","--show-current").decode().strip()==BRANCH
    assert not git("diff",START,"--",".",":(exclude)analysis/48_joint_experimental_constraint_reconstruction")
    originals={}
    for path in git("ls-tree","-r","--name-only",START).decode().splitlines():
        if not path.startswith("analysis/48_joint_experimental_constraint_reconstruction/"):
            assert (ROOT/path).read_bytes()==git("show",START+":"+path),path
            originals[path]=True
    manifest=json.loads((HERE/"output/executed_source_hashes.json").read_text())
    for path,digest in manifest.items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,path
    assert (HERE/"output/final_residuals.csv").read_bytes()==(HERE/"output/baseline_residuals.csv").read_bytes()
    assert (HERE/"output/final_residuals.json").read_bytes()==(HERE/"output/baseline_residuals.json").read_bytes()
    checkpoints=[]
    for path in sorted((HERE/"output/checkpoints").glob("48?.json")):
        rel=str(path.relative_to(ROOT))
        history=git("log","--format=%H","--",rel).decode().splitlines()
        # An uncommitted next receipt cannot yet have a commit identity.
        if not history: continue
        assert len(history)==1, "Published checkpoint receipt was rewritten"
        sha=history[0]
        record=json.loads(git("show",sha+":"+rel))
        assert git("show","-s","--format=%P",sha).decode().strip()==record["parent_sha"]
        for item in record["produced_files"]:
            assert hashlib.sha256(git("show",sha+":"+item["path"])).hexdigest()==item["sha256"],item["path"]
        checkpoints.append(dict(milestone=record["milestone"],sha=sha,parent_sha=record["parent_sha"],file_hashes_valid=True))
    rests=[json.loads(x.read_text()) for x in (HERE/"output/rests").glob("*.json")]
    trajectories=[json.loads(x.read_text()) for x in (HERE/"output/protocol_predictions").glob("*.json")]
    assert len(rests)==3 and len(trajectories)==12
    vector=json.loads((HERE/"output/parameter_vector.json").read_text())
    assert vector["shared_scientific_vector"]==json.loads((ROOT/"analysis/47_dynamic_potassium_recycling_reconstruction/input/active_parameters.json").read_text())
    expected=vector["shared_payload_sha256"]
    for record in rests+trajectories:
        assert record["shared_parameter_hash"]==expected
        assert record["positive_core"]
    for record in trajectories:
        root=next(x for x in rests if x["genotype"]==record["genotype"])
        assert record["initial_state"]==root["state"]
    result=dict(passed=True,preserved_original_files=len(originals),unchanged_executed_source_files=len(manifest),
        verified_published_checkpoint_receipts=checkpoints,positive_genotype_rests=3,positive_protocol_trajectories=12,
        shared_parameter_hash=expected,exact_genotype_initial_states=True,final_residuals_equal_baseline=True,
        new_scientific_solves_in_verification=0,
        final_receipt_note="A newly created 48F receipt is checked by its containing commit and remote verification after publication.")
    print(json.dumps(result,indent=2))
    return result


if __name__=="__main__": main()
