"""Read-only numerical postprocessing followed by an immutable Task 19 result freeze."""
from __future__ import annotations

from collections import Counter
import csv
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np

from . import task19_pump_coupling as task
from .run_calcium_fast_screen import write_json, write_rows
from .validation import sha256_file


def table(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def decoded_results():
    rows=table(task.OUT/"dynamic_results.csv")
    for row in rows:
        for k in ("ae4_expression","pump_scale","calcium_uM","Q_0_600_pL"):
            row[k]=float(row[k]) if row[k] else None
        for k in ("dynamic_gate_pass","new_trajectory"):
            assert row[k] in ("True","False")
            row[k]=row[k]=="True"
    return rows


def grouped_comparisons(comparisons):
    groups=[]
    for s in task.PUMPS:
        for family in ("ALL","AE4NA05","AE4NA20"):
            for ca in ("ALL",*task.CALCIUM):
                selected=[r for r in comparisons if r["pump_scale"]==s and
                    (family=="ALL" or r["routing_family"]==family) and
                    (ca=="ALL" or r["calcium_uM"]==ca)]
                row={"pump_scale":s,"routing_family":family,"calcium_uM":ca,
                    "expected_pairs":len(selected),"valid_pairs":sum(r["complete_factorial_comparison"] for r in selected)}
                for metric in task.METRICS:
                    values=[r[metric] for r in selected if r[metric] is not None]
                    row[metric+"_count"]=len(values)
                    for name,operation in (("min",min),("median",np.median),("max",max)):
                        row[metric+"_"+name]=float(operation(values)) if values else None
                groups.append(row)
    return groups


def state_flux_table(manifest,rows):
    output=[]
    for row in rows:
        common={k:row[k] for k in ("root_id","routing_family","ae4_expression","pump_scale","calcium_uM","status")}
        common["source_path"]=row.get("trajectory_path","")
        common["source_sha256"]=row.get("trajectory_sha256","")
        if not row.get("trajectory_path"):
            output.append({**common,"phase":"UNAVAILABLE","time_s":None})
            continue
        trace=task.read(task.REPO/row["trajectory_path"])["trajectory"]
        model=task.pump_model(manifest,row["root_id"],row["calcium_uM"],row["pump_scale"])
        times=np.asarray(trace["time_s"])
        states=np.asarray(trace["states"])
        for t in (0.,60.,180.,600.):
            indices=np.flatnonzero(times==t)
            if len(indices)!=1:
                output.append({**common,"phase":"UNAVAILABLE_TIMEPOINT","time_s":t})
                continue
            d=task.diagnostics(model,states[:,indices[0]],row["ae4_expression"],time_s=t)
            output.append({**common,"phase":"REST" if t==0 else "DYNAMIC",**d})
    return output


def continuation_failure_table(manifest):
    output=[]
    for root in sorted(manifest["roots"]):
        for s in task.PUMPS[1:]:
            saved=task.read(task.rest_path(root,s))
            for name in ("pump_leg","ae4_leg"):
                leg=saved[name]
                if leg["status"]=="COMPLETE":
                    continue
                row={"root_id":root,"routing_family":task.routing(manifest,root),
                    "target_pump_scale":s,"ae4_expression":1. if name=="pump_leg" else .05,
                    "leg":name,"status":leg["status"],"last_connected_value":leg.get("last_connected_value"),
                    "attempted_step_count":len(leg["steps"]),"source_path":task.relative(task.rest_path(root,s))}
                if leg["steps"]:
                    last=leg["steps"][-1]
                    attempt=last["attempts"][0]
                    row["last_rejected_value"]=last["proposed_value"]
                    row["last_optimizer_message"]=attempt["message"]
                    row["last_candidate_scaled_residual"]=attempt["max_abs_scaled_independent_rhs"]
                    near={k:v for k,v in attempt.get("boundary_proximity",{}).items() if v<=1e-5}
                    row["search_boundary_hits"]=json.dumps(near,sort_keys=True)
                    d=attempt.get("rejected_candidate_diagnostics",{})
                    for k in ("na_i_mM","k_i_mM","cl_i_mM","ph_i","volume_i_pL","max_abs_charge_fmol","max_abs_current_A","max_dimensionless_conservation_ratio"):
                        row["last_rejected_"+k]=d.get(k)
                output.append(row)
    return output


def freeze():
    c,digest=task.require_contract_commit()
    assert not (task.OUT/"results_frozen.json").exists(),"Results are immutable"
    manifest,_=task.load_freeze()
    rows=decoded_results()
    comparisons=task.factorial(rows,manifest)
    # Independently regenerate the already written primary comparison table.
    old=table(task.OUT/"factorial_comparisons.csv")
    assert len(old)==len(comparisons)==90
    for a,b in zip(old,comparisons):
        for key in task.METRICS:
            assert (float(a[key]) if a[key] else None)==b[key]
    groups=grouped_comparisons(comparisons)
    write_rows(task.OUT/"factorial_summary.csv",groups)
    write_rows(task.OUT/"state_flux_diagnostics.csv",state_flux_table(manifest,rows))
    failures=continuation_failure_table(manifest)
    write_rows(task.OUT/"continuation_failures.csv",failures)
    rests=table(task.OUT/"resting_states.csv")
    assert len(rests)==60
    numerical={"contract_sha256":digest,"condition_count":len(rows),
        "hash_valid_reused_trajectories":sum(r["status"]=="BASELINE_HASH_VALID_REUSE" for r in rows),
        "new_trajectories":sum(r["new_trajectory"] for r in rows),
        "valid_dynamic_conditions":sum(r["dynamic_gate_pass"] for r in rows),
        "dynamic_status_counts":dict(Counter(r["status"] for r in rows)),
        "rest_status_counts":dict(Counter(r["status"] for r in rests)),
        "by_pump":[],"comparison_groups":groups,
        "failures_are_retained":True,"missing_values_imputed":False,
        "experimental_data_used":False,"model_or_parameter_fit_performed":False}
    for s in task.PUMPS:
        r=[x for x in rows if x["pump_scale"]==s]
        st=[x for x in rests if float(x["pump_scale"])==s]
        numerical["by_pump"].append({"pump_scale":s,"resting_expected":len(st),
            "resting_valid":sum(x["resting_gate_pass"]=="True" for x in st),
            "dynamic_expected":len(r),"dynamic_valid":sum(x["dynamic_gate_pass"] for x in r),
            "paired_root_calcium_count":sum(x["complete_factorial_comparison"] for x in comparisons if x["pump_scale"]==s)})
    write_json(task.OUT/"numerical_summary.json",numerical)
    inputs=task.check_inputs()
    # Git trees prove that the production model and archive were not edited.
    protected={p:{"baseline_tree":task.git("rev-parse",task.BASE+":"+p),
        "working_changes":task.git("diff",task.BASE,"--",p)} for p in ("archive","model")}
    assert all(not v["working_changes"] for v in protected.values())
    for p,h in c["equation_source_sha256"].items():
        assert sha256_file(task.REPO/p)==h,p
    record={"created_utc":datetime.now(timezone.utc).isoformat(),"contract_sha256":digest,
        "contract_commit":"74159c95bedc739d63c0399c5f24166ad65de11d",
        "experimental_context_not_read_during_computation":True,
        "frozen_input_count":len(inputs["input_files"]),"all_187_frozen_inputs_unchanged":True,
        "production_equation_files_unchanged":True,"protected_paths":protected,
        "runner_sha256":sha256_file(Path(task.__file__)),"reporter_sha256":sha256_file(Path(__file__)),
        "hashes":{task.relative(p):sha256_file(p) for p in sorted(task.OUT.rglob("*")) if p.is_file()}}
    write_json(task.OUT/"results_frozen.json",record)
    print(json.dumps({k:v for k,v in numerical.items() if k!="comparison_groups"},indent=2))


def correct_branch_reporting():
    """Explicit reporting-only correction; preserve the first freeze and all Q hashes."""
    original=task.OUT/"results_frozen.json"
    prior=task.OUT/"results_frozen_v1.json"
    assert original.exists() and not prior.exists()
    prior.write_bytes(original.read_bytes())
    old=task.read(prior)
    manifest,_=task.load_freeze()
    write_rows(task.OUT/"resting_states.csv",task.resting_rows(manifest))
    write_rows(task.OUT/"state_flux_diagnostics.csv",state_flux_table(manifest,decoded_results()))
    changed={task.relative(task.OUT/n) for n in ("resting_states.csv","state_flux_diagnostics.csv")}
    for name,digest in old["hashes"].items():
        if name not in changed:
            assert sha256_file(task.REPO/name)==digest,name
    record={**old,"diagnostic_correction_utc":datetime.now(timezone.utc).isoformat(),
        "prior_freeze_path":task.relative(prior),"prior_freeze_sha256":sha256_file(prior),
        "diagnostic_correction":{
            "reason":"Production AE4 diagnostic j_na/j_k values precede expression scaling. Final top-level tables now label actual genotype-scaled branch currents and retain the original values in before_expression columns.",
            "scope":"Only the two top-level diagnostic CSV files changed. Original continuation logs retain the pre-expression branch fields under the original runner version.",
            "changed_files":sorted(changed),"new_resting_solves":0,"new_trajectories":0,
            "all_secretion_and_factorial_hashes_unchanged":True},
        "simulation_runner_sha256":old["runner_sha256"],
        "runner_sha256":sha256_file(Path(task.__file__)),
        "reporter_sha256":sha256_file(Path(__file__)),
        "hashes":{**old["hashes"],**{p:sha256_file(task.REPO/p) for p in changed},task.relative(prior):sha256_file(prior)}}
    write_json(original,record)
    print("Corrected branch-flux reporting; all numerical trajectories and comparison hashes unchanged")


if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("stage",choices=("freeze","correct-branch-reporting"),nargs="?",default="freeze")
    args=parser.parse_args()
    freeze() if args.stage=="freeze" else correct_branch_reporting()
