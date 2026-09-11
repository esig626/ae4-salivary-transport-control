"""Task 19 diagnostic runner; production equations and inherited inputs are immutable."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, replace
from datetime import datetime, timezone
import csv
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np

from .calibration import WTCalibrationSpec, CORE_INDEPENDENT_ROWS
from .genotype_evaluation import (
    _attempt_and_root, _cluster_roots, _core_coordinates, _regulatory_suffix,
    _state_from_core_coordinates, genotype_with_expression, simulate_genotype,
)
from .model import ModernFullModel, WT
from .nkcc_stimulation import StimulatedNkcc1Model
from .run_calcium_fast_screen import load_freeze, write_json, write_rows
from .task14_blind import build_model, routing, slug, trajectory_path, verify_contract
from .task14b_five_percent import _five_percent_core, _output_path, _finite_json
from .task16_wt_allocation import check_inputs
from .validation import (
    CONSERVATION_RESIDUAL_TOLERANCES, PRODUCTION_RADAU, attach_basal_regulation,
    sha256_file, sha256_object,
)

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / "results/19_ae4_loss_nak_pump_coupling"
BASE = "52a56cc8e8788cdb5cae60dcbaa75257de2878dd"
INTAKE = "631a77fd9aab879383080734c439cc7b4a049e5f"
BRANCH = "codex/task-19-ae4-loss-nak-pump-coupling"
EXPRESSIONS = (1.0, 0.05)
PUMPS = (1.0, 0.50, 0.10)
CALCIUM = (0.10, 0.25, 0.50)
METRICS = ("R_AE4", "R_pump_WT", "R_pump_AE4", "I", "R_coupled")


def read(path):
    return json.loads(Path(path).read_text())


def relative(path):
    return str(Path(path).relative_to(REPO))


def git(*args):
    return subprocess.check_output(["git", *args], cwd=REPO, text=True).strip()


def genotype(e):
    return WT if e == 1.0 else genotype_with_expression("AE4", e)


def pump_model(manifest, root_id, calcium, scale):
    """Scale total capacity; internal continuation scales are numerical only."""
    if not math.isfinite(scale) or not 0.10 <= scale <= 1.0:
        raise ValueError("Pump scale outside the contracted continuation interval")
    template = build_model(manifest, root_id, calcium)
    if scale == 1.0:
        return template
    membranes = replace(template.parameters.membranes,
        nak_capacity_fmol_s=template.parameters.membranes.nak_capacity_fmol_s * scale)
    return StimulatedNkcc1Model(ModernFullModel(
        parameters=replace(template.parameters, membranes=membranes),
        stimulus=template.stimulus, regulatory_model=template.regulatory_model,
        ae4_parameters=template.ae4_parameters, ae4_evaluator=template.ae4_evaluator))


def contract():
    path = OUT / "contract.json"
    return read(path), sha256_file(path)


def prepare():
    assert git("branch", "--show-current") == BRANCH
    assert not (OUT / "contract.json").exists(), "Contract is immutable"
    inherited, inherited_hash = verify_contract()
    manifest, _ = load_freeze()
    inputs = check_inputs()
    spec = WTCalibrationSpec()
    payload = {
        "contract_id": "TASK19_BASELINE_AE4_PUMP_FACTORIAL_V1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "branch": BRANCH, "baseline_commit": BASE, "intake_commit": INTAKE,
        "scientific_model_commit": inputs["baseline_commit"],
        "chloride_allocation": "inherited_baseline_only",
        "root_ids": sorted(manifest["roots"]),
        "routing_family_by_root": {r: routing(manifest, r) for r in sorted(manifest["roots"])},
        "routing_families": ["AE4NA05", "AE4NA20"],
        "ae4_expressions": list(EXPRESSIONS), "pump_scales": list(PUMPS),
        "calcium_uM": list(CALCIUM), "expected_resting_states": 60,
        "expected_dynamic_conditions": 180, "expected_new_dynamic_conditions": 120,
        "pump_scaling": {
            "changed_parameter": "membranes.nak_capacity_fmol_s",
            "rule": "target total capacity = inherited total capacity * s",
            "unchanged_partition": "membranes.apical_pump_fraction",
            "apical_capacity": "inherited total capacity * s * inherited apical fraction",
            "basolateral_capacity": "inherited total capacity * s * (1 - inherited apical fraction)",
            "stoichiometry": "3 Na cell to side : 2 K side to cell",
            "all_other_parameters_unchanged": True,
            "ae4_expression": "Genotype.ae4_expression only; complete conserved AE4 contribution"},
        "continuation": {
            "path": "For each target separately: (1,1) -> (1,s) -> (0.05,s).",
            "pump_maximum_absolute_step": 0.10, "expression_maximum_absolute_step": 0.05,
            "minimum_step_divisor": 1024, "maximum_attempted_steps_per_leg": 256,
            "retry": "Halve step after failed gate or disconnected root; double after success up to declared maximum.",
            "local_start_count": 3, "local_start_fraction": 0.001,
            "maximum_connected_step_fraction": 0.25, "max_nfev": 3000,
            "selection": "Nearest prior accepted root in max coordinate difference/span; tie by residual. Retain all local roots.",
            "coordinate_bounds": [asdict(b) for b in spec.coordinate_bounds],
            "bound_interpretation": "Inherited broad numerical search bounds, not measured physiological limits; failure is not proof of no physical root.",
            "failed_pump_parent": "Retain AE4-low arm as NOT_RUN_PARENT_REST_FAILURE; do not jump to another root.",
            "internal_points_are_scientific_conditions": False,
            "capacity_refitting": False},
        "resting_numerical_tolerances": inherited["continuation"]["numerical_tolerances"],
        "independent_rhs_scales": list(spec.independent_rhs_scales),
        "conservation_tolerances": CONSERVATION_RESIDUAL_TOLERANCES,
        "production_solver": asdict(PRODUCTION_RADAU),
        "protocol": inherited["protocol"], "time_grid_s": inherited["time_grid_s"],
        "dynamic_gates": "Complete production Radau grid; finite positive core; finite nonnegative flow; inherited conservation ratio <= 1; charge <= 1e-9 fmol.",
        "regulatory_member": inherited["regulatory_member"],
        "nkcc_normalisation_frozen": inherited["nkcc_normalisation_frozen"],
        "primary_comparisons": {
            "Q": "Trapezoidal integral of lumen outflow over the inherited 0-600 s grid, pL",
            "R_AE4": "Q(0.05,s) / Q(1,s)",
            "R_pump_WT": "Q(1,s) / Q(1,1)",
            "R_pump_AE4": "Q(0.05,s) / Q(0.05,1)",
            "I": "R_AE4(s) / R_AE4(1)",
            "R_coupled": "Q(0.05,s) / Q(1,1)",
            "reduction_percent": "100 * (1 - ratio)",
            "missing_cases": "Keep every condition; compute each ratio only when its own required inputs pass. Never impute missing Q.",
            "denominators": "Must be finite and strictly positive"},
        "interpretation": {
            "numerical_neutral_band": 0.0001,
            "material_interaction_fraction": 0.01,
            "material_rule": "Report exact ratios; 1% departure of I from 1 is a descriptive prespecified materiality threshold, not an experimental target.",
            "robust_success": "Both pump levels valid across all roots/calcium, matched AE4 ratio < 0.9999 and I < 0.99 at a tested level across all roots/calcium.",
            "coupled_reduction_alone_is_success": False,
            "experimental_phenotype_used_for_selection": False,
            "experimental_context_only_after_results_frozen": True,
            "no_pump_scales_added_after_results": True,
            "continuous_coupling_fit": False},
        "input_files": inputs["input_files"],
        "inherited_contract_sha256": inherited_hash,
        "equation_source_sha256": inherited["equation_source_sha256"],
        "root_payload_hashes": inherited["root_payload_hashes"],
        "reuse_rule": "Verify frozen file SHA256, exact parameter/state/genotype/protocol/solver/grid and flow integral; reuse all valid scale-1 files.",
        "no_archive_edits": True,
    }
    assert len(payload["root_ids"]) == 10
    write_json(OUT / "contract.json", payload)
    print("Prepared immutable contract for 60 resting and 180 dynamic conditions", flush=True)


def require_contract_commit():
    c, digest = contract()
    assert git("branch", "--show-current") == BRANCH
    assert git("rev-parse", "HEAD:" + relative(OUT / "contract.json")) == git("hash-object", relative(OUT / "contract.json"))
    assert c["pump_scales"] == list(PUMPS) and c["ae4_expressions"] == list(EXPRESSIONS)
    assert c["calcium_uM"] == list(CALCIUM) and c["baseline_commit"] == BASE
    check_inputs()
    return c, digest


def diagnostics(model, state, expression, time_s=0.0):
    ev = model.evaluate(time_s, state, genotype=genotype(expression))
    d = ev.diagnostics
    m, h, a, o = d.membranes, d.homeostasis, d.ae4, d.observables
    row = {"time_s": float(time_s), "ae4_expression": expression,
        "volume_i_pL": float(state[5]), "volume_l_pL": float(state[11]),
        "ph_i": o.cell_acid_base.ph, "ph_l": o.lumen_acid_base.ph,
        "v_apical_mV": 1000*m.v_apical_V, "v_basolateral_mV": 1000*m.v_basolateral_V,
        "v_transepithelial_mV": 1000*m.v_transepithelial_V,
        "pump_total_fmol_s": m.pump_apical_fmol_s + m.pump_basolateral_fmol_s,
        "pump_apical_fmol_s": m.pump_apical_fmol_s,
        "pump_basolateral_fmol_s": m.pump_basolateral_fmol_s,
        "pump_total_capacity_fmol_s": model.parameters.membranes.nak_capacity_fmol_s,
        "apical_pump_capacity_fraction": model.parameters.membranes.apical_pump_fraction,
        "ae4_na_branch_fmol_s": a.diagnostics["j_na_fmol_s"],
        "ae4_k_branch_fmol_s": a.diagnostics["j_k_fmol_s"],
        "ae4_net_cl_source_fmol_s": a.cl_cell_fmol_s,
        "ae4_na_source_fmol_s": a.na_cell_fmol_s, "ae4_k_source_fmol_s": a.k_cell_fmol_s,
        "ae4_hco3_source_fmol_s": a.hco3_cell_fmol_s,
        "nkcc1_inward_fmol_s": h.nkcc1_inward_fmol_s,
        "ae2_inward_fmol_s": h.ae2_inward_fmol_s,
        "nhe1_inward_fmol_s": h.nhe1_inward_fmol_s,
        "flow_pL_s": d.water.lumen_outflow_pL_s,
        "max_abs_charge_fmol": max(abs(v) for v in d.state_charge_fmol.values()),
        "max_abs_current_A": max(abs(v) for v in m.current_residuals_A.values()),
        "max_abs_scaled_independent_rhs": float(np.max(np.abs(ev.rhs[list(CORE_INDEPENDENT_ROWS)] / np.asarray(WTCalibrationSpec().independent_rhs_scales)))),
        "max_dimensionless_conservation_ratio": max(abs(v)/CONSERVATION_RESIDUAL_TOLERANCES[k] for k,v in d.conservation_residuals.items())}
    for side, concentrations in (("i", o.cell_concentrations_mM), ("l", o.lumen_concentrations_mM)):
        for name in ("na", "k", "cl", "tic", "hco3"):
            row[f"{name}_{side}_mM"] = concentrations[name]
    for k,v in d.conservation_residuals.items():
        row["conservation__" + k] = float(v)
    for k,v in m.currents_A.items():
        row["current_A__" + k] = float(v)
    for name,value in zip(model.state_names, ev.rhs):
        row["rhs__" + name] = float(value)
    return _finite_json(row)


def local_solve(model, previous, e):
    spec = WTCalibrationSpec()
    lower, upper = spec.coordinate_lower, spec.coordinate_upper
    span = upper-lower
    starts = [("previous_connected_root", previous)] + [(name, np.clip(previous + direction*.001*span,
        lower+1e-10*span, upper-1e-10*span)) for name,direction in (("local_minus",-1),("local_plus",1))]
    attempts, roots = [], []
    for name,start in starts:
        attempt, root = _attempt_and_root(model=model, transporter="AE4", expression=e,
            genotype=genotype(e), regulatory_suffix=_regulatory_suffix(model), spec=spec,
            start_id=name, start=start, lower=lower, upper=upper, max_nfev=3000)
        saved = asdict(attempt)
        # Preserve residuals at rejected candidates as well as accepted roots.
        if root is None:
            try:
                candidate = _state_from_core_coordinates(model, attempt.coordinates, _regulatory_suffix(model))
                saved["rejected_candidate_diagnostics"] = diagnostics(model, candidate, e)
                saved["boundary_proximity"] = {b.name: float(min((x-b.lower)/(b.upper-b.lower),
                    (b.upper-x)/(b.upper-b.lower))) for x,b in zip(attempt.coordinates, spec.coordinate_bounds)}
            except (ValueError, FloatingPointError, OverflowError, np.linalg.LinAlgError) as exc:
                saved["diagnostic_error"] = str(exc)
        attempts.append(saved)
        if root is not None:
            roots.append(root)
    roots = _cluster_roots(roots, expression=e, transporter="AE4", span=span,
        relative_tolerance=spec.cluster_relative_tolerance)
    selected = min(roots, key=lambda r:(float(np.max(np.abs(np.asarray(r.coordinates)-previous)/span)),
        r.max_abs_scaled_independent_rhs)) if roots else None
    distance = float(np.max(np.abs(np.asarray(selected.coordinates)-previous)/span)) if selected else None
    connected = selected is not None and distance <= .25
    if selected:
        d = diagnostics(model, attach_basal_regulation(model, selected.core_state), e)
        connected = connected and d["max_dimensionless_conservation_ratio"] <= 1
    return {"attempts": attempts, "roots": [asdict(r) for r in roots],
        "selected_root": asdict(selected) if selected else None,
        "distance_from_previous": distance, "connected": bool(connected)}


def continue_leg(manifest, root, core, target, fixed_pump=None):
    leg = "pump" if fixed_pump is None else "AE4"
    maximum = .10 if leg == "pump" else .05
    current, step = 1.0, maximum
    model = pump_model(manifest, root, .10, 1.0 if fixed_pump is None else fixed_pump)
    previous = _core_coordinates(model, core)
    history, selected = [], None
    while current > target + 1e-12 and len(history) < 256:
        proposed = max(target, round(current-step, 14))
        s, e = (proposed, 1.0) if leg == "pump" else (fixed_pump, proposed)
        model = pump_model(manifest, root, .10, s)
        report = local_solve(model, previous, e)
        report.update(leg=leg, previous_value=current, proposed_value=proposed,
            ae4_expression=e, pump_scale=s, numerical_point_only=True)
        history.append(report)
        if report["connected"]:
            current = proposed
            selected = report["selected_root"]
            previous = np.asarray(selected["coordinates"])
            step = min(maximum, step*2)
        else:
            step /= 2
            if step < maximum/1024:
                break
    reached = abs(current-target) <= 1e-12
    return {"leg": leg, "target": target, "fixed_pump": fixed_pump,
        "status": "COMPLETE" if reached else "REST_CONTINUATION_FAILURE",
        "reached_target": reached, "last_connected_value": current,
        "selected_root": selected if reached else None,
        "last_connected_root": selected, "steps": history}


def rest_path(root, s):
    return OUT / "resting" / f"{slug(root)}_pump{s:.2f}.json"


def rest_job(job):
    root,s = job
    c,digest = contract()
    path = rest_path(root,s)
    if path.exists():
        data = read(path)
        assert data["contract_sha256"] == digest
        return data
    manifest,_ = load_freeze()
    started = time.perf_counter()
    original = manifest["roots"][root]["core_state"]
    pump = continue_leg(manifest,root,original,s)
    near = continue_leg(manifest,root,pump["selected_root"]["core_state"],.05,fixed_pump=s) if pump["reached_target"] else {
        "status":"NOT_RUN_PARENT_REST_FAILURE", "reached_target":False, "selected_root":None, "steps":[]}
    result = {"root_id":root,"routing_family":routing(manifest,root),"pump_scale":s,
        "contract_sha256":digest,"pump_leg":pump,"ae4_leg":near,"wall_seconds":time.perf_counter()-started}
    write_json(path,_finite_json(result))
    return result


def baseline_core(manifest,root,e):
    return manifest["roots"][root]["core_state"] if e == 1.0 else _five_percent_core(root)


def resting_rows(manifest):
    rows=[]
    for root in sorted(manifest["roots"]):
        for s in PUMPS:
            saved = read(rest_path(root,s)) if s < 1 else None
            for e in EXPRESSIONS:
                leg = saved["pump_leg" if e == 1 else "ae4_leg"] if saved else None
                valid = leg["reached_target"] if leg else True
                core = (leg["selected_root"]["core_state"] if valid else None) if leg else baseline_core(manifest,root,e)
                row={"root_id":root,"routing_family":routing(manifest,root),"pump_scale":s,
                    "ae4_expression":e,"status":leg["status"] if leg else "BASELINE_HASH_VALID_REUSE",
                    "resting_gate_pass":valid,"core_state_sha256":sha256_object(core) if valid else None,
                    "resting_source_path":relative(rest_path(root,s)) if saved else "inherited_frozen_root",
                    "fixed_parameter_rank":leg["selected_root"]["normalized_jacobian_rank"] if leg and valid else (10 if not leg else None)}
                if valid:
                    model=pump_model(manifest,root,.10,s)
                    row.update(diagnostics(model,attach_basal_regulation(model,core),e))
                rows.append(row)
    return rows


def verify_reuse(manifest,root,e,ca):
    c,digest=contract()
    path=trajectory_path(root,ca,"WT","production_radau") if e == 1 else _output_path(root,ca)
    expected=c["input_files"][relative(path)]
    assert sha256_file(path)==expected["sha256"] and path.stat().st_size==expected["bytes"]
    saved=read(path)
    trace=saved["trajectory"]
    model=pump_model(manifest,root,ca,1.0)
    assert saved["root_id"]==root and saved["calcium_uM"]==ca and saved["status"]=="COMPLETE"
    assert saved["whole_cell_parameters_sha256"]==sha256_object(model.parameters)
    assert saved["ae4_parameters_sha256"]==sha256_object(model.ae4_parameters)
    assert trace["solver"]==asdict(PRODUCTION_RADAU) and trace["time_s"]==c["time_grid_s"]
    assert trace["numerical_gate_pass"] and trace["success"]
    assert np.array_equal(np.asarray(trace["states"])[:12,0], baseline_core(manifest,root,e))
    if e == 1:
        assert saved["genotype"]==trace["genotype_name"]=="WT"
    else:
        assert saved["genotype"]==asdict(genotype(e))
        assert saved["selected_core_state_sha256"]==sha256_object(baseline_core(manifest,root,e))
        assert saved["wt_denominator_sha256"]==sha256_file(trajectory_path(root,ca,"WT","production_radau"))
    assert math.isclose(trace["cumulative_flow_pL"][-1],float(np.trapezoid(trace["flow_pL_s"],trace["time_s"])),rel_tol=1e-12)
    return path,saved


def dynamic_path(root,e,s,ca):
    return OUT / "trajectories" / f"{slug(root)}_e{e:.2f}_pump{s:.2f}_ca{ca:.2f}.json"


def trace_summary(trace,model):
    states=np.asarray(trace["states"])
    charge=max(float(np.max(np.abs(states[0]+states[1]-states[2]-states[4]-model.parameters.geometry.fixed_cell_anion_equivalents_fmol))),
        float(np.max(np.abs(states[6]+states[7]-states[8]-states[10]))))
    c,_=contract()
    valid=bool(trace["success"] and trace["numerical_gate_pass"] and trace["time_s"]==c["time_grid_s"] and charge <= 1e-9)
    row={"dynamic_gate_pass":valid,"Q_0_600_pL":trace["cumulative_flow_pL"][-1] if valid else None,
        "max_abs_charge_fmol":charge,"max_dimensionless_conservation_ratio":trace["max_dimensionless_conservation_ratio"],
        "solver_success":trace["success"],"positive_core":trace["positive_core"],
        "all_flow_nonnegative":trace["all_flow_nonnegative"],"sample_count":len(trace["time_s"]),
        "solver_message":trace["message"]}
    row.update({"max_abs_conservation__"+k:v for k,v in trace["max_abs_conservation_residuals"].items()})
    return row


def dynamic_job(job):
    root,e,s,ca=job
    assert e in EXPRESSIONS and s in PUMPS and ca in CALCIUM
    manifest,_=load_freeze()
    c,digest=contract()
    model=pump_model(manifest,root,ca,s)
    row={"root_id":root,"routing_family":routing(manifest,root),"ae4_expression":e,
        "pump_scale":s,"calcium_uM":ca,"new_trajectory":s<1,"contract_sha256":digest}
    if s == 1:
        path,data=verify_reuse(manifest,root,e,ca)
        status="BASELINE_HASH_VALID_REUSE"
    else:
        rest=read(rest_path(root,s))["pump_leg" if e == 1 else "ae4_leg"]
        if not rest["reached_target"]:
            return {**row,"status":"NOT_RUN_REST_FAILURE","rest_status":rest["status"],
                "dynamic_gate_pass":False,"new_trajectory":False,"Q_0_600_pL":None}
        path=dynamic_path(root,e,s,ca)
        if path.exists():
            data=read(path)
            assert data["contract_sha256"]==digest and data["root_id"]==root
            assert data["whole_cell_parameters_sha256"]==sha256_object(model.parameters)
            assert data["core_state_sha256"]==sha256_object(rest["selected_root"]["core_state"])
        else:
            started=time.perf_counter()
            trace=simulate_genotype(model,rest["selected_root"]["core_state"],transporter="AE4",
                genotype=genotype(e),solver=PRODUCTION_RADAU,time_s=c["time_grid_s"])
            data={**row,"whole_cell_parameters_sha256":sha256_object(model.parameters),
                "ae4_parameters_sha256":sha256_object(model.ae4_parameters),"genotype":asdict(genotype(e)),
                "core_state_sha256":sha256_object(rest["selected_root"]["core_state"]),
                "trajectory":asdict(trace),"wall_seconds":time.perf_counter()-started}
            write_json(path,_finite_json(data))
            data=read(path)
        status="COMPLETE"
    row.update(trace_summary(data["trajectory"],model))
    row.update(status=status if row["dynamic_gate_pass"] else "DYNAMIC_NUMERICAL_FAILURE",
        trajectory_path=relative(path),trajectory_sha256=sha256_file(path))
    return _finite_json(row)


def ratios(q11,q051,q1s,q05s):
    def ratio(a,b):
        if a is None or b is None or not math.isfinite(a) or not math.isfinite(b) or b<=0:
            return None
        return a/b
    ra=ratio(q05s,q1s)
    baseline=ratio(q051,q11)
    return {"R_AE4":ra,"R_pump_WT":ratio(q1s,q11),"R_pump_AE4":ratio(q05s,q051),
        "I":ratio(ra,baseline),"R_coupled":ratio(q05s,q11)}


def factorial(rows,manifest):
    index={(r["root_id"],r["ae4_expression"],r["pump_scale"],r["calcium_uM"]):r for r in rows}
    expected={(r,e,s,c) for r in manifest["roots"] for e in EXPRESSIONS for s in PUMPS for c in CALCIUM}
    assert len(index)==len(rows)==180 and set(index)==expected
    output=[]
    for root in sorted(manifest["roots"]):
        for ca in CALCIUM:
            for s in PUMPS:
                arms=[index[root,e,scale,ca] for e,scale in ((1.,1.),(.05,1.),(1.,s),(.05,s))]
                q=[a.get("Q_0_600_pL") if a["dynamic_gate_pass"] else None for a in arms]
                row={"root_id":root,"routing_family":routing(manifest,root),"calcium_uM":ca,"pump_scale":s,
                    "complete_factorial_comparison":all(a["dynamic_gate_pass"] for a in arms),
                    **dict(zip(("Q_1_1","Q_005_1","Q_1_s","Q_005_s"),q)),**ratios(*q)}
                for label,arm in zip(("WT_normal","AE4_low_normal","WT_pump","AE4_low_pump"),arms):
                    row[label+"_status"]=arm["status"]
                for metric in METRICS:
                    row[metric+"_reduction_percent"]=100*(1-row[metric]) if row[metric] is not None else None
                output.append(row)
    return output


def pool(function,jobs,workers):
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures={executor.submit(function,j):j for j in jobs}
        for f in as_completed(futures):
            result=f.result()
            print(json.dumps({"job":futures[f],"status":result.get("status",result.get("pump_leg",{}).get("status")),
                "ae4_status":result.get("ae4_leg",{}).get("status")}),flush=True)
            yield result


def run_rests(workers):
    require_contract_commit()
    manifest,_=load_freeze()
    list(pool(rest_job,[(r,s) for r in sorted(manifest["roots"]) for s in PUMPS if s<1],workers))
    write_rows(OUT/"resting_states.csv",resting_rows(manifest))


def run_dynamics(workers):
    require_contract_commit()
    manifest,_=load_freeze()
    jobs=[(r,e,s,c) for r in sorted(manifest["roots"]) for s in PUMPS for e in EXPRESSIONS for c in CALCIUM]
    rows=list(pool(dynamic_job,jobs,workers))
    rows.sort(key=lambda r:(r["root_id"],-r["pump_scale"],-r["ae4_expression"],r["calcium_uM"]))
    write_rows(OUT/"dynamic_results.csv",rows)
    write_rows(OUT/"factorial_comparisons.csv",factorial(rows,manifest))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("stage",choices=("prepare","rest","dynamic"))
    parser.add_argument("--workers",type=int,default=9)
    args=parser.parse_args()
    if args.stage=="prepare": prepare()
    elif args.stage=="rest": run_rests(args.workers)
    else: run_dynamics(args.workers)


if __name__=="__main__":
    main()
