"""Authoritative Task 48C diagnostics, with reusable per-result checkpoints.

This does not optimise parameters. A structural contrast precedes numerical
diagnostics. Every established genotype has a solved stable resting state.
"""
from __future__ import annotations
import os
for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"
import csv
import json
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import root
from scipy.integrate import solve_ivp
from protocol_layer import *

OUT = HERE / "output"
KEEP = np.array([0, 1, 2, 3, 5, 6, 7, 8, 9, 11])


def dump(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + "\n")


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def observable(model, t, state, genotype):
    e = evaluate(model, t, state, genotype)
    d = e.diagnostics
    c = d.observables.cell_concentrations_mM
    return dict(time_s=float(t), cl_i_mM=float(c["cl"]), na_i_mM=float(c["na"]),
        k_i_mM=float(c["k"]), ph_i=float(d.observables.cell_acid_base.ph),
        volume_pL=float(state[5]), q_pL_s=float(d.water.lumen_outflow_pL_s),
        hco3_l_mM=float(d.observables.lumen_acid_base.hco3_mM),
        nkcc_cycle_fmol_s=float(d.homeostasis.nkcc1_inward_fmol_s),
        nhe_fmol_s=float(d.homeostasis.nhe1_inward_fmol_s),
        ae2_fmol_s=float(d.homeostasis.ae2_inward_fmol_s),
        ae4_fmol_s=float(d.ae4.cl_cell_fmol_s),
        nbc_fmol_s=float(d.regulatory["minimal_nbc_cycle_inward_fmol_s"]),
        dcl_dt_mM_s=float(concentration_derivative(state, e.rhs)),
        max_charge_fmol=float(max(abs(v) for v in d.state_charge_fmol.values())),
        max_current_residual_A=float(max(abs(v) for v in d.membranes.current_residuals_A.values())))


def resting_state(genotype_name):
    model = build_model("REST")
    genotype = GENOTYPES[genotype_name]
    key = rest_cache_key(model, genotype)
    path = OUT / "rests" / f"{genotype_name}.json"
    if path.exists():
        saved = json.loads(path.read_text())
        if saved["cache_key"] != key: raise RuntimeError("Rest cache inputs changed")
        return saved
    seed = original_seed()
    scales = seed[KEEP]
    fixed = model.parameters.geometry.fixed_cell_anion_equivalents_fmol

    def expand(z):
        y = np.zeros(13)
        y[KEEP] = scales * np.exp(z)
        y[4] = y[0] + y[1] - y[2] - fixed
        y[10] = y[6] + y[7] - y[8]
        return y

    def fun(z):
        return model.rhs(0., expand(z), genotype=genotype)[KEEP] / scales * 1000.

    fit = root(fun, np.zeros(10), method="hybr", options={"xtol": 1e-10, "maxfev": 1500})
    y = expand(fit.x)
    scaled_residual = float(np.max(np.abs(fun(fit.x))))
    if scaled_residual > 1e-7: raise RuntimeError(f"Rest not accepted: {genotype_name}, {fit.message}, {scaled_residual}")
    e = evaluate(model, 0, y, genotype)
    def log_rhs(z):
        v = expand(z)
        return model.rhs(0., v, genotype=genotype)[KEEP] / v[KEEP]
    jac = np.column_stack([(log_rhs(fit.x + np.eye(10)[k]*1e-5) - log_rhs(fit.x - np.eye(10)[k]*1e-5))/2e-5 for k in range(10)])
    eigen = np.linalg.eigvals(jac)
    result = dict(genotype=genotype_name, cache_key=key, shared_parameter_hash=model.task48_payload_hash,
        state=y.tolist(), state_names=list(model.state_names), seed_status="saved WT state is numerical seed only",
        nfev=int(fit.nfev), solver_success=bool(fit.success), solver_message=str(fit.message),
        scaled_rhs_max=scaled_residual, scaling="1000 s times RHS divided by seed amount or volume",
        full_rhs=e.rhs.tolist(), positive_core=bool(np.all(y[:12] > 0)),
        charge_residuals_fmol=dict(e.diagnostics.state_charge_fmol),
        core_eigenvalues_per_s=[[float(v.real), float(v.imag)] for v in eigen],
        regulatory_eigenvalue_per_s=-1./model.regulatory_model.ae4_regulatory_model.tau_activation_s,
        locally_stable=bool(np.max(eigen.real) < 0),
        observables=observable(model, 0, y, genotype), bath_projection=model.task48_bath_projection)
    if not result["positive_core"] or not result["locally_stable"]: raise RuntimeError("Rest positivity/stability gate failed")
    dump(path, result)
    print(json.dumps({"rest":genotype_name,"cl_mM":result["observables"]["cl_i_mM"],"ph":result["observables"]["ph_i"],"stable":True}), flush=True)
    return result


def structural_contrast():
    """Exact algebra plus independent runtime checks, not a restricted fit."""
    cch, both = build_model("CCH_ONLY"), build_model("CCH_IPR")
    y = original_seed(); y0 = y.copy(); y1 = y.copy()
    y0[-1] = 0.; y1[-1] = 1.
    g = GENOTYPES["AE4_KO"]
    difference = cch.rhs(1.,y0,genotype=g)[:12] - both.rhs(1.,y1,genotype=g)[:12]
    if np.max(np.abs(difference)) != 0: raise AssertionError("KO invariance failed")
    mean1, se1, mean2, se2 = 2.30, .10, .90, .09
    common = (mean1/se1**2 + mean2/se2**2)/(1/se1**2 + 1/se2**2)
    result=dict(
        hypothesis="Shared expected AE4 KO uptake for CCh and CCh plus IPR under a common source faithful observation rule",
        source="Peña Münzenmayer 2015 Table 1, AE4 KO CCh n=6 and combined n=9",
        predicted_contrast=0., experimental_contrast=mean1-mean2, units="1e-3 s^-1",
        exact_contrast_parameter_jacobian="zero for every admissible shared parameter, including AE4 gain and time constant",
        runtime_core_rhs_max_difference=float(np.max(np.abs(difference))),
        independent_mean_error_contrast_se=float(np.hypot(se1,se2)),
        independent_contrast_standard_errors=float((mean1-mean2)/np.hypot(se1,se2)),
        minimum_two_row_chi_square=float((mean1-mean2)**2/(se1**2+se2**2)),
        relaxed_optimal_common_mean=common,
        relaxed_two_row_standardised_residuals=[(common-mean1)/se1,(common-mean2)/se2],
        relaxed_mean_is_full_model_prediction=False,
        maximum_contrast_se_any_covariance=se1+se2,
        minimum_standard_errors_any_covariance=(mean1-mean2)/(se1+se2),
        additional_ipr_only_ko_prediction=0., ipr_only_observation=.20, ipr_only_se=.03,
        optical_qualification="Matched normalisation and recovery-window rule. Arbitrary cohort optical gains or differing windows are not introduced.",
        likelihood_qualification="Reported SE Gaussian working mean model; no global calibrated p value and no completed full joint likelihood.")
    dump(OUT/"structural_contrast.json", result)
    return result


def trajectory(genotype_name, arm):
    genotype = GENOTYPES[genotype_name]
    rest = resting_state(genotype_name)
    model = build_model(arm)
    name = genotype_name + "__" + arm
    path = OUT/"protocol_predictions"/(name+".json")
    if path.exists():
        saved=json.loads(path.read_text())
        if saved["shared_parameter_hash"] != model.task48_payload_hash or saved["rest_cache_key"] != rest["cache_key"]:
            raise RuntimeError("Trajectory cache mismatch")
        return saved
    seed = np.asarray(rest["state"])
    def rhs(t, value):
        e = model.evaluate(max(float(t),1e-9),value[:13],genotype=genotype)
        d=e.diagnostics
        q=d.water.lumen_outflow_pL_s
        return np.r_[e.rhs,q,q*d.observables.lumen_acid_base.hco3_mM]
    times=np.arange(0.,601.,10.)
    sol=solve_ivp(rhs,(0.,600.),np.r_[seed,0.,0.],method="BDF",rtol=1e-8,atol=1e-10,max_step=2.,t_eval=times)
    if not sol.success: raise RuntimeError(sol.message)
    rows=[]
    for i,t in enumerate(sol.t):
        row=observable(model,float(t),sol.y[:13,i],genotype)
        row.update({key:float(value) for key,value in zip(model.state_names,sol.y[:13,i])})
        row.update(water_integral_pL=float(sol.y[13,i]),bicarbonate_integral_fmol=float(sol.y[14,i]))
        rows.append(row)
    write_csv(path.with_suffix(".csv"),rows)
    amount=sol.y[13]; bicarbonate=sol.y[14]
    get=lambda t:float(amount[np.where(times==t)[0][0]])
    early_pH=rows[1]["ph_i"]-rows[0]["ph_i"]
    result=dict(genotype=genotype_name,protocol=arm,shared_parameter_hash=model.task48_payload_hash,
        rest_cache_key=rest["cache_key"],initial_state=seed.tolist(),final_state=sol.y[:13,-1].tolist(),
        solver=dict(method="BDF",rtol=1e-8,atol=1e-10,max_step_s=2.,nfev=int(sol.nfev)),
        duration_s=600, water_pL=float(amount[-1]),bicarbonate_fmol=float(bicarbonate[-1]),
        water_windows_pL={"0_120":get(120),"120_180":get(180)-get(120),"180_600":get(600)-get(180)},
        ph_change_first_10s=early_pH,ph_change_first_10s_status="qualitative numerical diagnostic, not BCECF amplitude or source slope window",
        maximum_charge_fmol=max(row["max_charge_fmol"] for row in rows),
        maximum_current_residual_A=max(row["max_current_residual_A"] for row in rows),
        positive_core=bool(np.all(sol.y[:12]>0)),
        ph_range=[min(row["ph_i"] for row in rows),max(row["ph_i"] for row in rows)],
        cl_range_mM=[min(row["cl_i_mM"] for row in rows),max(row["cl_i_mM"] for row in rows)],
        bath_projection=model.task48_bath_projection,
        prediction_status="unfitted conditional inherited model diagnostic")
    dump(path,result)
    print(json.dumps({"trajectory":name,"water_pL":result["water_pL"],"positive":result["positive_core"]}),flush=True)
    return result


def main():
    structural_contrast()
    for name in ["WT","AE4_KO","AE2_KO"]: resting_state(name)
    if "--rests-only" in sys.argv: return
    for name in ["WT","AE4_KO","AE2_KO"]:
        for arm in ["CCH_ONLY","IPR_ONLY","CCH_IPR","NHE_EIPA"]: trajectory(name,arm)


if __name__ == "__main__": main()
