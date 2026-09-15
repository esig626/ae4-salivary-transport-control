"""Construct the complete ledger without inventing unavailable predictions."""
import json
from pathlib import Path
import numpy as np
from run_baseline import HERE, OUT, dump, write_csv
from protocol_layer import shared_payload, payload_hash, gaussian_ratio_contrast

data=json.loads((HERE/"constraints.json").read_text())
rests={g:json.loads((OUT/"rests"/(g+".json")).read_text()) for g in ["WT","AE4_KO","AE2_KO"]}
pred={g:json.loads((OUT/"protocol_predictions"/(g+"__CCH_IPR.json")).read_text()) for g in rests}
rows=[]


def add(identifier, mean=None, se=None, n=None, *, genotype="", protocol="", units="", operator="", prediction=None, status="unavailable", reason="", residual=None, scored=False):
    if residual is None and prediction is not None and mean is not None and se:
        residual=(prediction-mean)/se
    rows.append(dict(constraint_id=identifier,experimental_mean=mean,SE=se,n=n,units=units,
        genotype=genotype,protocol=protocol,observation_operator=operator,prediction=prediction,
        standardised_residual=residual,status=status,reason=reason,scored=scored))


for family,groups in [("ae4_table1",["WT","AE4_KO"]),("ae2_table1",["control","AE2_KO"])]:
    for item in data[family]:
        obs=item["observable"]
        for group in groups:
            g="WT" if group=="control" else group
            mean,se,n=item[group]
            common=dict(genotype=group,protocol=obs,units=item.get("units","pH"))
            if obs in ["rest_Cl_i","rest_pH_i"]:
                key="cl_i_mM" if obs=="rest_Cl_i" else "ph_i"
                add(family+"/"+obs+"/"+group,mean,se,n,**common,
                    operator="genotype specific stationary concentration" if key=="cl_i_mM" else "stationary carbonate pH",
                    prediction=rests[g]["observables"][key],status="conditional_B_plus_projection",scored=True,
                    reason="one shared WT trajectory represents both strain controls; source bath projection qualification")
            elif obs=="initial_Cl_uptake_IPR" and group=="AE4_KO":
                add(family+"/"+obs+"/"+group,mean,se,n,**common,operator="slope of any common optical map on an exactly constant core state",
                    prediction=0.,status="analytic_invariance",scored=True,
                    reason="IPR has no core effect after AE4 deletion; stable genotype rest, no calibration magnitude needed")
            else:
                add(family+"/"+obs+"/"+group,mean,se,n,**common,operator="F0/F amplitude or source recovery regression slope",
                    reason="individual fluorescence calibration and source fitting windows unavailable; do not equate flux and signal")

w=pred["WT"]
deficit=100*(1-pred["AE4_KO"]["water_pL"]/w["water_pL"])
add("whole_gland/AE4_10min_deficit",35.,4.7,"6 WT; 6 KO",genotype="AE4_KO",protocol="CCH_IPR_600s",units="percent",
    operator="100*(1-integral(q_KO)/integral(q_WT))",prediction=deficit,status="conditional_gland_ratio",scored=True,
    reason="shared acinus to gland water scale; effective B+ bath; no fit")
for family,g,control in [("AE4","AE4_KO","WT"),("AE2","AE2_KO","control")]:
    item=data["salivary_HCO3_output"][family]
    for group in [control,g]:
        mean,se,n=item[group]
        add("HCO3/"+family+"/"+group,mean,se,n,genotype=group,protocol="CCH_IPR_600s",units="microequiv",
            operator="whole gland bicarbonate amount",status="retained_via_ratio",reason="no independently known acinus to gland scale; do not score twice")
    c,cs,_=item[control];k,ks,_=item[g]
    ratio=pred[g]["bicarbonate_fmol"]/w["bicarbonate_fmol"]
    observed_ratio=k/c
    approx_se=float(np.hypot(ks/c,k*cs/c**2))
    add("HCO3/"+family+"/ratio",observed_ratio,approx_se,genotype=g,protocol="CCH_IPR_600s",units="ratio",
        operator="ratio of integral(q*[HCO3]_l); Gaussian scale profile residual",
        prediction=ratio,status="conditional_gland_ratio",scored=True,
        residual=-gaussian_ratio_contrast(ratio,c,cs,k,ks),
        reason="independent group mean errors; SE column uses ratio delta approximation, residual profiles the unknown common scale; duct transfer unmodelled")

item=data["isolated_NKCC1"]["bumetanide_50uM_reduction_percent"]
add("isolated_NKCC1/bumetanide",item["mean"],item["SE"],item["n"],protocol="NKCC_ISOLATION_BUMETANIDE",units="percent",
    operator="SPQ uptake reduction",reason="CA inhibition and absent bicarbonate not representable; NKCC bath dependence absent")
ex=data["isolated_HCO3_dependent_exchanger"]
for g,item in ex["reduction_percent"].items():
    add("isolated_exchanger/reduction/"+g,item["mean"],item["SE"],item.get("n_KO"),genotype=g,protocol="EXCHANGER_ISOLATION",units="percent",
        operator="normalised SPQ slope reduction",reason="conditioning, optical map and exact slope windows unavailable")
for group in ["control_no_IPR","control_plus_IPR","AE4_KO_no_IPR","AE4_KO_plus_IPR"]:
    mean,se,n=ex["IPR_activation"][group]
    add("isolated_exchanger/IPR/"+group,mean,se,n,genotype=group,protocol="EXCHANGER_IPR",units="1e-3 s^-1",
        operator="normalised SPQ uptake",reason="isolated assay conditioning and exposure time unavailable; existing gain not fitted from saliva")
for group in ["Na_present","Na_free"]:
    mean,se,n=ex["AE4_Na_dependence_in_AE2_KO"][group]
    add("isolated_exchanger/"+group,mean,se,n,genotype="AE2_KO",protocol="EXCHANGER_NA_FREE",units="1e-3 s^-1",
        operator="normalised SPQ uptake",reason="exact sodium absent condition outside constitutive domain; no numerical epsilon substituted")
for item in data["additional_verified_primary_constraints"]:
    add("additional/"+item["observable"],item.get("mean"),item.get("SE"),protocol="source specified assay",units=item.get("units","qualitative"),
        operator="source fluorescence observation",reason=item.get("note",item.get("result",""))+"; no fabricated numerical mapping")

dump(OUT/"baseline_residuals.json",rows)
write_csv(OUT/"baseline_residuals.csv",rows)
contrast=json.loads((OUT/"structural_contrast.json").read_text())
qualitative={
    "AE4_timing":{name:pred["AE4_KO"]["water_windows_pL"][name]/w["water_windows_pL"][name] for name in w["water_windows_pL"]},
    "AE2_water_ratio":pred["AE2_KO"]["water_pL"]/w["water_pL"],
    "NKCC_isolation":"unavailable; no physiological compensation claim inferred from a different assay",
    "NHE":{g:{"combined_first10s_pH_change":pred[g]["ph_change_first_10s"],
        "EIPA_first10s_pH_change":json.loads((OUT/"protocol_predictions"/(g+"__NHE_EIPA.json")).read_text())["ph_change_first_10s"]} for g in rests},
    "NHE_qualification":"positive pH direction is insufficient: KO alkalinisation is much smaller and simultaneous EIPA does not reverse it; source inhibitor preincubation/window and BCECF calibration not identified, so this is a conditional diagnostic",
    "initial_Cl_exit":"fluorescence amplitude unavailable without calibration; no false quantitative pass",
    "bicarbonate":"conditional gland ratio retained; no claim about missing duct bicarbonate handling",
}
dump(OUT/"qualitative_checks.json",qualitative)
checks=dict(milestone="48C",accepted_rest_roots=3,locally_stable_roots=sum(x["locally_stable"] for x in rests.values()),
    completed_trajectories=12,parameter_optimisations=0,new_equation_mechanisms=0,
    baseline_residual_rows=len(rows),mapped_scored_rows=sum(x["scored"] for x in rows),
    structural_core_RHS_identity=contrast["runtime_core_rhs_max_difference"],
    minimum_two_row_chi_square=contrast["minimum_two_row_chi_square"],
    decision="inherited architecture rejected conditional on common source faithful observation rule; full numeric likelihood unavailable",
    max_rest_scaled_rhs=max(x["scaled_rhs_max"] for x in rests.values()),
    max_rest_charge_fmol=max(abs(v) for x in rests.values() for v in x["charge_residuals_fmol"].values()))
dump(OUT/"baseline_checks.json",checks);dump(OUT/"latest_checks.json",checks)

lines=["# Inherited architecture residuals", "", "The inherited architecture is not a defensible joint reconstruction under a common source faithful observation rule. This decision has a parameter independent component; it is not based on failure of a restricted optimiser.", "", "## Structural gate", "", f"AE4 KO CCh and CCh plus IPR have identical core dynamics for every shared admissible parameter vector. Beta input acts only on the deleted AE4 contribution. The observed uptake means are 2.30 ± 0.10 and 0.90 ± 0.09 (10⁻³ s⁻¹). The minimum two row Gaussian mean objective over any common prediction is {contrast['minimum_two_row_chi_square']:.8f}. This is a relaxed lower bound, not an optimised full model residual. Its contrast Jacobian is identically zero for all shared parameter directions. Existing AE4 activation gain and timing cannot change it.", "", "This statement assumes the same source faithful optical normalisation and recovery extraction rule. Exact per experiment windows and calibrations were not published. Arbitrary cohort optical gains or selectively chosen windows would weaken inference, but have no evidential basis and are not added here. No calibrated global p value is claimed.", "", "## Genotype resting states", "", "| Genotype | Rest Cl, mM | Rest pH | Local stability |", "| --- | ---: | ---: | --- |"]
for g,r in rests.items():lines.append(f"| {g} | {r['observables']['cl_i_mM']:.8f} | {r['observables']['ph_i']:.8f} | stable |")
lines += ["", "Each root uses its actual deletion, exact cell and lumen charge constraints and one shared inherited vector. The saved WT state was a numerical seed only. Local stability and residual checks do not establish global uniqueness or chronic biological adaptation.", "", "## Supported whole gland diagnostics", "", f"Under the source B+ projection, AE4 KO has a cumulative deficit of {deficit:.8f}% (a small increase in secretion), compared with the measured 35 ± 4.7%. AE2 KO/control water ratio is {qualitative['AE2_water_ratio']:.8f}. These are unfitted predictions from the genotype resting states, not Task 47 acute deletion trajectories.", "", "| Window, s | AE4 KO/WT water ratio |", "| --- | ---: |"]
for k,v in qualitative["AE4_timing"].items():lines.append(f"| {k.replace('_',' to ')} | {v:.8f} |")
lines += ["", "## Full data accounting and limits", "", f"The residual table contains {len(rows)} source and derived rows, with {sum(x['scored'] for x in rows)} quantitatively mapped diagnostic rows. Unavailable predictions and residuals are null in JSON and blank in CSV. Raw gland bicarbonate amounts are retained but only their ratios are scored. No target is silently dropped or assigned a fabricated fluorescence flux conversion.", "", "The independent group Gaussian interpretation uses SE as uncertainty of the reported mean. Correlations and individual animal clustering are unreported. If arbitrary covariance between the two contradictory means is allowed, their contrast SE is at most 0.19 and the discrepancy is still at least 7.3684 reported standard errors. This is not a finite sample p value.", "", "A full numerical observation Jacobian and joint optimum are not defined because essential protocols and observation maps are missing. The exact zero sensitivity of the contradictory contrast is an early rejection gate. The 33 provenance eligible parameters were not treated as an identifiable fit subspace. No optimisation or parameter correction was performed.", "", "NHE and EIPA pH diagnostics are retained, with their unverified timing and optical qualification. Isolated NKCC, isolated exchanger conditioning, sodium absent assays and exact SPQ values remain unavailable. The model therefore cannot claim to satisfy the isolated NKCC constraint or quantify compensation under it.", "", "All numerical trajectories use the qualified positive B+ reservoir projection declared at 48A: exact transported ion concentrations, fixed bath pH, 17 mM spectator osmoles and explicit divalent/carbonate charge bookkeeping. The omitted HEPES/titrant chemistry and duct/acinar assay topology remain limitations. No inherited intracellular equation was changed.", "", "Reproduce only in a separate review workspace using `run_baseline.py` then `summarise_baseline.py`. Existing per genotype and protocol records are reused after matching their input identities. No Tasks 46 or 47 calculation is rerun."]
(HERE/"BASELINE_RESIDUALS.md").write_text("\n".join(lines)+"\n")
print(json.dumps(checks))
