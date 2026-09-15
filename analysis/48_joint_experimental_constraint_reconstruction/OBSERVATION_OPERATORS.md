# Task 48 observation operator advice

Immutable input: `5c945c4a7857657b5f7c76d5bf9e3d3abc7d5721`.
Advisory worker only: no model evaluations, fits, canonical edits or publication.

## Main finding

The primary study uses SPQ **F0/F**, where F0 is fluorescence at the start of the experiment. It does not use F/F0. Table 1 exit is a finite fluorescence change, whereas uptake is a slope from a linear regression over the indicated recovery segment. Neither is an individual transporter flux. The primary Table 1 stimulation protocol is physiological bath followed by agonist, with initial Cl loss and recovery. Low Cl depletion followed by high Cl readdition belongs to the isolation assays and must not be substituted for Table 1.

Primary source: Peña-Münzenmayer et al. (2015), DOI 10.1074/jbc.M114.612895, [PMC4409235](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4409235/), methods, Table 1, Figures 1–5. Local source inspected: sibling `primary_data/jbc.html`. The source gives calibration concentrations 20, 50 and 80 mM, but no fluorescence readings, fitted quenching constant, or experimentally resolved calibration function. The exact per-experiment fitting windows and cross-observable covariances are not provided numerically. Resting Cl and pH are already calibrated measurements.

## Operators and mapping status

| Measurement | Correct model observable | Status and what is retained |
|---|---|---|
| Resting intracellular Cl | nCl_i / V_i at the genotype's own rest | Direct; compare each control and knockout mean and SE in mM |
| Resting pH | Model acid-base speciation pH at the genotype's own rest | Direct; compare each mean and SE; do not convert to H concentration before applying pH SE |
| SPQ initial exit | Z(t_min) − Z(t_pre), Z=F0/F | Conditional; signed change is negative. A Cl concentration change is a mechanistic diagnostic until calibration is fixed |
| SPQ initial uptake | OLS slope of Z(t) in the source recovery window | Conditional; use the same temporal window and net intracellular concentration trajectory |
| Isolated transporter uptake | OLS slope of Z after source depletion, inhibitor equilibration and readdition | Conditional on representable protocol and optical map; do not substitute isolated transporter J |
| NHE assay alkalinisation | BCECF excitation-ratio R(pH(t))/R(pH0) and its source-window slope | Exact optical amplitude blocked without BCECF calibration. Direction of pH change and EIPA reversal are qualitative mappings |
| Total saliva at 600 s | V_g=integral_0^600 q_out,g(t) dt; compare V_KO/V_control | Quantitatively mappable as a ratio subject to common effective acinar-to-gland transfer and no genotype-dependent duct water contribution |
| Secretion time course | Ratios of source-matched collection-bin integrals of q_out | Qualitative currently; exact means/SE require digitisation. Compare genotype to its own strain control |
| Salivary HCO3 amount | A_g=integral_0^600 q_out,g(t) [HCO3]_l,g(t) dt; compare A_KO/A_control | Conditional gland mapping; missing ducts may alter composition. Do not use luminal HCO3 concentration alone, total carbon, or total alkalinity |
| Reported percentage reductions | 100(1−O_KO/O_control), with O the correct assay observable | Conditional on O mapping; reported SE applies to the percentage, not individual means |

`q_out` is `diagnostics.water.lumen_outflow_pL_s`. Model HCO3 is `diagnostics.observables.lumen_acid_base.hco3_mM`. Their product is fmol/s. Integrating gives fmol; ratios remove unit and common gland scale. Model `outflow_sources_fmol_s['tic']` is total inorganic carbon, and `['alkalinity']` includes acid-base charge equivalents; neither equals a bicarbonate assay. The paper's methods say concentration was divided by saliva amount when calculating microequivalents, which is dimensionally inconsistent. Preserve the reported bicarbonate amounts without claiming to reproduce that wording as a physical operator; the concentration-times-volume amount relation is the dimensional definition.

## Net concentration, amount and fluorescence are distinct

The exact concentration derivative for model cell volume V is

    C = nCl / V
    dC/dt = (dnCl/dt − C dV/dt) / V.

Here dnCl/dt contains all represented chloride sources and sinks. The water correction cannot be dropped merely because two main transport pathways have been inhibited. `J_NKCC` alone may count cotransporter cycles rather than chloride molecules, with two Cl per cycle. Return the model's total Cl amount source before converting to concentration.

Under the explicitly conditional ideal quenching law F(C)=A/(1+K C),

    Z(t) = F(C0)/F(C(t)) = (1+K C(t))/(1+K C0),
    ΔZ = K ΔC/(1+K C0),
    slope(Z) = K slope(C)/(1+K C0).

This is a proposed measurement model, not a reconstructed source calibration. Source calibration concentrations alone do not determine K. F0 is starting fluorescence, not unquenched fluorescence. Positive Cl uptake produces positive Z slope, and Cl exit produces negative ΔZ, agreeing with the signs in Table 1.

For genotype ratio R = slope(Z_KO)/slope(Z_WT),

    R = slope(C_KO)/slope(C_WT) × (1+K C0_WT)/(1+K C0_KO).

Consequently, a simple ratio of concentration slopes is still conditional and is generally wrong if the starting concentrations differ. The AE4 resting concentrations demonstrably differ. Readdition does not guarantee identical depleted intracellular states or identical optical normalisation baselines.

Two permitted reduced comparisons are possible, with explicit information loss:

1. For a single trace, slope(Z)/abs(ΔZ) = slope(C)/abs(ΔC) when the same affine optical law applies to its two segments. This retains the recovery rate relative to initial loss, and loses the absolute exit amplitude and absolute uptake scale. Table 1 only supplies group means; a ratio of group means inherits a common-within-group-calibration assumption and is not the mean individual recovery ratio.
2. For two stimulus conditions with genuinely equal C0 and the same optical calibration, their Z slope ratio equals the corresponding concentration slope ratio. This retains relative stimulus response and loses absolute transport scale. Genotype comparisons must not silently borrow this same-C0 assumption.

Neither reduction can be advertised as satisfying all original independent quantitative constraints. Unmapped source values must remain in residual tables with `prediction=null`, `standardised_residual=null` and an explicit blocker.

## Timing and definition of initial uptake

For Table 1, the measurement is recovery following initial stimulus-induced exit. Instantaneous t=0+ chloride derivative is normally the exit phase and is not its initial uptake slope. Use an explicit source-matched interval following the minimum and record t_start, t_stop and regression convention. The source has drawn regression segments but does not supply their coordinates numerically. If the intervals are digitised, preserve the figure, axis mapping, points, script and digitisation uncertainty.

For isolation assays, determine the readdition time from the assay and fit the initial rising Z segment. Do not replace this with a fixed ad hoc 0–10 s or 0–30 s interval and call it measured. A preregistered algorithmic window may be presented as an assumption and should be tested for sensitivity, but cannot remove the missing source window.

For saliva, keep 0–120 s, 120–180 s and 180–600 s integrals as clearly labelled summary diagnostics. The reported initial 2–3 minute comparability and later deficit do not supply numeric error bands. A model in which the full deficit appears immediately fails the qualitative timing fingerprint. A formal acceptance tolerance requires the source figure errors or individual observations. The different AE4 versus AE2 control strain kinetics must not be fitted using an invented strain multiplier; compare within strain and report any limitation of one shared model trajectory.

## Uncertainty handling

All quoted 2015 means are means ± SE; n is number of experiments, with preparations from at least three mice per condition. SE directly quantifies uncertainty in the estimated population mean under the recommended Gaussian mean likelihood. Do not multiply or divide SE by sqrt(n) unless explicitly constructing a separate biological variance model; no individual-level covariance or mouse-clustering correction can be recovered from these summaries.

For y=mean_KO, x=mean_control, ratio r=y/x and covariance c between the means, the first-order ratio SE is

    SE(r)^2 = SE(y)^2/x^2 + y^2 SE(x)^2/x^4 − 2y c/x^3.

Independent genotype groups motivate c=0 for the simple gland ratio comparison, as an explicit assumption. Paired slope and exit measurements from the same Table 1 traces are correlated and their covariance is unknown. Report the chosen working independence assumption and do not give a formally calibrated global chi-square p-value. Multiple ratios sharing a control are correlated, even when original groups are independent:

    Cov(y_i/x, y_j/x) ≈ y_i y_j Var(x)/x^4.

Use one independent data representation per assay. Do not score both raw means and their derived ratio as independent observations. The published isolated exchanger reductions share a pooled control; their covariance cannot be recovered from the percentage SEs alone. A diagonal working objective is possible but must be labelled as such, with no claim of exact joint likelihood calibration.

For gland ratio prediction r_model, a Gaussian scale-profile residual avoids treating the denominator mean as exact:

    z = (mean_KO − r_model mean_control)
        / sqrt(SE_KO^2 + r_model^2 SE_control^2 − 2r_model c).

This is the signed profile least-squares contrast for a shared unknown scale (interior scale solution). It is preferable to applying an observed-ratio delta SE when the denominator is noisy. For final reporting give both original group means and the observed ratio with labelled approximate SE.

Derived independent-group ratio summaries from the seed, pending ledger verification:

| Assay | KO/control or stimulated/basal | Approximate ratio SE |
|---|---:|---:|
| AE4 HCO3 amount | 1.339623 | 0.384711 |
| AE2 HCO3 amount | 1.575000 | 0.476767 |
| Isolated exchanger control IPR/no IPR | 2.775641 | 0.674212 |
| Isolated exchanger AE4 KO IPR/no IPR | 1.136842 | 0.152628 |
| AE2 KO Na free/Na present | 0.377778 | 0.080306 |

These ratios retain dimensionless contrasts and discard absolute assay scale. For SPQ rows they remain **observed fluorescence ratios**, not automatically concentration-slope or carrier-flux ratios. Na dependence may involve paired experiments (both n=7); independence is a working assumption until pairing/covariance is established.

## Qualitative observations and numerical failure

A non-significant observation does not imply exact equality or define a numerical acceptable band. Keep NHE preservation, isolated NKCC similarity, AE2 secretion neutrality and initial exit preservation as qualitative checks unless means/SE have been obtained. Inhibiting the intended pathway to zero is a declared ideal inhibitor intervention, not a model claim that the measured residual response must vanish. The 4.4% NKCC assay residual, and 12.3% double exchanger KO residual, may include other sources or imperfect inhibition; do not add a generic residual current to make them fit.

If an assay cannot be represented (for example CA inhibition under instantaneous acid-base equilibrium, exact Na removal outside transporter domain, or experimental perfusion of isolated acini outside the lumen model), mark it blocked. Failure to construct a valid observation map is not proof that all possible shared parameter vectors are scientifically infeasible. A numerical failure for one parameter vector is not a general impossibility theorem.

## Inherited regulation audit affecting interpretation

`analysis/47_dynamic_potassium_recycling_reconstruction/common.py:load_model` already constructs `R1EffectiveActivation` inside `CompositeAe4NkccRegulation`, and passes it to `ModernFullModel`. `camp_pka.py` explicitly labels kinetic defaults as unmeasured sensitivity values. Thus the presence of AE4 cAMP/PKA regulation must be acknowledged. A gain re-identification is a parameter update, not discovery or insertion of an absent mechanism. Isolated endpoint activation data do not identify an acute regulatory time constant.

## Conditional architectural contradiction independent of calibration magnitude

If the audited AE4 KO equations give identical full trajectories for CCh only and CCh plus IPR, then a common deterministic observation operator must also give identical expected uptake. The unknown SPQ calibration magnitude does not rescue this invariance: both model trajectories share their genotype resting state, starting fluorescence normalisation and recovery-segment selection rule. The experimental uptake means differ by 2.30−0.90=1.40 in units of 10^-3 s^-1, with SE 0.10 and 0.09.

Under independent Gaussian mean errors, minimising the two residual squares over ANY common predicted uptake value yields

    min_m [(m−2.30)/0.10]^2 + [(m−0.90)/0.09]^2
      = (1.40)^2/(0.10^2+0.09^2)
      = 108.28729281767956.

This is a lower bound for this two-observation objective alone, not a computed joint optimum. With unknown covariance, SE of the difference cannot exceed 0.10+0.09=0.19. The discrepancy remains at least 7.368421 reported SE of the difference, or 54.293629 as a squared contrast. These bounds use reported SE estimates and are not finite-sample calibrated p-values.

This argument is conditional on a common source-faithful optical map and initial-slope extraction rule. Unevidenced systematic changes in optical calibration between cohorts, or deliberately different regression intervals applied to one identical nonlinear curve, could mathematically break measurement invariance. The paper provides no basis for adding such protocol-specific observation adjustments. State this condition explicitly, preserve missing calibration as a limitation, and do not pretend the whole undefined likelihood has been optimised. A correction acting only on AE4 cannot remove an exact invariant already present when AE4 is absent.


Orchestrator integration: reviewed and accepted for the 48A evidence freeze. This document supplies evidence and restrictions, not a completed fit.
