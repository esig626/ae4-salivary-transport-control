# 51A source and architecture freeze

## Scope and provenance

Operational start is the user's exact remote head `ad0e7344c8ba79552cbf54714b42e7131bcc1142`, descending from scientific parent `4bb2c89870fe3f4d4bb887cf9d40ba7de377a725`. The complete 2,151-file tree was reconstructed from authenticated Git objects and checked against tree `6889b99d839ce3b68808d0a17d532945dd3b2e83`. Task 50 archive independently remains `b16c30094b95f61a73d8f1977cd58e79c7bb50f6`. Five advisers read sources/code only; the orchestrator accepts decisions and performs all computation/writes/publication.

The mandatory ledger, phenotype convention, benchmark, evidence handoff, novelty/decision records, Tasks 49/48/39/40/41 reports, active source and prompt were read. RESULTS_LEDGER and analysis index were also consulted. Scoped text searches of prompts, analyses, results and source identify the earlier Task13B gate and Task49 conserved implementation. Task49 SOURCE_SEED already mentioned beta-NKCC literature: novelty is the authorised upstream equation/test, not first discovery of those papers. All 55 remote branch names were inspected; no claim is made to have freshly searched all historical branch contents. The prior exclusions remain binding.

## Independently checked source anchors

| Primary source | Verified information | Role and qualification |
|---|---|---|
| [Chaib et al. 1998](https://pubmed.ncbi.nlm.nih.gov/9880083/), DOI 10.1016/S0196-9781(98)00134-X | Rat submandibular IPR increases NKCC-mediated NH4 influx approximately 2.5-fold; bumetanide blocks; forskolin mimics; H8 inhibits | Prescribed central cross-species activity assumption; NH4/BCECF assay, not measured mouse Cl-loading gain |
| [Paulais and Turner 1992](https://www.jci.org/articles/view/115695), DOI 10.1172/JCI115695 | Rat parotid recovery approximately threefold; IPR K1/2 21.5 nM; beta1/cAMP/kinase evidence | Prescribed 3-fold sensitivity; correct first author is Paulais |
| [Tanimura et al. 1995](https://pubmed.ncbi.nlm.nih.gov/7559664/) | IPR-regulated phosphorylation; half-maximal approximately 20 nM; approximately sixfold activity in AlF4/IPR comparison context | Mechanistic support only; sixfold is not a rescue parameter |
| [Kurihara et al. 1999](https://pubmed.ncbi.nlm.nih.gov/10600770/) | Increased activity/phosphorylation and high-affinity bumetanide-binding sites, unchanged inhibitor sensitivity | Consistent with recruitment/activity; frozen concentration-response core retained |
| [Kurihara et al. 2002](https://pubmed.ncbi.nlm.nih.gov/11880270/) | cAMP mimics regulatory phosphorylation; PKA inhibition blocks | PKA participation, not proven simple direct purified-PKA phosphorylation |
| [Catalan et al. 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4343136/) | TMEM16A-independent IPR secretion, CFTR/ClC-2 controls, blocker sensitivity, swelling, outwardly rectifying Cl current | Independent mouse evidence; does not establish molecular LRRC8 identity or all-apical patch current |
| [Pena-Munzenmayer et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8887885/) | Beta/PKA activation of AE4, S173 dependence | Existing inherited AE4 arm retained once; no additional signalling state |

Quantitative rat gains are SOURCE FACTS in their assays and ASSUMPTION/SENSITIVITY when transferred to mouse. They are not mouse confidence bounds. Acute times and dose-response constants are not added to the model.

## Catalan supplement retrieval and direct-map decision

The main article was retrieved (HTTP 200; 192,871-byte HTML). The official [PMC supplement](https://pmc.ncbi.nlm.nih.gov/articles/instance/4343136/bin/pnas.201415739SI.pdf) returned HTTP 200 but 1,816-byte HTML download challenge, not PDF. The [PNAS supplement](https://www.pnas.org/doi/suppl/10.1073/pnas.1415739112/suppl_file/pnas.201415739si.pdf) and Europe PMC direct supplement returned 403. Europe PMC supplementaryFiles returned an XML non-open-access error; fullTextXML and PMC OA routes returned 404. No challenge was bypassed. **SI was not retrieved or inspected.** This is an access limitation, not a claim that SI lacks information.

The accessible primary article identifies nystatin-perforated whole-cell Figure5 recordings, quasi-physiological bicarbonate external solution, 5 uM IPR, 50 uM DCPIB, and -60 to +100 mV steps. The -23.6 mV ECl stated for Figure2 is a different CaCC protocol and cannot be reused for Figure5.

A direct map would require `g_V = I_apical / [beta*s_patch*(V_clamp-ECl_patch)]` with consistent signed current, known patch swelling, apical attribution and recorded/model-cell scale. These quantities are not jointly identified; rectification further qualifies a constant-ohmic conversion. No fitted nuisance factor is added. Direct mapping is **UNRESOLVED**, and the predeclared fallback is therefore used.

Reuse Task49 metrology unchanged: IPR current magnitude 208.59375 +/-35.15625 pA; DCPIB-blocked magnitude 107.8125 +/-17.1875 pA at -60 mV. Narrative swelling is 12.5+/-0.2%; graphical terminal proxy is 1.160311284 near450.92 s, raster envelope1.14786–1.17276 (not SEM). Preserve the discrepancy. Whole-gland outputs remain qualitative/ratio checks because model-to-gland scale is unidentified. Chronic TMEM16A deletion also alters NKCC1/AQP5; an acute mask is only an idealisation.

## Exact algebra and implementation contract

`M_NKCC = M_Ca + beta*(M_beta-1)`.

At beta zero the entire parent regulator is unchanged. At resting calcium and beta one it equals M_beta. At full inherited CCh it equals3.25 for central2.5 and3.75 for sensitivity3. The source-fixed Palk cycle is linear in this multiplier; reversal and 1Na:1K:2Cl remain exact. Implement a delegating adapter around CompositeAe4NkccRegulation, changing only its propagated NKCC multiplier. Both inherited NKCC and NBC wrappers consume that diagnostic, so active Ca cannot overwrite the beta increment. Retain the existing AE4 state/law; add no state.

Wrap the result in the unchanged Task49 VracModel: `G_V=g_V*beta*max(V/Vref-1,0)`. It contributes to both Cl amount balances and full current closure. The builder uses only Task48/40 parent machinery plus Task49 VRAC; no Task41 recruitment or Task50 model module enters this dependency chain. No direct expression argument enters G_V. Track1 uses the shared WT reference volume; Track2 uses each cached genotype reference volume, exactly as the declared rest convention requires.

At an exact AE4-null IPR rest with inward basal NKCC N0>0, the added source is deltaN*(1,1,2,0,0), deltaN=(M_beta-1)*N0. Initially V'=0, but under the inherited water law `V''=(Lb+La)*4*deltaN/Vref>0`. Thus finite positive g_V opens the swelling gate at sufficiently small positive time. This FORMAL DEDUCTION removes the old exact deadlock conditionally; it is not sustained physiological validation. The stationary identity remains `J_CaCC+J_VRAC=6P-H`.

## Fixed inference, observations and execution budget

Fallback fixes M_beta=2.5 and uses exactly one deterministic scalar conductance solve. Its sole residual is the WT Track2 IPR mean relative swelling over **300–600 s minus0.125**. This broad late-window amplitude is an explicit ASSUMPTION/SENSITIVITY, not a source-measured window or timing fit; it extends past the graphical trace. The narrative SEM does not quantify mapping/species uncertainty. No time point is selected after seeing output.

Conductance is nonnegative. Numerically bracket the one scalar residual at g=0 and1e-6 S (a solver bracket, not a biological prior or accepted production parameter), then Brent root with absolute g tolerance1e-13 S and at most40 iterations if bracketed. No grid, multistart or second fit. An unbracketed/inadmissible fallback is published honestly; it is not a global impossibility theorem. The fixed M_beta=3 sensitivity uses the SAME g, with no recalibration. The source figure discrepancy is reported, not an alternative calibration target.

All trajectories use the inherited600 s horizon, protocol onset at0+ (start1e-6 s), unchanged Ca0.25 uM when CCh and0.058 otherwise, beta0/1, one scientific process and one BLAS thread. Track1 uses inherited Radau rtol1e-7, amount/regulation atol1e-10, volume atol1e-12, maxstep2 s. Track2 and calibration use BDF rtol1e-8, same component atols, maxstep2 s. Samples at integer seconds plus the initial point; positivity/current/charge/carbon/water/speciation diagnostics are recorded. Numerical solver/conservation failure blocks that output; physiological departures are retained and flagged, not silently relaxed or discarded.

After51B, calibration calls plus three independent WT beta controls (standard IPR reused from root, TMEM16A-off, VRAC-off) are permitted; the primary fit and all attempts are logged. After51C, Track1 has3 central and3 predeclared M_beta=3 cases. Track2 has9 central protocol cases; WT IPR may reuse the calibration trajectory, and3 CCh cases may reuse hash-verified parent results if their complete diagnostic needs are met. The M_beta=3 Track2 IPR cases for WT and AE4KO check source sensitivity. No new resting solve: reuse exact cached states under beta-zero identity. No post-reveal computation changes parameters.

Freeze cumulative summaries at60,120,...600 s, mean60–600 s flow and endpoint flow. Report Cl concentration derivatives with dilution, initial direction, and fixed60–600 s concentration slope only as CONDITIONAL diagnostics. Absolute SPQ slopes are null without the unreported optical calibration/recovery window. Track1/Track2 outputs are never substituted for each other. Cached WT serves AE2 control with strain limitation explicit. Held-out contrasts are computed only after remote51D. Identification failure leads to explicit unavailable predictions, never a fabricated conductance or old Task49 conclusion.
