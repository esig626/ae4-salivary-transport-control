"""Task 21 WT physiology contract. This module performs no optimization."""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import subprocess

from .calibration import WTCalibrationSpec
from .run_calcium_fast_screen import load_freeze, write_json
from .task18_fixed_wt import now
from .validation import CONSERVATION_RESIDUAL_TOLERANCES, PRODUCTION_RADAU, sha256_file

REPO = Path(__file__).resolve().parents[2]
OUT = REPO/'results/21_pooled_ae4_wt_physiology'
ANALYSIS = REPO/'analysis/21_pooled_ae4_wt_physiology'
BRANCH = 'codex/task-21-pooled-ae4-wt-physiology-recalibration'
INTAKE = '0009581c6dd7f855570a5e1bc45d591be75e32dd'
BASE = '3f41579ef81f1e73cc37dc1ce4ab3380f98c4659'
CAPACITIES = ('AE4','NKCC1','AE2','NHE1','pump_apical','pump_basolateral',
              'K_apical','K_basolateral','CaCC','CO2_basolateral','CO2_apical',
              'paracellular_Na','paracellular_K','paracellular_Cl','paracellular_HCO3')
CALCIUM = (0.10,0.25,0.50)
PHYSIOLOGY = {
    'cl_i_mM':{'mean':50.10,'sem':1.50,'lower':47.10,'upper':53.10,'role':'measured WT target'},
    'ph_i':{'mean':6.91,'sem':0.07,'lower':6.77,'upper':7.05,'role':'measured WT target'},
    'na_i_mM':{'lower':2.0,'upper':60.0,'role':'broad physiology/provenance envelope',
        'reference_values_mM':[15.5,20.0,25.0],
        'source':'model/parameters.md Table 1 lineage: Grinstein and Foskett (1990), 20 mM; 2018 model 25 mM. Pena-Munzenmayer 2016 thermodynamic example: 15.5 mM, transferred/calculated from rat sublingual acini.',
        'rationale':'Round-number broad low-intracellular-Na envelope, spanning well below the cited reference values to severalfold above them; endpoints are modeling choices, not measured native SMG limits.',
        'direct_native_smg_interval_available':False},
    'k_i_mM':{'lower':60.0,'upper':200.0,'role':'broad physiology/provenance envelope',
        'reference_values_mM':[120.0,140.0,150.0],
        'source':'model/parameters.md Table 1 lineage: Pedersen and Petersen (1973), 120 mM; Pena-Munzenmayer 2016 Introduction/Discussion: intracellular K around 140-150 mM.',
        'rationale':'Retains a K-rich intracellular compartment but allows about half the 120 mM lineage reference and a generous upper margin over 140-150 mM. Endpoints are broad screening assumptions.',
        'direct_native_smg_interval_available':False},
    'volume_i_pL':{'lower':0.3,'upper':5.0,'role':'broad physiology/provenance envelope',
        'reference_value_pL':1.3,
        'source':'model/parameters.md and analysis/13B_modern_full_model/wt_calibration.md: Palk 2010/2018 published-model 1.3 pL convention, not a native SMG uncertainty estimate.',
        'rationale':'Allows roughly 0.23-3.85 times the lineage volume for an individual acinar cell; deliberately wider than the old 0.5-3 pL numerical box and removes the old narrow volume target.',
        'direct_native_smg_interval_available':False},
}


def git(*args):
    return subprocess.check_output(['git',*args],cwd=REPO,text=True).strip()


def prepare():
    if git('branch','--show-current') != BRANCH:
        raise AssertionError('Wrong Task 21 branch')
    path=OUT/'precalibration_contract.json'
    if path.exists():
        raise FileExistsError('Stage 0 contract is immutable')
    manifest,verification=load_freeze()
    spec=WTCalibrationSpec()
    sources=['model/parameters.md','analysis/13B_modern_full_model/evidence_ledger.md',
        'analysis/13B_modern_full_model/wt_calibration.md',
        'analysis/12_ae4_mechanism_reconstruction/experimental_evidence.md',
        'analysis/13_state_resolved_ae4/evidence_freeze.md',
        'src/modern_full_model/pooled_ae4.py', 'src/modern_full_model/model.py',
        'src/modern_full_model/acid_base.py','src/modern_full_model/membranes.py',
        'src/modern_full_model/parameters.py','src/modern_full_model/water.py',
        'src/modern_full_model/states.py','src/modern_full_model/camp_pka.py',
        'src/modern_full_model/nkcc_stimulation.py','src/modern_full_model/validation.py']
    write_json(path,{
        'contract_id':'TASK21_WT_PHYSIOLOGY_PREOPTIMIZATION_V1','created_utc':now(),
        'branch':BRANCH,'intake_commit':INTAKE,'inheritance_commit':BASE,
        'mechanism':'POOLED_CATION_112_NO_SLIP','mechanism_unchanged':True,
        'physiology':PHYSIOLOGY,
        'measured_screening_band_interpretation':'plus/minus two reported SEM; not population confidence intervals',
        'seed_root_ids':sorted(manifest['roots']),
        'seed_role':'Independent numerical starts and per-seed capacity references only; no legacy-state fit penalty',
        'calcium_uM':CALCIUM,'duration_s':600.0,'sample_step_s':1.0,
        'adjustable_capacities':CAPACITIES,
        'capacity_coordinates':'Independent positive apical/basolateral pump and K capacities; total/fractions derived',
        'capacity_domain':'Nonnegative physical capacities; log-fold objective on positive capacities. No inherited capacity box is a hard bound.',
        'positive_affinity_required':True,'strictly_positive_ae4_cl_loading_required':True,
        'ae4_chloride_share_target':None,'legacy_ae4_flux_matching':False,
        'objective_lexicographic':['sum of squared WT Cl and pH errors normalized by reported SEM',
            'minimize largest absolute log-fold capacity change',
            'minimize sum of squared log-fold capacity changes'],
        'state_reference_penalty':False,
        'deduplication':{'max_relative_difference':1e-6,'compare':['complete conserved state','full calibrated parameter vector'],
            'rule':'Componentwise difference <= 1e-6 times max(abs(first),abs(second),1e-12) in native units; different immutable parameter backgrounds remain distinct'},
        'retain':'Every distinct admissible calibrated solution passing the same WT-only gates; preserve failed seed decisions too',
        'fixed_parameters':['bath','acid-base constants','buffer amounts and bookkeeping','geometry parameters and equations',
            'water/permeability and outflow parameters and equations','beta/cAMP/PKA law and kinetics','calcium inputs',
            'pump stoichiometry','source signs','pooled AE4 law, equal selectivity and 1:1:2 working stoichiometry'],
        'numerical_tolerances':{'root_scaled':spec.root_scaled_tolerance,'omitted_row_raw':spec.omitted_row_raw_tolerance,
            'charge_fmol':spec.charge_tolerance_fmol,'current_A':spec.current_tolerance_A,
            'independent_rhs_scales':spec.independent_rhs_scales,
            'conservation':CONSERVATION_RESIDUAL_TOLERANCES},
        'solver':asdict(PRODUCTION_RADAU),
        'stage0_push_required_before_optimization':True,
        'wt_checkpoint_required_before_genotypes':True,
        'wt_checkpoint_contents':['all seed decisions','retained states and parameters','resting flux/conservation ledger',
            'all three WT dynamics and gates','complete candidate manifest'],
        'genotypes_after_wt_checkpoint':{'AE4_expression':0.05,'AE2_expression':0.0,
            'capacity_refits':0,'resting_method':'Connected production continuation, no compensation or refit'},
        'no_valid_wt_policy':'Stop before genotype evaluation; distinguish structural certificate from numerical nonconvergence',
        'experimental_ae4_deficit_used_in_contract_or_calibration':False,
        'genotype_results_used':False,'new_wt_optimization_results_inspected':False,
        'interpretation_firewall':'Experimental loss magnitude is context only after all genotype results are frozen',
        'source_hashes':{name:sha256_file(REPO/name) for name in sources},
        'inherited_wt_input_verification':verification,
        'seed_hashes':{r:{k:manifest['roots'][r][k] for k in ('core_state_sha256','whole_cell_parameters_sha256','ae4_parameters_sha256')} for r in manifest['roots']},
        'archive_tree_sha':git('rev-parse',INTAKE+':archive'),
    })
    print('Stage 0 declared; commit, push and remotely verify before optimization.')


if __name__=='__main__':
    prepare()
