"""Predeclared Task 22 domain and checkpoint checks; no optimization here."""
from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import subprocess

from .calibration import WTCalibrationSpec
from .run_calcium_fast_screen import load_freeze, write_json
from .task18_fixed_wt import now
from .task21_contract import CAPACITIES, PHYSIOLOGY as INTRACELLULAR
from .task21_calibration import WTProblem
from .validation import CONSERVATION_RESIDUAL_TOLERANCES, PRODUCTION_RADAU, sha256_file

REPO = Path(__file__).resolve().parents[2]
OUT = REPO / 'results/22_physiological_wt_domain'
ANALYSIS = REPO / 'analysis/22_physiological_wt_domain'
BRANCH = 'codex/task-22-physiological-wt-domain'
INTAKE = 'ecafebd3e3fc4237aeac3b063eb03deb67165955'
BASE = 'db734f52b8d71d670fcf793beab922967212bdd6'
CALCIUM = (0.10, 0.25, 0.50)
PHYSIOLOGY = {**INTRACELLULAR,
    'na_l_mM': {'lower': 100., 'upper': 200.},
    'k_l_mM': {'lower': 1., 'upper': 30.},
    'cl_l_mM': {'lower': 80., 'upper': 180.},
    'ph_l': {'lower': 6.8, 'upper': 8.0},
    'tic_l_mM': {'lower': 1., 'upper': 80.},
    'volume_l_pL': {'lower': .02, 'upper': .50}}
CAPACITY_FOLD = (.01, 100.)
OSMOTIC_RATIO = (.80, 1.20)
ANTI_CANCELLATION_MAX = 100.


def git(*args):
    return subprocess.check_output(['git', *args], cwd=REPO, text=True).strip()


def read(path):
    return json.loads(path.read_text())


def blob_sha(path):
    data = path.read_bytes()
    return hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest()


def check_receipt(receipt, path, *, repo=REPO):
    """Bind local bytes to the independently fetched remote commit and blob."""
    if receipt.get('remote_ref_verified') is not True:
        raise AssertionError('Remote ref was not verified')
    if receipt.get('remote_branch') != BRANCH:
        raise AssertionError('Wrong remote branch')
    if receipt.get('remote_blob_sha') != blob_sha(path):
        raise AssertionError('Checkpoint bytes differ from remote blob')
    checkpoint = receipt['checkpoint_sha']
    relative = str(path.relative_to(repo))
    if git('rev-parse', checkpoint + ':' + relative) != receipt['remote_blob_sha']:
        raise AssertionError('File absent from recorded checkpoint')
    git('merge-base', '--is-ancestor', checkpoint, 'HEAD')
    if receipt['remote_verified_utc'] < read(path)['created_utc']:
        raise AssertionError('Checkpoint verification predates its evidence')


def require_stage0(*, allow_frozen=False):
    if git('branch', '--show-current') != BRANCH:
        raise AssertionError('Wrong Task 22 branch')
    path = OUT / 'precalibration_contract.json'
    check_receipt(read(OUT / 'stage0_checkpoint.json'), path)
    contract = read(path)
    for name, digest in contract['source_hashes'].items():
        if sha256_file(REPO / name) != digest:
            raise AssertionError('Frozen input changed: ' + name)
    if git('rev-parse', 'HEAD:archive') != contract['archive_tree_sha']:
        raise AssertionError('Archive tree changed')
    if not allow_frozen and (OUT / 'wt_manifest.json').exists():
        raise AssertionError('WT evidence is frozen; further search prohibited')
    return contract


def authorize_genotypes():
    require_stage0(allow_frozen=True)
    path = OUT / 'wt_manifest.json'
    manifest = read(path)
    check_receipt(read(OUT / 'wt_checkpoint_receipt.json'), path)
    for name, digest in manifest['file_sha256'].items():
        if sha256_file(REPO / name) != digest:
            raise AssertionError('Frozen WT evidence changed: ' + name)
    if not manifest['search_complete']:
        raise AssertionError('WT search incomplete')
    if not manifest['retained_solution_ids']:
        raise RuntimeError('STOP_NO_ADMISSIBLE_WT: genotype evaluation prohibited')
    for cid in manifest['retained_solution_ids']:
        item = read(OUT / 'wt_solutions' / (cid + '.json'))
        if not item['resting_gate_pass']:
            raise AssertionError('Retained solution failed resting gates')
        panel = [r for r in read(OUT / 'wt_dynamic_results.json') if r['solution_id'] == cid]
        if {r['calcium_uM'] for r in panel} != set(CALCIUM) or len(panel) != 3:
            raise AssertionError('Incomplete WT dynamic panel')
        if not all(r['dynamic_gate_pass'] for r in panel):
            raise AssertionError('WT dynamics failed')
    return manifest


def prepare():
    if git('branch', '--show-current') != BRANCH:
        raise AssertionError('Wrong branch')
    path = OUT / 'precalibration_contract.json'
    if path.exists():
        raise FileExistsError('The contract is immutable')
    manifest, verification = load_freeze()
    refs = {root: WTProblem(manifest, root).reference for root in sorted(manifest['roots'])}
    if not all(value > 0 for row in refs.values() for value in row.values()):
        raise AssertionError('Unexpected zero reference; declare it before optimization')
    sources = ['AGENTS.md', 'prompts/22_physiologically_constrained_wt_recalibration.md',
        'model/parameters.md', 'analysis/13B_modern_full_model/evidence_ledger.md',
        'analysis/13B_modern_full_model/wt_calibration.md',
        'analysis/21_pooled_ae4_wt_physiology/wt_physiology_contract.md',
        'analysis/11_forensic_reconstruction/literature_crosscheck.md']
    sources += [str(p.relative_to(REPO)) for p in sorted((REPO / 'src/modern_full_model').glob('*.py'))
                if not p.name.startswith('task22')]
    witnesses = sorted((REPO / 'results/21_pooled_ae4_wt_physiology/wt_candidates').glob('*.json'))
    sources += [str(p.relative_to(REPO)) for p in witnesses]
    spec = WTCalibrationSpec()
    write_json(path, dict(contract_id='TASK22_PHYSIOLOGICAL_WT_DOMAIN_V1', created_utc=now(),
        branch=BRANCH, intake_commit=INTAKE, inheritance_commit=BASE,
        mechanism='POOLED_CATION_112_NO_SLIP', mechanism_unchanged=True,
        physiology=PHYSIOLOGY,
        luminal_rationale='Broad primary-acinar, near-isotonic/plasma-like physiology/model screens explicitly prescribed by Task 22; not matched mouse-SMG measurement uncertainty intervals.',
        physiology_sources=[
            dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC3517652/',
                 citation='Patterson et al. 2012, doi:10.1152/ajpgi.00364.2011',
                 basis='Primary acinar fluid is plasma-like and isotonic; subsequent ductal processing changes composition.'),
            dict(url='https://link.springer.com/article/10.1007/s11538-022-01041-3',
                 citation='Su et al. 2022, discussion of primary-fluid composition; model lineage already cited in repository literature_crosscheck.md',
                 basis='Reports primary-fluid ranges Na 140-150, K 5-15, Cl 80-130, HCO3 30-60 mM across cited animal/gland/stimulation contexts (Mangos 1973 and Young 1971). These are contextual ranges, not new native-SMG data or exact Task 22 targets.'),
            dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC5792321/',
                 citation='Vera-Siguenza et al. 2018, Table 1 and water-balance discussion; model/parameters.md',
                 basis='Lineage lumen Na 118.7, K 5.6, Cl 124.3 mM; model osmolarities bath/cell/lumen 292.6/296.6/297.4 mOsm. The inherited production model, not these old osmolarities, determines the actual ratios.')],
        luminal_ph_tic_volume_rationale='pH 6.8-8.0 and TIC 1-80 mM are broad prescribed acid-base envelopes (TIC is not bicarbonate alone). Volume .02-.50 pL brackets the .10 pL production reference and .026 pL old .02*1.3 geometry convention; neither reference is an uncertainty interval.',
        osmotic_ratio_bounds=OSMOTIC_RATIO,
        osmolarity_implementation='ModernFullModel.observables -> water.compartment_osmolarity_mOsm; cell finite buffer and other impermeant particles included exactly once; use production observables for acceptance.',
        adjustable_capacities=CAPACITIES, capacity_reference_by_root=refs,
        positive_capacity_fold_bounds=CAPACITY_FOLD, fractions_bounds=[0., 1.],
        capacity_coordinates='Task 21 fifteen existing coordinates; independent apical/basolateral pump and K capacities, with totals and fractions derived. Totals also remain within the fold box by positivity.',
        zero_reference_rule='No absolute nonzero exceptions declared. All zero-reference parameters remain exactly zero; all fifteen adjustable coordinates have positive references in all ten roots.',
        capacity_rationale='Predeclared broad two-orders-of-magnitude limits in each direction prevent million-fold optimizer escape; not inferred biological uncertainty.',
        anti_cancellation=dict(maximum_ratio=ANTI_CANCELLATION_MAX,
            positive_pool='max(J_AE4_Cl,0)+max(2*J_NKCC1,0)+max(J_AE2,0), evaluated at the candidate WT state using signed cell-source semantics from Tasks 18/20',
            fluxes='All individual transporter cycle fluxes and individual Na/K/Cl/TIC/alkalinity source components, channels converted from A to fmol/s, CO2 exchanges, paracellular sources, and lumen outflow sources. Test each before net cancellation. Water is reported separately in pL/s, never divided by an ion pool in fmol/s.',
            rationale='Broad anti-pathology and numerical-conditioning screen, not measured biology; no genotype fitting.'),
        seed_root_ids=sorted(manifest['roots']),
        task21_witness_files=[str(p.relative_to(REPO)) for p in witnesses],
        seed_rule='Every inherited root and every archived Task 21 resting witness is a numerical start; project infeasible coordinates into the declared box. Additional WT-only starts allowed; no state-distance objective.',
        objective_lexicographic=['((Cl_i-50.10)/1.50)^2+((pH_i-6.91)/.07)^2',
            'max absolute natural log of adjustable positive capacity/reference',
            'sum of squared natural logs of adjustable positive capacity/reference'],
        deduplication=dict(relative_tolerance=1e-6, floor_native_units=1e-12,
            fields=['complete_state', 'whole_cell_parameters', 'ae4_parameters'],
            retain='All distinct admissible solutions, sorted lexicographically; immutable backgrounds remain distinct.'),
        fixed='Bath chemistry, acid-base laws/constants, source signs, buffers and impermeant amounts, geometry equations/parameters, water/outflow equations/parameters, regulation, pathway inventory, working 1:1:2 stoichiometry and pooled mechanism are inherited unchanged.',
        no_chloride_share_target=True, no_legacy_ae4_flux_match=True,
        positive_affinity_required=True, strictly_positive_ae4_cl_loading_required=True,
        pooled_na_k_same_sign_required=True,
        resting_tolerances=dict(root_scaled=spec.root_scaled_tolerance,
            independent_rhs_scales=spec.independent_rhs_scales,
            omitted_rhs_raw=spec.omitted_row_raw_tolerance, charge_fmol=spec.charge_tolerance_fmol,
            current_A=spec.current_tolerance_A, conservation=CONSERVATION_RESIDUAL_TOLERANCES),
        dynamics=dict(calcium_uM=CALCIUM, duration_s=600., sample_step_s=1.,
            solver=asdict(PRODUCTION_RADAU), full_time_grid_required=True,
            all_existing_numerical_and_conservation_gates_required=True,
            state_charge_tolerance_fmol=spec.charge_tolerance_fmol, regulatory_fractions=[0.,1.],
            rest_physiology_screens_apply_to='Resting WT; dynamic physical/current/conservation/state-charge gates remain inherited.'),
        numerical_improvements='Scaling, continuation, precision and analytic Jacobians allowed; no equation, bound or tolerance changes.',
        stage0_push_and_remote_verification_before_new_optimization=True,
        complete_wt_push_and_remote_verification_before_any_genotype_or_historical_genotype_regression=True,
        empty_wt_policy='STOP before all genotype evaluation; no placeholder genotype comparison table.',
        failure_classes=['NUMERICAL_NONCONVERGENCE','PHYSIOLOGICAL_BOUND_CONFLICT',
            'CAPACITY_BOUND_CONFLICT','THERMODYNAMIC_CONFLICT','CONSERVATION_DYNAMIC_GATE_FAILURE',
            'PROVEN_STRUCTURAL_CONTRADICTION'],
        failure_interpretation='Local failure is not proof. Distinguish proved necessary-condition conflicts from observed active bounds or numerical residuals.',
        genotype_results_used=False, new_task22_optimization_run_or_inspected=False,
        inherited_input_verification=verification,
        source_hashes={name:sha256_file(REPO/name) for name in sources},
        archive_tree_sha=git('rev-parse', 'HEAD:archive')))
    print('Contract written. Commit, push, and remotely verify before optimization.')


if __name__ == '__main__':
    prepare()
