"""Checkpoint verification and reporting; no genotype simulation is hidden here."""
import hashlib
import math
from .task21_calibration import read, OUT, REPO, write_json, write_rows, inverse, git
from .validation import sha256_file


def checked_checkpoint():
    manifest_path=OUT/'wt_manifest.json';manifest=read(manifest_path)
    receipt=read(OUT/'wt_checkpoint_receipt.json')
    data=manifest_path.read_bytes()
    blob=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    if not receipt['remote_ref_verified'] or receipt['remote_manifest_blob_sha']!=blob:
        raise AssertionError('WT remote checkpoint is absent or changed')
    if not receipt['checkpoint_sha'] or receipt['remote_ref_verified_utc']<manifest['created_utc']:
        raise AssertionError('WT checkpoint ordering is invalid')
    if git('rev-parse',receipt['checkpoint_sha']+':results/21_pooled_ae4_wt_physiology/wt_manifest.json')!=blob:
        raise AssertionError('WT manifest is not in the recorded checkpoint commit')
    git('merge-base','--is-ancestor',receipt['checkpoint_sha'],'HEAD')
    for path,digest in manifest['file_sha256'].items():
        if sha256_file(REPO/path)!=digest:raise AssertionError(f'Frozen WT file changed: {path}')
    return manifest,receipt


def authorize_genotypes():
    manifest,receipt=checked_checkpoint()
    if not manifest['retained_solution_ids']:
        raise RuntimeError('STOP_NO_ADMISSIBLE_WT: genotype evaluation is prohibited')
    return manifest,receipt


def matched_ratio(wt,loss):
    if (wt['solution_id'],wt['calcium_uM'])!=(loss['solution_id'],loss['calcium_uM']):
        raise ValueError('Genotype and WT are not a matched pair')
    if not (wt['dynamic_gate_pass'] and loss['dynamic_gate_pass']):return None
    a,b=loss['Q_0_600_pL'],wt['Q_0_600_pL']
    if not (math.isfinite(a) and math.isfinite(b) and a>=0 and b>0):
        raise ValueError('Invalid secretion total')
    return a/b


def report_stopped_genotypes():
    manifest,receipt=checked_checkpoint()
    if manifest['retained_solution_ids']:
        raise AssertionError('This unavailable-result writer cannot replace eligible genotype tests')
    rows=[]
    for root in manifest['seed_root_ids']:
        for calcium in (.10,.25,.50):
            rows.append(dict(root_id=root,calcium_uM=calcium,status='NOT_RUN_NO_ADMISSIBLE_WT',
                WT_Q_0_600_pL=None,AE4_005_Q_0_600_pL=None,AE2_loss_Q_0_600_pL=None,
                AE4_005_over_WT=None,AE2_loss_over_WT=None,comparison_available=False))
    write_rows(OUT/'final_comparison.csv',rows)
    for filename,genotype in (('ae4_5pct_results.csv','AE4_0.05'),('ae2_loss_results.csv','AE2_loss')):
        write_rows(OUT/filename,[dict(root_id=r['root_id'],calcium_uM=r['calcium_uM'],genotype=genotype,
            status=r['status'],evaluated=False,capacity_refits=0,Q_0_600_pL=None) for r in rows])
    write_json(OUT/'final_summary.json',dict(completed_utc=inverse.now(),
        classification=manifest['classification'],retained_WT_count=0,
        rest_valid_candidate_count=len(manifest['rest_valid_candidate_ids']),
        tested_WT_candidate_count=len(manifest['tested_candidate_ids']),
        WT_dynamic_conditions=manifest['wt_dynamic_conditions'],WT_dynamic_valid=manifest['wt_dynamic_valid'],
        genotype_evaluations=0,genotype_specific_refits=0,experimental_loss_magnitude_compared=False,
        WT_checkpoint_sha=receipt['checkpoint_sha'],scientific_result='WT numerical calibration unresolved; no genotype conclusion',
        structural_infeasibility_proved=False,global_capacity_optimality_certified=False))


if __name__=='__main__':report_stopped_genotypes()
