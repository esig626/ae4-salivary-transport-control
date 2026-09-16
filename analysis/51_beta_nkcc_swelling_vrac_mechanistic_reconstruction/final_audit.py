"""File/Git-only final audit. This module never imports the scientific model."""
import hashlib
import json
import subprocess
from checkpoint import HERE, ROOT, OUT, START, BRANCH, git, integrity, dump, sha


def blob(commit, path):
    return subprocess.check_output(['git', 'show', commit + ':' + path], cwd=ROOT)


def main():
    destination = OUT / 'final_integrity_audit.json'
    if destination.exists():
        raise RuntimeError('Final audit already exists')
    assert git('branch', '--show-current') == BRANCH
    records = json.loads((OUT / 'publication_verification.json').read_text())
    assert [r['checkpoint'] for r in records] == ['51A', '51B', '51C', '51D', '51E']
    assert records[-1]['sha'] == git('rev-parse', 'HEAD')
    previous = START
    checks = []
    receipts = {}
    for record in records:
        assert record['verified']
        commit = record['sha']
        assert git('rev-parse', commit + '^') == previous
        assert git('rev-parse', commit + '^{tree}') == record['tree']
        path = str((OUT / 'checkpoints' / (record['checkpoint'] + '.json')).relative_to(ROOT))
        receipt = json.loads(blob(commit, path))
        receipts[record['checkpoint']] = receipt
        assert receipt['parent_sha'] == previous
        assert receipt['branch'] == BRANCH and receipt['requested_start_sha'] == START
        for artifact, expected in receipt['artifact_sha256'].items():
            actual = hashlib.sha256(blob(commit, artifact)).hexdigest()
            assert actual == expected, (record['checkpoint'], artifact)
        checks.append(dict(checkpoint=record['checkpoint'], sha=commit,
                           tree=record['tree'], receipt_artifacts_checked=len(receipt['artifact_sha256']), passed=True))
        previous = commit
    frozen_input_count = integrity()
    allowed_prefix = str(HERE.relative_to(ROOT)) + '/'
    differences = git('diff', '--name-only', START).splitlines()
    assert all(p.startswith(allowed_prefix) or p == 'docs/MANDATORY_RESEARCH_LEDGER.md' for p in differences)
    untracked = git('ls-files', '--others', '--exclude-standard').splitlines()
    assert all(p.startswith(allowed_prefix) for p in untracked)
    old_ledger = blob(START, 'docs/MANDATORY_RESEARCH_LEDGER.md')
    assert (ROOT / 'docs/MANDATORY_RESEARCH_LEDGER.md').read_bytes().startswith(old_ledger)
    protected = {
        '51B': ['task51_model.py'],
        '51C': ['trajectory.py', 'identify.py', 'output/parameter_freeze.json'],
        '51D': ['SOURCE_ARCHITECTURE_FREEZE.md', 'output/pre_AE4_predictions.json'],
        '51E': ['output/heldout_comparison.json'],
    }
    unchanged = {}
    for checkpoint, paths in protected.items():
        for relative in paths:
            path = HERE / relative
            digest = sha(path)
            assert digest == receipts[checkpoint]['artifact_sha256'][str(path.relative_to(ROOT))]
            unchanged[relative] = dict(since=checkpoint, sha256=digest)
    parameters = json.loads((OUT / 'parameter_freeze.json').read_text())
    predictions = json.loads((OUT / 'pre_AE4_predictions.json').read_text())
    comparison = json.loads((OUT / 'heldout_comparison.json').read_text())
    assert parameters['g_V_S'] is None and not parameters['heldout_data_used']
    assert not predictions['heldout_comparison_performed'] and not comparison['retuned']
    assert comparison['immutable_reveal_boundary']['sha'] == records[3]['sha']
    for relative, expected in comparison['post_reveal_frozen_hashes'].items():
        assert sha(HERE / relative) == expected
    cached_references = 0
    for row in predictions['predictions']:
        if row['status'] == 'EXACT_PARENT_IVP_CACHED_RESULT_REFERENCE':
            for suffix in ['json', 'csv']:
                assert sha(ROOT / row['parent_' + suffix]) == row['parent_' + suffix + '_sha256']
                cached_references += 1
    refs = json.loads((HERE / 'input/remote_refs_before_51F.json').read_text())
    assert refs['verified'] and refs['archive']['sha'] == 'b16c30094b95f61a73d8f1977cd58e79c7bb50f6'
    assert refs['task51']['sha'] == records[-1]['sha']
    verification = json.loads((OUT / 'implementation_verification.json').read_text())
    trajectories = [json.loads(p.read_text()) for p in sorted((OUT / 'trajectories').glob('*.json'))]
    assert len(trajectories) == 2 and all(t['status'] == 'COMPLETE' and t['physiologically_admissible'] for t in trajectories)
    counts = dict(distinct_software_tests=verification['unique_tests_passed'],
                  test_invocations=verification['test_invocations'],
                  verification_core_evaluations=verification['core_evaluations'],
                  calibration_trajectories=len(trajectories),
                  calibration_core_evaluations=sum(t['core_evaluations'] for t in trajectories),
                  accepted_parameter_fits=0, accepted_model_phenotype_trajectories=0,
                  resting_solves=0, post_reveal_model_evaluations=0,
                  new_model_evaluations_in_final_audit=0)
    counts['total_core_evaluations'] = counts['verification_core_evaluations'] + counts['calibration_core_evaluations']
    data = dict(passed=True, requested_start=START, branch=BRANCH,
                verified_dependency_head=records[-1]['sha'], checkpoint_audits=checks,
                total_historical_receipt_hashes_checked=sum(c['receipt_artifacts_checked'] for c in checks),
                original_frozen_input_hashes_checked=frozen_input_count,
                source_model_parameters_predictions_unchanged=unchanged,
                cached_parent_result_references_checked=cached_references,
                original_ledger_preserved_as_exact_prefix=True,
                original_ledger_prefix_sha256=hashlib.sha256(old_ledger).hexdigest(),
                outside_task_changes=['docs/MANDATORY_RESEARCH_LEDGER.md'],
                all_other_original_repository_files_unchanged=True,
                archive_remote_verification=refs['archive'],
                numerical_work_counts=counts,
                final_publication_rule='51F ref/commit/tree/receipt verified independently after publication; final SHA recorded outside the branch to avoid a seventh commit')
    dump(destination, data)
    print(json.dumps(dict(passed=True, historical_receipt_hashes=data['total_historical_receipt_hashes_checked'],
                          frozen_inputs=frozen_input_count, cached_parent_references=cached_references, counts=counts)))


if __name__ == '__main__':
    main()
