"""Verify final report artifacts against the independently replayed mathematics."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

P = Path(__file__).resolve().parent
ROOT = P.parents[1]
REL = str(P.relative_to(ROOT))
NUMERICAL_SHA = 'f783785a863440df459e9c5530beba4c7eb59110'
INITIAL_SHA = 'b7b775a0277d93263c4de339745806e5658ec9c2'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(name):
    return json.loads((P / name).read_text())


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)


def main():
    replay = read('output/clean_recomputation_verification.json')
    assert replay['status'] == 'passed' and not replay['failures']
    assert replay['numeric_values_compared'] == 331157
    assert replay['bit_identical_numeric_values'] == 331157
    assert replay['files_compared'] == 47
    assert sha(P / 'recompute_pinned_checkpoint.py') == replay['comparison_helper_sha256']
    numerical_scripts = ['core.py', 'run_analysis.py', 'verify.py', 'verify_exact.py',
                         'equilibrium_review.py', 'sensitivity_analysis.py',
                         'inverse_analysis.py', 'inverse_verify.py']
    for name in numerical_scripts:
        assert (P / name).read_bytes() == git('show', f'{NUMERICAL_SHA}:{REL}/{name}'), name

    # Report preparation cannot silently alter the verified numerical payloads.
    exclusions = {('output/completion_verification.json', 'analysed_commit'),
                  ('output/execution_summary.json', 'elapsed_s')}
    payload_hashes = {}
    for item in replay['files']:
        name = item['path']
        old = git('show', f'{NUMERICAL_SHA}:{REL}/{name}')
        new = (P / name).read_bytes()
        if name.endswith('.json'):
            left, right = json.loads(old), json.loads(new)
            for filename, key in exclusions:
                if name == filename:
                    if key == 'analysed_commit':
                        subprocess.run(['git', 'merge-base', '--is-ancestor',
                                        left[key], right[key]], cwd=ROOT, check=True)
                    left.pop(key); right.pop(key)
            assert left == right, name
        else:
            assert old == new, name
        payload_hashes[name] = sha(P / name)

    gates = read('output/sensitivity_full_gate_verification.json')
    clean_gates = read('output/clean_sensitivity_gate_verification.json')
    assert gates['status'] == clean_gates['status'] == 'passed'
    assert len(gates['cases']) == 384
    assert clean_gates['byte_identical'] and not clean_gates['failures']
    assert not clean_gates['excluded_fields']
    assert clean_gates['expected_json_sha256'] == sha(P / 'output/sensitivity_full_gate_verification.json')
    assert clean_gates['supplementary_script_sha256'] == sha(P / 'verify_sensitivity_gates.py')

    log = (P / 'build.log').read_text()
    assert not re.search(r'undefined|multiply defined|Overfull|LaTeX Warning|Package \w+ Warning', log)
    info = subprocess.check_output(['pdfinfo', str(P / 'AE4_full_system_mathematics.pdf')], text=True)
    pages = int(re.search(r'^Pages:\s+(\d+)', info, re.M).group(1))
    assert pages == 26
    review = read('output/report_independent_review.json')
    assert review['status'] == 'passed' and not review['unresolved_issues']
    assert review['covered_topic_count'] == 17
    assert review['pdf_sha256'] == sha(P / 'AE4_full_system_mathematics.pdf')
    report_files = sorted(P.glob('*.tex')) + [P / 'sources.bib', P / 'build.sh', P / 'build.log',
                                             P / 'AE4_full_system_mathematics.pdf']
    report_files += sorted((P / 'output').glob('*.tex'))
    report_files += sorted((P / 'output').glob('*.pdf')) + sorted((P / 'output').glob('*.png'))
    protected = git('diff', '--name-only', INITIAL_SHA, '--').decode().splitlines()
    assert all(name.startswith(REL + '/') for name in protected), protected
    summary = dict(status='passed', numerical_source_commit=NUMERICAL_SHA,
                   analysis_scripts_identical_to_numerical_checkpoint=numerical_scripts,
                   verified_numerical_payload_sha256=payload_hashes,
                   pdf_pages=pages, figures=5, references_resolved=True,
                   independently_reviewed_required_topics=17,
                   visual_review_pages=list(range(1, 27)),
                   visual_review_scope='Root inspected pages 1–8 and all five standalone figures; independent reviewers inspected pages 9–16 and 17–26, and all changed pages were reinspected.',
                   clean_replay_numeric_values=331157,
                   supplementary_gate_roots=384,
                   supplementary_clean_numeric_values=clean_gates['numeric_values_compared_exactly'],
                   report_artifact_sha256={str(path.relative_to(P)): sha(path) for path in report_files},
                   protected_path_differences=[],
                   scope='Artifact consistency and unchanged numerical payload checks. Visual and scientific review is recorded separately; no new scientific model solves.')
    (P / 'output/report_verification.json').write_text(json.dumps(summary, indent=2, sort_keys=True)+'\n')
    print(f'REPORT_VERIFICATION passed: {pages} pages; 47 numerical payloads unchanged; 384 supplementary gate roots.')


if __name__ == '__main__':
    main()
