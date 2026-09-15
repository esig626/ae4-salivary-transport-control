#!/usr/bin/env python3
"""Check the assembled report. This does not rerun the scientific simulations."""
from pathlib import Path
import hashlib
import json
import re
import subprocess

HERE = Path(__file__).resolve().parent

def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def load(p):
    return json.loads(p.read_text())

def main():
    manifest = load(HERE/'SOURCE_MANIFEST.json')
    bad = [r['report_path'] for r in manifest['files'] if digest(HERE/r['report_path']) != r['sha256']]
    if bad:
        raise RuntimeError('Changed source evidence: '+repr(bad))
    coverage = load(HERE/'generated/coverage.json')
    expected = {'parameter_inventory_records':137, 'parameter_active_flags':109,
                'equilibrium_cases':26, 'complete_jacobians':26, 'jacobian_entries':3146,
                'eigenvalues':286, 'state_entries':338, 'expression_derivative_points':4,
                'modal_rows':22, 'uncertain_parameters':16, 'parameter_metric_rows':208,
                'ordered_inverse_points':17, 'inverse_parameter_observable_pairs':16}
    for k,v in expected.items():
        if coverage[k] != v:
            raise RuntimeError(f'Coverage mismatch: {k}: {coverage[k]} != {v}')
    log=(HERE/'build/report.log').read_text(errors='replace')
    failures=[]
    patterns=[r'undefined citations',r'undefined references',r'Citation .+ undefined',
              r'Reference .+ undefined',r'Label\(s\) may have changed',r'LaTeX Error',
              r'Overfull \\[hv]box',r'Missing character:']
    for pattern in patterns:
        if re.search(pattern,log,re.I):
            failures.append(pattern)
    if failures:
        raise RuntimeError('Report log requires correction: '+repr(failures))
    pdf=HERE/'AE4_mathematical_analysis_report.pdf'
    info=subprocess.check_output(['pdfinfo',str(pdf)],text=True)
    pages=int(re.search(r'^Pages:\s+(\d+)',info,re.M).group(1))
    if pages<40:
        raise RuntimeError('Unexpectedly short report')
    text=subprocess.check_output(['pdftotext','-layout',str(pdf),'-'],text=True)
    if '??' in text:
        raise RuntimeError('Unresolved reference marker in PDF')
    topics=['carbonate','alkalinity','nondimensional','Jacobian','compensation','inverse',
            'provenance','thermodynamic','replay','26','137','331157','459','345','89.5']
    absent=[t for t in topics if t.lower() not in text.lower()]
    if absent:
        raise RuntimeError('Missing coverage terms: '+repr(absent))
    aux=(HERE/'build/report.aux').read_text()
    cites=set(k.strip() for group in re.findall(r'\\citation\{([^}]+)\}',aux) for k in group.split(','))
    bibkeys=set(re.findall(r'\\bibcite\{([^}]+)\}',aux))
    if cites-bibkeys:
        raise RuntimeError('Unresolved bibliography keys: '+repr(cites-bibkeys))
    report_inputs={p.relative_to(HERE).as_posix():digest(p) for p in sorted(HERE.glob('*.tex'))}
    result={
      'status':'REPORT_BUILD_AND_SOURCE_INTEGRITY_VERIFIED',
      'pdf':{'path':pdf.name,'pages':pages,'bytes':pdf.stat().st_size,'sha256':digest(pdf)},
      'pinned_inputs':manifest['pinned_inputs'],
      'unchanged_evidence_files':len(manifest['files']),
      'coverage':coverage,'resolved_citation_keys':len(cites),
      'report_inputs':report_inputs,
      'no_unresolved_citations_or_references':True,
      'no_overfull_boxes_or_missing_characters':True,
      'scientific_simulations_rerun_for_this_report':False,
      'scientific_replay_claims':'Transcribed from preserved Task 44 verification receipts, not newly executed here.',
      'visual_review':'Separate manual rendered page review required; this script does not assert visual review.',
      'assembly_script_sha256':digest(HERE/'assemble.py'),
      'verification_script_sha256':digest(Path(__file__))}
    (HERE/'BUILD_VERIFICATION.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'pages':pages,'source_files':len(manifest['files']),'citations':len(cites),'pdf_sha256':digest(pdf)},sort_keys=True))

if __name__=='__main__':
    main()
