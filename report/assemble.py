#!/usr/bin/env python3
"""Assemble a report from pinned evidence, without executing the model.

All numerical values are read from archived outputs. Formatting, unit display
conversions, table coverage checks and file hashing are the only calculations.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import re
import shutil
import subprocess
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PINS = {
    'task43': ('547f113d9123ab1976774639a69483faf40cef34', 'analysis/43_parameter_provenance'),
    'task44': ('e5fa9840bf96be3147c6118daee4942468c78f8b', 'analysis/44_full_system_mathematics'),
}
GEN = HERE / 'generated'
EVID = HERE / 'evidence'
COVERAGE: dict[str, object] = {}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def dump(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + '\n')


def esc(value: object) -> str:
    s = str(value)
    replacements = {'\\': r'\textbackslash{}', '&': r'\&', '%': r'\%', '$': r'\$', '#': r'\#',
                    '_': r'\_', '{': r'\{', '}': r'\}', '~': r'\textasciitilde{}', '^': r'\textasciicircum{}'}
    s = ''.join(replacements.get(c, c) for c in s)
    for a, b in [('−', '$-$'), ('≤', '$\\leq$'), ('≥', '$\\geq$'), ('→', '$\\to$'),
                 ('×', '$\\times$'), ('α', '$\\alpha$'), ('ρ', '$\\rho$'), ('τ', '$\\tau$'),
                 ('µ', '$\\mu$'), ('μ', '$\\mu$'), ('≈', '$\\simeq$')]:
        s = s.replace(a, b)
    return re.sub(r'\b[0-9a-f]{32,64}\b', lambda m: r'\nolinkurl{' + m.group(0) + '}', s)


def pathtex(value: object) -> str:
    return r'\nolinkurl{' + str(value) + '}'


def num(value: object, digits: int = 6) -> str:
    if value is None:
        return r'\textit{n/a}'
    if isinstance(value, bool):
        return 'Yes' if value else 'No'
    try:
        x = float(value)
    except (TypeError, ValueError):
        return esc(value)
    # Explicit e notation keeps dense scientific tables readable.
    return r'\texttt{' + format(x, f'.{digits}g') + '}'


def read(name: str, group: str = 'task44') -> object:
    return json.loads((EVID / group / name).read_text())


def csvrows(name: str, group: str = 'task44') -> list[dict[str, str]]:
    with (EVID / group / name).open(newline='') as handle:
        return list(csv.DictReader(handle))


def snapshot(name: str, supplied: str | None) -> list[dict[str, object]]:
    revision, original = PINS[name]
    dest = EVID / name
    dest.mkdir(parents=True, exist_ok=True)
    records = []
    if supplied:
        src = Path(supplied).resolve()
        if not src.is_dir():
            raise ValueError('Missing exact source directory: ' + str(src))
        files = [(p.relative_to(src).as_posix(), p.read_bytes()) for p in sorted(src.rglob('*')) if p.is_file()]
    else:
        data = subprocess.check_output(['git', 'archive', '--format=tar', revision, original], cwd=ROOT)
        files = []
        with tarfile.open(fileobj=io.BytesIO(data), mode='r:') as archive:
            for member in archive.getmembers():
                if not member.isfile():
                    continue
                prefix = original + '/'
                if not member.name.startswith(prefix):
                    raise ValueError('Unexpected archived path: ' + member.name)
                relative = member.name[len(prefix):]
                if '..' in Path(relative).parts:
                    raise ValueError('Unsafe archived path')
                stream = archive.extractfile(member)
                assert stream is not None
                files.append((relative, stream.read()))
    if not files:
        raise ValueError('No source files were retrieved for ' + name)
    for relative, data in files:
        target = dest / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
        assert target.read_bytes() == data
        records.append(dict(report_path=target.relative_to(HERE).as_posix(), source_commit=revision,
                            source_path=original + '/' + relative, size_bytes=len(data), sha256=sha(data)))
    return records


INPUT = re.compile(r'\\input\{([^{}]+)\}')
CITE = re.compile(r'(\\cite[a-zA-Z]*\*?(?:\[[^\]]*\]){0,2})\{([^{}]+)\}')
REF = re.compile(r'(\\(?:label|eqref|cref|Cref|ref|pageref))\{([^{}]+)\}')


def expand_inputs(text: str, folder: Path, chain: tuple[str, ...] = ()) -> str:
    def load(match: re.Match[str]) -> str:
        filename = match.group(1)
        p = folder / (filename if Path(filename).suffix else filename + '.tex')
        if str(p) in chain:
            raise ValueError('Circular TeX input')
        return '\n% Expanded source: ' + p.relative_to(HERE).as_posix() + '\n' + expand_inputs(p.read_text(), p.parent, chain + (str(p),)) + '\n'
    return INPUT.sub(load, text)


def body(group: str, prefix: str) -> None:
    folder = EVID / group
    original = (folder / 'report.tex').read_text()
    start = original.index(r'\section{')
    finish = original.index(r'\bibliographystyle')
    text = original[start:finish].replace(r'\appendix', '% Appendix integrated in combined section sequence')
    if group == 'task43':
        text = (folder / 'output/audit_counts.tex').read_text() + '\n' + text
    text = expand_inputs(text, folder)
    text = CITE.sub(lambda m: m.group(1) + '{' + ','.join(prefix + k.strip() for k in m.group(2).split(',')) + '}', text)
    text = REF.sub(lambda m: m.group(1) + '{' + ','.join(prefix + k.strip() for k in m.group(2).split(',')) + '}', text)
    text = re.sub(r'(\\includegraphics(?:\[[^\]]*\])?)\{([^{}]+)\}',
                  lambda m: m.group(1) + '{evidence/' + group + '/' + m.group(2) + '}', text)
    (GEN / (group + '_body.tex')).write_text(text)
    COVERAGE[group + '_full_exposition'] = dict(original_report_sha256=sha(original.encode()),
        all_sections_retained=True, original_section_count=len(re.findall(r'\\section\{', original[start:finish])),
        integrated_section_count=len(re.findall(r'\\section\{', text)), expanded_tex_sha256=sha(text.encode()))


def bibliography() -> None:
    combined = ['% Original bibliography records retained with source namespaces.\n']
    for group, prefix in [('task44', 'm44'), ('task43', 'p43')]:
        text = (EVID / group / 'sources.bib').read_text()
        text = re.sub(r'(@\w+\s*\{)([^,\s]+)(\s*,)', lambda m: m.group(1) + prefix + m.group(2) + m.group(3), text)
        combined.append('% ' + group + '\n' + text)
    (GEN / 'references.bib').write_text('\n'.join(combined))


def table(title: str, headers: list[str], rows: list[list[str]], spec: str | None = None, tiny: bool = False) -> str:
    if spec is None:
        spec = 'l' + 'r' * (len(headers) - 1)
    size = r'\scriptsize' if tiny else r'\footnotesize'
    head = ' & '.join(headers) + r' \\' + '\n'
    return ('\n\\begingroup' + size + '\n\\setlength{\\tabcolsep}{3pt}\n'
            + '\\begin{longtable}{' + spec + '}\n\\caption{' + title + '}\\\\\n\\toprule\n'
            + head + '\\midrule\\endfirsthead\n\\toprule\n' + head + '\\midrule\\endhead\n'
            + '\\bottomrule\\endfoot\n' + ''.join(' & '.join(r) + r' \\' + '\n' for r in rows)
            + '\\end{longtable}\n\\endgroup\n')


def source(name: str, group: str = 'task44') -> str:
    return '\n\\sourcefile{evidence/' + group + '/' + name + '}\n'


def numerical_appendix() -> None:
    roots = read('output/equilibria.json')['roots']
    review = read('output/equilibrium_review.json')
    rr = csvrows('output/equilibrium_summary.csv')
    assert len(roots) == len(rr) == len(review['equilibrium_rows']) == 26
    out = [r'\section{Complete equilibrium case tables}',
           'All 26 cases are printed, including the four physiological exclusions and the repeated null systems. '
           'Numeric entries use standard e notation where needed: 1.2e-5 means $1.2\\times10^{-5}$. '
           'Tables round for reading; every full precision value remains in the unchanged electronic source. '
           'The equilibrium evaluator time field is not an equilibration time.']
    def obs_table(title: str, fields: list[tuple[str, str]], tiny: bool = False) -> str:
        return table(title, ['Class', '$e_4$'] + [h for _, h in fields],
                     [[esc(r['mechanism']), num(r['ae4_expression'], 3)] + [num(r[k], 6 if not tiny else 5) for k, _ in fields] for r in rr], tiny=tiny)
    out += [obs_table('Stationary flow and intracellular ions. Flow is in pL/s; concentrations are in mM.',
                     [('q_pL_s', '$Q$'), ('na_i_mM', '$[Na]_i$'), ('k_i_mM', '$[K]_i$'), ('cl_i_mM', '$[Cl]_i$')]),
            obs_table('Intracellular carbon, alkalinity, bicarbonate and finite volumes. Concentrations are in mM and volumes in pL.',
                     [('T_i_mM', '$T_i$'), ('A_i_mM', '$A_i$'), ('hco3_i_mM', '$[HCO_3]_i$'), ('ph_i', '$pH_i$'), ('V_i_pL', '$V_i$'), ('V_l_pL', '$V_l$')], True),
            obs_table('Lumen composition and membrane potentials. Chloride is in mM, osmolarity in mOsm, and voltage in V.',
                     [('cl_l_mM', '$[Cl]_l$'), ('ph_l', '$pH_l$'), ('osm_lumen_mOsm', '$O_l$'), ('Va_V', '$V_a$'), ('Vb_V', '$V_b$')]),
            obs_table('Signed stationary fluxes in fmol/s. NKCC and pump columns are cycle fluxes.',
                     [('N', '$N$'), ('H', '$H$'), ('E', '$E$'), ('B', '$B$'), ('J4', '$J_4$'), ('P', '$P$'), ('CaCC', '$C_a$'), ('paraCl', '$j_{p,Cl}$')], True)]
    out.append(table('Separate physiological and thermodynamic classifications at every root.',
        ['Class', '$e_4$', 'Physiology', 'Failed gate', '$A_4/(RT)$', 'Flux compatibility'],
        [[esc(r['mechanism']), num(r['expression'], 3), 'Pass' if r['full_inherited_physiological'] else 'Fail',
          esc(', '.join(str(x.get('name', x.get('gate', x))) if isinstance(x, dict) else str(x) for x in r['failed_gates']) or 'None'), num(r['forward_affinity_over_RT']), esc(r['actual_flux_thermodynamics'])]
         for r in review['equilibrium_rows']], 'lrlp{.25\\linewidth}rl', True))
    out.append(obs_table('Root, state charge and current residual diagnostics. Independent residuals use the dimensionless vector field.',
        [('residual_max', '$\\|F\\|_\\infty$'), ('condition_number', '$\\kappa(J)$'),
         ('charge_max_fmol', 'Charge (fmol)'), ('current_residual_A', 'Current (A)'),
         ('max_amount_rhs', 'Amount RHS'), ('max_volume_rhs', 'Volume RHS')], True))
    out.append(source('output/equilibrium_summary.csv') + source('output/equilibrium_review.json'))
    out += [r'\section{All state vectors, Jacobians and eigenvalues}',
            r'The following matrices are the recorded dimensionless dynamic Jacobians $F_z$, in the eleven coordinate order '
            r'$(x_{Na},x_K,x_{Cl},x_T,v,y_{Na},y_K,y_{Cl},y_T,\ell,r)$. '
            r'The physical eigenvalues are in s$^{-1}$ after division by $t_0$. '
            r'The thirteen state entries are dimensional amounts and volumes followed by the regulatory fraction, in the production order specified above. '
            r'Full precision matrices and solver messages are retained in \texttt{equilibria.json}.', r'\begin{landscape}']
    for j, r in enumerate(roots):
        if j and j % 2 == 0:
            out.append(r'\clearpage')
        out.append(r'\subsection*{' + esc(r['mechanism']) + ', expression ' + str(r['ae4_expression']) + '}')
        out.append('{\\scriptsize State: ' + ', '.join(num(v, 7) for v in r['state']) + '.\\par\n')
        out.append('Solver flag: ' + esc(r['root_success']) + '; evaluations: ' + str(r['nfev'])
                   + '; residual: ' + num(r['residual_max']) + '; physiological: ' + esc(r['physiological']) + '.\\par}\n')
        out.append('\\begingroup\\scriptsize\\setlength{\\tabcolsep}{3pt}\\centering\n\\begin{tabular}{r' + 'r'*11 + '}\\toprule\n')
        out.append('Row/column & ' + ' & '.join(str(i) for i in range(1, 12)) + r' \\ \midrule' + '\n')
        for i, row in enumerate(r['jacobian_dimensionless']):
            assert len(row) == 11
            out.append(str(i+1) + ' & ' + ' & '.join(num(v, 5) for v in row) + r' \\' + '\n')
        out.append('\\bottomrule\\end{tabular}\\par\\endgroup\n')
        eigen = sorted(r['eigenvalues_per_s'], key=lambda x: x[0], reverse=True)
        def complex_text(v: list[float]) -> str:
            return num(v[0], 8) if abs(v[1]) < 1e-13 else num(v[0], 8) + ' + ' + num(v[1], 8) + '$\\mathrm i$'
        out.append('{\\scriptsize Eigenvalues: ' + '; '.join(complex_text(v) for v in eigen) + '.\\par}\n')
    out += [r'\end{landscape}', source('output/equilibria.json')]
    out += [r'\section{Complete expression sensitivities and modal accounting}']
    deriv = review['implicit_sensitivity_two_step_validation']
    out.append(table('Implicit expression derivatives at all four verified points. Flux derivatives are per unit expression.',
        ['$e_4$', '$Q_e$', '$N_e$', '$(J_4)_e$', '$P_e$', '$H_e$', '$(pH)_e$', '$C_N$'],
        [[num(r['expression'], 3)] + [num(v, 6) for v in r['implicit_derivatives']] + [num(r['chloride_compensation_gain'])] for r in deriv], tiny=True))
    out.append(table('Independent nearby root checks for each expression derivative and difference step.',
        ['$e_4$', 'Step', 'Rule', 'Max component error', 'State error', 'Root residual'],
        [[num(r['expression'], 3), num(t['step']), esc(t['rule']), num(max(t['component_relative_errors'])),
          num(t['state_derivative_relative_error']), num(t['nearby_root_max_residual'])] for r in deriv for t in r['trials']],
        'rrlrrr', True))
    for endpoint in review['modal_analysis']:
        out.append(table('All modal projections at expression ' + str(endpoint['expression']) + '. Projections use the declared scaled coordinates.',
            ['Mode', '$\\lambda$ (1/s)', '$\\tau$ (s)', 'Input', 'Output', 'Residue', 'Stationary contribution'],
            [[str(k+1), num(m['eigenvalue_per_s'][0]), num(m['decay_time_s']), num(m['expression_input_projection'][0]),
              num(m['secretion_output_projection'][0]), num(m['transfer_residue'][0]), num(m['stationary_dQ_de_contribution'][0])]
             for k, m in enumerate(endpoint['modes'])], tiny=True))
        # Every full right vector and inverse left row remains in the source JSON.
    attr = csvrows('output/coordinate_attribution.csv')
    out.append(table('Recorded coordinate chain rule attribution. This is not a unique causal allocation to individual transporters.',
        ['$e_4$', 'Coordinate', '$dQ/de_4$ contribution'],
        [[num(r['ae4_expression'], 3), pathtex(r['coordinate']), num(r['dQ_de'])] for r in attr], 'rlr'))
    out.append(source('output/equilibrium_review.json') + source('output/coordinate_attribution.csv'))
    sens = csvrows('output/sensitivity_extended.csv')
    metrics = ['WT_Q', 'WT_pH', 'null_stationary_deficit', 'NBC_alkalinity_share', 'noNBC_observed_support_share',
               'AE4_alkalinity_demand', 'NHE_physiological_upper', 'noNBC_AE4_capacity_upper',
               'noNBC_capacity_to_current_demand', 'noNBC_capacity_gap', 'chloride_compensation_gain',
               'WT_slowest_decay', 'null_slowest_decay']
    out += [r'\section{All extended parameter sensitivity outputs}',
            r'For each of the sixteen parameters, each recorded metric is listed as an unnormalised derivative with respect to $\log p$, '
            r'a log elasticity, and the recorded difference step disagreement. For pH the unnormalised derivative is the physically useful quantity; '
            r'the source also stores a formal log elasticity of the numerical pH value. It is retained here as metadata, not reinterpreted as a hydrogen concentration elasticity. '
            r'Near zero metrics can make a log elasticity undefined; n/a preserves that fact.']
    for r in sens:
        out.append(r'\subsection*{' + esc(r['parameter'].replace('_',' ')) + '}')
        out.append('Active setting: ' + pathtex(r['inventory_path']) + ' = ' + num(r['active_value'], 10) + '. '
                   + esc(r['interval_status']) + '\n')
        out.append(table('Complete metric derivatives for ' + esc(r['parameter'].replace('_',' ')) + '.',
            ['Metric', '$dY/d\\log p$', '$d\\log Y/d\\log p$', 'Step disagreement'],
            [[pathtex(k), num(r.get(k+'_per_log_parameter')), num(r.get(k+'_log_elasticity')), num(r.get(k+'_elasticity_step_disagreement'))] for k in metrics],
            'p{.44\\linewidth}rrr', True))
        out.append('Hydrogen concentration elasticity: ' + num(r['WT_H_concentration_log_elasticity'])
                   + '; deficit percentage points per log parameter: ' + num(r['null_deficit_percentage_points_per_log_parameter'])
                   + '; maximum IFT relative error: ' + num(r['max_IFT_relative_error']) + '.\\par\n')
    out.append(source('output/sensitivity_extended.csv') + source('output/sensitivity_extended.json'))
    out += [r'\section{Complete inverse scan and both crossing sensitivity observables}']
    scan = read('inverse_detail/ordered_family_scan.json')
    out.append(table('All seventeen ordered scan points in the prescribed inverse family. Deficits are percentages; output is pL.',
        ['$\\rho$', '$b$', '$I_A$', '$I_S$', '$100D_A$', '$100D_S$', 'Gate pass'],
        [[num(r['rho'], 9), num(r['b'], 8), num(r['cumulative_adaptive_pL'], 8), num(r['cumulative_sampled_pL'], 8),
          num(100*r['adaptive_deficit'], 7), num(100*r['sampled_deficit'], 7), num(r['physiological_600s'])] for r in scan], tiny=True))
    inv = read('inverse_detail/threshold_sensitivities.json')
    out.append(table('All sixteen parameter and observable pairs for the inverse threshold.',
        ['Parameter', 'Observable', '$d\\rho/d\\log p$', '$E_\\rho$', '$E_b$', 'Validation error'],
        [[esc(r['parameter'].replace('_',' ')), esc(r['observable']), num(r['d_rho_d_log_parameter']),
          num(r['crossing_log_elasticity']), num(r['residual_recruitment_log_elasticity']), num(r['maximum_validation_relative_error'])] for r in inv],
        'p{.25\\linewidth}lrrrr', True))
    out.append(source('inverse_detail/ordered_family_scan.json') + source('inverse_detail/threshold_sensitivities.json'))
    out += [r'\section{Complete scaling and capacity metadata}',
            'The saved group file also contains loose kinetic turnover bounds and the sharper pH conditioned bounds. '
            'They are distinct calculations, not contradictory versions of one bound. '
            'The fields named NHE global bound and AE4 without NBC global bound use the broad printed turnover bound; '
            'the physiological fields use the specified pH floor. Flux bound fields retain fmol/s even though they share the dimensionless group file.']
    groups = read('output/dimensionless_groups.json')
    additional = read('output/exact_verification.json')['dimensionless_additional_groups']
    for title, data in [('All saved scale and group scalars.', groups), ('Additional scale interpretation fields.', additional)]:
        rows = []
        for k, v in data.items():
            if isinstance(v, list):
                continue
            rows.append([pathtex(k), num(v, 10) if isinstance(v, (int,float,bool)) or v is None else esc(v)])
        out.append(table(title, ['Field', 'Recorded value'], rows, 'p{.46\\linewidth}p{.46\\linewidth}', True))
    out.append(source('output/dimensionless_groups.json') + source('output/exact_verification.json'))
    (GEN / 'numerical_appendix.tex').write_text('\n'.join(out))
    COVERAGE.update(equilibrium_cases=len(roots), complete_jacobians=len(roots), jacobian_entries=sum(len(x) for r in roots for x in r['jacobian_dimensionless']),
                    eigenvalues=sum(len(r['eigenvalues_per_s']) for r in roots), state_entries=sum(len(r['state']) for r in roots),
                    expression_derivative_points=len(deriv), modal_rows=sum(len(r['modes']) for r in review['modal_analysis']),
                    coordinate_attribution_rows=len(attr), uncertain_parameters=len(sens), parameter_metric_rows=len(sens)*len(metrics),
                    ordered_inverse_points=len(scan), inverse_parameter_observable_pairs=len(inv))


def parameter_appendix() -> None:
    records = read('output/parameter_inventory.json', 'task43')
    assert len(records) == 137
    out = [r'\section{Detailed provenance for every inventory record}',
           'All 137 records are included, not only the active or sensitive subset. The full original JSON and CSV additionally preserve '
           'all raw metadata and complete code location lists. Repeated generic caveats are not evidence of an independently measured constraint. '
           'Values, units, status and the following explanatory fields are transcribed from the audit without turning logical domains into physiological uncertainty intervals.']
    keys = [('equation_role','Equation role'), ('what_constrains_it','Constraint'),
            ('what_does_not_justify_it','Limit of justification'), ('logical_domain','Logical domain'),
            ('uncertainty_status','Uncertainty status'), ('manuscript_dependence','Dependence of claims'),
            ('constraint_experiment','Discriminating measurement')]
    for k, r in enumerate(records, 1):
        out.append(r'\subsubsection*{Record ' + str(k) + ': ' + esc(r['name'].split('.')[-1].replace('_',' ')) + '}')
        out.append('{\\small ' + pathtex(r['name']) + '\\par\n')
        value = num(r['value'], 16) if isinstance(r['value'], (int,float,bool)) else esc(r['value'])
        out.append('Value: ' + value + '; units: ' + esc(r.get('units','not stated'))
                   + '; active: ' + num(r['active']) + '; classification: ' + esc(r['classification']) + '.\\par\n')
        out.append('Default: ' + esc(r.get('default_value','not recorded')) + '; differs from default: '
                   + esc(r.get('differs_from_default','not recorded')) + '.\\par\n')
        for key, label in keys:
            if key in r:
                out.append('\\textbf{' + label + '.} ' + esc(r[key]) + '\\par\n')
        out.append('\\textbf{Source.} ' + pathtex(r.get('equation_source',r.get('source_path','not stated')))
                   + '. ' + esc(r.get('literature_source','')) + '\\par}\n')
    (GEN / 'parameter_details.tex').write_text('\n'.join(out))
    COVERAGE['parameter_inventory_records'] = len(records)
    COVERAGE['parameter_active_flags'] = sum(bool(r['active']) for r in records)


def verification_appendix() -> None:
    files = ['output/exact_verification.json','output/rescaling_verification.json','output/sensitivity_verification.json',
             'output/completion_verification.json','output/clean_recomputation_verification.json',
             'output/clean_sensitivity_gate_verification.json','output/report_verification.json',
             'inverse_detail/independent_verification.json']
    out = [r'\section{Recorded verification results and their scope}',
           'The tables below expose scalar receipts from the original verification. Nested trial arrays and per file comparisons are retained in full '
           'in the original JSON. This report has not rerun the numerical replay and does not claim that reproducibility alone supplies experimental validation.']
    def scalars(data: object, prefix: str = '', depth: int = 0) -> list[list[str]]:
        ans = []
        if isinstance(data, dict):
            for k, v in data.items():
                key = (prefix + '.' if prefix else '') + k
                if isinstance(v, dict) and depth < 2:
                    ans += scalars(v, key, depth+1)
                elif not isinstance(v, (list,dict)):
                    val = num(v,10) if isinstance(v,(int,float,bool)) or v is None else (pathtex(v) if isinstance(v,str) and '/' in v and ' ' not in v else esc(v))
                    ans.append([pathtex(key), val])
        return ans
    for name in files:
        out.append(r'\subsection*{' + esc(Path(name).stem.replace('_',' ')) + '}')
        out.append(table('Scalar receipt fields.', ['Field', 'Recorded value'], scalars(read(name)), 'p{.46\\linewidth}p{.46\\linewidth}', True))
        out.append(source(name))
    controls = read('output/numerical_controls.json','task43')
    out.append(r'\section{All audited numerical controls}')
    out.append(table('Production tolerances and inversion settings. These are not bounds on biological model error.',
                     ['Control','Recorded value'], scalars(controls), 'p{.49\\linewidth}p{.43\\linewidth}', True))
    out.append(source('output/numerical_controls.json','task43'))
    (GEN/'verification_appendix.tex').write_text('\n'.join(out))


def catalogue(manifest: list[dict[str,object]]) -> None:
    out = [r'\section{Complete electronic evidence catalogue}',
           'Every file below is included in the report folder as an unchanged copy from the pinned scientific commit. '
           'The full SHA256, source path, source commit and byte count are recorded in '+pathtex('SOURCE_MANIFEST.json')+'. '
           'No cache, missing historical seed or new experimental observation is inferred from this catalogue.']
    for group in PINS:
        rows = [[pathtex(str(r['report_path']).split(group+'/',1)[1]), num(r['size_bytes'],10)]
                for r in manifest if str(r['report_path']).startswith('evidence/'+group+'/')]
        out.append(table('Complete '+group+' snapshot.', ['File within snapshot','Bytes'], rows, 'p{.78\\linewidth}r', True))
    (GEN/'evidence_catalogue.tex').write_text('\n'.join(out))


def main() -> None:
    parser=argparse.ArgumentParser()
    parser.add_argument('--source43',help='Optional local exact Task 43 directory for offline assembly')
    parser.add_argument('--source44',help='Optional local exact Task 44 directory for offline assembly')
    args=parser.parse_args()
    GEN.mkdir(parents=True,exist_ok=True)
    manifest=[]
    for name, supplied in [('task43',args.source43),('task44',args.source44)]:
        if not supplied and (EVID/name/'report.tex').exists() and not (ROOT/'.git').exists():
            supplied=str(EVID/name)
        manifest += snapshot(name,supplied)
    dump(HERE/'SOURCE_MANIFEST.json',dict(pinned_inputs=PINS,files=manifest,scope='Unchanged copies of complete published Task 43 and Task 44 task directories.'))
    body('task44','m44');body('task43','p43');bibliography()
    numerical_appendix();parameter_appendix();verification_appendix();catalogue(manifest)
    COVERAGE['source_snapshot_files']=len(manifest)
    COVERAGE['scientific_recomputation_performed_by_report_assembly']=False
    dump(GEN/'coverage.json',COVERAGE)
    text='# Coverage of the consolidated technical report\n\n'
    text+='The full final source exposition is retained, and expanded tables are generated only from saved results.\n\n'
    for k,v in COVERAGE.items():
        text+='- '+k+': '+json.dumps(v,ensure_ascii=False)+'.\n'
    text+='\nAll original machine readable records, scripts, TeX fragments, verification receipts and the two original PDFs are preserved under `evidence/`, with byte hashes in `SOURCE_MANIFEST.json`.\n'
    text+='\nThe main report prints all state vectors, Jacobian matrices and eigenvalues. Full precision modal vectors, every nearby trial, diagnostic residual and original verification exclusion remain in the electronic evidence. No new scientific simulation is performed by report assembly.\n'
    (HERE/'COVERAGE.md').write_text(text)
    print(json.dumps({k:v for k,v in COVERAGE.items() if not isinstance(v,dict)},sort_keys=True))


if __name__=='__main__':
    main()
