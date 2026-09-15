"""Typeset a readable inventory from the exact active-object audit."""
from pathlib import Path
import json
D=Path(__file__).resolve().parent/'output'
a=json.loads((D/'parameter_inventory.json').read_text())
(D/'audit_counts.tex').write_text('\\newcommand{\\InventoryCount}{'+str(len(a))+'}\n\\newcommand{\\ActiveCount}{'+str(sum(r['active'] for r in a))+'}\n')
def esc(s):
    return str(s).replace('\\','\\textbackslash{}').replace('&','\\&').replace('%','\\%').replace('#','\\#').replace('_','\\_').replace('$','\\$')
def val(v):
    if isinstance(v,bool):return str(v)
    if isinstance(v,(int,float)):return f'{v:.9g}'
    return '\\path{'+str(v)+'}'
lines=[r'\setlength{\tabcolsep}{4pt}',r'\small',r'\begin{longtable}{p{.38\textwidth}p{.17\textwidth}p{.11\textwidth}p{.25\textwidth}}',
       r'\toprule Setting & Value & Used & Classification \\ \midrule\endfirsthead',
       r'\toprule Setting & Value & Used & Classification \\ \midrule\endhead',r'\bottomrule\endfoot']
for r in a:
    lines.append('\\path{'+r['name']+'} & '+val(r['value'])+' & '+('Yes' if r['active'] else 'No')+' & '+esc(r['classification'].replace('_',' '))+r' \\')
lines+=['\\end{longtable}','\\normalsize']
(D/'inventory_table.tex').write_text('\n'.join(lines)+'\n')
crosswalk=D/'task44_sensitivity_crosswalk.json'
if crosswalk.exists():
    cross=json.loads(crosswalk.read_text())
    with (D/'audit_counts.tex').open('a') as f:
        f.write('\\newcommand{\\SensitivitySourceCommit}{'+cross['source_commit']+'}\n')
    lines=[r'\begin{table}[tbp]',r'\centering\small',r'\begin{tabular}{p{.36\textwidth}rrrr}',
       r'\toprule Setting & $E_Q$ & $d\mathrm{pH}/d\log p$ & $E_D$ & $E_{\tau}$ \\ \midrule']
    for r in cross['parameters']:
        z=r['results'];values=[z[k] for k in ['WT_Q_log_elasticity','WT_pH_per_log_parameter','null_stationary_deficit_log_elasticity','WT_slowest_decay_log_elasticity']]
        lines.append(esc(r['parameter'].replace('_',' '))+' & '+' & '.join(f'{v:.5g}' for v in values)+r' \\')
    lines+=[r'\bottomrule\end{tabular}',r'\caption{Local stationary sensitivities imported from Task 44. No measured uncertainty intervals are implied.}',r'\end{table}']
    (D/'sensitivity_table.tex').write_text('\n'.join(lines)+'\n')
print('Full inventory table generated:',len(a),'records')
