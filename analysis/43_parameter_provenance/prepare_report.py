"""Typeset a readable inventory from the exact active-object audit."""
from pathlib import Path
import json
D=Path(__file__).resolve().parent/'output'
a=json.loads((D/'parameter_inventory.json').read_text())
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
print('Full inventory table generated:',len(a),'records')
