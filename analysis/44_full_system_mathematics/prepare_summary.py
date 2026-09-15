"""Report tables from the published numerical checkpoint; no model solves."""
import csv,json
from pathlib import Path
P=Path(__file__).resolve().parent;O=P/'output'
g=json.loads((O/'dimensionless_groups.json').read_text())
a=json.loads((O/'exact_verification.json').read_text())['dimensionless_additional_groups']
rows=[
(r'$n_{b,i}/(C_0V_0)$','Cell buffer amount',a['cell_buffer_scaled']),
(r'$n_{b,l}/(C_0L_0)$','Lumen buffer amount',a['lumen_buffer_scaled']),
(r'$n_{\rm fixed}/(C_0V_0)$','Fixed cell charge',a['fixed_cell_charge_scaled']),
(r'$n_{\rm other}/(C_0V_0)$','Other cell osmoles',a['cell_other_impermeant_scaled']),
(r'$K_1/C_0$','First carbonate dissociation',a['Ka1_over_C0']),
(r'$K_2/C_0$','Second carbonate dissociation',a['Ka2_over_C0']),
(r'$K_w/C_0^2$','Water dissociation',a['Kw_over_C0_squared']),
(r'$K_{b,i}/C_0$','Cell buffer dissociation',a['cell_buffer_K_over_C0']),
(r'$K_{b,l}/C_0$','Lumen buffer dissociation',a['lumen_buffer_K_over_C0']),
(r'$P_{{\rm CO}_2,b}C_0/J_0$','Basolateral carbon dioxide exchange',a['co2_basolateral']),
(r'$P_{{\rm CO}_2,a}C_0/J_0$','Apical carbon dioxide exchange',a['co2_apical']),
(r'$J_4(0^+)/J_0$','Initial AE4 to NKCC cycle ratio',a['AE4_to_NKCC_cycle_initial']),
(r'$J_4(0^+)/(2J_0)$','Initial chloride loading ratio',a['AE4_to_NKCC_chloride_initial']),
(r'$t_{\rm carrier}/t_0$','Inherited carrier relaxation ratio',a['carrier_QSS_time_ratio_initial'])]
lines=[r'\begin{tabular}{llr}',r'\toprule',r'Group & Meaning & Value \\',r'\midrule']
for symbol,meaning,value in rows:lines.append(symbol+' & '+meaning+' & '+f'{value:.6g}'+r' \\')
lines += [r'\bottomrule',r'\end{tabular}']
(O/'additional_groups_table.tex').write_text('\n'.join(lines)+'\n')
print('Additional dimensionless groups table generated from verified values.')
