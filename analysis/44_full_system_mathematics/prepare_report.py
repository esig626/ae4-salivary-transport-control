"""Generate report tables and figures solely from the new saved analysis."""
from pathlib import Path
import csv,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
D=Path(__file__).resolve().parent/'output'
def read(n):return json.loads((D/n).read_text())
def rows(n):return list(csv.DictReader((D/n).open()))
def table(name,header,body,align):
    text='\\begin{tabular}{'+align+'}\n\\toprule\n'+' & '.join(header)+' \\\\\n\\midrule\n'
    text+=''.join(' & '.join(map(str,r))+' \\\\\n' for r in body)
    text+='\\bottomrule\n\\end{tabular}\n';(D/name).write_text(text)
g=read('dimensionless_groups.json')
table('scale_table.tex',['Quantity','Value','Units'],[
 ['$C_0$',f"{g['C0_mM']:.0f}",'mM'],['$V_0$',f"{g['V0_pL']:.8f}",'pL'],['$L_0$',f"{g['L0_pL']:.8f}",'pL'],
 ['$J_0$',f"{g['J0_fmol_s']:.8f}",'fmol/s'],['$t_0$',f"{g['t0_s']:.5f}",'s'],
 ['$\\epsilon_L$',f"{g['epsilon_lumen']:.6f}",'1'],['$\\theta_4$',f"{g['theta_ae4']:.6f}",'1'],
 ['$G_B/J_0$',f"{g['nbc_capacity_ratio']:.6f}",'1'],['$G_2/J_0$',f"{g['ae2_capacity_ratio']:.6f}",'1'],
 ['$L_aC_0^2/J_0$',f"{g['water_apical']:.3f}",'1'],['$L_bC_0^2/J_0$',f"{g['water_basolateral']:.3f}",'1'],
 ['$L_tC_0^2/J_0$',f"{g['water_paracellular']:.3f}",'1']], 'lrl')
eq=read('equilibria.json')['roots'];wt=[r for r in eq if r['ae4_expression']==1.]
table('class_table.tex',['Class','$\\mathrm{pH}_i^*$','$[\\mathrm{Na}]_i^*$','$t_{\\mathrm{slow}}$ (s)','Failed criteria'],[
 [r['mechanism'],f"{r['observables']['ph_i']:.4f}",f"{r['observables']['na_i_mM']:.3f}",f"{r['slowest_decay_s']:.2f}",', '.join({'na_i_mM':'sodium','ph_i':'pH'}.get(k,k) for k in r['failed_gates']).replace('_','\\_') or 'None'] for r in wt], 'lrrrl')
p=rows('source_projections.csv');table('projection_table.tex',['Class','$\\nu_{Na}$','$\\nu_K$','$\\nu_{Cl}$','$\\nu_T$','$\\nu_A$','$\\kappa$','$\\mathcal F_4/(RT)$'],[
 [r['mechanism']]+[f"{float(r[k]):g}" for k in ['nu_Na','nu_K','nu_Cl','nu_T','nu_A','kappa']]+[f"{float(r['affinity_reference']):.3f}"] for r in p], 'lrrrrrrr')
s=rows('parameter_sensitivity.csv');table('sensitivity_table.tex',['Parameter','$\\partial\\log Q^*/\\partial\\log p$','$\\partial\\mathrm{pH}_i^*/\\partial\\log p$'],[
 [r['parameter'].replace('_',' '),f"{float(r['Q_log_elasticity']):.6f}",f"{float(r['pH_per_log_parameter']):.6f}"] for r in s], 'lrr')
base=next(r for r in eq if r['mechanism']=='C0' and r['ae4_expression']==1.)
null=next(r for r in eq if r['mechanism']=='C0' and r['ae4_expression']==0.)
sorted_lam=lambda r:sorted(r['eigenvalues_per_s'],key=lambda x:x[0],reverse=True)
fmt=lambda x: f"${x[0]:.6g}{x[1]:+.5g}\\mathrm{{i}}$" if abs(x[1])>1e-10 else f"${x[0]:.6g}$"
table('spectrum_table.tex',['Mode','WT (1/s)','Null (1/s)'],[[i+1,fmt(a),fmt(b)] for i,(a,b) in enumerate(zip(sorted_lam(base),sorted_lam(null)))], 'rrr')
nb=rows('nbc_continuation.csv');table('nbc_table.tex',['NBC capacity fraction','$J_4^*$ (fmol/s)','$J_B^*$ (fmol/s)','$\\mathrm{pH}_i^*$'],[
 [f"{float(r['NBC_fraction']):.2f}",f"{float(r['J4']):.6f}",f"{float(r['B']):.6f}",f"{float(r['ph_i']):.4f}"] for r in nb], 'rrrr')
def save(fig,name):
    fig.tight_layout();fig.savefig(D/(name+'.pdf'));fig.savefig(D/(name+'.png'),dpi=140);plt.close(fig)
c0=sorted([r for r in eq if r['mechanism']=='C0'],key=lambda r:r['ae4_expression'])
fig,ax=plt.subplots(figsize=(6.8,3.9));ax.plot([r['ae4_expression'] for r in c0],[r['observables']['q_pL_s']/base['observables']['q_pL_s'] for r in c0],marker='o')
ax.set(xlabel='Relative AE4 expression',ylabel='Stationary flow / WT stationary flow',title='Stationary secretion on the full neutral charge manifold');save(fig,'equilibrium_flow')
fig,ax=plt.subplots(figsize=(6.8,3.9));ax.plot([r['ae4_expression'] for r in c0],[r['slowest_decay_s'] for r in c0],marker='o')
ax.set(xlabel='Relative AE4 expression',ylabel='Slowest local decay time (s)',title='Slow modes are already present');save(fig,'slow_modes')
lr=rows('local_step_response.csv');fig,ax=plt.subplots(figsize=(6.8,3.9));ax.plot([float(r['time_s']) for r in lr],[float(r['linear_dQ_de'])*1e6 for r in lr],marker='o',label='Linear response')
ax.axhline(float(lr[0]['asymptotic_dQ_de'])*1e6,linestyle=':',label='Stationary sensitivity');ax.set(xlabel='Time after a small expression step (s)',ylabel='Flow sensitivity (10$^{-6}$ pL/s)',title='Transient sensitivity need not resemble stationary sensitivity');ax.legend();save(fig,'linear_response')
a=sorted([r for r in read('inverse_samples.json') if not r['tight']],key=lambda r:r['rho'])
fig,ax=plt.subplots(figsize=(6.8,3.9));ax.plot([r['rho'] for r in a],[100*r['adaptive_deficit'] for r in a],marker='o');ax.axhspan(30.3,39.7,alpha=.12,label='Reported phenotype comparison band')
ax.set(xlabel='AE4 dependent recruitment fraction, rho',ylabel='Cumulative null deficit at 600 s (%)',title='Inverse CaCC family: a numerical crossing, not a global bound');ax.legend(fontsize=8);save(fig,'inverse_threshold')
l=rows('inverse_long_time.csv');fig,ax=plt.subplots(figsize=(6.8,3.9));ax.plot([float(r['time_s'])/60 for r in l],[float(r['ph_i']) for r in l]);ax.axhline(7.3,linestyle=':',label='Inherited upper pH limit');ax.axvline(10,linestyle='--',label='Original 600 s horizon');ax.set(xlabel='Continued stimulation (min)',ylabel='Intracellular pH',title='The inverse construction is not valid indefinitely');ax.legend(fontsize=8);save(fig,'inverse_ph')
print('Report tables and five figures generated from saved values.')
