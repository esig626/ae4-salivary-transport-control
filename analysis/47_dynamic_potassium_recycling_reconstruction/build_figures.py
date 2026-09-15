"""Article figures from frozen CSV archives; no model import or simulation."""
from pathlib import Path
import csv
import gzip
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent
OUT=HERE/'figures'
OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':9,'axes.titlesize':10,'axes.labelsize':9,
                     'svg.fonttype':'none','svg.hashsalt':'task47-fixed-results',
                     'axes.spines.top':False,'axes.spines.right':False})

def rows(case):
    with gzip.open(HERE/'output'/case/'trajectory.csv.gz','rt') as f:
        return list(csv.DictReader(f))

def save(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=170,facecolor='white')
    fig.savefig(OUT/(name+'.svg'),metadata={'Date':None},facecolor='white')
    plt.close(fig)

def main():
    data={case:rows(case) for case in ['wt','ae4_null']}
    colours={'wt':'#224E70','ae4_null':'#B45F24'}
    labels={'wt':'WT','ae4_null':'AE4 null'}
    fig,axes=plt.subplots(3,2,figsize=(10.2,8.6))
    panels=[('q_out_pL_s',1e3,'Water output (fL/s)','A  Secretion'),
            ('nkcc1_cycle_inward_fmol_s',1,'NKCC cycles (fmol/s)','B  Chloride rescue recruits NKCC1'),
            ('k_total_efflux_fmol_s',1,'K efflux (fmol/s)','C  K channel recycling increases'),
            ('k_i_mM',1,'Intracellular K (mM)','D  K concentration changes modestly'),
            ('k_drive_basolateral_V',1e3,'Basolateral K driving voltage (mV)','E  A small extra driving voltage suffices'),
            ('pump_k_conditioned_utilisation',100,'Pump / available kinetic limit (%)','F  Pump saturation is not reached')]
    for ax,(field,scale,label,title) in zip(axes.flat,panels):
        for case,rr in data.items():
            if field=='k_drive_basolateral_V':
                rr=[r for r in rr if float(r['time_s'])>0]
            ax.plot([float(r['time_s']) for r in rr],[scale*float(r[field]) for r in rr],
                    color=colours[case],label=labels[case],lw=1.8)
        ax.set(xlim=(0,600),xlabel='Time (s)',ylabel=label,title=title)
        ax.grid(axis='y',alpha=.16)
        ax.set_xticks([0,150,300,450,600])
        if field=='pump_k_conditioned_utilisation':
            ax.set_ylim(0,105);ax.axhline(100,color='#777777',ls=':',lw=1)
    fig.legend(*axes[0,0].get_legend_handles_labels(),loc='upper center',ncol=2,
               frameon=False,bbox_to_anchor=(.5,.963))
    fig.suptitle('Frozen model: additional K recycling does not exhaust the pump',fontsize=13,y=.993)
    fig.subplots_adjust(left=.085,right=.98,bottom=.09,top=.90,wspace=.29,hspace=.53)
    fig.text(.085,.013,'Identical parameters and WT initial state. Pump limit is model based, not a measured physiological capacity.\nPanel E starts at the stimulated right limit; the instantaneous resting voltage jump is omitted.',fontsize=8,color='#444444')
    save(fig,'dynamic_recycling')

    summaries={case:json.loads((HERE/'output'/case/'summary.json').read_text())['integrals']['60_600'] for case in ['wt','ae4_null']}
    w,k=summaries['wt'],summaries['ae4_null']
    delta=lambda f:k[f]-w[f]
    kvals=[delta('k_nkcc_influx_fmol_s'),-delta('k_ae4_export_fmol_s'),delta('k_pump_influx_fmol_s'),-delta('k_total_efflux_fmol_s'),-delta('k_storage_fmol_s')]
    nvals=[delta('na_nkcc_entry_fmol_s'),-delta('na_ae4_export_fmol_s'),delta('na_nbc_entry_fmol_s'),delta('na_nhe_entry_fmol_s'),-delta('na_pump_clearance_fmol_s'),-delta('na_storage_fmol_s')]
    fig,axes=plt.subplots(2,1,figsize=(10.2,7.8))
    specs=[(kvals,['Extra NKCC\nentry','AE4 export\nremoved','Extra pump\nK entry','Extra channel\nefflux','Extra K\nstorage'],'A  Potassium: more channel efflux and more storage'),
           (nvals,['Extra NKCC\nentry','AE4 export\nremoved','Reduced NBC\nentry','Reduced NHE1\nentry','Extra pump\nNa clearance','Extra Na\nstorage'],'B  Sodium: loss of other entry routes offsets most of the added burden')]
    for ax,(values,names,title) in zip(axes,specs):
        assert abs(sum(values))<1e-8
        bars=ax.bar(np.arange(len(values)),values,color=['#B45F24' if v>0 else '#224E70' for v in values],width=.58)
        ax.axhline(0,color='#666666',lw=.7)
        ax.set_xticks(np.arange(len(values)),names)
        ax.set_ylabel('Null minus WT balance term (fmol)');ax.set_title(title,loc='left',pad=12)
        ax.set_ylim(-39,29);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
        for bar,v in zip(bars,values):
            ax.text(bar.get_x()+bar.get_width()/2,v+(1 if v>0 else -1),f'{v:+.2f}',ha='center',va='bottom' if v>0 else 'top',fontsize=9)
    fig.suptitle('Where the extra ion burden goes',fontsize=13,y=.985)
    fig.subplots_adjust(left=.09,right=.98,top=.91,bottom=.10,hspace=.42)
    fig.text(.09,.026,'Integrated model differences, 60 to 600 s. Positive terms add intracellular burden; negative terms offset it.\nStorage is placed on the clearance side of each balance. Each panel sums to zero.',fontsize=8,color='#444444')
    save(fig,'ion_balance_attribution')
    print('Built two figures from saved predictions; no model execution.')

if __name__=='__main__':main()
