"""Three predeclared paired/protocol diagnostics. No numerical NKCC tolerance."""
import argparse,json,subprocess
from pathlib import Path
import build_cbm as b
from analyse_region import compact_witness
from lp_tools import *

def run(which):
    out={'diagnostic':which,'input_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=b.AREA.parents[1],text=True).strip(),'builder_sha256':b.identity(),'numerical_NKCC_tolerance':None,'normal_stimulation_flux_equality_imposed':False}
    if which=='maximum':
        p=paired(fix(model(),{'J':1}),model('KO'));out['result']=solve(p)
        out['normalised_KO_maximum']=out['result']['objective']
        out['interpretation']='No quantitative cross-genotype stimulated NKCC relation follows from the distinct assay; paired region retains only shared capacities and WT J=1.'
    elif which=='replacement':
        w=compact_witness('WT');p=paired(fix(model(),w),fix(model('KO'),{'J':1}))
        out['WT_reference']=w;out['KO_NKCC_minimum']=solve(objective(p,{'KO:N':1},'min'));out['KO_NKCC_maximum']=solve(objective(p,{'KO:N':1},'max'))
        out['difference_range']=[out[k]['objective']-w['N'] for k in ('KO_NKCC_minimum','KO_NKCC_maximum')]
        out['missing_AE4_chloride_replacement_range']=[2*x/w['A'] for x in out['difference_range']]
        out['interpretation']='For this explicit illustrative WT witness, full output allows either no NKCC chloride replacement or complete replacement; negative values mean NKCC falls while AE2 changes. This is not an inferred WT partition or a phenotype fit.'
    else:
        blocked={'A':0,'E':0,'B':0,'J':0}
        p=paired(fix(model(),blocked),fix(model('KO'),blocked))
        out['steady_blocked_assay_NKCC_maximum']=solve(objective(p,{'WT:N':1,'KO:N':1}))
        out['assay_scope']='Idealised bicarbonate routes and CaCC blocked; no conversion of carbonic anhydrase inhibition to a numerical rate. Initial uptake requires storage.'
        out['initial_uptake_identity']={'NKCC_contribution_to_dNa_i':'N','NKCC_contribution_to_dK_i':'N','NKCC_contribution_to_dCl_i':'2N','NKCC_net_charge_storage':'N+N-2N=0'}
        out['interpretation']='Steady chloride balance would force N=0 in this blocked protocol. That mathematical obstruction shows why initial uptake cannot be represented by the normal steady CBM. Storage identity is assay-matched bookkeeping, not a new kinetic model or a quantitative fit. Equal assay uptake is compatible with shared capacities but neither predicted nor transferable to physiological stimulation.'
    return out
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--diagnostic',choices=['maximum','replacement','assay'],required=True);p.add_argument('--output',required=True);a=p.parse_args();d=run(a.diagnostic);Path(a.output).write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'diagnostic':a.diagnostic,'output':a.output}))
