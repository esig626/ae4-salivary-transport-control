"""One conserved, conditional transport network. No kinetic model is imported."""
from __future__ import annotations
import csv, hashlib, json, sys
from fractions import Fraction as F
from pathlib import Path
HERE=Path(__file__).resolve().parent
AREA=HERE.parent
sys.path.insert(0,str(AREA/'vendor'))
from carbonscope_flux_core import FluxModel, FluxMetabolite, FluxReaction, StoichiometricTerm, LinearObjective, ObjectiveTerm
# Coordinates are Na,K,Cl,TIC,TA in cell then lumen. External bath is a reservoir.
BIO_ROWS=['Na_i','K_i','Cl_i','TIC_i','TA_i','Na_l','K_l','Cl_l','TIC_l','TA_l']
REACTIONS=[
 ('N','NKCC1 bath to cell',[1,1,2,0,0],[0,0,0,0,0],False),
 ('A','AE4 equal routing chloride entry',[-.5,-.5,1,-2,-2],[0,0,0,0,0],False),
 ('E','AE2 chloride entry',[0,0,1,-1,-1],[0,0,0,0,0],True),
 ('H','NHE1 proton extrusion',[1,0,0,0,1],[0,0,0,0,0],False),
 ('B','NBC one Na and two bicarbonates entry',[1,0,0,2,2],[0,0,0,0,0],False),
 ('Pb','Basolateral Na/K ATPase',[-3,2,0,0,0],[0,0,0,0,0],False),
 ('Pa','Apical Na/K ATPase',[-3,2,0,0,0],[3,-2,0,0,0],False),
 ('Kb','Basolateral K efflux',[0,-1,0,0,0],[0,0,0,0,0],False),
 ('Ka','Apical K efflux',[0,-1,0,0,0],[0,1,0,0,0],False),
 ('J','Apical chloride export',[0,0,-1,0,0],[0,0,1,0,0],False),
 ('Db','CO2 bath to cell',[0,0,0,1,0],[0,0,0,0,0],True),
 ('Da','CO2 lumen to cell',[0,0,0,1,0],[0,0,0,-1,0],True),
 ('pNa','Paracellular Na lumen to bath',[0,0,0,0,0],[-1,0,0,0,0],True),
 ('pK','Paracellular K lumen to bath',[0,0,0,0,0],[0,-1,0,0,0],True),
 ('pCl','Paracellular Cl lumen to bath',[0,0,0,0,0],[0,0,-1,0,0],True),
 ('pHCO3','Paracellular bicarbonate lumen to bath',[0,0,0,0,0],[0,0,0,-1,-1],True),
 ('oNa','Luminal sodium outflow',[0,0,0,0,0],[-1,0,0,0,0],False),
 ('oK','Luminal potassium outflow',[0,0,0,0,0],[0,-1,0,0,0],False),
 ('oCl','Luminal chloride outflow',[0,0,0,0,0],[0,0,-1,0,0],False),
 ('oHCO3','Luminal bicarbonate outflow',[0,0,0,0,0],[0,0,0,-1,-1],False),
 ('oCO2','Luminal dissolved CO2 outflow',[0,0,0,0,0],[0,0,0,-1,0],False),
]
IDS=[r[0] for r in REACTIONS]

def rref(a):
    a=[[F(x) for x in row] for row in a]; piv=[]; r=0
    for c in range(len(a[0])):
        k=next((k for k in range(r,len(a)) if a[k][c]),None)
        if k is None: continue
        a[r],a[k]=a[k],a[r]; z=a[r][c]; a[r]=[x/z for x in a[r]]
        for k in range(len(a)):
            if k!=r and a[k][c]:
                z=a[k][c]; a[k]=[x-z*y for x,y in zip(a[k],a[r])]
        piv.append(c); r+=1
        if r==len(a): break
    return a,piv

def matrix():
    s=[[F((r[2]+r[3])[i]) for r in REACTIONS] for i in range(10)]
    # Cell/lumen net charge rows are audits, not individual transporter neutrality.
    for offset in (0,5): s.append([s[offset][j]+s[offset+1][j]-s[offset+2][j]-s[offset+4][j] for j in range(len(IDS))])
    neutral={'oNa':1,'oK':1,'oCl':-1,'oHCO3':-1}
    s.append([F(neutral.get(k,0)) for k in IDS])
    return s
ROWS=BIO_ROWS+['charge_i','charge_l','neutral_outflow']
S=matrix()
# Pivot columns of transposed S identify independent ORIGINAL biological rows.
_,INDEPENDENT=rref(list(map(list,zip(*S))))

def build(genotype='WT',M=100.0,objective=None,direction='max',fixed=None):
    if genotype not in ('WT','KO'): raise ValueError(genotype)
    if M<=0: raise ValueError('M must be positive')
    rx=[]; fixed=fixed or {}
    for j,(name,meaning,cell,lumen,reversible) in enumerate(REACTIONS):
        low=-M if reversible else 0.; high=1. if name=='J' else M
        if genotype=='KO' and name=='A': low=high=0.
        if name in fixed:
            value=float(fixed[name]); assert low<=value<=high
            low=high=value
        rx.append(FluxReaction(name,tuple(StoichiometricTerm(ROWS[i],float(S[i][j])) for i in INDEPENDENT if S[i][j]),low,high))
    obj={'J':1.} if objective is None else objective
    return FluxModel(tuple(FluxMetabolite(ROWS[i]) for i in INDEPENDENT),tuple(rx),LinearObjective(direction,tuple(ObjectiveTerm(k,v) for k,v in obj.items())))

def identity():
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

def validate(v,genotype='WT',M=100.):
    model=build(genotype,M); values=[float(v[k]) for k in IDS]
    residual=max(abs(sum(float(a)*x for a,x in zip(row,values))) for row in S)
    bound=max([0.]+[max(r.lower_bound-x,x-r.upper_bound,0.) for r,x in zip(model.reactions,values)])
    return {'max_full_balance_residual':residual,'max_bound_violation':bound,'valid':max(residual,bound)<1e-7}

def write_tables():
    def table(name,headers,rows):
        with (HERE/name).open('w',newline='') as f:
            w=csv.writer(f);w.writerow(headers);w.writerows(rows)
    table('stoichiometric_matrix.csv',['balance']+IDS,[[k]+[str(x) for x in row] for k,row in zip(ROWS,S)])
    meanings=['cell sodium','cell potassium','cell chloride','cell total inorganic carbon','cell total alkalinity','lumen sodium','lumen potassium','lumen chloride','lumen total inorganic carbon','lumen total alkalinity','cell charge audit Na+K-Cl-TA','lumen charge audit Na+K-Cl-TA','charge neutral luminal outflow']
    table('balance_ledger.csv',['id','meaning','retained_in_solver','storage_assumption'],[[k,m,i in INDEPENDENT,'zero sustained amount storage' if i<10 else 'exact charge constraint'] for i,(k,m) in enumerate(zip(ROWS,meanings))])
    table('reaction_ledger.csv',['id','meaning','cell_Na_K_Cl_TIC_TA','lumen_Na_K_Cl_TIC_TA','signed','source'],[[n,m,json.dumps(c),json.dumps(l),rev,'output/46h/evidence/stoichiometry.json; minimal bicarbonate/CO2 outflow split in CBM_PROTOCOL.md'] for n,m,c,l,rev in REACTIONS])
    table('bounds_provenance.csv',['id','lower_WT','upper_WT','lower_KO','upper_KO','classification','interpretation'],[[r.reaction_id,r.lower_bound,r.upper_bound,k.lower_bound,k.upper_bound,'WT computational big M; conditional screening sign; KO algebraic genotype zero' if r.reaction_id=='A' else ('computational normalisation' if r.reaction_id=='J' else 'computational big M; conditional screening sign'),'M=100, predeclared sensitivity M=1000; not measured capacity; source context output/46h/evidence/bounds.json'] for r,k in zip(build().reactions,build('KO').reactions)])
    rank=len(INDEPENDENT)
    data={'builder_sha256':identity(),'reactions':IDS,'full_rows':ROWS,'independent_rows':[ROWS[i] for i in INDEPENDENT],'full_rank':rank,'independent_rank':rank,'nullity':len(IDS)-rank,'charge_rows_redundant':all(i not in INDEPENDENT for i in (10,11)),'neutral_outflow_retained':12 in INDEPENDENT,'M':100,'scientific_optimisations':0}
    (HERE/'network_identity.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data))
if __name__=='__main__': write_tables()
