"""Task46-specific paired and affine constraints above the frozen solver core."""
from dataclasses import replace
import build_cbm as b
from carbonscope_flux_core import FluxModel,FluxMetabolite,FluxReaction,StoichiometricTerm,LinearObjective,ObjectiveTerm
from carbonscope_flux_core.solver import run_highs_fba

def model(genotype='WT',M=100,objective=None,direction='max',fixed=None):
    return b.build(genotype,M,objective,'maximise' if direction=='max' else 'minimise',fixed)

def objective(m,coeff,direction='max'):
    return replace(m,objective=LinearObjective('maximise' if direction=='max' else 'minimise',tuple(ObjectiveTerm(k,v) for k,v in coeff.items())))

def fix(m,values):
    rx=[]
    for r in m.reactions:
        if r.reaction_id in values:
            v=float(values[r.reaction_id]); assert r.lower_bound-1e-9<=v<=r.upper_bound+1e-9
            r=replace(r,lower_bound=v,upper_bound=v)
        rx.append(r)
    return replace(m,reactions=tuple(rx))

def paired(wt=None,ko=None):
    models=[('WT:',model() if wt is None else wt),('KO:',model('KO') if ko is None else ko)]
    ms=[];rs=[]
    for prefix,m in models:
        ms.extend(FluxMetabolite(prefix+x.metabolite_id,x.steady_state_balanced) for x in m.metabolites)
        rs.extend(FluxReaction(prefix+r.reaction_id,tuple(StoichiometricTerm(prefix+t.metabolite_id,t.coefficient) for t in r.stoichiometric_terms),r.lower_bound,r.upper_bound) for r in m.reactions)
    return FluxModel(tuple(ms),tuple(rs),LinearObjective('maximise',(ObjectiveTerm('KO:J',1),)))

def upper(m,coeff,cap,label):
    """Encode sum(a*v)<=cap exactly as one pseudo balance with nonnegative slack."""
    row='constraint:'+label; slack='slack:'+label; unit='unit:'+label
    ids={r.reaction_id for r in m.reactions};assert not {slack,unit}&ids and set(coeff)<=ids
    # A finite algebraic slack bound covers every vector in the original box.
    limit=abs(cap)+sum(abs(coeff.get(r.reaction_id,0))*max(abs(r.lower_bound),abs(r.upper_bound)) for r in m.reactions)
    rs=[replace(r,stoichiometric_terms=r.stoichiometric_terms+((StoichiometricTerm(row,coeff[r.reaction_id]),) if r.reaction_id in coeff else ())) for r in m.reactions]
    rs += [FluxReaction(slack,(StoichiometricTerm(row,1),),0,max(1.,limit)),FluxReaction(unit,(StoichiometricTerm(row,-cap),),1,1)]
    return replace(m,metabolites=m.metabolites+(FluxMetabolite(row),),reactions=tuple(rs))

def solve(m):
    r=run_highs_fba(m)
    return {'objective':r.objective_value,'fluxes':r.fluxes.to_dict(),'diagnostics':vars(r.diagnostics) if hasattr(r.diagnostics,'__dict__') else {k:getattr(r.diagnostics,k) for k in r.diagnostics.__slots__}}
