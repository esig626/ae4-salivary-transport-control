"""Analytical tests of paired indexing and affine resource constraints, not Task46G replay."""
import json,sys
from pathlib import Path
from lp_tools import *
m=FluxModel((FluxMetabolite('x'),),(FluxReaction('in',(StoichiometricTerm('x',1),),0,2),FluxReaction('out',(StoichiometricTerm('x',-1),),0,10)),LinearObjective('maximise',(ObjectiveTerm('out',1),)))
c=upper(m,{'in':1,'out':1},3,'sum');r=solve(c);assert abs(r['objective']-1.5)<1e-10
r2=solve(fix(m,{'out':1}));assert r2['fluxes']['in']==1
p=paired(m,m);p=objective(p,{'WT:out':1,'KO:out':1});r3=solve(p);assert r3['objective']==4
p=upper(p,{'KO:out':1,'WT:out':-1},-1,'difference');r4=solve(objective(p,{'KO:out':1}));assert r4['objective']==1
for result in (r,r2,r3,r4):assert max(result['diagnostics'].values())<1e-8
out={'passed':True,'analytical_tests':4,'expected':[1.5,1,4,1],'observed':[r['objective'],r2['objective'],r3['objective'],r4['objective']],'scope':'Task-specific block indexing, fixed bounds, affine cap, signed cross-block inequality; unchanged vendor'}
Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
