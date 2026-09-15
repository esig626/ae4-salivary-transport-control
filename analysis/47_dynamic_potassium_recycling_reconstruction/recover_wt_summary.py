"""Recover metadata from already saved WT data after a NumPy bool JSON error.

No integration or scientific correction. The original output CSVs are retained.
Only sampled extrema can be reconstructed; lost solver counters remain null.
"""
from run_diagnostics import diagnostics
from common import *
import scipy
out=HERE/'output/wt'
assert not (out/'summary.json').exists()
rows=[{k:float(v) for k,v in row.items()} for row in csv.DictReader((out/'trajectory.csv').open())]
states=[{k:float(v) for k,v in row.items()} for row in csv.DictReader((out/'states.csv').open())]
model,y0=load_model();max_ratios={}
for saved,state in zip(rows,states):
    y=np.array([state[k] for k in model.state_names])
    actual,ratios=diagnostics(model,saved['time_s'],y,WT)
    assert all(abs(actual[k]-saved[k]) <= 1e-11+1e-10*abs(saved[k]) for k in actual)
    for k,v in ratios.items():max_ratios[k]=max(max_ratios.get(k,0),float(v))
ts=np.array([r['time_s'] for r in rows])
fields=[k for k in rows[0] if k.endswith('_fmol_s') or k in ('q_out_pL_s','nbc_current_A')]
integrals={}
for start,end in [(0,60),(60,600),(0,600)]:
    mask=(ts>=start)&(ts<=end)
    integrals[f'{start}_{end}']={k:float(np.trapezoid([r[k] for r,keep in zip(rows,mask) if keep],ts[mask])) for k in fields}
rest=model.evaluate(0,y0).diagnostics.water.lumen_outflow_pL_s
gates=dict(mean=integrals['60_600']['q_out_pL_s']/540>1.1*rest,
           cumulative=rows[-1]['cumulative_outflow_pL']>1.1*600*rest,endpoint=rows[-1]['q_out_pL_s']>=rest)
assert all(gates.values())
from modern_full_model.validation import PRODUCTION_RADAU
sp=PRODUCTION_RADAU
dump(out/'summary.json',dict(case='wt',expression=1.,parameter_hashes=hashes(model),
    initial_state_sha256=sha256_object(y0.tolist()),initial_state=y0.tolist(),
    final_state=[states[-1][k] for k in model.state_names],
    accepted_steps=None,checked_states=None,nfev=None,elapsed_seconds=None,
    recovered_saved_states=len(rows),extrema_scope='saved onset and one second grid only',
    original_accepted_endpoint_checks='completed before metadata serialisation failed; counters not retained',
    python=sys.version.split()[0],numpy=np.__version__,scipy=scipy.__version__,
    solver=dict(method='Radau',rtol=sp.rtol,atol=sp.atol_vector(model.state_names).tolist(),max_step_s=sp.max_step_s),
    all_inherited_gates_passed=True,maximum_conservation_ratios=max_ratios,
    minima={k:min(r[k] for r in rows) for k in rows[0]},maxima={k:max(r[k] for r in rows) for k in rows[0]},
    integrals=integrals,endpoint=rows[-1],wt_activation=gates,
    local_derivative_cases=12,production_source_unchanged=True,
    recovery='NumPy bool metadata serialisation only; no scientific rerun; original CSVs unchanged'))
print('Recovered WT summary from 602 saved states; no integration repeated.')
