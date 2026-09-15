"""Task 45 completion: the prescribed Task 41 inverse family only.

All 600 s experiments retain the frozen WT conserved initial state and the
original Task 41 stimulus, source law, recruitment law and physiological gates.
Longer constant stimulation is a separately labelled diagnostic.  No production
files or historical results are written.  A cache stores compact new results;
--recompute ignores it for the clean checkout reproduction gate.
"""
from core import *
from scipy.optimize import brentq
from time import perf_counter
import argparse

HERE = Path(__file__).resolve().parent
OUT = HERE / 'output'
DETAIL = HERE / 'inverse_detail'
CACHE = DETAIL / 'trials'
GRID = np.r_[0., 1e-6, np.arange(1., 601.)]
TARGET = .303
SELECTED_RHO = 1. - .10511872843288446
SETTINGS = {'standard': {'rtol': 1e-8, 'max_step': 2.},
            'tight': {'rtol': 2e-9, 'max_step': 1.}}
PARAMETERS = ('NBC_capacity', 'AE4_carrier', 'NKCC_scale', 'NHE_carrier',
              'pump_capacity', 'CaCC_conductance', 'buffer_pool', 'water_apical')


def write_json(path, obj):
    """Atomic summaries allow independent parameter subprocesses to share cache."""
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.'+str(os.getpid())+'.tmp')
    temporary.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+'\n')
    temporary.replace(path)


def build_parameter(stim, name='baseline', h=0., constant=False):
    """One logarithmic perturbation; no resting solve and no retuning."""
    p = stim.parameters
    factor = float(np.exp(h))
    kw = {'constant': constant}
    if name == 'NBC_capacity':
        kw['nbc_parameters'] = replace(stim.nbc_parameters, capacity_fmol_s=stim.nbc_parameters.capacity_fmol_s*factor)
    elif name == 'AE4_carrier':
        kw['ae4_parameters'] = replace(stim.ae4_parameters, carrier_amount_fmol=stim.ae4_parameters.carrier_amount_fmol*factor)
    elif name == 'NKCC_scale':
        kw['nkcc1_kinetics'] = replace(stim.nkcc1_kinetics, alpha_eff_fmol_s=stim.nkcc1_kinetics.alpha_eff_fmol_s*factor)
    elif name == 'NHE_carrier':
        kw['parameters'] = replace(p, homeostasis=replace(p.homeostasis, nhe1_cha_carrier_amount_fmol=p.homeostasis.nhe1_cha_carrier_amount_fmol*factor))
    elif name == 'pump_capacity':
        kw['parameters'] = replace(p, membranes=replace(p.membranes, nak_capacity_fmol_s=p.membranes.nak_capacity_fmol_s*factor))
    elif name == 'CaCC_conductance':
        kw['parameters'] = replace(p, membranes=replace(p.membranes, g_cl_apical_S=p.membranes.g_cl_apical_S*factor))
    elif name == 'buffer_pool':
        kw['parameters'] = replace(p, geometry=replace(p.geometry, cell_buffer_total_fmol=p.geometry.cell_buffer_total_fmol*factor))
    elif name == 'water_apical':
        kw['parameters'] = replace(p, water=replace(p.water, apical_hydraulic_pL_s_mOsm=p.water.apical_hydraulic_pL_s_mOsm*factor))
    elif name != 'baseline':
        raise ValueError(name)
    return clone(stim, **kw)


class Experiment:
    def __init__(self, recompute=False):
        self.rest, self.stim, self.y0, self.task, audit = reference()
        self.recompute = recompute
        self.memory = {}
        self.qrest = observe(self.rest, 0., self.y0)['q_pL_s']
        for folder in (OUT, DETAIL, CACHE):
            folder.mkdir(parents=True, exist_ok=True)
        self.identity = dict(initial_state_sha256=sha256_object(self.y0.tolist()),
            production_parameter_hashes=self.task.parameter_hashes(self.stim),
            source_python_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted((ROOT/'src/modern_full_model').glob('*.py'))},
            gate_source='analysis/40_ae4_equal_cation_routing/validation_common.py::diagnose',
            original_grid=GRID.tolist(), target_deficit=TARGET,
            target_interpretation='Task comparison threshold 30.3 percent, one reported standard error below the experimental mean 35 percent (SE 4.7 percentage points); not a confidence interval or biological admissible range',
            solver_atol=1e-11, settings=SETTINGS, state_count=13,
            scope='Original Task 41 family; fixed inherited WT initial conserved state; local parameter changes only')
        write_json(DETAIL/'protocol.json', self.identity)

    def trial(self, rho, expression=0., precision='standard', parameter='baseline', h=0., end=600., constant=False, trajectory=False):
        args = dict(rho=float(rho), expression=float(expression), precision=precision,
            parameter=parameter, log_parameter_shift=float(h), end_s=float(end), constant=constant)
        numerical_identity={k:v for k,v in self.identity.items() if k!='target_interpretation'}
        key = hashlib.sha256(json.dumps([numerical_identity, args], sort_keys=True).encode()).hexdigest()
        path = CACHE/(key+'.json')
        if key in self.memory and not trajectory:
            return self.memory[key]
        if path.exists() and not self.recompute and not trajectory:
            rec = json.loads(path.read_text()); self.memory[key] = rec; return rec
        m = cacc(build_parameter(self.stim, parameter, h, constant), float(rho))
        genotype = replace(WT, ae4_expression=float(expression))
        cfg = dict(SETTINGS[precision])
        sol = integrate(m, self.y0, genotype, end=end, **cfg)
        sample_grid = GRID if end == 600 else np.r_[GRID, np.arange(601., end+1.)]
        monitor = np.unique(np.r_[sample_grid, sol.t])
        failures = {}; failure_details = {}; maxima = {}; minima = {}; conservation = {}
        grid_rows = {}; all_rows = []
        for t in monitor:
            y = sol.sol(t)[:-1]
            ev = m.evaluate(float(t), y, genotype=genotype)
            row, failed, ratios = self.task.diagnose(m, float(t), y, ev)
            row.update(cacc_recruitment_factor=float(ev.diagnostics.regulatory['ae4_dependent_cacc_recruitment_factor']),
                effective_cacc_conductance_S=float(ev.diagnostics.regulatory['effective_cacc_conductance_S']))
            for fail in failed:
                failures.setdefault(fail['gate'], float(t))
                failure_details.setdefault(fail['gate'], {'time_s': float(t), **fail})
            for field, value in row.items():
                minima[field] = min(minima.get(field, float('inf')), value)
                maxima[field] = max(maxima.get(field, -float('inf')), value)
            for field, ratio in ratios.items():
                conservation[field] = max(conservation.get(field, 0.), ratio)
            if t in sample_grid:
                grid_rows[float(t)] = row
            if trajectory:
                all_rows.append(row)
        rows = [grid_rows[float(t)] for t in sample_grid]
        q = np.array([r['q_out_pL_s'] for r in rows])
        time = np.array([r['time_s'] for r in rows])
        window = (time >= 60.) & (time <= 600.)
        first = time <= 600.
        sampled = float(np.trapezoid(q[first], time[first]))
        adaptive = float(sol.sol(600.)[-1])
        mean = float(np.trapezoid(q[window], time[window])/540.)
        activation = dict(mean_60_600_gt_1p10_rest=mean > 1.10*self.qrest,
            cumulative_0_600_gt_1p10_rest=sampled > 1.10*600*self.qrest,
            endpoint_at_least_rest=float(q[first][-1]) >= self.qrest)
        if expression == 1. and not all(activation.values()):
            failures['WT_stimulation_requirement'] = 600.
            failure_details['WT_stimulation_requirement'] = activation
        rec = dict(**args, b=1.-float(rho), tight=precision == 'tight',
            cumulative_adaptive_pL=adaptive, cumulative_sampled_pL=sampled,
            endpoint_cumulative_adaptive_pL=float(sol.y[-1,-1]),
            physiological=not failures, physiological_600s=not any(t <= 600. for t in failures.values()),
            failed_gates=failures, failure_details=failure_details,
            minimum=minima, maximum=maxima, ph_min=minima['ph_i'], ph_max=maxima['ph_i'],
            cl_min=minima['cl_i_mM'], max_conservation_ratio=max(conservation.values()),
            conservation_ratio_by_gate=conservation, nfev=int(sol.nfev), njev=int(sol.njev),
            monitored_points=len(monitor), solver_steps=len(sol.t),
            activation_checks=activation, activation_required=expression == 1.,
            initial_state_sha256=self.identity['initial_state_sha256'],
            final_state=sol.y[:-1,-1].tolist())
        if expression != 1.:
            wt = self.trial(0., 1., precision, parameter, h, 600., constant)
            rec['adaptive_deficit'] = 1.-adaptive/wt['cumulative_adaptive_pL']
            rec['sampled_deficit'] = 1.-sampled/wt['cumulative_sampled_pL']
            rec['WT_reference_physiological'] = wt['physiological']
            rec['WT_cumulative_adaptive_pL'] = wt['cumulative_adaptive_pL']
            rec['WT_cumulative_sampled_pL'] = wt['cumulative_sampled_pL']
        if trajectory:
            rec['precise_ph_gate_crossing_s'] = None
            for a, b in zip(all_rows, all_rows[1:]):
                if a['ph_i'] <= 7.3 < b['ph_i']:
                    rec['precise_ph_gate_crossing_s'] = float(brentq(lambda t: observe(m,t,sol.sol(t)[:-1],genotype)['ph_i']-7.3, a['time_s'], b['time_s'], xtol=1e-7))
                    break
            label = 'selected' if abs(rho-SELECTED_RHO)<1e-12 else 'crossing'
            write_csv(DETAIL/(label+'_long_time.csv'), rows)
            if label == 'crossing':
                compatible = [observe(m, t, sol.sol(t)[:-1], genotype) for t in np.arange(0.,end+1.,5.)]
                write_csv(OUT/'inverse_long_time.csv', compatible)
        write_json(path, rec); self.memory[key] = rec
        return rec

    def threshold(self, observable, precision='standard', parameter='baseline', h=0., bracket=(.89,.90)):
        fun = lambda rho: self.trial(rho, precision=precision, parameter=parameter, h=h)[observable+'_deficit']-TARGET
        a,b = bracket
        if fun(a)*fun(b) >= 0:
            a,b = .85,.95
        value = float(brentq(fun, a, b, xtol=2e-9))
        rec = self.trial(value, precision=precision, parameter=parameter, h=h)
        assert rec['physiological'] and rec['WT_reference_physiological'], rec['failed_gates']
        assert abs(rec[observable+'_deficit']-TARGET) < 1e-7
        return value, rec

    def baseline(self):
        for precision in SETTINGS:
            self.trial(0.,1.,precision)
        selected = {}
        for e in [1., .05, 0.]:
            rec = self.trial(SELECTED_RHO,e)
            case = {1.:'wt', .05:'ae4_5pct', 0.:'ae4_null'}[e]
            frozen = list(csv.DictReader((ROOT/'results/41_ae4_loss_algebraic_design/candidate_08'/(case+'_dense_timeseries.csv')).open()))
            recorded = float(frozen[-1]['cumulative_outflow_0_t_pL'])
            rec = dict(rec, frozen_task41_sampled_pL=recorded,
                reproduction_sampled_error_pL=rec['cumulative_sampled_pL']-recorded)
            assert abs(rec['reproduction_sampled_error_pL']) < 2e-8
            selected[case] = rec
        write_json(DETAIL/'selected_reproduction.json',selected)
        scan = [self.trial(r) for r in [0.,.1,.2,.3,.4,.5,.6,.7,.75,.8,.85,.875,.89,SELECTED_RHO,.90,.925,.95]]
        write_json(DETAIL/'ordered_family_scan.json',scan)
        crossings = {}
        for obs in ['adaptive','sampled']:
            standard, check = self.threshold(obs)
            tight, tight_check = self.threshold(obs,'tight')
            near = [self.trial(standard+delta) for delta in [-1e-4,1e-4]]
            crossings[obs] = dict(rho=standard,b=1-standard,ae4_dependent_recruitment_percent=100*standard,
                check=check,tight_rho=tight,tight_check=tight_check,
                solver_tolerance_rho_difference=tight-standard,nearby=near)
            assert near[0][obs+'_deficit'] < TARGET < near[1][obs+'_deficit']
            assert abs(tight-standard)<2e-6
        write_json(DETAIL/'crossings.json',crossings)
        all_base = [r for r in self.memory.values() if r['parameter']=='baseline' and r['expression']==0. and r['end_s']==600. and not r['constant']]
        write_json(OUT/'inverse_samples.json',all_base)
        a = crossings['adaptive']; s = crossings['sampled']
        write_json(OUT/'inverse_threshold.json',dict(target_lower_deficit=TARGET,target_upper_deficit=.397,
            target_interpretation=self.identity['target_interpretation'],
            numerical_first_crossing_rho=a['rho'],tight_check=self.trial(a['rho'],precision='tight'),
            tight_refined_rho=a['tight_rho'],nearby=a['nearby'],selected_task41_reproduction=selected['ae4_null'],
            monotone_at_sampled_points=bool(np.all(np.diff([r['adaptive_deficit'] for r in scan])>0)),
            scope='First numerical crossing in the ordered sampled prescribed family; monotonicity between samples and global exclusion are unproved.'))
        write_json(OUT/'inverse_legacy_threshold.json',dict(numerical_rho=s['rho'],check=self.trial(s['rho'],precision='tight'),
            tight_refined_rho=s['tight_rho'],quadrature='Original one-second grid including 0 and 1e-6 s'))
        print('INVERSE crossings', {k:round(v['rho'],9) for k,v in crossings.items()}, flush=True)
        return crossings

    def long_time(self):
        crossings = json.loads((DETAIL/'crossings.json').read_text())
        output = {}
        for label,rho in [('selected',SELECTED_RHO),('crossing',crossings['adaptive']['rho'])]:
            rec = self.trial(rho,end=3600.,constant=True,trajectory=True)
            rec['protocol_note'] = 'Separate constant full stimulation extension; gates monitored at accepted endpoints and integer seconds. Failure does not alter original 600 s Task 41 result.'
            model=cacc(build_parameter(self.stim,constant=True),rho)
            genotype=replace(WT,ae4_expression=0.)
            try:
                eq=stationary(model,Scaling(model,self.y0),np.array(rec['final_state']),genotype)
                state=np.array(eq['state']);ev=model.evaluate(1.,state,genotype=genotype)
                _,ff,rr=self.task.diagnose(model,1.,state,ev)
                eq['all_inherited_gate_failures']=ff
                eq['all_inherited_gates_pass']=not ff
                eq['maximum_inherited_conservation_ratio']=max(rr.values())
                rec['diagnostic_equilibrium']=eq
            except Exception as ex:
                rec['diagnostic_equilibrium']={'status':'root_not_established','error':str(ex)}
            output[label] = rec
            print('INVERSE long',label,rec['failed_gates'], flush=True)
        write_json(DETAIL/'long_time_summary.json',output)
        rec=output['crossing']
        write_json(OUT/'inverse_long_time_summary.json',dict(ph_gate_crossing_s=rec['precise_ph_gate_crossing_s'],
            end_ph=observe(cacc(build_parameter(self.stim,constant=True),rec['rho']),3600.,rec['final_state'],replace(WT,ae4_expression=0.))['ph_i'],
            end_s=3600.,rho=rec['rho'],all_failed_gates=rec['failed_gates'],all_gate_detail_file='inverse_detail/long_time_summary.json',
            protocol=rec['protocol_note']))

    def sensitivity(self, parameter=None):
        cross=json.loads((DETAIL/'crossings.json').read_text())
        rows=[]
        for name in ([parameter] if parameter else PARAMETERS):
            for obs in ['adaptive','sampled']:
                rho=cross[obs]['rho']; delta=2e-4
                d_rho=(self.trial(rho+delta)[obs+'_deficit']-self.trial(rho-delta)[obs+'_deficit'])/(2*delta)
                estimates=[]
                for h in [1e-3,5e-4]:
                    plus=self.trial(rho,parameter=name,h=h)
                    minus=self.trial(rho,parameter=name,h=-h)
                    d_log_p=(plus[obs+'_deficit']-minus[obs+'_deficit'])/(2*h)
                    implicit=-d_log_p/d_rho
                    # Nearby roots are independently solved after the derivative estimate.
                    rp,_=self.threshold(obs,parameter=name,h=h,bracket=(rho-.004,rho+.004))
                    rm,_=self.threshold(obs,parameter=name,h=-h,bracket=(rho-.004,rho+.004))
                    finite=(rp-rm)/(2*h)
                    error=abs(finite-implicit)/max(abs(finite),1e-8)
                    estimates.append(dict(log_step=h,d_rho_d_log_parameter=implicit,
                        finite_root_d_rho_d_log_parameter=finite,rho_plus=rp,rho_minus=rm,
                        derivative_relative_validation_error=error,
                        perturbed_gates_pass=bool(plus['physiological'] and minus['physiological'] and plus['WT_reference_physiological'] and minus['WT_reference_physiological'])))
                    assert error<2e-3,(name,obs,error)
                result=dict(parameter=name,observable=obs,rho=rho,
                    d_deficit_d_rho=d_rho,d_rho_d_log_parameter=estimates[-1]['d_rho_d_log_parameter'],
                    crossing_log_elasticity=estimates[-1]['d_rho_d_log_parameter']/rho,
                    residual_recruitment_log_elasticity=-estimates[-1]['d_rho_d_log_parameter']/(1.-rho),
                    finite_root_d_rho_d_log_parameter=estimates[-1]['finite_root_d_rho_d_log_parameter'],
                    difference_step_relative_error=abs(estimates[0]['d_rho_d_log_parameter']-estimates[1]['d_rho_d_log_parameter'])/max(abs(estimates[1]['d_rho_d_log_parameter']),1e-8),
                    maximum_validation_relative_error=max(x['derivative_relative_validation_error'] for x in estimates),
                    uncertainty_interval='Not established; local sensitivity only',
                    initialisation='Fixed inherited WT conserved state for WT and null; perturbed WT denominator; no rest solve',estimates=estimates)
                rows.append(result)
            write_json(DETAIL/('sensitivity_'+name+'.json'),rows[-2:])
            print('INVERSE sensitivity',name,round(rows[-2]['crossing_log_elasticity'],6),flush=True)
        self.combine_sensitivity()

    def combine_sensitivity(self):
        rows=[]
        for name in PARAMETERS:
            p=DETAIL/('sensitivity_'+name+'.json')
            if p.exists():rows.extend(json.loads(p.read_text()))
        for row in rows:
            row['residual_recruitment_log_elasticity']=-row['d_rho_d_log_parameter']/(1.-row['rho'])
        write_json(DETAIL/'threshold_sensitivities.json',rows)
        write_csv(DETAIL/'threshold_sensitivities.csv',[{k:v for k,v in r.items() if k!='estimates'} for r in rows])


def main():
    p=argparse.ArgumentParser();p.add_argument('--phase',choices=['baseline','long','sensitivity','combine','all'],default='all')
    p.add_argument('--parameter',choices=PARAMETERS);p.add_argument('--recompute',action='store_true')
    args=p.parse_args();start=perf_counter();e=Experiment(args.recompute)
    if args.phase in ('baseline','all'):e.baseline()
    if args.phase in ('long','all'):e.long_time()
    if args.phase in ('sensitivity','all'):e.sensitivity(args.parameter)
    if args.phase=='combine':e.combine_sensitivity()
    print('INVERSE phase complete',args.phase,round(perf_counter()-start,2),'s',flush=True)


if __name__=='__main__':main()
