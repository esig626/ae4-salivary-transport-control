"""Task52 paired chloride-reservoir experiment on byte-exact Task37 sources.

No model is constructed at import. The isolated package contains only the
historical pre-Palk dependency closure. All source replacements are explicit
electroneutral cycle vectors, after the final local NBC/homeostasis evaluation.
"""
from pathlib import Path
from dataclasses import replace
import importlib.util
import json
import math
import sys
import numpy as np

TASK = Path(__file__).resolve().parent
FROZEN = TASK / 'frozen_task37'
PACKAGE = 'task52_pre_palk'
if PACKAGE not in sys.modules:
    package_path = FROZEN / 'src/modern_full_model'
    spec = importlib.util.spec_from_file_location(
        PACKAGE, package_path/'__init__.py', submodule_search_locations=[str(package_path)])
    module = importlib.util.module_from_spec(spec)
    sys.modules[PACKAGE] = module
    spec.loader.exec_module(module)

from task52_pre_palk.model import ModernFullModel, WT, AE4_NULL
from task52_pre_palk.task31_nhe1_repair import build_model, load_background
from task52_pre_palk.nbc_minimal import MinimalNbcModel
from task52_pre_palk.nkcc_stimulation import StimulatedNkcc1Model
from task52_pre_palk.membranes import hill_activation, nernst_voltage_V, current_to_fmol_s
from task52_pre_palk.validation import StimulusArm, sha256_object

ALLOWED_G = (0.0, 2.32e-9, 4.49e-9)
ALLOWED_SUPPLY = (1.0, 1.05, 1.10)
CARRIER = 2.339370005697548e-05
CHECKPOINT = FROZEN / 'results/31_nhe1_mechanistic_repair/rest_checkpoints/185f5cee3bcec7a023cf82e4092872114f2292732ad9c4bab68ecf4e062f76e3.json'
COUNTS = {'parent_factories': 0, 'constructor_core_evaluations': 0,
          'local_evaluations': 0, 'paired_evaluations': 0}

def accepted_state():
    saved = json.loads(CHECKPOINT.read_text())
    assert saved['admissible'] and saved['background'] == 'R09'
    assert saved['ae4_expression'] == 1.0 and saved['carrier_amount_fmol'] == CARRIER
    assert sha256_object(saved['core_state']) == saved['core_state_sha256']
    y = np.r_[np.asarray(saved['core_state'], float), 0.0]
    assert sha256_object(y.tolist()) == '879a650761e9044691ab6b2a9b1bfc3be8f951ece176925117596cc2baa31abe'
    return y

def clone_parent(parent, *, parameters=None, protocol=None):
    """Rebuild identical historical wrappers; change only declared inputs."""
    core = parent.base_model.base_model
    stimulus = core.stimulus if protocol is None else protocol
    regulation = replace(core.regulatory_model, stimulus=stimulus)
    rebuilt = ModernFullModel(core.parameters if parameters is None else parameters,
        stimulus=stimulus, regulatory_model=regulation,
        ae4_parameters=core.ae4_parameters, ae4_evaluator=core.ae4_evaluator)
    COUNTS['constructor_core_evaluations'] += 1
    return MinimalNbcModel(StimulatedNkcc1Model(rebuilt), parent.nbc_parameters)

def make_parent(protocol='CCH_IPR'):
    if protocol not in ('CCH_IPR', 'CCH_ONLY', 'IPR_ONLY'):
        raise ValueError('Only the three frozen Task52 protocols are permitted')
    COUNTS['parent_factories'] += 1
    COUNTS['constructor_core_evaluations'] += 1
    parent = MinimalNbcModel(build_model(load_background('R09'),
                                        carrier_amount_fmol=CARRIER, stimulated=True))
    saved = json.loads(CHECKPOINT.read_text())
    assert sha256_object(parent.parameters) == saved['active_parameter_sha256']
    assert sha256_object(parent.ae4_parameters) == saved['ae4_parameter_sha256']
    if protocol != 'CCH_IPR':
        parent = clone_parent(parent, protocol=replace(parent.stimulus, arm=StimulusArm(protocol)))
    return parent

class SharedAuxiliaryModel:
    """Identical beta*G auxiliary current evaluated in full current closure.

    Protocol calcium is constant on its active interval. Adding G/gate to
    the stored maximal Cl conductance therefore adds exactly beta*G to the
    effective conductance. The original template is dispatched at beta zero.
    This reuses the complete unchanged NBC and IPR-only electrical closures.
    Diagnostic cl_apical is explicitly total; named components are separate.
    """
    def __init__(self, parent, G_aux_S):
        if G_aux_S not in ALLOWED_G:
            raise ValueError('G_aux must be one of the three predeclared values')
        self.parent = parent
        self.G_aux_S = float(G_aux_S)
        self.parameters = parent.parameters  # original, unchanged scientific parameters
        self.stimulus = parent.stimulus
        self.state_names = parent.state_names
        self.active = parent
        if G_aux_S and self.stimulus.contains_ipr:
            calcium = self.stimulus(1e-6).calcium_uM
            mp = parent.parameters.membranes
            gate = hill_activation(calcium, mp.calcium_half_uM, mp.calcium_hill)
            assert gate > 0
            effective = replace(parent.parameters, membranes=replace(mp,
                g_cl_apical_S=mp.g_cl_apical_S + self.G_aux_S/gate))
            self.active = clone_parent(parent, parameters=effective)

    def evaluate(self, t, y, *, genotype=WT):
        COUNTS['local_evaluations'] += 1
        stimulus = self.stimulus(float(t))
        beta = stimulus.beta_input
        if beta not in (0.0, 1.0):
            raise ValueError('Task52 has fixed zero/one beta inputs')
        model = self.active if beta else self.parent
        ev = model.evaluate(t, y, genotype=genotype)
        d = ev.diagnostics
        mp, c = self.parameters.membranes, self.parameters.constants
        ci, cl = d.observables.cell_concentrations_mM, d.observables.lumen_concentrations_mM
        ecl = nernst_voltage_V(ci['cl'], cl['cl'], valence=-1,
                              thermal_voltage_V=c.thermal_voltage_V)
        g_native = mp.g_cl_apical_S * hill_activation(stimulus.calcium_uM,
                                      mp.calcium_half_uM, mp.calcium_hill)
        g_aux = beta*self.G_aux_S
        drive = d.membranes.v_apical_V - ecl
        i_aux, i_native = g_aux*drive, g_native*drive
        currents = dict(d.membranes.currents_A)
        currents.update(cl_tmem16a=i_native, cl_auxiliary=i_aux,
                        cl_apical_total=currents['cl_apical'])
        regulatory = dict(d.regulatory)
        regulatory.update(task52_G_aux_S=self.G_aux_S, task52_g_aux_effective_S=g_aux,
            task52_g_tmem16a_effective_S=g_native, task52_E_Cl_V=ecl,
            task52_Cl_driving_force_V=-drive,
            task52_apical_current_component_residual_A=currents['cl_apical']-i_native-i_aux,
            task52_current_label='cl_apical and cl_apical_total are TOTAL; cl_tmem16a and cl_auxiliary are components')
        return replace(ev, diagnostics=replace(d,
            membranes=replace(d.membranes, currents_A=currents), regulatory=regulatory))

class PairedModel:
    def __init__(self, *, protocol='CCH_IPR', G_aux_S=2.32e-9, ko_supply_multiplier=1.0):
        if ko_supply_multiplier not in ALLOWED_SUPPLY:
            raise ValueError('Only matched, x1.05 and x1.10 KO supply are permitted')
        self.local = SharedAuxiliaryModel(make_parent(protocol), G_aux_S)
        self.ko_supply_multiplier = float(ko_supply_multiplier)
        self.parameters = self.local.parameters
        self.single_state_names = self.local.state_names
        self.state_names = tuple(g+'__'+n for g in ('WT','KO') for n in self.single_state_names)
        assert len(self.single_state_names) == 13

    def evaluate(self, t, y):
        COUNTS['paired_evaluations'] += 1
        y = np.asarray(y, float)
        if y.shape != (26,):
            raise ValueError('Paired state must have exactly 26 coordinates')
        wt = self.local.evaluate(t, y[:13], genotype=WT)
        ko = self.local.evaluate(t, y[13:], genotype=AE4_NULL)
        hw, hk = wt.diagnostics.homeostasis, ko.diagnostics.homeostasis
        N = self.ko_supply_multiplier * hw.nkcc1_inward_fmol_s
        E = self.ko_supply_multiplier * hw.ae2_inward_fmol_s
        dN, dE = N-hk.nkcc1_inward_fmol_s, E-hk.ae2_inward_fmol_s
        vector = np.asarray((dN, dN, 2*dN+dE, -dE, -dE))
        rhs = np.asarray(ko.rhs).copy()
        rhs[:5] += vector
        fields = ('na_cell_fmol_s','k_cell_fmol_s','cl_cell_fmol_s',
                  'tic_cell_fmol_s','alkalinity_cell_fmol_s')
        updated = {name:float(getattr(hk,name)+inc) for name,inc in zip(fields,vector)}
        qh = updated['na_cell_fmol_s']+updated['k_cell_fmol_s']-updated['cl_cell_fmol_s']-updated['alkalinity_cell_fmol_s']
        homeostasis = replace(hk, nkcc1_inward_fmol_s=N, ae2_inward_fmol_s=E,
                              charge_source_residual_fmol_s=qh, **updated)
        d = ko.diagnostics
        r = dict(d.regulatory)
        r.update(task52_ko_supply_multiplier=self.ko_supply_multiplier,
            task52_local_unimposed_nkcc1_cycles=hk.nkcc1_inward_fmol_s,
            task52_local_unimposed_ae2_cycles=hk.ae2_inward_fmol_s,
            task52_imposed_nkcc1_cycles=N, task52_imposed_ae2_cycles=E)
        residuals = dict(d.conservation_residuals)
        residuals['cell_bulk_charge_rate_fmol_s'] = float(rhs[0]+rhs[1]-rhs[2]-rhs[4])
        residuals['homeostasis_charge_fmol_s'] = float(qh)
        m, out = d.membranes, d.outflow_sources_fmol_s
        external_carbon = (homeostasis.tic_cell_fmol_s + d.ae4.tic_cell_fmol_s
            + m.cell_sources_fmol_s['tic'] + d.co2_fluxes_fmol_s['bath_to_cell']
            + m.lumen_sources_fmol_s['tic'] - out['tic'])
        residuals['carbon_accounting_fmol_s'] = float(rhs[3]+rhs[9]-external_carbon)
        ko = replace(ko, rhs=rhs, diagnostics=replace(d, homeostasis=homeostasis,
                     regulatory=r, conservation_residuals=residuals))
        return wt, ko

    def rhs(self, t, y):
        wt, ko = self.evaluate(t, y)
        return np.r_[wt.rhs, ko.rhs]

def dependency_audit():
    """Runtime module-origin audit, without model construction/evaluation."""
    forbidden = ('nkcc1_palk2010','ae4_cacc_recruitment','task41_selected',
                 'task50_effective_coupling','ae4_equal_cation_routing')
    loaded = {name:str(Path(m.__file__).resolve()) for name,m in sys.modules.items()
              if name.startswith(PACKAGE) and getattr(m,'__file__',None)}
    bad = [name for name in sys.modules if any(f in name for f in forbidden)]
    escaped = {name:path for name,path in loaded.items()
               if not Path(path).is_relative_to(FROZEN.resolve())}
    return {'pass':not bad and not escaped, 'forbidden_loaded':bad,
            'escaped_snapshot':escaped, 'loaded_scientific_modules':loaded}
