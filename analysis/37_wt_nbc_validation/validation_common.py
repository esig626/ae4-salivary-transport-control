"""Task 37 fixed inputs and diagnostics. No solve, search, or parameter changes."""
from pathlib import Path
from dataclasses import asdict
import json
import math
import os
import sys

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS", "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
from modern_full_model.model import WT
from modern_full_model.task31_nhe1_repair import build_model, load_background
from modern_full_model.nbc_minimal import MinimalNbcModel, MinimalNbcParameters
from modern_full_model.validation import CONSERVATION_RESIDUAL_TOLERANCES, sha256_object
from modern_full_model.membranes import current_to_fmol_s

ANALYSIS = ROOT / "analysis/37_wt_nbc_validation"
RESULTS = ROOT / "results/37_wt_nbc_validation"
PREPARED_HEAD = "d925785512a33b241ec24843660ff3ef88458549"
BRANCH = "codex/task-37-wt-nbc-validation"
CARRIER = 2.339370005697548e-05
Q_REST = 0.0010757073853493032
CHECKPOINT = ROOT / "results/31_nhe1_mechanistic_repair/rest_checkpoints/185f5cee3bcec7a023cf82e4092872114f2292732ad9c4bab68ecf4e062f76e3.json"
AMOUNT_ROWS = [0, 1, 2, 3, 4, 6, 7, 8, 9, 10]
VOLUME_ROWS = [5, 11]
EXTRA_SOURCES = {"minimal_nbc_charge_source_fmol_s", "minimal_nbc_carbon_source_fmol_s"}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def models_and_state():
    saved = json.loads(CHECKPOINT.read_text())
    assert saved["admissible"] and saved["background"] == "R09"
    assert saved["genotype"] == "WT" and saved["ae4_expression"] == 1.0
    assert saved["carrier_amount_fmol"] == CARRIER
    assert sha256_object(saved["core_state"]) == saved["core_state_sha256"]
    background = load_background("R09")
    rest = build_model(background, carrier_amount_fmol=CARRIER)
    stim_base = build_model(background, carrier_amount_fmol=CARRIER, stimulated=True)
    assert sha256_object(rest.parameters) == saved["active_parameter_sha256"]
    assert sha256_object(rest.ae4_parameters) == saved["ae4_parameter_sha256"]
    assert rest.parameters == stim_base.parameters
    assert rest.ae4_parameters == stim_base.ae4_parameters
    assert WT.ae4_expression == 1.0
    y0 = np.r_[np.asarray(saved["core_state"], dtype=float), rest.initial_state()[12:]]
    assert np.array_equal(y0[:12], saved["core_state"])
    assert np.array_equal(y0[12:], stim_base.initial_state()[12:])
    assert tuple(y0[12:]) == (0.0,)
    nbc = MinimalNbcParameters()
    assert nbc.capacity_fmol_s == 0.11570913197464398 and nbc.log_width == 2.0
    return saved, rest, MinimalNbcModel(rest), MinimalNbcModel(stim_base), y0


def row_for(model, t, y, evaluation=None):
    e = evaluation if evaluation is not None else model.evaluate(float(t), y, genotype=WT)
    d = e.diagnostics
    c = d.observables.cell_concentrations_mM
    m, h, r = d.membranes, d.homeostasis, d.regulatory
    f = model.parameters.constants.faraday_C_mol
    return {
        "time_s": float(t), "q_out_pL_s": float(d.water.lumen_outflow_pL_s),
        "na_i_mM": float(c["na"]), "k_i_mM": float(c["k"]), "cl_i_mM": float(c["cl"]),
        "tic_i_mM": float(c["tic"]), "hco3_i_mM": float(c["hco3"]),
        "ph_i": float(d.observables.cell_acid_base.ph), "volume_i_pL": float(y[5]),
        "osm_lumen_mOsm": float(d.observables.osmolarities_mOsm["lumen"]),
        "v_apical_V": float(m.v_apical_V), "v_basolateral_V": float(m.v_basolateral_V),
        "v_transepithelial_V": float(m.v_transepithelial_V),
        "nbc_activation": float(r["minimal_nbc_activation_fraction"]),
        "nbc_cycle_inward_fmol_s": float(r["minimal_nbc_cycle_inward_fmol_s"]),
        "nbc_current_A": float(r["minimal_nbc_current_A"]),
        "nhe1_multiplier": float(r["nhe1_stimulation_multiplier"]),
        "nhe1_inward_fmol_s": float(h.nhe1_inward_fmol_s),
        "nkcc1_cl_inward_fmol_s": float(2.0 * h.nkcc1_inward_fmol_s),
        "ae4_cl_inward_fmol_s": float(d.ae4.cl_cell_fmol_s),
        "ae2_cl_inward_signed_fmol_s": float(h.ae2_inward_fmol_s),
        "cacc_cl_export_fmol_s": float(current_to_fmol_s(m.currents_A["cl_apical"], valence=-1, faraday_C_mol=f)),
        "paracellular_cl_return_fmol_s": float(current_to_fmol_s(m.currents_A["para_cl"], valence=-1, faraday_C_mol=f)),
    }


def bookkeeping(model, e):
    """Independent conserved-row assembly from the declared transporter cycles."""
    d = e.diagnostics
    h, m, a, r = d.homeostasis, d.membranes, d.ae4, d.regulatory
    jn, jh, j2 = h.nkcc1_inward_fmol_s, h.nhe1_inward_fmol_s, h.ae2_inward_fmol_s
    jb = r["minimal_nbc_cycle_inward_fmol_s"]
    f = model.parameters.constants.faraday_C_mol
    jk_a = current_to_fmol_s(m.currents_A["k_apical"], valence=1, faraday_C_mol=f)
    jk_b = current_to_fmol_s(m.currents_A["k_basolateral"], valence=1, faraday_C_mol=f)
    jc = current_to_fmol_s(m.currents_A["cl_apical"], valence=-1, faraday_C_mol=f)
    pump = m.pump_apical_fmol_s + m.pump_basolateral_fmol_s
    co2 = sum(d.co2_fluxes_fmol_s.values())
    expected = np.asarray([
        jn + jh + jb + a.na_cell_fmol_s - 3 * pump,
        jn + a.k_cell_fmol_s + 2 * pump - jk_a - jk_b,
        2 * jn + j2 + a.cl_cell_fmol_s - jc,
        co2 + 2 * jb - j2 + a.tic_cell_fmol_s,
        jh + 2 * jb - j2 + a.alkalinity_cell_fmol_s,
        d.water.bath_to_cell_pL_s - d.water.cell_to_lumen_pL_s,
    ])
    isolated = np.asarray([
        m.cell_sources_fmol_s["na"] + 3 * pump,
        m.cell_sources_fmol_s["k"] - 2 * pump + jk_a + jk_b,
        m.cell_sources_fmol_s["cl"] + jc,
        m.cell_sources_fmol_s["tic"], m.cell_sources_fmol_s["alkalinity"],
    ])
    return {
        "cell_rhs_assembly_max_abs_fmol_s": float(np.max(np.abs(np.asarray(e.rhs[:5]) - expected[:5]))),
        "cell_water_assembly_abs_pL_s": float(abs(e.rhs[5] - expected[5])),
        "isolated_nbc_source": dict(zip(("na", "k", "cl", "tic", "ta"), map(float, isolated))),
        "isolated_nbc_signature_max_abs_fmol_s": float(np.max(np.abs(isolated - [jb, 0, 0, 2 * jb, 2 * jb]))),
        "nbc_source_plus_current_equivalents_fmol_s": float(jb - 2 * jb + r["minimal_nbc_current_A"] / f * 1e15),
    }


def all_numeric_finite(value):
    if isinstance(value, dict):
        return all(all_numeric_finite(v) for v in value.values())
    if isinstance(value, (tuple, list)):
        return all(all_numeric_finite(v) for v in value)
    if isinstance(value, (float, int, np.number)):
        return bool(np.isfinite(value))
    return True


def diagnose(model, t, y, e):
    row = row_for(model, t, y, e)
    failures = []
    if not (np.all(np.isfinite(y)) and np.all(np.isfinite(e.rhs))
            and all_numeric_finite(asdict(e.diagnostics))):
        failures.append({"gate": "finite_states_and_observables"})
    for k, v in enumerate(y[:12]):
        if v <= 0:
            failures.append({"gate": "positive_" + model.state_names[k], "value": float(v), "required": "> 0"})
    limits = (
        ("na_i_mM", None, 40.0, False), ("k_i_mM", 50.0, 200.0, True),
        ("cl_i_mM", 30.0, 80.0, True), ("ph_i", 6.6, 7.3, True),
        ("hco3_i_mM", 0.0, 100.0, False), ("tic_i_mM", 0.0, None, False),
        ("volume_i_pL", 0.0, 3.0, False),
    )
    for key, low, high, inclusive in limits:
        v = row[key]
        lower_bad = low is not None and (v < low if inclusive else v <= low)
        upper_bad = high is not None and (v > high if inclusive else v >= high)
        if lower_bad or upper_bad:
            failures.append({"gate": key, "value": v, "lower": low, "upper": high, "inclusive": inclusive})
    residuals = e.diagnostics.conservation_residuals
    if set(residuals) != set(CONSERVATION_RESIDUAL_TOLERANCES) | EXTRA_SOURCES:
        failures.append({"gate": "declared_conservation_diagnostic_keys"})
    ratios = {}
    for key, tol in CONSERVATION_RESIDUAL_TOLERANCES.items():
        ratios[key] = abs(float(residuals[key])) / tol
        if not math.isfinite(ratios[key]) or ratios[key] > 1:
            failures.append({"gate": key, "value": float(residuals[key]), "absolute_tolerance": tol})
    # These two new fields are signed physical NBC sources, not zero residuals.
    jb = row["nbc_cycle_inward_fmol_s"]
    for key, expected in (("minimal_nbc_charge_source_fmol_s", -jb),
                          ("minimal_nbc_carbon_source_fmol_s", 2 * jb)):
        if residuals[key] != expected:
            failures.append({"gate": key, "value": residuals[key], "expected": expected})
    return row, failures, ratios
