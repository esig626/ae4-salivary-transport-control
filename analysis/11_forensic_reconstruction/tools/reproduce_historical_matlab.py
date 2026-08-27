#!/usr/bin/env python3
"""Forensic evaluator for the archived salivary MATLAB equations.

This is an execution harness, not a clean scientific implementation.  It
transcribes the arithmetic expressions in the immutable archived MATLAB files
so that they can be evaluated when MATLAB/Octave is unavailable.  Source-line
provenance and the limitations of translated execution are documented in
``analysis/11_forensic_reconstruction/reproduction.md``.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
from pathlib import Path
from typing import Any, Callable

import numpy as np
from scipy.integrate import solve_ivp
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[3]
ARCHIVE_PROJECT = ROOT / "archive/legacy-2017/Ae4_Dynamics_Project"
RESULTS = ROOT / "results/11_forensic_reconstruction"
FORENSIC_COPY = ROOT / "tmp/forensics/historical_execution_copy"

STATE_NAMES_7 = ("Nal", "Kl", "Hi", "Na", "K", "Cl", "HCO3")
STATE_NAMES_5 = ("Hi", "Na", "K", "Cl", "HCO3")
SCENARIOS = ("wt", "ae2_knockout", "ae4_knockout")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def scalar(value: Any) -> float:
    return float(np.asarray(value).reshape(-1)[0])


def load_par_mat(path: Path = ARCHIVE_PROJECT / "Par.mat") -> tuple[dict[str, float], dict[str, str]]:
    raw = loadmat(path, squeeze_me=False, struct_as_record=False)
    struct = raw["par"][0, 0]
    values: dict[str, float] = {}
    dtypes: dict[str, str] = {}
    for name in struct._fieldnames:
        field = getattr(struct, name)
        values[name] = scalar(field)
        dtypes[name] = str(np.asarray(field).dtype)
    return values, dtypes


def original_calibration() -> dict[str, Any]:
    """Evaluate Original/Parameters.m and calibration lines 13--123 exactly."""
    # Original/Parameters.m, lines 3--89.
    Aa = 1300.0 * 2.0 / (28.7 * 2.51)
    Ab = 1.51 * Aa
    Hi0 = 28.7
    delta = (Aa + Ab) / 2.0
    KCaCC = 0.26 + 100.0
    KCaKC = 0.3 + 100.0
    eta1 = 1.49
    eta2 = 1.7
    r = 1.305e6
    alpha1 = 0.641
    R = 8.3144621
    T = 310.0
    F = 96485.3365
    RTF = 1e3 * R * T / F
    pHl, pHi, pHe = 6.81, 6.91, 7.41
    kn = 2.6e4 * 0.012
    kp = 11.0 * 0.012
    Cle, Nae, Ke = 102.6, 140.2, 5.3
    He = 1e3 * 10 ** (-pHe)
    HCO3e, CO2e = 40.0, 1.6
    Nal0, Kl0 = 118.7, 5.6
    Cll0 = Nal0 + Kl0
    Hl = 1e3 * 10 ** (-pHl)
    HCO3l = Kl0 + Nal0 - Cll0 + Hl
    CO2l = 11.6
    Na0, K0, Cl0 = 25.0, 120.0, 50.0
    H0 = 1e3 * 10 ** (-pHi)
    HCO30 = 10.0
    Va0, Vb0 = -50.24, -62.8
    Vt0 = Va0 - Vb0
    KCl, KB, KNa, KH = 5.6, 1e-4, 15.0, 4.5e-4
    k1 = 1.92e-2
    k2_raw = 1.341840733667157e-5
    alpha2 = 2.7e-2
    b1_pre = alpha2 * 0.01194 / (2 * (K0 + Na0 - Nal0 - Kl0 + H0))
    b2 = 7 * b1_pre
    b3 = 0.94 * b1_pre
    b1 = b1_pre * 74.4

    # Original/Saliva_Ae4.m, lines 13--123.
    Ca = 50e-3
    PCO2 = 0.197e4
    CO20 = (PCO2 * (CO2l + CO2e) - kn * HCO30 * H0) / (2 * PCO2 - kp)
    xl = (b2 / b1) * (
        2 * (Na0 + K0 + H0) + CO20 - (Nae + Ke + Cle + HCO3e)
    ) - (2 * (Nal0 + Kl0 - Na0 - K0 - H0) - CO20)
    s = 1e-3
    alpha_Nkcc1 = 0.02279 * 0.28

    VCaCC = RTF * math.log(Cll0 / Cl0)
    PCaCC = 1 / (1 + (KCaCC / Ca) ** eta1)
    vCaCC = PCaCC * (Va0 + VCaCC) / F
    VCaKC = RTF * math.log(Ke / K0)
    PCaKC = 1 / (1 + (KCaKC / Ca) ** eta2)
    vCaKC = PCaKC * (Vb0 - VCaKC) / F
    VtNa = RTF * math.log(Nal0 / Nae)
    vtNa = (Vt0 - VtNa) / F
    VtK = RTF * math.log(Kl0 / Ke)
    vtK = (Vt0 - VtK) / F

    qa = b1 * (2 * (Nal0 + Kl0 - Na0 - K0 - H0) - CO20 + xl)
    qb = b2 * (2 * (Na0 + K0 + H0) + CO20 - (Nae + Ke + Cle + HCO3e))
    qt = b3 * (2 * (Nal0 + Kl0) + xl - (Nae + Ke + Cle + HCO3e))
    qtot0 = qa + qt

    vNaK = r * (Ke * s) ** 2 * (Na0 * s) ** 3 / (
        (Ke * s) ** 2 + alpha1 * (Na0 * s) ** 3
    )
    a1, a2_raw, a3, a4_raw = 157.55, 2.0096e7, 1.0306, 1.3852e6
    Nkcc1 = alpha_Nkcc1 * (
        (a1 - a2_raw * (Na0 * s) * (K0 * s) * (Cl0 * s) ** 2)
        / (a3 + a4_raw * (Na0 * s) * (K0 * s) * (Cl0 * s) ** 2)
    )
    vAe4 = k1 * Cle * HCO30**2 * Na0 - (k2_raw * 118.7219) * Cl0 * 21**2 * Nae
    vAe2 = (Cle / (Cle + KCl)) * (HCO30 / (HCO30 + KB)) - (
        Cl0 / (Cl0 + KCl)
    ) * (21 / (21 + KB))
    vNhe1 = (Nae / (Nae + KNa)) * (H0 / (KH + H0)) - (
        Na0 / (Na0 + KNa)
    ) * (He / (KH + He))
    vBB = kp * CO20 - kn * HCO30 * H0
    JCO2 = PCO2 * (1.99997 * CO20 - CO2l - CO2e)

    GtNa = qtot0 * Nal0 / vtNa
    tNa = GtNa * vtNa
    GtK = qtot0 * Kl0 / vtK
    tK = GtK * vtK
    GCaCC = -qtot0 * Cll0 / vCaCC
    CaCC = GCaCC * vCaCC
    GBB = JCO2 / vBB
    BB = GBB * vBB
    GNhe1 = BB / vNhe1
    Nhe1 = GNhe1 * vNhe1
    alpha_NaK = ((tNa + tK) - Nkcc1) / (3 * vNaK)
    NaK = alpha_NaK * vNaK
    GCaKC = (Nkcc1 + 2 * (tNa + tK)) / (3 * vCaKC)
    CaKC = GCaKC * vCaKC
    GAe4 = (Nkcc1 - 3 * NaK + Nhe1) / vAe4
    Ae4 = GAe4 * vAe4
    GAe2 = (BB - 2 * Ae4) / vAe2
    Ae2 = vAe2 * GAe2

    dx = np.array(
        [
            tNa - qtot0 * Nal0,
            tK - qtot0 * Kl0,
            -CaCC - qtot0 * Cll0,
            qb - qa,
            (Nkcc1 - 3 * NaK + Nhe1 - Ae4 - (qb - qa) * Na0) / Hi0,
            (Nkcc1 + 2 * NaK - CaKC - (qb - qa) * K0) / Hi0,
            (2 * Nkcc1 + Ae2 + Ae4 + CaCC - (qb - qa) * Cl0) / Hi0,
            (BB - 2 * Ae4 - Ae2 - (qb - qa) * HCO30) / Hi0,
            (BB - Nhe1 - (qb - qa) * H0) / Hi0,
            (JCO2 - BB - (qb - qa) * CO20) / Hi0,
            -CaCC - (tNa + tK),
            -NaK - CaKC + (tNa + tK),
        ],
        dtype=float,
    )

    expected_par = {
        "gb": GBB,
        "pco2": PCO2,
        "kp": kp,
        "kn": kn,
        "Xi": Hi0 * (Na0 + K0 - HCO30 - Cl0 + H0),
        "a": PCO2 * 1.99997 - kp * GBB,
        "co2le": CO2l + CO2e,
        "aNaK": alpha_NaK,
        "r": r,
        "ke": Ke,
        "alpha1": alpha1,
        "s": s,
        "RTF": RTF,
        "KCaCC": KCaCC,
        "KCaKC": KCaKC,
        "eta1": eta1,
        "eta2": eta2,
        "F": F,
        "gtna": GtNa,
        "aNkcc1": alpha_Nkcc1,
        "gtk": GtK,
        "gcl": GCaCC,
        "gk": GCaKC,
        "g1": GNhe1,
        "g2": GAe2,
        "g4": GAe4,
        "b1": b1,
        "b2": b2,
        "b3": b3,
        "xl": xl,
        "k1": k1,
        "k2": k2_raw * 118.7219,
        "cle": Cle,
        "nae": Nae,
        "KCl": KCl,
        "KB": KB,
        "KNa": KNa,
        "KH": KH,
        "he": He,
        "a4": a4_raw * 1e-12,
        "a2": a2_raw * 1e-12,
        "a1": a1,
        "a3": a3,
        "Ie": Nae + Ke + Cle + HCO3e,
        "delta": delta * 1e-15,
    }
    fluxes = {
        "Ae4": Ae4,
        "Ae2": Ae2,
        "Nkcc1": Nkcc1,
        "NaK": NaK,
        "JCO2": JCO2,
        "Buffer": BB,
        "CaCC_outward": -CaCC,
        "CaKC": CaKC,
        "q_a": qa,
        "q_b": qb,
        "q_t": qt,
        "q_total": qtot0,
    }
    calibration = {
        "CO2_i": CO20,
        "H_i": H0,
        "pH_i": pHi,
        "Xi": expected_par["Xi"],
        "Psi_l": xl,
        "conductances": {
            "alpha_Nkcc1": alpha_Nkcc1,
            "PCO2": PCO2,
            "GtNa": GtNa,
            "GtK": GtK,
            "GCaCC": GCaCC,
            "GBB": GBB,
            "GNhe1": GNhe1,
            "alpha_NaK": alpha_NaK,
            "GCaKC": GCaKC,
            "GAe4": GAe4,
            "GAe2": GAe2,
        },
        "fluxes": fluxes,
        "residual_names": [f"dx{i}" for i in range(1, 13)],
        "residuals": dx.tolist(),
        "max_abs_residual": float(np.max(np.abs(dx))),
    }
    return {"calibration": calibration, "expected_par": expected_par}


def salivary7_diagnostics(
    x: np.ndarray, ca: float, par: dict[str, float], g4: float, g2: float
) -> tuple[np.ndarray, dict[str, float]]:
    """Transcription of Salivary.m (7-state outer model), before dx*10000."""
    nal, kl, hi, na, k, cl, hco = map(float, x)
    h = -(na + k - hco - cl - par["Xi"] / hi)
    co2 = (par["pco2"] * par["co2le"] - par["gb"] * par["kn"] * hco * h) / par["a"]
    cll = nal + kl
    jnak = par["aNaK"] * (
        par["r"] * par["ke"] ** 2 * na**3 * par["s"] ** 5
        / (par["s"] ** 2 * (par["ke"] ** 2 + par["alpha1"] * na**3 * par["s"]))
    )
    vcl = par["RTF"] * math.log(cll / cl)
    pcl = 1 / (1 + (par["KCaCC"] / ca) ** par["eta1"])
    vk = par["RTF"] * math.log(par["ke"] / k)
    pk = 1 / (1 + (par["KCaKC"] / ca) ** par["eta2"])
    vtna = par["RTF"] * math.log(nal / par["nae"])
    vtk = par["RTF"] * math.log(kl / par["ke"])
    denom = (par["gtna"] + par["gtk"]) * (
        par["gk"] * pk + par["gcl"] * pcl
    ) + (par["gk"] * pk) * (par["gcl"] * pcl)
    va = -(
        (par["gtna"] + par["gtk"]) * (par["F"] * jnak)
        - (par["gk"] * pk)
        * (par["gtna"] * vtna + par["gtk"] * vtk + (par["gtna"] + par["gtk"]) * vk)
        + (par["gcl"] * pcl * vcl)
        * ((par["gk"] * pk) + (par["gtna"] + par["gtk"]))
    ) / denom
    vb = -(
        ((par["gtna"] + par["gtk"]) + (par["gcl"] * pcl))
        * ((par["F"] * jnak) - (par["gk"] * pk) * vk)
        + (par["gcl"] * pcl) * (par["gtna"] * vtna + par["gtk"] * vtk)
        + (par["gtna"] + par["gtk"]) * (par["gcl"] * pcl * vcl)
    ) / denom
    vt = va - vb
    jcl = par["gcl"] * pcl * (va + vcl) / par["F"]
    jk = par["gk"] * pk * (vb - vk) / par["F"]
    jtna = par["gtna"] * (vt - vtna) / par["F"]
    jtk = par["gtk"] * (vt - vtk) / par["F"]
    qa = par["b1"] * (2 * (nal + kl - na - k - h) - co2 + par["xl"])
    qb = par["b2"] * (2 * (na + k + h) + co2 - par["Ie"])
    qt = par["b3"] * (2 * (nal + kl) + par["xl"] - par["Ie"])
    qtotal = qa + qt
    jnk1 = par["aNkcc1"] * (
        (par["a1"] - par["a2"] * na * k * cl**2)
        / (par["a3"] + par["a4"] * na * k * cl**2)
    )
    j4 = g4 * (par["k1"] * par["cle"] * hco**2 * na - par["k2"] * cl * 21**2 * par["nae"])
    j2 = g2 * (
        (par["cle"] / (par["cle"] + par["KCl"])) * (hco / (hco + par["KB"]))
        - (cl / (cl + par["KCl"])) * (21 / (21 + par["KB"]))
    )
    j1 = par["g1"] * (
        (par["nae"] / (par["nae"] + par["KNa"])) * (h / (par["KH"] + h))
        - (na / (na + par["KNa"])) * (par["he"] / (par["KH"] + par["he"]))
    )
    jb = par["gb"] * (par["kp"] * co2 - par["kn"] * hco * h)
    volume_rate = qb - qa
    raw = np.array(
        [
            jtna - qtotal * nal,
            jtk - qtotal * kl,
            volume_rate,
            (jnk1 - 3 * jnak + j1 - j4 - volume_rate * na) / hi,
            (jnk1 + 2 * jnak - jk - volume_rate * k) / hi,
            (2 * jnk1 + j2 + j4 + jcl - volume_rate * cl) / hi,
            (jb - 2 * j4 - j2 - volume_rate * hco) / hi,
        ],
        dtype=float,
    )
    diagnostics = {
        "Ca": ca,
        "H": h,
        "pH": math.log(1000 / h) / (math.log(2) + math.log(5)) if h > 0 else math.nan,
        "CO2": co2,
        "cell_volume_pL_from_delta_Hi": par["delta"] * hi * 1e12,
        "Cll": cll,
        "Va": va,
        "Vb": vb,
        "Vt": vt,
        "q_a": qa,
        "q_b": qb,
        "q_t": qt,
        "q_total": qtotal,
        "J_NaK": jnak,
        "J_NKCC1": jnk1,
        "J_Ae4": j4,
        "J_Ae2": j2,
        "J_NHE1": j1,
        "J_buffer": jb,
        "J_Cl_code": jcl,
        "J_K_code": jk,
        "J_tNa": jtna,
        "J_tK": jtk,
    }
    return raw, diagnostics


def salivary5_diagnostics(
    x: np.ndarray, ca: float, par: dict[str, float]
) -> tuple[np.ndarray, dict[str, float]]:
    """Transcription of Salivary2.m, before its dx*1000 multiplier."""
    hi, na, k, cl, hco = map(float, x)
    nal, kl = 118.7, 5.6
    cll = nal + kl
    h = -(na + k - hco - cl - par["Xi"] / hi)
    co2 = (par["pco2"] * par["co2le"] - par["gb"] * par["kn"] * hco * h) / par["a"]
    jnak = par["aNaK"] * (
        par["r"] * par["ke"] ** 2 * na**3 * par["s"] ** 5
        / (par["s"] ** 2 * (par["ke"] ** 2 + par["alpha1"] * na**3 * par["s"]))
    )
    vcl = par["RTF"] * math.log(cll / cl)
    pcl = 1 / (1 + (par["KCaCC"] / ca) ** par["eta1"])
    vk = par["RTF"] * math.log(par["ke"] / k)
    pk = 1 / (1 + (par["KCaKC"] / ca) ** par["eta2"])
    vtna = par["RTF"] * math.log(nal / par["nae"])
    vtk = par["RTF"] * math.log(kl / par["ke"])
    denom = (par["gtna"] + par["gtk"]) * (par["gk"] * pk + par["gcl"] * pcl) + (
        par["gk"] * pk
    ) * (par["gcl"] * pcl)
    va = -(
        (par["gtna"] + par["gtk"]) * (par["F"] * jnak)
        - (par["gk"] * pk)
        * (par["gtna"] * vtna + par["gtk"] * vtk + (par["gtna"] + par["gtk"]) * vk)
        + (par["gcl"] * pcl * vcl)
        * ((par["gk"] * pk) + (par["gtna"] + par["gtk"]))
    ) / denom
    vb = -(
        ((par["gtna"] + par["gtk"]) + (par["gcl"] * pcl))
        * ((par["F"] * jnak) - (par["gk"] * pk) * vk)
        + (par["gcl"] * pcl) * (par["gtna"] * vtna + par["gtk"] * vtk)
        + (par["gtna"] + par["gtk"]) * (par["gcl"] * pcl * vcl)
    ) / denom
    jcl = par["gcl"] * pcl * (va + vcl) / par["F"]
    jk = par["gk"] * pk * (vb - vk) / par["F"]
    qa = par["b1"] * (2 * (nal + kl - na - k - h) - co2 + par["xl"])
    qb = par["b2"] * (2 * (na + k + h) + co2 - par["Ie"])
    qt = par["b3"] * (2 * (nal + kl) + par["xl"] - par["Ie"])
    qtotal = qa + qt
    jnk1 = par["aNkcc1"] * (
        (par["a1"] - par["a2"] * na * k * cl**2)
        / (par["a3"] + par["a4"] * na * k * cl**2)
    )
    j4 = par["g4"] * (
        par["k1"] * par["cle"] * hco**2 * na - par["k2"] * cl * 21**2 * par["nae"]
    )
    j2 = par["g2"] * (
        (par["cle"] / (par["cle"] + par["KCl"])) * (hco / (hco + par["KB"]))
        - (cl / (cl + par["KCl"])) * (21 / (21 + par["KB"]))
    )
    j1 = par["g1"] * (
        (par["nae"] / (par["nae"] + par["KNa"])) * (h / (par["KH"] + h))
        - (na / (na + par["KNa"])) * (par["he"] / (par["KH"] + par["he"]))
    )
    jb = par["gb"] * (par["kp"] * co2 - par["kn"] * hco * h)
    volume_rate = qb - qa
    raw = np.array(
        [
            volume_rate,
            (jnk1 - 3 * jnak + j1 - j4 - volume_rate * na) / hi,
            (jnk1 + 2 * jnak - jk - volume_rate * k) / hi,
            (2 * jnk1 + j2 + j4 + jcl - volume_rate * cl) / hi,
            (jb - 2 * j4 - j2 - volume_rate * hco) / hi,
        ],
        dtype=float,
    )
    diagnostics = {
        "Ca": ca,
        "H": h,
        "pH": math.log(1000 / h) / (math.log(2) + math.log(5)) if h > 0 else math.nan,
        "CO2": co2,
        "cell_volume_pL_from_delta_Hi": par["delta"] * hi * 1e12,
        "Cll": cll,
        "Va": va,
        "Vb": vb,
        "Vt": va - vb,
        "q_a": qa,
        "q_b": qb,
        "q_t": qt,
        "q_total": qtotal,
        "J_NaK": jnak,
        "J_NKCC1": jnk1,
        "J_Ae4": j4,
        "J_Ae2": j2,
        "J_NHE1": j1,
        "J_buffer": jb,
        "J_Cl_code": jcl,
        "J_K_code": jk,
    }
    return raw, diagnostics


def scenario_par(par: dict[str, float], scenario: str) -> tuple[dict[str, float], float, float]:
    result = dict(par)
    g4, g2 = par["g4"], par["g2"]
    if scenario == "ae2_knockout":
        g2 = 0.0
        result["g2"] = 0.0
    elif scenario == "ae4_knockout":
        g4 = 0.0
        result["g4"] = 0.0
    elif scenario != "wt":
        raise ValueError(scenario)
    return result, g4, g2


def integrate_two_stage(
    raw_function: Callable[[np.ndarray, float], tuple[np.ndarray, dict[str, float]]],
    y0: np.ndarray,
    scale: float,
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    """Integrate the archived 0.05 -> 0.55 calcium step at t=100."""
    segments: list[tuple[np.ndarray, np.ndarray]] = []
    metadata: list[dict[str, Any]] = []
    current = np.asarray(y0, dtype=float)
    for index, (start, end, ca) in enumerate(((0.0, 100.0, 0.05), (100.0, 200.0, 0.55))):
        times = np.arange(start, end + 0.5, 1.0)

        def rhs(_t: float, state: np.ndarray) -> np.ndarray:
            raw, _ = raw_function(state, ca)
            return scale * raw

        solution = solve_ivp(
            rhs,
            (start, end),
            current,
            method="BDF",
            t_eval=times,
            rtol=1e-6,
            atol=1e-6,
        )
        if not solution.success:
            raise RuntimeError(f"BDF segment {index} failed: {solution.message}")
        segment_t, segment_y = solution.t, solution.y.T
        if index:
            segment_t, segment_y = segment_t[1:], segment_y[1:]
        segments.append((segment_t, segment_y))
        current = solution.y[:, -1]
        metadata.append(
            {
                "segment": index + 1,
                "t_span": [start, end],
                "calcium_numeric": ca,
                "nfev": solution.nfev,
                "njev": solution.njev,
                "nlu": solution.nlu,
                "message": solution.message,
            }
        )
    return (
        np.concatenate([part[0] for part in segments]),
        np.concatenate([part[1] for part in segments]),
        {"method": "scipy.solve_ivp/BDF", "rtol": 1e-6, "atol": 1e-6, "segments": metadata},
    )


def make_forensic_copy() -> dict[str, Any]:
    """Make a disposable read-only-source working copy outside archive/."""
    FORENSIC_COPY.mkdir(parents=True, exist_ok=True)
    copied: list[dict[str, Any]] = []
    for relative in (
        Path("Original/Parameters.m"),
        Path("Original/Saliva_Ae4.m"),
        Path("Par.m"),
        Path("Par.mat"),
        Path("Salivary.m"),
        Path("Salivary2.m"),
        Path("Salivary_ex.m"),
    ):
        source = ARCHIVE_PROJECT / relative
        target = FORENSIC_COPY / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            if sha256(target) != sha256(source):
                raise RuntimeError(f"existing forensic copy differs from archive source: {target}")
        else:
            shutil.copy2(source, target)
        copied.append(
            {
                "path": str(relative),
                "source_sha256": sha256(source),
                "copy_sha256": sha256(target),
                "bytes": source.stat().st_size,
            }
        )
    return {"path": str(FORENSIC_COPY.relative_to(ROOT)), "files": copied}


def write_outputs() -> dict[str, Any]:
    RESULTS.mkdir(parents=True, exist_ok=True)
    par, dtypes = load_par_mat()
    original = original_calibration()
    expected = original["expected_par"]
    comparisons: dict[str, dict[str, float]] = {}
    for field in sorted(par):
        expected_value = expected[field]
        absolute = par[field] - expected_value
        comparisons[field] = {
            "par_mat": par[field],
            "original_calibration_or_exact_transform": expected_value,
            "absolute_difference": absolute,
            "relative_difference": absolute / expected_value if expected_value else absolute,
        }

    copy_record = make_forensic_copy()
    trajectories: list[dict[str, Any]] = []
    residual_rows: list[dict[str, Any]] = []
    run_summaries: dict[str, Any] = {}
    initial_7 = np.array([118.7, 5.6, 28.7, 25.0, 120.0, 50.0, 10.0])
    initial_5 = np.array([28.7, 25.0, 120.0, 50.0, 10.0])

    for version in ("salivary_7_state", "salivary2_5_state"):
        state_names = STATE_NAMES_7 if version == "salivary_7_state" else STATE_NAMES_5
        y0 = initial_7 if version == "salivary_7_state" else initial_5
        scale = 10000.0 if version == "salivary_7_state" else 1000.0
        run_summaries[version] = {}
        for scenario in SCENARIOS:
            run_par, g4, g2 = scenario_par(par, scenario)
            if version == "salivary_7_state":
                evaluator = lambda state, ca, rp=run_par, v4=g4, v2=g2: salivary7_diagnostics(
                    state, ca, rp, v4, v2
                )
            else:
                evaluator = lambda state, ca, rp=run_par: salivary5_diagnostics(state, ca, rp)
            raw0, diag0 = evaluator(y0, 0.05)
            times, states, solver = integrate_two_stage(evaluator, y0, scale)
            diag_series: list[dict[str, float]] = []
            for t, state in zip(times, states, strict=True):
                ca = 0.05 if t <= 100 else 0.55
                _, diag = evaluator(state, ca)
                diag_series.append(diag)
                row: dict[str, Any] = {"version": version, "scenario": scenario, "time": float(t)}
                row.update({name: float(value) for name, value in zip(state_names, state, strict=True)})
                row.update({name: float(value) for name, value in diag.items()})
                trajectories.append(row)
            for index, value in enumerate(raw0):
                residual_rows.append(
                    {
                        "version": version,
                        "scenario": scenario,
                        "state_or_balance": state_names[index],
                        "raw_rhs_before_hidden_scale": float(value),
                        "rhs_after_hidden_scale": float(value * scale),
                    }
                )
            pre_index = int(np.where(times == 100.0)[0][0])
            end_index = -1
            pre_raw, _ = evaluator(states[pre_index], 0.05)
            end_raw, _ = evaluator(states[end_index], 0.55)
            q_values = np.array([entry["q_total"] for entry in diag_series])
            stimulated = times > 100
            peak_index = int(np.flatnonzero(stimulated)[np.argmax(q_values[stimulated])])
            run_summaries[version][scenario] = {
                "activities": {"g2": g2, "g4": g4},
                "hidden_rhs_scale": scale,
                "initial_state": dict(zip(state_names, map(float, y0), strict=True)),
                "initial_raw_rhs": dict(zip(state_names, map(float, raw0), strict=True)),
                "initial_scaled_rhs": dict(zip(state_names, map(float, raw0 * scale), strict=True)),
                "initial_diagnostics": diag0,
                "pre_stimulus_t100": {
                    "state": dict(zip(state_names, map(float, states[pre_index]), strict=True)),
                    "diagnostics": diag_series[pre_index],
                    "raw_rhs": dict(zip(state_names, map(float, pre_raw), strict=True)),
                    "max_abs_raw_rhs": float(np.max(np.abs(pre_raw))),
                    "max_abs_scaled_rhs": float(scale * np.max(np.abs(pre_raw))),
                },
                "stimulated_endpoint_t200": {
                    "state": dict(zip(state_names, map(float, states[end_index]), strict=True)),
                    "diagnostics": diag_series[end_index],
                    "raw_rhs": dict(zip(state_names, map(float, end_raw), strict=True)),
                    "max_abs_raw_rhs": float(np.max(np.abs(end_raw))),
                    "max_abs_scaled_rhs": float(scale * np.max(np.abs(end_raw))),
                },
                "stimulated_peak_sample": {
                    "time": float(times[peak_index]),
                    "q_total": float(q_values[peak_index]),
                },
                "solver": solver,
            }

    for version, version_runs in run_summaries.items():
        wt_end = version_runs["wt"]["stimulated_endpoint_t200"]["diagnostics"]["q_total"]
        wt_peak = version_runs["wt"]["stimulated_peak_sample"]["q_total"]
        for scenario, run in version_runs.items():
            end = run["stimulated_endpoint_t200"]["diagnostics"]["q_total"]
            peak = run["stimulated_peak_sample"]["q_total"]
            run["comparison_to_wt"] = {
                "endpoint_flow_ratio": end / wt_end,
                "endpoint_percent_change": 100 * (end / wt_end - 1),
                "peak_sample_flow_ratio": peak / wt_peak,
                "peak_sample_percent_change": 100 * (peak / wt_peak - 1),
            }

    summary: dict[str, Any] = {
        "execution_status": {
            "matlab_available": shutil.which("matlab") is not None,
            "octave_available": shutil.which("octave") is not None,
            "literal_matlab_execution": False,
            "reason": "Neither MATLAB nor Octave is installed; the only archived Salivary.m also has an arity mismatch with Original/Saliva_Ae4.m.",
            "equation_evaluation": "Direct Python transcription of immutable MATLAB arithmetic; scipy BDF substitutes for ode15s.",
        },
        "archive_working_copy": copy_record,
        "par_mat_sha256": sha256(ARCHIVE_PROJECT / "Par.mat"),
        "par_mat_fields": {
            name: {"value": par[name], "matlab_storage_dtype": dtypes[name]} for name in par
        },
        "original_calibration": original["calibration"],
        "original_to_par_concordance": {
            "fields_compared": len(comparisons),
            "max_abs_difference": max(abs(item["absolute_difference"]) for item in comparisons.values()),
            "max_abs_relative_difference": max(abs(item["relative_difference"]) for item in comparisons.values()),
            "fields": comparisons,
        },
        "runs": run_summaries,
        "published_reference_targets": {
            "resting_state": {
                "Na_i_mM": 25.0,
                "K_i_mM": 120.0,
                "Cl_i_mM": 50.0,
                "HCO3_i_mM": 12.1,
                "pH_i": 6.91,
                "CO2_i_mM": 6.6,
                "Na_l_mM": 118.7,
                "K_l_mM": 5.6,
                "Cl_l_mM": 124.3,
                "Va_mV": -50.24,
                "Vb_mV": -62.8,
                "cell_volume_pL": 1.3,
            },
            "phenotype": {
                "ae4_knockout_flow_percent_change_approx": -24.0,
                "ae2_knockout_flow": "no significant change",
            },
        },
    }

    summary_path = RESULTS / "reproduction_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    trajectory_path = RESULTS / "reproduction_trajectories.csv"
    fieldnames = sorted({key for row in trajectories for key in row})
    preferred = ["version", "scenario", "time"]
    fieldnames = preferred + [field for field in fieldnames if field not in preferred]
    with trajectory_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(trajectories)

    residual_path = RESULTS / "reproduction_residuals.csv"
    with residual_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(residual_rows[0]))
        writer.writeheader()
        writer.writerows(residual_rows)

    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="print key regression quantities")
    args = parser.parse_args()
    summary = write_outputs()
    if args.check:
        concise = {
            "original_max_abs_residual": summary["original_calibration"]["max_abs_residual"],
            "par_max_abs_difference": summary["original_to_par_concordance"]["max_abs_difference"],
            "seven_state_endpoint_flow_ratios": {
                scenario: run["comparison_to_wt"]["endpoint_flow_ratio"]
                for scenario, run in summary["runs"]["salivary_7_state"].items()
            },
            "five_state_endpoint_flow_ratios": {
                scenario: run["comparison_to_wt"]["endpoint_flow_ratio"]
                for scenario, run in summary["runs"]["salivary2_5_state"].items()
            },
        }
        print(json.dumps(concise, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
