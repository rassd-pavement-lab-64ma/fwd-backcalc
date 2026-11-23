#!/usr/bin/env python3
"""
fwd_backcalc.py - FWD backcalculation of pavement layer moduli.

Forward model: exact Burmister / Schiffman multilayer elastic solution
(Hankel-transform domain, per-layer A,B,C,D constants, solved as a
banded linear system per transform variable xi, numerically integrated
with Gauss-Legendre between Bessel zeros). Same theory class as
WESLEA/CHEVRON (the forward engines inside EVERCALC/CHEVDEF).

Backcalculation: bounded least-squares minimization of the sum of
squared RELATIVE deflection errors (EVERCALC's objective), solved in
log10(modulus) space with scipy TRF + multi-start, per the guidelines
of FHWA-HRT-16-011 Chapter 3.

Usage:
  python3 fwd_backcalc.py input.json

Input JSON:
{
  "load_kN": 40.0,
  "plate_radius_mm": 150.0,
  "sensors_mm": [0, 203, 305, 610, 914, 1219],
  "deflections_um": [467, 394, 350, 193, 127, 97],
  "layers": [
     {"thickness_mm": 135, "poisson": 0.35, "E_seed_MPa": 2700,
      "E_min_MPa": 690, "E_max_MPa": 13800},
     {"thickness_mm": 457, "poisson": 0.40, "E_seed_MPa": 170,
      "E_min_MPa": 35, "E_max_MPa": 3500},
     {"thickness_mm": null, "poisson": 0.45, "E_seed_MPa": 100,
      "E_min_MPa": 35, "E_max_MPa": 3500}
  ],
  "fix": []            // optional: layer indices (0-based) to hold at seed
}

Output JSON: backcalculated moduli, RMS%%, per-sensor fit, iterations.
"""

import json
import sys
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import j0, j1, jn_zeros
from scipy.optimize import least_squares


# ----------------------------------------------------------------------
# Forward model: multilayer elastic surface deflections
# ----------------------------------------------------------------------

def _solve_U0(xi, h, nu, G):
    """
    Batched solution of the layer constants for transform variable(s) xi
    with unit transformed surface pressure (p_bar = 1).

    h  : array (n-1,) finite layer thicknesses (same units as 1/xi)
    nu : array (n,) Poisson ratios
    G  : array (n,) shear moduli
    Returns U0(xi) : transformed vertical surface displacement, shape xi.shape
    """
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    m = xi.shape[0]
    n = len(nu)
    nvar = 4 * (n - 1) + 2
    A = np.zeros((m, nvar, nvar))
    b = np.zeros((m, nvar))

    def idx(i):
        # variable block for layer i (0-based); last layer has 2 vars
        return 4 * i if i < n - 1 else 4 * (n - 1)

    row = 0
    # --- surface boundary conditions (layer 0, z = 0) ---
    i = 0
    e = np.exp(-xi * h[0]) if n > 1 else None
    c0 = 2 if n > 1 else 1   # index of C0 (2-var block for infinite layer)
    A[:, row, 0] = xi * xi                      # xi * (xi*A0)
    A[:, row, c0] = -2.0 * nu[0] * xi           # xi * (-2 nu C0)
    if n > 1:
        A[:, row, 1] = xi * xi * e              # xi * (xi*B0') e^{-xi h}
        A[:, row, 3] = 2.0 * nu[0] * xi * e     # xi * (2 nu D0') e^{-xi h}
    row += 1
    # sigmaz(0) = -1 : xi^2[(xi A0 + (1-2nu) C0) + (-xi B0' + (1-2nu) D0') e^{-xi h}]
    A[:, row, 0] = xi ** 3
    A[:, row, c0] = xi ** 2 * (1 - 2 * nu[0])
    if n > 1:
        A[:, row, 1] = -xi ** 3 * e
        A[:, row, 3] = xi ** 2 * (1 - 2 * nu[0]) * e
    b[:, row] = -1.0
    row += 1

    # --- interface continuity at bottom of layer k (z = h_k) ---
    for k in range(n - 1):
        et = np.exp(-xi * h[k])          # e^{-xi h_k} for A,C part at z=h
        # top of layer k+1: z = 0
        if k + 1 < n - 1:
            eb = np.exp(-xi * h[k + 1])  # e^{-xi h_{k+1}} for B',D' part at z=0
        nuk, nuk1 = nu[k], nu[k + 1]
        Gk, Gk1 = G[k], G[k + 1]
        ik, ik1 = idx(k), idx(k + 1)
        c1 = 2 if k + 1 < n - 1 else 1   # C index of layer k+1 block

        # sigma_z continuity: layer k at z=h_k  minus  layer k+1 at z=0
        A[:, row, ik + 0] = xi ** 3 * et
        A[:, row, ik + 2] = xi ** 2 * ((1 - 2 * nuk) + xi * h[k]) * et
        A[:, row, ik + 1] = -xi ** 3
        A[:, row, ik + 3] = xi ** 2 * ((1 - 2 * nuk) - xi * h[k])
        A[:, row, ik1 + 0] = -xi ** 3
        A[:, row, ik1 + c1] = -xi ** 2 * (1 - 2 * nuk1)
        if k + 1 < n - 1:
            A[:, row, ik1 + 1] = xi ** 3 * eb
            A[:, row, ik1 + 3] = -xi ** 2 * (1 - 2 * nuk1) * eb
        row += 1

        # tau continuity
        A[:, row, ik + 0] = xi ** 2 * et
        A[:, row, ik + 2] = xi * (-2 * nuk + xi * h[k]) * et
        A[:, row, ik + 1] = xi ** 2
        A[:, row, ik + 3] = xi * (2 * nuk + xi * h[k])
        A[:, row, ik1 + 0] = -xi ** 2
        A[:, row, ik1 + c1] = -xi * (-2 * nuk1)
        if k + 1 < n - 1:
            A[:, row, ik1 + 1] = -xi ** 2 * eb
            A[:, row, ik1 + 3] = -xi * (2 * nuk1) * eb
        row += 1

        # 2G u_z continuity: (1/G_k) * U_k = (1/G_{k+1}) * U_{k+1}
        # U_k(h) = -xi[(xi A + (2-4nu)C + xi C h) et + (xi B' + (-(2-4nu))D' + xi D' h)]
        A[:, row, ik + 0] = -xi ** 2 * et / Gk
        A[:, row, ik + 2] = -xi * ((2 - 4 * nuk) + xi * h[k]) * et / Gk
        A[:, row, ik + 1] = -xi ** 2 / Gk
        A[:, row, ik + 3] = -xi * (-(2 - 4 * nuk) + xi * h[k]) / Gk
        A[:, row, ik1 + 0] = xi ** 2 / Gk1
        A[:, row, ik1 + c1] = xi * (2 - 4 * nuk1) / Gk1
        if k + 1 < n - 1:
            A[:, row, ik1 + 1] = xi ** 2 * eb / Gk1
            A[:, row, ik1 + 3] = xi * (-(2 - 4 * nuk1)) * eb / Gk1
        row += 1

        # 2G u_r continuity: (1/G_k) * V_k = (1/G_{k+1}) * V_{k+1}
        # V_k(h) = xi[(-xi A + C - xi C h) et + (xi B' + D' + xi D' h)]
        A[:, row, ik + 0] = -xi ** 2 * et / Gk
        A[:, row, ik + 2] = xi * (1 - xi * h[k]) * et / Gk
        A[:, row, ik + 1] = xi ** 2 / Gk
        A[:, row, ik + 3] = xi * (1 + xi * h[k]) / Gk
        A[:, row, ik1 + 0] = xi ** 2 / Gk1
        A[:, row, ik1 + c1] = -xi / Gk1
        if k + 1 < n - 1:
            A[:, row, ik1 + 1] = -xi ** 2 * eb / Gk1
            A[:, row, ik1 + 3] = -xi * eb / Gk1
        row += 1

    # row scaling for conditioning, then solve
    scale = np.abs(A).max(axis=2, keepdims=True)
    scale[scale == 0] = 1.0
    xs = np.linalg.solve(A / scale, (b / scale[:, :, 0])[..., None])[:, :, 0]

    # U0 = surface vertical displacement transform (unit p_bar):
    U0 = -xi * (xi * xs[:, 0] + (2 - 4 * nu[0]) * xs[:, c0])
    if n > 1:
        U0 = U0 - xi * (xi * xs[:, 1] - (2 - 4 * nu[0]) * xs[:, 3]) * e
    U0 = U0 / (2 * G[0])
    return U0


def surface_deflections(r, h, E, nu, q, a, t_max=220.0, n_gauss=8):
    """
    Surface deflections w(r) [same length unit as h, a] under uniform
    pressure q over circular plate radius a.
    r : array of radial offsets.  h : (n-1,) thicknesses, E, nu : (n,).
    Integral: w(r) = q * ∫ U0(xi) J1(xi a) J0(xi r) dxi
    Integrated in t = xi*a with Gauss-Legendre over sub-segments sized to
    the fastest Bessel oscillation.
    """
    h = np.asarray(h, float)
    E = np.asarray(E, float)
    nu = np.asarray(nu, float)
    G = E / (2 * (1 + nu))
    r = np.atleast_1d(np.asarray(r, float))
    w = np.zeros_like(r)

    xg, wg = leggauss(n_gauss)
    for j, rj in enumerate(r):
        # segment length: half-period of fastest oscillation / 2
        osc = max(1.0, rj / a)
        seg = np.pi / osc / 2.0
        nseg = int(np.ceil(t_max / seg))
        edges = np.linspace(0.0, nseg * seg, nseg + 1)
        tm = 0.5 * (edges[:-1] + edges[1:])   # (nseg,)
        th = 0.5 * (edges[1:] - edges[:-1])
        t = (tm[:, None] + th[:, None] * xg[None, :]).ravel()
        t = np.maximum(t, 1e-6)
        U = _solve_U0(t / a, h, nu, G)
        integrand = U * j1(t) * j0(t * rj / a)
        w[j] = q * float(np.sum(integrand.reshape(nseg, n_gauss) * (th[:, None] * wg[None, :])))
    return w


# ----------------------------------------------------------------------
# Backcalculation
# ----------------------------------------------------------------------

def backcalculate(defl_um, sensors_mm, h_mm, poisson, E_seed, E_min, E_max,
                  q_kPa, a_mm, fix=None, n_starts=5, verbose=False):
    """
    Minimize sum of squared relative deflection errors over log10(E).
    defl_um: measured deflections (microns)
    sensors_mm, h_mm: geometry (mm)
    E_seed/E_min/E_max: MPa per layer
    fix: indices of layers to hold at seed value
    Returns dict with moduli, RMS %, fit details.
    """
    defl_um = np.asarray(defl_um, float)
    sensors_mm = np.asarray(sensors_mm, float)
    h_mm = np.asarray(h_mm, float)
    poisson = np.asarray(poisson, float)
    E_seed = np.asarray(E_seed, float)
    E_min = np.asarray(E_min, float)
    E_max = np.asarray(E_max, float)
    n = len(poisson)
    fix = set(fix or [])

    free = [i for i in range(n) if i not in fix]
    q_Pa = q_kPa * 1e3
    a_m = a_mm * 1e-3
    h_m = h_mm * 1e-3
    r_m = sensors_mm * 1e-3

    def calc_defl(logE_free):
        E = E_seed.copy()
        for j, i in enumerate(free):
            E[i] = 10.0 ** logE_free[j]
        w_m = surface_deflections(r_m, h_m, E * 1e6, poisson, q_Pa, a_m)
        return w_m * 1e6  # microns

    def resid(logE_free):
        c = calc_defl(logE_free)
        return (c - defl_um) / defl_um

    lb = np.array([np.log10(E_min[i]) for i in free])
    ub = np.array([np.log10(E_max[i]) for i in free])
    seed = np.array([np.log10(E_seed[i]) for i in free])

    # multi-start: seed + low/high blends
    starts = [seed]
    rng = np.random.default_rng(20240117)
    for _ in range(n_starts - 1):
        starts.append(lb + (ub - lb) * rng.random(len(free)))

    best = None
    for s in starts:
        s = np.clip(s, lb, ub)
        try:
            sol = least_squares(resid, s, bounds=(lb, ub), method='trf',
                                x_scale='jac', ftol=1e-12, xtol=1e-12,
                                gtol=1e-12, max_nfev=400)
        except Exception:
            continue
        cost = 2 * sol.cost
        if best is None or cost < best[0]:
            best = (cost, sol)
    if best is None:
        raise RuntimeError("optimization failed for all starts")
    cost, sol = best
    E_out = E_seed.copy()
    for j, i in enumerate(free):
        E_out[i] = 10.0 ** sol.x[j]
    calc = calc_defl(sol.x)
    rms_pct = 100.0 * np.sqrt(np.mean(((calc - defl_um) / defl_um) ** 2))
    return {
        "E_MPa": list(np.round(E_out, 1)),
        "rms_percent": round(float(rms_pct), 3),
        "calc_defl_um": list(np.round(calc, 2)),
        "meas_defl_um": list(defl_um),
        "rel_error_percent": list(np.round(100 * (calc - defl_um) / defl_um, 2)),
        "n_function_evals": int(sol.nfev),
        "optimality": float(sol.optimality),
    }


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------

def main():
    cfg = json.load(open(sys.argv[1]))
    q_kPa = cfg["load_kN"] * 1e3 / (np.pi * (cfg["plate_radius_mm"] * 1e-3) ** 2) / 1e3
    layers = cfg["layers"]
    h_mm, poisson, E_seed, E_min, E_max = [], [], [], [], []
    for L in layers:
        th = L.get("thickness_mm")
        h_mm.append(th if th is not None else 0.0)
        poisson.append(L["poisson"])
        E_seed.append(L.get("E_seed_MPa", 0.5 * (L["E_min_MPa"] + L["E_max_MPa"])))
        E_min.append(L["E_min_MPa"])
        E_max.append(L["E_max_MPa"])
    h_mm = h_mm[:-1]  # last layer is semi-infinite

    res = backcalculate(cfg["deflections_um"], cfg["sensors_mm"], h_mm,
                        poisson, E_seed, E_min, E_max, q_kPa,
                        cfg["plate_radius_mm"], fix=cfg.get("fix"))
    res["load_kPa"] = round(q_kPa, 1)
    json.dump(res, sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
