#!/usr/bin/env python3
"""Checks for the theory of "When canonical covariance localisation fails".

Sections refer to manuscript.tex. Check names keep the labels of the 1 November 2026
draft: "Lemma 1" = Lemma 3.1, "Thm2" = Theorem 3.3, "Prop 4" = Proposition 4.1,
"Table 1" = Table 3.
  X2  inflation (Lemmas 3.5-3.6, Theorem 3.7): derivative vs finite differences, sign at
      tau*, existence / non-existence of a finite canonical minimiser for N=2, Stein
      identity at a*beta = 3, single sign change of Psi',
  X3  vector observations (Proposition 5.1): first variation vs central finite differences
      of the exact conditional risk (tensor quadrature); sign of Phi for two observed
      directions (Conjecture 5.4, numerical only),
  X4  W=I curvature ratio grid (Open problem 8.2, numerical only),
  X5  aligned-example curvature R''(0) against N (Section 4.2, Remark 4.2),
  MC  seeded Monte Carlo of the REAL ensemble procedure (actual ensembles, sample
      covariance, localisation in a rotated basis, update, physical loss),
  negative controls that must be detected.

No assert statements are used: every check appends to an explicit failure list, the
status is computed, results are written as JSON to a NEW directory given by --out,
and the process exits with status 1 on any failure.

Run:  uv run --with numpy==2.3.5 --with scipy==1.17.0 --with mpmath python3 -I -B theory_checks.py --out DIR
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import argparse
import hashlib
import json
import math
import platform
import sys
import time

import numpy as np
import scipy
from scipy import integrate, optimize
from scipy.special import gammaln, exp1
import mpmath

FAILURES = []
CHECKS = []


def check(name, ok, **info):
    ok = bool(ok)
    rec = {"name": name, "ok": ok}
    for k, v in info.items():
        rec[k] = _js(v)
    CHECKS.append(rec)
    if not ok:
        FAILURES.append(name)
    return ok


def _js(v):
    if isinstance(v, (np.floating,)):
        return float(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, np.ndarray):
        return [_js(x) for x in v.tolist()]
    if isinstance(v, (list, tuple)):
        return [_js(x) for x in v]
    if isinstance(v, dict):
        return {str(k): _js(x) for k, x in v.items()}
    if isinstance(v, float) and not math.isfinite(v):
        return str(v)
    return v


# ---------------------------------------------------------------- expectations
def pdf_gamma(d, nu):
    k = nu / 2.0
    return np.exp((k - 1) * np.log(d) - k * d + k * np.log(k) - gammaln(k))


def E1(f, nu, scale=1.0):
    """E f(D), D ~ chi2_nu/nu, adaptive quadrature. `scale` is a small length scale
    of f near 0 (e.g. 1/(beta*tau)), used only to place break points."""
    if nu == 1:
        # D = s^2, density of s on (0,inf): 2 (2 pi)^(-1/2) exp(-s^2/2)
        c = 2.0 / math.sqrt(2 * math.pi)
        g = lambda s: f(s * s) * c * math.exp(-0.5 * s * s)
        rs = math.sqrt(scale)
        pts = sorted(set([rs * 10.0 ** q for q in range(0, 40) if rs * 10.0 ** q < 0.5]
                         + [1.0, 3.0]))
        tot = 0.0
        lo = 0.0
        for p in pts + [12.0]:
            tot += integrate.quad(g, lo, p, limit=400, epsabs=1e-14, epsrel=1e-12)[0]
            lo = p
        return tot
    sd = math.sqrt(2.0 / nu)
    hi = 1 + 40 * sd + 40
    geo = [scale * 10.0 ** q for q in range(0, 40) if scale * 10.0 ** q < 0.5]
    pts = sorted(set([x for x in geo + [max(1 - 6 * sd, 1e-3), 1.0, 1 + 6 * sd]
                      if 0 < x < hi]))
    g = lambda d: f(d) * pdf_gamma(d, nu)
    tot = 0.0
    lo = 0.0
    for p in pts + [hi]:
        tot += integrate.quad(g, lo, p, limit=400, epsabs=1e-15, epsrel=1e-12)[0]
        lo = p
    return tot


_NODE_CACHE = {}


def exp_nodes(nu, h):
    """Trapezoidal rule in u = log d for E f(D), D ~ chi2_nu/nu.  The transformed
    integrand is analytic in the strip |Im u| < pi/2 for the rational integrands used
    here (poles of 1/(1+c d) sit at Im u = pi), so the error is O(exp(-pi^2/h))."""
    key = (nu, h)
    if key not in _NODE_CACHE:
        k = nu / 2.0
        lo = -(37.0 + abs(k * math.log(k) - gammaln(k))) / k
        u = np.arange(lo, 5.0 + h / 2, h)
        d = np.exp(u)
        w = h * np.exp(k * u - k * d + k * math.log(k) - gammaln(k))
        _NODE_CACHE[key] = (d, w)
    return _NODE_CACHE[key]


def grid(nu, h, dim):
    d, w = exp_nodes(nu, h)
    Ds = np.meshgrid(*([d] * dim), indexing="ij")
    Ws = np.meshgrid(*([w] * dim), indexing="ij")
    wt = np.ones_like(Ds[0])
    for x in Ws:
        wt = wt * x
    return [x.ravel() for x in Ds], wt.ravel()


# ------------------------------------------------- X2: inflation-robust obstruction
def Psi(tau, nu, a, beta):
    return E1(lambda d: (a + beta * (tau * d) ** 2) / (1 + beta * tau * d) ** 2, nu,
              1 / (beta * tau))


def dPsi(tau, nu, a, beta):
    return 2 * beta * E1(lambda d: d * (tau * d - a) / (1 + beta * tau * d) ** 3, nu,
                         1 / (beta * tau))


def gamma_infl(tau, nu, a, beta):
    return E1(lambda d: (tau * d - 1) * (tau * d - a) / (1 + beta * tau * d) ** 2, nu,
              1 / (beta * tau))


def gamma_lead_identity(tau, nu, a, beta):
    """E[psi (rho - rho(a))], psi = X(X-a)g^3, rho = (X-1)(1+beta X)/X."""
    ra = (a - 1) * (1 + beta * a) / a

    def f(d):
        X = tau * d
        g = 1 / (1 + beta * X)
        return X * (X - a) * g ** 3 * ((X - 1) * (1 + beta * X) / X - ra)
    return E1(f, nu, 1 / (beta * tau))


def tau_star(nu, a, beta):
    lo, hi = 1e-6, 1e6
    flo, fhi = dPsi(lo, nu, a, beta), dPsi(hi, nu, a, beta)
    if not (flo < 0 < fhi):
        return None
    return optimize.brentq(lambda t: dPsi(t, nu, a, beta), lo, hi, xtol=1e-13, rtol=1e-13)


def rot_risk_scalar(theta, nu, a, beta, W, tau=1.0, n=0.25):
    """Exact expected risk (coupling representation) for the 2-D scalar-observation
    model: rotate by theta, taper diag(tau,1) in rotated coords, update."""
    (D, Ev), wt = grid(nu, n, 2)
    c, s = math.cos(theta), math.sin(theta)
    l1, l2 = tau * D, Ev
    m11 = c * c * l1 + s * s * l2
    m21 = c * s * (l1 - l2)
    h = math.sqrt(beta)
    k1 = h * m11 / (1 + beta * m11)
    k2 = h * m21 / (1 + beta * m11)
    r = (a * (W[0, 0] + W[1, 1]) - 2 * a * h * (W[0, 0] * k1 + W[0, 1] * k2)
         + (a * beta + 1) * (W[0, 0] * k1 * k1 + 2 * W[0, 1] * k1 * k2 + W[1, 1] * k2 * k2))
    return float((wt * r).sum())


def x2_checks(res):
    out = {"grid": []}
    W = np.array([[2.0, 1.0], [1.0, 1.0]])
    nus = [1, 2, 3, 4, 9, 19, 49]
    betas = [0.05, 0.3, 1.0, 1.9, 2.5, 5.0, 20.0]
    nmiss = 0
    min_gamma = math.inf
    max_idgap = 0.0
    max_fd = 0.0
    for nu in nus:
        N = nu + 1
        for a in (1.0, 1.0 + 1.0 / N):
            for beta in betas:
                ts = tau_star(nu, a, beta)
                rec = {"N": N, "a": a, "beta": beta, "tau_star": ts}
                if ts is None:
                    rec["note"] = "no sign change of Psi' on [1e-6,1e6]"
                    nmiss += 1
                    # predicted only for N=2 with a*beta>3
                    check(f"X2 no-minimiser only when N=2,a*beta>3 (N={N},a={a:.3f},b={beta})",
                          nu == 1 and a * beta > 3, a_beta=a * beta)
                    out["grid"].append(rec)
                    continue
                if nu == 1:
                    check(f"X2 minimiser exists, N=2 a*beta<3 consistent (a={a},b={beta})",
                          a * beta < 3 + 1e-12, a_beta=a * beta)
                g = gamma_infl(ts, nu, a, beta)
                gid = gamma_lead_identity(ts, nu, a, beta)
                min_gamma = min(min_gamma, g)
                max_idgap = max(max_idgap, abs(g - gid) / max(1e-300, abs(g)))
                rec.update(gamma=g, gamma_identity=gid)
                # global minimum check on a log grid (finite minimiser beats grid)
                tg = np.logspace(-3, 4, 57)
                pv = [Psi(t, nu, a, beta) for t in tg]
                rec["Psi_star"] = Psi(ts, nu, a, beta)
                rec["Psi_grid_min"] = min(pv)
                rec["one_over_beta"] = 1 / beta
                check(f"X2 tau* is a global minimiser on grid (N={N},a={a:.3f},b={beta})",
                      rec["Psi_star"] <= min(pv) + 1e-10 and rec["Psi_star"] < 1 / beta,
                      **{k: rec[k] for k in ("Psi_star", "Psi_grid_min", "one_over_beta")})
                # derivative of rotated competitor vs central finite differences
                if nu <= 19 and beta <= 5:
                    hst = 1e-4
                    fd = (rot_risk_scalar(hst, nu, a, beta, W, ts)
                          - rot_risk_scalar(-hst, nu, a, beta, W, ts)) / (2 * hst)
                    an = 2 * beta * W[0, 1] * g
                    rec.update(fd=fd, analytic=an)
                    err = abs(fd - an)
                    max_fd = max(max_fd, err / max(1.0, abs(an)))
                    check(f"X2 derivative FD (N={N},a={a:.3f},b={beta})",
                          err <= 1e-6 * max(1.0, abs(an)), fd=fd, analytic=an)
                out["grid"].append(rec)
    check("X2 gamma(tau*) > 0 on whole grid", min_gamma > 0, min_gamma=min_gamma)
    check("X2 lead identity gamma = E[psi(rho-rho(a))] (rel 1e-7)", max_idgap < 1e-7,
          max_rel_gap=max_idgap)
    out["min_gamma"] = min_gamma
    out["max_rel_fd_err"] = max_fd
    # derivative also positive for tau away from tau*? (not claimed; recorded)
    # N = 2 large-tau asymptotics: sqrt(2 pi beta tau) J(tau) -> (a beta - 3) pi / 2
    asym = []
    for a in (1.0, 1.5):
        for beta in (0.5, 1.5, 2.5, 3.5, 6.0):
            c = a * beta
            phi = lambda y: (c - 1 - 2 * y) / (1 + y) ** 2
            vals = []
            for bt in (1e4, 1e6, 1e8):
                # sqrt(2 pi beta tau) J = 2 int phi(s^2) exp(-s^2/(2 beta tau)) ds
                L = math.sqrt(2 * bt) * 12
                pts = [0, 1, 10, 100, 1e3, 1e4, L]
                v = 0.0
                for lo, hi in zip(pts[:-1], pts[1:]):
                    if lo >= L:
                        break
                    v += integrate.quad(lambda s: 2 * phi(s * s) * math.exp(-s * s / (2 * bt)),
                                        lo, min(hi, L), limit=400, epsabs=1e-13)[0]
                vals.append(v)
            pred = (c - 3) * math.pi / 2
            asym.append({"a": a, "beta": beta, "scaled_J": vals, "pred": pred})
            check(f"X2 N=2 asymptotic constant (a*beta={c})", abs(vals[-1] - pred) < 2e-3,
                  scaled=vals[-1], pred=pred)
    out["N2_asymptotics"] = asym
    # N=2, a*beta>3: Psi(tau) > 1/beta on a wide grid (variation-diminishing claim)
    nd = []
    for a, beta in ((1.0, 3.5), (1.5, 2.5), (1.5, 6.0), (1.0, 10.0)):
        tg = np.logspace(-3, 7, 81)
        gap = min(Psi(t, 1, a, beta) - 1 / beta for t in tg)
        nd.append({"a": a, "beta": beta, "min_Psi_minus_1_over_beta": gap})
        check(f"X2 N=2 a*beta>3: Psi>1/beta on grid (a={a},b={beta})", gap > 0, gap=gap)
    out["N2_no_minimiser"] = nd
    # N=2, a*beta=3: Stein identity J(tau) = 2(1 - E[1/(1+y)])/(beta tau) > 0 (manuscript
    # proof of Lemma (best inflation factor), part 3); here J = E[(c-1-2y)/(1+y)^2], y = beta tau D.
    stein = []
    for beta in (0.7, 2.0, 3.0):
        a = 3.0 / beta
        for tau in (0.01, 0.3, 1.0, 7.0, 100.0):
            s = beta * tau
            J = E1(lambda d: (2.0 - 2.0 * s * d) / (1 + s * d) ** 2, 1, 1 / s)
            m1 = E1(lambda d: 1.0 / (1 + s * d), 1, 1 / s)
            pred = 2 * (1 - m1) / s
            stein.append({"beta": beta, "a": a, "tau": tau, "J": J, "pred": pred})
            check(f"X2 N=2 a*beta=3 Stein identity (beta={beta},tau={tau})",
                  abs(J - pred) <= 1e-9 * max(1.0, abs(pred)) and J > 0, J=J, pred=pred)
    out["N2_boundary_stein"] = stein
    # Monotone-likelihood-ratio claim: Psi' changes sign at most once (from - to +).
    nsc_max = 0
    for nu in (1, 2, 4, 19):
        for a in (1.0, 1.0 + 1.0 / (nu + 1)):
            for beta in (0.05, 1.0, 5.0):
                tg = np.logspace(-3, 5, 49)
                sg = [np.sign(dPsi(t, nu, a, beta)) for t in tg]
                sg = [x for x in sg if x != 0]
                changes = sum(1 for u, v in zip(sg[:-1], sg[1:]) if u != v)
                ok = changes <= 1 and (changes == 0 or sg[0] < 0)
                nsc_max = max(nsc_max, changes)
                check(f"X2 Psi' at most one sign change, - to + (nu={nu},a={a:.3f},b={beta})",
                      ok, changes=changes)
    out["max_sign_changes_dPsi"] = nsc_max
    res["X2"] = out


# ------------------------------------------------- X3: vector observations
def vec_risk(theta, plane, nu, a, G, W, n=None):
    m = G.shape[1]
    if n is None:
        n = {1: 0.3, 2: 0.25}.get(nu, 0.25) if m == 2 else {2: 0.45}.get(nu, 0.35)
    Dl, wt = grid(nu, n, m)
    i, j = plane
    Q = np.eye(m)
    c, s = math.cos(theta), math.sin(theta)
    Q[i, i] = c; Q[i, j] = -s; Q[j, i] = s; Q[j, j] = c
    Lam = np.zeros((len(wt), m, m))
    for q in range(m):
        Lam[:, q, q] = Dl[q]
    M = Q @ Lam @ Q.T
    k = G.shape[0]
    V = np.eye(k) + G @ M @ G.T
    K = np.linalg.solve(V, (M @ G.T).transpose(0, 2, 1)).transpose(0, 2, 1)
    IKG = np.eye(m) - K @ G
    r = (a * np.einsum("ab,nbc,nac->n", W, IKG, IKG)
         + np.einsum("ab,nbc,nac->n", W, K, K))
    return float((wt * r).sum())


def Phi(bi, bj, nu, a, n=0.2):
    (Di, Dj), wt = grid(nu, n, 2)
    qi = 1 + bi * Di
    qj = 1 + bj * Dj
    f = (Di - Dj) * (bi * (Di - a) / (qi ** 2 * qj) + bj * (Dj - a) / (qi * qj ** 2))
    return float((wt * f).sum())


def x3_checks(res):
    out = {"fd": []}
    b = np.array([2.0, 0.5, 0.0])
    R0 = np.array([[0.6, -0.8], [0.8, 0.6]])
    G = R0 @ np.hstack([np.diag(np.sqrt(b[:2])), np.zeros((2, 1))])
    check("X3 G^T G diagonal", np.allclose(G.T @ G, np.diag(b), atol=1e-14))
    W = np.array([[2.0, 0.6, 0.4], [0.6, 1.5, 0.3], [0.4, 0.3, 1.0]])
    R1 = np.array([[0.28, 0.96], [-0.96, 0.28]])
    G2 = R1 @ np.diag(np.sqrt([3.0, 0.4]))          # m=2, two observed directions
    W2 = np.array([[1.0, -0.7], [-0.7, 2.0]])
    check("X3 G2^T G2 diagonal", np.allclose(G2.T @ G2, np.diag([3.0, 0.4]), atol=1e-14))
    cases = [(nu, G2, W2, np.array([3.0, 0.4]), [(0, 1)]) for nu in (1, 2, 5)]
    cases += [(nu, G, W, b, [(0, 1), (0, 2), (1, 2)]) for nu in (2, 5)]
    for nu, GG, WW, bb, planes in cases:
        N = nu + 1
        for a in (1.0, 1 + 1 / N):
            for plane in planes:
                hst = 1e-4
                fd = (vec_risk(hst, plane, nu, a, GG, WW) - vec_risk(-hst, plane, nu, a, GG, WW)) / (2 * hst)
                i, j = plane
                an = 2 * WW[i, j] * Phi(bb[i], bb[j], nu, a)
                out["fd"].append({"N": N, "a": a, "m": GG.shape[1], "plane": plane, "fd": fd, "analytic": an})
                check(f"X3 first variation FD (m={GG.shape[1]},N={N},a={a:.3f},plane={plane})",
                      abs(fd - an) <= 2e-6 * max(1.0, abs(an)), fd=fd, analytic=an)
    # observed-unobserved reduces to beta*gamma (Theorem 2 constant)
    for nu in (1, 2, 9):
        for a in (1.0, 1 + 1 / (nu + 1)):
            for bi in (0.1, 1.0, 7.0):
                p = Phi(bi, 0.0, nu, a)
                g = E1(lambda d: (d - 1) * (d - a) / (1 + bi * d) ** 2, nu, 1 / bi)
                check(f"X3 Phi(b,0)=b*gamma (nu={nu},a={a:.3f},b={bi})",
                      abs(p - bi * g) < 1e-8 * max(1, abs(p)) and p > 0, Phi=p, bgamma=bi * g)
    # antisymmetry and tie
    check("X3 Phi antisymmetric", abs(Phi(2.0, 0.3, 2, 1.2) + Phi(0.3, 2.0, 2, 1.2)) < 1e-12)
    check("X3 Phi(b,b)=0", abs(Phi(1.3, 1.3, 3, 1.25)) < 1e-12)
    # numerical sign scan for two observed directions (NOT a proof)
    worst = None
    cnt = 0
    for nu in (1, 2, 3, 5, 10, 30):
        for a in (1.0, 1 + 1 / (nu + 1)):
            bs = [0.003, 0.03, 0.1, 0.3, 0.7, 1.0, 1.5, 3.0, 10.0, 30.0, 100.0]
            for x in bs:
                for y in bs:
                    if x <= y:
                        continue
                    p = Phi(x, y, nu, a, n=0.3)
                    cnt += 1
                    sc = p / (x - y)
                    if worst is None or sc < worst[0]:
                        worst = (sc, nu, a, x, y, p)
    out["sign_scan"] = {"count": cnt, "min_Phi_over_(bi-bj)": worst[0], "at": worst[1:]}
    check("X3 numerical: Phi(bi,bj)>0 whenever bi>bj>0 on scan", worst[0] > 0, worst=worst)
    res["X3"] = out


# ------------------------------------------------- X4 and X5: curvature constants
def AB(nu, a, beta):
    A = E1(lambda d: (d - a) * (d - 1) / (1 + beta * d) ** 3, nu, 1 / beta)
    B = E1(lambda d: ((d - 1) ** 2 + 2 / nu) / (1 + beta * d) ** 2, nu, 1 / beta)
    return A, B


def x4_x5_checks(res):
    out = {}
    worst = None
    rows = []
    Ns = [2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 100, 300, 1000]
    betas = np.logspace(-3, 3, 31)
    for N in Ns:
        nu = N - 1
        for a in (1.0, 1 + 1 / N):
            for beta in betas:
                A, B = AB(nu, a, beta)
                ratio = (a * beta + 1) * B / (2 * A)
                if A <= 0:
                    check(f"X4 A>0 (N={N},a={a},b={beta})", False, A=A)
                if worst is None or ratio < worst[0]:
                    worst = (ratio, N, a, float(beta))
            rows.append({"N": N, "a": a, "ratio_beta_1e-3": None})
    out["W_I_min_ratio"] = {"ratio": worst[0], "N": worst[1], "a": worst[2], "beta": worst[3],
                            "grid": {"N": Ns, "beta": [float(x) for x in betas]}}
    check("X4 numerical: 2A<(a beta+1)B on grid (W=I strict local min)", worst[0] > 1, worst=worst)
    # large-N limit ratio -> (1+beta)^2
    for beta in (0.5, 2.0):
        A, B = AB(999, 1 + 1 / 1000, beta)
        r = (1 + 1 / 1000 * beta + beta) * B / (2 * A)
        check(f"X4 large-N ratio ~ (1+beta)^2 (beta={beta})", abs(r / (1 + beta) ** 2 - 1) < 0.02,
              ratio=r, limit=(1 + beta) ** 2)
    # negative control: aligned W=diag(4,1) at N=3 is a local MAX: threshold must be < 4
    A, B = AB(2, 4 / 3, 1.0)
    thr = (4 / 3 + 1) * B / (2 * A)
    check("NEG X4 control: aligned example threshold < w1/w2=4 detected", thr < 4, threshold=thr)
    # X5: aligned example R''(0) = -16 A + 2(a+1) B (beta=1, w1=4, w2=1)
    tab = []
    for N in [2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 100, 200, 500, 1000]:
        nu = N - 1
        rec = {"N": N}
        for lab, a in (("known_mean", 1.0), ("ensemble_mean", 1 + 1 / N)):
            A, B = AB(nu, a, 1.0)
            rec[lab] = -16 * A + 2 * (a + 1) * B
            rec[lab + "_xN2"] = rec[lab] * N * N
        tab.append(rec)
    out["aligned_curvature"] = tab
    c = math.e * exp1(1.0)
    for lab, a in (("known_mean", 1.0), ("ensemble_mean", 4 / 3)):
        ex = 28 * a + 60 - (50 * a + 98) * c
        got = [r for r in tab if r["N"] == 3][0][lab]
        check(f"X5 N=3 closed form ({lab})", abs(ex - got) < 1e-10, exact=ex, quad=got)
    check("X5 N=3 bounds -56/773, -508/773 hold",
          [r for r in tab if r["N"] == 3][0]["known_mean"] <= -56 / 773
          and [r for r in tab if r["N"] == 3][0]["ensemble_mean"] <= -508 / 773)
    check("X5 ensemble mean: R''(0)<0 for all tabulated N",
          all(r["ensemble_mean"] < 0 for r in tab))
    check("X5 known mean: R''(0)>0 at N=2 and <0 for tabulated N>=3",
          tab[0]["known_mean"] > 0 and all(r["known_mean"] < 0 for r in tab[1:]))
    # negative control: claiming the N=3 curvature is <= -2 must be rejected
    check("NEG X5 wrong constant R3''<=-2 rejected",
          not ([r for r in tab if r["N"] == 3][0]["ensemble_mean"] <= -2.0))
    res["X4_X5"] = out


# ------------------------------------------------- Monte Carlo of the real procedure
def mc_real(L, N, thetas, taus, n_total, seed, known_mean=False, chunk=200000):
    """Physical-coordinate EnKF analysis, scalar observation H=(1,0), R=1.
    For each (theta,tau): transform w = Q^T L^{-1} x, sample covariance of the actual
    ensemble in w, keep diag, multiply the first variance by tau, update the ensemble
    mean (or known mean) with obs row H L Q, transform back, physical squared error.
    Returns per-sample losses for every (theta,tau) under common random numbers."""
    rng = np.random.default_rng(seed)
    m = 2
    Linv = np.linalg.inv(L)
    H = np.array([1.0, 0.0])
    out = {key: [] for key in zip(thetas, taus)}
    done = 0
    while done < n_total:
        n = min(chunk, n_total - done)
        Zens = rng.standard_normal((n, N, m))
        Xens = Zens @ L.T
        xt = rng.standard_normal((n, m)) @ L.T
        y = xt @ H + rng.standard_normal(n)
        xbar = Xens.mean(1)
        Xc = Xens - xbar[:, None, :]
        Sx = np.einsum("nri,nrj->nij", Xc, Xc) / (N - 1)
        prior = np.zeros((n, m)) if known_mean else xbar
        for th, tau in zip(thetas, taus):
            c, s = math.cos(th), math.sin(th)
            Q = np.array([[c, -s], [s, c]])
            T = Q.T @ Linv                    # w = T x
            Tinv = np.linalg.inv(T)
            Sw = np.einsum("ij,njk,lk->nil", T, Sx, T)
            dw = np.stack([tau * Sw[:, 0, 0], Sw[:, 1, 1]], 1)
            Hw = H @ Tinv                     # observation row in w coordinates
            wprior = prior @ T.T
            PHt = dw * Hw[None, :]
            denom = 1.0 + (Hw[None, :] ** 2 * dw).sum(1)
            Kw = PHt / denom[:, None]
            innov = y - wprior @ Hw
            wa = wprior + Kw * innov[:, None]
            xa = wa @ Tinv.T
            e = xa - xt
            out[(th, tau)].append((e * e).sum(1))
        done += n
    return {k: np.concatenate(v) for k, v in out.items()}


def mc_canon(W, h, N, configs, n_total, seed, chunk=200000):
    """Real ensemble procedure in canonical coordinates (scalar observation y = h z_1 + e).
    configs: list of (theta, t12). For theta != 0: complete localisation in the basis Q_theta
    (rotate the sample covariance S of the actual ensemble, keep the sample variances, map
    back); t12 is ignored. For theta == 0: canonical taper [[1,t12],[t12,1]] applied to S
    (t12 = 0 is complete canonical localisation). Then update the ensemble
    mean, and score (z_hat - z)^T W (z_hat - z), which equals physical squared error.
    Returns per-sample losses for every config under common random numbers."""
    rng = np.random.default_rng(seed)
    beta = h * h
    out = {c: [] for c in configs}
    done = 0
    while done < n_total:
        n = min(chunk, n_total - done)
        Z = rng.standard_normal((n, N, 2))
        zt = rng.standard_normal((n, 2))
        y = h * zt[:, 0] + rng.standard_normal(n)
        zbar = Z.mean(1)
        Zc = Z - zbar[:, None, :]
        S = np.einsum("nri,nrj->nij", Zc, Zc) / (N - 1)
        for th, t12 in configs:
            c, s_ = math.cos(th), math.sin(th)
            Q = np.array([[c, -s_], [s_, c]])
            Sw = np.einsum("ji,njk,kl->nil", Q, S, Q)
            Dw = np.zeros_like(Sw)
            Dw[:, 0, 0] = Sw[:, 0, 0]
            Dw[:, 1, 1] = Sw[:, 1, 1]
            M = np.einsum("ij,njk,lk->nil", Q, Dw, Q)
            if th == 0.0:
                M = S.copy()
                M[:, 0, 1] *= t12
                M[:, 1, 0] *= t12
            k = h * M[:, :, 0] / (1 + beta * M[:, 0, 0])[:, None]
            za = zbar + k * (y - h * zbar[:, 0])[:, None]
            e = za - zt
            out[(th, t12)].append(np.einsum("ni,ij,nj->n", e, W, e))
        done += n
    return {k: np.concatenate(v) for k, v in out.items()}


def mstat(x):
    return float(x.mean()), float(x.std(ddof=1) / math.sqrt(len(x)))


def mc_checks(res, n_mc):
    out = {}
    L = np.array([[1.0, 0.0], [1.0, 1.0]])
    W = L.T @ L
    N = 3
    nu = 2
    a = 1 + 1 / N
    beta = 1.0
    dlt = 0.04
    # (1) Theorem 2 derivative, tau=1 ; (2) inflation-robust derivative at tau*
    ts = tau_star(nu, a, beta)
    thetas = [-dlt, 0.0, dlt, -dlt, dlt, 0.15, -0.15]
    taus = [1.0, 1.0, 1.0, ts, ts, 1.0, 1.0]
    t0 = time.time()
    R = mc_real(L, N, thetas, taus, n_mc, seed=20261103)
    out["mc_seconds"] = time.time() - t0
    out["n_mc"] = n_mc
    g1 = E1(lambda d: (d - 1) * (d - a) / (1 + d) ** 2, nu)
    an1 = 2 * beta * W[0, 1] * g1
    fd1 = (R[(dlt, 1.0)] - R[(-dlt, 1.0)]) / (2 * dlt)
    m1, se1 = mstat(fd1)
    # O(dlt^2) bias of central difference, from the exact quadrature risk
    ex_fd1 = (rot_risk_scalar(dlt, nu, a, beta, W) - rot_risk_scalar(-dlt, nu, a, beta, W)) / (2 * dlt)
    z1 = (m1 - ex_fd1) / se1
    out["thm2"] = {"analytic": an1, "exact_central_difference": ex_fd1, "mc": m1, "se": se1, "z": z1}
    check("MC Thm2 derivative within 4 SE (real ensembles)", abs(z1) < 4, **out["thm2"])
    small = (rot_risk_scalar(1e-4, nu, a, beta, W) - rot_risk_scalar(-1e-4, nu, a, beta, W)) / 2e-4
    out["thm2"]["central_difference_bias_at_delta"] = ex_fd1 - an1
    check("MC-setting Thm2: exact small-step FD equals analytic (1e-6)",
          abs(small - an1) < 1e-6, fd=small, an=an1)
    g2 = gamma_infl(ts, nu, a, beta)
    an2 = 2 * beta * W[0, 1] * g2
    fd2 = (R[(dlt, ts)] - R[(-dlt, ts)]) / (2 * dlt)
    m2, se2 = mstat(fd2)
    ex_fd2 = (rot_risk_scalar(dlt, nu, a, beta, W, ts) - rot_risk_scalar(-dlt, nu, a, beta, W, ts)) / (2 * dlt)
    z2 = (m2 - ex_fd2) / se2
    out["inflation"] = {"tau_star": ts, "analytic": an2, "exact_central_difference": ex_fd2,
                        "mc": m2, "se": se2, "z": z2}
    check("MC inflation-robust derivative within 4 SE (real ensembles)", abs(z2) < 4, **out["inflation"])
    check("MC inflation-robust derivative positive (>5 SE)", m2 / se2 > 5, z_pos=m2 / se2)
    # canonical risk vs Table 1 (N=3: 2.831656)
    mc0, se0 = mstat(R[(0.0, 1.0)])
    out["canonical_risk"] = {"mc": mc0, "se": se0, "table": 2.831656}
    check("MC canonical risk vs Table 1 within 4 SE", abs(mc0 - 2.831656) < 4 * se0, mc=mc0, se=se0)
    # negative control: wrong rotation sign must increase risk; right sign decreases
    up, seu = mstat(R[(0.15, 1.0)] - R[(0.0, 1.0)])
    dn, sed = mstat(R[(-0.15, 1.0)] - R[(0.0, 1.0)])
    out["sign_control"] = {"R(+0.15)-R(0)": [up, seu], "R(-0.15)-R(0)": [dn, sed]}
    check("NEG MC wrong rotation sign increases risk (>4 SE)", up / seu > 4, z=up / seu)
    check("MC correct rotation sign decreases risk (>4 SE)", dn / sed < -4, z=dn / sed)
    # negative control: wrong constant (2x analytic) must be rejected
    zbad = (m1 - 2 * an1) / se1
    check("NEG MC deliberately wrong constant (2x gamma) rejected (|z|>8)", abs(zbad) > 8, z=zbad)
    # negative control: W=I (P=I) first variation is zero
    R2 = mc_real(np.eye(2), N, [-dlt, dlt], [1.0, 1.0], n_mc // 2, seed=20261104)
    m3, se3 = mstat((R2[(dlt, 1.0)] - R2[(-dlt, 1.0)]) / (2 * dlt))
    out["W_I_derivative"] = {"mc": m3, "se": se3}
    check("NEG MC W=I first variation zero (|z|<4)", abs(m3) < 4 * se3, mc=m3, se=se3)
    # Lemma 1 (canonical taper with T12 = 0.5), Prop 4 curvature (aligned example) and a
    # second Table 1 entry (rotated risk at theta_num, N = 3), all with real ensembles.
    n2 = n_mc // 2
    th_num = -2 * math.atan(1 / 8)
    Rc = mc_canon(W, 1.0, N, [(0.0, 0.0), (0.0, 0.5), (th_num, 0.0)], n2, seed=20261105)
    eta = float(mpmath.e * mpmath.e1(1))
    pred_l1 = beta * (a * beta + 1) / nu * (2 * eta - 1) * W[1, 1] * 0.25
    ml1, sel1 = mstat(Rc[(0.0, 0.5)] - Rc[(0.0, 0.0)])
    out["lemma1_taper_excess"] = {"mc": ml1, "se": sel1, "predicted": pred_l1}
    check("MC Lemma 1 taper excess within 4 SE (real ensembles)", abs(ml1 - pred_l1) < 4 * sel1,
          mc=ml1, se=sel1, predicted=pred_l1)
    mr, ser = mstat(Rc[(th_num, 0.0)])
    out["table1_rotated_N3"] = {"mc": mr, "se": ser, "table": 2.746111}
    check("MC rotated risk vs Table 1 (N=3) within 4 SE", abs(mr - 2.746111) < 4 * ser, mc=mr, se=ser)
    Wal = np.diag([4.0, 1.0])
    dc = 0.2
    Ra = mc_canon(Wal, 1.0, N, [(-dc, 0.0), (0.0, 0.0), (dc, 0.0)], n2, seed=20261106)
    m4, se4 = mstat((Ra[(dc, 0.0)] - 2 * Ra[(0.0, 0.0)] + Ra[(-dc, 0.0)]) / dc ** 2)
    ex4 = (rot_risk_scalar(dc, nu, a, beta, Wal) - 2 * rot_risk_scalar(0.0, nu, a, beta, Wal)
           + rot_risk_scalar(-dc, nu, a, beta, Wal)) / dc ** 2
    an4 = 28 * a + 60 - (50 * a + 98) * eta
    out["prop4_curvature"] = {"mc_second_difference": m4, "se": se4,
                              "exact_second_difference": ex4, "analytic": an4, "step": dc}
    check("MC Prop 4 curvature: second difference within 4 SE of exact (real ensembles)",
          abs(m4 - ex4) < 4 * se4, mc=m4, se=se4, exact=ex4)
    check("MC Prop 4 curvature negative (>4 SE)", m4 / se4 < -4, z=m4 / se4)
    small4 = (rot_risk_scalar(1e-3, nu, a, beta, Wal) - 2 * rot_risk_scalar(0.0, nu, a, beta, Wal)
              + rot_risk_scalar(-1e-3, nu, a, beta, Wal)) / 1e-6
    check("Prop 4: exact small-step second difference equals analytic (1e-4)",
          abs(small4 - an4) < 1e-4, fd=small4, analytic=an4)
    res["MC"] = out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-mc", type=int, default=4_000_000)
    args = ap.parse_args()
    outdir = os.path.abspath(args.out)
    if os.path.exists(outdir):
        print(f"refusing to overwrite existing path {outdir}", file=sys.stderr)
        raise SystemExit(2)
    os.makedirs(outdir)
    t0 = time.time()
    with open(os.path.abspath(__file__), "rb") as fh:
        sha = hashlib.sha256(fh.read()).hexdigest()
    res = {"script_sha256": sha, "python": sys.version, "numpy": np.__version__,
           "scipy": scipy.__version__, "mpmath": mpmath.__version__,
           "platform": platform.platform(), "optimize_flag": sys.flags.optimize}
    timings = {}
    import warnings
    wlog = []
    warnings.showwarning = lambda msg, cat, fn_, ln, file=None, line=None: wlog.append(f"{cat.__name__}: {msg}")
    for name, fn in (("X2", x2_checks), ("X3", x3_checks), ("X4_X5", x4_x5_checks)):
        t = time.time()
        fn(res)
        timings[name] = time.time() - t
        print(f"[{name}] done in {timings[name]:.1f}s; failures so far {len(FAILURES)}", flush=True)
    t = time.time()
    mc_checks(res, args.n_mc)
    timings["MC"] = time.time() - t
    print(f"[MC] done in {timings['MC']:.1f}s; failures so far {len(FAILURES)}", flush=True)
    res["timings_s"] = timings
    res["warnings"] = sorted(set(wlog))
    res["n_warnings"] = len(wlog)
    res["total_seconds"] = time.time() - t0
    res["checks"] = CHECKS
    res["n_checks"] = len(CHECKS)
    res["failures"] = FAILURES
    res["status"] = "PASS" if not FAILURES else "FAIL"
    with open(os.path.join(outdir, "theory_checks.json"), "w") as fh:
        json.dump(_js(res), fh, indent=1, sort_keys=True)
    print(f"status={res['status']} checks={len(CHECKS)} failures={len(FAILURES)} "
          f"total={res['total_seconds']:.1f}s")
    for f in FAILURES:
        print("FAILED:", f)
    if FAILURES:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
