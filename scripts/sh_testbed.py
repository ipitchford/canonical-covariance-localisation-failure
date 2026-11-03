#!/usr/bin/env python3
"""Reconstruction of the Snyder & Hakim (2022, JAMES, doi:10.1029/2021MS002937)
Section 4 one-dimensional testbed, evaluated with the leading-order basis
objective F(U) of the manuscript (Theorem 5) and an exact-conditional-risk
Monte Carlo of complete localisation in several whitening bases.

This is a RECONSTRUCTION of SH's setting, not a replication of their figures.

Setting (see NOTES.md for the justification of each choice):
  * Nx = 200 equally spaced points x_i = (i + 1/2)/Nx on (0,1), not periodic;
  * Matern nu = 5/2 correlation, C(r) = (1 + s + s^2/3) exp(-s), s = sqrt(5) r / rho,
    rho calibrated so that C(0.038) = 1/2 exactly (SH: "drops to 1/2 at r = 0.038");
  * unit prior variance; Ny = 50 observations of single grid values, every 4th point
    (indices 2, 6, ..., 198), independent unit-variance observation errors (R = I);
  * physical squared loss: W = L^T L in canonical coordinates, x = L z;
  * MSE normalised by the Kalman-filter (Bayes) MSE tr(WC), as in SH Figs 3 and 5.

Monte Carlo. For a fixed orthonormal basis U of the whitened state z ~ N(0, I),
complete localisation keeps the sample variances of w = U^T z. By rotation
invariance these are D_i iid chi2_nu / nu (nu = N - 1) and independent of the
ensemble mean, so the analysis-error covariance given D is exact:
  E(K) = a (I - K G)(I - K G)^T + K K^T,  K = M G^T (I + G M G^T)^{-1},
  M = U diag(D) U^T,  a = 1 + 1/N (ensemble mean) or a = 1 (known mean).
The same D draws are used for every basis (common random numbers).
The estimated-basis test does NOT use this shortcut: it simulates the forecast
perturbations in physical space, because an estimated whitening does not make
the transformed variables independent.

Usage:  python3 -I -B sh_testbed.py --out NEW_DIR [--reps R] [--quick]
No assert statements: failures are collected, a status is computed and the
process exits with status 1 when any check fails.
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import hashlib
import json
import platform
import sys
import time

import numpy as np
import scipy
import scipy.linalg as sl
from scipy.optimize import brentq

FAILURES = []
CHECKS = []


def check(name, ok, detail):
    """Record a named check; failures are collected, never raised."""
    CHECKS.append({"name": name, "ok": bool(ok), "detail": detail})
    if not ok:
        FAILURES.append(name)


# --------------------------------------------------------------------------
# Model construction
# --------------------------------------------------------------------------
def matern52(s):
    return (1.0 + s + s * s / 3.0) * np.exp(-s)


def calibrated_rho(r_half=0.038):
    x_half = brentq(lambda s: matern52(s) - 0.5, 0.5, 5.0, xtol=1e-15)
    return np.sqrt(5.0) * r_half / x_half, x_half


def build_model(nx=200, ny=50, rho=None, loss="physical"):
    if rho is None:
        rho, _ = calibrated_rho()
    xs = (np.arange(nx) + 0.5) / nx
    r = np.abs(xs[:, None] - xs[None, :])
    P = matern52(np.sqrt(5.0) * r / rho)
    ev, E = np.linalg.eigh(P)
    ev_min_raw = float(ev.min())
    ev = np.clip(ev, 1e-14, None)
    Phalf = (E * np.sqrt(ev)) @ E.T
    stride = nx // ny
    obs = stride * np.arange(ny) + stride // 2
    H = np.zeros((ny, nx))
    H[np.arange(ny), obs] = 1.0
    G = H @ Phalf                                    # R = I
    Uo, lam, Vt = np.linalg.svd(G, full_matrices=True)
    V = Vt.T
    # canonical observation operator: Uo^T G V = [diag(lam) 0]
    Gc = Uo.T @ G @ V
    off = Gc.copy()
    off[np.arange(ny), np.arange(ny)] = 0.0
    gc_offdiag = float(np.abs(off).max())
    if loss == "physical":
        W = V.T @ P @ V                              # W = L^T L with L = P^{1/2} V
    else:
        W = np.eye(nx)                               # prior-normalised loss
    W = 0.5 * (W + W.T)
    b = np.zeros(nx)
    b[:ny] = lam ** 2 / (1.0 + lam ** 2)
    c = 1.0 - b
    A = (c[:, None] * W) * c[None, :]
    A = 0.5 * (A + A.T)
    # alternative bases expressed as orthonormal bases of the canonical z-space
    U_eig = V.T @ E[:, ::-1]                         # SH Eq. 15 (eigenvectors of P)
    U_zca = V.T.copy()                               # symmetric (ZCA) whitening P^{-1/2}
    return dict(nx=nx, ny=ny, rho=rho, xs=xs, P=P, E=E, ev=ev, ev_min_raw=ev_min_raw,
                Phalf=Phalf, obs=obs, H=H, V=V, lam=lam, W=W, b=b, c=c, A=A,
                U_eig=U_eig, U_zca=U_zca, gc_offdiag=gc_offdiag)


# --------------------------------------------------------------------------
# Leading-order objective and its Riemannian (Cayley) descent
# --------------------------------------------------------------------------
def F_val(U, A, b):
    beta = np.einsum("ij,ij->j", U, b[:, None] * U)
    alpha = np.einsum("ij,ij->j", U, A @ U)
    return float(beta @ alpha)


def F_grad(U, A, b):
    BU = b[:, None] * U
    AU = A @ U
    beta = np.einsum("ij,ij->j", U, BU)
    alpha = np.einsum("ij,ij->j", U, AU)
    Ge = 2.0 * (BU * alpha + AU * beta)
    X = U.T @ Ge
    return float(beta @ alpha), 0.5 * (X - X.T)


def cayley(U, Om, t):
    m = U.shape[0]
    I = np.eye(m)
    return U @ np.linalg.solve(I + 0.5 * t * Om, I - 0.5 * t * Om)


def descend(U0, A, b, max_iter=4000, gtol=1e-7, label=""):
    """Armijo-backtracking Riemannian gradient descent with Cayley retraction."""
    U = U0.copy()
    f, Om = F_grad(U, A, b)
    f0 = f
    g0 = float(np.linalg.norm(Om))
    t = 1.0 / max(g0, 1e-300)
    hist = [(0, f, g0)]
    it = 0
    stalled = False
    for it in range(1, max_iter + 1):
        g2 = float(np.sum(Om * Om))
        g = np.sqrt(g2)
        if g <= gtol * max(f0, 1.0):
            break
        accepted = False
        for _ in range(60):
            Un = cayley(U, Om, t)
            fn = F_val(Un, A, b)
            if fn <= f - 1e-4 * t * g2:
                accepted = True
                break
            t *= 0.5
        if not accepted:
            stalled = True
            break
        # re-orthonormalise occasionally (Cayley is orthogonal up to round-off)
        if it % 200 == 0:
            Qr, Rr = np.linalg.qr(Un)
            Un = Qr * np.sign(np.diag(Rr))[None, :]
        U = Un
        f, Om = F_grad(U, A, b)
        t *= 2.0
        if it in (10, 30, 100, 300, 1000, 2000, 3000) or it % 1000 == 0:
            hist.append((it, f, float(np.linalg.norm(Om))))
    gfin = float(np.linalg.norm(Om))
    hist.append((it, f, gfin))
    orth = float(np.abs(U.T @ U - np.eye(U.shape[0])).max())
    return U, dict(label=label, F_start=f0, F_final=f, grad_norm_start=g0,
                   grad_norm_final=gfin, iterations=it, stalled=stalled,
                   converged=bool(gfin <= gtol * max(f0, 1.0)),
                   orthogonality_error=orth,
                   history=[{"iter": h[0], "F": h[1], "grad_norm": h[2]} for h in hist])


def haar(m, rng):
    Z = rng.standard_normal((m, m))
    Q, R = np.linalg.qr(Z)
    return Q * np.sign(np.diag(R))[None, :]


def haar_mean_exact(A, b):
    m = A.shape[0]
    trA, trB, trAB = np.trace(A), b.sum(), float(np.sum(np.diag(A) * b))
    # E[(u^T B u)(u^T A u)] = (trA trB + 2 tr(AB)) / (m (m + 2)) for u uniform on the sphere
    return float((trA * trB + 2.0 * trAB) / (m + 2.0)), float(trA * trB / m)


# --------------------------------------------------------------------------
# Exact conditional risk, canonical z-space (fixed oracle basis)
# --------------------------------------------------------------------------
class ZBasis:
    def __init__(self, U, mdl):
        ny = mdl["ny"]
        lam = mdl["lam"]
        self.U = U
        self.GU = lam[:, None] * U[:ny, :]                 # G U, ny x m
        WU = mdl["W"] @ U
        self.GWU = lam[:, None] * WU[:ny, :]                # rows of G W U
        self.Wu = U.T @ WU                                  # U^T W U
        self.trW = float(np.trace(mdl["W"]))
        self.lam2 = lam ** 2


def z_risk(basis, D):
    """Exact conditional risk parts for a batch D (R x m): risk(a) = a*X + Y."""
    GU = basis.GU
    GUD = GU[None, :, :] * D[:, None, :]                    # R x ny x m
    S = GUD @ GU.T                                          # R x ny x ny
    ny = GU.shape[0]
    S[:, np.arange(ny), np.arange(ny)] += 1.0
    Z = np.linalg.solve(S, GUD)                             # K^T = Z U^T
    T1 = Z @ basis.Wu                                       # R x ny x m
    kwk = np.einsum("rij,rij->ri", T1, Z)                   # diag(K^T W K)
    cross = np.einsum("ij,rij->r", basis.GWU, Z)            # tr(G W K)
    X = basis.trW - 2.0 * cross + kwk @ basis.lam2
    Y = kwk.sum(axis=1)
    return X, Y


def z_risk_direct(U, mdl, d, a):
    """Unvectorised reference formula, used only as an implementation check."""
    m, ny = mdl["nx"], mdl["ny"]
    G = np.zeros((ny, m))
    G[np.arange(ny), np.arange(ny)] = mdl["lam"]
    M = (U * d) @ U.T
    K = M @ G.T @ np.linalg.inv(np.eye(ny) + G @ M @ G.T)
    IKG = np.eye(m) - K @ G
    Ecov = a * IKG @ IKG.T + K @ K.T
    return float(np.sum(mdl["W"] * Ecov))


def mc_z(bases, mdl, Ns, reps, seed, chunk=100):
    """CRN Monte Carlo over D for each basis and N; returns per-rep risks."""
    m = mdl["nx"]
    out = {}
    for N in Ns:
        nu = N - 1
        rng = np.random.default_rng([seed, N])
        Xs = {k: [] for k in bases}
        Ys = {k: [] for k in bases}
        done = 0
        while done < reps:
            n = min(chunk, reps - done)
            D = rng.chisquare(nu, size=(n, m)) / nu
            for k, bs in bases.items():
                X, Y = z_risk(bs, D)
                Xs[k].append(X)
                Ys[k].append(Y)
            done += n
        out[N] = {k: (np.concatenate(Xs[k]), np.concatenate(Ys[k])) for k in bases}
    return out


def summarise(mc, kf, Ns, ref="canonical"):
    res = {}
    for N in Ns:
        rN = {}
        for conv, a in (("ensemble_mean", 1.0 + 1.0 / N), ("known_mean", 1.0)):
            vals = {k: (a * X + Y) / kf for k, (X, Y) in mc[N].items()}
            n = len(next(iter(vals.values())))
            ent = {}
            for k, v in vals.items():
                e = {"mean": float(v.mean()), "se": float(v.std(ddof=1) / np.sqrt(n))}
                if k != ref and ref in vals:
                    dlt = v - vals[ref]
                    e["diff_vs_canonical"] = float(dlt.mean())
                    e["se_diff"] = float(dlt.std(ddof=1) / np.sqrt(n))
                    exc_ref = vals[ref].mean() - 1.0
                    e["excess_reduction_vs_canonical"] = float(-dlt.mean() / exc_ref) if exc_ref > 0 else None
                ent[k] = e
            rN[conv] = ent
        res[str(N)] = rN
    return res


def leading(mdl, Fdict, Ns):
    W, c = mdl["W"], mdl["c"]
    trWC = float(np.sum(np.diag(W) * c))
    trWC2 = float(np.sum(np.diag(W) * c * c))
    out = {}
    for N in Ns:
        out[str(N)] = {
            k: {"ensemble_mean": 1.0 + (trWC2 / N + 2.0 * F / (N - 1)) / trWC,
                "known_mean": 1.0 + (2.0 * F / (N - 1)) / trWC}
            for k, F in Fdict.items()}
    return out


# --------------------------------------------------------------------------
# Physical-space simulation for estimated bases
# --------------------------------------------------------------------------
def gaspari_cohn(r, c):
    z = r / c
    out = np.zeros_like(z)
    m1 = z <= 1.0
    m2 = (z > 1.0) & (z < 2.0)
    x = z[m1]
    out[m1] = -0.25 * x ** 5 + 0.5 * x ** 4 + 0.625 * x ** 3 - 5.0 / 3.0 * x ** 2 + 1.0
    x = z[m2]
    out[m2] = (x ** 5 / 12.0 - 0.5 * x ** 4 + 0.625 * x ** 3 + 5.0 / 3.0 * x ** 2
               - 5.0 * x + 4.0 - 2.0 / (3.0 * x))
    return out


def sym_pow(S, p):
    ev, E = np.linalg.eigh(0.5 * (S + S.T))
    ev = np.maximum(ev, 1e-13 * ev.max())
    return (E * ev ** p) @ E.T


def phys_risk_batch(Tz, Tinv, mdl, Zc, a_list):
    """Exact conditional physical risk given forecast perturbations.

    Zc: R x m x nu standard-normal contrasts, forecast perturbations P^{1/2} Zc.
    Tz = T P^{1/2} maps Zc to transformed perturbations; Tinv = T^{-1}.
    """
    P, obs = mdl["P"], mdl["obs"]
    nu = Zc.shape[2]
    Y = Tz[None, :, :] @ Zc
    d = np.sum(Y * Y, axis=2) / nu                          # R x m
    Tobs = Tinv[obs, :]                                     # ny x m
    MHt = (Tinv[None, :, :] * d[:, None, :]) @ Tobs.T       # R x m x ny
    HMHt = MHt[:, obs, :]
    ny = len(obs)
    HMHt[:, np.arange(ny), np.arange(ny)] += 1.0
    K = np.swapaxes(np.linalg.solve(HMHt, np.swapaxes(MHt, 1, 2)), 1, 2)   # R x m x ny
    Pobs = P[np.ix_(obs, obs)]
    PcolObs = P[:, obs]
    trP = float(np.trace(P))
    KtK = np.swapaxes(K, 1, 2) @ K
    X = trP - 2.0 * np.einsum("rij,ij->r", K, PcolObs) + np.einsum("rij,ij->r", KtK, Pobs)
    Yv = np.einsum("rij,rij->r", K, K)
    return {a: a * X + Yv for a in a_list}


def estimated_bases(mdl, Np, rng, gc_half, ridge, descent_iter):
    """Bases built from a GC-localised pilot sample covariance of size Np."""
    m, ny, H = mdl["nx"], mdl["ny"], mdl["H"]
    Zp = rng.standard_normal((m, Np))
    Xp = mdl["Phalf"] @ Zp
    Xp = Xp - Xp.mean(axis=1, keepdims=True)
    S = Xp @ Xp.T / (Np - 1)
    xs = mdl["xs"]
    Tgc = gaspari_cohn(np.abs(xs[:, None] - xs[None, :]), gc_half) if gc_half is not None else 1.0
    Ph = Tgc * S + ridge * np.eye(m)
    Ph_half = sym_pow(Ph, 0.5)
    Ph_mhalf = sym_pow(Ph, -0.5)
    Gh = H @ Ph_half
    _, lh, Vht = np.linalg.svd(Gh, full_matrices=True)
    Vh = Vht.T
    Wh = Vh.T @ Ph @ Vh
    Wh = 0.5 * (Wh + Wh.T)
    bh = np.zeros(m)
    bh[:ny] = lh ** 2 / (1.0 + lh ** 2)
    ch = 1.0 - bh
    Ah = (ch[:, None] * Wh) * ch[None, :]
    Uh, info = descend(np.eye(m), 0.5 * (Ah + Ah.T), bh, max_iter=descent_iter, label=f"pilot Np={Np}")
    rel_err = float(np.linalg.norm(Ph - mdl["P"]) / np.linalg.norm(mdl["P"]))
    P_half = mdl["Phalf"]
    out = {
        "est_canonical": (Vh.T @ Ph_mhalf @ P_half, Ph_half @ Vh),
        "est_zca": (Ph_mhalf @ P_half, Ph_half),
        "est_F_opt": (Uh.T @ Vh.T @ Ph_mhalf @ P_half, Ph_half @ Vh @ Uh),
    }
    # value of the TRUE leading objective for each estimated transform (rows normalised)
    return out, dict(pilot_rel_frobenius_error=rel_err, descent=info)


def rms_offdiag(Tz):
    """RMS off-diagonal TRUE correlation of the transformed variables w = T x.

    Tz = T P^{1/2}, so cov(w) = Tz Tz^T. Zero for an oracle whitening basis.
    """
    Q = Tz @ Tz.T
    dq = np.sqrt(np.diag(Q))
    Cn = Q / dq[:, None] / dq[None, :]
    off = Cn - np.diag(np.diag(Cn))
    return float(np.sqrt(np.sum(off ** 2) / Q.shape[0]))


# --------------------------------------------------------------------------
# Fourier / homogeneous corollary of Theorem 6
# --------------------------------------------------------------------------
def fourier_check(m, rho, r, rng, descent_iter):
    xs = np.arange(m) / m
    dd = np.abs(xs[:, None] - xs[None, :])
    dd = np.minimum(dd, 1.0 - dd)                         # periodic distance
    P = matern52(np.sqrt(5.0) * dd / rho)
    # real orthonormal Fourier basis
    cols = [np.ones(m) / np.sqrt(m)]
    for k in range(1, m // 2):
        cols.append(np.sqrt(2.0 / m) * np.cos(2 * np.pi * k * np.arange(m) / m))
        cols.append(np.sqrt(2.0 / m) * np.sin(2 * np.pi * k * np.arange(m) / m))
    cols.append(np.cos(np.pi * np.arange(m)) / np.sqrt(m))
    Phi = np.array(cols).T
    Dg = Phi.T @ P @ Phi
    p = np.diag(Dg).copy()
    diag_err = float(np.abs(Dg - np.diag(p)).max())
    pos = bool(p.min() > 0)
    # H = I, R = r I: canonical basis = Fourier basis, W = diag(p), c = r/(r+p)
    w = p
    c = r / (r + p)
    a = c * c * w
    b = 1.0 - c
    A = np.diag(a)
    F_can = float(np.sum(a * b))
    tol = 1e-12 * max(1.0, p.max())
    best = None
    best_below = None
    nviol = 0
    nviol_below_r = 0
    npairs_below_r = 0
    for i in range(m):
        for j in range(i + 1, m):
            if abs(p[i] - p[j]) <= tol:
                continue
            prod = (a[i] - a[j]) * (b[i] - b[j])
            below = (p[i] <= r) and (p[j] <= r)
            npairs_below_r += int(below)
            if prod > 0:
                nviol += 1
                nviol_below_r += int(below)
                if best is None or prod > best[0]:
                    best = (prod, i, j)
                if below and (best_below is None or prod > best_below[0]):
                    best_below = (prod, i, j)
    res = dict(m=m, rho=rho, r=r, fourier_diag_error=diag_err, spectrum_positive=pos,
               p_min=float(p.min()), p_max=float(p.max()),
               n_modes_below_r=int(np.sum(p <= r)), F_canonical=F_can,
               distinct_pairs_both_below_r=npairs_below_r,
               violating_pairs=nviol, violating_pairs_both_below_r=nviol_below_r)
    for tag, bp in (("best_pair", best), ("best_pair_both_below_r", best_below)):
        if bp is None:
            continue
        prod, i, j = bp
        U = np.eye(m)
        s = 1.0 / np.sqrt(2.0)
        U[np.ix_([i, j], [i, j])] = np.array([[s, -s], [s, s]])
        F_mix = F_val(U, A, b)
        res[tag] = dict(modes=[int(i), int(j)], p=[float(p[i]), float(p[j])],
                        predicted_45deg_reduction=float(prod / 2.0),
                        numerical_45deg_reduction=float(F_can - F_mix))
    if best is not None:
        # descent from a small random perturbation of the canonical (Fourier) basis
        Xi = rng.standard_normal((m, m))
        Xi = 1e-3 * (Xi - Xi.T)
        U0 = cayley(np.eye(m), Xi, 1.0)
        U_opt, info = descend(U0, A, b, max_iter=descent_iter, label="fourier")
        res["F_after_descent_from_perturbed_canonical"] = info["F_final"]
        res["descent"] = {k: info[k] for k in ("iterations", "converged", "grad_norm_final")}
    # all-pairs statement: every distinct pair with p_j < p_i <= r violates the ordering
    allviol = True
    for i in range(m):
        for j in range(m):
            if p[j] < p[i] - tol and p[i] <= r:
                if not ((a[i] - a[j]) * (b[i] - b[j]) > 0):
                    allviol = False
    res["every_distinct_pair_below_r_violates"] = bool(allviol)
    # worked example from the domain review: r = 1, p = (0.2, 0.5)
    pe = np.array([0.2, 0.5])
    ce = 1.0 / (1.0 + pe)
    ae, be = ce * ce * pe, 1 - ce
    res["worked_example_r1_p02_p05"] = dict(F_canonical=float(np.sum(ae * be)),
                                            reduction_45deg=float((ae[1] - ae[0]) * (be[1] - be[0]) / 2))
    return res


# --------------------------------------------------------------------------
def basis_stats(U, mdl):
    b, A = mdl["b"], mdl["A"]
    beta = np.einsum("ij,ij->j", U, b[:, None] * U)
    alpha = np.einsum("ij,ij->j", U, A @ U)
    eff = float(beta.sum() ** 2 / np.sum(beta ** 2))
    return dict(sum_b=float(beta.sum()), max_b=float(beta.max()),
                n_b_above_0p05=int(np.sum(beta > 0.05)),
                effective_number_of_observed_directions=eff,
                corr_alpha_beta=float(np.corrcoef(alpha, beta)[0, 1]),
                top5_b=[float(v) for v in np.sort(beta)[::-1][:5]])


def config_analysis(mdl, rng, descent_iter, haar_samples):
    m = mdl["nx"]
    A, b, W, c = mdl["A"], mdl["b"], mdl["W"], mdl["c"]
    Cm = np.diag(c)
    comm = float(np.linalg.norm(W @ Cm - Cm @ W) / np.linalg.norm(W @ Cm))
    I = np.eye(m)
    F = {"canonical": F_val(I, A, b), "P_eigen": F_val(mdl["U_eig"], A, b),
         "ZCA": F_val(mdl["U_zca"], A, b)}
    hm_exact, hm_lead = haar_mean_exact(A, b)
    hs = np.array([F_val(haar(m, rng), A, b) for _ in range(haar_samples)])
    hs_mean = float(hs.mean()) if len(hs) else None
    hs_se = float(hs.std(ddof=1) / np.sqrt(len(hs))) if len(hs) > 1 else None
    U_opt, info_can = descend(I, A, b, max_iter=descent_iter, label="from canonical")
    U_opt2, info_zca = descend(mdl["U_zca"], A, b, max_iter=descent_iter, label="from ZCA")
    F["F_opt_from_canonical"] = info_can["F_final"]
    F["F_opt_from_ZCA"] = info_zca["F_final"]
    trWC = float(np.sum(np.diag(W) * c))
    trWC2 = float(np.sum(np.diag(W) * c * c))
    WC = W * c[None, :]
    # relative sd of a single realisation of the squared error at the KF (Gaussian error):
    # sd(e^T W e) / E(e^T W e) = sqrt(2 tr((WC)^2)) / tr(WC)
    single_sd = float(np.sqrt(2.0 * np.sum(WC * WC.T)) / trWC)
    return dict(commutator_rel=comm, single_realisation_relative_sd_at_KF=single_sd,
                implied_se_of_MSE_ratio_800_unpaired_realisations=single_sd / np.sqrt(800.0), F=F, haar_mean_exact=hm_exact,
                haar_mean_trA_trB_over_m=hm_lead,
                haar_sampled_mean=hs_mean, haar_sampled_se=hs_se,
                haar_samples=int(len(hs)), trWC=trWC, trWC2=trWC2, trAB_canonical=float(np.sum(np.diag(A) * b)),
                descent_from_canonical=info_can, descent_from_ZCA=info_zca), U_opt, U_opt2


def eq_ensemble_size(Ngrid, exc_opt, target):
    """Smallest real N with excess_opt(N) = target by log-log interpolation."""
    lN = np.log(np.array(Ngrid, float))
    le = np.log(np.array(exc_opt))
    if target >= exc_opt[0]:
        return float(Ngrid[0]), "at_or_below_grid_minimum"
    if target <= exc_opt[-1]:
        return None, "beyond_grid"
    for k in range(len(Ngrid) - 1):
        if exc_opt[k] >= target >= exc_opt[k + 1]:
            t = (np.log(target) - le[k]) / (le[k + 1] - le[k])
            return float(np.exp(lN[k] + t * (lN[k + 1] - lN[k]))), "interpolated"
    return None, "non_monotone"


def make_figure(outdir, Ns, summ, lead, kf_name):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "svg.hashsalt": "sh-testbed", "pdf.fonttype": 42})
    sty = {"canonical": ("Canonical (SH Eq. 8)", "#b2182b", "o", "-", "full"),
           "P_eigen": ("Eigenvectors of P (SH Eq. 15)", "#f4a582", "s", ":", "none"),
           "ZCA": ("Symmetric (ZCA) whitening", "#2166ac", "^", ":", "none"),
           "F_opt": ("Best basis found for F", "#1b7837", "D", "-", "full")}
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.1), constrained_layout=True)
    Nn = np.array(Ns, float)
    offs = {"canonical": 0.96, "P_eigen": 1.0, "ZCA": 1.0, "F_opt": 1.04}
    for k, (lab, col, mk, ls, fill) in sty.items():
        mu = np.array([summ[str(N)]["ensemble_mean"][k]["mean"] for N in Ns])
        se = np.array([summ[str(N)]["ensemble_mean"][k]["se"] for N in Ns])
        ld = np.array([lead[str(N)][k]["ensemble_mean"] for N in Ns])
        xo = Nn * offs[k]
        kw = dict(color=col, marker=mk, ms=4.5, capsize=2, mfc=(col if fill == "full" else "white"), mew=1.0)
        axs[0].errorbar(xo, mu, yerr=2 * se, lw=1.2, ls=ls, label=lab, **kw)
        axs[1].errorbar(xo, mu - 1.0, yerr=2 * se, lw=0, elinewidth=1.0, **kw)
        if k in ("canonical", "F_opt"):
            axs[1].plot(Nn, ld - 1.0, color=col, lw=1.0, ls="--")
    for ax in axs:
        ax.set_xscale("log")
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.set_xticks(Ns)
        ax.set_xticklabels([str(N) for N in Ns])
        ax.set_xlabel("Ensemble size N")
    axs[0].set_ylabel("Analysis MSE / Kalman-filter MSE")
    axs[0].legend(frameon=False, fontsize=7.5)
    axs[0].set_title("(a) Normalised analysis error", fontsize=9, loc="left")
    axs[1].set_yscale("log")
    axs[1].set_ylabel("MSE / Kalman-filter MSE $-$ 1")
    axs[1].set_title("(b) Excess; dashed: leading-order formula", fontsize=9, loc="left")
    for ext in ("pdf", "png"):
        meta = {"CreationDate": None, "Producer": None} if ext == "pdf" else {"Software": None}
        fig.savefig(os.path.join(outdir, f"sh_testbed_mse.{ext}"), dpi=200, metadata=meta)
    plt.close(fig)


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True, help="new output directory (must not exist)")
    ap.add_argument("--reps", type=int, default=4000, help="Monte Carlo replicates per N (main table)")
    ap.add_argument("--seed", type=int, default=20261103)
    ap.add_argument("--quick", action="store_true", help="small smoke-test settings")
    args = ap.parse_args()

    if os.path.exists(args.out):
        print(f"refusing to overwrite existing path: {args.out}", file=sys.stderr)
        raise SystemExit(1)
    os.makedirs(args.out)
    t_start = time.time()
    timing = {}
    q = args.quick
    reps = 200 if q else args.reps
    dit = 300 if q else 4000
    dit_s = 200 if q else 2500
    seed = args.seed
    rng = np.random.default_rng([seed, 0])

    with open(os.path.abspath(__file__), "rb") as fh:
        script_sha = hashlib.sha256(fh.read()).hexdigest()

    rho0, x_half = calibrated_rho()
    results = {
        "description": "Reconstruction of the Snyder-Hakim (2022) Section 4 testbed; not a replication of their figures.",
        "script_sha256": script_sha,
        "python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__,
        "platform": platform.platform(), "optimize_flag": sys.flags.optimize,
        "settings": dict(reps=reps, seed=seed, quick=q, descent_max_iter=dit, descent_max_iter_sensitivity=dit_s),
        "setup": dict(nx=200, ny=50, domain="(0,1), x_i=(i+1/2)/200, not periodic",
                      matern_nu=2.5, matern_form="(1+s+s^2/3)exp(-s), s=sqrt(5) r/rho",
                      rho=rho0, s_at_half_correlation=x_half, r_half=0.038,
                      r_half_if_rho_sqrt2_over_40=float(x_half * np.sqrt(2) / 40 / np.sqrt(5)),
                      r_half_if_rho_1p5_over_40=float(x_half * 1.5 / 40 / np.sqrt(5)),
                      prior_variance=1.0, obs_error_variance=1.0,
                      obs_indices="4k+2, k=0..49 (point values)", loss="physical sum of squares"),
    }

    # ---------------- base configuration ----------------
    t0 = time.time()
    mdl = build_model(200, 50, rho0, "physical")
    check("prior_covariance_positive_definite", mdl["ev_min_raw"] > 0,
          f"min eigenvalue {mdl['ev_min_raw']:.3e}")
    check("canonical_obs_operator_diagonal", mdl["gc_offdiag"] < 1e-10,
          f"max off-diagonal {mdl['gc_offdiag']:.2e}")
    base, U_opt, U_opt_z = config_analysis(mdl, rng, dit, 0 if q else 40)
    if q:
        hs = np.array([F_val(haar(200, rng), mdl["A"], mdl["b"]) for _ in range(5)])
        base["haar_sampled_mean"] = float(hs.mean())
        base["haar_sampled_se"] = float(hs.std(ddof=1) / np.sqrt(5))
        base["haar_samples"] = 5
    if base["F"]["F_opt_from_ZCA"] < base["F"]["F_opt_from_canonical"]:
        U_best, best_from = U_opt_z, "ZCA"
    else:
        U_best, best_from = U_opt, "canonical"
    base["F"]["F_opt"] = min(base["F"]["F_opt_from_canonical"], base["F"]["F_opt_from_ZCA"])
    base["F_opt_basis_used"] = f"descent from {best_from}"
    base["canonical_lam_range"] = [float(mdl["lam"].min()), float(mdl["lam"].max())]
    U_haar = haar(200, np.random.default_rng([seed, 99]))
    base["F"]["Haar_sample_basis"] = F_val(U_haar, mdl["A"], mdl["b"])
    base["basis_stats"] = {k: basis_stats(U, mdl) for k, U in
                           (("canonical", np.eye(200)), ("P_eigen", mdl["U_eig"]),
                            ("ZCA", mdl["U_zca"]), ("F_opt", U_best))}
    # where does the F-optimal basis put its weight?  fraction of each optimal vector in the observed block
    obs_frac = np.sum(U_best[:50, :] ** 2, axis=0)
    base["F_opt_vectors_observed_block_fraction"] = dict(
        min=float(obs_frac.min()), median=float(np.median(obs_frac)), max=float(obs_frac.max()),
        mean=float(obs_frac.mean()))
    check("F_canonical_equals_trAB", abs(base["F"]["canonical"] - base["trAB_canonical"]) < 1e-9 * base["trAB_canonical"],
          f"{base['F']['canonical']:.12g} vs {base['trAB_canonical']:.12g}")
    check("commutator_nonzero_base", base["commutator_rel"] > 1e-3, f"{base['commutator_rel']:.4f}")
    check("F_opt_below_canonical", base["F"]["F_opt"] < 0.9 * base["F"]["canonical"],
          f"{base['F']['F_opt']:.4f} vs {base['F']['canonical']:.4f}")
    check("F_opt_not_above_ZCA", base["F"]["F_opt"] <= base["F"]["ZCA"] + 1e-9,
          f"{base['F']['F_opt']:.4f} vs ZCA {base['F']['ZCA']:.4f}")
    if base["haar_samples"] >= 2:
        z = abs(base["haar_sampled_mean"] - base["haar_mean_exact"]) / max(base["haar_sampled_se"], 1e-12)
        check("haar_sampled_mean_matches_exact", z < 5.0, f"z = {z:.2f}")
    # Riemannian gradient finite-difference check at canonical
    fg, Om = F_grad(np.eye(200), mdl["A"], mdl["b"])
    Xi = np.random.default_rng([seed, 7]).standard_normal((200, 200))
    Xi = 0.5 * (Xi - Xi.T)
    h = 1e-6
    fd = (F_val(cayley(np.eye(200), Xi, -h), mdl["A"], mdl["b"])
          - F_val(cayley(np.eye(200), Xi, h), mdl["A"], mdl["b"])) / (2 * h)
    an = float(np.sum(Om * Xi))
    check("riemannian_gradient_fd", abs(fd - an) < 1e-5 * max(1.0, abs(an)), f"fd {fd:.8g} analytic {an:.8g}")
    timing["base_F_and_descent_s"] = time.time() - t0
    results["base"] = base

    # ---------------- main Monte Carlo ----------------
    t0 = time.time()
    Ns = [4, 8, 16, 32, 64]
    bases = {"canonical": ZBasis(np.eye(200), mdl), "P_eigen": ZBasis(mdl["U_eig"], mdl),
             "ZCA": ZBasis(mdl["U_zca"], mdl), "F_opt": ZBasis(U_best, mdl),
             "Haar_sample": ZBasis(U_haar, mdl)}
    # implementation check of the vectorised formula against the direct formula
    rchk = np.random.default_rng([seed, 5])
    Dchk = rchk.chisquare(7, size=(3, 200)) / 7
    worst = 0.0
    for k, bs in bases.items():
        X, Y = z_risk(bs, Dchk)
        for i in range(3):
            for a in (1.0, 1.125):
                dr = z_risk_direct(bs.U, mdl, Dchk[i], a)
                worst = max(worst, abs(a * X[i] + Y[i] - dr) / dr)
    check("vectorised_risk_matches_direct_formula", worst < 1e-10, f"max rel diff {worst:.2e}")
    kf = base["trWC"]
    mc = mc_z(bases, mdl, Ns, reps, seed + 1)
    summ = summarise(mc, kf, Ns)
    Fd = {"canonical": base["F"]["canonical"], "P_eigen": base["F"]["P_eigen"], "ZCA": base["F"]["ZCA"],
          "F_opt": base["F"]["F_opt"], "Haar_sample": base["F"]["Haar_sample_basis"]}
    lead = leading(mdl, Fd, Ns)
    results["mc_main"] = dict(Ns=Ns, reps=reps, kf_mse=kf, summary=summ, leading_order=lead)
    for N in Ns:
        s = summ[str(N)]["ensemble_mean"]
        dz = s["F_opt"]["diff_vs_canonical"] / s["F_opt"]["se_diff"]
        check(f"F_opt_beats_canonical_N{N}", dz < -3.0, f"diff {s['F_opt']['diff_vs_canonical']:.5f}, z {dz:.1f}")
    # leading-order agreement at the largest N (excess within 15%)
    s64 = summ["64"]["ensemble_mean"]
    for k in ("canonical", "F_opt"):
        ratio = (s64[k]["mean"] - 1.0) / (lead["64"][k]["ensemble_mean"] - 1.0)
        check(f"leading_order_excess_agrees_N64_{k}", abs(ratio - 1.0) < 0.15, f"MC/leading excess ratio {ratio:.4f}")
    timing["mc_main_s"] = time.time() - t0

    # ---------------- equivalent ensemble size ----------------
    t0 = time.time()
    Ngrid = [2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64]
    reps_eq = 100 if q else 1500
    mc_eq = mc_z({"F_opt": bases["F_opt"], "ZCA": bases["ZCA"]}, mdl, Ngrid, reps_eq, seed + 2)
    eqres = {}
    for kb in ("F_opt", "ZCA"):
        exc = [float(np.mean(((1 + 1 / N) * mc_eq[N][kb][0] + mc_eq[N][kb][1]) / kf) - 1.0) for N in Ngrid]
        rowd = {}
        for N in (8, 16, 32, 64):
            tgt = summ[str(N)]["ensemble_mean"]["canonical"]["mean"] - 1.0
            neq, how = eq_ensemble_size(Ngrid, exc, tgt)
            # leading-order version: tr(WC^2)/n + 2F/(n-1) = canonical leading excess
            tl = lead[str(N)]["canonical"]["ensemble_mean"] - 1.0
            Fk = Fd[kb]
            f = lambda n: (base["trWC2"] / n + 2 * Fk / (n - 1)) / kf - tl
            nl = brentq(f, 1.0001, 1e6) if f(1.0001) > 0 else None
            rowd[str(N)] = {"N_equivalent_MC": neq, "how": how, "N_equivalent_leading_order": nl}
        eqres[kb] = {"grid": Ngrid, "excess_on_grid": exc, "equivalent": rowd}
    results["equivalent_ensemble_size"] = dict(reps=reps_eq, results=eqres,
                                               note="N_eq such that the basis at N_eq matches canonical at N (ensemble-mean update)")
    timing["equivalent_N_s"] = time.time() - t0

    # ---------------- prior-normalised loss W = I ----------------
    t0 = time.time()
    mdlI = build_model(200, 50, rho0, "prior")
    AI, bI = mdlI["A"], mdlI["b"]
    FI_can = F_val(np.eye(200), AI, bI)
    fgI, OmI = F_grad(np.eye(200), AI, bI)
    starts = {"ZCA": mdlI["U_zca"], "P_eigen": mdlI["U_eig"],
              "Haar1": haar(200, np.random.default_rng([seed, 11]))}
    if not q:
        starts["Haar2"] = haar(200, np.random.default_rng([seed, 12]))
    dres = {}
    for k, U0 in starts.items():
        _, inf = descend(U0, AI, bI, max_iter=dit_s, label=f"W=I from {k}")
        dres[k] = {kk: inf[kk] for kk in ("F_start", "F_final", "iterations", "converged", "grad_norm_final")}
    minF = min(v["F_final"] for v in dres.values())
    check("W_eq_I_canonical_stationary", float(np.linalg.norm(OmI)) < 1e-10 * max(FI_can, 1.0),
          f"grad norm {float(np.linalg.norm(OmI)):.2e}")
    check("W_eq_I_no_basis_below_canonical", minF >= FI_can - 1e-8 * FI_can,
          f"min descent F {minF:.10g} vs canonical {FI_can:.10g}")
    kfI = float(np.sum(mdlI["c"]))
    basesI = {"canonical": ZBasis(np.eye(200), mdlI), "P_eigen": ZBasis(mdlI["U_eig"], mdlI),
              "ZCA": ZBasis(mdlI["U_zca"], mdlI)}
    NsI = [8, 32]
    mcI = mc_z(basesI, mdlI, NsI, 200 if q else 2000, seed + 3)
    summI = summarise(mcI, kfI, NsI)
    results["prior_normalised"] = dict(F_canonical=FI_can, F_ZCA=F_val(mdlI["U_zca"], AI, bI),
                                       F_P_eigen=F_val(mdlI["U_eig"], AI, bI),
                                       canonical_grad_norm=float(np.linalg.norm(OmI)),
                                       descents=dres, min_descent_F=minF, kf_mse=kfI, mc=summI)
    for N in NsI:
        s = summI[str(N)]["ensemble_mean"]["ZCA"]
        check(f"W_eq_I_ZCA_not_better_N{N}", s["diff_vs_canonical"] > -3 * s["se_diff"],
              f"diff {s['diff_vs_canonical']:.5f} se {s['se_diff']:.5f}")
    timing["prior_normalised_s"] = time.time() - t0

    # ---------------- sensitivity ----------------
    t0 = time.time()
    sens = {}
    cfgs = [("rho_x0.5", 0.5 * rho0, 50), ("rho_x1.5", 1.5 * rho0, 50), ("Ny25", rho0, 25), ("Ny100", rho0, 100)]
    for name, rho, ny in cfgs:
        md = build_model(200, ny, rho, "physical")
        ca, Uo1, Uo2 = config_analysis(md, rng, dit_s, 0)
        Ub = Uo2 if ca["F"]["F_opt_from_ZCA"] < ca["F"]["F_opt_from_canonical"] else Uo1
        ca["F"]["F_opt"] = min(ca["F"]["F_opt_from_canonical"], ca["F"]["F_opt_from_ZCA"])
        for kk in ("descent_from_canonical", "descent_from_ZCA"):
            ca[kk] = {x: ca[kk][x] for x in ("F_start", "F_final", "iterations", "converged", "grad_norm_final")}
        bs = {"canonical": ZBasis(np.eye(200), md), "P_eigen": ZBasis(md["U_eig"], md),
              "ZCA": ZBasis(md["U_zca"], md), "F_opt": ZBasis(Ub, md)}
        Nss = [8, 32]
        ms = mc_z(bs, md, Nss, 150 if q else 2000, seed + 4)
        ca["mc"] = summarise(ms, ca["trWC"], Nss)
        ca["leading_order"] = leading(md, {k: ca["F"][k] for k in ("canonical", "P_eigen", "ZCA", "F_opt")}, Nss)
        ca["rho"] = rho
        ca["ny"] = ny
        ca["prior_min_eigenvalue"] = md["ev_min_raw"]
        sens[name] = ca
        check(f"sensitivity_{name}_F_opt_beats_canonical_N8",
              ca["mc"]["8"]["ensemble_mean"]["F_opt"]["diff_vs_canonical"] < 0,
              f"{ca['mc']['8']['ensemble_mean']['F_opt']['diff_vs_canonical']:.5f}")
    results["sensitivity"] = sens
    timing["sensitivity_s"] = time.time() - t0

    # ---------------- estimated (pilot) basis ----------------
    t0 = time.time()
    # pilot configurations: (label, Np, GC half-width or None, ridge); one pilot draw per Np,
    # shared across regularisations (paired); all independent of the analysis ensembles
    pilot_cfgs = [("Np50_gc0.2_ridge1e-6", 50, 0.2, 1e-6), ("Np50_gc0.2_ridge1e-4", 50, 0.2, 1e-4),
                  ("Np50_gc0.1_ridge1e-3", 50, 0.1, 1e-3),
                  ("Np200_gc0.2_ridge1e-6", 200, 0.2, 1e-6), ("Np200_gc0.2_ridge1e-4", 200, 0.2, 1e-4),
                  ("Np200_gc0.1_ridge1e-3", 200, 0.1, 1e-3),
                  ("Np2000_raw_noridge", 2000, None, 0.0)]
    if q:
        pilot_cfgs = [pilot_cfgs[0], pilot_cfgs[3]]
    Ns_est = [8, 16]
    reps_est = 100 if q else 800
    oracle = {"canonical": (mdl["V"].T, mdl["Phalf"] @ mdl["V"]),
              "ZCA": (np.eye(200), mdl["Phalf"]),
              "F_opt": (U_best.T @ mdl["V"].T, mdl["Phalf"] @ mdl["V"] @ U_best)}
    # Tz = T P^{1/2} acts on the contrasts; oracle canonical T = V^T P^{-1/2} gives Tz = V^T
    est = {}
    est_info = {}
    for label, Np, gch, ridge in pilot_cfgs:
        prng = np.random.default_rng([seed, 1000 + Np])
        bset, info = estimated_bases(mdl, Np, prng, gch, ridge, 200 if q else 800)
        for k, v in bset.items():
            est[f"{k}|{label}"] = v
            info[f"{k}_rms_offdiag_true_correlation"] = rms_offdiag(v[0])
        info["descent"] = {x: info["descent"][x] for x in ("F_start", "F_final", "iterations", "converged")}
        est_info[label] = info
    allb = dict(oracle)
    allb.update(est)
    est_res = {}
    for N in Ns_est:
        nu = N - 1
        erng = np.random.default_rng([seed, 2000 + N])
        a_list = (1.0 + 1.0 / N,)
        vals = {k: [] for k in allb}
        done = 0
        while done < reps_est:
            n = min(50, reps_est - done)
            Zc = erng.standard_normal((n, 200, nu))
            for k, (Tz, Tinv) in allb.items():
                vals[k].append(phys_risk_batch(Tz, Tinv, mdl, Zc, a_list)[a_list[0]])
            done += n
        vals = {k: np.concatenate(v) / kf for k, v in vals.items()}
        n = reps_est
        ref = vals["canonical"]
        gain_oracle = float(ref.mean() - vals["F_opt"].mean())
        row = {}
        for k, v in vals.items():
            e = {"mean": float(v.mean()), "se": float(v.std(ddof=1) / np.sqrt(n))}
            if k != "canonical":
                d = v - ref
                e["diff_vs_oracle_canonical"] = float(d.mean())
                e["se_diff"] = float(d.std(ddof=1) / np.sqrt(n))
                e["fraction_of_oracle_F_opt_gain"] = float(-d.mean() / gain_oracle) if gain_oracle > 0 else None
            row[k] = e
        est_res[str(N)] = row
        # consistency of the physical-space route with the z-space shortcut (oracle bases)
        for kb in ("canonical", "F_opt", "ZCA"):
            zc = summ[str(N)]["ensemble_mean"][kb]
            zz = abs(row[kb]["mean"] - zc["mean"]) / np.hypot(row[kb]["se"], zc["se"])
            check(f"physical_route_matches_z_route_{kb}_N{N}", zz < 4.0, f"z = {zz:.2f}")
    results["estimated_basis"] = dict(configs=[{"label": c[0], "Np": c[1], "gc_half_width": c[2], "ridge": c[3]}
                                               for c in pilot_cfgs],
                                      Ns=Ns_est, reps=reps_est, pilot_info=est_info, results=est_res,
                                      note=("one pilot sample per Np, independent of the analysis ensembles; "
                                            "P_hat = GC o S_pilot + ridge I (GC half-width in domain units, support twice that); "
                                            "est_canonical, est_zca, est_F_opt are the canonical, symmetric and F-descended "
                                            "whitening transforms of P_hat; analysis ensembles simulated in physical space; "
                                            "risks normalised by the true KF MSE"))
    timing["estimated_basis_s"] = time.time() - t0

    # ---------------- Fourier corollary ----------------
    t0 = time.time()
    four = fourier_check(64, rho0, 1.0, rng, 300 if q else 2000)
    # negative control: r below the smallest prior mode variance -> no violating pair
    rctl = 0.5 * four["p_min"]
    four_ctrl = fourier_check(64, rho0, rctl, rng, 50)
    results["fourier_corollary"] = dict(main=four, control_r_below_all_modes=four_ctrl)
    check("fourier_basis_diagonalises_circulant", four["fourier_diag_error"] < 1e-10, f"{four['fourier_diag_error']:.2e}")
    check("fourier_spectrum_positive", four["spectrum_positive"], f"p_min {four['p_min']:.3e}")
    check("fourier_violations_exist_with_modes_below_r", four["violating_pairs_both_below_r"] > 0,
          f"{four['violating_pairs_both_below_r']} of {four['distinct_pairs_both_below_r']}")
    check("fourier_every_distinct_pair_below_r_violates", four["every_distinct_pair_below_r_violates"], "")
    for tag in ("best_pair", "best_pair_both_below_r"):
        if tag in four:
            dd = abs(four[tag]["predicted_45deg_reduction"] - four[tag]["numerical_45deg_reduction"])
            check(f"fourier_45deg_reduction_formula_{tag}", dd < 1e-12 * max(1.0, four["F_canonical"]), f"diff {dd:.2e}")
        else:
            check(f"fourier_45deg_reduction_formula_{tag}", False, "no violating pair found")
    check("fourier_descent_lowers_F", four.get("F_after_descent_from_perturbed_canonical", np.inf) < four["F_canonical"] - 1e-6,
          f"{four.get('F_after_descent_from_perturbed_canonical')} vs {four['F_canonical']}")
    check("fourier_control_no_violation", four_ctrl["violating_pairs"] == 0, f"{four_ctrl['violating_pairs']}")
    we = four["worked_example_r1_p02_p05"]
    check("fourier_worked_example", abs(we["reduction_45deg"] - 1.0 / 144.0) < 1e-14, f"{we['reduction_45deg']:.6f}")
    timing["fourier_s"] = time.time() - t0

    # ---------------- figure ----------------
    t0 = time.time()
    make_figure(args.out, Ns, summ, lead, "kf")
    timing["figure_s"] = time.time() - t0

    timing["total_s"] = time.time() - t_start
    results["timing_s"] = timing
    results["checks"] = CHECKS
    results["failures"] = FAILURES
    results["status"] = "PASS" if not FAILURES else "FAIL"
    with open(os.path.join(args.out, "results.json"), "w") as fh:
        json.dump(results, fh, indent=1, sort_keys=True)
    # timing-free copy for byte comparison between runs
    det = {k: v for k, v in results.items() if k not in ("timing_s", "optimize_flag", "platform")}
    with open(os.path.join(args.out, "results_deterministic.json"), "w") as fh:
        json.dump(det, fh, indent=1, sort_keys=True)
    print(json.dumps({"status": results["status"], "failures": FAILURES, "timing_s": timing}, indent=1))
    if FAILURES:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
