#!/usr/bin/env python3
"""Deterministic numerical checks for the two-dimensional example of the paper.

Model (manuscript Eq. (ex)): P = [[1,1],[1,2]], H = (1,0), R = 1, physical
squared loss, so W = [[2,1],[1,1]] in canonical coordinates and the Bayes risk
is 2.  Two independently coded expectations are evaluated: a one-dimensional
Laplace-transform integral and a tensor generalised Gauss--Laguerre rule over
the two independent sample variances.  Both rest on the same conditional-risk
identity and coupling representation, so their agreement checks the quadrature
and its implementation, not the model; the model itself is checked by the
Monte Carlo of the real ensemble procedure in theory_checks.py.

These are floating-point checks, not interval enclosures, and not the proof of
any theorem.

No assert statements are used.  Every check is recorded, the status is computed
from the checks, and the script exits with code 1 on any failure.  All outputs
(JSON and the figure) are written to a NEW directory given by --out; the script
refuses to overwrite an existing path and leaves no temporary directory behind.

Usage:  python3 -I -B numerical_checks.py --out DIR
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import sys
import tempfile
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import gammaln, roots_genlaguerre

W_PHYS = np.array([[2.0, 1.0], [1.0, 1.0]])
THETA_STAR = -math.atan(2.0) / 4
THETA_RATIONAL = -2 * math.atan(1 / 8)

CHECKS: list[dict] = []
FAILURES: list[str] = []


def check(name: str, ok: bool, **info) -> bool:
    ok = bool(ok)
    rec = {"name": name, "ok": ok}
    rec.update({k: (float(v) if isinstance(v, (np.floating,)) else v) for k, v in info.items()})
    CHECKS.append(rec)
    if not ok:
        FAILURES.append(name)
    return ok


def risk_laplace(n: int, theta: float, known_mean: bool = False, W=W_PHYS) -> float:
    """Expected risk via one-dimensional Laplace integrals (beta = 1)."""
    shape = (n - 1) / 2
    a = 1.0 if known_mean else 1 + 1 / n
    c, s = math.cos(theta), math.sin(theta)
    cc, ss, cs = c * c, s * s, c * s

    def integrand(t: float, mode: int) -> float:
        ra = 1 / (1 + cc * t / shape)
        rb = 1 / (1 + ss * t / shape)
        base = math.exp(-t - shape * (math.log1p(cc * t / shape)
                                      + math.log1p(ss * t / shape)))
        if mode == 0:
            return base
        if mode == 1:
            return t * base
        if mode == 2:
            return cs * base * (ra - rb)
        if mode == 3:
            return cs * t * base * (ra - rb)
        return cs * cs * t * base * ((ra - rb) ** 2 + (ra * ra + rb * rb) / shape)

    vals = [quad(integrand, 0, np.inf, args=(j,), epsabs=3e-13,
                 epsrel=3e-13, limit=250)[0] for j in range(5)]
    a0, b0, c0, c1, v2 = vals
    ek1, ek2 = 1 - a0, c0
    ek11, ek12, ek22 = 1 - 2 * a0 + b0, c0 - c1, v2
    return float(a * np.trace(W) - 2 * a * (W[0, 0] * ek1 + W[0, 1] * ek2)
                 + (a + 1) * (W[0, 0] * ek11 + 2 * W[0, 1] * ek12 + W[1, 1] * ek22))


def risk_laguerre(n: int, theta: float, order: int = 160, W=W_PHYS) -> float:
    """Direct two-dimensional expectation by a tensor generalised Gauss--Laguerre rule."""
    shape = (n - 1) / 2
    nodes, weights = roots_genlaguerre(order, shape - 1)
    weights = weights / weights.sum()
    d, e = nodes[:, None] / shape, nodes[None, :] / shape
    c, s = math.cos(theta), math.sin(theta)
    u = c * c * d + s * s * e
    v = c * s * (d - e)
    k1, k2 = u / (1 + u), v / (1 + u)
    a = 1 + 1 / n
    loss = (a * np.trace(W) - 2 * a * (W[0, 0] * k1 + W[0, 1] * k2)
            + (a + 1) * (W[0, 0] * k1 * k1 + 2 * W[0, 1] * k1 * k2
                         + W[1, 1] * k2 * k2))
    return float(weights @ loss @ weights)


def derivative_integral(n: int, beta: float = 1, known_mean: bool = False) -> float:
    """2 beta E[(D-1)(D-a)/(1+beta D)^2] (Theorem 2 with W_12 = 1)."""
    shape = (n - 1) / 2
    a = 1.0 if known_mean else 1 + 1 / n

    def f(d: float) -> float:
        if d == 0:
            return 0.0
        density = math.exp(shape * math.log(shape) - gammaln(shape)
                           + (shape - 1) * math.log(d) - shape * d)
        return 2 * beta * (d - 1) * (d - a) / (1 + beta * d) ** 2 * density

    return sum(quad(f, lo, hi, epsabs=2e-12, epsrel=2e-12, limit=300)[0]
               for lo, hi in [(0, 1), (1, np.inf)])


def f_leading(theta: float) -> float:
    """Leading objective F(theta) = (5 - cos 4 theta + 2 sin 4 theta)/16."""
    return (5 - math.cos(4 * theta) + 2 * math.sin(4 * theta)) / 16


def matrix_checks() -> dict:
    l = np.array([[1.0, 0.0], [1.0, 1.0]])
    p = np.array([[1.0, 1.0], [1.0, 2.0]])
    h = np.array([[1.0, 0.0]])
    check("L L^T = P", np.array_equal(l @ l.T, p))
    check("L^T L = W", np.array_equal(l.T @ l, W_PHYS))
    check("H L = e_1^T", np.array_equal(h @ l, h))
    check("det P = 1", abs(np.linalg.det(p) - 1.0) < 1e-15)
    c = np.diag([0.5, 1.0])
    b = np.eye(2) - c
    aa = c @ W_PHYS @ c
    max_error = 0.0
    for theta in np.linspace(-math.pi, math.pi, 73):
        co, si = math.cos(theta), math.sin(theta)
        u = np.array([[co, -si], [si, co]])
        direct = np.dot(np.diag(u.T @ b @ u), np.diag(u.T @ aa @ u))
        max_error = max(max_error, abs(direct - f_leading(theta)))
    check("leading-angle formula matches direct matrix evaluation (1e-14)", max_error < 1e-14,
          max_error=max_error)
    return {"leading_angle_formula_max_error": max_error}


def make_figure(rows: list[dict], outdir: Path) -> None:
    with tempfile.TemporaryDirectory(prefix="mplcfg-") as mpl_dir:
        os.environ["MPLCONFIGDIR"] = mpl_dir
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                             "axes.spines.top": False, "axes.spines.right": False,
                             "svg.hashsalt": "canon", "pdf.compression": 6})
        fig, axes = plt.subplots(2, 1, figsize=(6.6, 7.0), layout="constrained")
        colours = {2: "#A66532", 10: "#7C3F73", 100: "#277B78"}
        degrees = np.linspace(-40, 20, 181)
        for n, colour in colours.items():
            r0 = risk_laplace(n, 0)
            ys = [100 * (risk_laplace(n, math.radians(x)) - r0) / (r0 - 2) for x in degrees]
            axes[0].plot(degrees, ys, color=colour, label=f"N = {n}", linewidth=2)
        axes[0].axhline(0, color="#615A50", linewidth=0.8)
        axes[0].axvline(math.degrees(THETA_STAR), color="#B8A582", linestyle=":", linewidth=1.4)
        axes[0].set(xlabel="Rotation from canonical coordinates (degrees)",
                    ylabel="Change in excess risk over Bayes (%)",
                    title="(a) A negative rotation lowers expected error", ylim=(-16, 22))
        axes[0].legend(frameon=False)
        ns = np.array([r["N"] for r in rows])
        reduction = [r["rational_rotation_excess_risk_reduction_percent"] for r in rows]
        axes[1].plot(ns, reduction, "o-", color="#6C526F", linewidth=2,
                     label="Fixed rational rotation, −14.25°")
        # In this example the leading excess is 3/(2N) + 2F/(N-1) and 3/2 + 2F(0) = 2,
        # so the relative leading excess reduction equals 100 (F(0) - F(theta)).
        asymp = 100 * (f_leading(0) - f_leading(THETA_RATIONAL))
        axes[1].axhline(asymp, color="#A66532", linestyle="--",
                        label=f"Leading-order limit for this rotation: {asymp:.2f}%")
        axes[1].set(xscale="log", xlabel="Ensemble size N",
                    ylabel="Reduction in excess risk over Bayes (%)",
                    title="(b) The relative excess-risk gain persists as N grows")
        axes[1].legend(frameon=False, fontsize=9)
        fig.savefig(outdir / "risk_rotation.png", dpi=220, metadata={"Software": None})
        fig.savefig(outdir / "risk_rotation.pdf",
                    metadata={"Creator": None, "Producer": None, "CreationDate": None})
        plt.close(fig)
    os.environ.pop("MPLCONFIGDIR", None)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True, help="new output directory (must not exist)")
    args = ap.parse_args()
    outdir = Path(os.path.abspath(args.out))
    if outdir.exists():
        print(f"refusing to overwrite existing path {outdir}", file=sys.stderr)
        raise SystemExit(2)
    script_sha = hashlib.sha256(Path(os.path.abspath(__file__)).read_bytes()).hexdigest()

    rows = []
    crosschecks = []
    for n in [2, 3, 5, 10, 20, 50, 100, 200]:
        canonical = risk_laplace(n, 0)
        chosen = risk_laplace(n, THETA_RATIONAL)
        asymptotic_angle = risk_laplace(n, THETA_STAR)
        opt = minimize_scalar(lambda th: risk_laplace(n, th),
                              bounds=(-math.pi / 4, 0), method="bounded",
                              options={"xatol": 2e-12})
        derivative = derivative_integral(n)
        step = 1e-5
        finite_diff = (risk_laplace(n, step) - risk_laplace(n, -step)) / (2 * step)
        check(f"N={n}: analytic derivative matches central difference (3e-8)",
              abs(derivative - finite_diff) < 3e-8, analytic=derivative, fd=finite_diff)
        check(f"N={n}: rational rotation beats canonical", chosen < canonical,
              canonical=canonical, rotated=chosen)
        rows.append({"N": n, "canonical_risk": canonical,
                     "rational_rotation_risk": chosen,
                     "MSE_over_KF_canonical": canonical / 2,
                     "MSE_over_KF_rotated": chosen / 2,
                     "rational_rotation_total_risk_reduction_percent": 100 * (canonical - chosen) / canonical,
                     "rational_rotation_excess_risk_reduction_percent": 100 * (canonical - chosen) / (canonical - 2),
                     "asymptotic_angle_risk": asymptotic_angle,
                     "numerical_minimiser_degrees": math.degrees(opt.x),
                     "numerical_minimum_risk": opt.fun,
                     "analytic_derivative_numeric_integral": derivative,
                     "central_difference_derivative": finite_diff})
        for theta in [0.0, THETA_RATIONAL, THETA_STAR, 0.2]:
            r1 = risk_laplace(n, theta)
            r2 = risk_laguerre(n, theta, 160)
            err = abs(r1 - r2)
            check(f"N={n}, theta={theta:.6f}: Laplace vs Laguerre (2e-8)", err < 2e-8,
                  laplace=r1, laguerre=r2)
            crosschecks.append({"N": n, "theta": theta, "absolute_difference": err})

    beta_checks = [{"N": n, "beta": beta, "derivative": derivative_integral(n, beta)}
                   for n in [2, 3, 10, 100] for beta in [0.01, 0.1, 1.0, 10.0, 100.0]]
    for r in beta_checks:
        check(f"N={r['N']}, beta={r['beta']}: derivative positive", r["derivative"] > 0,
              value=r["derivative"])
    convergence = []
    for n in [100, 1000, 10000, 100000]:
        for theta in [0, THETA_STAR]:
            observed = n * (risk_laplace(n, theta) - 2)
            predicted = 1.5 + 2 * f_leading(theta)
            convergence.append({"N": n, "theta": theta,
                                "scaled_excess": observed, "limit": predicted})
    last = [c for c in convergence if c["N"] == 100000]
    for c in last:
        check(f"N=1e5, theta={c['theta']:.4f}: scaled excess within 1e-4 of limit",
              abs(c["scaled_excess"] - c["limit"]) < 1e-4, observed=c["scaled_excess"],
              limit=c["limit"])

    # Negative controls.
    neg = {}
    for n in (3, 10):
        up = risk_laplace(n, abs(THETA_RATIONAL))
        neg[f"wrong_sign_rotation_N{n}"] = {"canonical": risk_laplace(n, 0), "wrong_sign": up}
        check(f"NEG N={n}: rotation of the wrong sign increases risk", up > risk_laplace(n, 0),
              canonical=risk_laplace(n, 0), wrong_sign=up)
    eye = np.eye(2)
    for n in (3, 10):
        fd_eye = (risk_laplace(n, 1e-5, W=eye) - risk_laplace(n, -1e-5, W=eye)) / 2e-5
        neg[f"W_identity_first_variation_N{n}"] = fd_eye
        check(f"NEG N={n}: W=I first variation is zero (|fd|<1e-8)", abs(fd_eye) < 1e-8,
              fd=fd_eye)
    bad = derivative_integral(3) * 2
    fd3 = rows[1]["central_difference_derivative"]
    check("NEG doubled derivative constant is rejected", abs(bad - fd3) > 1e-3, claimed=bad, fd=fd3)

    report = {"script": "numerical_checks.py", "script_sha256": script_sha,
              "python": sys.version.split()[0], "numpy": np.__version__,
              "scipy": scipy.__version__, "platform": platform.platform(),
              "optimize_flag": sys.flags.optimize,
              "model": {"P": [[1, 1], [1, 2]], "H": [[1, 0]], "R": 1,
                        "W": W_PHYS.tolist(), "Bayes_risk": 2},
              "quadrature_scope": "floating-point deterministic checks, not rigorous interval enclosures",
              "theta_star_degrees": math.degrees(THETA_STAR),
              "rational_rotation": {"cos": "63/65", "sin": "-16/65",
                                    "theta_degrees": math.degrees(THETA_RATIONAL)},
              "identities": matrix_checks(), "risk_rows": rows,
              "independent_quadrature_comparisons": crosschecks,
              "independent_quadrature_max_difference": max(x["absolute_difference"] for x in crosschecks),
              "positive_derivative_checks": beta_checks,
              "large_N_convergence": convergence,
              "negative_controls": neg}
    outdir.mkdir(parents=True)
    make_figure(rows, outdir)
    for name in ("risk_rotation.png", "risk_rotation.pdf"):
        check(f"figure {name} written", (outdir / name).stat().st_size > 1000)
    report["checks"] = CHECKS
    report["n_checks"] = len(CHECKS)
    report["failures"] = FAILURES
    report["status"] = "PASS" if not FAILURES else "FAIL"
    (outdir / "numerical_checks.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "n_checks": len(CHECKS),
                      "failures": FAILURES,
                      "independent_quadrature_max_difference": report["independent_quadrature_max_difference"],
                      "rows": [(r["N"], round(r["canonical_risk"], 6), round(r["rational_rotation_risk"], 6),
                                round(r["numerical_minimiser_degrees"], 3)) for r in rows]}, indent=1))
    if FAILURES:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
