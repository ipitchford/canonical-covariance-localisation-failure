#!/usr/bin/env python3
"""Exact rational-arithmetic certificates for the canonical-localisation paper.

The script checks, with fractions.Fraction (no floating point), the finite
arithmetic used in the proofs:

* the polynomial-square lower bound c >= 460/773 for c = E[1/(1+D)], D ~ Exp(1);
* the aligned-example curvature bounds -56/773 (known mean) and -508/773
  (ensemble mean, N = 3);
* the global second-derivative bound H_0 = 18*66 + 45*488 = 23148;
* the rational finite-angle certificate t_nu = g_nu/(4 H_0) and the risk-gap
  bound -g_nu^2/(8 H_0), with exact orthogonality of the Cayley rotation;
* negative controls: deliberately wrong constants must be rejected.

It does not formalise the probabilistic or differentiation arguments; those
are proved in the manuscript.

No assert statements are used. Every check is recorded; the status is computed
from the recorded checks; the script exits with code 1 if any check fails.
Results are written to a NEW directory given by --out (the script refuses to
overwrite an existing path). Python 3 standard library only.

Usage:  python3 -I -B exact_certificates.py --out DIR
"""

import argparse
import hashlib
import json
import os
import platform
import sys
from fractions import Fraction as F
from math import factorial

CHECKS = []
FAILURES = []


def check(name, ok, **info):
    ok = bool(ok)
    rec = {"name": name, "ok": ok}
    rec.update({k: str(v) for k, v in info.items()})
    CHECKS.append(rec)
    if not ok:
        FAILURES.append(name)
    return ok


def exponential_integral(poly):
    """Integral_0^infinity poly(x) exp(-x) dx, exactly (coefficients low to high)."""
    return sum((coef * factorial(i) for i, coef in enumerate(poly)), F(0))


def multiply(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def hessian_bound(moments):
    """18 E[T(1+T)^2] + 45 E[T^2(1+T)^2] from bounds on E T^k, k = 0..4."""
    e1 = moments[1] + 2 * moments[2] + moments[3]
    e2 = moments[2] + 2 * moments[3] + moments[4]
    return e1, e2, 18 * e1 + 45 * e2


def certificate(nu, h0):
    g = F(2 * nu, 2 * nu * nu + 11 * nu + 12)
    t = g / (4 * h0)
    gap = -g * g / (8 * h0)
    cosine = (1 - t * t) / (1 + t * t)
    sine = -2 * t / (1 + t * t)
    return g, t, gap, cosine, sine


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True, help="new output directory (must not exist)")
    args = ap.parse_args()
    outdir = os.path.abspath(args.out)
    if os.path.exists(outdir):
        print(f"refusing to overwrite existing path {outdir}", file=sys.stderr)
        raise SystemExit(2)

    with open(os.path.abspath(__file__), "rb") as fh:
        script_sha = hashlib.sha256(fh.read()).hexdigest()

    # 1. Polynomial-square certificate for c = E[1/(1+D)], D ~ Exp(1).
    p = [F(713, 773), F(-413, 773), F(113, 773), F(-13, 773), F(1, 1546)]
    integral_p = exponential_integral(p)
    integral_weighted_square = exponential_integral(multiply([F(1), F(1)], multiply(p, p)))
    lower_c = 2 * integral_p - integral_weighted_square
    check("integral of p e^{-x} equals 460/773", integral_p == F(460, 773), value=integral_p)
    check("integral of (1+x) p^2 e^{-x} equals 460/773",
          integral_weighted_square == F(460, 773), value=integral_weighted_square)
    check("lower bound c >= 460/773", lower_c == F(460, 773), value=lower_c)

    # 2. Aligned-example curvature: R''(0) = 28a + 60 - (50a + 98)c, decreasing in c.
    curvature = {}
    for name, a, expected in [("known_mean", F(1), F(-56, 773)),
                              ("ensemble_mean_N3", F(4, 3), F(-508, 773))]:
        upper = 28 * a + 60 - (50 * a + 98) * lower_c
        curvature[name] = str(upper)
        check(f"curvature upper bound {name} is negative", upper < 0, value=upper)
        check(f"curvature upper bound {name} equals {expected}", upper == expected, value=upper)

    # 3. Global second-derivative bound from moment bounds of T = D1 + D2.
    moment_bounds = [F(1), F(2), F(8), F(48), F(384)]
    e1, e2, h0 = hessian_bound(moment_bounds)
    check("E[T(1+T)^2] bound equals 66", e1 == 66, value=e1)
    check("E[T^2(1+T)^2] bound equals 488", e2 == 488, value=e2)
    check("H_0 equals 23148", h0 == 23148, value=h0)
    # The moment bounds dominate the exact gamma moments for nu = 1..200.
    # T ~ Gamma(shape nu, scale 2/nu): E T^k = (2/nu)^k nu (nu+1) ... (nu+k-1).
    worst = F(0)
    for nu in range(1, 201):
        for k in range(1, 5):
            exact = F(2, nu) ** k
            for i in range(k):
                exact *= nu + i
            worst = max(worst, exact / moment_bounds[k])
    check("gamma moment bounds 2,8,48,384 hold for nu = 1..200", worst <= 1, max_ratio=worst)

    # 4. Rational finite-angle certificate.
    finite_angle = {}
    for nu in [1, 2, 3, 9, 99]:
        g, t, gap, cosine, sine = certificate(nu, h0)
        ok_angle = 0 < t < 1
        ok_gap = gap < 0
        ok_orth = cosine * cosine + sine * sine == 1
        check(f"nu={nu}: 0 < t < 1", ok_angle, t=t)
        check(f"nu={nu}: risk gap bound negative", ok_gap, gap=gap)
        check(f"nu={nu}: Cayley rotation exactly orthogonal", ok_orth)
        # Taylor step: -g*delta + H0*delta^2/2 <= -g^2/(8 H0) at the extreme
        # deltas t and 2t allowed by t <= delta <= 2t.
        worst_taylor = max(-g * d + h0 * d * d / 2 for d in (t, 2 * t))
        check(f"nu={nu}: Taylor bound at delta in {{t,2t}} below gap", worst_taylor <= gap,
              worst=worst_taylor, gap=gap)
        finite_angle[str(nu + 1)] = {"gradient_lower": str(g), "cayley_parameter": str(t),
                                     "risk_gap_upper": str(gap),
                                     "orthogonality_identity": "exact" if ok_orth else "FAILED"}
    check("N=3 Cayley parameter equals 1/972216",
          finite_angle["3"]["cayley_parameter"] == "1/972216",
          value=finite_angle["3"]["cayley_parameter"])
    check("N=3 risk gap bound equals -1/20416536",
          finite_angle["3"]["risk_gap_upper"] == "-1/20416536",
          value=finite_angle["3"]["risk_gap_upper"])

    # 5. Negative controls: wrong constants must be rejected by the same tests.
    neg = {}
    p_bad = list(p)
    p_bad[0] = F(714, 773)
    bad_c = 2 * exponential_integral(p_bad) - exponential_integral(
        multiply([F(1), F(1)], multiply(p_bad, p_bad)))
    neg["perturbed_polynomial_lower_c"] = str(bad_c)
    check("NEG perturbed polynomial does not reproduce 460/773", bad_c != F(460, 773), value=bad_c)
    # Halving the fourth-moment bound (384 -> 192) would shrink H_0; the exact
    # gamma moment at nu = 1 (E T^4 = 384) must expose that bound as false.
    bad_bounds = [F(1), F(2), F(8), F(48), F(192)]
    exact_t4_nu1 = F(2) ** 4 * 1 * 2 * 3 * 4
    neg["halved_fourth_moment_H0"] = str(hessian_bound(bad_bounds)[2])
    check("NEG halved fourth-moment bound (smaller H_0) is violated at nu=1",
          exact_t4_nu1 > bad_bounds[4], exact=exact_t4_nu1, claimed=bad_bounds[4])
    g, t, gap, _, _ = certificate(2, h0)
    worst_double = max(-g * d + h0 * d * d / 2 for d in (t, 2 * t))
    check("NEG a gap claimed twice as large is rejected", worst_double > 2 * gap,
          worst=worst_double, claimed=2 * gap)
    wrong_curv = 28 * F(4, 3) + 60 - (50 * F(4, 3) + 98) * F(459, 773)
    check("NEG weaker lower bound 459/773 does not give -508/773",
          wrong_curv != F(-508, 773), value=wrong_curv)

    status = "PASS" if not FAILURES else "FAIL"
    report = {
        "script": "exact_certificates.py",
        "script_sha256": script_sha,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "optimize_flag": sys.flags.optimize,
        "scope": ("exact rational arithmetic for the finite constants used in the proofs; "
                  "the probabilistic and differentiation arguments are proved in the manuscript"),
        "polynomial_coefficients": [str(x) for x in p],
        "integral_p": str(integral_p),
        "integral_weighted_square": str(integral_weighted_square),
        "lower_bound_c": str(lower_c),
        "curvature_upper_bounds": curvature,
        "hessian_bound": str(h0),
        "finite_angle_examples_by_ensemble_size": finite_angle,
        "negative_controls": neg,
        "checks": CHECKS,
        "n_checks": len(CHECKS),
        "failures": FAILURES,
        "status": status,
    }
    os.makedirs(outdir)
    with open(os.path.join(outdir, "exact_certificates.json"), "w") as fh:
        json.dump(report, fh, indent=2)
        fh.write("\n")
    print(f"status={status} checks={len(CHECKS)} failures={len(FAILURES)}")
    for name in FAILURES:
        print("FAILED:", name)
    if FAILURES:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
