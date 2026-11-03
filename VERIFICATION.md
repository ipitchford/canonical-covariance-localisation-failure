# Verification ledger

Authoritative statement: `manuscript.tex` and its build, `canonical_localisation.pdf` (3 November 2026). Numbering below refers to that PDF.

Evidence types:
- **Proof**: a written argument in the paper. No proof is machine-verified.
- **Exact**: rational arithmetic in `scripts/exact_certificates.py`.
- **Numerical**: floating-point quadrature or seeded Monte Carlo, not interval-certified.

Recorded outputs live in `results/<script>/`. Every script records its own SHA-256, computes its status from explicit checks and exits with code 1 on failure.

Check names inside `theory_checks.json` keep the labels of the 1 November draft:

| Check-name label | Current paper |
|---|---|
| "Lemma 1" | Lemma 3.1 |
| "Thm2" | Theorem 3.3 |
| "Prop 4" | Proposition 4.1 |
| "Table 1" | Table 3 |

## Claims, proofs and checks

| Claim | Where proved | Check (script → check name) | Recorded result |
|---|---|---|---|
| Exact conditional risk (2.2)–(2.3), both mean conventions | Section 2 | `theory_checks` MC: canonical risk vs Table 3 (real ensembles) | MC 2.82784 ± 0.00177 vs 2.831656 (within 4 SE) |
| Lemma 3.1: complete localisation is the best canonical taper with T₁₁ = 1 (PSD not needed) | Lemma 3.1 | `theory_checks` "MC Lemma 1 taper excess" (T₁₂ = 0.5, N = 3) | MC 0.05581 ± 0.00051 vs exact 0.05620 |
| Remark 3.2: covariance-error criteria cannot rank whitened bases | One-line computation | — | Proof only |
| Theorem 3.3: R′_N(0) = 2βW₁ⱼγ, γ > 0 for every N ≥ 2 | Theorem 3.3 | `numerical_checks`: analytic derivative vs central difference at 8 N (tol 3e-8); positive sign at 20 (N, β) pairs; `theory_checks` MC Thm 2 | All pass; MC 0.5921 ± 0.0024 vs exact central difference 0.5935 |
| Corollary 3.4: failure whenever Hᵀ is not an eigenvector of P | Corollary 3.4 | — | Proof only |
| Lemma 3.5: risk of inflated canonical tapers (any T₁₁ = τ > 0) | Lemma 3.5 | `theory_checks` X2 grid | Pass |
| Lemma 3.6(1): Ψ′ has at most one sign change (monotone likelihood ratio) | Lemma 3.6 | `theory_checks` "X2 Psi' at most one sign change" (24 cases) | Pass |
| Lemma 3.6(2)–(3): minimiser exists for N ≥ 3; for N = 2 iff aβ < 3 | Lemma 3.6 (Beta integrals; Stein identity for aβ ≥ 3) | `theory_checks` X2: existence pattern on grid (14 cases N = 2), asymptotic constant π(c−3)/2 (10), Stein identity J = 2(1−m₁)/s (15), Ψ > 1/β for aβ > 3 (4) | All pass |
| Theorem 3.7: a rotation with the **same** rescaling diag(τ₀,1,…,1) beats every canonical taper with any T₁₁ > 0; γ(τ*) > 0 at any stationary τ*. An unrescaled rotation need not (N = 3: 2.7453 vs infimum 2.7408; Remark 3.8) | Theorem 3.7 | `theory_checks` X2: γ(τ*) > 0 on grid, identity γ = E[ψ(ϱ−ϱ(a))] (rel. 1e-7), derivative vs FD (69 cases); MC of real ensembles | All pass; MC 0.5256 ± 0.0023 vs 0.5260; τ* = 2.37206 |
| Rational certificate (4.2), gap ≤ −1/20416536 at N = 3 | Section 4.1, Appendix A | `exact_certificates`: g_ν, t_ν, H₀ = 23148, Taylor bound, exact orthogonality, gamma-moment bounds ν ≤ 200 | 37 checks pass, including 4 negative controls |
| Proposition 4.1: exact curvature formula | Proposition 4.1 | `theory_checks` "Prop 4: exact small-step second difference equals analytic"; MC second difference (real ensembles) | Pass; MC −0.879 ± 0.018 vs exact −0.878 at step 0.2 (analytic −0.8652) |
| Aligned example: strict local maximum at N = 3, bounds −56/773, −508/773 | Section 4.2 | `exact_certificates` (460/773 certificate, both bounds); `theory_checks` X5 closed form | Pass |
| Remark 4.2: aligned curvature at other N | Numerical remark | `theory_checks` X5 (N = 2…1000) | Positive at N = 2 (known mean), negative otherwise; N²R″ → −8, −12 |
| Proposition 5.1: vector-observation first variation | Proposition 5.1 | `theory_checks` "X3 first variation FD" (18 cases) | Pass |
| Corollaries 5.2–5.3 | Corollaries | Follow from Proposition 5.1 | Proof only |
| Conjecture 5.4: Φ(x, y) > 0 for x > y > 0 | **Not proved** | `theory_checks` X3 Φ grid | Positive in all 660 cases; min Φ/(x−y) = 1.7e-8 |
| Lemmas 6.2–6.3, Theorem 6.1: uniform expansion with remainder O(ν^{−3/2}) | Section 6 (Schur-multiplier norm bound proved in Lemma 6.2) | `sh_testbed` "leading_order_excess_agrees_N64" | Pass. O(ν⁻²) is a sketch only (Remark 6.4). |
| Theorem 7.1: sharp leading-order criterion | Theorem 7.1 | `numerical_checks` leading-angle formula (73 angles, 1e-14); `sh_testbed` F_canonical = tr(AB), commutator, descent | Pass |
| Lemma 7.4, Corollary 7.6: large-N failure against every unit-diagonal PSD canonical taper | Section 7.1 | — | Proof only; N₀ not explicit |
| Remark 7.5: unit-diagonal PSD tapers in any whitened basis have leading term ≥ 2F(U) (Wishart moment + AM–GM) | Remark 7.5 | — | Proof only; the identity and inequality were checked independently in the proof audit |
| Remark 7.2: P = [[3,1],[1,3]], H = I, R = I: canonical leading-order optimal under physical loss (ordering product −28/3375) | Remark | Direct computation (numpy: minimum of F over angles equals tr(AB) = 0.276148) | Pass |
| Proposition 8.1 (scalar observation, W = I): strict local minimum modulo completion if ζ > 0, not a local minimum if ζ < 0; ζ = 0 unresolved | Proposition 8.1 | `numerical_checks` and `theory_checks` "W=I first variation zero" | Pass |
| Open problem 8.2: ζ > 0 | **Not proved** | `theory_checks` X4 grid (14 N, 31 β, both means) | Ratio > 1 everywhere; minimum 1.0020 |
| Section 9: closed-form angle θ* = −¼ arctan 2, 30.90 % and 7.73 % | Section 9 | `numerical_checks` | Pass |
| Table 3 (risks for the 2-D example) | Numerical | `numerical_checks`: Laplace vs Gauss–Laguerre (32 comparisons); MC rotated risk at N = 3 | Max difference 4.0e-14; MC 2.7445 ± 0.0024 vs 2.746111 |
| Section 10, Tables 5–6, Figure 2: testbed reconstruction | Numerical | `sh_testbed` (39 checks) | PASS; see `results/sh_testbed/results.json` |
| Corollary 10.1: Fourier bases under physical loss | Corollary 10.1 | `sh_testbed` fourier_* checks | Pass; 1058 + 671 violating pairs; 45° reduction matches to 1e-15 |

## Negative controls (must fail, and did)

**`exact_certificates`:**
- a perturbed polynomial does not reproduce 460/773;
- a halved fourth-moment bound is violated by the exact gamma moment;
- a gap claimed to be twice as large is rejected by the Taylor bound;
- a weaker lower bound does not give −508/773.

**`numerical_checks`:**
- a rotation of the wrong sign raises the risk (N = 3 and 10);
- W = I gives a zero first variation;
- a doubled derivative constant is rejected.

**`theory_checks`:**
- a wrong rotation sign raises the risk in the real-ensemble Monte Carlo (z > 4);
- a 2×γ constant is rejected;
- W = I gives a zero Monte Carlo derivative;
- the aligned threshold is detected below 4;
- a false curvature bound is rejected.

## Regeneration and stability

Every script was run from the package root into a new `results/<script>/` directory, with numpy 2.3.5, scipy 1.17.0 and matplotlib 3.10.8. Each was then rerun with `python -O` into a scratch directory. The two sets of outputs agree apart from timings, platform strings and the recorded optimisation flag; the check, failure and status fields are identical. This confirms that no check depends on `assert`.

Against the 1 November 2026 recorded `numerical_checks.json`, every risk, derivative and quadrature value agrees to a relative difference of at most 7.5 × 10⁻¹¹. Only the optimiser's minimising angles differ, by at most 3 × 10⁻⁷ relative.

The testbed numbers equal those of the development run made with numpy 2.4.6 and scipy 1.17.1, to 1e-12 relative.

## What is not verified

- No proof is machine-checked.
- Floating-point values are not interval enclosures.
- The O(ν⁻²) remainder is a sketch only.
- Conjecture 5.4 and Open problem 8.2 are supported only numerically.
- The Furrer–Bengtsson (2007) full text could not be retrieved during this revision. It is cited only for what its title and record state; its expansions were not compared with Theorem 6.1.
