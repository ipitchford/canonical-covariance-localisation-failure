# AI index: when canonical covariance localisation fails (v0.2.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with
[README.md](README.md) and the paper [canonical_localisation.pdf](canonical_localisation.pdf).
DOI: https://doi.org/10.5281/zenodo.23941881 (concept DOI 10.5281/zenodo.23926883). Supersedes 0.1.0-candidate (10.5281/zenodo.23926884) by a citation correction only; see [CORRECTIONS.md](CORRECTIONS.md).
Release page: https://evidencepress.org/releases/canonical-covariance-localisation-failure/

## Question

Snyder and Hakim (2022, JAMES, doi:10.1029/2021MS002937) conjectured that canonical coordinates (whiten the prior
covariance P and observation-error covariance R, then rotate by the singular value decomposition of the whitened
observation operator) give the smallest mean squared analysis error of an ensemble Kalman filter after covariance
localisation. This package answers: it depends on the loss used to score the analysis.

## Normalisation

One analysis step; state x = L z with P = L L^T; ensemble of size N, i.i.d. Gaussian, independent of the truth; nu = N - 1;
oracle (true) P used to build every basis; deterministic symmetric Schur-product tapers; loss (x_hat - x)^T Omega (x_hat - x)
with whitened weight W = L^T Omega L. Physical loss: Omega = I. Prior-normalised (Mahalanobis) loss: Omega = P^{-1}, W = I.
Complete localisation keeps only sample variances in the chosen orthonormal basis U = [u_1 ... u_m] of whitened space.
C = (I + G^T G)^{-1}, B = I - C, A = C W C. Section 2 and Table 2 of [manuscript.tex](manuscript.tex) fix the notation.

## Claims and where they are proved

Numbering follows [canonical_localisation.pdf](canonical_localisation.pdf). Full ledger: [VERIFICATION.md](VERIFICATION.md).

| ID | Claim | Location | Evidence type |
|---|---|---|---|
| C1 | Lemma 3.1: complete localisation is the best canonical taper with unit observed diagonal (zero-correlation case of Flowerdew 2015) | Section 3.1 | proof; Monte Carlo check |
| C2 | Theorem 3.3 / Corollary 3.4: scalar observation, physical loss, any N >= 2: if H^T is not an eigenvector of P, a small rotation plus complete localisation beats every canonical taper with unit observed diagonal | Section 3.2 | proof; derivative checks (quadrature and Monte Carlo of the ensemble procedure) |
| C3 | Lemmas 3.5-3.6, Theorem 3.7: with the same deterministic rescaling diag(tau_0,1,...,1), a rotation beats every canonical taper with any observed-variance factor tau > 0; an unrescaled rotation need not (Remark 3.8) | Section 3.3 | proof; numerical checks |
| C4 | Two-dimensional example: explicit rational rotation, risk gap <= -1/20416536 at N = 3 | Section 4.1, Appendix A | proof + exact rational certificate |
| C5 | Proposition 4.1: aligned example, canonical basis a strict local maximum at N = 3 (bounds -56/773, -508/773) | Section 4.2 | proof + exact; other N numerical (Remark 4.2) |
| C6 | Proposition 5.1, Corollaries 5.2-5.3: several observations, any N: not stationary if the loss couples an observed and an unobserved canonical direction; stationary if [W, C] = 0 (comparison with rotated complete localisations only) | Section 5 | proof; finite-difference checks |
| C7 | Conjecture 5.4: sign property for two observed directions | Section 5 | CONJECTURE; numerical only (660 cases) |
| C8 | Theorem 6.1 (Lemmas 6.2-6.3): uniform large-ensemble expansion, leading basis-dependent term 2F(U)/nu with F(U) = sum_i (u_i^T B u_i)(u_i^T A u_i), remainder O(nu^{-3/2}) | Section 6 | proof; testbed agreement at N = 64 |
| C9 | Theorem 7.1: canonical coordinates minimise F iff [A, B] = 0 and the eigenvalues are oppositely ordered; uses the classical antimonotone trace inequality (Bourin 1999; Fujii 2005), not claimed | Section 7 | proof |
| C10 | Lemma 7.4, Corollary 7.6: if the criterion fails, canonical tapers are beaten for all large N | Section 7.1 | proof (N_0 not explicit) |
| C11 | Prior-normalised loss: canonical coordinates leading-order optimal among complete localisations and unit-diagonal tapers in whitened bases (Theorem 7.1, Remark 7.5) | Sections 7-8 | proof |
| C12 | Proposition 8.1: scalar observation, W = I, finite N: strict local minimum modulo completion if zeta > 0, not a local minimum if zeta < 0 | Section 8 | proof |
| C13 | Open problem 8.2: zeta > 0 | Section 8 | OPEN; numerical only (minimum ratio 1.0020) |
| C14 | Section 9: closed-form rotation angle theta* = -(1/4) arctan 2 for the 2-D example; 30.90% and 7.73% figures | Section 9 | proof; numerical checks |
| C15 | Section 10: reconstruction (not replication) of the Snyder-Hakim 200-variable testbed: for N = 8 to 64 the best basis found removes 42-54% and ZCA whitening 41-54% of canonical excess error (N = 8: canonical 1.137, best found 1.063 times Kalman MSE); canonical best found under prior-normalised loss; pilot-estimated bases keep much less of the gain | Section 10, Tables 4-6, Figure 2 | NUMERICAL ILLUSTRATION; seeded Monte Carlo |
| C16 | Corollary 10.1: Fourier bases under physical loss are not leading-order optimal when two modes have prior variance below the observation-error variance | Section 10.3 | proof; numerical check with control |

## Verification artefacts

- [scripts/exact_certificates.py](scripts/exact_certificates.py): exact rational arithmetic (37 checks, 4 negative controls). Standard library only.
- [scripts/numerical_checks.py](scripts/numerical_checks.py): deterministic quadrature for the 2-D example (82 checks, 5 negative controls).
- [scripts/theory_checks.py](scripts/theory_checks.py): inflation, vector-observation, prior-normalised and aligned-example checks plus a seeded
  Monte Carlo of the actual ensemble procedure with 4 x 10^6 ensembles and negative controls (296 checks).
- [scripts/sh_testbed.py](scripts/sh_testbed.py): testbed reconstruction (39 checks), seed 20261103.
- Recorded outputs: [results/](results/). Every script records its own SHA-256, writes to a new `--out` directory, computes its status
  from explicit checks (no `assert`) and exits 1 on failure; normal and `-O` runs agree apart from timings and version strings.
- Floating-point results are illustrations and consistency checks, not interval certificates.

## Exclusions

Cycling, model error, sample-adaptive localisation, shrinkage, hybrid covariances, tuned distance-based localisation,
operational data and m >> N asymptotics are not covered. With several observations the finite-N result compares only with
rotated complete localisations; the comparison with every canonical taper is asymptotic. The O(nu^{-2}) remainder is a sketch.

## Review and provenance

Internal AI-assisted review only: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). No human peer review, independent human reproduction
or machine-checked proof. Pre-revision audit notes (1 November 2026): [development_notes/](development_notes/). Prior-art scrutiny:
[RESEARCH_GATES.json](RESEARCH_GATES.json). Bibliography generated from Crossref records: [refs.bib](refs.bib),
[tools/make_refs.py](tools/make_refs.py). Licences: [LICENSE.md](LICENSE.md).
