# When canonical covariance localisation fails

**Exact finite-ensemble counterexamples, a sharp leading-order criterion, and the Snyder–Hakim test case**

Author: Anonymous · Date: 3 November 2026 (version 0.2.0-candidate: 4 November 2026) · Text CC0 1.0 · Code MIT (see `LICENSE.md`)

Version 0.2.0-candidate (unrefereed research candidate) · DOI [10.5281/zenodo.23941881](https://doi.org/10.5281/zenodo.23941881) · Supersedes 0.1.0-candidate ([10.5281/zenodo.23926884](https://doi.org/10.5281/zenodo.23926884)); see `CORRECTIONS.md` · Concept DOI (all versions) [10.5281/zenodo.23926883](https://doi.org/10.5281/zenodo.23926883) · Release page <https://evidencepress.org/releases/canonical-covariance-localisation-failure/>

## In plain English

Weather and ocean forecasts are corrected with new observations using an *ensemble Kalman filter*. The filter estimates forecast uncertainty from a small ensemble of forecasts, so the estimate is noisy. *Localisation* removes the noisiest parts of that estimate. How well this works depends on the variables in which it is done.

Snyder and Hakim (2022) proposed a special set of variables, *canonical coordinates*. In these variables the exact update splits into independent one-dimensional pieces. They conjectured that this choice is the best one to make before localising.

This paper shows that the answer depends on how errors are scored.

- **Errors scored in the ordinary physical variables.** Snyder and Hakim used this score in their experiments. Here the conjecture is false for every ensemble size. A slightly rotated set of variables does strictly better than every canonical localisation that keeps the observed variance as estimated. If the canonical approach may also rescale that variance, a rotation still wins provided it applies the same rescaling. Physical scoring does not always rule canonical coordinates out: the paper's criterion (Theorem 7.1) decides when they are best, and gives a small example where they are.
- **Errors scored relative to the forecast uncertainty (Mahalanobis error).** Here the conjecture is true for large ensembles, within the class of localisations studied. For small ensembles and a single observation, every case we computed supports it, but it is not proved.
- **Snyder and Hakim's own 200-variable test case,** rebuilt from their paper and scored their way. Ordinary symmetric whitening removes about half of the canonical method's excess error at small ensemble sizes. That gain is worth roughly doubling the ensemble. The gain is fragile, however, when the prior covariance itself has to be estimated.

The proofs are written out in full but are not machine-verified. Numerical results are reported as illustrations, kept separate from exact-arithmetic certificates. The paper has had internal, AI-assisted review only. It has not been externally peer reviewed.

## Files

| Path | Contents |
|---|---|
| `canonical_localisation.pdf` | The paper (31 pages). |
| `manuscript.tex`, `refs.bib` | LaTeX source and bibliography. |
| `tools/make_refs.py` | Builds `refs.bib` (28 entries) from saved Crossref records, copying every title, author, year, venue and DOI. Usage: `python3 -I tools/make_refs.py refs.bib RECORDS_DIR [MORE_DIRS]`. |
| `figures/` | `risk_rotation.pdf/png` (Figure 1) and `sh_testbed_mse.pdf/png` (Figure 2), copied from `results/`. |
| `scripts/exact_certificates.py` | Exact rational arithmetic for the constants used in the proofs. Python standard library only. |
| `scripts/numerical_checks.py` | Deterministic quadrature for the two-dimensional example (Table 3, Figure 1). |
| `scripts/theory_checks.py` | Checks for the inflation, vector-observation, prior-normalised-loss and aligned-example results. Includes a seeded Monte Carlo of the actual ensemble procedure (4 × 10⁶ ensembles) and negative controls. |
| `scripts/sh_testbed.py` | Reconstruction of the Snyder–Hakim one-dimensional testbed (Section 10, Tables 4–6, Figure 2). |
| `results/<script>/` | Fresh outputs of each script, each written to a new directory. |
| `VERIFICATION.md` | Each claim, its proof location, the check that supports it and the recorded result. |
| `REVIEW_SUMMARY.md` | Objections raised in internal review and how each was resolved. |
| `development_notes/` | Pre-revision (1 November 2026) audit notes, kept for the record. |
| `AI_INDEX.md` | Machine-oriented map of claims, evidence, scope and limits. |
| `RESEARCH_GATES.json` | Dated record of the prior-art scrutiny, contribution map and limitation triage. |
| `CITATION.cff`, `.zenodo.json` | Citation and archive metadata. |
| `requirements.txt`, `.github/workflows/replay.yml` | Pinned packages and the CI replay (manifest check; all four scripts in normal and `-O` modes). |
| `CORRECTIONS.md` | Correction history: what changed from 0.1.0-candidate and why. |
| `MANIFEST.sha256` | SHA-256 hashes of every other file in the package. |

## How to reproduce

Requirements:
- [uv](https://docs.astral.sh/uv/), which fetches Python and the pinned packages;
- a TeX distribution with `pdflatex` and `bibtex`, for the PDF only.

Run every command from this directory.

Each script writes into the directory given by `--out`. It **refuses to overwrite** an existing path. Use a new directory name, and never point it at the recorded `results/` subdirectories.

```bash
RUN="uv run --with numpy==2.3.5 --with scipy==1.17.0 --with matplotlib==3.10.8 --with mpmath python3 -I -B"
$RUN scripts/exact_certificates.py --out rerun/exact_certificates   # < 1 s
$RUN scripts/numerical_checks.py   --out rerun/numerical_checks     # about 20 s
$RUN scripts/theory_checks.py      --out rerun/theory_checks        # about 40 s
$RUN scripts/sh_testbed.py         --out rerun/sh_testbed           # about 90 s
```

The commands above are for bash. In zsh, write each command out in full instead of using `$RUN`.

Each script behaves the same way when it runs:
- it is single-threaded;
- it needs less than 1 GB of memory;
- it records its own SHA-256 hash, its Python and library versions, and every individual check;
- it prints `status=PASS` or `"status": "PASS"`;
- it exits with code 1 if any check fails.

The scripts contain no `assert` statements, so `python -O` does not switch the checks off. The recorded runs and `-O` reruns gave identical results apart from timings, platform strings and the recorded `-O` flag.

`exact_certificates.py` needs only the standard library. You can also run it with a plain `python3 -I -B scripts/exact_certificates.py --out DIR`.

To rebuild the PDF, run the following and then copy `manuscript.pdf` to `canonical_localisation.pdf`. With the supplied files the build gives no warnings and no overfull boxes.

```bash
pdflatex manuscript && bibtex manuscript && pdflatex manuscript && pdflatex manuscript
```

To check the integrity of the package before you regenerate anything:

```bash
shasum -a 256 -c MANIFEST.sha256
```

You can use `sha256sum -c` on Linux instead. A rebuilt PDF will have a different hash because of its embedded metadata.

## Recorded results (this package)

| Script | Checks | Status |
|---|---|---|
| `exact_certificates.py` | 37 (including 4 negative controls) | PASS |
| `numerical_checks.py` | 82 (including 5 negative controls) | PASS |
| `theory_checks.py` | 296 (including negative controls) | PASS |
| `sh_testbed.py` | 39 | PASS |

`numerical_checks.py` was compared with the 1 November 2026 recorded output. All risks, derivatives and quadrature comparisons agree to a relative difference of at most 7.5 × 10⁻¹¹. The only larger differences are the optimiser's minimising angles, which differ by up to 3 × 10⁻⁷ in relative terms because the objective is flat there; the paper does not use those angles beyond one decimal place.

## Scope

What the paper covers:
- one analysis step with an independent Gaussian ensemble;
- an oracle prior covariance, used to build every basis;
- deterministic tapers;
- squared-error loss, either physical or prior-normalised.

What it does not cover:
- cycling, model error, sample-adaptive localisation, shrinkage or hybrid covariances;
- tuned distance-based localisation;
- operational data.

Limits on the results:
- **Several observations.** The finite-ensemble result compares only against rotated complete localisations. The comparison with every canonical taper holds only for large ensembles.
- **Section 10.** The testbed is a reconstruction from Snyder and Hakim's text, not a replication of their figures.
- **Open questions.** Two statements remain open, with numerical support:
  - the sign conjecture for two observed directions (Conjecture 5.4);
  - the prior-normalised inequality ζ > 0 (Open problem 8.2).
- **Literature search.** The search was bounded and is dated 3 November 2026. It does not establish priority.

## How to cite

Anonymous (2026). *When canonical covariance localisation fails: Exact finite-ensemble counterexamples, a sharp leading-order criterion, and the Snyder–Hakim test case* (version 0.2.0-candidate) [unrefereed research candidate]. Zenodo. https://doi.org/10.5281/zenodo.23941881

To cite all versions, use the concept DOI https://doi.org/10.5281/zenodo.23926883. Machine-readable metadata are in `CITATION.cff`.
