# Summary of the internal review and the revision

The 1 November 2026 version of this paper went through one internal review round on 3 November 2026. Five reviewer roles took part: an editor, a methodology reviewer, a domain reviewer, a cross-disciplinary reviewer and a devil's advocate. A separate adversarial check was also run. All reviews were AI-assisted. This was **not** external peer review.

This file lists the reviewers' objections and how the revision answers each one. It contains no scores or grades.

## Overall verdict on correctness

No reviewer found an error in a stated result. Each load-bearing formula and constant was independently recomputed:
- Lemma 1, Theorem 2, the 460/773 certificate, H₀ = 23148, t = 1/972216 and the −1/20416536 gap;
- the Theorem 5 expansion, Theorem 6 and the two-dimensional closed form.

Reviewers also ran Monte Carlo simulations of the actual ensemble procedure. These confirmed the derivative, the curvature and the table entries.

All reviewers agreed that reading Snyder and Hakim's conjecture with physical mean squared error is faithful to their experiments.

## Objections and resolutions

| # | Objection | Resolution in this version |
|---|---|---|
| 1 | The failure looks like a small effect in a 2 × 2 toy problem. The criterion was never applied to Snyder and Hakim's own test case. | New Section 10 reconstructs their 200-variable testbed. For N = 8 to 64 the best basis found removes 42–54% and symmetric (ZCA) whitening 41–54% of canonical excess error, roughly equivalent to doubling the ensemble. Sensitivity, prior-normalised loss and estimated bases are also reported. It is labelled a reconstruction, not a replication. |
| 2 | The canonical side was not allowed to inflate the observed variance, so the baseline could be called a strawman. | Lemma 3.5 covers every symmetric taper with any T₁₁ = τ > 0, positive semidefinite or not. Lemma 3.6 identifies the best τ and when it exists. Theorem 3.7 proves that a rotation still wins at every N when it applies the same rescaling in the rotated basis; Remark 3.8 shows that an unrescaled rotation need not. |
| 3 | The headline depends on the loss, and this was buried. | The loss split now leads the paper, in the subtitle, the abstract and Section 1.3. Under physical loss the conjecture fails at every N. Under prior-normalised loss it holds at leading order (Theorem 7.1). Its finite-N status is reduced to a scalar inequality (Proposition 8.1, Open problem 8.2). |
| 4 | The fidelity discussion did not engage with Snyder and Hakim's own evidence. | Section 1.2 quotes their definition and restatement, with section numbers. It notes that their Eq. 15 basis is a rotation in our family. It states their convention of localising both covariances, and our convention for several observations. It cites their remark that the transformation is an oracle target. |
| 5 | The prose read like an internal record and lacked intuition. | The author line is now "Anonymous". The search and agent notes are gone from the text. The paper now has a background section for non-specialists and a plain-words mechanism before Theorem 3.3. It also has a results-at-a-glance table, a roadmap, a notation table and a short reproducibility and assurance statement. |
| 6 | Literature coverage was thin. | Section 1.5 relates the work to optimal localisation (Flowerdew; Ménétrier et al.; Periáñez et al.; Anderson; Morzfeld and Hodyss; Atkinson et al.). It also covers gain bias and sampling-error expansions (van Leeuwen; Hodyss and Morzfeld; Sacher and Bartello; Furrer and Bengtsson) and transformed or spectral localisation (Buehner and Charron; Kepert; Lorenc; Nadeem and Potthast; Shlyaeva et al.). On the statistics side it cites Kessy et al. and Ledoit and Wolf. It also cites Snyder (2025) and Bourin. All 28 references (Fujii 2005 added in the verification round) are generated from Crossref records. Lemma 3.1 is credited as the zero-correlation case of Flowerdew's optimal localisation. |
| 7 | The remainder argument in Theorem 5 was only sketched. | Lemmas 6.2–6.3 now write it out. They give explicit gain-derivative bounds, the split of second-order terms, the moment and tail bounds, and a short proof of the Schur-multiplier norm bound. The O(ν^{−3/2}) rate is not claimed sharp; O(ν^{−2}) is given only as a sketch. |
| 8 | With several observations the comparison class was narrower than in the scalar case, and the paper did not say so. | Several additions address this. Lemma 7.4 extends the leading-order taper comparison to several observations, and Corollary 7.6 gives the large-N failure against every unit-diagonal PSD canonical taper, with its exact scope in Remark 7.7. At finite N, Proposition 5.1 gives the exact first variation for any observation matrix. Corollary 5.2 shows failure whenever the loss couples an observed and an unobserved direction, and Corollary 5.3 shows stationarity when [W, C] = 0. The two-observed-direction sign property is stated as Conjecture 5.4, with numerical support. |
| 9 | Finite-N behaviour under prior-normalised loss was left open. | For a scalar observation, Proposition 8.1 proves a strict local minimum modulo completion if ζ > 0 and no local minimum if ζ < 0; ζ = 0 is unresolved. The inequality ζ > 0 holds in every computed case (minimum ratio 1.0020) but is posed as an open problem, not claimed. |
| 10 | The aligned example was proved only at N = 3. | Remark 4.2 tabulates the curvature from N = 2 to 1000. It is positive at N = 2 with a known mean and negative otherwise. The decay is about −8/N² and −12/N², which corrects a rougher figure in one review. N = 3 remains the proved case. |
| 11 | Script hygiene: results were written in place, checks used `assert`, the status was hard-coded, a temporary directory was left behind, and nothing simulated the real procedure. | All four scripts now use `--out` to write to a new directory and refuse to overwrite. They collect failures explicitly, compute their status, exit with code 1 on failure, record their own SHA-256 and leave no temporary files. `theory_checks.py` adds a seeded Monte Carlo of the actual ensemble procedure, covering Lemma 3.1, Theorem 3.3, Theorem 3.7, Proposition 4.1 and two table entries, plus negative controls. `-O` reruns give the same results. |
| 12 | Notation clashed: a and aᵢ; scalar A, B and matrix A, B; c; K. | The eigenvalues are now λᵢ, μᵢ. The curvature constants are κ, ω and ζ. The integral is η. The two-dimensional coefficients are F̄, f_c and f_s. Remainder constants are ϰ. |
| 13 | The fixed-dimension asymptotics do not match the operational regime m ≫ N. | Remark 6.5 presents F as a basis-ranking heuristic outside fixed dimension. The testbed checks how accurate the leading order is. |
| 14 | The homogeneous (Fourier) case that Snyder and Hakim highlight was not examined. | Corollary 10.1 shows that under physical loss the Fourier basis is not leading-order optimal whenever two modes have prior variance below the observation-error variance. This is checked numerically, with a control. |
| 15 | The bases are oracles; would an estimated basis keep the gain? | Section 10 tests pilot-estimated bases. The gain survives only with accurate estimation of small prior variances, and heavy regularisation reverses it. Estimated canonical bases are robust. Data-driven estimation is listed as an open problem. |
| 16 | Several minor points. | All are addressed. Corollary 3.4 now proves that a canonical square root exists and does not depend on the completion. The orientation of the contrast columns is specified. The observation-noise variance is defined in Section 2. Section 9 explains why the tabulated rotation (−14.25°) differs from θ* (−15.86°) and reports the per-N optimal angles. Table 3 adds a Kalman-filter-normalised column. A Schur–Horn reading of the criterion is given. Appendix B notes that the two quadratures share a model, which is why the Monte Carlo is needed. |

## Points not adopted, or adopted only in part

- **Tuned distance-based localisation and cycling.** Neither was added to the testbed. Both are listed as limits and open problems.
- **References without a verified record or a real role in the text.** These were not cited. Examples are Bhatia's *Matrix Analysis*, the Mandel and Kwiatkowski convergence papers and some localisation papers suggested in review.
- **Optional schematic figure.** It was not added. The mechanism is explained in words instead (Section 3.2).
- **Furrer and Bengtsson (2007).** The accepted manuscript could not be retrieved during the revision, because the repository was unavailable. The paper cites it only for its subject and says explicitly that its expansions were not compared term by term with Theorem 6.1.
- **Two references the draft cited from memory.** The draft theory fragment cited Karlin's *Total Positivity* and a specific theorem number in Horn and Johnson; both have been dropped. They are replaced by self-contained arguments: a monotone-likelihood-ratio proof, a Stein-identity proof and a direct Schur-multiplier proof. Horn and Johnson is cited only generally.

## Reviewer computations reconciled

- **Testbed numbers.** Two reviewers' independent reconstructions of the testbed agree with Section 10 at the Monte Carlo level.
- **One exploratory figure.** It had used an unintended correlation scale, set by a command-line argument. Rerun with the intended scale, it agrees with the other reconstructions.

## Verification round (3 November 2026)

A proof audit, an adversarial audit and an editor's re-review of the revised version found no error in any proof. They asked for wording fixes, all made:

| Point | Change |
|---|---|
| Inflation result overstated in the abstract and Table 1 | Both now say that the rotated competitor uses the same deterministic rescaling. Remark 3.8 records that an unrescaled rotation does not beat the best inflated canonical taper (N = 3: 2.7453 vs 2.7408). |
| Several observations | The abstract and Table 1 now say that the finite-N result compares only with canonical complete localisation. |
| Prior-normalised loss | Scalar observation only. Strict local minimum modulo completion if ζ > 0, not a local minimum if ζ < 0, ζ = 0 unresolved. The "iff" is removed. The leading-order comparison class is stated, and Remark 7.5 extends it to unit-diagonal tapers in any whitened basis (Wishart moment plus AM–GM). |
| Discussion implied that physical loss alone condemns canonical coordinates | Remark 7.2 gives a counterexample (P = [[3,1],[1,3]], H = I, R = I), checked by hand and numerically. The loss-concentration mechanism is restricted to the testbed, and practical advice is limited to a known prior covariance. |
| "Optimum" used for an optimiser's best value | Replaced by "best value found" throughout, including the figure legend (script relabelled and rerun; numbers unchanged). |
| Testbed details | Table 1 states each method's range. The omitted estimated-basis setting (N_p = 200, δ = 10⁻⁴) is reported. The 1.136 vs 1.137 difference is explained. Remark 6.5 is corrected. |
| Figure 1 label 7.64% vs 7.73% | The caption explains both. |
| Remark 7.7 claim that the optimal factor tends to one | Removed. |
| Mostajeran–Sepulchre lemma number | The published version could not be opened, so the number is dropped. The inequality is credited to Bourin and Fujii, as those authors do; Fujii (2005) is added from its Crossref record. |
| Abstract hard for general readers | A plain-language first paragraph was added, and the matrix criterion is now described in words. |
