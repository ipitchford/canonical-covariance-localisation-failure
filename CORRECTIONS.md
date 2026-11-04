# Corrections

## 0.2.0-candidate (4 November 2026; doi:10.5281/zenodo.23941881)

Supersedes 0.1.0-candidate (3 November 2026; doi:10.5281/zenodo.23926884). That record and its tag `v0.1.0-candidate` are unchanged.

**What was wrong.** Version 0.1.0 cited Furrer and Bengtsson (2007, *J. Multivariate Anal.* 98:227-255, doi:10.1016/j.jmva.2006.08.003) only for its subject. It said that their expansions had not been compared with ours.

**What changed.** Their paper has now been compared with this one in full. Section 1.5 of the paper now states the following.

- Their entrywise optimal taper (their Eq. 20) has the same Q²/(1+Q²) form as Flowerdew's factor and vanishes where the true covariance does. It predates Flowerdew (2015).
- For a diagonal prior with H = R = I, their second-order approximation to the propagated-analysis covariance error (their Eq. 17) is minimised with zero off-diagonal taper entries. This is a covariance-error analogue of Lemmas 3.1 and 7.4.
- Their Jensen argument that the sample gain is biased low (their Eq. 6 and Corollaries 2-3) is credited alongside van Leeuwen and Hodyss and Morzfeld as the mechanism behind Theorem 3.3. Their gain-bias inflation rule (their Eq. 26) is contrasted with Lemma 3.6.
- Their expansion (their Proposition 6, Eqs. 14-15) differs from Theorem 6.1. It is the Frobenius error of the propagated analysis covariance, in the fixed prior eigenbasis, for the untapered estimator, valid on a convergence event, with a remainder order stated but not proved. It contains no basis objective like F(U), and they never optimise over bases.

**What did not change.** No theorem, proof, script, recorded output or figure has changed. The contribution and priority assessment is unchanged; the comparison found no earlier statement of the basis question, the criterion of Theorem 7.1 or the finite-N obstructions.

Files changed: `manuscript.tex` and the PDF, `README.md`, `VERIFICATION.md`, `REVIEW_SUMMARY.md`, `RESEARCH_GATES.json`, `AI_INDEX.md`, `CITATION.cff`, `.zenodo.json`, `development_notes/novelty_audit.md` (dated addendum), this file and `MANIFEST.sha256`.
