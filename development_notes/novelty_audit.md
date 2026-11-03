> Development note: this file describes the 1 November 2026 pre-revision version of the manuscript; theorem numbers and some statements differ from the current paper.

# Independent novelty audit: canonical transformation and localisation risk

Audit date: 1 November 2026. This records a bounded literature search and inspected primary sources. It does not establish historical priority.

## Proposed contribution being checked

A finite-sample counterexample to the universal mean-squared-analysis-error optimality conjectured by Snyder and Hakim (2022) for canonical whitening/SVD coordinates followed by covariance localisation. The main witness uses P=[[1,1],[1,2]], H=[1,0], R=1, and physical Euclidean squared-error loss. In whitened coordinates x=Lz with L=[[1,0],[1,1]], the observation is y=z1+epsilon and the loss matrix is W=L^T L=[[2,1],[1,1]]. Rotating whitened coordinates before retaining only sample covariance diagonals gives a nonzero exact risk derivative at the canonical basis for every finite sample covariance degrees of freedom. This replaces an earlier shear witness; the mathematical derivative is unchanged because W12=1 in both models, but risk levels and optimal angles differ.

The strengthened claim also covers every state dimension m>=2 and every scalar observation: if H^T is not an eigenvector of P, canonical complete localisation is nonstationary for every N>=2, and an improving rotation beats every deterministic canonical correlation-matrix Schur taper. A second, separately verified example P=diag(4,1), H=[1,0], R=4, N=3 has an aligned observation but canonical complete localisation is a strict local maximum. These are exact finite-sample statements, distinct from the asymptotic objective below.

The mathematical scope must state the estimator convention (known prior mean or sampled ensemble mean), the sample covariance normalisation, and the allowed canonical taper class. A deterministic Schur correlation taper conventionally has unit diagonal; allowing variance shrinkage or arbitrary covariance estimators is a different comparison class.

## Original claim and fidelity check

Snyder, C., & Hakim, G. J. (2022). An optimal linear transformation for data assimilation. *Journal of Advances in Modeling Earth Systems, 14*(6), e2021MS002937. https://doi.org/10.1029/2021MS002937

Full text inspected at https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021MS002937 .

Decision-relevant passages:

- Section 2, discussion of Figure 1: mean-squared error over all state elements of the updated ensemble mean. This supports original-state physical mean-squared loss.
- Section 3, opening paragraphs: the best transformation is the one that yields smallest average mean-squared analysis error after an ensemble Kalman filter using covariance localisation; localisation may depend on transformation. Only linearity and invertibility constrain transformations. The statement is expressly a conjecture, supported by intuition and numerical comparisons rather than a proof.
- Section 3.2: whitening the state and observation errors, then applying the singular value decomposition of the scaled observation operator, makes the population covariance and Kalman gain diagonal.
- Section 4: the proposed implementation transforms ensemble members and observations, sets all off-diagonal entries of the transformed state sample covariance to zero, then updates. Numerical comparisons involve one specific covariance/observation model, with loss in original variables.

The counterexample therefore addresses the paper's broad conjectured finite-sample risk optimality. It does not invalidate the exact population diagonalisation, canonical correlation decomposition, or its information-theoretic identities.

## Searches performed

Both general web search systems were used, with primary sources opened where material. Query families included:

1. Snyder Hakim 2022 covariance localization optimal linear transformation canonical ensemble Kalman.
2. Exact title and DOI (2021MS002937), with and without exclusions for duplicate publisher/index sites.
3. Snyder + Hakim + conjecture + localization; Snyder + Hakim + optimal transformation + counterexample.
4. Canonical observation operators + localization, including 2025/2026 date terms.
5. Covariance localization + loss function + transformation; covariance localization + risk + basis.
6. Ensemble Kalman + rotation + mean squared + optimal transformation.
7. Localisation/localization + canonical + counterexample or optimal Kalman risk.
8. The titles of adjacent optimal-localisation and sampling-error papers.
9. Snyder and Hakim (2022) citations and subsequent author conference presentations.

No exact published counterexample, proof of this conjecture, or formula-level match for the two-dimensional rotation derivative was found in these queries. Search engines do not index all literature, and an unindexed preprint or earlier different terminology remains possible.

## Inspected comparisons

### Original-source follow-up

Snyder, C. (2024). *Sampling error in the ensemble Kalman filter for small ensembles and high-dimensional states* [Conference presentation]. https://ncda-fs.web.norce.cloud/WS2024/snyder.pdf

The presentation was inspected, including ensemble conventions, original and canonical-coordinate error covariance formulas, its use of the canonical transformation, and asymptotic random-matrix approximations. It uses the coordinates as an analytical tool. It does not demonstrate global localisation-risk optimality and does not give the proposed rotated-basis counterexample.

Snyder, C. (2025). Sampling error in the ensemble Kalman filter for small ensembles and high-dimensional states. *EGU General Assembly 2025*, EGU25-14641. https://doi.org/10.5194/egusphere-egu25-14641

Full conference abstract inspected at https://meetingorganizer.copernicus.org/EGU25/EGU25-14641.html . It uses the 2022 canonical coordinates with random-matrix approximations to study sample gain, ensemble covariance and analysis error. The abstract does not claim a resolution of the localisation optimality conjecture. The full presentation linked from the EGU page was unavailable to unregistered visitors; a related 2024 presentation was inspected instead.

### Adjacent optimal-localisation theory

Morzfeld, M., & Hodyss, D. (2023). A theory for why even simple covariance localization is so useful in ensemble data assimilation. *Monthly Weather Review, 151*(3), 717–736. https://doi.org/10.1175/MWR-D-22-0255.1

Publisher search excerpts and authors' conference poster were inspected. Full publisher text returned HTTP 403 to direct retrieval in this audit. Poster: https://ncda-fs.web.norce.cloud/WS2024/morzfeld_poster.pdf . The work distinguishes optimal covariance estimation from optimal gain approximation, derives fixed-coordinate taper rules, and analyses their similarities. This is prior art for the general point that optimality depends on the criterion. It does not, in the inspected material, optimise over whitening rotations or present the proposed exact finite-N failure of the canonical basis. Do not present loss-dependence of localisation in general as newly discovered.

Atkinson, R. S., Flowerdew, J., Hughes, S., & Roulstone, I. (2024). Comparison of traditional and hybrid forms of optimal localisation for mitigation of sampling error in ensemble Kalman filters. *Tellus A: Dynamic Meteorology and Oceanography, 76*(1), 57–73. https://doi.org/10.16993/tellusa.35

Full publisher PDF inspected: https://tellusjournal.org/articles/35/files/662f688d72830.pdf . It compares fixed-coordinate single-observation gain and covariance localisation rules, including hybrid damping towards climatology, and studies failures when the single-observation derivation is applied to denser observations. It does not cite Snyder–Hakim and does not optimise basis rotations. Its existence limits novelty claims to the specific canonical-transform conjecture and exact rotation theorem.

Al Ghattas, O., & Sanz-Alonso, D. (2024). Non-asymptotic analysis of ensemble Kalman updates: Effective dimension and localization. *Information and Inference: A Journal of the IMA, 13*(1), iaad043. https://doi.org/10.1093/imaiai/iaad043

ArXiv v3 abstract and HTML sections/references inspected: https://arxiv.org/html/2208.03246v3 . The work establishes sample-complexity and covariance-estimation bounds for ensemble updates, including localisation. The inspected sections do not treat exact optimality among all canonical rotations. This is neighbouring finite-sample theory, not a collision with the proposed rotation derivative.

Hodyss, D., & Morzfeld, M. (2023). How sampling errors in covariance estimates cause bias in the Kalman gain and impact ensemble data assimilation. *Monthly Weather Review, 151*(9), 2413–2426. https://doi.org/10.1175/MWR-D-23-0029.1

Publisher search excerpts and bibliographic record inspected; direct full-text retrieval returned HTTP 403. The paper analyses gain bias from covariance sampling error. This is prior art for the mechanism of nonlinear gain bias; its full text should be inspected before any broad originality claim concerning that mechanism. No basis-rotation result appeared in inspected excerpts.

Furrer, R., & Bengtsson, T. (2007). Estimation of high-dimensional prior and posterior covariance matrices in Kalman filter variants. *Journal of Multivariate Analysis, 98*(2), 227–255. https://doi.org/10.1016/j.jmva.2006.08.003

Publisher abstract, author-hosted bibliographic record/abstract, and ZORA repository metadata were inspected on 1 November 2026. The publisher full-text request returned HTTP 403; the repository lists an accepted manuscript but requires login for download. The abstract explicitly describes mean-squared error expressions, covariance expansions and optimal taper rules. Thus broad novelty claims for sampling-error expansions or optimal tapering are inappropriate. No exact rotated-basis result was visible in the inspected material; the full paper's absence from this audit remains an access boundary, not a claim that its formulas cannot imply some part of the asymptotic reduction. Author record: https://user.math.uzh.ch/furrer/research/publications.shtml . Repository: https://doi.org/10.5167/uzh-21542 .

## Assessment by contribution unit

| Unit | Bounded assessment |
|---|---|
| Canonical population whitening/SVD diagonalisation | Known; explicitly acknowledged by Snyder–Hakim. |
| Criterion-dependence of optimal localisation | Known in fixed-coordinate literature; no novelty claim warranted. |
| Sampling error creates gain bias | Known; cite Hodyss–Morzfeld and earlier references. |
| Exact finite-N rotated-basis counterexample for P=[[1,1],[1,2]], H=[1,0], R=1 | No collision found within this search. |
| Every finite-N strict first-variation obstruction, including sampled mean | No collision found within this search; independent algebra below supports correctness, not novelty. |
| General higher-dimensional rank-one rotational nonstationarity criterion | No exact collision found in the searches recorded below. |
| Aligned exact finite-N local-maximum example | No exact collision found in the searches recorded below. |
| General leading EnKF risk objective F(U) | No exact EnKF collision found in formula and alias searches; bounded uncertainty remains. |
| Abstract antimonotone-pair sufficiency for minimising F | COLLISION with a classical matrix Chebyshev inequality; retain as a standard tool applied to the EnKF objective. |

Permitted wording: the result refutes the specified mathematical optimality claim under the stated estimator and loss conventions, with no earlier resolution found in this recorded search.

Avoid: first-ever discovery that localisation depends on basis or loss; universal superiority in practical geophysical systems; exhaustive literature search; guarantee of historical priority.

## Independent algebra sanity check

Let Q(theta) have columns (cos(theta),sin(theta)) and (-sin(theta),cos(theta)). Let A,B be independent chi-square_nu/nu variables. Rotational invariance lets the diagonals of the rotated sample covariance be represented by this fixed pair in distribution. After rotating back,

B_theta=Q(theta)diag(A,B)Q(theta)^T,
k_theta=B_theta e1/(1+e1^T B_theta e1).

Then k_0=(A/(1+A),0)^T and k'_0=(0,(A-B)/(1+A))^T.

For known prior mean, the expected physical squared-error derivative is

R'(0)=2 W12 E[(A-1)^2/(1+A)^2],

strictly positive if W12>0 and nu is finite.

For an independently sampled ensemble mean with variance cI, c=1/N and nu=N-1, the derivative becomes

R'(0)=2 W12 E[(A-1-c)(A-1)/(1+A)^2].

The first term is strictly positive. Gamma integration by parts gives

E[(A-1)/(1+A)^2]=-(4/nu)E[A/(1+A)^3]<0,

so the second term also increases positivity. This is an algebraic check, not a numerical test or formal proof-assistant verification. Rotation-sign conventions may reverse the derivative sign while preserving the nonstationarity conclusion.

## Formula-first rerun after the general leading-risk reduction

### Frozen object and normalisation

For fixed finite m,k,G,W with W positive definite, use prior-whitened state z~N(0,I), noise-whitened observations y=Gz+epsilon, independent unbiased Gaussian sample covariance with nu=N-1, and complete localisation in an orthogonal basis U=[u_i]. Put

C=(I+G^T G)^(-1), B=I-C, A=CWC.

The proposed expansion is

R_N(U)=tr(WC)+N^(-1)tr(WC^2)+2nu^(-1)F(U)+O(nu^(-3/2)),
F(U)=sum_i (u_i^T B u_i)(u_i^T A u_i),

uniform in U for fixed dimensions. The intended new unit is the EnKF risk identification and its consequences; no novelty is asserted for the abstract quadratic-product functional, Taylor expansion, orthogonal group, diagonal pinching, or rearrangement inequality.

Let P_i=u_i u_i^T and Delta_U(X)=sum_i P_i X P_i. Then the same object is

F(U)=tr(Delta_U(A)Delta_U(B))=sum_i tr(A P_i B P_i).

This is an entrywise-verified identity, not merely an analogy. In a common eigenbasis for commuting A,B, let M_ji=U_ji^2. Then M is doubly stochastic and

F(U)=a^T M M^T b,

where a,b are the jointly indexed eigenvalue vectors. The usual antimonotone-pair condition is (a_i-a_j)(b_i-b_j)<=0 for all i,j. These aliases make an independent matrix-inequality novelty claim inappropriate.

### Exact fingerprints

In the main model, A=[[1/2,1/2],[1/2,1]], B=diag(1/2,0), and

F(theta)=(5-cos(4theta)+2sin(4theta))/16.

Three exact values are F(0)=1/4, F(-pi/8)=3/16 and F(pi/4)=3/8. Here tr(A)=3/2, det(A)=1/4, its characteristic polynomial is lambda^2-(3/2)lambda+1/4; tr(B)=1/2 and det(B)=0. The commutator AB-BA=[[0,-1/4],[1/4,0]]. Sign changes and permutations of U's columns leave F invariant. Embedding this 2x2 pair in higher dimensions with A extended by identities and B by zeros leaves the rotating-block objective unchanged. No nontrivial coefficient sequence or generating function is claimed; the chi-square Laplace transform (1+2t/nu)^(-nu/2) is standard.

The three physical-prior embeddings P_2=[[1,1],[1,2]], P_3=P_2 direct-sum[1], and P_4=P_2 direct-sum I_2 have determinants one, traces 3,4,5, and characteristic polynomials (lambda^2-3lambda+1)(lambda-1)^(m-2). These elementary embedding fingerprints were used to distinguish the model from a new matrix-family claim.

### Alias map and actual extra queries (1 November 2026)

Probability/statistics: Wishart diagonal variance estimation; nonlinear gain bias; decision-weighted covariance risk. Optimisation: orthogonal basis selection; products of Rayleigh quotients; joint diagonalisation. Matrix analysis/mathematical physics: diagonal pinching, dephasing, conditional expectation, antimonotone operator pairs. Combinatorics: orthostochastic/doubly stochastic matrices and rearrangement. Special-function language: gamma or chi-square resolvent moments; no new special function is claimed.

Both web search systems were used. Additional literal query families were:

- ensemble Kalman localisation "rotation" "risk" "basis"
- "ensemble Kalman" "optimal" "basis" "mean squared" rotation
- "ensemble Kalman" "basis" "asymptotic risk"
- minimize "product of diagonal" "two symmetric matrices" orthogonal
- "sum" "Rayleigh quotients" "orthonormal" product minimization
- "u_i^T A u_i" "u_i^T B u_i"
- "localization" "CWC" Kalman
- "Kalman" "pinching"
- "covariance" "orthostochastic" "risk"
- "pinching" "oppositely ordered" matrices trace
- matrix "antimonotone pair" "trace"
- "monotone pair" "pinching" trace
- "Chebyshev" "diagonal" matrices trace inequality
- "Bourin" "monotone pairs" inequalities trace

The direct EnKF and exact-formula searches gave no matching expansion or canonical-basis classification. Many broad results concern ensemble-member rotations, which preserve ensemble mean and covariance, rather than state-coordinate rotations before nonlinear localisation; these are different operations. Search negatives remain corpus-limited.

### Verified classical antecedent: do not claim new matrix inequality

Mostajeran, C., & Sepulchre, R. (2018). Ordering positive definite matrices. *Information Geometry, 1*, 287–313. https://doi.org/10.1007/s41884-018-0003-7 . Full publisher text was inspected, particularly Lemma 3, equation (54), with its attribution to Bourin (1999) and Fujii (2005).

Lemma 3 states the trace inequality tr(A X B X)>=tr(AB X^2) for Hermitian X when A,B form an antimonotone functional-calculus pair. Set X=P_i=u_i u_i^T and sum over the rank-one projections. Since P_i^2=P_i and sum_i P_i=I, the result is exactly F(U)>=tr(AB), with equality at a joint eigenbasis. This verifies a formula-level collision for the sufficient abstract matrix condition. The elementary necessity argument using plane rotations does not justify calling the resulting abstract criterion a newly discovered matrix principle.

Bourin, J.-C. (1999). Some inequalities for norms on matrices and operators. *Linear Algebra and its Applications, 292*(1–3), 139–154. https://doi.org/10.1016/S0024-3795(99)00024-5 . Publisher bibliographic record and abstract inspected, not its full proof. The same rank-one inequality was verified directly in Bourin's *Compressions, dilations and matrix inequalities*, Chapter 4, Lemma 1.2(2), https://rgmia.org/papers/monographs/cdmi.pdf . The document states (u^T A u)(u^T B u)>=u^T ABu for an antimonotone pair. Only the relevant chapter/lemma was inspected, not the entire 93-page document. An exact publication date for this hosted monograph was not established, so the dated 2018 primary article is the safe explicit reference.

Outcome: the proposed EnKF risk identification and exact counterexamples remain **BOUNDED UNCERTAINTY / clear within the recorded search boundary**. The abstract sufficiency theorem is **COLLISION**, while its specific use to classify canonical leading risk may be a useful application/bridge. No author enquiry was made, no private databases were searched, and there is no priority certificate. A fresh stronger or structurally different formula requires another gate.
