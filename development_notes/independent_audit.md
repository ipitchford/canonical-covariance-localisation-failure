> Development note: this file describes the 1 November 2026 pre-revision version of the manuscript; theorem numbers and some statements differ from the current paper.

# Independent mathematical audit of the canonical-localisation counterexample

Audit date: 1 November 2026. This audit was derived separately from the proposed argument by reconstructing the estimation risk and differentiating it. It is an informal mathematical audit, not proof-assistant verification.

## 1. Model and source matching

The inspected primary source is Snyder, C., & Hakim, G. J. (2022), *An optimal linear transformation for data assimilation*, Journal of Advances in Modeling Earth Systems, 14, e2021MS002937, https://doi.org/10.1029/2021MS002937. Full text: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2021MS002937.

Section 3 defines optimality through the smallest averaged mean-squared analysis error, allows the localisation to depend on the chosen transformation, and places only linearity and invertibility restrictions on transforms. Sections 3.1–3.2 identify whitening followed by the observation-operator singular-value decomposition; Section 4 uses complete localisation in these coordinates. Section 7 reiterates the conjecture without proof. The counterexample targets mean-square error in the *original physical coordinates*. A transformed-coordinate unweighted loss would be a different problem.

The present theorem compares all **deterministic correlation-matrix Schur tapers** in the canonical coordinates, as well as complete localisation after a competing rotation. It does not claim dominance over arbitrary sample-adaptive taper functions, covariance inflation, or other modifications of the estimator. The source does not explicitly formalise those wider choices, so this boundary should be stated.

## 2. Reconstructed risk identity

Let z~N(0,I_p), let y=h z_1+epsilon with epsilon~N(0,1), and let the independent ensemble contain N>=2 draws. Write nu=N-1 and beta=h^2>0. Let the ensemble mean m and unbiased covariance S be computed in the usual way. Gaussian sampling makes m independent of S and Cov(m-z)=a I_p, where a=1+1/N. For the known prior mean, replace m by 0 and a by 1.

For any covariance-derived gain k, the estimate is m+k(y-h m_1). Under weighted squared error with W positive definite, averaging over z, epsilon, and m conditional on k gives

F(k)=a tr(W)-2a h e_1^T W k+(a beta+1)k^T W k.

This identity is exact. It accounts for the finite-ensemble mean and observation noise separately. It holds for N=2, despite the rank-one sample covariance, because no inversion of S occurs and the gain denominator is strictly positive.

## 3. Rotation derivative

For a two-dimensional rotation Q_theta=[[cos theta,-sin theta],[sin theta,cos theta]], rotational invariance of the Gaussian ensemble allows the two retained diagonal sample variances in the rotating frame to be represented by independent D_1,D_2~chi-square_nu/nu for every theta. This is a coupling of their *marginal law* across theta, which is legitimate when differentiating the expected risk; it does not assert that diagonal variances of one fixed ensemble remain constant while its coordinates rotate.

Put B_theta=Q_theta diag(D_1,D_2) Q_theta^T and k_theta=h B_theta e_1/(1+beta e_1^T B_theta e_1). At zero,

k_0=(h D_1/(1+beta D_1),0)^T,
k'_0=(0,h(D_1-D_2)/(1+beta D_1))^T.

Substitution into the exact risk identity yields

R'(0)=2 beta W_12 E[(D_1-a)(D_1-D_2)/(1+beta D_1)^2]
      =2 beta W_12 E[(D-1)(D-a)/(1+beta D)^2].

The factor is strictly positive for every finite nu>=1 and a>=1. To see this, write it as

E[(D-1)^2/(1+beta D)^2]-(a-1)E[(D-1)/(1+beta D)^2].

The first term is positive because D is nondegenerate. The second expectation is Cov(D,(1+beta D)^(-2))<0, because the second function is strictly decreasing. A direct proof uses independent D,D': twice the covariance equals E[(D-D')(f(D)-f(D'))], which is strictly negative. Hence the full expression is positive.

Conclusion: if W_12>0, sufficiently small negative theta strictly improves risk. The sign reverses when W_12<0. The proposed sign and all normalisations check out.

## 4. Canonical localisation optimised over all deterministic tapers

In dimension two, an admissible correlation taper has diagonal entries one and off-diagonal t, with |t|<=1. Canonical localisation gives

k=(h S_11/(1+beta S_11),h t S_12/(1+beta S_11)).

Conditional on D=S_11, S_12 is centred normal with variance D/nu. Terms linear in S_12 therefore vanish after expectation. The exact excess risk is

R_can(t)-R_can(0)
 = beta(a beta+1) W_22 t^2 E[S_12^2/(1+beta S_11)^2]
 = beta(a beta+1) W_22 t^2/nu E[D/(1+beta D)^2].

In the direct-observation counterexample, beta=1 and W_22=1, so the prefactor reduces to a+1. In general W_22>0 because W is positive definite, and this excess risk is positive whenever t!=0. Thus complete canonical localisation is indeed the best deterministic canonical correlation taper, and the rotated candidate beats that entire class.

In dimension p, let C be any deterministic correlation matrix and retain the first-column coefficients t_j=C_j1. Conditional on S_11=D, (S_21,...,S_p1) is a centred Gaussian vector with covariance (D/nu)I_(p-1). Therefore

R_can(C)-R_can(I)
 = beta(a beta+1)/nu E[D/(1+beta D)^2] sum_(j>1) W_jj C_j1^2.

The minimisers have C_j1=0 for every j>1. They need not have a diagonal nullspace block, so the minimum is not unique in dimensions above two. This qualification prevents an unnecessarily strong uniqueness claim.

## 5. Generic higher-dimensional obstruction

A rotation between canonical axes 1 and j has derivative

2 beta W_1j E[(D-1)(D-a)/(1+beta D)^2].

Thus a necessary condition for local optimality of the canonical transformation under complete localisation is W e_1 parallel to e_1. This is only a necessary condition; vanishing first derivatives do not establish optimality.

For ordinary Euclidean loss in physical variables x=Lz, W=L^T L. If H L=h e_1^T, then W e_1 parallel to e_1 is equivalent to P H^T parallel to H^T, where P=LL^T. Proof: L^T H^T=h e_1, so P H^T=h L e_1. If W e_1=lambda e_1, multiplying by L^(-T) gives L e_1=(lambda/h)H^T, hence P H^T=lambda H^T; the converse follows by reversing the calculation. Therefore canonical localisation generally fails first-order stationarity whenever a scalar observation is not aligned with an eigenvector of P.

## 6. Explicit rational finite-angle certificate

The existence argument already proves the theorem. The following additional certificate supplies a fully specified invertible transform with no numerical integration. It is intentionally conservative.

Use the direct-observation physical model L=[[1,0],[1,1]], P=[[1,1],[1,2]], H=[1,0], with beta=1 and W=L^T L=[[2,1],[1,1]]. Take either a=1 or a=1+1/N. The following bounds also apply to the earlier transposed-shear example, because its W has the same norm and W_12=1. Let nu=N-1. Gamma moments give

E(D-1)^2=2/nu,
E[D(D-1)^2]=(2/nu)(1+4/nu),
E[D^2(D-1)^2]=(2/nu)(1+14/nu+24/nu^2).

Cauchy–Schwarz applied to |D-1|/(1+D) and |D-1|(1+D) gives

E[(D-1)^2/(1+D)^2]
 >= [E(D-1)^2]^2 / E[(D-1)^2(1+D)^2]
 = nu/(2nu^2+11nu+12).

Since a>=1 only increases the derivative factor,

R'(0)>=g_nu:=2nu/(2nu^2+11nu+12)>0.

Here is a uniform second-derivative bound. Write T=D_1+D_2 and Delta=|D_1-D_2|. Since B_theta has eigenvalues D_1,D_2,

||B||<=T, ||B'||=Delta, ||B''||=2Delta,

and because 1+B_11>=1, differentiating k=B e_1/(1+B_11) gives

||k||<=T,
||k'||<=T(1+T),
||k''||<=2T(1+T)^2.

The norm of W is (3+sqrt(5))/2<3, and a<=3/2. Therefore

|R''(theta)|
 <= 18 E[T(1+T)^2]+45 E[T^2(1+T)^2]
 <= 18*66+45*488=23148=:H_0.

The last bounds follow from T~Gamma(shape=nu,scale=2/nu), so for nu>=1 the first four moments are at most 2,8,48,384. These moment bounds also dominate the differentiated integrand and justify differentiation under expectation.

Set

t_nu=g_nu/(4H_0),
Q_nu=1/(1+t_nu^2) [[1-t_nu^2,2t_nu],[-2t_nu,1-t_nu^2]].

This matrix has rational entries. It is the rotation by theta=-delta, where delta=2 arctan(t_nu). As 0<t_nu<1, t_nu<=delta<=2t_nu. Taylor's theorem then yields

R(Q_nu)-R(I)
 <= -g_nu delta+(H_0/2)delta^2
 <= -g_nu t_nu+2H_0 t_nu^2
 = -g_nu^2/(8H_0)<0.

For N=3, nu=2, this simplifies to

g_nu=2/21,
t_nu=1/972216,
R(Q_nu)-R(I)<=-1/20416536.

Those fractions were recalculated with Python's exact Fraction class; H_0 was also recalculated as 23148. This is an exact conservative upper bound on the risk difference, not an approximate Monte Carlo result.

## 7. Remaining limitations and recommendations

* State the physical loss and the estimator precisely; a loss changed with coordinates would change the claim.
* Define the taper class as deterministic positive-semidefinite correlation matrices. Do not imply coverage of arbitrary adaptive covariance corrections.
* Keep the unknown-mean and known-mean cases separate. Both are proved, with a=1+1/N and a=1 respectively.
* Spell out the rotational-invariance coupling; otherwise readers may mistakenly think the same realised covariance was held diagonal while rotating.
* The rational witness proves strict improvement, but its angle and gap are deliberately tiny because of a loose global bound. Numerical examples with useful angles can illustrate magnitude independently; they are not needed for correctness.
* The failure of stationarity is generic in the rank-one observation class, but does not solve the remaining optimisation over transforms or disprove every narrower invariant-loss conjecture.
* This audit did not independently complete a current-literature novelty search; its source inspection establishes the original conjecture, not universal priority.

## 8. Additional audit of the complete proof draft

After reconstructing the core independently, I read `research_enkf/proof_draft.md` in full, including its general large-ensemble extension. Equations (1)–(4), the explicit physical model, and the coordinate-free eigenvector condition agree with the derivations above. Equation (2) correctly includes the factor beta(a beta+1), and correctly limits uniqueness to the observed–unobserved taper coefficients.

I also checked the leading objective for arbitrary observation matrix G. With C=(I+G^T G)^(-1), B=I-C, and A=CWC, differentiating the Kalman gain gives dK(E)=C E G^T(I+GG^T)^(-1). Completing the quadratic error expression around the true gain then gives the leading covariance-estimation contribution (2/nu) sum_i (u_i^T B u_i)(u_i^T A u_i). The independent sample mean contributes tr(WC^2)/N. These factors and matrix orders are correct.

For B=diag(b1,b2), A=[[p,q],[q,r]], direct expansion gives the trigonometric coefficients claimed in the draft: c=(b1-b2)(p-r)/4, d=(b1-b2)q/2, and constant tr(B)tr(A)/2+c. For the chosen example this yields F(theta)=(5-cos(4theta)+2sin(4theta))/16, with theta*=-(1/4) arctan(2), minimum (5-sqrt(5))/16, and canonical value 1/4. Bayes risk is 2; the leading sample-mean term is 3/(2N). The stated 30.90% reduction concerns the covariance-estimation contribution, while 7.73% concerns total excess over Bayes risk asymptotically. Neither is a reduction of the entire physical risk by that percentage.

The Taylor-remainder justification in the draft uses bounded dimension and fixed G,W. Those restrictions are essential: the stated o(1/nu) does not automatically remain valid when dimension grows with N. Within the fixed-dimensional scope the moment and exponentially small tail argument is sufficient.

## 9. Global criterion for the leading objective

I independently checked the proposed stronger characterisation of canonical optimality for F(U). In a canonical basis B is diagonal, with entries b_i. The derivatives 2(b_i-b_j)A_ij vanish for every plane exactly when A and B commute. If they fail to commute, a nonzero first derivative produces a decreasing rotation.

If A and B commute, simultaneously diagonalise them, including diagonalising A within repeated-eigenvalue blocks of B. Let their paired eigenvalues be a_i,b_i. A pair with (a_i-a_j)(b_i-b_j)>0 admits a rotation that decreases F; this follows from the two-dimensional cosine coefficient. Conversely, suppose the pairs can be ordered so a_i are nondecreasing and b_i nonincreasing. For any orthogonal basis matrix U, the entrywise squares M_ji=U_ji^2 form a doubly stochastic matrix, and F(U)=a^T(MM^T)b. The matrix MM^T is also doubly stochastic. The rearrangement inequality, or its immediate convex-combination extension to doubly stochastic matrices, gives F(U)>=a^T b=tr(AB). Thus the commuting/opposite-order criterion is necessary and sufficient for global optimality of the canonical basis for this leading objective.

Degenerate eigenvalues require the phrase 'there exists a simultaneous eigenvalue ordering with opposite order'. A need not already be diagonal within equal-b_i blocks in the arbitrary canonical basis supplied by a singular-value decomposition. Such block rotations leave the canonical risk objective unchanged.

This characterisation concerns the leading large-N objective. It must not be advertised as an exact finite-N global optimisation theorem. The separate aligned finite-N counterexample in `aligned_curvature.md` shows that the obstruction survives beyond the first derivative and is not repaired by simple covariance–observation commutation.

## Audit verdict

The proposed rank-one finite-ensemble counterexample is mathematically sound under the assumptions stated above. The derivative, strict sign, and optimality within the entire deterministic canonical Schur-taper class survived independent reconstruction. The high-dimensional generic obstruction and the explicit rational finite-angle certificate strengthen the result. No fatal mathematical defect was found in the inspected claim.
