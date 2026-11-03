> Development note: this file describes the 1 November 2026 pre-revision version of the manuscript; theorem numbers and some statements differ from the current paper.

# A second exact counterexample: canonical localisation can be a local maximum even with aligned observations

This supplements the first-order obstruction. It uses diagonal physical covariance and a direct observation of an eigenvector. All calculations below are exact; no floating-point quadrature is used in the proof.

## General two-dimensional curvature

Use the same whitened Gaussian ensemble model and notation as the main proof, with W=diag(w1,w2), beta=h^2>0, nu=N-1, and a=1+1/N (or a=1 if the true prior mean is used). Let D,E be independent chi-square_nu/nu random variables.

For the rotation Q_theta, write u=(B_theta)_11 and v=(B_theta)_21. Then

u=D+(E-D)theta^2+O(theta^4),
v=(D-E)theta+O(theta^3).

The exact conditional risk is

w1 (a+beta u^2)/(1+beta u)^2
 +w2 [a+beta(a beta+1)v^2/(1+beta u)^2].

This expression follows from the conditional quadratic risk identity in the main proof; in particular, it includes the theta dependence of the gain denominator. The derivative of the first rational factor with respect to u is 2beta(u-a)/(1+beta u)^3. Consequently

R'(0)=0,

R''(0)=-4beta w1 A_(nu,a,beta)+2beta(a beta+1)w2 B_(nu,beta),

where

A_(nu,a,beta)=E[(D-a)(D-1)/(1+beta D)^3],
B_(nu,beta)=E[((D-1)^2+2/nu)/(1+beta D)^2].

Both A and B are positive for nu>=1, beta>0, a>=1. For A use

A=E[(D-1)^2/(1+beta D)^3]-(a-1)Cov(D,(1+beta D)^(-3)),

whose first term is positive and covariance is negative. B is manifestly positive. Thus the exact threshold is

tau_(nu,a,beta)=(a beta+1) B_(nu,beta)/(2 A_(nu,a,beta)).

If w1/w2>tau, the canonical basis is a strict local maximum of risk along this rotation family. If w1/w2<tau it is a strict local minimum along this family. Equality leaves higher-order terms to determine the local behaviour. No claim of global optimality follows from the positive-curvature case.

## A modest fixed physical model

Take

P=diag(4,1), L=diag(2,1), H=[1,0], R_noise=4,

and N=3, so nu=2. Then z=L^(-1)x has identity covariance. Dividing the observation by 2 gives y'=z1+epsilon' with noise variance one, so beta=1 and W=P=diag(4,1). The observation direction H^T is an eigenvector of P. The source's canonical transform is z, up to sign changes.

Let D be an exponential random variable of mean one and put

c=E[(1+D)^(-1)]=integral_0^infinity e^(-x)/(1+x) dx.

Integration by parts gives

E[(1+D)^(-2)]=1-c,
E[(1+D)^(-3)]=c/2.

Substituting into the curvature factors yields

A=(2a+5)c-a-3,
B=6-9c,

and therefore

R''(0)=28a+60-(50a+98)c.

## Elementary exact lower bound on c

Define the rational polynomial

p(x)=(713-413x+113x^2-13x^3+x^4/2)/773.

Since the square is nonnegative,

0 <= integral_0^infinity (1+x)[1/(1+x)-p(x)]^2 e^(-x) dx
   =c-2 integral_0^infinity p(x)e^(-x)dx
      +integral_0^infinity (1+x)p(x)^2 e^(-x)dx.

Using integral_0^infinity x^k e^(-x)dx=k!, direct rational arithmetic gives

integral p(x)e^(-x)dx=460/773,
integral (1+x)p(x)^2e^(-x)dx=460/773.

Hence c>=460/773. Equality is impossible because 1/(1+x) is not this polynomial almost everywhere, but a non-strict lower bound already suffices.

For the known prior mean, a=1, giving

R''(0)=88-148c <= 88-148*(460/773)=-56/773<0.

For the ensemble-estimated mean, a=4/3, giving

R''(0)=292/3-(494/3)c <= -508/773<0.

Thus canonical complete localisation is a strict local maximum in both mean conventions. Every sufficiently small nonzero positive or negative rotation has lower risk. By the canonical-taper theorem in the main proof, such a rotation also beats every deterministic canonical correlation taper.

## Meaning and limitations

This example strengthens the failure result: physical loss need not contain cross-coordinate terms, and the observation can be an eigenvector of the physical prior covariance. Here rotation changes the way independent variance estimates are pooled. The improvement therefore survives a natural proposed repair of the conjecture based solely on covariance–observation alignment.

The exact second-derivative threshold is valid for all finite N, a>=1, beta>0. The explicit ratio-four example above is certified at N=3. Numerical investigation suggests other small uniform ratios may work for every N, but no such universal small-ratio theorem is claimed here.

The square-polynomial lower bound and both curvature bounds were recalculated with Python Fraction arithmetic on 1 November 2026. The verifier needs only the standard library factorial and Fraction implementations. No machine proof of the surrounding unrestricted theorem was performed.
