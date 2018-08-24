\name{mfmr_em_alg}
\alias{mfmr_em_alg}
\alias{mfmr_em_alg_hom}

\title{ECM Algorithm for MFMR}
\description{Internal function to perform the MFMR ECM algorithm}

\usage{
mfmr_em_alg(
  Yb, Y, G, X, K, trace,
  N, Q, S, P, B,
  pvec, alpha, beta, Sigma,
  maxit=1e3, tol=rep(1e-3,2)
)
mfmr_em_alg_hom( Yb, Y, X, trace, maxit=1e3, tol=1e-3 )
}

\arguments{
  \item{Yb}{NxB matrix of binary traits}
  \item{Y}{NxP matrix of quantitative traits}
  \item{G}{NxS matrix of heterogeneous covariates}
  \item{X}{NxQ matrix of homogeneous covariates}
  \item{K}{Number of subtypes}
  \item{trace}{Whether to print updates at each ECM iteration}
  \item{N, Q, S, P, B}{Data dimensions}
  \item{pvec, alpha, beta, Sigma}{Initial parameters}
  \item{maxit}{Maximum ECM iterations before termination}
  \item{tol}{Convergence tolerance (relative changes in alpha and beta)}
}
\value{
  \item{pvec, alpha, beta, Sigma}{Final parameters}
  \item{ll}{Final log-likelihood}
  \item{niter}{Number ECM iterations}
  \item{conv}{Convergence indicator}
}

\author{Andy Dahl}
\details{These functions should rarely be called directly. mfmr_em_alg_hom is the K=1 version and reduces to multi-trait regression (with mixed binary and quantitative output variables).}
%\examples{
%  mfmr( Y, G )
%}
