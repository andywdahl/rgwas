\name{fit_z}
\alias{fit_z}

\title{Multitrait Finite Mixture of Regressions}
\description{Returns the fitted posterior responsibilities given data and MFMR parameters.}

\usage{
fit_z( Yb, Yq, G, X, out )
}

\arguments{
  \item{Yb}{NxB matrix of binary traits}
  \item{Yq}{NxP matrix of quantitative traits}
  \item{G}{NxS matrix of heterogeneous covariates}
  \item{X}{NxQ matrix of homogeneous covariates}
  \item{out}{MFMR fit as output from mfrm_em_alg}
}
\value{
  \item{pmat}{NxK matrix of fitted cluster probabilities}
}

\details{This function can be used to predict subtypes out of sample (or in cross-validation).}

\author{Andy Dahl}
