\name{score_K}
\alias{score_K}

\title{Cross-validation for MFMR}
\description{Fits cross-validated likelihood and entropy for MFMR to choose K.}
\usage{
score_K( Yb, Yq, G, X=NULL, K, n.folds=10, ... )
}

\arguments{
  \item{Yb}{NxB matrix of binary traits}
  \item{Yq}{NxP matrix of quantitative traits}
  \item{G}{NxS matrix of heterogeneous covariates}
  \item{X}{NxQ matrix of homogeneous covariates}
  \item{K}{Number of subtypes, >=1}
  \item{n.folds}{Number of folds to split data into}
  \item{...}{Arguments passsed to mfrm_em_alg, the workhorse function}
}
\value{An n.folds x 2 matrix with the out-of-sample log-likelihood (column 1) and entropy (column 2) for each fold.}

\author{Andy Dahl}
\details{We recommend scoring various values of K, starting from 1, and to only study K less than or equal to the cross-validated likelihood maximizer, though other metrics could be used.}
\examples{
N <- 2000
P <- 27
B <- 3
z <- rep( 1:2, each=N/2 )

G <- matrix( rnorm( N*3 ), N, 3 )
beta1 <- matrix( rnorm(3*P), 3, P )
beta2 <- matrix( rnorm(3*P), 3, P )
Yq  <- matrix( NA, N, P )
Yq[z==1,] <- G[z==1,] \%*\% beta1 + matrix( rnorm( N/2*P ), N/2, P )
Yq[z==2,] <- G[z==1,] \%*\% beta2 + matrix( rnorm( N/2*P ), N/2, P )
Yq <- scale(Yq)

Yb  <- matrix( NA, N, B )
Yb[z==1,] <- matrix( rbinom( N/2*B, 1, .25 ), N/2, B )
Yb[z==2,] <- matrix( rbinom( N/2*B, 1, .75 ), N/2, B )

# average out-of-sample likelihood for K=1,2,3
mean( score_K( Yb, Yq, G, K=1         )[,1] )
mean( score_K( Yb, Yq, G, K=2, nrun=1 )[,1] )
mean( score_K( Yb, Yq, G, K=3, nrun=1 )[,1] )
}
