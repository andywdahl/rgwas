\name{score_K}
\alias{score_K}
\alias{score_K_ps}

\title{Cross-validation for MFMR}
\description{Fits cross-validated likelihood and entropy for MFMR to choose K.}
\usage{
score_K( Yb, Yq, G, X=NULL, K, n.folds=10, mc.cores=1, nrun=10, ... )
score_K_ps( Yb, Yq, G, X=NULL, K, n.folds=10, mc.cores=1, nrun=10, tau=0, ... )
}

\arguments{
  \item{Yb}{NxB matrix of binary traits}
  \item{Yq}{NxP matrix of quantitative traits}
  \item{G}{NxS matrix of heterogeneous covariates}
  \item{X}{NxQ matrix of homogeneous covariates}
  \item{K}{Number of subtypes, >=1}
  \item{n.folds}{Number of folds to split data into}
  \item{tau}{Hard-threshold for cluster probabilities, needed to compute prediction strength}
  \item{nrun}{Number of random initializations to run}
  \item{mc.cores}{Number of cores to use in mclapply}
  \item{...}{Arguments passsed to mfrm_em_alg, the workhorse function}
}
\value{score_K returns an n.folds x 2 matrix with the out-of-sample log-likelihood (column 1) and entropy (column 2) for each fold.

score_K_ps returns an n.folds-length vector with estimated prediction strengths using a 50/50 train/test split.}

\author{Andy Dahl}
\details{We recommend scoring various values of K, starting from 1, and to only study K less than or equal to the maximizer of the cross-validated likelihood or prediction strength. However, other metrics could be used, and prediction strength can be rather conservative for smaller sample sizes.}
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
