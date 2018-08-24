\name{mfmr}
\alias{mfmr}

\title{Multitrait Finite Mixture of Regressions}
\description{The MFMR EM algorithm to cluster samples into subtypes using multi-trait data and large-effect covariates.}

\usage{
mfmr(
  Yb, Yq, G, X=NULL, K,
  nrun=10, init_sd=1e-2, seeds=1:nrun+898,
  mc.cores=1,
  trace=F,
  ...
)
}

\arguments{
  \item{Yb}{NxB matrix of binary traits}
  \item{Yq}{NxP matrix of quantitative traits}
  \item{G}{NxS matrix of heterogeneous covariates}
  \item{X}{NxQ matrix of homogeneous covariates}
  \item{K}{Number of subtypes}
  \item{nrun}{Number of random initializations to run}
  \item{init_sd}{Standard deviation for the randomly cluster means}
  \item{seeds}{Seeds for the random restarts}
  \item{mc.cores}{Number of cores to use in mclapply}
  \item{trace}{Whether to print updates at each EM iteration}
  \item{...}{Arguments passsed to mfrm_em_alg, the workhorse function}
}
\value{
  \item{pmat}{NxK matrix of MFMR posterior cluster responsibilities}
  \item{out}{The MFMR EM algorithm output for the initialization giving the highest likliheood}
  \item{ll}{The MFMR log-likelihood for the initialization giving the highest likliheood}
  \item{lls}{The MFMR log-likelihoods for each random initialization}
}

\author{Andy Dahl}
\note{Currently, both Yb and Yq must be provided.}
\details{The covariates (X and G) should contain variables that may confound the primary subtype structure, e.g. sex or population structure. The traits (Yb and Yq) should contain variables that differ in distribution between subtypes. G should usually contain a constant column for main subtype effects. X is mostly used for internal purposes (in droptest) and should usually be omitted.}
\examples{
N <- 5000
P <- 27
B <- 3
z <- rep( 1:2, each=N/2 ) # subtypes

G <- matrix( rnorm( N*3 ), N, 3 ) # covariates
beta1 <- matrix( rnorm(3*P), 3, P )
beta2 <- matrix( rnorm(3*P), 3, P )

# quantitative traits
Yq  <- matrix( NA, N, P )
Yq[z==1,] <- G[z==1,] \%*\% beta1 + matrix( rnorm( N/2*P ), N/2, P )
Yq[z==2,] <- G[z==1,] \%*\% beta2 + matrix( rnorm( N/2*P ), N/2, P )
Yq <- scale(Yq)

# binary traits
Yb  <- matrix( NA, N, B )
Yb[z==1,] <- matrix( rbinom( N/2*B, 1, .25 ), N/2, B )
Yb[z==2,] <- matrix( rbinom( N/2*B, 1, .75 ), N/2, B )

out <- mfmr( Yb, Yq, G, K=2 )

# accuracy measures:
round( cor( out$pmat[,1], z )^2 , 2 )
round( cor( as.numeric( beta1 ), as.numeric( out$out$beta[1,,-(1:B)] ) ), 2 )
round( cor( as.numeric( beta1 ), as.numeric( out$out$beta[2,,-(1:B)] ) ), 2 )
}
