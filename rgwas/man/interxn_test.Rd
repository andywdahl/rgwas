\name{interxn_test}
\alias{interxn_test}

\title{Covariate interaction test}
\description{Tests the interaction of g and pmat either with ordinary linear regression (if bin=FALSE) or logistic regression (if bin=TRUE).}

\usage{
interxn_test( X, y, g, pmat, bin )
}

\arguments{
  \item{X}{NxQ matrix of background covariates}
  \item{y}{Trait vector with length N}
  \item{g}{Focal covariate with length N}
  \item{pmat}{Nx(K-1) matrix of environmental covariates (e.g. subytpe proportions)}
  \item{bin}{Whether y is binary}
}
\value{
  \item{pvals}{Vector of homogeneous, heterogeneous, and global tests for g}
  \item{fit}{Output from the glm with interaction effects between g and pmat}
}

\author{Andy Dahl}
\note{The environments, pmat, should be full rank. If pmat corresponds to K subtype proportions, the columns should contatin (any) K-1 of the proportions, meaning pmat has non-negative entries and row sums in [0,1].}
\details{Fits the homogeneous, heterogeneous, and global tests for the effect of g on y, conditional on X. Heterogeneity is defined by the interaction between g and pmat.}
\examples{

N <- 5000
P <- 27
B <- 3
z <- rep( 1:2, each=N/2 ) # subtypes

G <- matrix( rnorm( N*3 ), N, 3 ) # large-effect covariates
g <- scale(rbinom(N,2,.5)) 	      # small-effect covariates
beta1 <- matrix( rnorm(3*P), 3, P )
beta2 <- matrix( rnorm(3*P), 3, P )

# quantitative traits
Yq  <- matrix( NA, N, P )
Yq[z==1,] <- G[z==1,] \%*\% beta1 + matrix( rnorm( N/2*P ), N/2, P )
Yq[z==2,] <- G[z==1,] \%*\% beta2 + matrix( rnorm( N/2*P ), N/2, P )

# add small g effect to first two quantitative traits
Yq <- scale(Yq)
Yq[z==1,1]<- Yq[z==1,1] + g[z==1] * rnorm( 1, sd=1e-1*sqrt(2) )
Yq[    ,2]<- Yq[    ,2] + g       * rnorm( 1, sd=1e-1 )
Yq <- scale(Yq)

# binary traits
Yb  <- matrix( NA, N, B )
Yb[z==1,] <- matrix( rbinom( N/2*B, 1, .25 ), N/2, B )
Yb[z==2,] <- matrix( rbinom( N/2*B, 1, .75 ), N/2, B )

out <- mfmr( Yb, Yq, G, K=2, nrun=2 )
round( cor( out$pmat[,1], z )^2 , 2 ) # MFMR estimates z well

# MFMRX test for small-effect g
# Khatri-Rao product allows G to have heterogeneoues effects
X <- t(sapply( 1:nrow(G), function(i) G[i,,drop=FALSE] \%x\% out$pmat[i,] ))

pmat1 <- out$pmat[,-1,drop=FALSE] # full-rank version of pmat
X <- cbind( X, pmat1 )

interxn_test( X=X, y=Yq[,1], g=g, pmat=pmat1, bin=FALSE )$pvals # truly het
interxn_test( X=X, y=Yq[,2], g=g, pmat=pmat1, bin=FALSE )$pvals # truly hom
interxn_test( X=X, y=Yq[,3], g=g, pmat=pmat1, bin=FALSE )$pvals # truly null

}
%\keyword{algebra}
\seealso{\code{\link{droptest}}}
