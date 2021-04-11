\name{droptest}
\alias{droptest}
\alias{rgwas}

\title{RGWAS Test for Large-Effect Covariates}
\description{Performs the MFMR test used in RGWAS for large-effect covariates.}

\usage{
droptest(
  Yb, Yq, G, X=NULL, test_inds, K,
  mc.cores=1,
  mfmrx=FALSE, pmat=NULL,
  ...
)
}

\arguments{
  \item{Yb}{NxB matrix of binary traits}
  \item{Yq}{NxP matrix of quantitative traits}
  \item{G}{NxS matrix of heterogeneous covariates}
  \item{X}{NxQ matrix of homogeneous covariates}
  \item{test_inds}{Columns of G to run the droptest on}
  \item{K}{Number of subtypes}
  \item{mc.cores}{Number of cores to use in mclapply}
  \item{mfmrx}{Whether to refit MFMR for each column of G}
  \item{pmat}{If mfmrx=TRUE, pmat can be given to skip the initial MFMR run used to find pmat}
  \item{...}{Arguments passsed to mfmr}
}
\value{
  \item{betas}{A K x ncol(G) x (ncol(Yb)+ncol(Yq))-dimensional array of effect sizes in each subtype for each covariate-trait pair. If subtype==TRUE, only one trait is tested and the array reduces to K x ncol(G).}
  \item{ses}{Array where each entry gives the standard error for the corresponding entry in betas.}
  \item{pvals}{A 3 x ncol(G) x (ncol(Yb)+ncol(Yq))-dimensional array of p-values for the hom., het. and global tests for each covariate-trait pair. If subtype==TRUE, only the heterogeneity test is performed, and only for one trait, giving an ncol(G)-length vector.}
}

\author{Andy Dahl}
%\references{}
\details{This tests columns of G in turn with the large-effect RGWAS test. For small-effect covariates, interxn_test can be used directly.}
\note{G (usually) should have a column of 1s for main subtype effects.}
\examples{
N  <- 2000
P  <- 27
B  <- 3
z  <- rep( 1:2, each=N/2 )

# heterogeneous for first 25 traits, then homogeneous, then null
G     <- cbind( 1, matrix( rnorm( N*2 ), N, 2 ) )
beta1 <- cbind( matrix( rnorm(3*(P-2)), 3, P-2 ), 1, 0 )
beta2 <- cbind( matrix( rnorm(3*(P-2)), 3, P-2 ), 1, 0 )
Yq    <- matrix( NA, N, P )
Yq[z==1,] <- G[z==1,] \%*\% beta1 + matrix( rnorm( N/2*P ), N/2, P )
Yq[z==2,] <- G[z==1,] \%*\% beta2 + matrix( rnorm( N/2*P ), N/2, P )

# scale to help mfmr initialization
Yq <- scale(Yq)

Yb <- matrix( NA, N, B )
Yb[z==1,]  <- matrix( rbinom( N/2*B, 1, .25 ), N/2, B )
Yb[z==2,]  <- matrix( rbinom( N/2*B, 1, .75 ), N/2, B )

out  <- droptest( Yb, Yq, G, test_inds=2:2, K=2 )

# mostly significant for first 9 quant. traits, not last one
round( out$pvals['Hom',1,B+1:P], 3 )
# mostly significant for first 8 quant. traits, not last two
round( out$pvals['Het',1,B+1:P], 3 )
}

\seealso{\code{\link{interxn_test}}}
