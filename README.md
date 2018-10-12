# RGWAS

RGWAS infers subtypes with MFMR and then tests covariates downstream for heterogeneity across subtypes.

## Installation
R CMD INSTALL rgwas_1.0.tar.gz

## MFMR

mfmr learns latent sample subtypes using quantitative and binary traits and covariates. Ideally, the traits share subtypes, meaning groups of samples differ in distribution across many of the traits. The covariates can have different effects in different subtypes, which can be leveraged to improve subtype estimates.

## Large-effect covariate test 

droptest performs the RGWAS test for large-effect covariates. It fits MFMR for each covariate in turn, not allowing the tested covariate to have an heterogeneous effect across subtypes in MFMR (which overfits, see Supplementary Figure 1 of our RGWAS paper).

## Small-effect covariate test (MFMRX)

interxn_test performs the RGWAS test for small-effect covariates, which uses a single fit of MFMR that ignores the tested covariates. This is computationally efficient and enables GWAS-scale tests.

# Example analysis of simple simulated data

## Simulating data

RGWAS relies on several different data matrices. First, the covariates are partitioned based on whether their effects are homogeneoues (`X`), heterogeneoues (`G`), or negligible (`snps`). In practice, we usually take `X` to be null, because including homogeneoues covariates in `G` does not cause problems in our simulations (Dahl et al. 2018).

We generate the heterogeneoues covariates and SNPs in one of the simplest ways possible:
```R
N <- 3e3 # sample size

K <- 2   # number subtypes
z <- sample( K, N, replace=TRUE ) 

X0<- matrix( rnorm(N), N, 1 ) # null covar
X <- matrix( rnorm(N), N, 1 ) # hom covar
G <- matrix( rnorm(N), N, 1 ) # het covar

S <- 1e4 # number SNPs
snps <- matrix( rbinom(N*S,2,.1), N, S )
```
In the last line, I defined the subtypes `Z`, which are used to simulate the data. RGWAS aims to recover `Z`.

The second type of data RGWAS uses are phenotype matrices: `Y` for quantitative traits and `Yb` for binary traits. First, I simulate quantitative traits, again in a very simple (and computationally inefficient!) way:
```R
P <- 20  # number binary+quantitative traits
Y0  <- matrix( NA, N, P )
alpha <- rnorm(P) # homogeneoues effects
beta <- matrix( rnorm(K*P), K, P )    # heterogeneoues effects
submeans <- matrix( rnorm(K*P), K, P )
for( i in 1:N ) # very very slow
  Y0[i,] <- submeans[z[i],] + 1/Q0*X[i,] %*% alpha + 1/Q1*G[i,] %*% beta[z[i],] + rnorm(P)
round(rhomat <- cor(Y0),2)
round(rhomat1 <- cor(Y0[z==1,]),2)
round(rhomat2 <- cor(Y0[z==2,]),2)
mean(abs(rhomat[upper.tri(rhomat)]))   # smaller
mean(abs(rhomat1[upper.tri(rhomat1)])) # larger
mean(abs(rhomat2[upper.tri(rhomat2)])) # also larger
```
Note that I added a mean subtype effect

To make binary traits, I'll just threshold some columns of `Y0`, implicitly treating them as liabilities:
```R
bphens <- 1:(P/2)
Yb  <- apply( Y0[,bphens], 2, function(y) as.numeric( y > quantile(y,.8) ) )
qphens <- P/2+1:(P/2)
Y   <- Y0[,qphens]
```

## Running MFMR

Now I run MFMR on the traits and covariates. Imaginging that I don't know which columns in `X` and `G` are null, homogeneoues or heterogeneoues, I combine all into the putatitively heterogeneoues covariates inside MFMR (the `G` argument). I also add an intercept column to capture mean subtype effects, which are extremely helpful in practice--this is the entire signal driving covariate-unaware methods, like Gaussian mixture models of k-means.

```R
out <- mfmr( Yb, Y, cbind(1,X0,X,G), K=2 )
```
In this extremely simple simulation, MFMR seems to converge to the same likelihood mode for each of the `nrun=10` random restarts. In practice, however, random restarts is usually important: it does not guarantee global likelihood maximization, but it dramatically reduces the probability of obtaining a practically useless mode.

### Assessing subtype estimates

To see whether MFMR provides a reasonable estimate of the subtypes, I calculate the R-squared between true and estimated subtype probabilities:
```R
cor( out$pmat[,1], z )^2
```
Though intuitive for our example with `K=2` subtypes, simple correlation is not generally a useful metric for assessing similarity between proportions.

### Assessing regression estimates

I can also see whether MFMR estimated the regression coefficients accurately. First, the homogeneoues effects should resemble the subtype-specific effects estimates in both groups:
```R
# truly-hom column of cbind(1,X0,X,G)
qalphahat1 <- out$out$beta[1,3,qphens]
qalphahat2 <- out$out$beta[2,3,qphens]
cor( alpha[qphens], qalphahat1 )
cor( alpha[qphens], qalphahat2 )
```

To assess the heterogeneoues effects requires matching up the labels in MFMR to the true, simulated labels. I.e. `beta[1,]` may correspond to the estimated `betahat[2,]` because of label swapping. So I just compute the error metrics for each labelling option (because there are only 2 when `K=2`):
```R
# truly-het column of cbind(1,X0,X,G)
qbetahat <- out$out$beta[1,4,qphens]
cor( beta[1,qphens], qbetahat )^2
cor( beta[2,qphens], qbetahat )^2
```

## Testing large-effect covariates

Covariates with broad phenotypic effects are difficult to test because they can perturb subtype estimates if improperly modelled. Our test for these covariates involves treating each, in turn, as homogeneoues, ,....
```R
dropout  <- droptest( Yb, Y, cbind(1,X0,X,G), test_inds=2:4, K=2 )

round( dropout$pvals['Hom',2 ,qphens], 3 ) # all truly null
round( dropout$pvals['Hom',3:4,qphens], 3 ) # all truly alternate

round( dropout$pvals['Het',2:3,qphens], 3 ) # all truly null
round( dropout$pvals['Het',4  ,qphens], 3 ) # all truly alternate
```

## Testing small-effect covariates


