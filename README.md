# RGWAS

RGWAS is a two-step approach to find and then validate novel subtypes in multi-trait datasets, and is described in detail in [Dahl et al 2019 Plos Genetics](https://doi.org/10.1371/journal.pgen.1008009)

There are two key steps: Step 1 is to fit subtypes with `mfmr`; Step 2 is to test covariates for subtype-specific effects with `droptest` or `interxn_test`.

We also have unpublished tests for Step 2 when Step 1 is applied only in cases, which allows disease subtyping: `droptest_mlr` and `interxn_test_mlr`. In this setting, all controls are grouped together.

## Installation
R CMD INSTALL rgwas_1.0.tar.gz

## MFMR

`mfmr` takes covariates and traits (binary and/or continuous) as input and produces subtypes as output.

The covariates can have different effects in different subtypes, which can be leveraged to improve subtype estimates:
  ...
  ,,,,

## Large-effect covariate test 

`droptest` performs the RGWAS test for large-effect covariates. It fits MFMR for each covariate in turn, not allowing the tested covariate to have an heterogeneous effect across subtypes in MFMR (which overfits, see Supplementary Figure 1 of our RGWAS paper).

## Small-effect covariate test (MFMRX)

`interxn_test` performs the RGWAS test for small-effect covariates, which uses a single fit of MFMR that ignores the tested covariates. This is computationally efficient and enables GWAS-scale tests.

# Example analysis of simple simulated data

## Simulating data

RGWAS relies on several different data matrices. First, the covariates are partitioned based on whether their effects are homogeneoues (`X`), heterogeneoues (`G`), or negligible (`snps`). In practice, we usually take `X` to be null, because including homogeneoues covariates in `G` does not cause problems in our simulations (Dahl et al. 2018).

We generate the heterogeneoues covariates and SNPs in one of the simplest ways possible:
```R
library(rgwas)

N <- 3e3 # sample size
S <- 1e2 # number SNPs
K <- 2   # number subtypes

z <- sample( K, N, replace=TRUE ) # subtypes
G0<- matrix( rnorm(N*5), N, 5 )     # null covar
X <- matrix( rnorm(N), N, 1 )     # hom covar
G <- matrix( rnorm(N), N, 1 )     # het covar
snps <- matrix( rbinom(N*S,2,.1), N, S )
```
RGWAS aims to recover `z`, the true subtypes used to simulate the data.

The second type of data RGWAS uses are phenotype matrices: `Y` for quantitative traits and `Yb` for binary traits. First, I simulate quantitative traits, again in a very simple (and computationally inefficient!) way:
```R
P   <- 20  # number binary+quantitative traits
Y0  <- matrix( NA, N, P )
alpha <- rnorm(P) # homogeneoues effects
beta  <- matrix( rnorm(K*P), K, P )    # heterogeneoues effects
mus   <- matrix( rnorm(K*P), K, P )
for( i in 1:N )
  Y0[i,] <- X[i,] %*% alpha + G[i,] %*% beta[z[i],] + rnorm(P)
```
I added a mean subtype effect which (a) is usually realistic and (b) makes subtyping much easier.

To make binary traits, I'll treat some columns of `Y0` as liabilities and threshold them:
```R
bphens <- 1:3
qphens <- 4:P
Yb  <- apply( Y0[,bphens], 2, function(y) as.numeric( y > quantile(y,.8) ) )
Y   <- scale(Y0[,qphens])
```

## Running MFMR

Now I run MFMR on the traits and covariates. Imaginging that I don't know which columns in `X` and `G` are null, homogeneoues or heterogeneoues, I combine all into `covars` and treat them as putatitively heterogeneoues inside MFMR (the `G` argument). I also add an intercept column to capture mean subtype effects, which are extremely helpful in practice--this is the entire signal driving covariate-unaware methods, like Gaussian mixture models of k-means.

```R
covars <- cbind(1,G0,X,G)
out    <- mfmr( Yb, Y, covars, K=2 )
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
# truly-hom column of covars
qalphahat1 <- out$out$beta[1,3,qphens]
qalphahat2 <- out$out$beta[2,3,qphens]
cor( alpha[qphens], qalphahat1 )
cor( alpha[qphens], qalphahat2 )
```

To assess the heterogeneoues effects requires matching up the labels in MFMR to the true, simulated labels. I.e. `beta[1,]` may correspond to the estimated `betahat[2,]` because of label swapping. So I just compute the error metrics for each labelling option (because there are only 2 when `K=2`):
```R
# truly-het column of covars
qbetahat <- out$out$beta[1,4,qphens]
cor( beta[1,qphens], qbetahat )^2
cor( beta[2,qphens], qbetahat )^2
```

## Testing large-effect covariates

Covariates with broad phenotypic effects are difficult to test because they can perturb subtype estimates if improperly modelled. Our proposal is to refit MFMR for each tested large-effect covariate in turn, treating the tested covariate as homogeneoues within MFMR to balance over- and under-fitting the covariate effect (see our paper for details). In practice, this means the test is best performed with a specialized for loop, which I implemented in `droptest`:
```R
dropout  <- droptest( Yb, Y, covars, test_inds=2:4, K=2 )

round( dropout$pvals['Hom',2  ,qphens], 3 ) # all truly null
round( dropout$pvals['Hom',3:4,qphens], 3 ) # all truly alternate

round( dropout$pvals['Het',2:3,qphens], 3 ) # all truly null
round( dropout$pvals['Het',4  ,qphens], 3 ) # all truly alternate
```

## Testing small-effect covariates

Covariates with negligible phenotypic effects, on the other hand, can be ignored when fitting MFMR. This is very similar to fitting linear mixed models with variance components learned assuming each individual SNP has roughly zero effect. In this scenario, testing can be done independently of fitting MFMR. So I take the original MFMR fit, treating all covariates in `covars` as heterogeneoues, i.e. `out`. Then, using these subtypes, I perform standard fixed effect tests for SNP-subtype interaction:
```R 
# condition on all covariate-subtype interactions with Khatri-Rao product:
covars_x_E  <- t(sapply( 1:N, function(i) covars[i,,drop=FALSE] %x% out$pmat[i,,drop=FALSE] ))

pmat   <- out$pmat[,-1,drop=FALSE] # full-rank version of pmat
pvalsq <- apply( snps, 2, function(g) interxn_test( covars_x_E, Y[,1], g, pmat, bin=FALSE )$pvals )
ks.test( pvalsq['Hom',] )
ks.test( pvalsq['Het',] )

pvalsb <- apply( snps, 2, function(g) interxn_test( covars_x_E, Yb[,1], g, pmat, bin=TRUE )$pvals )
ks.test( pvalsb['Hom',] )
ks.test( pvalsb['Het',] )
```
This can be performed for all traits with an `apply` function (e.g. mclapply from the `parallel` R package) or for loop. Because MFMR does not need to be refit, testing with `interxn_tst` just amounts to performing t/F-tests for linear/logistic regression, for which `interxn_tst` is essentially just a (hopefully) helpful interface.

To simulate true positive results for SNP heterogeneity, I add a small effects of SNP 1: a homogeneoues effect on quantitative trait 2, and a heterogeneoues effect on quantitative trait 1:
```R 
# scale so effect sizes are more easily interpretable
Y    <- scale(Y) 
snps <- scale(snps)

# add small g effect to first two quantitative traits
Y[z==1,1]<- Y[z==1,1] + snps[z==1,1] * .05 * sqrt(2) # het
Y[    ,2]<- Y[    ,2] + snps[    ,1] * .05           # hom

# refit pmat with new, perturbed phens
out <- mfmr( Yb, Y, covars, K=2 )

# test with new pmat and phens
pmat    <- out$pmat[,-1,drop=FALSE] # full-rank version of pmat
pvalsq1 <- apply( snps, 2, function(g) interxn_test( covars_x_E, Y[,1], g, pmat, bin=FALSE )$pvals )
pvalsq2 <- apply( snps, 2, function(g) interxn_test( covars_x_E, Y[,2], g, pmat, bin=FALSE )$pvals )

pvalsq1[,1] # signif for Hom and Het
pvalsq2[,1] # signif for Hom only

# everything other than binary and quantitative trait 1 is (rightly) null
ks.test( pvalsq1['Hom',-1], 'punif' )$p
ks.test( pvalsq1['Het',-1], 'punif' )$p
ks.test( pvalsq2['Hom',-1], 'punif' )$p
ks.test( pvalsq2['Het',-1], 'punif' )$p
```

## Choosing the number of subtypes (K)

Choosing K is not generally straightforward. Ultimately, we heavily rely on simulations in Dahl et al 2019 Plos Gen that suggest downstream inference is calibrated for any (reasonble) choice of K, and our goal in choosing K is to maximize power. We strongly suggest sensitivity analyses in practice, meaning evaluating results after increasing or decreasing K by 1.

For modest sample sizes (eg <10K), we recommend choosing K to maximize out-of-sample likelihood. This can be accomplished with `score_K`, which performs `n.fold`-fold cross-validation:
```R 
meanll  <- numeric(3)
for( K in 1:3 )
        meanll[K]       <- median( score_K( G=cbind(G,G0), X=X, Yb=Yb, Yq=Y, K=K, n.folds=3 )[,1] ) ### n.folds=3 just for illustration
meanll # maximized at K=2, which is true in this simple simulation
[1] -22046.82 -21375.95 -21549.32
```
Note this liable to fit K that is too large, as the penalty for superfluous clusters is low (eg there is no likelihood cost to adding a cluster with weight 0). We prefer to err on the side of conservatism, meaning choosing lower values of K when multiple choices of K give similar likelihoods.

Robust lower bounds on K can be established with prediction strength metrics at large sample sizes. However, this comes at the cost of splitting the sample into halves, which can be prohibitive for small sample sizes.
