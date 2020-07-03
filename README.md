# RGWAS

RGWAS is a two-step approach to find and validate novel subtypes in multi-trait datasets:
1. Find subtypes by clustering traits with `mfmr`
2. Validate subtypes by testing for subtype-specific covariate effects with `droptest` or `interxn_test`

RGWAS is described fully in [Dahl et al 2019 Plos Genetics](https://doi.org/10.1371/journal.pgen.1008009). We also have unpublished tests for Step 2 when Step 1 is applied only in cases, which allows disease subtyping (Scenario 3 below).

### Installation
R CMD INSTALL rgwas_1.0.tar.gz

# Step 1: Learning subtypes with MFMR

`mfmr` takes covariates and traits (binary and/or continuous) as input and produces subtypes as output.

The covariates can have different effects in different subtypes, which can be leveraged to improve subtype estimates:
 - Traits are assumed to differ in distribution between subtypes.
 - Covariates are assumed to differ in effects between subtypes (but not necessarily in distribution).

Intuitively, the goal is to learn structure in the traits after accounting for covariates. This removes confounding structure from the data (eg sex, age, or ethnicity). This adjustment can both (a) solve false positive subtypes and (b) enable discovery of subtle true positive subtypes.
 
# Step 2: Testing subtypes for heterogeneous covariates

## Scenario 1: Covariates with small effects

Covariates with small effects can be ignored in the clustering step to infer subtypes. Importantly, this enables genome-wide testing.

In this scenario, subtypes from Step 1 are fixed and directly input to `interxn_test`, along with the tested covariate.

## Scenario 2: Covariates with large effects

When covariate effects are large, they cannot be ignored during clustering and can bias inferred subtypes due to overfitting.

In this scenario, `droptest` should be used, which internally re-fits the subtypes for each covariate to avoid overfitting with a score-like test. 

Note: these tests do not directly correspond to the subtypes inferred in Step 1, which would cause overfitting. Nonetheless, this provides a calibrated test for covariate heterogeneity that can be used for the subtypes which are inferred in Step 1. This technical complexity is essential for calibrated inference, but is not conceptually important.

## Scenario 3: Case-only subtypes

In case/control datasets with multiple measured secondary traits, it is possible to specifically learn subtypes of only the cases in Step 1. This is especially useful for case-only secondary traits (e.g. drug response, or tumor features).

In this scenario, the above interaction tests are meaningless, as subtype status is perfectly collinear with case/control status. Instead, we propose tests based on multinomial logistic regression. This enables testing whether a covariate differentially pushes people towards a particular disease subtype. For example, BMI is likely a more salient risk factor for type 2 diabetes than type 1 diabetes.

In this case, the above `droptest` and `interxn_test` should be replaced by `droptest_mlr` and `interxn_test_mlr`, respectively.

Note: This approach is susceptible to detecting disease-irrelevant differences between subtypes, e.g. [hair color-based subtypes of T2D](https://www.nature.com/articles/ng.3751).

# Example analysis of simple simulated data

## Simulating data

I simulate covariates that are null (`G0`), homogeneoues (`X`), heterogeneoues (`G`), or negligible (`snps`):
```R
set.seed(1234)
library(rgwas)

N <- 5e3 # sample size
S <- 1e2 # number SNPs
K <- 2   # number subtypes

z <- sample( K, N, replace=TRUE ) # subtypes
G0<- matrix( rnorm(N*5), N, 5 )   # null covar
X <- matrix( rnorm(N), N, 1 )     # hom covar
G <- matrix( rnorm(N), N, 1 )     # het covar
snps <- matrix( rbinom(N*S,2,.1), N, S )
```

I now simulate quantitative phenotypes using these covariates:
```R
P   <- 30  # number binary+quantitative traits
Y0  <- matrix( NA, N, P )
alpha <- rnorm(P) # homogeneoues effects
beta  <- matrix( rnorm(K*P), K, P )    # heterogeneoues effects
mus   <- matrix( rnorm(K*P)/2, K, P )
for( i in 1:N )
  Y0[i,] <- mus[z[i],] + X[i,] %*% alpha + G[i,] %*% beta[z[i],] + rnorm(P)
```
The goal of RGWAS is to recover the true subtypes `z` and  effects sizes `alpha` and `beta` in Step 1. In Step 2, the goal of RGWAS is to test whether the covariate effects are truly heterogeneous, i.e. to distinguish `alpha` and `beta`.

Note: I added subtype-specific means (`mus`), which is usually realistic and makes subtyping much easier. In particular, methods like Gaussian Mixture Models and k-means effectively use only the signal in `mus` to learn subtypes.

To make binary traits, I'll treat some columns of `Y0` as liabilities and threshold them:
```R
bphens <- 1:3
qphens <- 4:P
Yb  <- apply( Y0[,bphens], 2, function(y) as.numeric( y > quantile(y,.8) ) )
Yq  <- scale(Y0[,qphens])
```

## Step 1: Running MFMR

Now I run MFMR on the traits and covariates. All covariates are lumped together into `covars` and treated as heterogeneous by MFMR--in practice, it is not known which covariates are null/homogeneous/heterogeneous a priori:
```R
covars <- cbind(1,X,G,G0)
out    <- mfmr( Yb, Yq, covars, K=2, nrun=3 )
```
Note: In this simple simulation, MFMR converges to the same likelihood for each of the `nrun=3` random restarts. In practice, however, random restarts are important to consistently obtain useful solutions.

### Assessing the MFMR fit

The R-squared between true subtypes and the estimated subtype probabilities is:
```R
cor( out$pmat[,1], z )^2                   # [1] 0.9748958
```

To check that MFMR estimates the regression coefficients accurately:
```R
bhats   <- out$out$beta

# 1st component of bhats estimates intercepts:
cor( mus[1,qphens], bhats[1,1,qphens] )^2  # [1] 0.5339208

# 2nd component of bhats estimates truly homogenous effects:
cor( alpha[qphens], bhats[1,2,qphens] )^2  # [1] 0.9557434

# 3rd component of bhats estimates truly heterogenous effects:
cor( beta[1,qphens], bhats[1,3,qphens] )^2 # [1] 0.9725267
cor( beta[1,qphens], bhats[2,3,qphens] )^2 # [1] 0.04442622
```
Some notes:
- Though intuitive for `K=2` subtypes, correlation is not generally appropriate for comparing proportions.
- The correlation is calculated specifically on the quantitative phenotypes to avoid worrying about the liability scale transformation, but similar results are obtained by examining `bphens` instead.
- The mean effects in `mus` are correlated with `bhats` for either subtype because the phenotypes are centered.
- The homogeneous effects in `alpha` are correlated with `bhats` for either subtype.
- The heterogeneous effects in `beta` are only correlated with `bhats` for one subtype. Due to label switching, the true subtype labels may not match inferred subtype labels--i.e. `beta[1,]` may match `bhats[1,3,]` or `bhats[2,3,]`.

## Step 2: Testing covariate heterogeneity

### Testing large-effect covariates

Covariates with large effects are difficult to test because they perturb subtype estimates. To balance over- and under-fitting this effect, RGWAS refits MFMR while treating the tested covariate as homogeneoues. This is implemented in `droptest`:
```R
dropout  <- droptest( Yb, Yq, covars, test_inds=2:4, K=2 )

### Effects for X and G truly exist and test should be non-null:
quantile( dropout$pvals['Hom',2:3,qphens] )
#           0%           25%           50%           75%          100% 
# 0.000000e+00  0.000000e+00 4.062414e-278  3.367239e-73  6.501983e-01 

### Effects for X and G truly do not exist and test should be null:
quantile( dropout$pvals['Hom',4. ,qphens] )
#         0%        25%        50%        75%       100% 
# 0.03665779 0.22945576 0.52966551 0.71049620 0.98079611 

### Effects for X are truly Hom and so Het test should be null:
quantile( dropout$pvals['Het',2,qphens] )
#         0%        25%        50%        75%       100% 
# 0.01771635 0.27058785 0.48521272 0.71587194 0.98876730 

### Effects for X are truly Het and so Het test should be non-null:
quantile( dropout$pvals['Het',3,qphens] ) 
#           0%          25%          50%          75%         100% 
# 3.802443e-14 1.232082e-08 2.247578e-04 2.369354e-01 6.692128e-01 
```

### Testing small-effect covariates

Covariates with small effects do not meaningfully perturb subtype estimates and can be ignored when fitting MFMR. [EMMAX](https://www.nature.com/articles/ng.548) leverages a similar idea to expedite linear mixed models for GWAS. Therefore, RGWAS uses the original MFMR fit and performs standard fixed effect tests for SNP-subtype interaction:
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

Note: In this example, I adjust for all MFMR covariates in `covars` and conservatively treat them as heterogeneoues.


To simulate true positive results for SNP heterogeneity, I add a small effects of SNP 1: a homogeneoues effect on quantitative trait 2, and a heterogeneoues effect on quantitative trait 1:
```R 
# scale so effect sizes are more easily interpretable
Yq   <- scale(Yq) 
snps <- scale(snps)

# add small g effect to first two quantitative traits
Yq[z==1,1]<- Yq[z==1,1] + snps[z==1,1] * .05 * sqrt(2) # het
Yq[    ,2]<- Yq[    ,2] + snps[    ,1] * .05           # hom

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

Choosing K is not generally straightforward. Ultimately, we heavily rely on simulations in Dahl et al 2019 Plos Gen that suggest downstream inference is calibrated for any (reasonble) choice of K, and we view choosing K as maximizing power. We prefer to err on the side of conservatism, meaning choosing lower values of K when multiple choices of K give similar likelihoods. Further, we strongly suggest sensitivity analyses in practice, meaning evaluating results after increasing or decreasing K by 1.

For modest sample sizes (eg <10K), we recommend choosing K to maximize cross-validated likelihood with `score_K`:
```R 
mean_ll  <- numeric(3)
for( K in 1:3 )
  mean_ll[K]  <- mean( score_K( G=covars, Yb=Yb, Yq=Yq, K=K, n.folds=3, nrun=3 )[,'oos_ll'] ) ### n.folds=3 and nrun=3 are just for illustration
mean_ll # [1] -51080.62 -48649.89 -48727.61
```
Note: This liable to overfit K as the penalty for superfluous clusters is low (eg there is no likelihood cost to adding a cluster with weight 0). 

Robust lower bounds on K can be established with prediction strength metrics. However, this can be conservative because it requires splitting the data into halves. But, for larger sample sizes (or simple simulations like this), prediction strength can be a powerful way of robustly demonstrating that at least `K` subtypes exist and is implemented in `score_K_ps`:
```R 
mean_ps  <- numeric(3)
for( K in 2:3 )
  mean_ps[K]  <- mean( score_K_ps( G=covars, Yb=Yb, Yq=Yq, K=K, n.folds=3, nrun=3 ) )
mean_ps # [1] 0.0000000 0.9890182 0.7707128
```
