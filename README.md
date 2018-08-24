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
