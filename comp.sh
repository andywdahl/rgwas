R_LIBS=$HOME/R/libs 
export R_LIBS
R CMD build rgwas
R CMD INSTALL rgwas
#R CMD check 
R CMD check rgwas_1.0.tar.gz --no-manual

