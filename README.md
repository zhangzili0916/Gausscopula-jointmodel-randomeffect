# Gausscopula-jointmodel-randomeffect
R code for "Joint modelling of longitudinal measurements and survival times via a copula approach" This repository includes R code for the simulation study and real data analysis. All the studies were based on RStudio of version 1.3.1073. The file "simulationstudy.R" includes the R codes to perform the simulation study and the file "realdatastudy.R" includes the R codes for really data application. Function "Y.lmmsurcondicopulajoineRcorrNOPLOTni" generate joint data of longitudinal measurements and survival times via the multivariate Gaussian. Function "Y.lmmsurcondicopulajoineRGHObest" perform the estimation process. The log-likelihhood functions "apploglikObmatrixni", "apploglikObuncormatrixni" and "apploglikObLongmatrixni" are required in the estimation process.
