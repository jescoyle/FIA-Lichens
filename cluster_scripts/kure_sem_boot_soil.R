# This script runs bootstrap parameter estimation for the Soil SEM model

library(lavaan,  lib.loc="/nas02/home/j/r/jrcoyle/Rlibs/")

load('soil_model.Rdata')

endfit_std = bootstrapLavaan(endfit, R=10000, FUN=function(x) c(parameterEstimates(x)$est,standardizedSolution(x)$est.std))

save.image('soil_model.Rdata')