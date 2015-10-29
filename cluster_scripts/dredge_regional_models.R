# This script dredges regional richness spatial error models

library(sp)
library(MuMIn)
library(spdep)

options(na.action='na.fail') # So that models are not fit to different subsets of data



load('regmods_objects_FIA.RData')


ROreg_dredge = dredge(RO_regmod, subset=dc('pseas_reg_mean', 'pseas_reg_mean2'))
RHreg_dredge = dredge(RH_regmod, subset=dc('regS_tree','regS_tree2'))


save(ROreg_dredge, RHreg_dredge, file='regmods_dredged.RData')




