# This script dredges regional richness spatial error models

library(sp)
library(MuMIn)
library(spdep)

options(na.action='na.fail') # So that models are not fit to different subsets of data



load('regmods_objects_FIA.RData')


R_dredge = dredge(R_regmod, subset=dc('pseas_reg_mean', 'pseas_reg_mean2')&dc('regS_tree','regS_tree2'))

save(R_dredge, file='regmod_full_dredged.RData')




