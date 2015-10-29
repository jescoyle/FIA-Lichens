## This script dredges models on the computing cluster

options(na.action='na.fail') # So that models are not fit to different subsets of data


library(MASS) # glm.nb
library(MuMIn) # dredge
library(snow) # makeCluster
library(parallel)

use_data_test = read.csv('localmods_data.csv', row.names=1)

# Read in previously saved tables
modcompare = read.csv('Univariate model shapes GLM-NB.csv', row.names=1)

# Which variables have AIC supported concave-down relationships?
sq_vars = rownames(subset(modcompare, concavity=='down'&type=='quadratic'))

# Define sets of predictors to be used in models
climvars = c('mat','iso','pseas','wetness','rain_lowRH')
LOvars = c(climvars, 'radiation','bark_moist_pct.ba','wood_SG.ba','light.mean', 'PC1','bigTrees')
LHvars = c('bark_moist_pct.rao.ba','wood_SG.rao.ba','lightDist.mean','PIE.ba.tree','propDead','diamDiversity')
ROvars = paste(climvars, 'reg_mean', sep='_')
RHvars = c(paste(climvars, 'reg_var', sep='_'), 'regS_tree')
Rvars = c(RHvars, ROvars)
Lvars = c(LHvars, LOvars)
Hvars = c(RHvars, LHvars)
Ovars = c(ROvars, LOvars)
allvars = c(Rvars, Lvars)

### Local vs. regional variance partitioning of Local richness
# Define data for including squared terms
sqdata = use_data_test[,sq_vars]
sqdata = sqdata^2
colnames(sqdata) = paste(colnames(sqdata), 2, sep='')

full_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness','reg',allvars)], sqdata[,paste(sq_vars[sq_vars %in% allvars],2, sep='')]))

RregS_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness','reg',Rvars)], sqdata[,paste(sq_vars[sq_vars %in% Rvars],2, sep='')]))
R_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness',Rvars)], sqdata[,paste(sq_vars[sq_vars %in% Rvars],2, sep='')]))
regS_mod = glm.nb(richness~reg, data=use_data_test)
L_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness', Lvars)], sqdata[,paste(sq_vars[sq_vars %in% Lvars],2,sep='')]))

RH_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness',RHvars)], sqdata[,paste(sq_vars[sq_vars %in% RHvars],2, sep='')]), link='log')
RO_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness',ROvars)], sqdata[,paste(sq_vars[sq_vars %in% ROvars],2, sep='')]), link='log')

LH_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness',LHvars)], sqdata[,paste(sq_vars[sq_vars %in% LHvars],2, sep='')]), link='log')
LO_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness',LOvars)], sqdata[,paste(sq_vars[sq_vars %in% LOvars],2, sep='')]), link='log')

H_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness',Hvars)], sqdata[,paste(sq_vars[sq_vars %in% Hvars],2, sep='')]), link='log')
O_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness',Ovars)], sqdata[,paste(sq_vars[sq_vars %in% Ovars],2, sep='')]), link='log')

save(use_data_test, H_mod, Hvars, sqdata, sq_vars, file='cluster_data_H.RData')

# Set up cluster for parallel computing
c1 = makeCluster(6, 'SOCK')
clusterEvalQ(c1, library(MASS))
clusterEvalQ(c1, load('cluster_data_H.RData'))


## Dredge models

Hdredge = pdredge(H_mod, cluster=c1, beta=T, subset=dc('propDead','propDead2')&dc('wood_SG.rao.ba','wood_SG.rao.ba2')&dc('bark_moist_pct.rao.ba','bark_moist_pct.rao.ba2')&dc('pseas_reg_var', 'pseas_reg_var2')&dc('regS_tree', 'regS_tree2'), m.min=2)

save(Hdredge, file='local_rich_het_model_dredge.RData')

stopCluster(c1)


# Subset to best models
best_weight = Hdredge$weight[1] 
keep_models = which(Hdredge$weight/best_weight>= 0.05)
Hdredge_best = Hdredge[keep_models,]

save(Hdredge_best, file='local_rich_het_models_best.Rdata')




