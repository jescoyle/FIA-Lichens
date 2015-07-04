# This script is used to fit SEM models of lichen FIA data

source('./UNC/Projects/FIA Lichen/GitHub/FIA-Lichens/load_data.R')
source('./GitHub/FIA-Lichens/fia_lichen_analysis_functions.R')


################################################################################
### Path Analysis ###

#library(lavaan)
#library(semPlot)
#library(semTools)
#library(corrplot)


# Read in table of predictor variable types
predtypes = read.csv('predictors.csv', row.names=1)


#################################################################
### Determine variables to be used in SEM

library(MuMIn)
library(MASS)

# Read in data where variables are one their original scale, but some have been transformed to reduce skew (trans_data: see make_data_tables_for_analysis.R)
use_data = trans_data

# Define data sets for choosing variables (fit) and testing models (test)
test_data = use_data[testplots$yrplot.id,] 
fit_data = use_data[fitplots$yrplot.id,] 

clim_vars = c('ap','rh','mat','pseas')

## Define heterogeneity and optimality variables from PCA
RHclimvars = paste(clim_vars, 'reg_cv', sep='_')
ROvars = paste(clim_vars, 'reg_mean', sep='_')

RHPCA = prcomp(use_data[,RHclimvars], scale=T, center=T)
summary(RHPCA); RHPCA # all vars load positively on PC1 and it explains 62% of variance
RH_clim = predict(RHPCA)[,'PC1']

RH_for = use_data$regS_tree

ROPCA = prcomp(use_data[,ROvars], scale=T, center=T)
summary(ROPCA); ROPCA # PC1 explains 58% of variance and corresponds to increasing ap, rh and decreasing pseas, but mat not strongly correlated
RO_clim1 = predict(ROPCA)[,'PC1']
RO_clim2 = predict(ROPCA)[,'PC2']

LOclimvars = clim_vars
LOforvars = c('bigTrees', 'bark_moist_pct.ba', 'wood_SG.ba', 'light.mean', 'PC1')       
LHforvars = c('diamDiversity','bark_moist_pct.rao.ba','wood_SG.rao.ba','lightDist.mean','PIE.ba.tree','propDead')

LOclimPCA = prcomp(use_data[,LOclimvars], center=T, scale=T)
LOclimPCA; summary(LOclimPCA)
LO_clim1 = predict(LOclimPCA)[,'PC1']
LO_clim2 = predict(LOclimPCA)[,'PC2']

LOforPCA = prcomp(use_data[,LOforvars], center=T, scale=T)
summary(LOforPCA); LOforPCA # PC1 and 2 explain 50% of variance
LO_for1 = predict(LOforPCA)[,'PC1']
LO_for2 = predict(LOforPCA)[,'PC2'] 

LHforPCA = prcomp(use_data[,LHforvars], center=T, scale=T)
summary(LHforPCA); LHforPCA # PC1 and 2 explain 55% of variance
LH_for1 = -1*predict(LHforPCA)[,'PC1'] # increasing bark/wood/tree diversity and propDead
LH_for2 = predict(LHforPCA)[,'PC2'] # all increasing but most strongly size diversity and light diversity


# Make a dataframe for the SEM (including scaled richness/abundance for use in sem
sem_data = data.frame(model_data[c('lichen.rich','tot_abun_log', 'regS')], RH_for, RH_clim, RO_clim1, RO_clim2, LH_for1, LH_for2, LO_clim1, LO_clim2, LO_for1, LO_for2)
write.csv(sem_data, './SEM models/piecewise_sem_data.csv', row.names=T)

sem_data_test = sem_data[testplots$yrplot.id,]
write.csv(sem_data_test, './SEM models/piecewise_sem_data_test.csv', row.names=T)

#################################################################
### Piecewise SEM

sem_data = read.csv('./SEM models/piecewise_sem_data.csv', row.names=1)
sem_data_test = read.csv('./SEM models/piecewise_sem_data_test.csv', row.names=1)
sem_data_fit = sem_data[fitplots$yrplot.id,]

# This might not be an entirely appropriate approach since I think I have a model where the errors are correlated: 

library(devtools)
install_github("jslefche/piecewiseSEM")
library(piecewiseSEM)
library(sp)
library(spdep)

## Regional richness model is ML fit spatial error model fit by errorsarlm in spdep and assumes normal errors

# Define spatial dataframe for spatial analysis
spdata = master[rownames(sem_data_fit),c('LAT','LON')]
coordinates(spdata) = c('LON','LAT'); proj4string(spdata) = CRS("+proj=longlat")

# Write out test data set for Jon
play_data = cbind(spdata, sem_data_fit)
play_data = play_data[1:50, ]

write.csv(play_data, 'play_data.csv')

# Create a neighbors object with observations as neighbors if they are within 1000km (overlap in cells used to calculate reg variables)
reg_nb = dnearneigh(spdata, 0, 1000, row.names=rownames(sem_data_fit))
summary(reg_nb)

# Create a listw object with weights based on overlapping areas of 500km circles
# The weight is the percentage of the area shared with overlapping circle from another observation
A_intersect = function(d, R){ 2*(R^2)*acos(d/(2*R)) - 0.5*d*sqrt(4*(R^2) - (d^2))}
area_dist_func = function(d) A_intersect(d, 500)/(pi*(500^2))
reg_weights = lapply(nbdists(reg_nb,spdata), area_dist_func) 
hist(unlist(reg_weights)) # distribution of weights
reg_listw = nb2listw(reg_nb, glist=reg_weights, style='W')

# Re-scale variables to similar variances
# PCA based variables are already on a similar enough scale, need only fix regS and regS_tree
use_data = sem_data_fit
use_data$regS = use_data$regS/100
use_data$RH_for = use_data$RH_for/10


# Spatial error model
regS_mod_sp = errorsarlm(regS ~ RO_clim1 + RO_clim2 + I(RO_clim1^2) + I(RO_clim2^2) + RH_clim + RH_for, 
	data=use_data, listw=reg_listw)
# gives warning if variables are not re-scaled to be on similar scale
summary(regS_mod_sp)
res = resid(regS_mod_sp)
plot(res~use_data$regS)
hist(res)
moran.test(res, reg_listw)

# Regional forest ~ climate model
RHfor_mod_sp = errorsarlm(RH_for ~ RH_clim + RO_clim1 + RO_clim2 + I(RO_clim1^2) + I(RO_clim2^2), data=use_data, listw=reg_listw)
summary(RHfor_mod_sp)
res = resid(RHfor_mod_sp)
plot(res~use_data$RH_for)
hist(res)
moran.test(res, reg_listw)

# Non-spatial regional models
regS_mod = lm(regS ~ RO_clim1 + RO_clim2 + I(RO_clim1^2) + I(RO_clim2^2) + RH_clim + RH_for, data=use_data)
RHfor_mod = lm(RH_for ~ RH_clim + RO_clim1 + RO_clim2, data=use_data)
summary(regS_mod); summary(RHfor_mod)


# Local richness model - GLM with negative binomial error distribution
S_mod = glm.nb(lichen.rich ~ LO_clim1 + I(LO_clim1^2) + LO_clim2 + I(LO_clim2^2) + LO_for1 + I(LO_for1^2)+ LO_for2 + I(LO_for2^2) + LH_for1 + LH_for2 + regS, link='log', data=use_data)
summary(S_mod)
res = resid(S_mod)
plot(res~use_data$lichen.rich)
hist(res)

# Local forest ~ climate models (including regional effects)
LOfor1_mod = lm(LO_for1 ~ LO_clim1 + I(LO_clim1^2) + LO_clim2 + I(LO_clim2^2) , data=use_data)
LOfor2_mod = lm(LO_for2 ~ LO_clim1 + I(LO_clim1^2) + LO_clim2 + I(LO_clim2^2) , data=use_data)
LHfor1_mod = lm(LH_for1 ~ LO_clim1 + I(LO_clim1^2) + LO_clim2 + I(LO_clim2^2) , data=use_data)
LHfor2_mod = lm(LH_for2 ~ LO_clim1 + I(LO_clim1^2) + LO_clim2 + I(LO_clim2^2), data=use_data)
summary(LOfor1_mod);summary(LOfor2_mod);summary(LHfor1_mod);summary(LHfor2_mod)

# Models of regional climate effects on local climate
LOclim1_mod = lm(LO_clim1 ~ RO_clim1 + RO_clim2 + RH_clim, data=use_data)
LOclim2_mod = lm(LO_clim2 ~ RO_clim1 + RO_clim2 + RH_clim, data=use_data)
summary(LOclim1_mod); summary(LOclim2_mod)


# Fit piecewise SEM using piecewiseSEM package
full_modlist = list(regS_mod, RHfor_mod, S_mod, LOfor1_mod, LOfor2_mod, LHfor1_mod, LHfor2_mod, LOclim1_mod, LOclim2_mod)
full_sem = sem.fit(full_modlist, use_data)

subset(full_sem$missing.paths, p.value >0.05)



## Testing by hand in order to use spatial models

full_modlist = read.csv('./SEM models/Piecewise/full_modlist.csv')
reg2rich_modlist = read.csv('./SEM models/Piecewise/reg2rich_modlist.csv')


# A function to parse and fit a model in a model list table
make_mod = function(yvar, testvar, condvars, form, mod_data, listw=NA){
	require(stringr)

	if(is.na(condvars)) condvars=''
	covars = str_trim(unlist(strsplit(condvars, ';')))
	covars = c(covars, testvar)
	
	yvals = mod_data[,yvar]
	data_df = data.frame(mod_data[,covars])

	if(form=='lm'){
		mod = lm(yvals ~ ., data=data_df)
	}
	
	if(form=='errorsarlm'){
		mod = errorsarlm(yvals ~ ., data = data_df, listw=listw)
	}

	if(form=='glm.nb'){
		mod = glm.nb(yvals~ ., data = data_df, link='log')
	}

	mod
}


reg2rich_mods = sapply(1:nrow(reg2rich_modlist), function(i){
	x = reg2rich_modlist
	make_mod(yvar=x[i,'Var2'], testvar=x[i,'Var1'], condvars=x[i,'Z'], form=x[i,'Model'], mod_data=use_data, listw=reg_listw)
})


# A function that extracts the p-value of a given coefficient
reg2_pvals = sapply(1:nrow(reg2rich_modlist), function(i){
	x = reg2rich_modlist
	yvar=x[i,'Var2']
	testvar=x[i,'Var1']
	form = x[i,'Model']

	this_mod = reg2rich_mods[[i]]

	if(form %in% c('lm', 'glm.nb')){
		mod_tab = anova(this_mod)
		p = mod_tab[testvar, 'Pr(>F)']
	}

	if(form=='errorsarlm'){
		mod_tab = summary(this_mod)$Coef
		
		p = mod_tab[nrow(mod_tab), 'Pr(>|z|)']	
	}
	p
})


reg2rich_modlist$P = reg2_pvals
subset(reg2rich_modlist, P < 0.05)

# Calculate C statistic
log(reg2_pvals)


##################################################################
## Read in tables of parameter estimates and effects

load('./SEM models/Oct2014 No Measurement Error/regToRich_nopol_Fric_testdata_output.RData')

## TABLES FOR FRIC DID NOT GENERATE PROPERLY IN KURE DUE TO AN ERROR. WILL NEED TO RE-OUTPUT THESE LATER.


## noabun tables are from the reg2_noabun model

# All parameter estimares
allsp_ests = read.csv('./SEM models/Oct2014 No Measurement Error/nopol_AllSp_testdata_parameterEstimates.csv')
noabun_ests = read.csv('./SEM models/Oct2014 No Measurement Error/noabun_regToRich_nopol_AllSp_testdata_parameterEstimates.csv')
reg2rich_ests = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_AllSp_testdata_parameterEstimates.csv')
soil_ests = read.csv('./SEM models/Oct2014 No Measurement Error/soilmod_regTorich_nopol_AllSp_testdata_parameterEstimates.csv')
recip_ests = read.csv('./SEM models/Oct2014 No Measurement Error/nopol_recip_AllSp_testdata_parameterEstimates.csv')
reg2_recip_ests = read.csv('./SEM models/Oct2014sursement Error/regTorich_nopol_recip_AllSp_testdata_parameterEstimates.csv')

# Total effects
allsp = read.csv('./SEM models/Oct2014 No Measurement Error/nopol_AllSp_testdata_totaleffects.csv', row.names=1)
fric = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_fric_testdata_totaleffects.csv', row.names=1) # tables did not generate in KURE, need to do this separately
reg2 = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_AllSp_testdata_totaleffects.csv', row.names=1)
noabun = read.csv('./SEM models/Oct2014 No Measurement Error/noabun_regTorich_nopol_AllSp_testdata_totaleffects.csv', row.names=1)
soil = read.csv('./SEM models/Oct2014 No Measurement Error/soilmod_regTorich_nopol_AllSp_testdata_totaleffects.csv', row.names=1)
recip = read.csv('./SEM models/Oct2014 No Measurement Error/nopol_recip_AllSp_testdata_totaleffects.csv', row.names=1)
reg2_recip = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_recip_AllSp_testdata_totaleffects.csv', row.names=1)

# Direct effects
allsp_d = read.csv('./SEM models/Oct2014 No Measurement Error/nopol_AllSp_testdata_directeffects_richness.csv', row.names=1)
fric_d = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Fric_testdata_directeffects_richness.csv', row.names=1)
reg2_d = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_AllSp_testdata_directeffects_richness.csv', row.names=1)
noabun_d = read.csv('./SEM models/Oct2014 No Measurement Error/noabun_regTorich_nopol_AllSp_testdata_directeffects_richness.csv', row.names=1)
soil_d = read.csv('./SEM models/Oct2014 No Measurement Error/soilmod_regTorich_nopol_AllSp_testdata_directeffects_richness.csv', row.names=1)
recip_d = read.csv('./SEM models/Oct2014 No Measurement Error/nopol_recip_AllSp_testdata_directeffects_richness.csv', row.names=1)
reg2_recip_d = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_recip_AllSp_testdata_directeffects_richness.csv', row.names=1)

# Direct effect on abundance
#allsp_da = read.csv('./SEM models/No Pollution/nopol_AllSp_testdata_directeffects_abundance.csv', row.names=1)
#parm_da = read.csv('./SEM models/finalmod_Parm_testdata_directeffects_abundance.csv', row.names=1)
#phys_da = read.csv('./SEM models/finalmod_Phys_testdata_directeffects_abundance.csv', row.names=1)

# Direct effects on regional richness
allsp_dr = read.csv('./SEM models/Oct2014 No Measurement Error/nopol_AllSp_testdata_directeffects_regS.csv', row.names=1)
reg2_dr = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_AllSp_testdata_directeffects_regS.csv', row.names=1)
noabun_dr = read.csv('./SEM models/Oct2014 No Measurement Error/noabun_regTorich_nopol_AllSp_testdata_directeffects_regS.csv', row.names=1)
soil_dr = read.csv('./SEM models/Oct2014 No Measurement Error/soilmod_regTorich_nopol_AllSp_testdata_directeffects_regS.csv', row.names=1)
recip_dr = read.csv('./SEM models/Oct2014 No Measurement Error/nopol_recip_AllSp_testdata_directeffects_regS.csv', row.names=1)
reg2_recip_dr = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_recip_AllSp_testdata_directeffects_regS.csv', row.names=1)


# Indirect effects via abundance
allsp_i = read.csv('./SEM models/Oct2014 No Measurement Error/nopol_AllSp_testdata_indirecteffects_via_abundance.csv', row.names=1)
fric_i = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Fric_testdata_indirecteffects_via_abundance.csv', row.names=1)
reg2_i = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_AllSp_testdata_indirecteffects_via_abundance.csv', row.names=1)
soil_i = read.csv('./SEM models/Oct2014 No Measurement Error/soilmod_regTorich_nopol_AllSp_testdata_indirecteffects_via_abundance.csv', row.names=1)
recip_i = read.csv('./SEM models/Oct2014 No Measurement Error/nopol_recip_AllSp_testdata_indirecteffects_via_abundance.csv', row.names=1)
reg2_recip_i = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_recip_AllSp_testdata_indirecteffects_via_abundance.csv', row.names=1)


# Indirect effects of regional scale predictors
allsp_ir = read.csv('./SEM models/Oct2014 No Measurement Error/nopol_AllSp_testdata_regionalvars_indirecteffects.csv')
fric_ir = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Fric_testdata_regionalvars_indirecteffects.csv')
reg2_ir = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_AllSp_testdata_regionalvars_indirecteffects.csv')
noabun_ir = read.csv('./SEM models/Oct2014 No Measurement Error/noabun_regTorich_nopol_AllSp_testdata_regionalvars_indirecteffects.csv')
soil_ir = read.csv('./SEM models/Oct2014 No Measurement Error/soil_regTorich_nopol_AllSp_testdata_regionalvars_indirecteffects.csv')
recip_ir = read.csv('./SEM models/Oct2014 No Measurement Error/nopol_recip_AllSp_testdata_regionalvars_indirecteffects.csv')

# Indirect effects via forest structure
allsp_if = read.csv('./SEM models/Oct2014 No Measurement Error/nopol_AllSp_testdata_indirecteffects_via_forest.csv')
parm_if = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Parm_testdata_indirecteffects_via_forest.csv')
phys_if = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Phys_testdata_indirecteffects_via_forest.csv')
fric_if = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Fric_testdata_indirecteffects_via_forest.csv')


##################################################################
### Results ###

# How many predictors have significant total effects?
names(which(apply(allsp[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0))
names(which(apply(parm[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0))
names(which(apply(phys[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0))
names(which(apply(reg2[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0))
names(which(apply(soil[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0))
names(which(apply(recip[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0))

tot_sig = allsp[which(apply(allsp[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0),]
tot_sig = reg2[which(apply(reg2[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0),]
tot_sig = soil[which(apply(soil[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0),]
tot_sig = rev[which(apply(rev[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0),]

tot_sig = tot_sig[order(abs(tot_sig$std.all)),]

predtypes[rownames(tot_sig),]

# How many local-scale predictors have significant direct effects on local richness?
dir_sig = allsp_d[which(apply(allsp_d[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0),]
dir_sig = reg2_d[which(apply(reg2_d[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0),]
dir_sig = soil_d[which(apply(soil_d[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0),]
dir_sig = rev_d[which(apply(rev_d[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0),]

dir_sig = dir_sig[order(abs(dir_sig$std.all)),]
predtypes[rownames(dir_sig),]

# How many local-scale predictors have significant indirect effects on local richness?
indir_sig = reg2_i[which(apply(reg2_i[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0),]
indir_sig = soil_i[which(apply(soil_i[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0),]

indir_sig = indir_sig[order(abs(indir_sig$std.all)),]
predtypes[rownames(indir_sig),]

# How many regional-scale predictors have significant direct effects on regional richness?
dir_sig = reg2_dr[which(apply(reg2_dr[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0),]
dir_sig = dir_sig[order(abs(dir_sig$std.all)),]
predtypes[rownames(dir_sig),]

# Compare significance of total effects across models
tot_sig = data.frame(Base = apply(allsp[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0,
	Recip = apply(recip[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0,
	Reg2 = apply(reg2[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0,	
	Reg2_Recip = apply(reg2_recip[,c('std.ci.lower','std.ci.upper')], 1, prod)>0)
tot_vals = data.frame(Var = rownames(allsp), Base = allsp$std.all,
	Recip = recip$std.all,
	Reg2 = reg2$std.all,
	Reg2_Recip = reg2_recip$std.all)
plot(tot_vals$Base, tot_vals$Recip); abline(0,1); abline(h=0,v=0)
plot(tot_vals$Reg2, tot_vals$Reg2_Recip); abline(0,1); abline(h=0,v=0)
plot(tot_vals$Base, tot_vals$Reg2); abline(0,1); abline(h=0,v=0)


allsp_d[order(abs(allsp_d$std.all)),]

# Compare significance of total effects between model with and without abundance
rownames(reg2)==rownames(noabun)
tot_sig = data.frame(Reg2_sig = apply(reg2[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0,
	NoAbun_sig = apply(noabun[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0)
data.frame(predictor = rownames(reg2), Reg2 = reg2$std.all, NoAbun = noabun$std.all, tot_sig)

###########################################################################
### Assess model fit

# Load model specifications
source('./GitHub/FIA-Lichens/sem_model_specs.R')

# Load previously fit models (from Kure)

load('./SEM models/Oct2014 No Measurement Error/nopol_AllSp_testdata_output.RData')
load('./SEM models/Oct2014 No Measurement Error/noabun_nopol_AllSp_testdata_output.RData')
load('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_AllSp_testdata_output.RData')
load('./SEM models/Oct2014 No Measurement Error/regToRich_nopol_recip_AllSp_testdata_output.RData')

nopol_fit
noabun_nopol_fit
regTorich_nopol_fit
regTorich_nopol_recip_fit

AIC(nopol_fit, noabun_nopol_fit, regTorich_nopol_fit, regTorich_nopol_recip_fit)
anova(regTorich_nopol_fit, nopol_fit) # including direct paths to local richness makes the model significantly better.
anova(regTorich_nopol_recip_fit, regTorich_nopol_fit) # including reciprocal effects of local richness on regional richness significant improves model.

# Assess model fit
examine_mod = regTorich_nopol_fit 

summary(examine_mod, standardized=T, rsq=TRUE, fit.measures=T)
fitMeasures(examine_mod)
logLik(examine_mod)
hist(resid(examine_mod)$cov)
qqnorm(resid(examine_mod, type='standardized')$cov)

res = resid(examine_mod)$cov
which(abs(res)>.15, arr.ind=T)

res_order = order(abs(res), decreasing=T)
resdf = data.frame(residual=res[res_order], 
	var1=matrix(rownames(res), nrow=nrow(res), ncol=ncol(res))[res_order],
	var2=matrix(colnames(res), nrow=nrow(res), ncol=ncol(res), byrow=T)[res_order])

subset(resdf, abs(residual)>0.1)[]
subset(resdf, var1=='lichen.rich_log')
subset(resdf, var1=='regS'&abs(residual)>0.1)

subset(resdf, var1=='regS'& )


# Residual correlation figure: only plot residual correlations greater than .1
cortab = res
cortabsig = 1-abs(cortab)
png('./Figures/nopol regTorich residual correlations.png', height=1700, width=1700, type='cairo')
corrplot(cortab, method='square', type='upper', diag=F, 
	order='original', hclust.method='complete', p.mat=cortabsig,
	sig.level=.9, insig='blank', tl.cex=1.5, tl.col=1, cl.cex=2, mar=c(1,1,4,1))
dev.off()

# What are the residuals? Obs - Exp covariances
obscov = inspectSampleCov(path_regTorich_nopol, data=working_data_test)$cov #Make sure to change to test or fit data where appropriate
modcov = fitted(regTorich_nopol_fit)$cov 
covdiff =  obscov - modcov 
solid = abs(covdiff)>0.1
mypch = c(1,16)

richVar = 'lichen.rich_log'#'Phys_log'#
regVar = 'regS' #'regPhys'#

svg('./Figures/nopol regTorich recip covariance residuals.svg', height=6, width=6)
plot(obscov[covdiff!=0], modcov[covdiff!=0], xlab='Observed Covariances', 
	ylab='Model Covariances', las=1, pch = mypch[solid[covdiff!=0]+1])
abline(0,1)
abline(h=0,v=0)
points(obscov[richVar,], modcov[richVar,], pch=mypch[solid[richVar,]+1], col=2, cex=1.5)
points(obscov[regVar,], modcov[regVar,], pch=mypch[solid[regVar,]+1], col='blue', cex=1.1)
legend('topleft', c('Local richness','Regional richness'), col=c('red','blue'), 
	pch=16, pt.cex=c(1.5,1.1), bty='n')
mtext('Direct Regional Paths + Reciprocal Richness Model', 3, 0.5) 
dev.off()

cbind(obscov['lichen.rich_log',], modcov['lichen.rich_log',])

## Compare model predictions to observed data.
# Model covariance and mean structure
examine_mod = regTorich_nopol_fit
mu = fitted(examine_mod)$mean
Sigma= matrix(fitted(examine_mod)$cov, nrow=length(mu), ncol=length(mu))

# Random variables that follow model structure
nvars = length(mu); nobs = nrow(working_data_test)
rvars = matrix(rnorm(nvars*nobs, mean=mu), nrow=nvars, ncol=nobs)
L = chol(Sigma)
predvars = t(L) %*% matrix(rnorm(nvars*nobs), nrow=nvars, ncol=nobs)
predvars = data.frame(t(predvars)); colnames(predvars) = names(mu)
 
## Fit model to randomly generated data. What are the fit statistics?
examine_path = path_regTorich_nopol
random_fit = sem(examine_path, data=predvars, fixed.x=T, estimator='ML', se='robust.sem')
summary(random_fit, standardized=T, rsq=TRUE, fit.measures=T)
obsrand = inspectSampleCov(examine_path, data=predvars)$cov
modrand = fitted(random_fit)$cov 

svg('./Figures/No Measurement Error/compare regTorich random vs model residuals.svg', height=6, width=12)
par(mfrow=c(1,2))
plot(obsrand[covdiff!=0], modrand[covdiff!=0], xlab='Observed Covariances', 
	ylab='Model Covariances', las=1)
abline(0,1)
abline(h=0,v=0)
#points(obsrand['lichen.rich_log',], modrand['lichen.rich_log',], pch=16, col=2, cex=1.5)
#points(obsrand['regS',], modrand['regS',], pch=16, col='blue', cex=1.1)
#legend('topleft', c('Local richness','Regional richness'), col=c('red','blue'), 
#	pch=16, pt.cex=c(1.5,1.1), bty='n')
mtext('Random Data', 3, 0.5) 

plot(obscov[covdiff!=0], modcov[covdiff!=0], xlab='Observed Covariances', 
	ylab='Model Covariances', las=1)
abline(0,1)
abline(h=0,v=0)
#points(obscov['lichen.rich_log',], modcov['lichen.rich_log',], pch=16, col=2, cex=1.5)
#points(obscov['regS',], modcov['regS',], pch=16, col='blue', cex=1.1)
#legend('topleft', c('Local richness','Regional richness'), col=c('red','blue'), 
#	pch=16, pt.cex=c(1.5,1.1), bty='n')
mtext('Observed Data', 3, 0.5) 
dev.off()


## Residuals from Base Model with Pollution
pol_fit = sem(path_pol, data=working_data_test, fixed.x=T, estimator='ML', se='robust.sem')
obscov = inspectSampleCov(path_pol, data=working_data_test)$cov 
modcov = fitted(pol_fit)$cov 
covdiff =  obscov - modcov 
solid = abs(covdiff)>0.1
mypch = c(1,16)

svg('./Figures/finalmod covariance residuals.svg', height=6, width=6)
plot(obscov[covdiff!=0], modcov[covdiff!=0], xlab='Observed Covariances', 
	ylab='Model Covariances', las=1, pch = mypch[solid[covdiff!=0]+1])
abline(0,1)
abline(h=0,v=0)
points(obscov[richVar,], modcov[richVar,], pch=mypch[solid[richVar,]+1], col=2, cex=1.5)
points(obscov[regVar,], modcov[regVar,], pch=mypch[solid[regVar,]+1], col='blue', cex=1.1)
legend('topleft', c('Local richness','Regional richness'), col=c('red','blue'), 
	pch=16, pt.cex=c(1.5,1.1), bty='n')
points(obscov[c('totalNS','totalNS_reg'),], modcov[c('totalNS','totalNS_reg'),], pch=16, col='green', cex=1)
legend('bottomright', 'Pollution', col='green', pch=16, bty='n')
dev.off()










#################################################################################
### Variation partitioning of SEMs ###

path_local_noabun = "

	# Latent variables
	lichen_rich =~ sqrt(0.75)*lichen.rich_log
	
	# Local environment/climate effects on forest structure
	bark_moist_pct.rao.ba ~ fh1cm1*wetness + fh1cm2*rain_lowRH + fh1cm3*iso + fh1cm4*pseas + fh1cm5*mat + fh1cm6*radiation
	wood_SG.rao.ba ~ fh2cm1*wetness + fh2cm2*rain_lowRH + fh2cm3*iso + fh2cm4*pseas + fh2cm5*mat + fh2cm6*radiation
	LogSeed.rao.ba ~ fh3cm1*wetness + fh3cm2*rain_lowRH + fh3cm3*iso + fh3cm4*pseas + fh3cm5*mat + fh3cm6*radiation
	PIE.ba.tree ~ fh4cm1*wetness + fh4cm2*rain_lowRH + fh4cm3*iso + fh4cm4*pseas + fh4cm5*mat + fh4cm6*radiation
	propDead ~ fh5cm1*wetness + fh5cm2*rain_lowRH + fh5cm3*iso + fh5cm4*pseas + fh5cm5*mat + fh5cm6*radiation
	lightDist.mean ~ fh6cm1*wetness + fh6cm2*rain_lowRH + fh6cm3*iso + fh6cm4*pseas + fh6cm5*mat + fh6cm6*radiation
	diamDiversity ~ fh7cm1*wetness + fh7cm2*rain_lowRH + fh7cm3*iso + fh7cm4*pseas + fh7cm5*mat + fh7cm6*radiation
	bark_moist_pct.ba ~ fm1cm1*wetness + fm1cm2*rain_lowRH + fm1cm3*iso + fm1cm4*pseas + fm1cm5*mat + fm1cm6*radiation
	wood_SG.ba ~ fm2cm1*wetness + fm2cm2*rain_lowRH + fm2cm3*iso + fm2cm4*pseas + fm2cm5*mat + fm2cm6*radiation
	LogSeed.ba ~ fm3cm1*wetness + fm3cm2*rain_lowRH + fm3cm3*iso + fm3cm4*pseas + fm3cm5*mat + fm3cm6*radiation
	bigTrees ~ fm4cm1*wetness + fm4cm2*rain_lowRH + fm4cm3*iso + fm4cm4*pseas + fm4cm5*mat + fm4cm6*radiation
	light.mean  ~ fm5cm1*wetness + fm5cm2*rain_lowRH + fm5cm3*iso + fm5cm4*pseas + fm5cm5*mat + fm5cm6*radiation
	PC1 ~ fm6cm1*wetness + fm6cm2*rain_lowRH + fm6cm3*iso + fm6cm4*pseas + fm6cm5*mat + fm6cm6*radiation

	# Local climate effects on pollution
	totalNS ~ p1cm1*wetness + p1cm2*rain_lowRH

	# Effects on local lichen richness
	lichen_rich ~ r1fh1*bark_moist_pct.rao.ba + r1fh2*wood_SG.rao.ba + r1fh3*LogSeed.rao.ba +
		r1fh4*PIE.ba.tree + r1fh5*propDead + r1fh6*lightDist.mean + r1fh7*diamDiversity +
		r1fm1*bark_moist_pct.ba + r1fm2*wood_SG.ba + r1fm3*LogSeed.ba + r1fm4*bigTrees + r1fm5*light.mean + r1fm6*PC1 +
		r1cm1*wetness + r1cm2*rain_lowRH + r1cm3*iso + r1cm4*pseas + r1cm5*mat + r1cm6*radiation + r1R1*regS + r1p1*totalNS
		
	## Covariances among endogenous predictors in the same group
	# Dont need to specificy for Climate LO b/c these are exogenous and are calculated automatically in lavaan
	
	# Forest LH
	PIE.ba.tree ~~ bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba + propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba ~~ propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba ~~ wood_SG.rao.ba + LogSeed.rao.ba
	wood_SG.rao.ba ~~ LogSeed.rao.ba
	propDead ~~ lightDist.mean + diamDiversity
	lightDist.mean ~~ diamDiversity

	# Forest LO
	bark_moist_pct.ba + wood_SG.ba + LogSeed.ba ~~ bigTrees + light.mean + PC1
	bark_moist_pct.ba ~~ wood_SG.ba + LogSeed.ba 
	wood_SG.ba ~~ LogSeed.ba 
	bigTrees ~~ light.mean + PC1
	light.mean ~~ PC1
"
path_local = "

	# Latent variables
	lichen_rich =~ sqrt(0.75)*lichen.rich_log
	
	# Local environment/climate effects on forest structure
	bark_moist_pct.rao.ba ~ fh1cm1*wetness + fh1cm2*rain_lowRH + fh1cm3*iso + fh1cm4*pseas + fh1cm5*mat + fh1cm6*radiation
	wood_SG.rao.ba ~ fh2cm1*wetness + fh2cm2*rain_lowRH + fh2cm3*iso + fh2cm4*pseas + fh2cm5*mat + fh2cm6*radiation
	LogSeed.rao.ba ~ fh3cm1*wetness + fh3cm2*rain_lowRH + fh3cm3*iso + fh3cm4*pseas + fh3cm5*mat + fh3cm6*radiation
	PIE.ba.tree ~ fh4cm1*wetness + fh4cm2*rain_lowRH + fh4cm3*iso + fh4cm4*pseas + fh4cm5*mat + fh4cm6*radiation
	propDead ~ fh5cm1*wetness + fh5cm2*rain_lowRH + fh5cm3*iso + fh5cm4*pseas + fh5cm5*mat + fh5cm6*radiation
	lightDist.mean ~ fh6cm1*wetness + fh6cm2*rain_lowRH + fh6cm3*iso + fh6cm4*pseas + fh6cm5*mat + fh6cm6*radiation
	diamDiversity ~ fh7cm1*wetness + fh7cm2*rain_lowRH + fh7cm3*iso + fh7cm4*pseas + fh7cm5*mat + fh7cm6*radiation
	bark_moist_pct.ba ~ fm1cm1*wetness + fm1cm2*rain_lowRH + fm1cm3*iso + fm1cm4*pseas + fm1cm5*mat + fm1cm6*radiation
	wood_SG.ba ~ fm2cm1*wetness + fm2cm2*rain_lowRH + fm2cm3*iso + fm2cm4*pseas + fm2cm5*mat + fm2cm6*radiation
	LogSeed.ba ~ fm3cm1*wetness + fm3cm2*rain_lowRH + fm3cm3*iso + fm3cm4*pseas + fm3cm5*mat + fm3cm6*radiation
	bigTrees ~ fm4cm1*wetness + fm4cm2*rain_lowRH + fm4cm3*iso + fm4cm4*pseas + fm4cm5*mat + fm4cm6*radiation
	light.mean  ~ fm5cm1*wetness + fm5cm2*rain_lowRH + fm5cm3*iso + fm5cm4*pseas + fm5cm5*mat + fm5cm6*radiation
	PC1 ~ fm6cm1*wetness + fm6cm2*rain_lowRH + fm6cm3*iso + fm6cm4*pseas + fm6cm5*mat + fm6cm6*radiation

	# Local climate effects on pollution
	totalNS ~ p1cm1*wetness + p1cm2*rain_lowRH

	# Effects on local lichen richness
	lichen_rich ~ r1fh1*bark_moist_pct.rao.ba + r1fh2*wood_SG.rao.ba + r1fh3*LogSeed.rao.ba +
		r1fh4*PIE.ba.tree + r1fh5*propDead + r1fh6*lightDist.mean + r1fh7*diamDiversity +
		r1fm1*bark_moist_pct.ba + r1fm2*wood_SG.ba + r1fm3*LogSeed.ba + r1fm4*bigTrees + r1fm5*light.mean + r1fm6*PC1 +
		r1cm1*wetness + r1cm2*rain_lowRH + r1cm3*iso + r1cm4*pseas + r1cm5*mat + r1cm6*radiation + r1p1*totalNS +
		r1a1*tot_abun_log
	
	# Effects on local lichen abundance
	tot_abun_log ~ a1fh1*bark_moist_pct.rao.ba + a1fh2*wood_SG.rao.ba + a1fh3*LogSeed.rao.ba +
		a1fh4*PIE.ba.tree + a1fh5*propDead + a1fh6*lightDist.mean + a1fh7*diamDiversity +
		a1fm1*bark_moist_pct.ba + a1fm2*wood_SG.ba + a1fm3*LogSeed.ba + a1fm4*bigTrees + a1fm5*light.mean + a1fm6*PC1 +
		a1cm1*wetness + a1cm2*rain_lowRH + a1cm3*iso + a1cm4*pseas + a1cm5*mat + a1cm6*radiation + a1p1*totalNS 

	## Covariances among endogenous predictors in the same group
	# Dont need to specificy for Climate LO b/c these are exogenous and are calculated automatically in lavaan
	
	# Forest LH
	PIE.ba.tree ~~ bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba + propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba ~~ propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba ~~ wood_SG.rao.ba + LogSeed.rao.ba
	wood_SG.rao.ba ~~ LogSeed.rao.ba
	propDead ~~ lightDist.mean + diamDiversity
	lightDist.mean ~~ diamDiversity

	# Forest LO
	bark_moist_pct.ba + wood_SG.ba + LogSeed.ba ~~ bigTrees + light.mean + PC1
	bark_moist_pct.ba ~~ wood_SG.ba + LogSeed.ba 
	wood_SG.ba ~~ LogSeed.ba 
	bigTrees ~~ light.mean + PC1
	light.mean ~~ PC1
"

path_regional = "

	# Latent variables
	lichen_rich =~ sqrt(0.75)*lichen.rich_log
	
	# Regional climate, pollution, and regional forest heterogeneity effects on regional richness
	regS ~ R1CM1*wetness_reg_mean + R1CM2*rain_lowRH_reg_mean + R1CM3*iso_reg_mean + R1CM4*pseas_reg_mean + R1CM5*mat_reg_mean +
		R1CH1*wetness_reg_var + R1CH2*rain_lowRH_reg_var + R1CH3*iso_reg_var + R1CH4*pseas_reg_var + R1CH5*mat_reg_var +
		R1FH1*regS_tree + R1P1*totalNS_reg

	# Regional climate effects on regional forest heterogeneity
	regS_tree ~ FH1CM1*wetness_reg_mean + FH1CM2*rain_lowRH_reg_mean + FH1CM3*iso_reg_mean + FH1CM4*pseas_reg_mean + FH1CM5*mat_reg_mean +
		FH1CH1*wetness_reg_var + FH1CH2*rain_lowRH_reg_var + FH1CH3*iso_reg_var + FH1CH4*pseas_reg_var + FH1CH5*mat_reg_var

	# Regional climate effects on pollution
	totalNS_reg ~ P1CM1*wetness_reg_mean + P1CM2*rain_lowRH_reg_mean

	# Effects on local lichen richness
	lichen_rich ~ r1R1*regS 

"

# Includes lichen_rich ~ tot_abun_log b/c this is such a strong correlation.
path_regional_regTorich = "

	# Latent variables
	lichen_rich =~ sqrt(0.75)*lichen.rich_log
	
	# Regional climate, pollution, and regional forest heterogeneity effects on regional richness
	regS ~ R1CM1*wetness_reg_mean + R1CM2*rain_lowRH_reg_mean + R1CM3*iso_reg_mean + R1CM4*pseas_reg_mean + R1CM5*mat_reg_mean +
		R1CH1*wetness_reg_var + R1CH2*rain_lowRH_reg_var + R1CH3*iso_reg_var + R1CH4*pseas_reg_var + R1CH5*mat_reg_var +
		R1FH1*regS_tree + R1P1*totalNS_reg

	# Regional climate effects on regional forest heterogeneity
	regS_tree ~ FH1CM1*wetness_reg_mean + FH1CM2*rain_lowRH_reg_mean + FH1CM3*iso_reg_mean + FH1CM4*pseas_reg_mean + FH1CM5*mat_reg_mean +
		FH1CH1*wetness_reg_var + FH1CH2*rain_lowRH_reg_var + FH1CH3*iso_reg_var + FH1CH4*pseas_reg_var + FH1CH5*mat_reg_var

	# Regional climate effects on pollution
	totalNS_reg ~ P1CM1*wetness_reg_mean + P1CM2*rain_lowRH_reg_mean

	# Effects on local lichen richness
	lichen_rich ~ r1R1*regS + 
		r1CM1*wetness_reg_mean + r1CM3*rain_lowRH_reg_mean + r1CM3*iso_reg_mean + r1CM4*pseas_reg_mean + r1CM5*mat_reg_mean +
		r1CH1*wetness_reg_var + r1CH2*rain_lowRH_reg_var + r1CH3*iso_reg_var + r1CH4*pseas_reg_var + r1CH5*mat_reg_var +
		r1P1*totalNS_reg + r1FH1*regS_tree + r1a1*tot_abun_log
"

path_regional_regTorich_noabun = "

	# Latent variables
	lichen_rich =~ sqrt(0.75)*lichen.rich_log
	
	# Regional climate, pollution, and regional forest heterogeneity effects on regional richness
	regS ~ R1CM1*wetness_reg_mean + R1CM2*rain_lowRH_reg_mean + R1CM3*iso_reg_mean + R1CM4*pseas_reg_mean + R1CM5*mat_reg_mean +
		R1CH1*wetness_reg_var + R1CH2*rain_lowRH_reg_var + R1CH3*iso_reg_var + R1CH4*pseas_reg_var + R1CH5*mat_reg_var +
		R1FH1*regS_tree + R1P1*totalNS_reg

	# Regional climate effects on regional forest heterogeneity
	regS_tree ~ FH1CM1*wetness_reg_mean + FH1CM2*rain_lowRH_reg_mean + FH1CM3*iso_reg_mean + FH1CM4*pseas_reg_mean + FH1CM5*mat_reg_mean +
		FH1CH1*wetness_reg_var + FH1CH2*rain_lowRH_reg_var + FH1CH3*iso_reg_var + FH1CH4*pseas_reg_var + FH1CH5*mat_reg_var

	# Regional climate effects on pollution
	totalNS_reg ~ P1CM1*wetness_reg_mean + P1CM2*rain_lowRH_reg_mean

	# Effects on local lichen richness
	lichen_rich ~ r1R1*regS + 
		r1CM1*wetness_reg_mean + r1CM3*rain_lowRH_reg_mean + r1CM3*iso_reg_mean + r1CM4*pseas_reg_mean + r1CM5*mat_reg_mean +
		r1CH1*wetness_reg_var + r1CH2*rain_lowRH_reg_var + r1CH3*iso_reg_var + r1CH4*pseas_reg_var + r1CH5*mat_reg_var +
		r1P1*totalNS_reg + r1FH1*regS_tree 
"
# Fit models
# Change data set to test when ready to report final results

regional_fit =  sem(path_regional, data=working_data_fit, fixed.x=T, estimator='ML', se='robust.sem')
regional_regTorich_fit =  sem(path_regional_regTorich, data=working_data_fit, fixed.x=T, estimator='ML', se='robust.sem')
regional_regTorich_noabun_fit = sem(path_regional_regTorich_noabun, data=working_data_fit, fixed.x=T, estimator='ML', se='robust.sem')
local_fit = sem(path_local, data=working_data_fit, fixed.x=T, estimator='ML', se='robust.sem')
local_noabun_fit = sem(path_local_noabun, data=working_data_fit, fixed.x=T, estimator='ML', se='robust.sem')

regTorich_fit = sem(path_regTorich, data=working_data_fit, fixed.x=T, estimator='ML', se='robust.sem')
noabun_fit = sem(path_noabun, data=working_data_fit, fixed.x=T, estimator='ML', se='robust.sem')
finalmod_fit = sem(path_finalmod, data=working_data_fit, fixed.x=T, estimator='ML', se='robust.sem')

nopol_fit = sem(path_nopol, data=working_data_fit, fixed.x=T, estimator='ML', se='robust.sem')




## Comparisons:

## Variance partitionings here are not accurate because the full models include covariances 
## and paths between local and regional variables that are not in either of the sub models.

# path_local_noabun vs path_regional vs path_noabun: varpart R2
# formula for r2 = 1 - residual variance(lichen_rich)
# note: this differs from formula in Kline (see discussion: https://groups.google.com/forum/#!searchin/lavaan/r$20squared/lavaan/VrPYcTb_194/YxZxqW974AUJ)
summary(local_noabun_fit, rsq=T); r1=0.522
summary(regional_fit, rsq=T); r2=0.174
summary(noabun_fit, rsq=T); r3=0.651
rs = c(r1, r2, r3); names(rs) = c('Local','Regional','Full')
partvar2(rs)
# greater contribution of local predictors (b/c regional predictors mediated by regional richness)

# path_local vs path_regional_regTorich vs path_regTorich: varpart R2 # local abundance will show up in both models
summary(local_fit, rsq=T); r1=0.805
summary(regional_regTorich_fit, rsq=T); r2=.762
summary(regTorich_fit, rsq=T); r3=0.856
rs = c(r1, r2, r3); names(rs) = c('Local','Regional','Full')
partvar2(rs)
# Not very informative.

# path_local_noabun vs path_regional_regTorich_noabun vs path_noabun: varpart R2
summary(local_noabun_fit, rsq=T); r1=0.522
summary(regional_regTorich_noabun_fit, rsq=T); r2=.472
summary(noabun_fit, rsq=T); r3=0.651
rs = c(r1, r2, r3); names(rs) = c('Local','Regional','Full')
partvar2(rs)
# equal contributions of local and regional predictors


# path_local vs path_noabun_local: better to have abundance in local model
anova(local_fit, local_noabun_fit, test=T)

# path_local_noabun vs path_noabun: local model better than full
anova(noabun_fit, local_noabun_fit, test=T) # performs the Chi-square difference test which is not the LR test

# path_regional vs path_noabun: regional model better than full
anova(noabun_fit, regional_fit, test=T) # performs the Chi-square difference test which is not the LR test

# path_local vs path_regTorich: local model better
anova(regTorich_fit, local_fit, test=T) # performs the Chi-square difference test which is not the LR test

# path_regional_regTorich vs path_regTorich: regional model better than full
anova(regTorich_fit, regional_regTorich_fit, test=T) # performs the Chi-square difference test which is not the LR test

# path_local vs path_finalmod: local model better
anova(finalmod_fit, local_fit, test=T)

# finalmod vs regTorich: better to have direct paths in model
anova(regTorich_fit, finalmod_fit, test=T)


AIC(nopol_fit, regTorich_fit, finalmod_fit)

anova(regTorich_nopol_fit, nopol_fit)


#######################################################################################
### Figures ###



####################################
### Compare total vs direct effects in regTorich model


######################################
## Compare total effects
library(lattice)

# Color scheme
mycols = matrix(c('#b3b5ffff','#6b6dd7ff','#8dff94ff','#38af4fff'), nrow=2)
mycolsbw = c('grey80','white')
colnames(mycols) = c('regional','local')
rownames(mycols) = c('het','opt')
names(mycolsbw) = c('regional','local')
#mycols_trans = paste(mycols, '50', sep='')
#names(mycols_trans) = c('regional','local')
#plot(1:length(mycols),1:length(mycols), type='n'); text(1:length(mycols),1:length(mycols),labels=names(mycols), cex=2, col=mycols) # Check colors
mypcols = c('white','grey30')
mypch = c(22,23)
myadj=.15
#mytypes = expression('C'['H'],'C'['O'],'F'['H'],'F'['O'],'P','','') # symbols used in plot to denote variable types
#names(mytypes)=c('CH','CM','FH','FM','P','R','A')
#myshade = c('55','99','')
#names(myshade) = c('het','opt','')
mytypes = expression('C','F','P','','') # symbols used in plot to denote variable types
names(mytypes)=c('C','F','P','R','A')


# Make df of total effects including direct effects of regS and abundance

# Order variables from lowest to highest total effects
total = reg2 # This changes based on what response variable is being analyzed
total = rbind(total, reg2_d[c('regS','tot_abun_log'),]) #

ordered_vars = rownames(total[order(total$std.all),])

# Put tables in same order
use_total = total[ordered_vars,]

# Define range limits that will include 95% confidence intervals
myrange = range(use_total[,c('std.ci.lower','std.ci.upper')], na.rm=T)+c(-.04, .04)
myrange[1] = -0.85 #-1

# Make plot
svg('./Figures/No Measurement Error/Standardized total effects on AllSp richness nopol regTorich.svg', height=20, width=19)
dotplot(as.numeric(factor(rownames(use_total), levels = ordered_vars))~std.all, data=use_total, 
	xlab=list('Standardized Effect',cex=3), ylab='',
	main='',cex.lab=3,aspect=5/3, xlim=myrange,
	panel=function(x,y){
	
		# Add horizontal boxes
		colorcombos = predtypes[ordered_vars, c('mode','scale')]
		colorcombos['regS','mode'] = 'opt'
		colororder = apply(colorcombos, 1, function(x) mycols[x[1],x[2]])

		panel.rect(myrange[1]-0.01,1:length(ordered_vars)-.5, myrange[2]+0.01, 1:length(ordered_vars)+.5,
			col=colororder, border='grey50')
	
		# Add vertical line at 0
		panel.abline(v=c(-.2,0,.2), col='grey30', lty=2, lwd=2)		
		
		# Add null distribution segments for total effects
		panel.segments(use_total$std.ci.lower, y,
			use_total$std.ci.upper, y, 
			col='black', lwd=4.5, lend=1)
		# Add points for total estimated effects
		panel.points(x, y, col='black', fill=mypcols[1], pch=mypch[1], cex=3, lwd=3) 
	
		# Add text labeling the variable type
		vartypes =  sapply(predtypes[ordered_vars,'label'], function(x) toupper(substr(x, 1, 1))) 
		panel.text(myrange[1]+0.05, y, labels=mytypes[vartypes], cex=2)
		
	},
	scales=list(y=list(labels=varnames[ordered_vars,'midName'], 
		cex=3, col='black'),
		x=list(cex=3, tick.number=8))
)
dev.off()

### Direct effects on regional richness

# Define data set to use
direct_reg = reg2_recip_dr

# Order variables from lowest to highest direct effects
ordered_vars = rownames(direct_reg[order(direct_reg$std.all),])

# Put tables in same order
use_df = direct_reg[ordered_vars,]

# Define range limits that will include 95% confidence intervals
myrange = range(use_df[,c('std.ci.lower','std.ci.upper')], na.rm=T)+c(-.04, .04)
myrange[1] = -2

# Make plot
svg('./Figures/No Measurement Error/Standardized direct effects on AllSp regional richness regTorich nopol recip.svg', height=9, width=19)
dotplot(as.numeric(factor(rownames(use_df), levels = ordered_vars))~std.all, data=use_df, 
	xlab=list('Standardized Effect',cex=3), ylab='',
	main='',cex.lab=3,aspect=4/5, xlim=myrange,
	panel=function(x,y){
	
		# Add horizontal boxes
		colorcombos = predtypes[ordered_vars, c('mode','scale')]
		colorcombos['lichen.rich_log','mode'] = 'opt' # Use with reciprocal richness models
		colororder = apply(colorcombos, 1, function(x) mycols[x[1],x[2]])
		panel.rect(myrange[1]-0.01,1:length(ordered_vars)-.5, myrange[2]+0.01, 1:length(ordered_vars)+.5,
			col=colororder, border='grey50')
	

		# Add vertical line at 0
		panel.abline(v=c(-.2,0,.2), col='grey30', lty=2, lwd=2)		
		
		# Add null distribution segments for direct effects
		panel.segments(use_df$std.ci.lower, y,
			use_df$std.ci.upper, y, 
			col='black', lwd=4.5, lend=1)
		# Add points for direct estimated effects
		panel.points(x, y, col='black', fill=mypcols[1], pch=mypch[1], cex=3, lwd=3) 
	
		# Add text labeling the variable type
		#vartypes =  sapply(predtypes[ordered_vars,'label'], function(x) toupper(substr(x, 1, 1)))
		#panel.text(myrange[1]+0.1, y, labels=mytypes[vartypes], cex=2)
		
	},
	scales=list(y=list(labels=varnames[ordered_vars,'midName'], 
		cex=3, col='black'),
		x=list(cex=3, tick.number=8))
)
dev.off()

### Direct and indirect effects on local richness

# Define datasets to use
direct = reg2_recip_d
indirect = reg2_recip_i
direct = direct[c(rownames(indirect),'tot_abun_log'),]
#direct = direct[rownames(indirect),] # Use for plotting noabun model

# Order variables from lowest to highest direct effects
ordered_vars = rownames(direct[order(direct$std.all),])

# Put tables in same order
use_direct = direct[ordered_vars,]
use_indirect = indirect[ordered_vars,]

# Define range limits that will include 95% confidence intervals
myrange = range(rbind(use_direct,use_indirect)[,c('std.ci.lower','std.ci.upper')], na.rm=T)+c(-.04, .04)
myrange[1] = -.34

jitter = 0

# Make plot
svg('./Figures/No Measurement Error/Standardized direct and indirect effects on AllSp richness regTorich nopol recip.svg', height=13, width=19)
dotplot(as.numeric(factor(rownames(use_direct), levels = ordered_vars))~std.all, data=use_direct, 
	xlab=list('Standardized Effect',cex=3), ylab='',
	main='',cex.lab=3,aspect=5/4, xlim=myrange,
	panel=function(x,y){
	
		# Add horizontal boxes
		colorcombos = predtypes[ordered_vars, c('mode','scale')]
		colororder = apply(colorcombos, 1, function(x) mycols[x[1],x[2]])

		panel.rect(myrange[1]-0.01,1:length(ordered_vars)-.5, myrange[2]+0.01, 1:length(ordered_vars)+.5,
			col=colororder, border='grey50')

		# Add vertical line at 0
		panel.abline(v=c(-.2,0,.2), col='grey30', lty=2, lwd=2)		
		
		# Add null distribution segments for indirect effects
		panel.segments(use_indirect$std.ci.lower, y-jitter,
			use_indirect$std.ci.upper, y-jitter, 
			col='black', lwd=4.5, lend=1)
		# Add points for indirect estimated effects
		panel.points(use_indirect$std.all, y-jitter, col='black', fill=mypcols[2], pch=mypch[2], cex=3, lwd=3) 

		# Add null distribution segments for direct effects
		panel.segments(use_direct$std.ci.lower, y+jitter,
			use_direct$std.ci.upper, y+jitter, 
			col='black', lwd=4.5, lend=1)
		# Add points for direct estimated effects
		panel.points(x, y+jitter, col='black', fill=mypcols[1], pch=mypch[1], cex=3, lwd=3) 
	

		# Add text labeling the variable type
		#vartypes =  sapply(predtypes[ordered_vars,'label'], function(x) toupper(substr(x, 1, 1)))
		#panel.text(myrange[1]+.03, y, labels=mytypes[vartypes], cex=2)
		
	},
	scales=list(y=list(labels=varnames[ordered_vars,'midName'], 
		cex=3, col='black'),
		x=list(cex=3, tick.number=8)),
	key=list(x=.5, y=1, corner=c(.5,0), lines=list(type='o', pch=mypch, fill=mypcols, lwd=3, pt.lwd=3),
		text=list(c('Direct Effect','Indirect Effect')),
		background='white', cex=3, divide=1, padding.text=5, border=F, columns=2)
)
dev.off()



### Indirect effects via abundance on local richness

# Define datasets to use
indirect = reg2_i

# Order variables from lowest to highest direct effects
ordered_vars = rownames(indirect[order(indirect$std.all),])

# Put tables in same order
use_indirect = indirect[ordered_vars,]

# Define range limits that will include 95% confidence intervals
myrange = range(use_indirect[,c('std.ci.lower','std.ci.upper')], na.rm=T)+c(-.04, .04)
myrange[1] = -.34

jitter = 0

# Make plot
svg('./Figures/Standardized indirect effects on AllSp richness regTorich nopol.svg', height=13, width=19)
dotplot(as.numeric(factor(rownames(use_indirect), levels = ordered_vars))~std.all, data=use_indirect, 
	xlab=list('Standardized Effect',cex=3), ylab='',
	main='',cex.lab=3,aspect=5/4, xlim=myrange,
	panel=function(x,y){
	
		# Add horizontal boxes
		combos = predtypes[ordered_vars, c('mode','scale')]
		colororder = apply(combos, 1, function(x) mycols[x[1],x[2]])

		panel.rect(myrange[1]-0.01,1:length(ordered_vars)-.5, myrange[2]+0.01, 1:length(ordered_vars)+.5,
			col=colororder, border='grey50')

		# Add vertical line at 0
		panel.abline(v=c(-.2,0,.2), col='grey30', lty=2, lwd=2)		
		
		# Add null distribution segments for indirect effects
		panel.segments(use_indirect$std.ci.lower, y-jitter,
			use_indirect$std.ci.upper, y-jitter, 
			col='black', lwd=4.5, lend=1)
		# Add points for indirect estimated effects
		panel.points(use_indirect$std.all, y-jitter, col='black', fill=mypcols[2], pch=mypch[2], cex=3, lwd=3) 

		# Add text labeling the variable type
		vartypes =  sapply(predtypes[ordered_vars,'label'], function(x) toupper(substr(x, 1, 1)))
		panel.text(myrange[1]+.03, y, labels=mytypes[vartypes], cex=2)
		
	},
	scales=list(y=list(labels=varnames[ordered_vars,'midName'], 
		cex=3, col='black'),
		x=list(cex=3, tick.number=8))
)
dev.off()


###########################################
### Are totaleffects at regional scale larger than total effects at local scale?

total = reg2# This changes based on what response variable is being analyzed
total = rbind(total, reg2_d[c('regS','tot_abun_log'),1:4])

cvars = c('wetness','rain_lowRH','pseas','iso','mat')

lvars = c(cvars, 'PIE.ba.tree')
rvars = c(paste(cvars, 'reg_mean', sep='_'), 'regS_tree')

svg('./Figures/No Measurement Error/compare local-regional total effects nopol regTorich.svg', height=5.5, width=5.5)
par(mar=c(5,5,1,1))
par(lwd=2)
par(lend=1)
par(cex.lab=1.2)
plot(total[lvars, 'std.all'], total[rvars, 'std.all'], 	
	xlab='Local Effect', ylab='Regional Effect', type='n', las=1,
	xlim=c(-0.9,0.9), ylim=c(-0.9,0.9))
usr=par('usr')
polygon(usr[c(1,1,2,2)],usr[c(3,4,4,3)], col='#8dff94fb')
polygon(c(usr[1],0,usr[1],usr[2],0,usr[2],usr[1]),
	c(usr[3],0,usr[4],usr[4],0,usr[3],usr[3]), col='#b3b5ffff')
abline(h=0,v=0,lty=2, col='black')
arrows(total[lvars,'std.ci.lower'], total[rvars,'std.all'],
	total[lvars,'std.ci.upper'], total[rvars,'std.all'], col='black', code=3, angle=90, length=0.05)
arrows(total[lvars,'std.all'], total[rvars,'std.ci.lower'],
	total[lvars,'std.all'], total[rvars,'std.ci.upper'], col='black', code=3, angle=90, length=0.05)
points(total[lvars, 'std.all'], total[rvars, 'std.all'], pch=15, col='black')

# Labels for regTorich
text(total[lvars,'std.all'], total[rvars,'std.all'], varnames[lvars,'midName'],
	pos=c(2,2,4,4,2,4), offset=1.5, col='black')

# Labels for base model
#text(total[lvars,'std.all'], total[rvars,'std.all'], varnames[lvars,'midName'],
#	pos=c(3,2,4,3,1,4), offset=.9)

dev.off()





###########################################
## Make table of indirect climate effects

## Regional Scale Climate variables
indirect_ir = allsp_ir

order_clim = c('wetness','rain_lowRH','pseas','mat','iso','radiation')
order_clim_reg = c(paste(order_clim[1:5],'reg_mean', sep='_'), paste(order_clim[1:5], 'reg_var', sep='_'), 'regS_tree')

climReg_IE_tab = data.frame(var = order_clim_reg, predictor = varnames[order_clim_reg,'displayName'])
use_loc = subset(indirect_ir, IEvar1 %in% c('clim_loc','PIE.ba.tree')); rownames(use_loc) = use_loc$predictor
climReg_IE_tab$loc = use_loc[order_clim_reg, 'std.all']
use_regS = subset(indirect_ir, IEvar1=='regS'); rownames(use_regS) = use_regS$predictor
climReg_IE_tab$regS = use_regS[order_clim_reg,'std.all']
use_FH = subset(indirect_ir, IEvar1=='regS_tree'&is.na(IEvar2)); rownames(use_FH) = use_FH$predictor
climReg_IE_tab$FH = use_FH[order_clim_reg, 'std.all']

climReg_IE_tab$loc_sig = apply(use_loc[order_clim_reg,c('std.ci.lower','std.ci.upper')], 1, prod)>0
climReg_IE_tab$regS_sig = apply(use_regS[order_clim_reg,c('std.ci.lower','std.ci.upper')], 1, prod)>0
climReg_IE_tab$FH_sig = apply(use_FH[order_clim_reg,c('std.ci.lower','std.ci.upper')], 1, prod)>0

sub_vars_FH = subset(indirect_ir, (IEvar1=='regS_tree')&!is.na(IEvar2)); colnames(sub_vars_FH)[1] = 'var'; colnames(sub_vars_FH)[colnames(sub_vars_FH)=='std.all'] <- 'subpath_FH'
sub_vars_FH$subFH_sig = apply(sub_vars_FH[,c('std.ci.lower','std.ci.upper')], 1, prod)>0

climReg_IE_tab = merge(climReg_IE_tab, sub_vars_FH[,c('var','IEvar2','subpath_FH','subFH_sig')])
climReg_IE_tab = climReg_IE_tab = climReg_IE_tab[,c('var','IEvar2','predictor','loc','regS','FH','subpath_FH',
	'loc_sig','regS_sig','FH_sig','subFH_sig')]

write.csv(climReg_IE_tab, './Tables/No Measurement Error/Indirect regional climate effects nopol.csv', row.names=F)

# Table for regTorich model which includes direct effect
indirect_ir = reg2_ir
direct = reg2_d

order_clim = c('wetness','rain_lowRH','pseas','mat','iso','radiation')
order_clim_reg = c(paste(order_clim[1:5],'reg_mean', sep='_'), paste(order_clim[1:5], 'reg_var', sep='_'), 'regS_tree')

climReg_IE_tab = data.frame(var = order_clim_reg, predictor = varnames[order_clim_reg,'displayName'])
use_dir = direct[order_clim_reg,]
climReg_IE_tab$dir = use_dir$std.all
use_loc = subset(indirect_ir, IEvar1 %in% c('clim_loc','PIE.ba.tree')); rownames(use_loc) = use_loc$predictor
climReg_IE_tab$loc = use_loc[order_clim_reg, 'std.all']
use_regS = subset(indirect_ir, IEvar1=='regS'); rownames(use_regS) = use_regS$predictor
climReg_IE_tab$regS = use_regS[order_clim_reg,'std.all']
use_FH = subset(indirect_ir, IEvar1=='regS_tree'&is.na(IEvar2)); rownames(use_FH) = use_FH$predictor
climReg_IE_tab$FH = use_FH[order_clim_reg, 'std.all']

climReg_IE_tab$dir_sig = apply(use_dir[,c('std.ci.lower','std.ci.upper')], 1, prod)>0
climReg_IE_tab$loc_sig = apply(use_loc[order_clim_reg,c('std.ci.lower','std.ci.upper')], 1, prod)>0
climReg_IE_tab$regS_sig = apply(use_regS[order_clim_reg,c('std.ci.lower','std.ci.upper')], 1, prod)>0
climReg_IE_tab$FH_sig = apply(use_FH[order_clim_reg,c('std.ci.lower','std.ci.upper')], 1, prod)>0

sub_vars_FH = subset(indirect_ir, (IEvar1=='regS_tree')&!is.na(IEvar2)); colnames(sub_vars_FH)[1] = 'var'; colnames(sub_vars_FH)[colnames(sub_vars_FH)=='std.all'] <- 'subpath_FH'
sub_vars_FH$subFH_sig = apply(sub_vars_FH[,c('std.ci.lower','std.ci.upper')], 1, prod)>0

climReg_IE_tab = merge(climReg_IE_tab, sub_vars_FH[,c('var','IEvar2','subpath_FH','subFH_sig')])
climReg_IE_tab = climReg_IE_tab = climReg_IE_tab[,c('var','IEvar2','predictor','dir','loc','regS','FH','subpath_FH',
	'dir_sig','loc_sig','regS_sig','FH_sig','subFH_sig')]

write.csv(climReg_IE_tab, './Tables/No Measurement Error/Indirect regional climate effects regTorich nopol.csv', row.names=F)


################################################
## Plot direct vs. indirect effects via abundance for local predictors
use_vars = subset(predtypes, scale=='local'&label!=''&type%in%c('env','forest'))
use_vars = use_vars[-grep('soil', rownames(use_vars)),]

use_d = reg2_d[rownames(use_vars),]
use_i = reg2_i[rownames(use_vars),]

mypch = c(22,23)
mycol=c('white','grey30')

svg('./Figures/No Measurement Error/compare richness abundance effects local vars regTorich.svg', height=5, width=5 )
par(mar=c(4,4,1,1))
plot(use_d$std.all, use_i$std.all, type='n', 
	xlim=c(-.3, .3), ylim=c(-.3,.3), las=1, ylab='Indirect Effect via Abundance',
	xlab='Direct Effect on Richness', cex.axis=1, cex.lab=1)
usr=par('usr')
polygon(c(usr[1],0,usr[1],usr[2],0,usr[2],usr[1]),
	c(usr[3],0,usr[4],usr[4],0,usr[3],usr[3]), col='grey80')
abline(h=0,v=0)
arrows(use_d$std.ci.lower, use_i$std.all,
	use_d$std.ci.upper, use_i$std.all,
	code=3, angle=90, lwd=2, length=0.05)
arrows(use_d$std.all, use_i$std.ci.lower,
	use_d$std.all, use_i$std.ci.upper,
	code=3, angle=90, lwd=2, length=0.05)
points(use_d$std.all, use_i$std.all, 
	pch=mypch[ifelse(use_vars$type=='env',1,2)], bg=mycol[ifelse(use_vars$mode=='het',1,2)], 
	lwd=2, col='black', cex=2)

#legend('topright',c('Heterogeneity','Optimality'), pch=mypch, pt.bg=mycol, 
#	pt.lwd=2, bg='white', box.lwd=1, pt.cex=2)

#text(-.09,.28,'Significant Effect\non Abundance', font=2, adj=1)
#sigvars_a = names(which(apply(allsp_da[use_vars,c('std.ci.lower','std.ci.upper')], 1, prod)>0))
#text(-.09, allsp_da[sigvars_a,'std.all'], labels=varnames[sigvars_a, 'midName'],
#	adj=1)

#text(.09,.28,'Significant Effect\non Richness', font=2, adj=0)
#sigvars_r = names(which(apply(allsp_d[use_vars,c('std.ci.lower','std.ci.upper')],1,prod)>0))
#text(0.09, allsp_da[sigvars_r,'std.all'], labels=varnames[sigvars_r,'midName'], adj=0)

dev.off()

# Plot legend
svg('./Figures/direct indirect legend.svg', height=3, width=3)
par(mar=c(0,0,0,0))
plot.new()
points(c(0,.1), rep(.1,2), pch=mypch, bg=mycol[2], cex=2, lwd=2)
points(.1, 0, pch=mypch[2], bg=mycol[1], cex=2, lwd=2)
text(rep(.15,2), c(.1,0), labels = c('Optimality','Heterogeneity'), adj=0)
text(c(0,.1), rep(0.15,2), labels=c('Climate','Forest Structure'), srt=45, adj=0)
dev.off()

compare_id = data.frame(direct = use_d$std.all, indirect = use_i$std.all)
rownames(compare_id) = rownames(use_d)
compare_id$sigD  = apply(use_d[,c('std.ci.lower','std.ci.upper')], 1, prod) >0
compare_id$sigI  = apply(use_i[,c('std.ci.lower','std.ci.upper')], 1, prod) >0
compare_id[order(compare_id$direct, decreasing=T),]




##################################################
### Draw path diagram for significant paths

use_pred = subset(predtypes, !(label %in% c('', 'r1')))
use_pred = rbind(use_pred, predtypes['lichen.rich_log',]) # This changes based on which local richness variable is used in the model e.g. fric


## Going to make diagram in two modules: focused on local and regional effects

## Local module:

# Only look at estimates corresponding to paths and correlations between local-regional variables
localvars = rownames(subset(use_pred, scale=='local'))
paths = subset(allsp_ests, op=='~'&(rhs %in% localvars)&(lhs %in% localvars)) # Change based on response variables of interest
paths = rbind(paths, subset(allsp_ests, label=='r1:R1'))

# Calculate which paths are significant
sigpaths = paths[which(apply(paths[,c('std.ci.lower','std.ci.upper')],1,prod)>0),]

# Histogram of significant paths
hist(sigpaths$std.all)
sigpaths$sigcat = cut(abs(sigpaths$std.all), c(0,.2,.4,.6,.8,1))

# Find all variables in the model
allvars = unique(c(paths$rhs,paths$lhs))
allvars[allvars=='regS'] <-'reg' # Change based on response variables of interest
allvars[allvars=='tot_abun_log'] <- 'abun_log' # Change based on response variables of interest

# Replace varnames in sigpaths to generic
sigpaths$lhs[sigpaths$lhs=='tot_abun_log']<-'abun_log'
sigpaths$lhs[sigpaths$lhs=='regS']<-'reg'
sigpaths$rhs[sigpaths$rhs=='tot_abun_log']<-'abun_log'
sigpaths$rhs[sigpaths$rhs=='regS']<-'reg'

## Create a matrix of variable locations
var_locs = matrix(0, nrow=length(allvars), ncol=2, byrow=T)
colnames(var_locs) = c('X','Y')
rownames(var_locs) = allvars

unit = 1

# Start with richness and abundance in center
var_locs['abun_log',] = c(-1,0) 
var_locs['lichen.rich_log',] = c(1,0)

# Add pollution and regS to corners
#var_locs['totalNS',] = c(-2.5,2.5)
var_locs['reg',] = c(0,-2.5) 

# Add climate and local environment vars to top and bottom
c_vars = rownames(use_pred[grep('cm', use_pred$label),])
var_locs[c_vars,] = cbind(3/(length(c_vars)-1)*(0:(length(c_vars)-1))-1.5,rep(3,length(c_vars)))

# Add forest vars to sides
fh_vars = rownames(use_pred[grep('fh', use_pred$label),])
fm_vars = rownames(use_pred[grep('fm', use_pred$label),])

var_locs[fh_vars,]=cbind(rep(3,length(fh_vars)), 3.5/(length(fh_vars)-1)*(0:(length(fh_vars)-1))-2)
var_locs[fm_vars,]=cbind(rep(-3,length(fm_vars)), 3.5/(length(fm_vars)-1)*(0:(length(fm_vars)-1))-2)

var_locs=data.frame(var_locs)

# Make plot
mycex=1.5

svg('./Figures/path diagram nopol model local module.svg', height=5, width=8.5) # Change based on response variable of interest
par(mar=c(0,0,0,0))
par(lend="butt")
plot(var_locs, xlim=c(-6,6), ylim=c(-2.5,6), type='n', axes=F, xlab='', ylab='')
for(i in order(abs(sigpaths$std.all))){
	from = var_locs[sigpaths[i,'rhs'],]
	to = var_locs[sigpaths[i,'lhs'],]
	thickness = as.numeric(sigpaths[i,'sigcat'])*mycex
	col=c('red','black')[(sigpaths[i,'std.all']>0)+1]

	if(as.numeric(sigpaths[i,'sigcat'])>1) arrows(from$X, from$Y, to$X, to$Y, length=0.15, lwd=thickness, col=col)
}
text(var_locs[fm_vars,], labels=varnames[fm_vars,'midName'], pos=2, offset=.25)
text(var_locs[fh_vars,], labels=varnames[fh_vars,'midName'], pos=4, offset=.25)
text(var_locs[c_vars,], labels=varnames[c_vars,'midName'],adj=-.1, srt=45)
text(var_locs[c('lichen.rich_log','abun_log'),], labels=c('Local\nrichness','Abundance'),pos=1, offset=1)
text(var_locs['totalNS',], labels=varnames['totalNS','midName'], pos=2, offset=.25)
text(var_locs['reg',], labels='Regional richness', pos=1, offset=.25)

#legend(x=-5.5, y=5, xjust=0, yjust=1, c('0.0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1.0'), lwd=mycex*(1:5), 
#	bty='n')
#legend(x=-6, y=5, xjust=0, yjust=1, rep('', 5), lwd=mycex*(1:5), col='red', bty='n')

legend(x=-5.5, y=5, xjust=0, yjust=1, c('0.2-0.4','0.4-0.6','0.6-0.8','0.8-1.0'), lwd=mycex*(2:5), 
	bty='n')
legend(x=-6, y=5, xjust=0, yjust=1, rep('', 4), lwd=mycex*(2:5), col='red', bty='n')
text(-5.5, 5, '+', cex=2, adj=-1)
text(-6, 5, '-', cex=2, adj=-2.3)

dev.off()


## Regional module:

# Only look at estimates corresponding to paths and correlations between local-regional variables
regvars = rownames(subset(use_pred, scale=='regional'))
paths = subset(allsp_ests, op=='~'&(rhs %in% regvars)&(lhs %in% regvars)) # paths between regional variables
cors = subset(allsp_ests, op=='~~'&(lhs %in% c_vars)&(rhs %in% regvars)) # correlations between local and regional variables
cors = rbind(cors, subset(allsp_ests, label %in%c('p1:P1','fh4:FH1')))

# Calculate which paths are significant
sigpaths = paths[which(apply(paths[,c('std.ci.lower','std.ci.upper')],1,prod)>0),]
sigcors = cors[which(apply(cors[,c('std.ci.lower','std.ci.upper')],1,prod)>0),]

# Histogram of significant paths
hist(sigpaths$std.all)
sigpaths$sigcat = cut(abs(sigpaths$std.all), seq(0,2,.2))
sigcors$sigcat = cut(abs(sigcors$std.all), seq(0,2,.2))

# Find all variables in the model
allvars = unique(c(paths$rhs,paths$lhs,cors$lhs))
allvars[allvars=='regS'] <-'reg' # Change based on response variables of interest

# Replace varnames in sigpaths to generic
sigpaths$lhs[sigpaths$lhs=='regS']<-'reg'
sigpaths$rhs[sigpaths$rhs=='regS']<-'reg'

## Create a matrix of variable locations
var_locs = matrix(0, nrow=length(allvars), ncol=2, byrow=T)
colnames(var_locs) = c('X','Y')
rownames(var_locs) = allvars

unit = 1

# Start with regional richness in center
var_locs['reg',] = c(0,-1) 

# Add regional pollution and forest to lower corners
#var_locs['totalNS_reg',] = c(1.5,-1.5)
var_locs['regS_tree',] = c(-1.5,-1.5)
#var_locs['totalNS',] = c(2.5,-2.5)
var_locs['PIE.ba.tree',] = c(-2.5,-2.5)

# Add regional climate vars diagonally to top
CH_vars = rownames(use_pred[grep('CH', use_pred$label),])
CM_vars = rownames(use_pred[grep('CM', use_pred$label),])

xpos = seq(1, 3, length.out=length(CH_vars))
slope = 0.75; yint = 1.5
var_locs[CH_vars,] = cbind(xpos, (-slope)*xpos+yint)
var_locs[CM_vars,] = cbind(-xpos, slope*(-xpos)+yint)

# Add local climate vars to top
var_locs[c_vars,]= cbind(3/(length(c_vars)-1)*(0:(length(c_vars)-1))-1.5,rep(2.5,length(c_vars)))

var_locs=data.frame(var_locs)

# Make plot
mycex=.75

svg('./Figures/nopol model large paths regional module.svg', height=5, width=8.5) # Change based on response variable of interest
par(mar=c(0,0,0,0))
par(lend="butt")
plot(var_locs, xlim=c(-7.1,7.1), ylim=c(-2.5,5), type='n', axes=F, xlab='', ylab='')
for(i in order(abs(sigpaths$std.all))){
	from = var_locs[sigpaths[i,'rhs'],]
	to = var_locs[sigpaths[i,'lhs'],]
	thickness = as.numeric(sigpaths[i,'sigcat'])*mycex
	col=c('red','black')[(sigpaths[i,'std.all']>0)+1]

	if(abs(sigpaths[i,'std.all'])>0.2) arrows(from$X, from$Y, to$X, to$Y, length=0.15, lwd=thickness, col=col)
}
for(i in order(abs(sigcors$std.all))){
	from = var_locs[sigcors[i,'rhs'],]
	to = var_locs[sigcors[i,'lhs'],]
	thickness = as.numeric(sigcors[i,'sigcat'])*mycex
	col=c('red','black')[(sigcors[i,'std.all']>0)+1]

	if(abs(sigcors[i,'std.all'])>0.2) arrows(from$X, from$Y, to$X, to$Y, length=0.15, lwd=thickness, col=col, code=3)
}
text(var_locs[c_vars,], labels=varnames[c_vars,'midName'], adj=-.1, srt=45)
text(3, var_locs[CH_vars,'Y'], labels=varnames[CH_vars,'midName'], pos=4, offset=.25)
text(-3, var_locs[CM_vars,'Y'], labels=varnames[CM_vars,'midName'],pos=2, offset=.25)
#text(var_locs['totalNS',], labels=varnames['totalNS','midName'], pos=2, offset=.25)
#text(var_locs['totalNS_reg',], labels=varnames['totalNS_reg','midName'], pos=2, offset=.6)
text(var_locs['regS_tree',], labels=varnames['regS_tree','midName'], pos=2, offset=.6)
text(var_locs['PIE.ba.tree',], labels=varnames['PIE.ba.tree','midName'], pos=2, offset=.25)
text(var_locs['reg',], labels='Regional\nrichness', pos=1, offset=.6)

#legend(x=-5.5, y=4, xjust=0, yjust=1, c('0.0-0.4','0.4-0.8','0.8-1.2','1.2-1.6','1.6-2.0'), 
#	lwd=mycex*seq(2,10,2), bty='n')
#legend(x=-6, y=4, xjust=0, yjust=1, rep('', 5), lwd=mycex*seq(2,10,2), col='red', bty='n')
legend(x=-5.5, y=4, xjust=0, yjust=1, c('0.2-0.4','0.4-0.8','0.8-1.2','1.2-1.6','1.6-2.0'), 
	lwd=mycex*seq(2,10,2), bty='n')
legend(x=-6, y=4, xjust=0, yjust=1, rep('', 5), lwd=mycex*seq(2,10,2), col='red', bty='n')
text(-5.5, 4, '+', cex=2, adj=-1)
text(-6, 4, '-', cex=2, adj=-2.3)

dev.off()


## Regional climate modules: one for each variable
mycex=1.5
unit=2
var_locs = matrix(unit*c(0.25,-1,.5,0,-.5,0,0,1,-.5,-1), nrow=5, ncol=2, byrow=T)
colnames(var_locs) = c('X','Y')
rownames(var_locs) = c('lichen.rich_log','regS','lvar','rvar','tot_abun_log')
var_locs = data.frame(var_locs)

# Model with direct paths from regional climate to local richness
use_ests = reg2rich_ests

regS_eff = subset(use_ests, lhs=='lichen.rich_log'&rhs=='regS')[,c('std.ci.lower','std.all','std.ci.upper')]
abun_eff = subset(use_ests, lhs=='lichen.rich_log'&rhs=='tot_abun_log')[,c('std.ci.lower','std.all','std.ci.upper')]
use_vars = c('mat','iso','pseas','wetness','rain_lowRH')
use_breaks = seq(0,2,length.out=11)

svg('./Figures/No Measurement Error/regional climate paths nopol regTorich model.svg', height=14, width=7)
par(mfrow=c(5,1))
par(mar=c(0,0,0,0))
for(i in use_vars){

	lvar = i
	rvar = paste(i, 'reg_mean', sep='_')

	d_eff = subset(use_ests, lhs=='lichen.rich_log'&rhs==rvar)[,c('std.ci.lower','std.all','std.ci.upper')]
	rl_cor = subset(use_ests, lhs==lvar&rhs==rvar)[,c('std.ci.lower','std.all','std.ci.upper')]
	l_eff = subset(use_ests, lhs=='lichen.rich_log'&rhs==lvar)[,c('std.ci.lower','std.all','std.ci.upper')]
	r_eff = subset(use_ests, lhs=='regS'&rhs==rvar)[,c('std.ci.lower','std.all','std.ci.upper')]
	la_eff = subset(use_ests, lhs=='tot_abun_log'&rhs==lvar)[,c('std.ci.lower','std.all','std.ci.upper')]

	effs = rbind(d_eff, l_eff, r_eff, regS_eff, rl_cor, la_eff, abun_eff)
	effs$from = c('rvar','lvar','rvar','regS','rvar','lvar','tot_abun_log')
	effs$to = c('lichen.rich_log','lichen.rich_log','regS','lichen.rich_log','lvar','tot_abun_log','lichen.rich_log')
	effs$code = c(2,2,2,2,3,2,2)
	effs$col = c('#FF0000','#000000')[as.numeric(effs$std.all>0)+1]
	effs$lty = c(3,1)[as.numeric(apply(effs[,c('std.ci.lower','std.ci.upper')], 1, prod)>0)+1]
	
	#effs$lowCat = as.numeric(cut(abs(effs$std.ci.lower),use_breaks))
	#effs$hiCat = as.numeric(cut(abs(effs$std.ci.upper),use_breaks))
	effs$estCat = as.numeric(cut(abs(effs$std.all),use_breaks))
	
	plot(var_locs, xlim=c(-(unit+1),unit+1), ylim=c(-(unit+1),unit+1), type='n', axes=F, xlab='', ylab='')
	for(a in 1:nrow(effs)){
		arrows(var_locs[effs$from[a],'X'], var_locs[effs$from[a],'Y'],
			var_locs[effs$to[a],'X'], var_locs[effs$to[a],'Y'], col=effs$col[a], 
			lwd=mycex*effs$estCat[a], code=effs$code[a], length=0.15, lty=effs$lty[a])
		midpoint = apply(var_locs[c(effs$from[a], effs$to[a]),], 2, mean)
		midlab = format(effs[a,'std.all'], digits=2)
		confint = paste('(',format(effs[a,'std.ci.lower'], digits=2), ',',
			format(effs[a, 'std.ci.upper'], digits=2), ')', sep='')
	
		points(midpoint['X'], midpoint['Y'], cex=10, col='white', pch=15)
		text(midpoint['X'], midpoint['Y'], paste(midlab, confint, sep='\n'))
	}
	text(var_locs['rvar',], labels=varnames[rvar,'midName'], pos=3, offset=.5)
	text(var_locs['lvar',], labels=varnames[lvar,'midName'], pos=2, offset=.5)
	text(var_locs['lichen.rich_log',], labels='Local richness', pos=1, offset=.5)
	text(var_locs['regS',], labels='Regional richness', pos=4, offset=.5)
	text(var_locs['tot_abun_log',], labels='Lichen abundance', pos=1, offset=.5)

}
dev.off()

# Model with direct paths from regional climate to local richness
use_ests = allsp_ests

regS_eff = subset(use_ests, lhs=='lichen.rich_log'&rhs=='regS')[,c('std.ci.lower','std.all','std.ci.upper')]
abun_eff = subset(use_ests, lhs=='lichen.rich_log'&rhs=='tot_abun_log')[,c('std.ci.lower','std.all','std.ci.upper')]
use_vars = c('mat','iso','pseas','wetness','rain_lowRH')
use_breaks = seq(0,2,length.out=11)

svg('./Figures/No Measurement Error/regional climate paths nopol model.svg', height=14, width=7)
par(mfrow=c(5,1))
par(mar=c(0,0,0,0))
for(i in use_vars){

	lvar = i
	rvar = paste(i, 'reg_mean', sep='_')

	rl_cor = subset(use_ests, lhs==lvar&rhs==rvar)[,c('std.ci.lower','std.all','std.ci.upper')]
	l_eff = subset(use_ests, lhs=='lichen.rich_log'&rhs==lvar)[,c('std.ci.lower','std.all','std.ci.upper')]
	r_eff = subset(use_ests, lhs=='regS'&rhs==rvar)[,c('std.ci.lower','std.all','std.ci.upper')]
	la_eff = subset(use_ests, lhs=='tot_abun_log'&rhs==lvar)[,c('std.ci.lower','std.all','std.ci.upper')]

	effs = rbind(l_eff, r_eff, regS_eff, rl_cor, la_eff, abun_eff)
	effs$from = c('lvar','rvar','regS','rvar','lvar','tot_abun_log')
	effs$to = c('lichen.rich_log','regS','lichen.rich_log','lvar','tot_abun_log','lichen.rich_log')
	effs$code = c(2,2,2,3,2,2)
	effs$col = c('#FF0000','#000000')[as.numeric(effs$std.all>0)+1]
	effs$lty = c(3,1)[as.numeric(apply(effs[,c('std.ci.lower','std.ci.upper')], 1, prod)>0)+1]
	
	#effs$lowCat = as.numeric(cut(abs(effs$std.ci.lower),use_breaks))
	#effs$hiCat = as.numeric(cut(abs(effs$std.ci.upper),use_breaks))
	effs$estCat = as.numeric(cut(abs(effs$std.all),use_breaks))
	
	plot(var_locs, xlim=c(-(unit+1),unit+1), ylim=c(-(unit+1),unit+1), type='n', axes=F, xlab='', ylab='')
	for(a in 1:nrow(effs)){
		arrows(var_locs[effs$from[a],'X'], var_locs[effs$from[a],'Y'],
			var_locs[effs$to[a],'X'], var_locs[effs$to[a],'Y'], col=effs$col[a], 
			lwd=mycex*effs$estCat[a], code=effs$code[a], length=0.15, lty=effs$lty[a])
		midpoint = apply(var_locs[c(effs$from[a], effs$to[a]),], 2, mean)
		midlab = format(effs[a,'std.all'], digits=2)
		confint = paste('(',format(effs[a,'std.ci.lower'], digits=2), ',',
			format(effs[a, 'std.ci.upper'], digits=2), ')', sep='')
	
		points(midpoint['X'], midpoint['Y'], cex=10, col='white', pch=15)
		text(midpoint['X'], midpoint['Y'], paste(midlab, confint, sep='\n'))
	}
	text(var_locs['rvar',], labels=varnames[rvar,'midName'], pos=3, offset=.5)
	text(var_locs['lvar',], labels=varnames[lvar,'midName'], pos=2, offset=.5)
	text(var_locs['lichen.rich_log',], labels='Local richness', pos=1, offset=.5)
	text(var_locs['regS',], labels='Regional richness', pos=4, offset=.5)
	text(var_locs['tot_abun_log',], labels='Lichen abundance', pos=1, offset=.5)

}
dev.off()



## Only show variables with significant paths to richness
## Variables to remove have to be identified by hand

## THIS DOES NOT REFLECT MOST RECENT ANALYSES
# Variables to remove
subset(sigpaths, rhs=='bark_moist_pct.rao.ba') # Check variables that are difficult to distinguish on path diagram
remove_vars = c('diamDiversity','bigTrees','LogSeed.ba','bark_moist_pct.ba','lightDist.mean','LogSeed.rao.ba','wood_SG.rao.ba')
keep_vars = sigvars[!(sigvars %in% remove_vars)]

sigpaths2 = subset(sigpaths, (lhs %in% keep_vars)&(rhs %in% keep_vars))

svg('./Figures/full model significant effects on richness.svg', height=5, width=8.5)
par(mar=c(0,0,0,0))
par(lend="butt")
plot(var_locs, xlim=c(-6,6), ylim=c(-2.5,6), type='n', axes=F, xlab='', ylab='')
for(i in order(abs(sigpaths2$std.all))){
	from = var_locs[sigpaths2[i,'rhs'],]
	to = var_locs[sigpaths2[i,'lhs'],]
	thickness = as.numeric(sigpaths2[i,'sigcat'])*mycex
	col=c('red','black')[(sigpaths2[i,'std.all']>0)+1]

	arrows(from$X, from$Y, to$X, to$Y, length=0.15, lwd=thickness, col=col)
}
text(var_locs[fm_vars,], labels=varnames[fm_vars,'midName'], pos=2, offset=.25)
text(var_locs[fh_vars,], labels=varnames[fh_vars,'midName'], pos=4, offset=.25)
text(var_locs[c_vars,], labels=varnames[c_vars,'midName'],adj=-.1, srt=45)
text(var_locs[c('lichen_rich','tot_abun_log'),], labels=c('Local\nrichness','Abundance'),pos=1, offset=1)
text(var_locs['radiation',], labels='Solar radiation', pos=1, offset=.25)
text(var_locs['totalNS',], labels=varnames['totalNS','midName'], pos=2, offset=.25)
text(var_locs['regS',], labels=varnames['regS','midName'], pos=4, offset=.25)

legend(x=-5.5, y=5, xjust=0, yjust=1, c('0.0-0.2','0.2-0.4','0.4-0.6','0.6-0.8'), lwd=mycex*(1:4), 
	bty='n')
legend(x=-6, y=5, xjust=0, yjust=1, rep('', 4), lwd=mycex*(1:4), col='red', bty='n')
text(-5.5, 5, '+', cex=2, adj=-1)
text(-6, 5, '-', cex=2, adj=-2.3)
dev.off()


##################################################
## Functional Richness

## A table comparing direct effects of local variables on species richness vs functional richness
reg2_d$sig = apply(reg2_d[,c('std.ci.lower','std.ci.upper')], 1, prod)>0
fric_d$sig = apply(fric_d[,c('std.ci.lower','std.ci.upper')], 1, prod)>0
use_cols = c('std.all','std.se','sig')
compare_tab = cbind(predictor=varnames[rownames(reg2_d),'midName'],
	reg2_d[,use_cols], fric_d[rownames(reg2_d),use_cols])
compare_tab[order(compare_tab[,2], decreasing=T),]

use_vars = subset(predtypes, scale=='local' & label!='' & type%in%c('env','forest','abundance'))
use_vars = use_vars[-grep('soil', rownames(use_vars)),]

use_reg2 = reg2_d[rownames(use_vars), ]
use_fric = fric_d[rownames(use_vars),]

# order by largest difference
ordered_vars = rownames(use_fric)[order(abs(use_fric$std.all) - abs(use_reg2$std.all))]



order(abs(use_fric$std.all - use_reg2$std.all)*ifelse(abs(use_fric$std.all)>abs(use_reg2$std.all), 1, -1))]
use_reg2 = use_reg2[ordered_vars,]
use_fric = use_fric[ordered_vars,]

myrange = c(-.3,.9)
myshade = c('55','99','')
names(myshade) = c('het','opt','')

# Comparing direct effects on richness and fric
svg('./Figures/No Measurement Error/Standardized direct effects on Fric vs richness regTorich nopol.svg', height=13, width=16)
dotplot(as.numeric(factor(rownames(use_reg2), levels = ordered_vars))~std.all, data=use_reg2, 
	xlab=list('Standardized Direct Effect',cex=3), ylab='',
	main='',cex.lab=3,aspect=9/10, xlim=myrange,
	panel=function(x,y){
	
		# Add horizontal boxes
		shading = mycols[predtypes[ordered_vars, 'mode'],'local']
		
		panel.rect(myrange[1]-0.1,1:length(ordered_vars)-.5, myrange[2]+0.1, 1:length(ordered_vars)+.5,
			col=shading, border='grey50')
	
		# Add horizontal boxes
		#panel.rect(-2,1:length(ordered_vars)-.5, 2, 1:length(ordered_vars)+.5,
		#	col='white', border='grey50')

		# Add vertical line at 0
		panel.abline(v=0, col='grey30', lty=2, lwd=2)		
		
		# Add null distribution segments for fric
		panel.segments(use_fric$std.ci.lower, y+myadj,
			use_fric$std.ci.upper, y+myadj, 
			col='black', lwd=4.5, lend=1)
		# Add points for direct estimated effects
		panel.points(use_fric$std.all, y+myadj, col='black', fill=mypcols[2], pch=mypch[2], cex=3, lwd=3) 
		
		# Add null distribution segments for richness
		panel.segments(use_reg2$std.ci.lower, y,
			use_reg2$std.ci.upper, y, 
			col='black', lwd=4.5, lend=1)
		# Add points for total estimated effects
		panel.points(x, y, col='black', fill=mypcols[1], pch=mypch[1], cex=3, lwd=3) 
	
		# Add text labeling the variable type
		#vartypes =  sapply(predtypes[ordered_vars,'label'], function(x) toupper(substr(x, 1, 1)))
		#panel.text(-.24, y, labels=mytypes[vartypes], cex=2)
		
	},
	scales=list(y=list(labels=varnames[ordered_vars,'midName'], 
		cex=3, col='black'),
		x=list(cex=3, tick.number=8)),
	key=list(x=.5, y=1, corner=c(.5,0), lines=list(type='o', pch=mypch, fill=mypcols, lwd=3, pt.lwd=3),
		text=list(c('Species Richness','Functional Richness')),
		background='white', cex=3, divide=1, padding.text=5, border=F, columns=2)
)
dev.off()


order(abs(use_fric$std.all) - abs(use_reg2$std.all))


######################################################################
### Effect of removing abundance from model


## Plot effects of No Abundance model vs Direct Regional Paths model (Reg2Rich)

# Set to use types of effects interested in
use_x = subset(reg2_ir, IEvar1=='regS')
use_y = subset(noabun_ir, IEvar1=='regS')

rownames(use_x) = use_x$predictor
rownames(use_y) = use_y$predictor

# Make sure tables are in same order
rownames(use_x) == rownames(use_y)
use_x = use_x[rownames(use_y),]


# Examine differences in effects (in order)
diff = abs(use_x$std.all - use_y$std.all)
names(diff) = rownames(use_x)
diff = diff[order(diff, decreasing=T)]
diff_tab = data.frame(diff, reg2 = use_x[names(diff), 'std.all'], noabun= use_y[names(diff),'std.all'])
diff_tab$signChange = diff_tab$reg2*diff_tab$noabun<0

use_vars = rownames(subset(diff_tab, diff > 0.1))
#use_vars = rownames(subset(diff_tab, signChange))
use_x = use_x[use_vars,]
use_y = use_y[use_vars,]

use_pch = c(22,23,21)[factor(predtypes[use_vars,'mode'], levels=c('het','opt',''))]
use_col = apply(predtypes[use_vars,], 1, function(x){
	if(x['mode']==''){ 'grey30' } else { mycols[x['mode'],x['scale']] }
})
use_lim = range(c(use_x[,c('std.ci.lower','std.ci.upper')],use_y[,c('std.ci.lower','std.ci.upper')]))#+c(-.05, .05)

svg('./Figures/No Measurement Error/compare reg2 vs noabun indirect effects via regS.svg', height=5, width=5)
par(mar=c(4,4,1,1))
plot(use_x$std.all, use_y$std.all, type='n', 
	xlim=use_lim, ylim=use_lim, las=1, ylab='Indirect Effect (No Abundance)',
	xlab='Indirect Effect (Direct Regional Paths)', cex.axis=1, cex.lab=1)
abline(h=0,v=0)
abline(0,1)
abline(0,-1)
arrows(use_x$std.ci.lower, use_y$std.all,
	use_x$std.ci.upper, use_y$std.all,
	code=3, angle=90, lwd=2, length=0.05)
arrows(use_x$std.all, use_y$std.ci.lower,
	use_x$std.all, use_y$std.ci.upper,
	code=3, angle=90, lwd=2, length=0.05)
points(use_x$std.all, use_y$std.all, bg=use_col, pch=use_pch, lwd=2, cex=1.5)#

dev.off()




## Plot total effects for reg2 and no abun

# Order variables from lowest to highest total effects
total = reg2 # This changes based on what response variable is being analyzed
total = rbind(total, reg2_d[c('regS','tot_abun_log'),]) #

total_noabun = noabun
total_noabun = rbind(total_noabun, noabun_d['regS',]) #

ordered_vars = rownames(total[order(total$std.all),])

# Put tables in same order
use_total = total[ordered_vars,]
use_noabun = total_noabun[ordered_vars,]

# Define range limits that will include 95% confidence intervals
myrange = range(c(use_total[,c('std.ci.lower','std.ci.upper')], use_noabun[,c('std.ci.lower','std.ci.upper')]), na.rm=T)+c(-.04, .04)
jitter = 0.05

mypcols = c('white','grey30')
mypch = c(22,23)

# Make plot
svg('./Figures/No Measurement Error/Standardized total effects on AllSp richness nopol regTorich and noabun.svg', height=20, width=22)
dotplot(as.numeric(factor(rownames(use_total), levels = ordered_vars))~std.all, data=use_total, 
	xlab=list('Standardized Effect',cex=3), ylab='',
	main='',cex.lab=3,aspect=5/3, xlim=myrange,
	panel=function(x,y){
	
		# Add horizontal boxes
		colorcombos = predtypes[ordered_vars, c('mode','scale')]
		colorcombos['regS','mode'] = 'opt'
		colororder = apply(colorcombos, 1, function(x) mycols[x[1],x[2]])

		panel.rect(myrange[1]-0.01,1:length(ordered_vars)-.5, myrange[2]+0.01, 1:length(ordered_vars)+.5,
			col=colororder, border='grey50')
	
		# Add vertical line at 0
		panel.abline(v=c(-.2,0,.2), col='grey30', lty=2, lwd=2)

		# Add null distribution segments for total effects of noabun model
		panel.segments(use_noabun$std.ci.lower, y+jitter,
			use_noabun$std.ci.upper, y+jitter, 
			col='black', lwd=4.5, lend=1)
		# Add points for total estimated effects of noabun model
		panel.points(use_noabun$std.all, y+jitter, col='black', fill=mypcols[2], pch=mypch[2], cex=3, lwd=3) 
	
		# Add null distribution segments for total effects or reg2 model
		panel.segments(use_total$std.ci.lower, y-jitter,
			use_total$std.ci.upper, y-jitter, 
			col='black', lwd=4.5, lend=1)
		# Add points for total estimated effects of reg2 model
		panel.points(x, y-jitter, col='black', fill=mypcols[1], pch=mypch[1], cex=3, lwd=3) 
	
},
	scales=list(y=list(labels=varnames[ordered_vars,'midName'], 
		cex=3, col='black'),
		x=list(cex=3, tick.number=8)),
	key=list(x=.5, y=1, corner=c(.5,0), lines=list(type='o', pch=mypch, fill=mypcols, lwd=3, pt.lwd=3),
		text=list(c('Direct Regional Paths model','No Abundance model')),
		background='white', cex=3, divide=1, padding.text=5, border=F, columns=2)
)
dev.off()












#############################################################################
### Soil Models Comparison






##########################################################################
### Spatial Autocorrelation

## Need to read in and fit models from cluster scripts


# Predict values

reg2_pred = predict(regTorich_nopol_fit) #This doesn't work because of a bug in lavaan 0.5-16













