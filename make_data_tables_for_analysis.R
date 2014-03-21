## This script compiles data tables for FIA lichen plots into:
## master : a data table with all variables from 2228 plots
## model_data : data from plots used in models from all plots without any data missing, includes PCA variables
## trans_data : environmental data log or sqrt transformed to reduce skew
## working_data_unstd :  trans_data scaled by linear factor (usually 10) to put variables on similar range
## working_data : trans_data scaled to have mean 0 and std dev 1
## Outliers are analyzed and removed from all data sets and working data sets are divided into equal sized test ('_test') and fitting ('_fit') sets.


setwd('C://Users/jrcoyle/Documents/UNC/Projects/FIA Lichen')
options(stringsAsFactors=F)
varnames=read.csv('varnames.csv', row.names=1)
source('./GitHub/FIA-Lichens/fia_lichen_analysis_functions.R')

########################################################################
### Combine data files into a master data file and working data file ###

# FIA plot data
plot_locs = read.csv('./Data/fia_lichen_plot_locations.csv')
county_data = read.csv('./Data/fia_lichen_county_data.csv')
plot_data = merge(plot_locs, county_data[,c('COUNTYNM','state.abbr','STATE','COUNTY','POP2000','POP2010','SQMI')])

# Lichen richness data
rich_current = read.csv('./Data/lichen_richness_current.csv')
rich_legacy = read.csv('./Data/lichen_richness_legacy.csv')
rich_current = rich_current[,names(rich_legacy)]
rich_data = rbind(rich_current, rich_legacy)

# Lichen abundance data
abun_data = read.csv('./Data/lichen abundance based on tree occupancy.csv') # Not available for all plots b/c originally calculated after subsetting.

# FIA tree data
tree_data = read.csv('./Data/TreeData/master_data_forest.csv')[,c('yrplot.id','S.tree','D.abun.tree',
	'D.area.tree','maxDiam','numTreesBig','numTreesSm','propDead','numDead','numCut',
	'PIE.stems.tree','PIE.ba.tree','wood_moist_pct.rao.pres','bark_moist_pct.rao.pres',
	'wood_SG.rao.pres','bark_SG.rao.pres','LogSeed.rao.pres','wood_moist_pct.rao.stems',
	'bark_moist_pct.rao.stems','wood_SG.rao.stems','bark_SG.rao.stems','LogSeed.rao.stems',
	'wood_moist_pct.stems','bark_moist_pct.stems','wood_SG.stems','bark_SG.stems','LogSeed.stems',
	'wood_moist_pct.rao.ba','bark_moist_pct.rao.ba','wood_SG.rao.ba','bark_SG.rao.ba','LogSeed.rao.ba',
	'wood_moist_pct.ba','bark_moist_pct.ba','wood_SG.ba','bark_SG.ba','LogSeed.ba','diamDist.mean',
	'n.stems','basal.area','light.mean','light.var','lightDist.mean','totalCirc')]
tree_pca = read.csv('./Data/TreeData/tree_funcgrp_pca1-3.csv')

# Regional tree richness
regS_tree = read.csv('./Data/TreeData/Regional tree diversity/fia_lichen_tree_regS.csv')

# Lichen functional diversity data
fd_data = read.csv('./Data/LichenTraits/fia_lichen_LIAS_means_diversity.csv')

# Lichen regional richness data
reg_data = read.csv('./Data/Regional Richness/fia_lichen_reg_richness.csv')

# Environmental data
env_data = read.csv('./Data/fia_lichen_env_data_points.csv')
env_plot_data = read.csv('./Data/fia_lichen_env_data_plots.csv')
env_reg_data = read.csv('./Data/fia_lichen_env_data_regional.csv')

# Merge to make master data file
master = merge(plot_data, env_plot_data, all.x=T, all.y=F)
master = merge(master, rich_data, all.x=T, all.y=F)
master = merge(master, abun_data, all.x=T, all.y=F)
master = merge(master, tree_data, all.x=T, all.y=F)
master = merge(master, tree_pca, all.x=T, all.y=F)
master = merge(master, regS_tree, all.x=T, all.y=F)
master = merge(master, fd_data, all.x=T, all.y=F)
master = merge(master, reg_data, all.x=T, all.y=F)
master = merge(master, env_data, all.x=T, all.y=F)
master = merge(master, env_reg_data, all.x=T, all.y=F)

# Save data
write.csv(master, './Data/fia_lichen_master_data.csv', row.names=F)

rownames(master) = master$yrplot.id
###############################################################################
### Data Subsetting ###
master = read.csv('./Data/fia_lichen_master_data.csv', row.names=1)

# Use recent plots after plot design had been standardized
model_data = subset(master, MEASYEAR>=1997)

# Remove plots with only one large tree (heterogeneity measurements are NA)
model_data = subset(model_data, numTreesBig>1) # removes 44 plots widely distributed across US

## Define variables potentially used in analysis
predictors = read.csv('predictors.csv')

# Subset by predictors that are included in model_data (not derived PCs)
measured_pred = subset(predictors, pred %in% colnames(model_data))

# Plot histograms of predictors
pdf('./Figures/Predictor variable histograms.pdf', height=6, width=6)
for(p in measured_pred$pred){
	hist(model_data[,p], main=varnames[p,'displayName'])
	mtext(paste('# Missing =',sum(is.na(model_data[,p]))), side=3, line=0, adj=1)
}
dev.off()

# Remove records that are missing data in these variables
model_data = model_data[rowSums(is.na(model_data[,measured_pred$pred]))==0,] # 1927 plots

## Test for correlations among variables

# Standardize variables
use_response = 'lichen.rich'
use_pred = measured_pred$pred
use_data = model_data[,c(use_response,use_pred)]

# Transform skewed variables (except proportions)
logTrans_vars = c('totalCirc', 'ap','ap_reg_var','pseas_reg_var','rh_reg_var',
	'wetness_reg_var','rain_lowRH_reg_var')
sqrtTrans_vars = c('bark_SG.rao.ba','bark_moist_pct.rao.ba','wood_SG.rao.ba',
	'wood_moist_pct.rao.ba','LogSeed.rao.ba', 'S.tree')
for(v in logTrans_vars){
	use_data[,v] = log10(use_data[,v])
}
for(v in sqrtTrans_vars){
	use_data[,v] = sqrt(use_data[,v])
}
# Center and scale data
use_data[,2:ncol(use_data)] = scale(use_data[,2:ncol(use_data)] )

# Pairwise correlations
library('corrplot')

# Record correlations between variables to determine whether to get rid of some.
cortab = cor(use_data, use='complete.obs')
which(cortab>0.7&cortab<1, arr.ind=T) 
write.csv(cortab, 'correlation matrix std vars.csv', row.names=T)

cortabsig = 1-abs(cortab)

useorder = c(measured_pred$pred[order(measured_pred$type, measured_pred$scale, measured_pred$mode)], 'lichen.rich')
cortab = cortab[useorder,useorder]
cortabsig = cortabsig[useorder,useorder]

png('./Figures/correlation matrix std vars.png', height=1200, width=1200, type='cairo')
corrplot(cortab[2:nrow(cortab),2:ncol(cortab)], method='square', type='upper', diag=F, 
	order='original', hclust.method='complete', p.mat=cortabsig[2:nrow(cortab),2:ncol(cortab)],
	sig.level=.6, insig='blank', tl.cex=1.5, tl.col=1, cl.cex=2, mar=c(1,1,4,1))
dev.off()

## Define new PCA variable pairs for intrinsically correlated variables
newvars = data.frame(yrplot.id=rownames(model_data))

# max tree size and tree size range
diam_pca = prcomp(na.omit(use_data[,c('diamDist.mean','maxDiam')]))
diam_vars = data.frame(predict(diam_pca))
names(diam_vars) = c('bigTrees','diamDiversity')
diam_vars$diamDiversity = -1*diam_vars$diamDiversity
diam_vars$yrplot.id = rownames(diam_vars)
newvars = merge(newvars, diam_vars, all.x=T)

## Create data set with variables used for modeling
myvars = c('lichen.rich','Parmeliaceae','Physciaceae','fric','fdiv','raoQ','wetness','rain_lowRH',
	'mat','iso','pseas','totalNS','radiation','wetness_reg_mean','rain_lowRH_reg_mean',
	'mat_reg_mean','iso_reg_mean','pseas_reg_mean','wetness_reg_var','rain_lowRH_reg_var',
	'mat_reg_var','iso_reg_var','pseas_reg_var','totalNS_reg','regS_tree',
	'bark_moist_pct.ba','bark_moist_pct.rao.ba','wood_SG.ba','wood_SG.rao.ba','PC1',
	'LogSeed.ba','LogSeed.rao.ba','PIE.ba.tree','propDead','light.mean','lightDist.mean',
	'regS','regParm','regPhys','tot_abun_log','parm_abun_log','phys_abun_log'
)

model_data = cbind(model_data[,myvars], newvars[,2:ncol(newvars)])

# Subset predictor table by variables to be used in subsequent models
model_pred = subset(predictors, pred %in% colnames(model_data))

## Create scaled and transformed datasets

trans_data = model_data
logTrans_vars = c('pseas_reg_var','wetness_reg_var','rain_lowRH_reg_var')
sqrtTrans_vars = c('bark_moist_pct.rao.ba','wood_SG.rao.ba', 'LogSeed.rao.ba')

for(v in logTrans_vars){
	trans_data[,v] = log10(trans_data[,v])
}
for(v in sqrtTrans_vars){
	trans_data[,v] = sqrt(trans_data[,v])
}

hist(trans_data$bigTrees) # Not much I can do about transforming this, so I won't

working_data = trans_data

# Make transformation of richness response used in models
working_data$lichen.rich_log = log(working_data$lichen.rich+1)
working_data$Parm_log = log(working_data$Parmeliaceae+1)
working_data$Phys_log = log(working_data$Physciaceae+1)

# Rescale by mean and stddev for standardized data
# Note: this scales the response variable (lichen richness), which may not be what we want to do
working_data = data.frame(scale(working_data, center=T, scale=T))

# Plot correlation matrix of rescaled and transformed data
cortab = cor(working_data[,model_pred$pred], use='complete.obs')
cortabsig = 1-abs(cortab)

png('./Figures/correlation matrix working vars.png', height=900, width=900, type='cairo')
corrplot(cortab, method='square', type='upper', diag=F, 
	order='hclust', hclust.method='complete', p.mat=cortabsig,
	sig.level=0.6, insig='blank', tl.cex=1.5, tl.col=1, cl.cex=2, mar=c(1,1,4,1))
dev.off()

## Determine which points are outliers and remove.

outliers = working_data[,model_pred$pred]
outliers[,]<-NA
for(v in model_pred$pred){
	ols = lm(working_data[,'lichen.rich_log']~working_data[,v])
	cd = cooks.distance(ols)
	outs = which(cd >4/nrow(working_data))

	outliers[outs,v]<-cd[outs]
}
outliers = outliers[apply(outliers, 1, function(x) sum(!is.na(x))>0),]
outliers = data.frame(outliers)
outliers$numOut = apply(outliers, 1, function(x) sum(!is.na(x)))
write.csv(outliers, 'cooks D outliers.csv', row.names=T)

# All 35 plots with 1 species are outliers
rownames(subset(model_data, lichen.rich==1)) %in% rownames(outliers)[which(outliers$numOut>12)]
outliers[rownames(subset(model_data, lichen.rich==1)),]

# Take out all plots with 1-2 species from outliers so that they will be ignored when assessing other outliers
outliers = subset(outliers, rownames(outliers) %in% rownames(model_data[model_data$lichen.rich>2,]))

# 18 plots with more that 2 species are outliers in 1/4 of the predictor variables.
dim(subset(model_data, lichen.rich>2&rownames(model_data) %in% rownames(outliers)[which(outliers$numOut>8)]))

# Used to check outliers in each variable
i=model_pred$pred[33]
ols = lm(working_data$lichen.rich_log~working_data[,i])
opar <- par(mfrow = c(2, 2), oma = c(0, 0, 1.1, 0))
plot(ols, las = 1)
cd = cooks.distance(ols)
which(cd > 4/nrow(working_data))
outliers[order(outliers[,i], decreasing=T),][1:20,]
subset(model_data, rownames(model_data) %in% names(which(cd>0.01)))

# mat - none
# iso - none
# pseas - none
# radiation - none
# totalNS - none
# bark_moist_pct.ba - none
# bark_moist_pct.rao.ba: 1998_17_43_6379
# wood_SG - none
# wood_SG.rao.ba - none
# LogSeed.ba : 1999_41_25_7306
# LogSeed.rao.ba - none
# PIE.ba.tree - none
# propDead - none
# light.mean - none: 1997_56_7_6475, 2007_4_9_87353 are outliers with 2-3 trees
# lightDist.mean : 2007_4_19_83376
# totalCirc - none
# regS - none
# regParm - none
# regPhys - none
# tot_abun_log : several appear to be outlier b/c there is one species
# parm_abun_log - none
# phys_abun_log - none
# bigTrees - none
# diamDiversity : 2004_16_49_85627, two trees diams 5.5, 27.1 leads to large diameter difference for relatively small trees, PCA exacerbates this
# wetness - none
# rain_lowRH - none
# PC1 - none
# mat_reg_mean, iso_reg_mean, pseas_reg_mean, wetness_reg_mean, rain_lowRH_reg_mean - none
# mat_reg_var, iso_reg_var, pseas_reg_var, wetness_reg_var, rain_lowRH_reg_var - none\
# regS_tree - none

# Remove outliers
remove_plots = c('2004_16_49_85627','2007_4_19_83376','1999_41_25_7306','1998_17_43_6379')
working_data = subset(working_data, !(rownames(working_data) %in% remove_plots))
model_data = subset(model_data, !(rownames(model_data) %in% remove_plots))
trans_data = subset(trans_data, !(rownames(trans_data) %in% remove_plots))


## Save data sets
write.csv(working_data, './Data/fia_lichen_working_data.csv', row.names=T)
write.csv(trans_data, './Data/fia_lichen_trans_data.csv', row.names=T)
write.csv(model_data, './Data/fia_lichen_model_data.csv', row.names=T)

## Divide data into fitting and testing data sets
allplots = rownames(model_data)
usedata = master[allplots, c('state.abbr', 'yrplot.id')]

# Only run once !!! 
#unlist(tapply(usedata$yrplot.id, usedata$state.abbr, function(x){
#	n = ifelse(runif(1)>0.5, ceiling(length(x)/2), floor(length(x)/2))
#	sample(x, n)
#}))->fitplots
#length(fitplots) # stopped at 961
#testplots = allplots[!(allplots %in% fitplots)]

## Write out list of test and fit plots
#write.csv(data.frame(yrplot.id=testplots), './Data/model test plots.csv')
#write.csv(data.frame(yrplot.id=fitplots), './Data/model fit plots.csv')



################# OLD CODE ###############

# Make a list of predictors of each type
climate = c('ap','mat','iso','pseas','rh')
local_env = c('radiation')
pollution = c('totalNS')
forest_het = c('bark_SG.rao.ba', 'bark_moist_pct.rao.ba', 'wood_SG.rao.ba', 'wood_moist_pct.rao.ba', 
	'LogSeed.rao.ba','lightDist.mean','PIE.ba.tree','S.tree', 'propDead', 'diamDist.mean')
forest_hab = c('bark_SG.ba', 'bark_moist_pct.ba', 'wood_SG.ba', 'wood_moist_pct.ba', 
	'LogSeed.ba','light.mean','totalCirc', 'PC1')
forest_time = c('maxDiam')
region = c('regS')

predictors = data.frame(pred = c(climate, pollution, local_env, forest_het, forest_hab, forest_time, region),
	type = c(rep('climate',length(climate)+length(pollution)+length(local_env)),rep('forest',length(forest_het)+length(forest_hab)+length(forest_time)), rep('region',length(region))),
	shape = c(rep(2,length(climate)), rep(1, length(pollution)), rep(2, length(local_env)), rep(1, length(forest_het)), rep(2, length(forest_hab)), rep(1, length(forest_time)), rep(1,length(region))),
	hyp = c(rep('resource',length(climate)+length(pollution)+length(local_env)), rep('niche', length(forest_het)), rep('resource', length(forest_hab)), rep('time', length(forest_time)), rep('region', length(region)))
)
predictors[predictors$pred=='totalCirc','shape']<-1


#working_data_unstd = trans_data

# For unstd data: rescale by constant so that variances in path analysis will be of similar scale
#working_data_unstd$mat = working_data_unstd$mat/10
#working_data_unstd$pseas = working_data_unstd$pseas/10
#working_data_unstd$radiation = working_data_unstd$radiation/1000000
#working_data_unstd$totalNS = working_data_unstd$totalNS/100
#working_data_unstd$bark_moist_pct.ba = working_data_unstd$bark_moist_pct.ba/10
#working_data_unstd$bark_moist_pct.rao.ba = working_data_unstd$bark_moist_pct.rao.ba*10
#working_data_unstd$wood_SG.rao.ba = working_data_unstd$wood_SG.rao.ba*10
#working_data_unstd$wood_SG.ba = working_data_unstd$wood_SG.ba*10
#working_data_unstd$LogSeed.rao.ba = working_data_unstd$LogSeed.rao.ba*10
#working_data_unstd$PIE.ba.tree = working_data_unstd$PIE.ba.tree*10
#working_data_unstd$propDead = working_data_unstd$propDead*10
#working_data_unstd$light.mean = working_data_unstd$light.mean/10
#working_data_unstd$lightDist.mean = working_data_unstd$lightDist.mean/10
#working_data_unstd$regS = working_data_unstd$regS/10
#working_data_unstd$regParm = working_data_unstd$regParm/10
#working_data_unstd$regPhys = working_data_unstd$regPhys/10
#working_data_unstd$bigTrees = working_data_unstd$bigTrees/10
#working_data_unstd$PC1 = working_data_unstd$PC1*10


#working_data_unstd$lichen.rich_log = log(working_data_unstd$lichen.rich+1)
#working_data_unstd$lichen.rich = working_data_unstd$lichen.rich/10



