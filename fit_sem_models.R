# This script is used to fit SEM models of lichen FIA data

source('./UNC/Projects/FIA Lichen/GitHub/FIA-Lichens/load_data.R')
source('./GitHub/FIA-Lichens/fia_lichen_analysis_functions.R')

save.image('sem_analysis.RData')
load('sem_analysis.RData')

################################################################################
### Path Analysis ###

library(lavaan)
library(semPlot)
library(semTools)


# Read in table of predictor variable types
predtypes = read.csv('predictors.csv', row.names=1)

# Full path model with measurement error model on lichen.rich_log
#3/24/2014: Decided to add a bunch of regional scale predictors 
#1/3/2014: Decided to use standardized data because the model will not fit to unscaled data and I did not want to apply an arbitrary scaling.
#           Standardized data is in working_data dataframe. Note that the response variables are also standardized.     
# I played with using log richness and adding abundace with this model before fitting the final model below.

## Define models 
# Rules for coefficient names:
# dependent variable first, independent variable second
# Capital = regional, lowercase = local
# R = richness, A = abundance, F = forest, C = climate/environment, P = pollution
# H = heterogeneity, M = optimality (only used for F and C)
# number indicates which variable in the category: e.g. fh3 is LogSeed.rao.ba
# labels are not actually used- I relabel stuff later when creating the summary datatables from bootstrapped parameter estimates.

path_finalmod = "

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

	# Regional climate, pollution, and regional forest heterogeneity effects on regional richness
	regS ~ R1CM1*wetness_reg_mean + R1CM2*rain_lowRH_reg_mean + R1CM3*iso_reg_mean + R1CM4*pseas_reg_mean + R1CM5*mat_reg_mean +
		R1CH1*wetness_reg_var + R1CH2*rain_lowRH_reg_var + R1CH3*iso_reg_var + R1CH4*pseas_reg_var + R1CH5*mat_reg_var +
		R1FH1*regS_tree + R1P1*totalNS_reg

	# Regional climate effects on regional forest heterogeneity
	regS_tree ~ FH1CM1*wetness_reg_mean + FH1CM2*rain_lowRH_reg_mean + FH1CM3*iso_reg_mean + FH1CM4*pseas_reg_mean + FH1CM5*mat_reg_mean +
		FH1CH1*wetness_reg_var + FH1CH2*rain_lowRH_reg_var + FH1CH3*iso_reg_var + FH1CH4*pseas_reg_var + FH1CH5*mat_reg_var

	# Local climate effects on pollution
	totalNS ~ p1cm1*wetness + p1cm2*rain_lowRH

	# Regional climate effects on pollution
	totalNS_reg ~ P1CM1*wetness_reg_mean + P1CM2*rain_lowRH_reg_mean

	# Effects on local lichen richness
	lichen_rich ~ r1fh1*bark_moist_pct.rao.ba + r1fh2*wood_SG.rao.ba + r1fh3*LogSeed.rao.ba +
		r1fh4*PIE.ba.tree + r1fh5*propDead + r1fh6*lightDist.mean + r1fh7*diamDiversity +
		r1fm1*bark_moist_pct.ba + r1fm2*wood_SG.ba + r1fm3*LogSeed.ba + r1fm4*bigTrees + r1fm5*light.mean + r1fm6*PC1 +
		r1cm1*wetness + r1cm2*rain_lowRH + r1cm3*iso + r1cm4*pseas + r1cm5*mat + r1cm6*radiation + r1R1*regS + r1p1*totalNS +
		r1a1*tot_abun_log
	
	# Effects on local lichen abundance
	tot_abun_log ~ a1fh1*bark_moist_pct.rao.ba + a1fh2*wood_SG.rao.ba + a1fh3*LogSeed.rao.ba +
		a1fh4*PIE.ba.tree + a1fh5*propDead + a1fh6*lightDist.mean + a1fh7*diamDiversity +
		a1fm1*bark_moist_pct.ba + a1fm2*wood_SG.ba + a1fm3*LogSeed.ba + a1fm4*bigTrees + a1fm5*light.mean + a1fm6*PC1 +
		a1cm1*wetness + a1cm2*rain_lowRH + a1cm3*iso + a1cm4*pseas + a1cm5*mat + a1cm6*radiation + a1p1*totalNS 

	## Covariances between exogenous and endogenous variables
	# Dont need to specificy for Climate L-R b/c these are exogenous and are calculated automatically in lavaan

	# Local-regional pollution
	totalNS ~~ p1P1*totalNS_reg
	
	# Local-regional forest structure
	regS_tree ~~ FH1fh4*PIE.ba.tree

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

# Fit model
fit_fixed = sem(path_finalmod, data=working_data_fit, fixed.x=T, estimator='ML', se='robust')
fit_free = sem(path_finalmod, data=working_data_fit, fixed.x=F, estimator='ML', se='robust')

ests_fixed = parameterEstimates(fit_fixed, standardized=T)
ests_free = parameterEstimates(fit_free, standardized=T)

summary(fit_fixed, standardized=T, rsq=TRUE, fit.measures=T)
fitMeasures(fit_fixed)

mi1 = modindices(fit_fixed)
mi = miPowerFit(fit_fixed)
problem_parms = subset(mi, decision %in% c('M', 'EPC:M'))
problem_parms[order(problem_parms$mi),]
res = resid(fit_fixed)

library(corrplot)
cortab = res$cov
cortabsig = 1-abs(cortab)

# only plot residual correlations greater than .1
png('./Figures/finalmod residual correlations.png', height=1700, width=1700, type='cairo')
corrplot(cortab, method='square', type='upper', diag=F, 
	order='original', hclust.method='complete', p.mat=cortabsig,
	sig.level=.9, insig='blank', tl.cex=1.5, tl.col=1, cl.cex=2, mar=c(1,1,4,1))
dev.off()

which(abs(res$cov[lower.tri(res$cov)])>.2, arr.ind=T)


## Model without abundance
# See sem_boot_noabun_allsp.R script
path_noabun = "

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

	# Regional climate, pollution, and regional forest heterogeneity effects on regional richness
	regS ~ R1CM1*wetness_reg_mean + R1CM2*rain_lowRH_reg_mean + R1CM3*iso_reg_mean + R1CM4*pseas_reg_mean + R1CM5*mat_reg_mean +
		R1CH1*wetness_reg_var + R1CH2*rain_lowRH_reg_var + R1CH3*iso_reg_var + R1CH4*pseas_reg_var + R1CH5*mat_reg_var +
		R1FH1*regS_tree + R1P1*totalNS_reg

	# Regional climate effects on regional forest heterogeneity
	regS_tree ~ FH1CM1*wetness_reg_mean + FH1CM2*rain_lowRH_reg_mean + FH1CM3*iso_reg_mean + FH1CM4*pseas_reg_mean + FH1CM5*mat_reg_mean +
		FH1CH1*wetness_reg_var + FH1CH2*rain_lowRH_reg_var + FH1CH3*iso_reg_var + FH1CH4*pseas_reg_var + FH1CH5*mat_reg_var

	# Local climate effects on pollution
	totalNS ~ p1cm1*wetness + p1cm2*rain_lowRH

	# Regional climate effects on pollution
	totalNS_reg ~ P1CM1*wetness_reg_mean + P1CM2*rain_lowRH_reg_mean

	# Effects on local lichen richness
	lichen_rich ~ r1fh1*bark_moist_pct.rao.ba + r1fh2*wood_SG.rao.ba + r1fh3*LogSeed.rao.ba +
		r1fh4*PIE.ba.tree + r1fh5*propDead + r1fh6*lightDist.mean + r1fh7*diamDiversity +
		r1fm1*bark_moist_pct.ba + r1fm2*wood_SG.ba + r1fm3*LogSeed.ba + r1fm4*bigTrees + r1fm5*light.mean + r1fm6*PC1 +
		r1cm1*wetness + r1cm2*rain_lowRH + r1cm3*iso + r1cm4*pseas + r1cm5*mat + r1cm6*radiation + r1R1*regS + r1p1*totalNS
		
	## Covariances between exogenous and endogenous variables
	# Dont need to specificy for Climate L-R b/c these are exogenous and are calculated automatically in lavaan

	# Local-regional pollution
	totalNS ~~ p1P1*totalNS_reg
	
	# Local-regional forest structure
	regS_tree ~~ FH1fh4*PIE.ba.tree

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

## Model with paths from regional variables to local richness
path_regTorich = "

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

	# Regional climate, pollution, and regional forest heterogeneity effects on regional richness
	regS ~ R1CM1*wetness_reg_mean + R1CM2*rain_lowRH_reg_mean + R1CM3*iso_reg_mean + R1CM4*pseas_reg_mean + R1CM5*mat_reg_mean +
		R1CH1*wetness_reg_var + R1CH2*rain_lowRH_reg_var + R1CH3*iso_reg_var + R1CH4*pseas_reg_var + R1CH5*mat_reg_var +
		R1FH1*regS_tree + R1P1*totalNS_reg

	# Regional climate effects on regional forest heterogeneity
	regS_tree ~ FH1CM1*wetness_reg_mean + FH1CM2*rain_lowRH_reg_mean + FH1CM3*iso_reg_mean + FH1CM4*pseas_reg_mean + FH1CM5*mat_reg_mean +
		FH1CH1*wetness_reg_var + FH1CH2*rain_lowRH_reg_var + FH1CH3*iso_reg_var + FH1CH4*pseas_reg_var + FH1CH5*mat_reg_var

	# Local climate effects on pollution
	totalNS ~ p1cm1*wetness + p1cm2*rain_lowRH

	# Regional climate effects on pollution
	totalNS_reg ~ P1CM1*wetness_reg_mean + P1CM2*rain_lowRH_reg_mean

	# Effects on local lichen richness
	lichen_rich ~ r1fh1*bark_moist_pct.rao.ba + r1fh2*wood_SG.rao.ba + r1fh3*LogSeed.rao.ba +
		r1fh4*PIE.ba.tree + r1fh5*propDead + r1fh6*lightDist.mean + r1fh7*diamDiversity +
		r1fm1*bark_moist_pct.ba + r1fm2*wood_SG.ba + r1fm3*LogSeed.ba + r1fm4*bigTrees + r1fm5*light.mean + r1fm6*PC1 +
		r1cm1*wetness + r1cm2*rain_lowRH + r1cm3*iso + r1cm4*pseas + r1cm5*mat + r1cm6*radiation + r1R1*regS + r1p1*totalNS +
		r1CM1*wetness_reg_mean + r1CM3*rain_lowRH_reg_mean + r1CM3*iso_reg_mean + r1CM4*pseas_reg_mean + r1CM5*mat_reg_mean +
		r1CH1*wetness_reg_var + r1CH2*rain_lowRH_reg_var + r1CH3*iso_reg_var + r1CH4*pseas_reg_var + r1CH5*mat_reg_var +
		r1P1*totalNS_reg + r1FH1*regS_tree + r1a1*tot_abun_log
	
	# Effects on local lichen abundance
	tot_abun_log ~ a1fh1*bark_moist_pct.rao.ba + a1fh2*wood_SG.rao.ba + a1fh3*LogSeed.rao.ba +
		a1fh4*PIE.ba.tree + a1fh5*propDead + a1fh6*lightDist.mean + a1fh7*diamDiversity +
		a1fm1*bark_moist_pct.ba + a1fm2*wood_SG.ba + a1fm3*LogSeed.ba + a1fm4*bigTrees + a1fm5*light.mean + a1fm6*PC1 +
		a1cm1*wetness + a1cm2*rain_lowRH + a1cm3*iso + a1cm4*pseas + a1cm5*mat + a1cm6*radiation + a1p1*totalNS 

	## Covariances between exogenous and endogenous variables
	# Dont need to specificy for Climate L-R b/c these are exogenous and are calculated automatically in lavaan

	# Local-regional pollution
	totalNS ~~ p1P1*totalNS_reg
	
	# Local-regional forest structure
	regS_tree ~~ FH1fh4*PIE.ba.tree

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

regTorich_fit =  sem(path_regTorich, data=working_data_fit, fixed.x=T, estimator='ML', se='robust.sem')
summary(regTorich_fit, standardized=T, rsq=TRUE, fit.measures=T)

# Examine model residuals
regTorich_res = resid(regTorich_fit)
which(abs(regTorich_res$cov)>0.2, arr.ind=T)

# Examine significance of model paths
regTorich_ests = parameterEstimates(regTorich_fit, standardized=T, ci=T, level=0.95)
subset(regTorich_ests , pvalue>0.05)
regTorich_paths = subset(regTorich_ests , op=='~')
regTorich_paths[order(regTorich_paths$std.all, decreasing=T),]

# Write out a data set to be used on the cluster
write.csv(working_data[testplots$yrplot.id,], './SEM models/standardized_test_dataset.csv', row.names=T)

## Parameter estimates for models were bootstrapped using scripts run on the Kure comuting cluster
# see: sem_boot_regTorich_allsp.R, sem_boot_finalmod_allsp.R, sem_boot_noabun_allsp.R

path_nopol = "

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

	# Regional climate, pollution, and regional forest heterogeneity effects on regional richness
	regS ~ R1CM1*wetness_reg_mean + R1CM2*rain_lowRH_reg_mean + R1CM3*iso_reg_mean + R1CM4*pseas_reg_mean + R1CM5*mat_reg_mean +
		R1CH1*wetness_reg_var + R1CH2*rain_lowRH_reg_var + R1CH3*iso_reg_var + R1CH4*pseas_reg_var + R1CH5*mat_reg_var +
		R1FH1*regS_tree

	# Regional climate effects on regional forest heterogeneity
	regS_tree ~ FH1CM1*wetness_reg_mean + FH1CM2*rain_lowRH_reg_mean + FH1CM3*iso_reg_mean + FH1CM4*pseas_reg_mean + FH1CM5*mat_reg_mean +
		FH1CH1*wetness_reg_var + FH1CH2*rain_lowRH_reg_var + FH1CH3*iso_reg_var + FH1CH4*pseas_reg_var + FH1CH5*mat_reg_var

	# Effects on local lichen richness
	lichen_rich ~ r1fh1*bark_moist_pct.rao.ba + r1fh2*wood_SG.rao.ba + r1fh3*LogSeed.rao.ba +
		r1fh4*PIE.ba.tree + r1fh5*propDead + r1fh6*lightDist.mean + r1fh7*diamDiversity +
		r1fm1*bark_moist_pct.ba + r1fm2*wood_SG.ba + r1fm3*LogSeed.ba + r1fm4*bigTrees + r1fm5*light.mean + r1fm6*PC1 +
		r1cm1*wetness + r1cm2*rain_lowRH + r1cm3*iso + r1cm4*pseas + r1cm5*mat + r1cm6*radiation + r1R1*regS + r1a1*tot_abun_log
	
	# Effects on local lichen abundance
	tot_abun_log ~ a1fh1*bark_moist_pct.rao.ba + a1fh2*wood_SG.rao.ba + a1fh3*LogSeed.rao.ba +
		a1fh4*PIE.ba.tree + a1fh5*propDead + a1fh6*lightDist.mean + a1fh7*diamDiversity +
		a1fm1*bark_moist_pct.ba + a1fm2*wood_SG.ba + a1fm3*LogSeed.ba + a1fm4*bigTrees + a1fm5*light.mean + a1fm6*PC1 +
		a1cm1*wetness + a1cm2*rain_lowRH + a1cm3*iso + a1cm4*pseas + a1cm5*mat + a1cm6*radiation + a1p1*totalNS 

	## Covariances between exogenous and endogenous variables
	# Dont need to specificy for Climate L-R b/c these are exogenous and are calculated automatically in lavaan

	# Local-regional forest structure
	regS_tree ~~ FH1fh4*PIE.ba.tree

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



##################################################################
## Read in tables of parameter estimates and effects
# All parameter estimares
allsp_ests = read.csv('./SEM models/No Pollution/nopol_AllSp_testdata_parameterEstimates.csv')
#parm_ests = read.csv('./SEM models/finalmod_Parm_testdata_parameterEstimates.csv')
#phys_ests = read.csv('./SEM models/finalmod_Phys_testdata_parameterEstimates.csv')
#raoq_ests = read.csv('./SEM models/finalmod_RaoQ_testdata_parameterEstimates.csv')
noabun_ests = read.csv('./SEM models/No Pollution/noabun_nopol_AllSp_testdata_parameterEstimates.csv')
reg2rich_ests = read.csv('./SEM models/No Pollution/regTorich_nopol_AllSp_testdata_parameterEstimates.csv')

# Total effects
allsp = read.csv('./SEM models/No Pollution/nopol_AllSp_testdata_totaleffects.csv', row.names=1)
#parm = read.csv('./SEM models/finalmod_Parm_testdata_totaleffects.csv', row.names=1)
#phys = read.csv('./SEM models/finalmod_Phys_testdata_totaleffects.csv', row.names=1)
#raoq = read.csv('./SEM models/finalmod_RaoQ_testdata_totaleffects.csv', row.names=1)
#fric = read.csv('./SEM models/finalmod_Fric_testdata_totaleffects.csv', row.names=1)
reg2 = read.csv('./SEM models/No Pollution/regTorich_nopol_AllSp_testdata_totaleffects.csv', row.names=1)
noabun = read.csv('./SEM models/No Pollution/noabun_nopol_AllSp_testdata_totaleffects.csv', row.names=1)

# Direct effects
allsp_d = read.csv('./SEM models/No Pollution/nopol_AllSp_testdata_directeffects_richness.csv', row.names=1)
#parm_d = read.csv('./SEM models/finalmod_Parm_testdata_directeffects_richness.csv', row.names=1)
#phys_d = read.csv('./SEM models/finalmod_Phys_testdata_directeffects_richness.csv', row.names=1)
#fric_d = read.csv('./SEM models/finalmod_Fric_testdata_directeffects_richness.csv', row.names=1)
#raoq_d = read.csv('./SEM models/finalmod_RaoQ_testdata_directeffects_richness.csv', row.names=1)
reg2_d = read.csv('./SEM models/No Pollution/regTorich_nopol_AllSp_testdata_directeffects_richness.csv', row.names=1)
noabun_d = read.csv('./SEM models/No Pollution/noabun_nopol_AllSp_testdata_directeffects_richness.csv', row.names=1)

# Direct effect on abundance
allsp_da = read.csv('./SEM models/No Pollution/nopol_AllSp_testdata_directeffects_abundance.csv', row.names=1)
#parm_da = read.csv('./SEM models/finalmod_Parm_testdata_directeffects_abundance.csv', row.names=1)
#phys_da = read.csv('./SEM models/finalmod_Phys_testdata_directeffects_abundance.csv', row.names=1)

# Direct effects on regional richness
allsp_dr = read.csv('./SEM models/No Pollution/nopol_AllSp_testdata_directeffects_regS.csv', row.names=1)
reg2_dr = read.csv('./SEM models/No Pollution/regTorich_nopol_AllSp_testdata_directeffects_regS.csv', row.names=1)

# Indirect effects via abundance
allsp_i = read.csv('./SEM models/No Pollution/nopol_AllSp_testdata_indirecteffects_via_abundance.csv', row.names=1)
#parm_i = read.csv('./SEM models/finalmod_Parm_testdata_indirecteffects_via_abundance.csv', row.names=1)
#phys_i = read.csv('./SEM models/finalmod_Phys_testdata_indirecteffects_via_abundance.csv', row.names=1)
#fric_i = read.csv('./SEM models/finalmod_Fric_testdata_indirecteffects_via_abundance.csv', row.names=1)
#raoq_i = read.csv('./SEM models/finalmod_RaoQ_testdata_indirecteffects_via_abundance.csv', row.names=1)
reg2_i = read.csv('./SEM models/No Pollution/regTorich_nopol_AllSp_testdata_indirecteffects_via_abundance.csv', row.names=1)

# Indirect effects of regional scale predictors
allsp_ir = read.csv('./SEM models/No Pollution/nopol_AllSp_testdata_regionalvars_indirecteffects.csv')
#parm_ir = read.csv('./SEM models/finalmod_Parm_testdata_regionalvars_indirecteffects.csv')
#phys_ir = read.csv('./SEM models/finalmod_Phys_testdata_regionalvars_indirecteffects.csv')
#fric_ir = read.csv('./SEM models/finalmod_Fric_testdata_regionalvars_indirecteffects.csv')
#raoq_ir = read.csv('./SEM models/finalmod_RaoQ_testdata_regionalvars_indirecteffects.csv')
reg2_ir = read.csv('./SEM models/No Pollution/regTorich_nopol_AllSp_testdata_regionalvars_indirecteffects.csv')

# Indirect effects via forest structure
allsp_if = read.csv('./SEM models/No Pollution/nopol_AllSp_testdata_indirecteffects_via_forest.csv')
#parm_if = read.csv('./SEM models/finalmod_Parm_testdata_indirecteffects_via_forest.csv')
#phys_if = read.csv('./SEM models/finalmod_Phys_testdata_indirecteffects_via_forest.csv')
#fric_if = read.csv('./SEM models/finalmod_Fric_testdata_indirecteffects_via_forest.csv')
#raoq_if = read.csv('./SEM models/finalmod_RaoQ_testdata_indirecteffects_via_forest.csv')

# Indirect effects via local pollution
#allsp_ip = read.csv('./SEM models/finalmod_AllSp_testdata_indirecteffects_via_pollution.csv')
#parm_ip = read.csv('./SEM models/finalmod_Parm_testdata_indirecteffects_via_pollution.csv')
#phys_ip = read.csv('./SEM models/finalmod_Phys_testdata_indirecteffects_via_pollution.csv')
#fric_ip = read.csv('./SEM models/finalmod_Fric_testdata_indirecteffects_via_pollution.csv')
#raoq_ip = read.csv('./SEM models/finalmod_RaoQ_testdata_indirecteffects_via_pollution.csv')

##################################################################
### Results ###

# How many predictors have significant total effects?
names(which(apply(allsp[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0))
names(which(apply(parm[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0))
names(which(apply(phys[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0))
names(which(apply(reg2[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0))

tot_sig = allsp[which(apply(allsp[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0),]
tot_sig = tot_sig[order(abs(tot_sig$std.all)),]

tot_sig = reg2[which(apply(reg2[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0),]
tot_sig = tot_sig[order(abs(tot_sig$std.all)),]

# Compare significance of total effects across taxa
tot_sig = data.frame(AllSp = apply(allsp[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0,
	Parm = apply(parm[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0,
	Phys = apply(phys[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0)
data.frame(AllSp = allsp$std.all,
	Parm = parm$std.all,
	Phys = phys$std.all)

# Compare significance of direct effects across taxa
dir_sig = data.frame(AllSp = apply(allsp_d[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0,
	Parm = apply(parm_d[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0,
	Phys = apply(phys_d[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0)
colSums(dir_sig)

allsp_d[order(abs(allsp_d$std.all)),]

# Examine order of effects among forest structure variables
allsp_f = subset(allsp, type %in% c('FH','FM'))
allsp_df = subset(allsp_d, type %in% c('FH','FM'))

allsp_f[order(abs(allsp_f$std.all)),]
allsp_df[order(abs(allsp_df$std.all)),]


parm_f = subset(parm, type %in% c('FH','FM'))
parm_df = subset(parm_d, type %in% c('FH','FM'))
parm_f[order(abs(parm_f$std.all)),]
parm_df[order(abs(parm_df$std.all)),]

phys_f = subset(phys, type %in% c('FH','FM'))
phys_df = subset(phys_d, type %in% c('FH','FM'))
phys_f[order(abs(phys_f$std.all)),]
phys_df[order(abs(phys_df$std.all)),]




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

# Load previously fit models (from Kure)

load('./SEM models/No Pollution/nopol_AllSp_testdata_output.RData')
load('./SEM models/No Pollution/noabun_nopol_AllSp_testdata_output.RData')
load('./SEM models/No Pollution/regTorich_nopol_AllSp_testdata_output.RData')

nopol_fit
noabun_nopol_fit
regTorich_nopol_fit

AIC(nopol_fit, noabun_nopol_fit, regTorich_nopol_fit)
anova(regTorich_nopol_fit, nopol_fit) # including direct paths to local richness makes the model significantly better.

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

# Assess model fit
examine_mod = nopol_fit 

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

subset(resdf, abs(residual)>0.1)
subset(resdf, var1=='lichen.rich_log')
subset(resdf, var1=='regS')


# Residual correlation figure: only plot residual correlations greater than .1
cortab = res
cortabsig = 1-abs(cortab)
png('./Figures/nopol noabun residual correlations.png', height=1700, width=1700, type='cairo')
corrplot(cortab, method='square', type='upper', diag=F, 
	order='original', hclust.method='complete', p.mat=cortabsig,
	sig.level=.9, insig='blank', tl.cex=1.5, tl.col=1, cl.cex=2, mar=c(1,1,4,1))
dev.off()

# What are the residuals? Obs - Exp covariances
obscov = inspectSampleCov(path_noabun_nopol, data=working_data_test)$cov #Make sure to change to test or fit data where appropriate
modcov = fitted(examine_mod)$cov 
covdiff =  obscov - modcov 

svg('./Figures/nopol noabun covariance residuals.svg', height=6, width=6)
plot(obscov[covdiff!=0], modcov[covdiff!=0], xlab='Observed Covariances', 
	ylab='Model Covariances', las=1)
abline(0,1)
abline(h=0,v=0)
points(obscov['lichen.rich_log',], modcov['lichen.rich_log',], pch=16, col=2, cex=1.5)
points(obscov['regS',], modcov['regS',], pch=16, col='blue', cex=1.1)
legend('topleft', c('Local richness','Regional richness'), col=c('red','blue'), 
	pch=16, pt.cex=c(1.5,1.1), bty='n')
#points(obscov[c('totalNS','totalNS_reg'),], modcov[c('totalNS','totalNS_reg'),], pch=16, col='green', cex=1)
#legend('bottomright', 'Pollution', col='green', pch=16, bty='n')
mtext('No Abundance Model', 3, 0.5) 
dev.off()

cbind(obscov['lichen.rich_log',], modcov['lichen.rich_log',])

## Compare model predictions to observed data.
# Model covariance and mean structure
mu = fitted(examine_mod)$mean
Sigma= matrix(fitted(examine_mod)$cov, nrow=length(mu), ncol=length(mu))

# Random variables that follow model structure
nvars = length(mu); nobs = nrow(working_data_test)
rvars = matrix(rnorm(nvars*nobs, mean=mu), nrow=nvars, ncol=nobs)
L = chol(Sigma)
predvars = t(L) %*% matrix(rnorm(nvars*nobs), nrow=nvars, ncol=nobs)
predvars = data.frame(t(predvars)); colnames(predvars) = names(mu)
 
## Fit model to randomly generated data. What are the fit statistics?
random_fit = sem(path_noabun_nopol, data=predvars, fixed.x=T, estimator='ML', se='robust.sem')
summary(random_fit, standardized=T, rsq=TRUE, fit.measures=T)
obsrand = inspectSampleCov(path_noabun_nopol, data=predvars)$cov
modrand = fitted(random_fit)$cov 

svg('./Figures/compare noabun_nopol random vs model residuals.svg', height=6, width=12)
par(mfrow=c(1,2))
plot(obsrand[covdiff!=0], modrand[covdiff!=0], xlab='Observed Covariances', 
	ylab='Model Covariances', las=1)
abline(0,1)
abline(h=0,v=0)
points(obsrand['lichen.rich_log',], modrand['lichen.rich_log',], pch=16, col=2, cex=1.5)
points(obsrand['regS',], modrand['regS',], pch=16, col='blue', cex=1.1)
legend('topleft', c('Local richness','Regional richness'), col=c('red','blue'), 
	pch=16, pt.cex=c(1.5,1.1), bty='n')
mtext('No Abundance Model: Random Data', 3, 0.5) 

plot(obscov[covdiff!=0], modcov[covdiff!=0], xlab='Observed Covariances', 
	ylab='Model Covariances', las=1)
abline(0,1)
abline(h=0,v=0)
points(obscov['lichen.rich_log',], modcov['lichen.rich_log',], pch=16, col=2, cex=1.5)
points(obscov['regS',], modcov['regS',], pch=16, col='blue', cex=1.1)
legend('topleft', c('Local richness','Regional richness'), col=c('red','blue'), 
	pch=16, pt.cex=c(1.5,1.1), bty='n')
mtext('No Abundance Model: Observed Data', 3, 0.5) 
dev.off()


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

# Make df of total effects including direct effects of regS and abundance

# Order variables from lowest to highest total effects
total = noabun # This changes based on what response variable is being analyzed
#total = rbind(total, reg2_d[c('regS','tot_abun_log'),])
total = rbind(total, noabun_d['regS',])

ordered_vars = rownames(total[order(total$std.all),])

# Put tables in same order
use_total = total[ordered_vars,]

library(lattice)

# Define range limits that will include 95% confidence intervals
myrange = range(use_total[,c('std.ci.lower','std.ci.upper')], na.rm=T)+c(-.04, .04)
myrange[1] = -1 #-.8

# Color scheme: http://colorschemedesigner.com/#3341SsYrGvyw0
mycols = c('#2415B0','#00BF32')
mycolsbw = c('grey80','white')
names(mycols) = c('regional','local')
names(mycolsbw) = c('regional','local')
mycols_trans = paste(mycols, '50', sep='')
names(mycols_trans) = c('regional','local')
plot(1:length(mycols),1:length(mycols), type='n'); text(1:length(mycols),1:length(mycols),labels=names(mycols), cex=2, col=mycols) # Check colors
mypcols = c('white','grey30')
mypch = c(22,23)
myadj=.15
#mytypes = expression('C'['H'],'C'['O'],'F'['H'],'F'['O'],'P','','') # symbols used in plot to denote variable types
#names(mytypes)=c('CH','CM','FH','FM','P','R','A')
myshade = c('55','99','')
names(myshade) = c('het','opt','')
mytypes = expression('C','F','P','','') # symbols used in plot to denote variable types
names(mytypes)=c('C','F','P','R','A')

# Make plot
svg('./Figures/Standardized total effects on AllSp richness noabun nopol.svg', height=20, width=19)
dotplot(as.numeric(factor(rownames(use_total), levels = ordered_vars))~std.all, data=use_total, 
	xlab=list('Standardized Effect',cex=3), ylab='',
	main='',cex.lab=3,aspect=5/3, xlim=myrange,
	panel=function(x,y){
	
		# Add horizontal boxes
		shading = myshade[predtypes[ordered_vars, 'mode']]
		shading[is.na(shading)]<- '99'
		panel.rect(myrange[1]-0.1,1:length(ordered_vars)-.5, myrange[2]+0.1, 1:length(ordered_vars)+.5,
			col=paste(mycols[predtypes[ordered_vars,'scale']], shading, sep=''), border='grey50')
	
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
		panel.text(myrange[1]-0.05, y, labels=mytypes[vartypes], cex=2)
		
	},
	scales=list(y=list(labels=varnames[ordered_vars,'midName'], 
		cex=3, col='black'),
		x=list(cex=3, tick.number=8))
)
dev.off()

### Direct effects on regional richness

# Define data set to use
direct_reg = allsp_dr

# Order variables from lowest to highest direct effects
ordered_vars = rownames(direct_reg[order(direct_reg$std.all),])

# Put tables in same order
use_df = direct_reg[ordered_vars,]

# Define range limits that will include 95% confidence intervals
myrange = range(use_df[,c('std.ci.lower','std.ci.upper')], na.rm=T)+c(-.04, .04)
myrange[1] = -1.8

myshade = c('55','99','')
names(myshade) = c('het','opt','')
mytypes = expression('C','F','P','','') # symbols used in plot to denote variable types
names(mytypes)=c('C','F','P','R','A')

# Make plot
svg('./Figures/Standardized direct effects on AllSp regional richness nopol.svg', height=9, width=19)
dotplot(as.numeric(factor(rownames(use_df), levels = ordered_vars))~std.all, data=use_df, 
	xlab=list('Standardized Effect',cex=3), ylab='',
	main='',cex.lab=3,aspect=4/5, xlim=myrange,
	panel=function(x,y){
	
		# Add horizontal boxes
		panel.rect(myrange[1]-0.1,1:length(ordered_vars)-.5, myrange[2]+0.1, 1:length(ordered_vars)+.5,
			col=paste(mycols[predtypes[ordered_vars,'scale']], myshade[predtypes[ordered_vars, 'mode']], sep=''),
			border='grey50')

		# Add vertical line at 0
		panel.abline(v=c(-.2,0,.2), col='grey30', lty=2, lwd=2)		
		
		# Add null distribution segments for direct effects
		panel.segments(use_df$std.ci.lower, y,
			use_df$std.ci.upper, y, 
			col='black', lwd=4.5, lend=1)
		# Add points for direct estimated effects
		panel.points(x, y, col='black', fill=mypcols[1], pch=mypch[1], cex=3, lwd=3) 
	
		# Add text labeling the variable type
		vartypes =  sapply(predtypes[ordered_vars,'label'], function(x) toupper(substr(x, 1, 1)))
		panel.text(-1.7, y, labels=mytypes[vartypes], cex=2)
		
	},
	scales=list(y=list(labels=varnames[ordered_vars,'midName'], 
		cex=3, col='black'),
		x=list(cex=3, tick.number=8))
)
dev.off()

### Direct and indirect effects on local richness

# Define datasets to use
direct = allsp_d
indirect = allsp_i

direct = subset(direct, rownames(direct)!='regS')
indirect = subset(indirect, rownames(indirect)!='regS')

# Order variables from lowest to highest direct effects
ordered_vars = rownames(direct[order(direct$std.all),])

# Put tables in same order
use_direct = direct[ordered_vars,]
use_indirect = indirect[ordered_vars,]

# Define range limits that will include 95% confidence intervals
myrange = range(c(use_direct[,c('std.ci.lower','std.ci.upper')],use_indirect[,c('std.ci.lower','std.ci.upper')]), na.rm=T)+c(-.04, .04)
myrange[1] = -.38

jitter = 0.15

# Make plot
svg('./Figures/Standardized direct indirect effects on AllSp richness nopol.svg', height=13, width=19)
dotplot(as.numeric(factor(rownames(use_direct), levels = ordered_vars))~std.all, data=use_direct, 
	xlab=list('Standardized Effect',cex=3), ylab='',
	main='',cex.lab=3,aspect=5/4, xlim=myrange,
	panel=function(x,y){
	
		# Add horizontal boxes
		shading = myshade[predtypes[ordered_vars, 'mode']]
		shading[is.na(shading)]<- '99'
		panel.rect(myrange[1]-0.1,1:length(ordered_vars)-.5, myrange[2]+0.1, 1:length(ordered_vars)+.5,
			col=paste(mycols[predtypes[ordered_vars,'scale']], shading, sep=''), border='grey50')

		# Add vertical line at 0
		panel.abline(v=c(-.2,0,.2), col='grey30', lty=2, lwd=2)		
		
		# Add null distribution segments for direct effects
		panel.segments(use_direct$std.ci.lower, y+jitter,
			use_direct$std.ci.upper, y+jitter, 
			col='black', lwd=4.5, lend=1)
		# Add points for direct estimated effects
		panel.points(x, y+jitter, col='black', fill=mypcols[1], pch=mypch[1], cex=3, lwd=3) 

		# Add null distribution segments for indirect effects
		panel.segments(use_indirect$std.ci.lower, y-jitter,
			use_indirect$std.ci.upper, y-jitter, 
			col='black', lwd=4.5, lend=1)
		# Add points for indirect estimated effects
		panel.points(use_indirect$std.all, y-jitter, col='black', fill=mypcols[2], pch=mypch[2], cex=3, lwd=3) 
	
		# Add text labeling the variable type
		vartypes =  sapply(predtypes[ordered_vars,'label'], function(x) toupper(substr(x, 1, 1)))
		panel.text(-.34, y, labels=mytypes[vartypes], cex=2)
		
	},
	scales=list(y=list(labels=varnames[ordered_vars,'midName'], 
		cex=3, col='black'),
		x=list(cex=3, tick.number=8)),
	key=list(x=.5, y=1, corner=c(.5,0), lines=list(type='o', pch=mypch, fill=mypcols, lwd=3, pt.lwd=3),
		text=list(c('Direct Effect','Indirect Effect')),
		background='white', cex=3, divide=1, padding.text=5, border=F, columns=2)
)
dev.off()

###########################################
### Are total effects at regional scale larger than total effects at local scale?
cvars = c('wetness','rain_lowRH','pseas','iso','mat')

lvars = c(cvars, 'PIE.ba.tree')
rvars = c(paste(cvars, 'reg_mean', sep='_'), 'regS_tree')

# keep track of which model is currently saved as 'total'
svg('./Figures/compare local-regional total effects regTorich nopol.svg', height=5.5, width=5.5)
par(mar=c(5,5,1,1))
plot(total[lvars, 'std.all'], total[rvars, 'std.all'], 	
	xlab='Local Effect', ylab='Regional Effect', type='n', las=1,
	xlim=c(-1,1), ylim=c(-1,1))
usr=par('usr')
polygon(c(usr[1],0,usr[1],usr[2],0,usr[2],usr[1]),
	c(usr[3],0,usr[4],usr[4],0,usr[3],usr[3]), col='grey80')
abline(h=0,v=0, lwd=1, col='black')
arrows(total[lvars,'std.ci.lower'], total[rvars,'std.all'],
	total[lvars,'std.ci.upper'], total[rvars,'std.all'], code=3, angle=90, length=0.05)
arrows(total[lvars,'std.all'], total[rvars,'std.ci.lower'],
	total[lvars,'std.all'], total[rvars,'std.ci.upper'], code=3, angle=90, length=0.05)

# Labels for regTorich
text(total[lvars,'std.all'], total[rvars,'std.all'], varnames[lvars,'midName'],
	pos=c(2,2,4,4,2,4), offset=1.3)

# Labels for finalmod
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

write.csv(climReg_IE_tab, 'Indirect regional climate effects nopol.csv', row.names=F)


indirect = phys_i
indirectR = phys_ir
indirectF = phys_if
total = phys
direct = phys_d

use_total_sub = subset(total, type %in% c('C','L'))
use_direct_sub = subset(direct, type %in% c('C','L'))
use_indirect_sub = subset(indirect, type %in% c('C','L'))


use_total_sub = use_total_sub[order_clim,]
use_direct_sub = use_direct_sub[order_clim,]
use_indirect_sub = use_indirect_sub[order_clim,]
use_indirectR_sub = indirectR[order_clim,]
use_indirectF_sub = indirectF[c(paste(order_clim,'FM', sep='_'),paste(order_clim,'FH', sep='_')),]


climEff_tab = data.frame(Predictor=varnames[use_direct_sub$predictor,'displayName'], 
	direct = use_direct_sub$std.all,
	indirectR = use_indirectR_sub$std.all,
	indirectA = use_indirect_sub$std.all,
	indirectFH = subset(use_indirectF_sub, Ftype=='FH')$std.all,
	indirectFM = subset(use_indirectF_sub, Ftype=='FM')$std.all,
	directSig = apply(use_direct_sub[,c('std.ci.lower','std.ci.upper')], 1, prod)>0,
	indirectRSig = apply(use_indirectR_sub[,c('std.ci.lower','std.ci.upper')], 1, prod)>0,
	indirectASig = apply(use_indirect_sub[,c('std.ci.lower','std.ci.upper')], 1, prod)>0,
	indirectFHSig = apply(subset(use_indirectF_sub, Ftype=='FH')[,c('std.ci.lower','std.ci.upper')], 1, prod)>0,
	indirectFMSig = apply(subset(use_indirectF_sub, Ftype=='FM')[,c('std.ci.lower','std.ci.upper')], 1, prod)>0
)

write.csv(climEff_tab, './SEM models/Compare effects climate variables Phys.csv', row.names=F)


################################################
## Plot direct vs. indirect effects via abundance for local predictors
use_vars = rownames(predtypes)[predtypes$scale=='L']
use_vars = use_vars[use_vars!='abun_log']
use_vars = use_vars[-grep('soil', use_vars)]

mypch = c(22,23)
mycol=c('white','grey30')
myfact = ifelse(allsp_d[use_vars,'type']=='FH', 1, 2)

svg('./Figures/compare richness abundance effects forest vars.svg', height=5, width=5 )
par(mar=c(4,4,1,1))
plot(allsp_d[use_vars,'std.all'], allsp_da[use_vars,'std.all'], type='n', 
	xlim=c(-.3, .3), ylim=c(-.3,.3), las=1, ylab='Direct Effect on Abundance',
	xlab='Direct Effect on Richness', cex.axis=1, cex.lab=1)
usr=par('usr')
polygon(c(usr[1],0,usr[1],usr[2],0,usr[2],usr[1]),
	c(usr[3],0,usr[4],usr[4],0,usr[3],usr[3]), col='grey80')
abline(h=0,v=0)
arrows(allsp_d[use_vars,'std.ci.lower'], allsp_da[use_vars,'std.all'],
	allsp_d[use_vars,'std.ci.upper'], allsp_da[use_vars,'std.all'],
	code=3, angle=90, lwd=2, length=0.05)
arrows(allsp_d[use_vars,'std.all'], allsp_da[use_vars,'std.ci.lower'],
	allsp_d[use_vars,'std.all'], allsp_da[use_vars,'std.ci.upper'],
	code=3, angle=90, lwd=2, length=0.05)
points(allsp_d[use_vars,'std.all'], allsp_da[use_vars,'std.all'], 
	pch=mypch[myfact], bg=mycol[myfact], lwd=2, col='black', cex=2)
legend('topright',c('Heterogeneity','Optimality'), pch=mypch, pt.bg=mycol, 
	pt.lwd=2, bg='white', box.lwd=1, pt.cex=2)

#text(-.09,.28,'Significant Effect\non Abundance', font=2, adj=1)
#sigvars_a = names(which(apply(allsp_da[use_vars,c('std.ci.lower','std.ci.upper')], 1, prod)>0))
#text(-.09, allsp_da[sigvars_a,'std.all'], labels=varnames[sigvars_a, 'midName'],
#	adj=1)

#text(.09,.28,'Significant Effect\non Richness', font=2, adj=0)
#sigvars_r = names(which(apply(allsp_d[use_vars,c('std.ci.lower','std.ci.upper')],1,prod)>0))
#text(0.09, allsp_da[sigvars_r,'std.all'], labels=varnames[sigvars_r,'midName'], adj=0)

dev.off()



##################################################
### Draw path diagram for significant paths

use_pred = subset(predtypes, label !='')
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
var_locs['lichen_rich',] = c(1,0)

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

svg('./Figures/nopol model large paths local module.svg', height=5, width=8.5) # Change based on response variable of interest
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
text(var_locs[c('lichen_rich','abun_log'),], labels=c('Local\nrichness','Abundance'),pos=1, offset=1)
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
unit=2
var_locs = matrix(unit*c(0.25,-1,.5,0,-.5,0,0,1,-.5,-1), nrow=5, ncol=2, byrow=T)
colnames(var_locs) = c('X','Y')
rownames(var_locs) = c('lichen_rich','regS','lvar','rvar','tot_abun_log')
var_locs = data.frame(var_locs)

# Model with direct paths from regional climate to local richness
regS_eff = subset(reg2rich_ests, lhs=='lichen_rich'&rhs=='regS')[,c('std.ci.lower','std.all','std.ci.upper')]
abun_eff = subset(reg2rich_ests, lhs=='lichen_rich'&rhs=='tot_abun_log')[,c('std.ci.lower','std.all','std.ci.upper')]
use_vars = c('mat','iso','pseas','wetness','rain_lowRH')
use_breaks = seq(0,2,length.out=11)

svg('./Figures/regional climate paths nopol regTorich model.svg', height=14, width=7)
par(mfrow=c(5,1))
par(mar=c(0,0,0,0))
for(i in use_vars){

	lvar = i
	rvar = paste(i, 'reg_mean', sep='_')

	d_eff = subset(reg2rich_ests, lhs=='lichen_rich'&rhs==rvar)[,c('std.ci.lower','std.all','std.ci.upper')]
	rl_cor = subset(reg2rich_ests, lhs==lvar&rhs==rvar)[,c('std.ci.lower','std.all','std.ci.upper')]
	l_eff = subset(reg2rich_ests, lhs=='lichen_rich'&rhs==lvar)[,c('std.ci.lower','std.all','std.ci.upper')]
	r_eff = subset(reg2rich_ests, lhs=='regS'&rhs==rvar)[,c('std.ci.lower','std.all','std.ci.upper')]
	la_eff = subset(reg2rich_ests, lhs=='tot_abun_log'&rhs==lvar)[,c('std.ci.lower','std.all','std.ci.upper')]

	effs = rbind(d_eff, l_eff, r_eff, regS_eff, rl_cor, la_eff, abun_eff)
	effs$from = c('rvar','lvar','rvar','regS','rvar','lvar','tot_abun_log')
	effs$to = c('lichen_rich','lichen_rich','regS','lichen_rich','lvar','tot_abun_log','lichen_rich')
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
	text(var_locs['lichen_rich',], labels='Local richness', pos=1, offset=.5)
	text(var_locs['regS',], labels='Regional richness', pos=4, offset=.5)
	text(var_locs['tot_abun_log',], labels='Lichen abundance', pos=1, offset=.5)

}
dev.off()

# Model without direct paths from regional climate to local richness
regS_eff = subset(allsp_ests, lhs=='lichen_rich'&rhs=='regS')[,c('std.ci.lower','std.all','std.ci.upper')]
abun_eff = subset(allsp_ests, lhs=='lichen_rich'&rhs=='tot_abun_log')[,c('std.ci.lower','std.all','std.ci.upper')]
use_vars = c('mat','iso','pseas','wetness','rain_lowRH')
use_breaks = seq(0,2,length.out=11)

svg('./Figures/regional climate paths nopol model.svg', height=14, width=7)
par(mfrow=c(5,1))
par(mar=c(0,0,0,0))
for(i in use_vars){

	lvar = i
	rvar = paste(i, 'reg_mean', sep='_')

	rl_cor = subset(allsp_ests, lhs==lvar&rhs==rvar)[,c('std.ci.lower','std.all','std.ci.upper')]
	l_eff = subset(allsp_ests, lhs=='lichen_rich'&rhs==lvar)[,c('std.ci.lower','std.all','std.ci.upper')]
	r_eff = subset(allsp_ests, lhs=='regS'&rhs==rvar)[,c('std.ci.lower','std.all','std.ci.upper')]
	la_eff = subset(reg2rich_ests, lhs=='tot_abun_log'&rhs==lvar)[,c('std.ci.lower','std.all','std.ci.upper')]

	effs = rbind(l_eff, r_eff, regS_eff, rl_cor, la_eff, abun_eff)
	effs$from = c('lvar','rvar','regS','rvar','lvar','tot_abun_log')
	effs$to = c('lichen_rich','regS','lichen_rich','lvar','tot_abun_log','lichen_rich')
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
	text(var_locs['lichen_rich',], labels='Local richness', pos=1, offset=.5)
	text(var_locs['regS',], labels='Regional richness', pos=4, offset=.5)
	text(var_locs['tot_abun_log',], labels='Lichen abundance', pos=1, offset=.5)

}
dev.off()







## Only show variables with significant paths to richness
## Variables to remove have to be identified by hand

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
### Compare Parmeliaceae and Physciaceae


## Compare paths

physpaths = subset(phys_ests, op=='~') 
parmpaths = subset(parm_ests, op=='~') 

# Which paths differ by more than 0.1?
diff_cutoff = .1
physpaths$diff = abs(parmpaths$std.all - physpaths$std.all)>diff_cutoff
parmpaths$diff = abs(parmpaths$std.all - physpaths$std.all)>diff_cutoff

physpaths$sigcat = cut(abs(physpaths$std.all), c(0,.2,.4,.6,.8,1))
parmpaths$sigcat = cut(abs(parmpaths$std.all), c(0,.2,.4,.6,.8,1))

# Find all variables in the model
allvars = unique(c(parmpaths$rhs,parmpaths$lhs))
allvars[allvars=='regParm'] <-'reg' # Change based on response variables of interest
allvars[allvars=='parm_abun_log'] <- 'abun_log' # Change based on response variables of interest

# Subset paths to only those that are significant
parmpaths = parmpaths[which(apply(parmpaths[,c('std.ci.lower','std.ci.upper')],1,prod)>0),]
physpaths = physpaths[which(apply(physpaths[,c('std.ci.lower','std.ci.upper')],1,prod)>0),]

# Replace varnames in sigpaths to generic
physpaths$lhs[physpaths$lhs=='phys_abun_log']<-'abun_log'
physpaths$lhs[physpaths$lhs=='regPhys']<-'reg'
physpaths$rhs[physpaths$rhs=='phys_abun_log']<-'abun_log'
physpaths$rhs[physpaths$rhs=='regPhys']<-'reg'
parmpaths$lhs[parmpaths$lhs=='parm_abun_log']<-'abun_log'
parmpaths$lhs[parmpaths$lhs=='regParm']<-'reg'
parmpaths$rhs[parmpaths$rhs=='parm_abun_log']<-'abun_log'
parmpaths$rhs[parmpaths$rhs=='regParm']<-'reg'

# make sure var_locs is read-in from above.

# Make plots
mycex=1.5

svg('./Figures/phys parm model sig paths 0.1 diff.svg', height=10, width=8.5) # Change based on response variable of interest
par(mar=c(0,0,1.5,0))
par(mfrow=c(2,1))
par(lend="butt")
plot(var_locs, xlim=c(-6,6), ylim=c(-2.5,6), type='n', axes=F, xlab='', ylab='')
for(i in order(abs(parmpaths$std.all))){
	if(parmpaths$diff[i]){
		from = var_locs[parmpaths[i,'rhs'],]
		to = var_locs[parmpaths[i,'lhs'],]
		thickness = as.numeric(parmpaths[i,'sigcat'])*mycex
		col=c('red','black')[(parmpaths[i,'std.all']>0)+1]

		arrows(from$X, from$Y, to$X, to$Y, length=0.15, lwd=thickness, col=col)
	}
}
text(var_locs[fm_vars,], labels=varnames[fm_vars,'midName'], pos=2, offset=.25)
text(var_locs[fh_vars,], labels=varnames[fh_vars,'midName'], pos=4, offset=.25)
text(var_locs[c_vars,], labels=varnames[c_vars,'midName'],adj=-.1, srt=45)
text(var_locs[c('lichen_rich','abun_log'),], labels=c('Local\nrichness','Abundance'),pos=1, offset=1)
text(var_locs['radiation',], labels='Solar radiation', pos=1, offset=.25)
text(var_locs['totalNS',], labels=varnames['totalNS','midName'], pos=2, offset=.25)
text(var_locs['reg',], labels='Regional richness', pos=4, offset=.25)

mtext('Parmeliaceae', 3, 0, cex=1.2, font=2, adj=0)

legend(x=-5.5, y=6, xjust=0, yjust=1, c('0.0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1.0'), lwd=mycex*(1:5), 
	bty='n')
legend(x=-6, y=6, xjust=0, yjust=1, rep('', 5), lwd=mycex*(1:5), col='red', bty='n')
text(-5.5, 6, '+', cex=2, adj=-1)
text(-6, 6, '-', font=2, adj=-2.3)

plot(var_locs, xlim=c(-6,6), ylim=c(-2.5,6), type='n', axes=F, xlab='', ylab='')
for(i in order(abs(physpaths$std.all))){
	if(physpaths$diff[i]){
		from = var_locs[physpaths[i,'rhs'],]
		to = var_locs[physpaths[i,'lhs'],]
		thickness = as.numeric(physpaths[i,'sigcat'])*mycex
		col=c('red','black')[(physpaths[i,'std.all']>0)+1]

		arrows(from$X, from$Y, to$X, to$Y, length=0.15, lwd=thickness, col=col)
	}
}
text(var_locs[fm_vars,], labels=varnames[fm_vars,'midName'], pos=2, offset=.25)
text(var_locs[fh_vars,], labels=varnames[fh_vars,'midName'], pos=4, offset=.25)
text(var_locs[c_vars,], labels=varnames[c_vars,'midName'],adj=-.1, srt=45)
text(var_locs[c('lichen_rich','abun_log'),], labels=c('Local\nrichness','Abundance'),pos=1, offset=1)
text(var_locs['radiation',], labels='Solar radiation', pos=1, offset=.25)
text(var_locs['totalNS',], labels=varnames['totalNS','midName'], pos=2, offset=.25)
text(var_locs['reg',], labels='Regional richness', pos=4, offset=.25)

mtext('Physciaceae', 3, 0, cex=1.2, font=2, adj=0)
dev.off()


## Plot difference in total effects ordered by magnitude of difference

# Order variables from lowest to highest total effects
ordered_vars = rownames(parm)[order(abs(parm$std.all - phys$std.all))]
ordered_vars = ordered_vars[!(ordered_vars %in% c('FH','FM'))] # Drop effects of FM and FH categories (they may be non-sensical)

# Put tables in same order
use_parm = parm[ordered_vars,]
use_phys = phys[ordered_vars,]

# Define range limits that will include 95% confidence intervals
myrange = range(c(use_parm[,c('std.ci.lower','std.ci.upper')],
	use_phys[,c('std.ci.lower','std.ci.upper')]), na.rm=T)+c(-.04, .04)
myrange[1] = -.8

mycols = c('#2415B0','#00BF32')
mycolsbw = c('grey80','white')
names(mycols) = c('R','L')
names(mycolsbw) = c('R','L')
mycols_trans = paste(mycols, '50', sep='')
names(mycols_trans) = c('R','L')
mypch = c(22,23)
mypcols = c('white','grey30')
myadj=.15
mytypes = expression('C'['R'],'C'['L'],'F'['H'],'F'['O'],'C'['R'],'R','') # symbols used in plot to denote variable types
names(mytypes)=c('C','L','FH','FM','P','R','A')

# Total and direct standardized effects on same graph
svg('./Figures/Standardized total effects on Parm Phys richness bw.svg', height=12, width=19)
dotplot(as.numeric(factor(rownames(use_parm), levels = ordered_vars))~std.all, data=use_phys, 
	xlab=list('Standardized Total Effect',cex=3), ylab='',
	main='',cex.lab=3,aspect=9/10, xlim=myrange,
	panel=function(x,y){
	
		# Add horizontal boxes
		panel.rect(-2,1:length(ordered_vars)-.5, 2, 1:length(ordered_vars)+.5,
			col=mycolsbw[predtypes[ordered_vars,'scale']], border='grey50')
		
		# Add vertical line at 0
		panel.abline(v=0, col='grey30', lty=2, lwd=2)		
		
		# Add null distribution segments for Parmeliaceae
		panel.segments(use_parm$std.ci.lower, y+myadj,
			use_parm$std.ci.upper, y+myadj, 
			col='black', lwd=4.5, lend=1)
		# Add points for direct estimated effects
		panel.points(use_parm$std.all, y+myadj, col='black', fill=mypcols[2], pch=mypch[2], cex=3, lwd=3) 
		
		# Add null distribution segments for Physciaceae
		panel.segments(use_phys$std.ci.lower, y,
			use_phys$std.ci.upper, y, 
			col='black', lwd=4.5, lend=1)
		# Add points for total estimated effects
		panel.points(x, y, col='black', fill=mypcols[1], pch=mypch[1], cex=3, lwd=3) 
	
		# Add text labeling the variable type
		panel.text(-.75, y, labels=mytypes[predtypes[ordered_vars,'type']], cex=2)
		
	},
	scales=list(y=list(labels=varnames[ordered_vars,'midName'], 
		cex=3, col='black'),
		x=list(cex=3, tick.number=8)),
	key=list(x=1, y=0, corner=c(1,0), lines=list(type='o', pch=mypch, fill=mypcols, lwd=3, pt.lwd=3),
		text=list(c('Physciaceae','Parmaliaceae')),
		background='#FFFFFFCC', cex=3, divide=1, padding.text=5, border='black')
)
dev.off()

## Plot direct effects on richness vs abundance for Parm and Phys side-by-side

use_vars = rownames(predtypes)[predtypes$scale=='L']
use_vars = use_vars[use_vars!='abun_log']
use_vars = use_vars[-grep('soil', use_vars)]

mypch = c(22,23)
mycol=c('white','grey30')
myfact = ifelse(parm_d[use_vars,'type']=='FH', 1, 2)

svg('./Figures/compare richness abundance effects forest vars parm vs phys.svg', height=5, width=5 )
par(mar=c(4,4,1,1))
par(mfrow=c(1,2))

sigvars_a = names(which(apply(parm_da[use_vars,c('std.ci.lower','std.ci.upper')], 1, prod)>0))
sigvars_r = names(which(apply(parm_d[use_vars,c('std.ci.lower','std.ci.upper')],1,prod)>0))
sigvars = unique(c(sigvars_a,sigvars_r))
sigvars_L = sigvars[parm_d[sigvars,'std.all']<0]
sigvars_R = sigvars[parm_d[sigvars,'std.all']>0]

plot(parm_d[use_vars,'std.all'], parm_da[use_vars,'std.all'], type='n', 
	xlim=c(-.3, .3), ylim=c(-.4,.4), las=1, ylab='Direct Effect on Abundance',
	xlab='Direct Effect on Richness', cex.axis=1, cex.lab=1)
usr=par('usr')
polygon(c(usr[1],0,usr[1],usr[2],0,usr[2],usr[1]),
	c(usr[3],0,usr[4],usr[4],0,usr[3],usr[3]), col='grey80')
abline(h=0,v=0)
arrows(parm_d[use_vars,'std.ci.lower'], parm_da[use_vars,'std.all'],
	parm_d[use_vars,'std.ci.upper'], parm_da[use_vars,'std.all'],
	code=3, angle=90, lwd=2, length=0.05)
arrows(parm_d[use_vars,'std.all'], parm_da[use_vars,'std.ci.lower'],
	parm_d[use_vars,'std.all'], parm_da[use_vars,'std.ci.upper'],
	code=3, angle=90, lwd=2, length=0.05)
points(parm_d[use_vars,'std.all'], parm_da[use_vars,'std.all'], 
	pch=mypch[myfact], bg=mycol[myfact], lwd=2, col='black', cex=2)
legend('topright',c('Heterogeneity','Optimality'), pch=mypch, pt.bg=mycol, 
	pt.lwd=2, bg='white', box.lwd=1, pt.cex=2)

text(parm_d[sigvars_R, 'std.ci.upper']+.02,parm_da[sigvars_R,'std.all'], labels=varnames[sigvars_R, 'midName'], adj=0)
text(parm_d[sigvars_L, 'std.ci.lower']-.02,parm_da[sigvars_L,'std.all'], labels=varnames[sigvars_L, 'midName'], adj=1)

sigvars_a = names(which(apply(phys_da[use_vars,c('std.ci.lower','std.ci.upper')], 1, prod)>0))
sigvars_r = names(which(apply(phys_d[use_vars,c('std.ci.lower','std.ci.upper')],1,prod)>0))
sigvars = unique(c(sigvars_a,sigvars_r))
sigvars_L = sigvars[parm_d[sigvars,'std.all']<0]
sigvars_R = sigvars[parm_d[sigvars,'std.all']>0]

plot(phys_d[use_vars,'std.all'], phys_da[use_vars,'std.all'], type='n', 
	xlim=c(-.3, .3), ylim=c(-.3,.4), las=1, ylab='Direct Effect on Abundance',
	xlab='Direct Effect on Richness', cex.axis=1, cex.lab=1)
usr=par('usr')

xvals = seq(usr[1],usr[2],length.out=11)
fun1 = function(x) x
fun2 = function(x) -x
polygon(xvals, pmax(fun1(xvals), fun2(xvals)), col='grey80')
polygon(xvals, pmin(fun1(xvals), fun2(xvals)), col='grey80')

abline(h=0,v=0)
arrows(phys_d[use_vars,'std.ci.lower'], phys_da[use_vars,'std.all'],
	phys_d[use_vars,'std.ci.upper'], phys_da[use_vars,'std.all'],
	code=3, angle=90, lwd=2, length=0.05)
arrows(phys_d[use_vars,'std.all'], phys_da[use_vars,'std.ci.lower'],
	phys_d[use_vars,'std.all'], phys_da[use_vars,'std.ci.upper'],
	code=3, angle=90, lwd=2, length=0.05)
points(phys_d[use_vars,'std.all'], phys_da[use_vars,'std.all'], 
	pch=mypch[myfact], bg=mycol[myfact], lwd=2, col='black', cex=2)
legend('topright',c('Heterogeneity','Optimality'), pch=mypch, pt.bg=mycol, 
	pt.lwd=2, bg='white', box.lwd=1, pt.cex=2)

text(phys_d[sigvars_R, 'std.ci.upper']+.02,phys_da[sigvars_R,'std.all'], labels=varnames[sigvars_R, 'midName'], adj=0)
text(phys_d[sigvars_L, 'std.ci.lower']-.02,phys_da[sigvars_L,'std.all'], labels=varnames[sigvars_L, 'midName'], adj=1)



dev.off()


##################################################
## Functional Richness

## A table comparing direct effects of local variables on species richness vs functional richness
allsp_d$sig = apply(allsp_d[,c('std.ci.lower','std.ci.upper')], 1, prod)>0
fric_d$sig = apply(fric_d[,c('std.ci.lower','std.ci.upper')], 1, prod)>0
use_cols = c('std.all','std.se','sig')
compare_tab = cbind(predictor=varnames[rownames(allsp_d),'midName'],
	allsp_d[,use_cols], fric_d[rownames(allsp_d),use_cols])
compare_tab[order(compare_tab[,2], decreasing=T),]

use_allsp = subset(allsp_d, type %in% c('FH','FM','L'))
use_fric = fric_d[rownames(use_allsp),]
ordered_vars = rownames(use_allsp)[order(use_allsp$type,use_allsp$std.all)]
use_allsp = use_allsp[ordered_vars,]
use_fric = use_fric[ordered_vars,]

# Comparing direct effects on richness and fric
svg('./Figures/Standardized direct effects on Fric vs richness.svg', height=13, width=19)
dotplot(as.numeric(factor(rownames(use_allsp), levels = ordered_vars))~std.all, data=use_allsp, 
	xlab=list('Standardized Direct Effect',cex=3), ylab='',
	main='',cex.lab=3,aspect=9/10, xlim=c(-.2,.2),
	panel=function(x,y){
	
		# Add horizontal boxes
		panel.rect(-2,1:length(ordered_vars)-.5, 2, 1:length(ordered_vars)+.5,
			col='white', border='grey50')

		# Add vertical line at 0
		panel.abline(v=0, col='grey30', lty=2, lwd=2)		
		
		# Add null distribution segments for fric
		panel.segments(use_fric$std.ci.lower, y+myadj,
			use_fric$std.ci.upper, y+myadj, 
			col='black', lwd=4.5, lend=1)
		# Add points for direct estimated effects
		panel.points(use_fric$std.all, y+myadj, col='black', fill=mypcols[2], pch=mypch[2], cex=3, lwd=3) 
		
		# Add null distribution segments for richness
		panel.segments(use_allsp$std.ci.lower, y,
			use_allsp$std.ci.upper, y, 
			col='black', lwd=4.5, lend=1)
		# Add points for total estimated effects
		panel.points(x, y, col='black', fill=mypcols[1], pch=mypch[1], cex=3, lwd=3) 
	
		# Add text labeling the variable type
		panel.text(-.19, y, labels=mytypes[use_allsp$type], cex=2)
		
	},
	scales=list(y=list(labels=varnames[ordered_vars,'midName'], 
		cex=3, col='black'),
		x=list(cex=3, tick.number=8)),
	key=list(x=.5, y=1, corner=c(.5,0), lines=list(type='o', pch=mypch, fill=mypcols, lwd=3, pt.lwd=3),
		text=list(c('Species Richness','Functional Richness')),
		background='white', cex=3, divide=1, padding.text=5, border=F, columns=2)
)
dev.off()









