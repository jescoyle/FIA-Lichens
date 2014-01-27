# This script is used to fit SEM models of lichen FIA data

source('./UNC/Projects/FIA Lichen/GitHub/FIA-Lichens/load_data.R')
source('./GitHub/FIA-Lichens/fia_lichen_analysis_functions.R')

################################################################################
### Path Analysis ###

library(lavaan)
library(semPlot)
library(semTools)


# Read in table of predictor variable types
predtypes = read.csv('./SEM models/var_types.csv', row.names=1)
predtypes = subset(predtypes, !(rownames(predtypes) %in% c('FM','FH')))

# Full path model with measurement error model on lichen.rich_log
# 1/3/2014: Decided to use standardized data because the model will not fit to unscaled data and I did not want to apply an arbitrary scaling.
#           Standardized data is in working_data dataframe. Note that the response variables are also standardized.     
# I modified covariances and played with using log richness and adding abundace with this model before fitting the final model below.
path_measerr = '

	# Latent variables
	lichen_rich =~ sqrt(0.75)*lichen.rich_log

	# Climate effects on forest: both regional climate and local env (radiation) impact forest
	wood_SG.ba + wood_SG.rao.ba + bark_moist_pct.rao.ba + LogSeed.rao.ba + PIE.ba.tree + propDead + lightDist.mean + diamDiversity + bark_moist_pct.ba + LogSeed.ba + totalCirc + bigTrees + light.mean + PC1 ~
		wetness + rain_lowRH + iso + pseas + mat + radiation

	# Climate effects on regional richness: local env does not impact regional richness
	regS ~ wetness + rain_lowRH + iso + pseas + mat

	# Climate effects on pollution
	totalNS ~ wetness + rain_lowRH

	# Covariances
	PIE.ba.tree ~~ bark_moist_pct.rao.ba + LogSeed.rao.ba + wood_SG.rao.ba
	PIE.ba.tree ~~ diamDiversity + bigTrees + totalCirc + PC1
	bark_moist_pct.rao.ba ~~ LogSeed.rao.ba + wood_SG.rao.ba
	LogSeed.rao.ba ~~ wood_SG.rao.ba
	LogSeed.ba ~~ wood_SG.ba + bark_moist_pct.ba
	wood_SG.ba ~~ bark_moist_pct.ba
	totalCirc ~~ diamDiversity + bigTrees
	light.mean ~~ diamDiversity + totalCirc
	totalNS ~~ iso + pseas + mat
	regS ~~ totalNS
	rain_lowRH + wetness ~~ iso + pseas + mat + radiation
	mat ~~ iso + pseas + radiation
	iso ~~ pseas + radiation
	pseas ~~ radiation

	lichen_rich ~ bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba +
		PIE.ba.tree + propDead + lightDist.mean + diamDiversity +
		bark_moist_pct.ba + wood_SG.ba + LogSeed.ba + totalCirc + bigTrees + light.mean + PC1 +
		wetness + rain_lowRH + iso + pseas + mat + radiation + regS + totalNS +
		tot_abun_log
	
	tot_abun_log ~ bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba +
		PIE.ba.tree + propDead + lightDist.mean + diamDiversity +
		bark_moist_pct.ba + wood_SG.ba + LogSeed.ba + totalCirc + bigTrees + light.mean + PC1 +
		wetness + rain_lowRH + iso + pseas + mat + radiation + totalNS
'

# Fit model
fit_measerr = sem(path_measerr, data=working_data_fit, fixed.x=F)

summary(fit_measerr, standardized=T, rsq=TRUE, fit.measures=T)
fitMeasures(fit_measerr)
measerr_ests = parameterEstimates(fit_measerr, standardized=T)

modindices(fit_measerr)
mi = miPowerFit(fit_measerr)
problem_parms = subset(mi, decision %in% c('M', 'EPC:M'))
problem_parms[order(problem_parms$mi),]
res = resid(fit_measerr)
which(abs(res$cov)>.2, arr.ind=T)

## Define model that will calculate indirect and total effects
path_measerr = '

	# Latent variables
	lichen_rich =~ sqrt(0.75)*lichen.rich_log

	# Climate and local environment effects on forest
	bark_moist_pct.rao.ba ~ FH1a*wetness + FH1b*rain_lowRH + FH1c*iso + FH1d*pseas + FH1e*mat + FH1f*radiation
	wood_SG.rao.ba ~ FH2a*wetness + FH2b*rain_lowRH + FH2c*iso + FH2d*pseas + FH2e*mat +  FH2f*radiation
	LogSeed.rao.ba ~ FH3a*wetness + FH3b*rain_lowRH + FH3c*iso + FH3d*pseas + FH3e*mat + FH3f*radiation
	PIE.ba.tree ~ FH4a*wetness + FH4b*rain_lowRH + FH4c*iso + FH4d*pseas + FH4e*mat + FH4f*radiation
	propDead ~ FH5a*wetness + FH5b*rain_lowRH + FH5c*iso + FH5d*pseas + FH5e*mat + FH5f*radiation
	lightDist.mean ~ FH6a*wetness + FH6b*rain_lowRH + FH6c*iso + FH6d*pseas + FH6e*mat + FH6f*radiation
	diamDiversity ~ FH7a*wetness + FH7b*rain_lowRH + FH7c*iso + FH7d*pseas + FH7e*mat + FH7f*radiation
	bark_moist_pct.ba ~ FM1a*wetness + FM1b*rain_lowRH + FM1c*iso + FM1d*pseas + FM1e*mat + FM1f*radiation
	wood_SG.ba ~ FM2a*wetness + FM2b*rain_lowRH + FM2c*iso + FM2d*pseas + FM2e*mat + FM2f*radiation
	LogSeed.ba ~ FM3a*wetness + FM3b*rain_lowRH + FM3c*iso + FM3d*pseas + FM3e*mat + FM3f*radiation
	totalCirc ~ FM4a*wetness + FM4b*rain_lowRH + FM4c*iso + FM4d*pseas + FM4e*mat + FM4f*radiation
	bigTrees ~ FM5a*wetness + FM5b*rain_lowRH + FM5c*iso + FM5d*pseas + FM5e*mat + FM5f*radiation
	light.mean  ~ FM6a*wetness + FM6b*rain_lowRH + FM6c*iso + FM6d*pseas + FM6e*mat + FM6f*radiation
	PC1 ~ FM7a*wetness + FM7b*rain_lowRH + FM7c*iso + FM7d*pseas + FM7e*mat + FM7f*radiation

	# Climate effects on regional richness
	regS ~ R1*wetness + R2*rain_lowRH + R3*iso + R4*pseas + R5*mat

	# Climate effects on pollution
	totalNS ~ P1*wetness + P2*rain_lowRH

	# Covariances
	# NO LONGER USING fixed.x=F B/C IT IS BETTER TO KEEP EXOGENOUS VARIABLES AS FIXED (NON-RANDOM)
	# COVARIANCES AUTOMATICALLY CALCULATED
	#PIE.ba.tree ~~ bark_moist_pct.rao.ba + LogSeed.rao.ba + wood_SG.rao.ba
	#PIE.ba.tree ~~ diamDiversity + bigTrees + totalCirc + PC1
	#bark_moist_pct.rao.ba ~~ LogSeed.rao.ba + wood_SG.rao.ba
	#LogSeed.rao.ba ~~ wood_SG.rao.ba
	#LogSeed.ba ~~ wood_SG.ba + bark_moist_pct.ba
	#wood_SG.ba ~~ bark_moist_pct.ba
	#totalCirc ~~ diamDiversity + bigTrees
	#light.mean ~~ diamDiversity + totalCirc
	#totalNS ~~ iso + pseas + mat
	#regS ~~ totalNS
	#rain_lowRH + wetness ~~ iso + pseas + mat + radiation
	#mat ~~ iso + pseas + radiation
	#iso ~~ pseas + radiation
	#pseas ~~ radiation
	
	# Lichen richness regression
	lichen_rich ~ FH1*bark_moist_pct.rao.ba + FH2*wood_SG.rao.ba + FH3*LogSeed.rao.ba +
		FH4*PIE.ba.tree + FH5*propDead + FH6*lightDist.mean + FH7*diamDiversity +
		FM1*bark_moist_pct.ba + FM2*wood_SG.ba + FM3*LogSeed.ba + FM4*totalCirc + FM5*bigTrees + FM6*light.mean + FM7*PC1 +
		C1*wetness + C2*rain_lowRH + C3*iso + C4*pseas + C5*mat + C6*radiation + R*regS + P*totalNS +
		A*tot_abun_log

	# Lichen abundance regression
	tot_abun_log ~ AFH1*bark_moist_pct.rao.ba + AFH2*wood_SG.rao.ba + AFH3*LogSeed.rao.ba +
		AFH4*PIE.ba.tree + AFH5*propDead + AFH6*lightDist.mean + AFH7*diamDiversity +
		AFM1*bark_moist_pct.ba + AFM2*wood_SG.ba + AFM3*LogSeed.ba + AFM4*totalCirc + AFM5*bigTrees + AFM6*light.mean + AFM7*PC1 +
		AC1*wetness + AC2*rain_lowRH + AC3*iso + AC4*pseas + AC5*mat + AC6*radiation + AP*totalNS
	
	# Indirect effects of forest structure variables on richness via abundance
	IE_bark_moist_pct.rao.ba := AFH1*A
	IE_wood_SG.rao.ba := AFH2*A
	IE_LogSeed.rao.ba := AFH3*A
	IE_PIE.ba.tree := AFH4*A
	IE_propDead := AFH5*A
	IE_lightDist.mean := AFH6*A
	IE_diamDiversity := AFH7*A
	IE_bark_moist_pct.ba := AFM1*A
	IE_wood_SG.ba := AFM2*A
	IE_LogSeed.ba := AFM3*A
	IE_totalCirc := AFM4*A
	IE_bigTrees := AFM5*A
	IE_light.mean := AFM6*A
	IE_PC1 := AFM7*A

	# Indirect effect of pollution via abundance
	IE_totalNS := AP*A
	
	# Indirect effects of climate and local environment on lichen richness
	IE_wetness_FH := FH1a*(FH1 + IE_bark_moist_pct.rao.ba) + FH2a*(FH2 + IE_wood_SG.rao.ba) +
		FH3a*(FH3 + IE_LogSeed.rao.ba) + FH4a*(FH4 + IE_PIE.ba.tree) + FH5a*(FH5 + IE_propDead) + 
		FH6a*(FH6 + IE_lightDist.mean) + FH7a*(FH7 + IE_diamDiversity)
	IE_wetness_FM := FM1a*(FM1 + IE_bark_moist_pct.ba) + FM2a*(FM2 + IE_wood_SG.ba) + 
		FM3a*(FM3 + IE_LogSeed.ba) + FM4a*(FM4 + IE_totalCirc) + FM5a*(FM5 + IE_bigTrees) + 
		FM6a*(FM6 + IE_light.mean) + FM7a*(FM7 + IE_PC1)
	IE_wetness_R := R1*R
	IE_wetness_P := P1*(P + IE_totalNS)
	IE_wetness_A := AC1*A
	IE_wetness := IE_wetness_FH + IE_wetness_FM + IE_wetness_R + IE_wetness_P + IE_wetness_A

	IE_rain_lowRH_FH := FH1b*(FH1 + IE_bark_moist_pct.rao.ba) + FH2b*(FH2 + IE_wood_SG.rao.ba) +
		FH3b*(FH3 + IE_LogSeed.rao.ba) + FH4b*(FH4 + IE_PIE.ba.tree) + FH5b*(FH5 + IE_propDead) + 
		FH6b*(FH6 + IE_lightDist.mean) + FH7b*(FH7 + IE_diamDiversity)
	IE_rain_lowRH_FM := FM1b*(FM1 + IE_bark_moist_pct.ba) + FM2b*(FM2 + IE_wood_SG.ba) + 
		FM3b*(FM3 + IE_LogSeed.ba) + FM4b*(FM4 + IE_totalCirc) + FM5b*(FM5 + IE_bigTrees) + 
		FM6b*(FM6 + IE_light.mean) + FM7b*(FM7 + IE_PC1)
	IE_rain_lowRH_R := R2*R
	IE_rain_lowRH_P := P2*(P + IE_totalNS)
	IE_rain_lowRH_A := AC2*A
	IE_rain_lowRH := IE_rain_lowRH_FH + IE_rain_lowRH_FM + IE_rain_lowRH_R + IE_rain_lowRH_P + IE_rain_lowRH_A

	IE_iso_FH := FH1c*(FH1 + IE_bark_moist_pct.rao.ba) + FH2c*(FH2 + IE_wood_SG.rao.ba) +
		FH3c*(FH3 + IE_LogSeed.rao.ba) + FH4c*(FH4 + IE_PIE.ba.tree) + FH5c*(FH5 + IE_propDead) + 
		FH6c*(FH6 + IE_lightDist.mean) + FH7c*(FH7 + IE_diamDiversity)
	IE_iso_FM := FM1c*(FM1 + IE_bark_moist_pct.ba) + FM2c*(FM2 + IE_wood_SG.ba) + 
		FM3c*(FM3 + IE_LogSeed.ba) + FM4c*(FM4 + IE_totalCirc) + FM5c*(FM5 + IE_bigTrees) +
		FM6c*(FM6 + IE_light.mean) +  FM7c*(FM7 + IE_PC1)
	IE_iso_R := R3*R
	IE_iso_A := AC3*A
	IE_iso := IE_iso_FH + IE_iso_FM + IE_iso_R + IE_iso_A

	IE_pseas_FH := FH1d*(FH1 + IE_bark_moist_pct.rao.ba) + FH2d*(FH2 + IE_wood_SG.rao.ba) +
		FH3d*(FH3 + IE_LogSeed.rao.ba) + FH4d*(FH4 + IE_PIE.ba.tree) + FH5d*(FH5 + IE_propDead) + 
		FH6d*(FH6 + IE_lightDist.mean) + FH7d*(FH7 + IE_diamDiversity)
	IE_pseas_FM := FM1d*(FM1 + IE_bark_moist_pct.ba) + FM2d*(FM2 + IE_wood_SG.ba) + 
		FM3d*(FM3 + IE_LogSeed.ba) + FM4d*(FM4 + IE_totalCirc) + FM5d*(FM5 + IE_bigTrees) + 
		FM6d*(FM6 + IE_light.mean) + FM7d*(FM7 + IE_PC1)
	IE_pseas_R := R4*R
	IE_pseas_A := AC4*A
	IE_pseas := IE_pseas_FH + IE_pseas_FM + IE_pseas_R + IE_pseas_A
	
	IE_mat_FH := FH1e*(FH1 + IE_bark_moist_pct.rao.ba) + FH2e*(FH2 + IE_wood_SG.rao.ba) +
		FH3e*(FH3 + IE_LogSeed.rao.ba) + FH4e*(FH4 + IE_PIE.ba.tree) + FH5e*(FH5 + IE_propDead) + 
		FH6e*(FH6 + IE_lightDist.mean) + FH7e*(FH7 + IE_diamDiversity)
	IE_mat_FM := FM1e*(FM1 + IE_bark_moist_pct.ba) + FM2e*(FM2 + IE_wood_SG.ba) + 
		FM3e*(FM3 + IE_LogSeed.ba) + FM4e*(FM4 + IE_totalCirc) + FM5e*(FM5 + IE_bigTrees) + 
		FM6e*(FM6 + IE_light.mean) + FM7e*(FM7 + IE_PC1)
	IE_mat_R := R5*R
	IE_mat_A := AC5*A
	IE_mat := IE_mat_FH + IE_mat_FM + IE_mat_R + IE_mat_A

	IE_radiation_FH := FH1f*(FH1 + IE_bark_moist_pct.rao.ba) + FH2f*(FH2 + IE_wood_SG.rao.ba) +
		FH3f*(FH3 + IE_LogSeed.rao.ba) + FH4f*(FH4 + IE_PIE.ba.tree) + FH5f*(FH5 + IE_propDead) + 
		FH6f*(FH6 + IE_lightDist.mean) + FH7f*(FH7 + IE_diamDiversity)
	IE_radiation_FM := FM1f*(FM1 + IE_bark_moist_pct.ba) + FM2f*(FM2 + IE_wood_SG.ba) + 
		FM3f*(FM3 + IE_LogSeed.ba) + FM4f*(FM4 + IE_totalCirc) + FM5f*(FM5 + IE_bigTrees) + 
		FM6f*(FM6 + IE_light.mean) + FM7f*(FM7 + IE_PC1)
	IE_radiation_A := AC6*A
	IE_radiation := IE_radiation_FH + IE_radiation_FM + IE_radiation_A

	# Total effects on lichen richness
	TE_wetness := IE_wetness + C1	
	TE_rain_lowRH := IE_rain_lowRH + C2
	TE_iso := IE_iso + C3
	TE_pseas := IE_pseas + C4
	TE_mat := IE_mat + C5
	TE_radiation := IE_radiation + C6

	TE_bark_moist_pct.rao.ba := IE_bark_moist_pct.rao.ba + FH1
	TE_wood_SG.rao.ba := IE_wood_SG.rao.ba + FH2
	TE_LogSeed.rao.ba := IE_LogSeed.rao.ba + FH3
	TE_PIE.ba.tree := IE_PIE.ba.tree + FH4
	TE_propDead := IE_propDead + FH5
	TE_lightDist.mean := IE_lightDist.mean + FH6
	TE_diamDiversity := IE_diamDiversity + FH7
	TE_bark_moist_pct.ba := IE_bark_moist_pct.ba + FM1
	TE_wood_SG.ba := IE_wood_SG.ba + FM2
	TE_LogSeed.ba := IE_LogSeed.ba + FM3
	TE_totalCirc := IE_totalCirc + FM4
	TE_bigTrees := IE_bigTrees + FM5
	TE_light.mean := IE_light.mean + FM6
	TE_PC1 := IE_PC1 +  FM7	

	TE_FH := TE_bark_moist_pct.rao.ba +	TE_wood_SG.rao.ba + TE_LogSeed.rao.ba +
		TE_PIE.ba.tree + TE_propDead + TE_lightDist.mean + TE_diamDiversity
	TE_FM := TE_bark_moist_pct.ba + TE_wood_SG.ba + TE_LogSeed.ba + TE_totalCirc +
		TE_bigTrees + TE_light.mean + TE_PC1

	TE_totalNS := IE_totalNS + P

'
endfit =  sem(path_measerr, data=working_data_test, fixed.x=T, estimator='ML', se='robust.sem')
summary(endfit, standardized=T, rsq=TRUE, fit.measures=T)

# Examine model residuals
res = resid(endfit)
which(abs(res$cov)>0.2, arr.ind=T)

# Save model output
save(endfit, file='./SEM models/path_measerr_finalmod_testdata.Rdata')

# Examine significance of model paths
endfit_ests = parameterEstimates(endfit, standardized=T, ci=T, level=0.95)
subset(endfit_ests, pvalue>0.05)
endfit_paths = subset(endfit_ests, op=='~')
endfit_paths[order(endfit_paths$std.all, decreasing=T),]

write.csv(endfit_ests, './SEM models/finalmod_testdata_estimates.csv', row.names=F)

# Get standardized coeficient estimates
stdsol = standardizedSolution(endfit)

# Compare standardized coefs
cbind(stdsol$se, endfit_ests$se)

# Set up data for bootstrapping standardized coefficients
# This will be done on the cluster
endfit_std = bootstrapLavaan(endfit, R=3, FUN=function(x) c(parameterEstimates(x)$est,standardizedSolution(x)$est.std))

# Make a data frame that has all response variables in it
new_response = master[rownames(working_data),c('fric','fdiv','raoQ')]
new_response$Parm_log = log(master[rownames(working_data),'Parmeliaceae']+1)
new_response$Phys_log = log(master[rownames(working_data),'Physciaceae']+1)
new_response = scale(new_response, center=T, scale=T)
alldata=data.frame(working_data,new_response)

write.csv(alldata[testplots$yrplot.id,], './SEM models/Bootstrap/standardized_test_dataset.csv', row.names=T)

alldata = read.csv('./SEM models/Bootstrap/standardized_test_dataset.csv', row.names=1)
endfit =  sem(path_measerr, data=alldata, fixed.x=F, estimator='ML', se='robust.sem')

### See scripts for bootstrap calculations done on the cluster
# kure_sem_boot_[allsp, parm, phys, fric, raoQ].R


## Read in tables of parameter estimates and effects
# All parameter estimares
allsp_ests = read.csv('./SEM models/AllSp_testdata_parameterEstimates.csv')

# Total effects
allsp = read.csv('./SEM models/AllSp_testdata_totaleffects.csv', row.names=1)
parm = read.csv('./SEM models/Parm_testdata_totaleffects.csv', row.names=1)
phys = read.csv('./SEM models/Phys_testdata_totaleffects.csv', row.names=1)
fric = read.csv('./SEM models/Fric_testdata_totaleffects.csv', row.names=1)
raoq = read.csv('./SEM models/RaoQ_testdata_totaleffects.csv', row.names=1)

# Direct effecst
allsp_d = read.csv('./SEM models/AllSp_testdata_directeffects_richness.csv', row.names=1)
parm_d = read.csv('./SEM models/Parm_testdata_directeffects_richness.csv', row.names=1)
phys_d = read.csv('./SEM models/Phys_testdata_directeffects_richness.csv', row.names=1)
fric_d = read.csv('./SEM models/Fric_testdata_directeffects_richness.csv', row.names=1)
raoq_d = read.csv('./SEM models/RaoQ_testdata_directeffects_richness.csv', row.names=1)

# Direct effect on abundance
allsp_da = read.csv('./SEM models/AllSp_testdata_directeffects_abundance.csv', row.names=1)

# Indirect effects via abundance
allsp_i = read.csv('./SEM models/AllSp_testdata_indirecteffects_via_abundance.csv', row.names=1)
parm_i = read.csv('./SEM models/Parm_testdata_indirecteffects_via_abundance.csv', row.names=1)
phys_i = read.csv('./SEM models/Phys_testdata_indirecteffects_via_abundance.csv', row.names=1)
fric_i = read.csv('./SEM models/Fric_testdata_indirecteffects_via_abundance.csv', row.names=1)
raoq_i = read.csv('./SEM models/RaoQ_testdata_indirecteffects_via_abundance.csv', row.names=1)

# Indirect effects via regional richness
allsp_ir = read.csv('./SEM models/AllSp_testdata_indirecteffects_via_regS.csv', row.names=1)
parm_ir = read.csv('./SEM models/Parm_testdata_indirecteffects_via_regS.csv', row.names=1)
phys_ir = read.csv('./SEM models/Phys_testdata_indirecteffects_via_regS.csv', row.names=1)
fric_ir = read.csv('./SEM models/Fric_testdata_indirecteffects_via_regS.csv', row.names=1)
raoq_ir = read.csv('./SEM models/RaoQ_testdata_indirecteffects_via_regS.csv', row.names=1)

# Indirect effects via forest structure
allsp_if = read.csv('./SEM models/AllSp_testdata_indirecteffects_via_forest.csv', row.names=1)
parm_if = read.csv('./SEM models/Parm_testdata_indirecteffects_via_forest.csv', row.names=1)
phys_if = read.csv('./SEM models/Phys_testdata_indirecteffects_via_forest.csv', row.names=1)
fric_if = read.csv('./SEM models/Fric_testdata_indirecteffects_via_forest.csv', row.names=1)
raoq_if = read.csv('./SEM models/RaoQ_testdata_indirecteffects_via_forest.csv', row.names=1)

##################################################################
### Results ###

# How many predictors have significant total effects?
names(which(apply(allsp[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0))
tot_sig = allsp[which(apply(allsp[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0),]

tot_sig = tot_sig[order(abs(tot_sig$std.all)),]


# Examine order of effects among forest structure variables
allsp_f = subset(allsp, type %in% c('FH','FM'))
allsp_df = subset(allsp_d, type %in% c('FH','FM'))

allsp_f[order(abs(allsp_f$std.all)),]
allsp_df[order(abs(allsp_df$std.all)),]


#######################################################################################
### Figures ###


######################################
## Compare direct effects on richness vs abundance

## Compare direct effects vs. indirect effects via abundance vs. total effects

# Order variables from lowest to highest total effects
total = fric # This changes based on what response variable is being analyzed
ordered_vars = rownames(total[order(total$std.all),])
ordered_vars = ordered_vars[!(ordered_vars %in% c('FH','FM'))] # Drop effects of FM and FH categories (they may be non-sensical)

# Put tables in same order
use_direct = fric_d[ordered_vars,] # This changes based on what response variable is being analyzed
use_total = total[ordered_vars,]

library(lattice)

# Define range limits that will include 95% confidence intervals
myrange = range(c(use_total[,c('std.ci.lower','std.ci.upper')],
	use_direct[,c('std.ci.lower','std.ci.upper')]), na.rm=T)+c(-.04, .04)
myrange[1] = -.8

# Color scheme: http://colorschemedesigner.com/#3341SsYrGvyw0
mycols = c('#000000','#07395D','#066788','#688e60','#05633B','#902E07','#DD8615')
names(mycols) = c('A','C','L','FH','FM','P','R')
mycols_trans = paste(mycols, '50', sep='')
plot(1:7,1:7, type='n'); text(1:7,1:7,labels=names(mycols), cex=2, col=mycols) # Check colors
names(mycols_trans) = c('A','C','L','FH','FM','P','R')
mypch = c(22,23)
mypcols = c('white','grey30')
myadj=.15
mytypes = expression('C'['R'],'C'['L'],'F'['H'],'F'['O'],'C'['R'],'R','') # symbols used in plot to denote variable types
names(mytypes)=c('C','L','FH','FM','P','R','A')

# Total and direct standardized effects on same graph
svg('./Figures/New Coordinates/Standardized direct total effects on Fric.svg', height=12, width=19)
dotplot(as.numeric(factor(rownames(use_total), levels = ordered_vars))~std.all, data=use_total, 
	xlab=list('Standardized Effect',cex=3), ylab='',
	main='',cex.lab=3,aspect=9/10, xlim=myrange,
	panel=function(x,y){
	
		# Add horizontal boxes
		panel.rect(-2,1:length(ordered_vars)-.5, 2, 1:length(ordered_vars)+.5,
			col='white', border='grey50')
		
		# Add vertical line at 0
		panel.abline(v=0, col='grey30', lty=2, lwd=2)		
		
		# Add null distribution segments for direct effects
		panel.segments(use_direct$std.ci.lower, y+myadj,
			use_direct$std.ci.upper, y+myadj, 
			col='black', lwd=4.5, lend=1)
		# Add points for direct estimated effects
		panel.points(use_direct$std.all, y+myadj, col='black', fill=mypcols[2], pch=mypch[2], cex=3, lwd=3) 
		
		# Add null distribution segments for total effects
		panel.segments(use_total$std.ci.lower, y,
			use_total$std.ci.upper, y, 
			col='black', lwd=4.5, lend=1)
		# Add points for total estimated effects
		panel.points(x, y, col='black', fill=mypcols[1], pch=mypch[1], cex=3, lwd=3) 
	
		# Add text labeling the variable type
		panel.text(-.75, y, labels=mytypes[use_total$type], cex=2)
		
	},
	scales=list(y=list(labels=varnames[ordered_vars,'midName'], 
		cex=3, col='black'),
		x=list(cex=3, tick.number=8)),
	key=list(x=1, y=0, corner=c(1,0), lines=list(type='o', pch=mypch, fill=mypcols, lwd=3, pt.lwd=3),
		text=list(c('Total effect','Direct effect')),
		background='white', cex=3, divide=1, padding.text=5, border='black')
)
dev.off()

## Larger figure including more indirect effects via panels
myrange = range(c(use_total[,c('std.ci.lower','std.ci.upper')],
	use_direct[,c('std.ci.lower','std.ci.upper')],
	use_indirect[,c('std.ci.lower','std.ci.upper')],
	use_indirectF[,c('std.ci.lower','std.ci.upper')],
	use_indirectR[,c('std.ci.lower','std.ci.upper')]), na.rm=T)+c(-.04, .04)

# Climate: direct/total
use_total_sub = subset(use_total, type %in% c('C','L'))
use_direct_sub = subset(use_direct, type %in% c('C','L'))
use_indirect_sub = subset(use_indirect, type %in% c('C','L'))
order_clim = use_direct_sub[order(use_direct_sub$std.all),'predictor']

use_total_sub = use_total_sub[order_clim,]
use_direct_sub = use_direct_sub[order_clim,]
use_indirect_sub = use_indirect_sub[order_clim,]
use_indirectR_sub = use_indirectR[order_clim[order_clim %in% rownames(use_indirectR)],]
use_indirectF_sub = use_indirectF_sub[c(paste(order_clim,'FM', sep='_'),paste(order_clim,'FH', sep='_')),]

svg('./Figures/New Coordinates/Effects figure P1 climate direct total.svg', height=2.25, width=7, bg='transparent', pointsize=8)
dotplot(as.numeric(factor(use_total_sub$predictor, levels = order_clim))~std.all, data=use_total_sub, 
	xlab=list('Standardized Effect',cex=1), ylab='',
	main='',cex.lab=1,aspect=5/12, xlim=myrange,
	panel=function(x,y){
		
		# Add horizontal boxes
		#panel.rect(-2,1:length(order_clim)-.5, 2, 1:length(order_clim)+.5,
		#	col=mycols_trans[use_total_sub[order_clim,'type']], border='transparent' )

		# Add vertical line at 0
		panel.abline(v=0, col='black', lty=1, lwd=1)		
		
		# Add 95% CI segments for total effects
		panel.segments(use_total_sub$std.ci.lower, y,
			use_total_sub$std.ci.upper, y, 
			col=mypcols[2], lwd=2, lend=1)
		# Add points for total estimated effects
		panel.points(x, y, col=mypcols[2], pch=mypch[1], cex=1, lwd=2) 
		
		# Add 95% CI for direct effects
		panel.segments(use_direct_sub$std.ci.lower, y+myadj,
			use_direct_sub$std.ci.upper, y+myadj, 
			col=mypcols[1], lwd=2, lend=1)
		# Add points for direct estimated effects
		panel.points(use_direct_sub$std.all, y+myadj, col=mypcols[1], pch=mypch[3], cex=1, lwd=2) 
			
	},
	scales=list(y=list(labels=varnames[order_clim,'midName'], 
		cex=1, col=mycols[use_total_sub[order_clim,'type']]),
		x=list(cex=1, tick.number=8)),
	key=list(x=1, y=.48, corner=c(0,.5), lines=list(type='o', pch=mypch[1:2], col=mypcols[1:2], lwd=1, size=2),
		text=list(c('Direct effect','Total effect')),
		background='#ffffff00', cex=1, divide=1, padding.text=2, between=1)
)
dev.off()

###########################################
## Make table of indirect climate effects

climEff_tab = data.frame(varnames[use_direct_sub$predictor,'displayName'], 
	direct = use_direct_sub$std.all,
	directSig = apply(use_direct_sub[,c('std.ci.lower','std.ci.upper')], 1, prod)>0,
	indirectA = use_indirect_sub$std.all,
	indirectASig = apply(use_indirect_sub[,c('std.ci.lower','std.ci.upper')], 1, prod)>0,
	indirectFH = subset(use_indirectF_sub, Ftype=='FH')$std.all,
	indirectFHSig = apply(subset(use_indirectF_sub, Ftype=='FH')[,c('std.ci.lower','std.ci.upper')], 1, prod)>0,
	indirectFM = subset(use_indirectF_sub, Ftype=='FM')$std.all,
	indirectFMSig = apply(subset(use_indirectF_sub, Ftype=='FM')[,c('std.ci.lower','std.ci.upper')], 1, prod)>0,
	indirectR = use_indirectR[order_clim,]$std.all,
	indirectRSig = apply(use_indirectR[order_clim,c('std.ci.lower','std.ci.upper')], 1, prod)>0
)

write.csv(climEff_tab, './SEM models/Compare effects climate variables.csv', row.names=F)


################################################
## Plot direct vs. indirect effects via abundance for local predictors
use_vars = rownames(predtypes)[predtypes$scale=='L']
use_vars = use_vars[use_vars!='abun_log']
use_vars = use_vars[-grep('soil', use_vars)]

mypch = c(22,23)
mycol=c('white','grey30')
myfact = ifelse(allsp_d[use_vars,'type']=='FH', 1, 2)

svg('./Figures/compare richness abundance effects forest vars.svg', height=6, width=6 )
par(mar=c(4,4,1,1))
plot(allsp_d[use_vars,'std.all'], allsp_da[use_vars,'std.all'], type='n', 
	xlim=c(-.3, .3), ylim=c(-.25,.25), las=1, ylab='Direct Effect on Abundance',
	xlab='Direct Effect on Richness', cex.axis=1.2, cex.lab=1.2)
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
legend('bottomright',c('Heterogeneity','Optimality'), pch=mypch, pt.bg=mycol, 
	pt.lwd=2, bg='white', box.lwd=1)

text(-.07,.22,'Significant Effect\non Abundance', font=2, adj=1)
sigvars_a = names(which(apply(allsp_da[use_vars,c('std.ci.lower','std.ci.upper')], 1, prod)>0))
text(-.07, allsp_da[sigvars_a,'std.all'], labels=varnames[sigvars_a, 'midName'],
	adj=1)

text(.08,.22,'Significant Effect\non Richness', font=2, adj=0)
sigvars_r = names(which(apply(allsp_d[use_vars,c('std.ci.lower','std.ci.upper')],1,prod)>0))
text(0.08, allsp_da[sigvars_r,'std.all'], labels=varnames[sigvars_r,'midName'], adj=0)

dev.off()



##################################################
### Draw path diagram for significant paths

# Only look at estimates corresponding to paths
paths = subset(allsp_ests, op=='~')

# Calculate which paths are significant
sigpaths = paths[which(apply(paths[,c('std.ci.lower','std.ci.upper')],1,prod)>0),]

# Histogram of significant paths
hist(sigpaths$std.all)
sigpaths$sigcat = cut(abs(sigpaths$std.all), c(0,.2,.4,.6,.8,1))

# Find all variables involved in significant relationships
sigvars  = unique(c(sigpaths$lhs, sigpaths$rhs))

## Create a matrix of variable locations
var_locs = matrix(0, nrow=length(sigvars), ncol=2, byrow=T)
colnames(var_locs) = c('X','Y')
rownames(var_locs) = sigvars

unit = 1

# Start with richness and abundance in center
var_locs['tot_abun_log',] = c(-1,0)
var_locs['lichen_rich',] = c(1,0)

# Add pollution and regS to corners
var_locs['totalNS',] = c(-2.5,2.5)
var_locs['regS',] = c(2.5,2.5)

# Add climate and local environment vars to top and bottom
c_vars = rownames(subset(predtypes[sigvars,], type=='C'))
var_locs[c_vars,] = cbind(3/(length(c_vars)-1)*(0:(length(c_vars)-1))-1.5,rep(3,length(c_vars)))
var_locs['radiation',] = c(0,-2.5)

# Add forest vars to sides
fh_vars = rownames(subset(predtypes[sigvars,], type=='FH'))
fm_vars = rownames(subset(predtypes[sigvars,], type=='FM'))

var_locs[fh_vars,]=cbind(rep(3,length(fh_vars)), 3.5/(length(fh_vars)-1)*(0:(length(fh_vars)-1))-2)
var_locs[fm_vars,]=cbind(rep(-3,length(fh_vars)), 3.5/(length(fh_vars)-1)*(0:(length(fh_vars)-1))-2)

var_locs=data.frame(var_locs)

# Make plot
mycex=1.5

svg('./Figures/full model significant paths.svg', height=5, width=8.5)
par(mar=c(0,0,0,0))
par(lend="butt")
plot(var_locs, xlim=c(-6,6), ylim=c(-2.5,6), type='n', axes=F, xlab='', ylab='')
for(i in order(abs(sigpaths$std.all))){
	from = var_locs[sigpaths[i,'rhs'],]
	to = var_locs[sigpaths[i,'lhs'],]
	thickness = as.numeric(sigpaths[i,'sigcat'])*mycex
	col=c('red','black')[(sigpaths[i,'std.all']>0)+1]

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

## Only show variables with significant paths to richness

# Variables to remove
subset(sigpaths, rhs=='PC1') # Check variables that are difficult to distinguish on path diagram
remove_vars = c('lightDist.mean','propDead','LogSeed.rao.ba','radiation','bark_moist_pct.ba')
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
svg('./Figures/New Coordinates/Standardized direct effects on Fric vs richness.svg', height=13, width=19)
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












############ OLD CODE ####################

## Original estimates smear graph with direct, indirect and total effects
# Color scheme: http://colorschemedesigner.com/#3341SsYrGvyw0
mycols = c('#000000','#07395D','#066788','#688e60','#05633B','#902E07','#DD8615')
names(mycols) = c('A','C','L','FH','FM','P','R')
mycols_trans = paste(mycols, '50', sep='')
names(mycols_trans) = c('A','C','L','FH','FM','P','R')
mypch = c(3,3,3) #c(0,4,5)#
mypcols = c('black','grey30','white')
myadj=.15

# Order variables from lowest to highest total effects
total = allsp
ordered_vars = rownames(total[order(total$std.all),])
ordered_vars = ordered_vars[!(ordered_vars %in% c('FH','FM'))] # Drop effects of FM and FH categories (they may be non-sensical)

# Put tables in same order
use_direct = allsp_d[ordered_vars,]
use_indirect = allsp_i[ordered_vars,] # regS and tot_abun will be missing
use_total = total[ordered_vars,]
use_indirectF = allsp_if
use_indirectR = allsp_ir

library(lattice)

# Define range limits that will include 95% confidence intervals
myrange = range(c(use_total[,c('std.ci.lower','std.ci.upper')],
	use_direct[,c('std.ci.lower','std.ci.upper')],
	use_indirect[,c('std.ci.lower','std.ci.upper')]), na.rm=T)+c(-.04, .04)

# Total and direct standardized effects on same graph
png('./Figures/New Coordinates/Standardized direct total indirect via abundance effects on AllSp richness.png', height=1000, width=1400, type="cairo-png")
dotplot(as.numeric(factor(rownames(use_total), levels = ordered_vars))~std.all, data=use_total, 
	xlab=list('Standardized Effect',cex=3), ylab='',
	main='',cex.lab=3,aspect=9/10, xlim=myrange,
	panel=function(x,y){
		
		# Add horizontal boxes
		panel.rect(-2,1:length(ordered_vars)-.5, 2, 1:length(ordered_vars)+.5,
			col=mycols_trans[use_total[ordered_vars,'type']], border='transparent' )

		# Add vertical line at 0
		panel.abline(v=0, col='grey30', lty=2, lwd=2)		
		
		# Add null distribution segments for direct effects
		panel.segments(use_direct$std.ci.lower, y+myadj,
			use_direct$std.ci.upper, y+myadj, 
			col=mypcols[3], lwd=4.5, lend=1)
		# Add points for direct estimated effects
		panel.points(use_direct$std.all, y+myadj, col=mypcols[3], pch=mypch[3], cex=3, lwd=3) 
		
		# Add null distribution segments for indirect effects		
		panel.segments(use_indirect$std.ci.lower, y-myadj,
			use_indirect$std.ci.upper, y-myadj, 
			col=mypcols[2], lwd=4.5, lend=1)	
		# Add points for indirect estimated effects
		panel.points(use_indirect$std.all, y-myadj, col=mypcols[2], pch=mypch[2], cex=3, lwd=3)
		
		# Add null distribution segments for total effects
		panel.segments(use_total$std.ci.lower, y,
			use_total$std.ci.upper, y, 
			col=mypcols[1], lwd=4.5, lend=1)
		# Add points for total estimated effects
		panel.points(x, y, col=mypcols[1], pch=mypch[1], cex=3, lwd=3) 
		
	},
	scales=list(y=list(labels=varnames[ordered_vars,'midName'], 
		cex=3, col=mycols[use_total[ordered_vars,'type']]),
		x=list(cex=3, tick.number=8)),
	key=list(x=1, y=0, corner=c(1,0), lines=list(type='o', pch=mypch, col=mypcols, lwd=3),
		text=list(c('Total effect','Indirect effect\nvia abundance','Direct effect')),
		background='#bbbbbbff', cex=3, divide=1, padding.text=5, border='black')
)
dev.off()







## Compare direct and indirect effects ###
# in order to test whether indirect effects are higher for mean condition variables

# Define data set to plot
use_direct = allsp_d[rownames(allsp_i),]
use_indirect = allsp_i

# Remove points where direct and indirect effects are non-significant
sig_d = sign(use_direct$std.ci.upper)==sign(use_direct$std.ci.lower)
sig_i = sign(use_indirect$std.ci.upper)==sign(use_indirect$std.ci.lower)

use_direct = use_direct[sig_d+sig_i>0,]
use_indirect = use_indirect[sig_d+sig_i>0,]

# Color scheme: http://colorschemedesigner.com/#3341SsYrGvyw0
mycols = c('#000000','#155C8D','#DD4E15','#0E975D','#711392','#DD8615')
names(mycols) = c('A','C','FH','FM','P','R')
mypch = 15:18

pdf('./Figures/New Coordinates/Compare direct and indirect effects richness.pdf', height=8, width=8)
par(mar=c(5,5,1,1))
plot(use_indirect$std.all,use_direct$std.all, xlim=c(-.55, .66), ylim=c(-.55, .66), 
	type='n', xlab='Indirect effect via abundance', ylab='Direct effect', las=1,
	cex.axis=1.5, cex.lab=1.5)
usr = par('usr')
polygon(c(usr[1],0,usr[2],usr[2],0,usr[1]), c(usr[1],0,-usr[2],usr[2],0,-usr[1]),
     col = "grey90", border = NA)
abline(h=0,v=0, lty=1, lwd=2, col='grey30')
abline(0,1, lty=2, lwd=2, col='grey30')
abline(0,-1, lty=2, lwd=2, col='grey30')
arrows(use_indirect$std.all,use_direct$std.ci.lower,use_indirect$std.all,
	use_direct$std.ci.upper,length=.05, code=3, angle=90, lwd=2, col=mycols[use_direct$type])
arrows(use_indirect$std.ci.lower,use_direct$std.all,use_indirect$std.ci.upper,
	use_direct$std.all,length=.05, code=3, angle=90, lwd=2, col=mycols[use_direct$type])
points(use_indirect$std.all, use_direct$std.all, col=mycols[use_direct$type], 
	pch=mypch[factor(use_direct$type)], cex=2)
box()

# Label outliers
outs = c('totalNS','iso','wetness')
text(use_indirect[outs,'std.all'],use_direct[outs, 'std.all'],varnames[outs,'shortName'], 
	pos=c(1,1,1), offset=2, col=mycols[use_direct[outs,'type']], cex=1.5)

dev.off()




## Make plot comparing Total richness, Parmeliaceae and Physciaceae

## Effects table
rownames(parm_d)==rownames(phys_d) #checking
parm_phys_effectstab = data.frame(parm_d[,c('std.all','std.ci.lower','std.ci.upper','type')], phys_d[,c('std.all','std.ci.lower','std.ci.upper','pval_parmphys')])
parm_phys_effectstab$efftype = 'direct'

rownames(parm_i)==rownames(phys_i) #checking
parm_phys_effectstab = rbind(parm_phys_effectstab,
	data.frame(parm_i[,c('std.all','std.ci.lower','std.ci.upper','type')], phys_i[,c('std.all','std.ci.lower','std.ci.upper','pval_parmphys')], efftype='indirect')
)
rownames(parm)==rownames(phys) #checking
parm_phys_effectstab = rbind(parm_phys_effectstab,
	data.frame(parm[,c('std.all','std.ci.lower','std.ci.upper','type')], phys[,c('std.all','std.ci.lower','std.ci.upper','pval_parmphys')], efftype='total')
)

myorder = order(parm_phys_effectstab$type, rownames(parm_phys_effectstab))
parm_phys_effectstab = parm_phys_effectstab[myorder,]
parm_phys_effectstab$varname = varnames[rownames(parm_phys_effectstab), 'displayName']

write.csv(parm_phys_effectstab, './SEM Models/Compare effects Parmeliaceae vs Physciaceae.csv', row.names=F)

options(scipen = 10)
parm_phys_effectstab$parm_ci = paste('(',signif(parm_phys_effectstab$std.ci.lower,digits=2),', ', signif(parm_phys_effectstab$std.ci.upper,digits=2),')', sep='')
parm_phys_effectstab$phys_ci = paste('(',signif(parm_phys_effectstab$std.ci.lower.1,digits=2), ', ', signif(parm_phys_effectstab$std.ci.upper.1,digits=2),')', sep='')

write.csv(parm_phys_effectstab[,c('parm_ci','phys_ci')], './SEM Models/Compare effects Parmeliaceae vs Physciaceae ci.csv', row.names=F)


ordered_vars = rownames(allsp[order(allsp$std.all),])
ordered_vars = ordered_vars[!(ordered_vars %in% c('FH','FM'))] # Drop effects of FM and FH categories (they may be non-sensical)

# put tables in same order
parm = parm[ordered_vars,]
phys = phys[ordered_vars,]
allsp = allsp[ordered_vars,]

mypch = c(3,4,5)#c('I','I','I') #
mypcols = c('black','grey30','white')

myrange = range(c(allsp[,c('std.ci.lower','std.ci.upper')],
	parm[,c('std.ci.lower','std.ci.upper')],
	phys[,c('std.ci.lower','std.ci.upper')]), na.rm=T)+c(-.03, .03)

png('./Figures/New Coordinates/Standardized total effects on richness across families.png', height=800, width=1600, type="cairo-png")
dotplot(1:nrow(allsp)~std.all, data=allsp, 
	xlab=list('Standardized coefficient',cex=3), ylab='',
	main='',cex.lab=3,aspect=4/5, xlim=myrange,
	panel=function(x,y){
		
		# Add horizontal boxes
		panel.rect(-2,1:length(ordered_vars)-.5, 2, 1:length(ordered_vars)+.5,
			col=mycols_trans[factor(allsp[ordered_vars,'type'])], border='transparent' )

		# Add vertical line at 0
		panel.abline(v=0, col='grey30', lty=2, lwd=2)		
		
		# Add null distribution segments for all species
		panel.segments(allsp$std.ci.lower, y,
			allsp$std.ci.upper, y, 
			col=mypcols[1], lwd=4.5, lend=1)
		# Add points for total estimated effects
		panel.points(x, y, col=mypcols[1], pch=mypch[1], cex=3, lwd=3) 
		
		# Add null distribution segments for Physciaceae
		panel.segments(phys$std.ci.lower, y+myadj,
			phys$std.ci.upper, y+myadj, 
			col=mypcols[3], lwd=4.5, lend=1)
		# Add points for Physciaceae
		panel.points(phys$std.all, y+myadj, col=mypcols[3], pch=mypch[3], cex=3, lwd=3) 
		
		# Add null distribution segments for Parmeliaceae		
		panel.segments(parm$std.ci.lower, y-myadj,
			parm$std.ci.upper, y-myadj, 
			col=mypcols[2], lwd=4.5, lend=1)	
		# Add points for indirect estimated effects
		panel.points(parm$std.all, y-myadj, col=mypcols[2], pch=mypch[2], cex=3, lwd=3)
	},
	scales=list(y=list(labels=varnames[ordered_vars,'displayName'], 
		cex=3, col=mycols[factor(allsp[ordered_vars,'type'])]),
		x=list(cex=3, tick.number=8)),
	key=list(x=1, y=.48, corner=c(1,.5), lines=list(type='o', pch=mypch, col=mypcols, lwd=3),
		text=list(c('All Species','Parmeliaceae','Physciaceae')),
		background='#00000033', cex=3, divide=1, padding.text=5)
)
dev.off()







