# This script is used to fit SEM models of lichen FIA data



source('./UNC/Projects/FIA Lichen/GitHub/FIA-Lichens/load_data.R')
source('./GitHub/FIA-Lichens/fia_lichen_analysis_functions.R')

################################################################################
### Path Analysis ###

library(lavaan)
library(semPlot)
library(semTools)

# Full path model with measurement error model on lichen.rich_log
# Previously decided to use unstandardized data
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
	rain_lowRH + wetness ~~ iso + pseas + mat
	mat ~~ iso + pseas
	iso ~~ pseas

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
fit_measerr = sem(path_measerr, data=working_data_unstd_fit, fixed.x=F)

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

	# Climate effects on forest
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
	rain_lowRH + wetness ~~ iso + pseas + mat
	mat ~~ iso + pseas
	iso ~~ pseas
	
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
	
	# Indirect effects of climate on lichen richness
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
endfit =  sem(path_measerr, data=working_data_unstd_test, fixed.x=F, se='robust')
summary(endfit, standardized=T, rsq=TRUE, fit.measures=T)

# Examine model residuals
res = resid(endfit)
which(abs(res$cov)>0.2, arr.ind=T)

save(endfit, file='./SEM models/path_measerr_finalmod_testdata.Rdata')

endfit_ests = parameterEstimates(endfit, standardized=T, boot.ci.type='perc')
subset(endfit_ests, pvalue>0.05)
endfit_paths = subset(endfit_ests, op=='~')
endfit_paths[order(endfit_paths$std.all, decreasing=T),]

write.csv(endfit_ests, './SEM models/finalmod_testdata_estimates.csv', row.names=F)

# Set up data for bootstrapping standardized coefficients
# This will be done on the cluster
endfit_std = bootstrapLavaan(endfit, R=3, FUN=function(x) c(coef(x),standardizedSolution(x)$est.std))

alldata=data.frame(working_data_unstd_test, 
	master[rownames(working_data_unstd_test),c('fric','fdiv','raoQ','Parmeliaceae','Physciaceae')])
write.csv(alldata, './SEM models/Bootstrap/unstandardized_test_dataset.csv', row.names=T)

alldata = read.csv('./SEM models/Bootstrap/unstandardized_test_dataset.csv', row.names=1)
endfit =  sem(path_measerr, data=alldata, fixed.x=F, se='robust')


### See scripts for bootstrap calculations done on the cluster


## Calculate confidence intervals from bootstrapped results
## Uses Rdata files saved from Kure cluster runs
library(lavaan)
load('./SEM Models/Bootstrap/bootstrap_lichen_richness_model.Rdata')

coefrange = (length(coef(endfit))+1):(nrow(standardizedSolution(endfit))+length(coef(endfit)))

#phys_boot = endfit_std[,coefrange] # save results to use later
#parm_boot = endfit_std[,coefrange] # save results to use later

## Test whether parameter values are different between Parm and Phys
# Calculate difference between Parm and Phys estimates for each bootsample
diffs = parm_boot - phys_boot

# Calculate the probability of getting 0 from the distribution of bootsample differences
pvals = apply(diffs, 2, function(x) min(sum(x<0), sum(x>0))/1000*2 ) 
parameterEstimates(endfit)[which(pvals<0.05),'label']

## Create tables of estimates

summary(endfit, standardized=T, rsq=TRUE, fit.measures=T)
endfit_ests = parameterEstimates(endfit)
endfit_ests$std.all = apply(endfit_std[,coefrange], 2, median)
endfit_ests$std.ci.lower = apply(endfit_std[,coefrange], 2, function(x) quantile(x, p=0.025))
endfit_ests$std.ci.upper = apply(endfit_std[,coefrange], 2, function(x) quantile(x, p=0.975))
#endfit_ests$pval_parmphys = pvals
sum(endfit_ests$std.all > endfit_ests$std.ci.upper) # Checking
sum(endfit_ests$std.all < endfit_ests$std.ci.lower) # Checking


# Check whether bootstrapped standardized solutions match the regular standardized solution
plot(apply(endfit_std[,coefrange], 2, median)~standardizedSolution(endfit)$est.std)
abline(0,1,col=2)
hist(apply(endfit_std[,coefrange], 2, median)-standardizedSolution(endfit)$est.std)
hist(apply(endfit_std[,coefrange], 2, mean)-standardizedSolution(endfit)$est.std)


## Make a table of direct effects on lichen richness
## Keep in mind that se for non-standardized coefficients were caluclated using the robust method, not bootstrapped
direct_rich = subset(endfit_ests, (lhs=='lichen_rich')&(op=='~'))
direct_rich = direct_rich[,c('rhs','est','se','z','pvalue','ci.lower','ci.upper','std.all','std.ci.lower','std.ci.upper')] #, 'pval_parmphys', 'pval_parmphys'
names(direct_rich)[1] = 'predictor'
rownames(direct_rich) = direct_rich$predictor
rownames(direct_rich)[rownames(direct_rich)=='regS'] = 'reg' # These change with each dataset
rownames(direct_rich)[rownames(direct_rich)=='tot_abun_log'] = 'abun_log' # These change with each dataset

## Make a table of total effects on lichen richness
total = endfit_ests[grep('TE',endfit_ests$label),]
total = total[,c('lhs','est','se','z','pvalue','ci.lower','ci.upper','std.all','std.ci.lower','std.ci.upper')]#, 'pval_parmphys', 'pval_parmphys'
names(total)[1] = 'predictor'
total$predictor = substring(total$predictor, first=4)
addrows = subset(direct_rich, !(predictor %in% total$predictor))
colnames(addrows) = colnames(total)
total = rbind(total, addrows)
rownames(total) = total$predictor
rownames(total)[rownames(total)=='regS'] = 'reg' # These change with each dataset
rownames(total)[rownames(total)=='tot_abun_log'] = 'abun_log' # These change with each dataset; tot, parm, phys

## Make a table of indirect effects on lichen richness via abundance
c(paste('IE',c('bark_moist_pct.rao.ba','wood_SG.rao.ba','LogSeed.rao.ba','PIE.ba.tree','propDead','lightDist.mean','diamDiversity',
		'bark_moist_pct.ba','wood_SG.ba','LogSeed.ba','totalCirc','bigTrees','light.mean','totalNS'), sep='_'),
	paste('IE',c('wetness','rain_lowRH','iso','pseas','mat'),'A', sep='_'))->indir_vars
indirect = subset(endfit_ests, label %in% indir_vars)
indirect = indirect[,c('lhs','est','se','z','pvalue','ci.lower','ci.upper','std.all','std.ci.lower','std.ci.upper')]#, 'pval_parmphys', 'pval_parmphys'
names(indirect)[1] = 'predictor'
indirect$predictor = substring(indirect$predictor, first=4)
indirect$predictor = sapply( indirect$predictor, function(x) ifelse(substr(x, nchar(x), nchar(x))=='A', substr(x, 1, nchar(x)-2), x))	
rownames(indirect) = indirect$predictor

## Make a table of direct effects on lichen abundance
direct_abun = subset(endfit_ests, (lhs=='tot_abun_log')&(op=='~')) # This changes with each model; tot, parm, phys
direct_abun = direct_abun[,c('rhs','est','se','z','pvalue','ci.lower','ci.upper','std.all','std.ci.lower','std.ci.upper')]#,'pval_parmphys' ,'pval_parmphys'
names(direct_abun)[1] = 'predictor'
rownames(direct_abun) = direct_abun$predictor

## Append variable types column
vartypes = read.csv('./SEM models/var_types.csv', row.names=1)

total$type = vartypes[rownames(total),'type']
direct_rich$type = vartypes[rownames(direct_rich),'type']
direct_abun$type = vartypes[rownames(direct_abun),'type']
indirect$type = vartypes[rownames(indirect),'type']

## Save tables- make sure to changes names to appropriate dataset
write.csv(total, './SEM models/AllSp_testdata_totaleffects.csv', row.names=T)
write.csv(direct_rich, './SEM models/AllSp_testdata_directeffects_richness.csv', row.names=T)
write.csv(direct_abun, './SEM models/AllSp_testdata_directeffects_abundance.csv', row.names=T)
write.csv(indirect, './SEM models/AllSp_testdata_indirecteffects_via_abundance.csv', row.names=T)

### Figures ###

## Compare direct effects on richness vs abundance
use_direct = direct_rich[rownames(direct_abun),]

mycols = c('dodgerblue','darkred','forestgreen','purple')
mypch = 15:18

png('./Figures/New Coordinates/Compare direct effects on Physciaceae richness vs abunance.png', height=500, width=500)
par(mar=c(5,5,1,1))
plot(as.numeric(as.matrix(direct_abun[,c('std.ci.lower','std.ci.upper')])),
	as.numeric(as.matrix(use_direct[,c('std.ci.lower','std.ci.upper')])), 
	type='n', xlab='Effect on lichen abundance', ylab='Effect on lichen richness', las=1)
abline(h=0,v=0, lty=3, lwd=3, col='grey30')
arrows(direct_abun$std.all,use_direct$std.ci.lower,direct_abun$std.all,
	use_direct$std.ci.upper,length=.05, code=3, angle=90, lwd=2, col=mycols[factor(use_direct$type)])
arrows(direct_abun$std.ci.lower,use_direct$std.all,direct_abun$std.ci.upper,
	use_direct$std.all,length=.05, code=3, angle=90, lwd=2, col=mycols[factor(use_direct$type)])
points(direct_abun$std.all, use_direct$std.all, col=mycols[factor(use_direct$type)], 
	pch=mypch[factor(use_direct$type)], cex=2)
dev.off()

## Compare direct effects vs. indirect effects via abundance vs. total effects
ordered_vars = rownames(total[order(total$std.all),])
ordered_vars = ordered_vars[!(ordered_vars %in% c('FH','FM'))] # Drop effects of FM and FH categories (they may be non-sensical)

# put tables in same order
use_direct = direct_rich[ordered_vars,]
use_indirect = indirect[ordered_vars,]
use_total = total[ordered_vars,]

library(lattice)
#library(gplots)

# Color scheme: http://colorschemedesigner.com/#3341SsYrGvyw0
mycols = c('#000000','#07395D','#902E07','#05633B','#711392','#DD8615')
names(mycols) = c('A','C','FH','FM','P','R')
mycols_trans = paste(mycols, '66', sep='')
mypch = c('I','I','I') #c(0,4,5)#
mypcols = c('black','grey30','white')
myadj=.15

myrange = range(c(use_total[,c('std.ci.lower','std.ci.upper')],
	use_direct[,c('std.ci.lower','std.ci.upper')],
	use_indirect[,c('std.ci.lower','std.ci.upper')]), na.rm=T)+c(-.04, .04)


# Total and direct standardized effects on same graph
png('./Figures/New Coordinates/Standardized direct total indirect via abundance effects on Physciaceae richness.png', height=800, width=1400, type="cairo-png")
dotplot(as.numeric(factor(rownames(use_total), levels = ordered_vars))~std.all, data=use_total, 
	xlab=list('Standardized coefficient',cex=3), ylab='',
	main='',cex.lab=3,aspect=4/5, xlim=myrange,
	panel=function(x,y){
		
		# Add horizontal boxes
		panel.rect(-2,1:length(ordered_vars)-.5, 2, 1:length(ordered_vars)+.5,
			col=mycols_trans[factor(total[ordered_vars,'type'])], border='transparent' )

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
		cex=3, col=mycols[factor(total[ordered_vars,'type'])]),
		x=list(cex=3, tick.number=8)),
	key=list(x=1, y=.48, corner=c(1,.5), lines=list(type='o', pch=mypch, col=mypcols, lwd=3),
		text=list(c('Total effect','Indirect effect\nvia abundance','Direct effect')),
		background='#00000033', cex=3, divide=1, padding.text=5)
)
dev.off()

## Make plot comparing Total richness, Parmeliaceae and Physciaceae

# Using total effects
parm = read.csv('./SEM models/Parmeliaceae_testdata_totaleffects.csv', row.names=1)
phys = read.csv('./SEM models/Physciaceae_testdata_totaleffects.csv', row.names=1)
allsp = read.csv('./SEM models/AllSp_testdata_totaleffects.csv', row.names=1)

parm_i = read.csv('./SEM models/Parmeliaceae_testdata_indirecteffects_via_abundance.csv', row.names=1)
phys_i = read.csv('./SEM models/Physciaceae_testdata_indirecteffects_via_abundance.csv', row.names=1)
allsp_i = read.csv('./SEM models/AllSp_testdata_indirecteffects_via_abundance.csv', row.names=1)

parm_d = read.csv('./SEM models/Parmeliaceae_testdata_directeffects_richness.csv', row.names=1)
phys_d = read.csv('./SEM models/Physciaceae_testdata_directeffects_richness.csv', row.names=1)
allsp_d = read.csv('./SEM models/AllSp_testdata_directeffects_richness.csv', row.names=1)

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


### Compare effects on functional diversity versus richness ###
load('./SEM Models/Bootstrap/bootstrap_fric_model.Rdata')
load('./SEM Models/Bootstrap/bootstrap_raoQ_model.Rdata')
load('./SEM Models/Bootstrap/bootstrap_lichen_richness_model.Rdata')

coefrange = (length(coef(endfit))+1):(nrow(standardizedSolution(endfit))+length(coef(endfit)))

# Make table of all parameter estimates
endfit_ests = parameterEstimates(endfit)
endfit_ests$bootrow = 1:nrow(endfit_ests)
endfit_ests$std.all = apply(endfit_std[,coefrange], 2, median)
endfit_ests$std.ci.lower = apply(endfit_std[,coefrange], 2, function(x) quantile(x, p=0.025)) 
endfit_ests$std.ci.upper = apply(endfit_std[,coefrange], 2, function(x) quantile(x, p=0.975)) 
sum(endfit_ests$std.all > endfit_ests$std.ci.upper) # Checking
sum(endfit_ests$std.all < endfit_ests$std.ci.lower) # Checking

fric_est = endfit_ests
raoq_est = endfit_ests
rich_est = endfit_ests

fric_boot = endfit_std[,coefrange] # save results to use later
raoq_boot = endfit_std[,coefrange] # save results to use later
rich_boot = endfit_std[,coefrange]

# Make a table of direct effects on lichen richness
# Keep in mind that se for non-standardized coefficients were caluclated using the robust method, not bootstrapped
# Change lhs== to 'fric', 'raoQ', or 'lichen_rich' as necessary
direct_rich = subset(endfit_ests, (lhs=='lichen_rich')&(op=='~'))
direct_rich = direct_rich[,c('rhs','est','se','z','pvalue','ci.lower','ci.upper','std.all','std.ci.lower','std.ci.upper','bootrow')]
names(direct_rich)[1] = 'predictor'
rownames(direct_rich) = direct_rich$predictor
rownames(direct_rich)[rownames(direct_rich)=='regS'] = 'reg' # These change with each dataset
rownames(direct_rich)[rownames(direct_rich)=='tot_abun_log'] = 'abun_log' # These change with each dataset

# Make a table of total effects on lichen richness
total = endfit_ests[grep('TE',endfit_ests$label),]
total = total[,c('lhs','est','se','z','pvalue','ci.lower','ci.upper','std.all','std.ci.lower','std.ci.upper','bootrow')]
names(total)[1] = 'predictor'
total$predictor = substring(total$predictor, first=4)
addrows = subset(direct_rich, !(predictor %in% total$predictor))
colnames(addrows) = colnames(total)
total = rbind(total, addrows)
rownames(total) = total$predictor
rownames(total)[rownames(total)=='regS'] = 'reg' # These change with each dataset
rownames(total)[rownames(total)=='tot_abun_log'] = 'abun_log' # These change with each dataset

# Make a table of indirect effects on lichen richness via abundance
c(paste('IE',c('bark_moist_pct.rao.ba','wood_SG.rao.ba','LogSeed.rao.ba','PIE.ba.tree','propDead','lightDist.mean','diamDiversity',
		'bark_moist_pct.ba','wood_SG.ba','LogSeed.ba','totalCirc','bigTrees','light.mean','totalNS'), sep='_'),
	paste('IE',c('wetness','rain_lowRH','iso','pseas','mat'),'A', sep='_'))->indir_vars
indirect = subset(endfit_ests, label %in% indir_vars)
indirect = indirect[,c('lhs','est','se','z','pvalue','ci.lower','ci.upper','std.all','std.ci.lower','std.ci.upper','bootrow')]
names(indirect)[1] = 'predictor'
indirect$predictor = substring(indirect$predictor, first=4)
indirect$predictor = sapply( indirect$predictor, function(x) ifelse(substr(x, nchar(x), nchar(x))=='A', substr(x, 1, nchar(x)-2), x))	
rownames(indirect) = indirect$predictor

# Make a table of direct effects on lichen abundance
direct_abun = subset(endfit_ests, (lhs=='tot_abun_log')&(op=='~'))
direct_abun = direct_abun[,c('rhs','est','se','z','pvalue','ci.lower','ci.upper','std.all','std.ci.lower','std.ci.upper','bootrow')]
names(direct_abun)[1] = 'predictor'
rownames(direct_abun) = direct_abun$predictor

# Add variable indicating type
vartypes = read.csv('./SEM models/var_types.csv', row.names=1)
total$type = vartypes[rownames(total),'type']
direct_rich$type = vartypes[rownames(direct_rich),'type']
direct_abun$type = vartypes[rownames(direct_abun),'type']
indirect$type = vartypes[rownames(indirect),'type']

# Save tables
write.csv(total, './SEM models/raoQ_testdata_totaleffects.csv', row.names=T)
write.csv(direct_rich, './SEM models/raoQ_testdata_directeffects_richness.csv', row.names=T)
write.csv(direct_abun, './SEM models/raoQ_testdata_directeffects_abundance.csv', row.names=T)
write.csv(indirect, './SEM models/raoQ_testdata_indirecteffects_via_abundance.csv', row.names=T)

# Rename for future use
fric = total
fric_i = indirect
fric_d = direct_rich
raoq = total
raoq_i = indirect
raoq_d = direct_rich
allsp = total
allsp_i = indirect
allsp_d = direct_rich

## Calculate P-values for difference with lichen richness model

# Make one effects table for each metric
fric$efftype = 'total'
fric_d$efftype = 'direct'
fric_i$efftype = 'indirect'
fric_all = rbind(fric, fric_d, fric_i)
raoq$efftype = 'total'
raoq_d$efftype = 'direct'
raoq_i$efftype = 'indirect'
raoq_all = rbind(raoq, raoq_d, raoq_i)
allsp$efftype = 'total'
allsp_d$efftype = 'direct'
allsp_i$efftype = 'indirect'
allsp_all = rbind(allsp, allsp_d, allsp_i)

# Put all tables in same order
fric_all = fric_all[order(fric_all$predictor, fric_all$efftype),]
raoq_all = raoq_all[order(raoq_all$predictor, raoq_all$efftype),]
allsp_all = allsp_all[order(allsp_all$predictor, allsp_all$efftype),]
cbind(fric_all[,c('predictor','efftype')], raoq_all[,c('predictor','efftype')], allsp_all[,c('predictor','efftype')])

# Calculate difference between fric or raoq and lichen richness estimates in bootsamples
diffs_fric = fric_boot[,fric_all$bootrow] - rich_boot[,allsp_all$bootrow]
diffs_raoq = raoq_boot[,fric_all$bootrow] - rich_boot[,allsp_all$bootrow]

# Calculate the probability of getting 0 from the distribution of bootsample differences
fric_all$pval_rich = apply(diffs_fric, 2, function(x) min(sum(x<0), sum(x>0))/1000*2 ) 
raoq_all$pval_rich = apply(diffs_raoq, 2, function(x) min(sum(x<0), sum(x>0))/1000*2 ) 

# Save
write.csv(fric_all, './SEM Models/Compare effects fric vs lichen richness.csv', row.names=F)
write.csv(raoq_all, './SEM Models/Compare effects raoQ vs lichen richness.csv', row.names=F)

## Examine differences
i = which(fric_all$pval_rich<0.05)
fric_all[i, c('predictor','std.all','std.ci.lower','std.ci.upper','efftype','pval_rich')]
allsp_all[i, c('predictor','std.all','std.ci.lower','std.ci.upper','efftype')]
# See notebook for results
subset(fric_all, type=='FH')

i = which(raoq_all$pval_rich<0.05)
raoq_all[i, c('predictor','std.all','std.ci.lower','std.ci.upper','efftype','pval_rich')]
allsp_all[i, c('predictor','std.all','std.ci.lower','std.ci.upper','efftype')]
# See notebook for results
subset(fric_all, type=='FH')


### Compare direct and indirect effects ###
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



