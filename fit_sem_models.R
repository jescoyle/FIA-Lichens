# This script is used to fit SEM models of lichen FIA data

source('./UNC/Projects/FIA Lichen/GitHub/FIA-Lichens/load_data.R')
source('./GitHub/FIA-Lichens/fia_lichen_analysis_functions.R')

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
path_measerr = "

	# Latent variables
	lichen_rich =~ sqrt(0.75)*lichen.rich_log

	# Local environment/climate effects on forest structure
	wood_SG.ba + wood_SG.rao.ba + bark_moist_pct.rao.ba + LogSeed.rao.ba + PIE.ba.tree + propDead + lightDist.mean + diamDiversity + bark_moist_pct.ba + LogSeed.ba + bigTrees + light.mean + PC1 ~
		wetness + rain_lowRH + iso + pseas + mat + radiation

	# Regional climate, pollution, and regional forest heterogeneity effects on regional richness
	regS ~ wetness_reg_mean + rain_lowRH_reg_mean + iso_reg_mean + pseas_reg_mean + mat_reg_mean +
		wetness_reg_var + rain_lowRH_reg_var + iso_reg_var + pseas_reg_var + mat_reg_var +
		regS_tree + totalNS_reg

	# Regional climate effects on regional forest heterogeneity
	regS_tree ~ wetness_reg_mean + rain_lowRH_reg_mean + iso_reg_mean + pseas_reg_mean + mat_reg_mean +
		wetness_reg_var + rain_lowRH_reg_var + iso_reg_var + pseas_reg_var + mat_reg_var

	# Local climate effects on pollution
	totalNS ~ wetness + rain_lowRH

	# Regional climate effects on pollution
	totalNS_reg ~ wetness_reg_mean + rain_lowRH_reg_mean

	# Effects on local lichen richness
	lichen_rich ~ bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba +
		PIE.ba.tree + propDead + lightDist.mean + diamDiversity +
		bark_moist_pct.ba + wood_SG.ba + LogSeed.ba + bigTrees + light.mean + PC1 +
		wetness + rain_lowRH + iso + pseas + mat + radiation + regS + totalNS +
		wetness_reg_mean + rain_lowRH_reg_mean + iso_reg_mean + pseas_reg_mean + mat_reg_mean +
		wetness_reg_var + rain_lowRH_reg_var + iso_reg_var + pseas_reg_var + mat_reg_var +
		totalNS_reg + regS_tree + tot_abun_log
	
	# Effects on local lichen abundance
	tot_abun_log ~ bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba +
		PIE.ba.tree + propDead + lightDist.mean + diamDiversity +
		bark_moist_pct.ba + wood_SG.ba + LogSeed.ba + bigTrees + light.mean + PC1 +
		wetness + rain_lowRH + iso + pseas + mat + radiation + totalNS +
		wetness_reg_mean + rain_lowRH_reg_mean + iso_reg_mean + pseas_reg_mean + mat_reg_mean +
		wetness_reg_var + rain_lowRH_reg_var + iso_reg_var + pseas_reg_var + mat_reg_var +
		totalNS_reg + regS_tree

	## Covariances between exogenous and endogenous variables
	# Don't need to specificy for Climate L-R b/c these are exogenous and are calculated automatically in lavaan

	# Local-regional pollution
	totalNS ~~ totalNS_reg
	
	# Local-regional forest structure
	regS_tree ~~ PIE.ba.tree

	## Covariances among endogenous predictors in the same group
	# Don't need to specificy for Climate LO b/c these are exogenous and are calculated automatically in lavaan
	
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
fit_fixed = sem(path_measerr, data=working_data_fit, fixed.x=T, estimator='ML', se='robust')
fit_free = sem(path_measerr, data=working_data_fit, fixed.x=F, estimator='ML', se='robust')

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

png('./Figures/sem residual correlations.png', height=1700, width=1700, type='cairo')
corrplot(cortab, method='square', type='upper', diag=F, 
	order='original', hclust.method='complete', p.mat=cortabsig,
	sig.level=.99, insig='blank', tl.cex=1.5, tl.col=1, cl.cex=2, mar=c(1,1,4,1))
dev.off()

which(abs(res$cov[lower.tri(res$cov)])>.2, arr.ind=T)

semPaths(fit_measerr, 'std')

## Model without abundance
# See sem_boot_noabun_allsp.R script

## Define models 
# Rules for coefficient names:
# dependent variable first, independent variable second
# Capital = regional, lowercase = local
# R = richness, A = abundance, F = forest, C = climate/environment, P = pollution
# H = heterogeneity, M = optimality (only used for F and C)
# number indicates which variable in the category: e.g. fh3 is LogSeed.rao.ba


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


### TABLE NAMES HAVE CHANGED!!! ###


## Read in tables of parameter estimates and effects
# All parameter estimares
allsp_ests = read.csv('./SEM models/AllSp_testdata_parameterEstimates.csv')
parm_ests = read.csv('./SEM models/Parm_testdata_parameterEstimates.csv')
phys_ests = read.csv('./SEM models/Phys_testdata_parameterEstimates.csv')

# Total effects
allsp = read.csv('./SEM models/AllSp_testdata_totaleffects.csv', row.names=1)
parm = read.csv('./SEM models/Parm_testdata_totaleffects.csv', row.names=1)
phys = read.csv('./SEM models/Phys_testdata_totaleffects.csv', row.names=1)
fric = read.csv('./SEM models/Fric_testdata_totaleffects.csv', row.names=1)
raoq = read.csv('./SEM models/RaoQ_testdata_totaleffects.csv', row.names=1)

# Direct effects
allsp_d = read.csv('./SEM models/AllSp_testdata_directeffects_richness.csv', row.names=1)
parm_d = read.csv('./SEM models/Parm_testdata_directeffects_richness.csv', row.names=1)
phys_d = read.csv('./SEM models/Phys_testdata_directeffects_richness.csv', row.names=1)
fric_d = read.csv('./SEM models/Fric_testdata_directeffects_richness.csv', row.names=1)
raoq_d = read.csv('./SEM models/RaoQ_testdata_directeffects_richness.csv', row.names=1)

# Direct effect on abundance
allsp_da = read.csv('./SEM models/AllSp_testdata_directeffects_abundance.csv', row.names=1)
parm_da = read.csv('./SEM models/Parm_testdata_directeffects_abundance.csv', row.names=1)
phys_da = read.csv('./SEM models/Phys_testdata_directeffects_abundance.csv', row.names=1)

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
names(which(apply(parm[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0))
names(which(apply(phys[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0))

tot_sig = allsp[which(apply(allsp[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0),]
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

#######################################################################################
### Figures ###


######################################
## Compare direct effects on richness vs abundance

## Compare direct effects vs. total effects

# Order variables from lowest to highest total effects
total = allsp # This changes based on what response variable is being analyzed
ordered_vars = rownames(total[order(total$std.all),])
ordered_vars = ordered_vars[!(ordered_vars %in% c('FH','FM'))] # Drop effects of FM and FH categories (they may be non-sensical)

# Put tables in same order
use_direct = allsp_d[ordered_vars,] # This changes based on what response variable is being analyzed
use_total = total[ordered_vars,]


library(lattice)

# Define range limits that will include 95% confidence intervals
myrange = range(c(use_total[,c('std.ci.lower','std.ci.upper')],
	use_direct[,c('std.ci.lower','std.ci.upper')]), na.rm=T)+c(-.04, .04)
myrange[1] = -.8

# Color scheme: http://colorschemedesigner.com/#3341SsYrGvyw0
#mycols = c('#000000','#07395D','#066788','#688e60','#05633B','#902E07','#DD8615')
#names(mycols) = c('A','C','L','FH','FM','P','R')
mycols = c('#2415B0','#00BF32')
mycolsbw = c('grey80','white')
names(mycols) = c('R','L')
names(mycolsbw) = c('R','L')
mycols_trans = paste(mycols, '50', sep='')
names(mycols_trans) = c('R','L')
plot(1:length(mycols),1:length(mycols), type='n'); text(1:length(mycols),1:length(mycols),labels=names(mycols), cex=2, col=mycols) # Check colors
#names(mycols_trans) = c('A','C','L','FH','FM','P','R')
mypch = c(22,23)
mypcols = c('white','grey30')
myadj=.15
mytypes = expression('C'['R'],'C'['L'],'F'['H'],'F'['O'],'C'['R'],'R','') # symbols used in plot to denote variable types
names(mytypes)=c('C','L','FH','FM','P','R','A')

# Total and direct standardized effects on same graph
svg('./Figures/Standardized direct total effects on AllSp richness bw.svg', height=12, width=19)
dotplot(as.numeric(factor(rownames(use_total), levels = ordered_vars))~std.all, data=use_total, 
	xlab=list('Standardized Effect',cex=3), ylab='',
	main='',cex.lab=3,aspect=9/10, xlim=myrange,
	panel=function(x,y){
	
		# Add horizontal boxes
		panel.rect(-2,1:length(ordered_vars)-.5, 2, 1:length(ordered_vars)+.5,
			col=mycolsbw[predtypes[ordered_vars,'scale']], border='grey50')
		
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
		panel.text(-.75, y, labels=mytypes[predtypes[ordered_vars,'type']], cex=2)
		
	},
	scales=list(y=list(labels=varnames[ordered_vars,'midName'], 
		cex=3, col='black'),
		x=list(cex=3, tick.number=8)),
	key=list(x=1, y=0, corner=c(1,0), lines=list(type='o', pch=mypch, fill=mypcols, lwd=3, pt.lwd=3),
		text=list(c('Total effect','Direct effect')),
		background='#FFFFFFCC', cex=3, divide=1, padding.text=5, border='black')
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

indirect = phys_i
indirectR = phys_ir
indirectF = phys_if
total = phys
direct = phys_d

use_total_sub = subset(total, type %in% c('C','L'))
use_direct_sub = subset(direct, type %in% c('C','L'))
use_indirect_sub = subset(indirect, type %in% c('C','L'))

order_clim = c('iso','rain_lowRH','radiation','mat','pseas','wetness')

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

# Only look at estimates corresponding to paths
paths = subset(phys_ests, op=='~') # Change based on response variables of interest

# Calculate which paths are significant
sigpaths = paths[which(apply(paths[,c('std.ci.lower','std.ci.upper')],1,prod)>0),]

# Histogram of significant paths
hist(sigpaths$std.all)
sigpaths$sigcat = cut(abs(sigpaths$std.all), c(0,.2,.4,.6,.8,1))

# Find all variables in the model
allvars = unique(c(paths$rhs,paths$lhs))
allvars[allvars=='regPhys'] <-'reg' # Change based on response variables of interest
allvars[allvars=='phys_abun_log'] <- 'abun_log' # Change based on response variables of interest

# Replace varnames in sigpaths to generic
sigpaths$lhs[sigpaths$lhs=='phys_abun_log']<-'abun_log'
sigpaths$lhs[sigpaths$lhs=='regPhys']<-'reg'
sigpaths$rhs[sigpaths$rhs=='phys_abun_log']<-'abun_log'
sigpaths$rhs[sigpaths$rhs=='regPhys']<-'reg'

# Find all variables involved in significant relationships
sigvars  = unique(c(sigpaths$lhs, sigpaths$rhs))

## Create a matrix of variable locations
var_locs = matrix(0, nrow=length(allvars), ncol=2, byrow=T)
colnames(var_locs) = c('X','Y')
rownames(var_locs) = allvars

unit = 1

# Start with richness and abundance in center
var_locs['abun_log',] = c(-1,0) 
var_locs['lichen_rich',] = c(1,0)

# Add pollution and regS to corners
var_locs['totalNS',] = c(-2.5,2.5)
var_locs['reg',] = c(2.5,2.5) 

# Add climate and local environment vars to top and bottom
c_vars = rownames(subset(predtypes[allvars,], type=='C'))
var_locs[c_vars,] = cbind(3/(length(c_vars)-1)*(0:(length(c_vars)-1))-1.5,rep(3,length(c_vars)))
var_locs['radiation',] = c(0,-2.5)

# Add forest vars to sides
fh_vars = rownames(subset(predtypes[allvars,], type=='FH'))
fm_vars = rownames(subset(predtypes[allvars,], type=='FM'))

var_locs[fh_vars,]=cbind(rep(3,length(fh_vars)), 3.5/(length(fh_vars)-1)*(0:(length(fh_vars)-1))-2)
var_locs[fm_vars,]=cbind(rep(-3,length(fm_vars)), 3.5/(length(fm_vars)-1)*(0:(length(fm_vars)-1))-2)

var_locs=data.frame(var_locs)

# Make plot
mycex=1.5

svg('./Figures/full model significant paths Phys.svg', height=5, width=8.5) # Change based on response variable of interest
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
text(var_locs[c('lichen_rich','abun_log'),], labels=c('Local\nrichness','Abundance'),pos=1, offset=1)
text(var_locs['radiation',], labels='Solar radiation', pos=1, offset=.25)
text(var_locs['totalNS',], labels=varnames['totalNS','midName'], pos=2, offset=.25)
text(var_locs['reg',], labels='Regional richness', pos=4, offset=.25)

legend(x=-5.5, y=5, xjust=0, yjust=1, c('0.0-0.2','0.2-0.4','0.4-0.6','0.6-0.8','0.8-1.0'), lwd=mycex*(1:5), 
	bty='n')
legend(x=-6, y=5, xjust=0, yjust=1, rep('', 5), lwd=mycex*(1:5), col='red', bty='n')
text(-5.5, 5, '+', cex=2, adj=-1)
text(-6, 5, '-', cex=2, adj=-2.3)

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








	## Indirect Effects	
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
	#IE_totalCirc := AFM4*A
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
		FM3a*(FM3 + IE_LogSeed.ba) + FM5a*(FM5 + IE_bigTrees) + 
		FM6a*(FM6 + IE_light.mean) + FM7a*(FM7 + IE_PC1)
	IE_wetness_R := R1*R
	IE_wetness_P := P1*(P + IE_totalNS)
	IE_wetness_A := AC1*A
	IE_wetness := IE_wetness_FH + IE_wetness_FM + IE_wetness_R + IE_wetness_P + IE_wetness_A

	IE_rain_lowRH_FH := FH1b*(FH1 + IE_bark_moist_pct.rao.ba) + FH2b*(FH2 + IE_wood_SG.rao.ba) +
		FH3b*(FH3 + IE_LogSeed.rao.ba) + FH4b*(FH4 + IE_PIE.ba.tree) + FH5b*(FH5 + IE_propDead) + 
		FH6b*(FH6 + IE_lightDist.mean) + FH7b*(FH7 + IE_diamDiversity)
	IE_rain_lowRH_FM := FM1b*(FM1 + IE_bark_moist_pct.ba) + FM2b*(FM2 + IE_wood_SG.ba) + 
		FM3b*(FM3 + IE_LogSeed.ba) + FM5b*(FM5 + IE_bigTrees) + 
		FM6b*(FM6 + IE_light.mean) + FM7b*(FM7 + IE_PC1)
	IE_rain_lowRH_R := R2*R
	IE_rain_lowRH_P := P2*(P + IE_totalNS)
	IE_rain_lowRH_A := AC2*A
	IE_rain_lowRH := IE_rain_lowRH_FH + IE_rain_lowRH_FM + IE_rain_lowRH_R + IE_rain_lowRH_P + IE_rain_lowRH_A

	IE_iso_FH := FH1c*(FH1 + IE_bark_moist_pct.rao.ba) + FH2c*(FH2 + IE_wood_SG.rao.ba) +
		FH3c*(FH3 + IE_LogSeed.rao.ba) + FH4c*(FH4 + IE_PIE.ba.tree) + FH5c*(FH5 + IE_propDead) + 
		FH6c*(FH6 + IE_lightDist.mean) + FH7c*(FH7 + IE_diamDiversity)
	IE_iso_FM := FM1c*(FM1 + IE_bark_moist_pct.ba) + FM2c*(FM2 + IE_wood_SG.ba) + 
		FM3c*(FM3 + IE_LogSeed.ba) + FM5c*(FM5 + IE_bigTrees) +
		FM6c*(FM6 + IE_light.mean) +  FM7c*(FM7 + IE_PC1)
	IE_iso_R := R3*R
	IE_iso_A := AC3*A
	IE_iso := IE_iso_FH + IE_iso_FM + IE_iso_R + IE_iso_A

	IE_pseas_FH := FH1d*(FH1 + IE_bark_moist_pct.rao.ba) + FH2d*(FH2 + IE_wood_SG.rao.ba) +
		FH3d*(FH3 + IE_LogSeed.rao.ba) + FH4d*(FH4 + IE_PIE.ba.tree) + FH5d*(FH5 + IE_propDead) + 
		FH6d*(FH6 + IE_lightDist.mean) + FH7d*(FH7 + IE_diamDiversity)
	IE_pseas_FM := FM1d*(FM1 + IE_bark_moist_pct.ba) + FM2d*(FM2 + IE_wood_SG.ba) + 
		FM3d*(FM3 + IE_LogSeed.ba) + FM5d*(FM5 + IE_bigTrees) + 
		FM6d*(FM6 + IE_light.mean) + FM7d*(FM7 + IE_PC1)
	IE_pseas_R := R4*R
	IE_pseas_A := AC4*A
	IE_pseas := IE_pseas_FH + IE_pseas_FM + IE_pseas_R + IE_pseas_A
	
	IE_mat_FH := FH1e*(FH1 + IE_bark_moist_pct.rao.ba) + FH2e*(FH2 + IE_wood_SG.rao.ba) +
		FH3e*(FH3 + IE_LogSeed.rao.ba) + FH4e*(FH4 + IE_PIE.ba.tree) + FH5e*(FH5 + IE_propDead) + 
		FH6e*(FH6 + IE_lightDist.mean) + FH7e*(FH7 + IE_diamDiversity)
	IE_mat_FM := FM1e*(FM1 + IE_bark_moist_pct.ba) + FM2e*(FM2 + IE_wood_SG.ba) + 
		FM3e*(FM3 + IE_LogSeed.ba) + FM5e*(FM5 + IE_bigTrees) + 
		FM6e*(FM6 + IE_light.mean) + FM7e*(FM7 + IE_PC1)
	IE_mat_R := R5*R
	IE_mat_A := AC5*A
	IE_mat := IE_mat_FH + IE_mat_FM + IE_mat_R + IE_mat_A

	IE_radiation_FH := FH1f*(FH1 + IE_bark_moist_pct.rao.ba) + FH2f*(FH2 + IE_wood_SG.rao.ba) +
		FH3f*(FH3 + IE_LogSeed.rao.ba) + FH4f*(FH4 + IE_PIE.ba.tree) + FH5f*(FH5 + IE_propDead) + 
		FH6f*(FH6 + IE_lightDist.mean) + FH7f*(FH7 + IE_diamDiversity)
	IE_radiation_FM := FM1f*(FM1 + IE_bark_moist_pct.ba) + FM2f*(FM2 + IE_wood_SG.ba) + 
		FM3f*(FM3 + IE_LogSeed.ba) + FM5f*(FM5 + IE_bigTrees) + 
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
	#TE_totalCirc := IE_totalCirc + FM4
	TE_bigTrees := IE_bigTrees + FM5
	TE_light.mean := IE_light.mean + FM6
	TE_PC1 := IE_PC1 +  FM7	

	TE_FH := TE_bark_moist_pct.rao.ba +	TE_wood_SG.rao.ba + TE_LogSeed.rao.ba +
		TE_PIE.ba.tree + TE_propDead + TE_lightDist.mean + TE_diamDiversity
	TE_FM := TE_bark_moist_pct.ba + TE_wood_SG.ba + TE_LogSeed.ba +
		TE_bigTrees + TE_light.mean + TE_PC1

	TE_totalNS := IE_totalNS + P



