# This script bootstraps parameter estimates for the finalmod model fit to all species richness.
# No direct paths from regional variables to local richness


library(lavaan, lib.loc="/nas02/home/j/r/jrcoyle/Rlibs/")

options(stringsAsFactors=F)

working_soil_test = read.csv('working_soil_test.csv', row.names=1)
predtypes = read.csv('predictors.csv', row.names=1)

## Model with paths from regional variables to local richness
path_soilmod_nopol = "


	# Local environment/climate effects on forest structure
	bark_moist_pct.rao.ba ~ fh1cm1*wetness + fh1cm2*rain_lowRH + fh1cm3*iso + fh1cm4*pseas + fh1cm5*mat + fh1cm6*radiation + fh1sm1*soilPC1 + fh1sm2*soilPC2
	wood_SG.rao.ba ~ fh2cm1*wetness + fh2cm2*rain_lowRH + fh2cm3*iso + fh2cm4*pseas + fh2cm5*mat + fh2cm6*radiation + fh2sm1*soilPC1 + fh2sm2*soilPC2
	LogSeed.rao.ba ~ fh3cm1*wetness + fh3cm2*rain_lowRH + fh3cm3*iso + fh3cm4*pseas + fh3cm5*mat + fh3cm6*radiation + fh3sm1*soilPC1 + fh3sm2*soilPC2
	PIE.ba.tree ~ fh4cm1*wetness + fh4cm2*rain_lowRH + fh4cm3*iso + fh4cm4*pseas + fh4cm5*mat + fh4cm6*radiation + fh4sm1*soilPC1 + fh4sm2*soilPC2
	propDead ~ fh5cm1*wetness + fh5cm2*rain_lowRH + fh5cm3*iso + fh5cm4*pseas + fh5cm5*mat + fh5cm6*radiation + fh5sm1*soilPC1 + fh5sm2*soilPC2
	lightDist.mean ~ fh6cm1*wetness + fh6cm2*rain_lowRH + fh6cm3*iso + fh6cm4*pseas + fh6cm5*mat + fh6cm6*radiation + fh6sm1*soilPC1 + fh6sm2*soilPC2
	diamDiversity ~ fh7cm1*wetness + fh7cm2*rain_lowRH + fh7cm3*iso + fh7cm4*pseas + fh7cm5*mat + fh7cm6*radiation + fh7sm1*soilPC1 + fh7sm2*soilPC2
	bark_moist_pct.ba ~ fm1cm1*wetness + fm1cm2*rain_lowRH + fm1cm3*iso + fm1cm4*pseas + fm1cm5*mat + fm1cm6*radiation + fm1sm1*soilPC1 + fm1sm2*soilPC2
	wood_SG.ba ~ fm2cm1*wetness + fm2cm2*rain_lowRH + fm2cm3*iso + fm2cm4*pseas + fm2cm5*mat + fm2cm6*radiation + fm2sm1*soilPC1 + fm2sm2*soilPC2
	LogSeed.ba ~ fm3cm1*wetness + fm3cm2*rain_lowRH + fm3cm3*iso + fm3cm4*pseas + fm3cm5*mat + fm3cm6*radiation + fm3sm1*soilPC1 + fm3sm2*soilPC2
	bigTrees ~ fm4cm1*wetness + fm4cm2*rain_lowRH + fm4cm3*iso + fm4cm4*pseas + fm4cm5*mat + fm4cm6*radiation + fm4sm1*soilPC1 + fm4sm2*soilPC2
	light.mean  ~ fm5cm1*wetness + fm5cm2*rain_lowRH + fm5cm3*iso + fm5cm4*pseas + fm5cm5*mat + fm5cm6*radiation + fm5sm1*soilPC1 + fm5sm2*soilPC2
	PC1 ~ fm6cm1*wetness + fm6cm2*rain_lowRH + fm6cm3*iso + fm6cm4*pseas + fm6cm5*mat + fm6cm6*radiation + fm6sm1*soilPC1 + fm6sm2*soilPC2

	# Regional climate, pollution, and regional forest heterogeneity effects on regional richness
	regS ~ R1CM1*wetness_reg_mean + R1CM2*rain_lowRH_reg_mean + R1CM3*iso_reg_mean + R1CM4*pseas_reg_mean + R1CM5*mat_reg_mean +
		R1CH1*wetness_reg_var + R1CH2*rain_lowRH_reg_var + R1CH3*iso_reg_var + R1CH4*pseas_reg_var + R1CH5*mat_reg_var +
		R1FH1*regS_tree 

	# Regional climate effects on regional forest heterogeneity
	regS_tree ~ FH1CM1*wetness_reg_mean + FH1CM2*rain_lowRH_reg_mean + FH1CM3*iso_reg_mean + FH1CM4*pseas_reg_mean + FH1CM5*mat_reg_mean +
		FH1CH1*wetness_reg_var + FH1CH2*rain_lowRH_reg_var + FH1CH3*iso_reg_var + FH1CH4*pseas_reg_var + FH1CH5*mat_reg_var

	# Effects on local lichen richness
	lichen.rich_log ~ r1fh1*bark_moist_pct.rao.ba + r1fh2*wood_SG.rao.ba + r1fh3*LogSeed.rao.ba +
		r1fh4*PIE.ba.tree + r1fh5*propDead + r1fh6*lightDist.mean + r1fh7*diamDiversity +
		r1fm1*bark_moist_pct.ba + r1fm2*wood_SG.ba + r1fm3*LogSeed.ba + r1fm4*bigTrees + r1fm5*light.mean + r1fm6*PC1 +
		r1cm1*wetness + r1cm2*rain_lowRH + r1cm3*iso + r1cm4*pseas + r1cm5*mat + r1cm6*radiation + r1sm1*soilPC1 + r1sm2*soilPC2 +
		r1R1*regS + r1a1*tot_abun_log + r1CM1*wetness_reg_mean + r1CM3*rain_lowRH_reg_mean + r1CM3*iso_reg_mean + r1CM4*pseas_reg_mean + r1CM5*mat_reg_mean +
		r1CH1*wetness_reg_var + r1CH2*rain_lowRH_reg_var + r1CH3*iso_reg_var + r1CH4*pseas_reg_var + r1CH5*mat_reg_var +
		r1FH1*regS_tree
	
	# Effects on local lichen abundance
	tot_abun_log ~ a1fh1*bark_moist_pct.rao.ba + a1fh2*wood_SG.rao.ba + a1fh3*LogSeed.rao.ba +
		a1fh4*PIE.ba.tree + a1fh5*propDead + a1fh6*lightDist.mean + a1fh7*diamDiversity +
		a1fm1*bark_moist_pct.ba + a1fm2*wood_SG.ba + a1fm3*LogSeed.ba + a1fm4*bigTrees + a1fm5*light.mean + a1fm6*PC1 +
		a1cm1*wetness + a1cm2*rain_lowRH + a1cm3*iso + a1cm4*pseas + a1cm5*mat + a1cm6*radiation + 
		a1sm1*soilPC1 + a1sm2*soilPC2 

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

soilmod_nopol_fit =  sem(path_soilmod_nopol, data=working_soil_test, fixed.x=T, estimator='ML', se='robust.sem')
use_fit = soilmod_nopol_fit

mod_boot = bootstrapLavaan(use_fit, R=10000, FUN=function(x) c(parameterEstimates(x)$est,standardizedSolution(x)$est.std))
soilmod_nopol_boot = mod_boot

# Used to re-calculate tables outside of Kure
#load('soilmod_nopol_AllSp_testdata_output.RData') 
#mod_boot = soilmod_nopol_boot
#use_fit = soilmod_nopol_fit 

# Save raw bootstrap output and models
response = 'AllSp'
modform = 'soilmod_regTorich_nopol'
save(soilmod_nopol_fit, soilmod_nopol_boot, path_soilmod_nopol, file=paste(modform,response,'testdata_output.RData', sep='_'))

## Calculate table of bootstrapped parameter estimates
ests = parameterEstimates(use_fit)[,c('label','lhs','op','rhs')]
nEst = ncol(mod_boot)/2 # number of parameters
estlabs = paste(predtypes[ests$lhs,'label'],predtypes[ests$rhs,'label'], sep=':')
cbind(estlabs, ests[,c('lhs','op','rhs')]) #checking
stdmat = mod_boot[,(nEst+1):(nEst*2)]; colnames(stdmat) = estlabs
ests$label = estlabs
ests$std.all = apply(mod_boot[,(nEst+1):(2*nEst)], 2, mean)
ests$std.se = apply(mod_boot[,(nEst+1):(2*nEst)], 2, function(x) sqrt(var(x)))
ests$std.ci.lower = apply(mod_boot[,(nEst+1):(2*nEst)], 2, function(x) quantile(x, p=0.025))
ests$std.ci.upper = apply(mod_boot[,(nEst+1):(2*nEst)], 2, function(x) quantile(x, p=0.975))

# Calculate indirect effects via abundance for all variables
use_ests = ests$lhs=='tot_abun_log'&ests$op=='~'
IE_abun_mat = stdmat[,use_ests]*stdmat[,'r1:a1'] # multiply paths: var -> abun -> rich
colnames(IE_abun_mat) = ests[use_ests,c('rhs')]
indirect_abun = data.frame(predictor=colnames(IE_abun_mat))
indirect_abun$std.all = apply(IE_abun_mat, 2, mean)
indirect_abun$std.se = apply(IE_abun_mat, 2, function(x) sqrt(var(x)))
indirect_abun$std.ci.lower = apply(IE_abun_mat, 2, function(x) quantile(x, p=0.025))
indirect_abun$std.ci.upper = apply(IE_abun_mat, 2, function(x) quantile(x, p=0.975))

# Calculate direct effects on local richness
use_ests = ests$lhs=='lichen.rich_log'&ests$op=='~'
DE_mat = stdmat[,use_ests]
colnames(DE_mat) = ests[use_ests,c('rhs')]
direct_rich = data.frame(predictor=colnames(DE_mat))
direct_rich$std.all = apply(DE_mat, 2, mean)
direct_rich$std.se = apply(DE_mat, 2, function(x) sqrt(var(x)))
direct_rich$std.ci.lower = apply(DE_mat, 2, function(x) quantile(x, p=0.025))
direct_rich$std.ci.upper = apply(DE_mat, 2, function(x) quantile(x, p=0.975))

# Calculate direct effects on abundance
use_ests = ests$lhs=='tot_abun_log'&ests$op=='~'
DE_abun_mat = stdmat[,use_ests]
colnames(DE_abun_mat) = ests[use_ests,'rhs']
direct_abun = data.frame(predictor=colnames(DE_abun_mat))
direct_abun$std.all = apply(DE_abun_mat, 2, mean)
direct_abun$std.se = apply(DE_abun_mat, 2, function(x) sqrt(var(x)))
direct_abun$std.ci.lower = apply(DE_abun_mat, 2, function(x) quantile(x, p=0.025))
direct_abun$std.ci.upper = apply(DE_abun_mat, 2, function(x) quantile(x, p=0.975))

# Calculate direct effects on regional richness
use_ests = ests$lhs=='regS'&ests$op=='~'
DE_reg_mat = stdmat[,use_ests]
colnames(DE_reg_mat) = ests[use_ests,c('rhs')]
direct_reg = data.frame(predictor=colnames(DE_reg_mat))
direct_reg$std.all = apply(DE_reg_mat, 2, mean)
direct_reg$std.se = apply(DE_reg_mat, 2, function(x) sqrt(var(x)))
direct_reg$std.ci.lower = apply(DE_reg_mat, 2, function(x) quantile(x, p=0.025))
direct_reg$std.ci.upper = apply(DE_reg_mat, 2, function(x) quantile(x, p=0.975))

# Calculate total effects for local forest structure variables
fvars = rownames(predtypes[grep('f',predtypes$label),])
TE_f_mat = IE_abun_mat[,fvars] + DE_mat[,fvars]
colnames(TE_f_mat) = fvars

# Calculate indirect effects via forest structure for local climate
cvars = rownames(predtypes[grep('c',predtypes$label),])

# Make an array that holds bootstrap replicates of paths c -> f
cf_mat = array(NA, dim=c(nrow(stdmat),length(fvars),length(cvars)), dimnames=list(1:nrow(stdmat),fvars, cvars))
combos = expand.grid(fvars, cvars, stringsAsFactors=F)
for(i in 1:nrow(combos)){
	yvar = combos[i,1]
	xvar = combos[i,2]
	this_lab = paste(predtypes[yvar,'label'],predtypes[xvar,'label'], sep=':')
	cf_mat[,yvar, xvar] = stdmat[,this_lab]
}
# multiply c -> f paths by total effects of each f variable
cm_IE_for_mat = array(apply(cf_mat, 3, function(x) x*TE_f_mat),dim=c(nrow(stdmat),length(fvars),length(cvars)), dimnames=list(1:nrow(stdmat),fvars, cvars)) 

# sum effects across fh and fm categories for each c variable
cm_IE_fh_mat = apply(cm_IE_for_mat[,fvars[predtypes[fvars,'mode']=='het'],], c(1,3), function(x) sum(x))
cm_IE_fm_mat = apply(cm_IE_for_mat[,fvars[predtypes[fvars,'mode']=='opt'],], c(1,3), function(x) sum(x))

indirect_for = data.frame(predictor=rep(cvars, 2), Ftype = rep(c('het','opt'), each=length(cvars)))
indirect_for$std.all = apply(cbind(cm_IE_fh_mat, cm_IE_fm_mat), 2, mean)
indirect_for$std.se = apply(cbind(cm_IE_fh_mat, cm_IE_fm_mat), 2, function(x) sqrt(var(x)))
indirect_for$std.ci.lower = apply(cbind(cm_IE_fh_mat, cm_IE_fm_mat), 2, function(x) quantile(x, p=0.025))
indirect_for$std.ci.upper = apply(cbind(cm_IE_fh_mat, cm_IE_fm_mat), 2, function(x) quantile(x, p=0.975))

# Calculate total effects of local climate
TE_cm_mat = DE_mat[,cvars]+IE_abun_mat[,cvars]+cm_IE_fm_mat[,cvars]+cm_IE_fh_mat[,cvars]

# Calculate indirect effects via regional richness (for all regional scale variables)
use_ests = ests$lhs=='regS'&ests$op=='~'
IE_regS_mat = stdmat[,use_ests]*DE_mat[,'regS'] # var -> regS * regS -> lichen_rich
colnames(IE_regS_mat) = ests[use_ests, 'rhs']

# Calculate indirect path of regional forest via local forest
# ests[grep(':FH1', ests$label),] # only correlated with fh4 (PIE.ba.tree)
FH_IE_for_mat = stdmat[,'fh4:FH1']*TE_f_mat[,'PIE.ba.tree']
FH_IE_for_mat = matrix(FH_IE_for_mat, nrow=nrow(stdmat)); colnames(FH_IE_for_mat) = 'regS_tree'

# Calculate total effect/path of regional forest (effect does not include indirect path via local forest)
TP_FH_mat = DE_mat[,'regS_tree'] + IE_regS_mat[,'regS_tree'] + FH_IE_for_mat
TP_FH_mat = matrix(TP_FH_mat, nrow=nrow(stdmat)); colnames(TP_FH_mat) = 'regS_tree'
TE_FH_mat = DE_mat[,'regS_tree'] + IE_regS_mat[,'regS_tree']
TE_FH_mat = matrix(TE_FH_mat, nrow=nrow(stdmat)); colnames(TE_FH_mat) = 'regS_tree'

## Calculate indirect effects of regional climate via regional versus local paths 
# e.g. do or do not go through local predictors
Cvars = rownames(predtypes)[grep('C', predtypes$label)]

# Forest (via forest's regional, local, and direct effects)
use_labs = paste('FH1',predtypes[Cvars, 'label'], sep=':')
C_IE_FH_mat_reg = stdmat[,use_labs]*IE_regS_mat[,'regS_tree']
colnames(C_IE_FH_mat_reg) = Cvars
C_IE_FH_mat_loc = stdmat[,use_labs]*as.numeric(FH_IE_for_mat)
colnames(C_IE_FH_mat_loc) = Cvars
C_IE_FH_mat_dir = stdmat[,use_labs]*DE_mat[,'regS_tree']
colnames(C_IE_FH_mat_dir) = Cvars
C_IE_FH_mat = C_IE_FH_mat_reg + C_IE_FH_mat_loc + C_IE_FH_mat_dir

# RegS
#IE_regS_mat[,Cvars]

# Local climate: Make an array that holds bootstrap replicates of paths CH/M -> cm
Cc_mat = array(NA, dim=c(nrow(stdmat),length(cvars),length(Cvars)), dimnames=list(1:nrow(stdmat), cvars, Cvars))
combos = expand.grid(cvars, Cvars, stringsAsFactors=F)
for(i in 1:nrow(combos)){
	yvar = combos[i,1]
	xvar = combos[i,2]
	this_lab = paste(predtypes[yvar,'label'],predtypes[xvar,'label'], sep=':')
	Cc_mat[,yvar, xvar] = stdmat[,this_lab]
}
# multiply C -> c paths by total effects of each c variable
C_IE_cm_mat_ind = array(apply(Cc_mat, 3, function(x) x*TE_cm_mat),dim=c(nrow(stdmat),length(cvars),length(Cvars)), dimnames=list(1:nrow(stdmat),cvars, Cvars)) 

# sum across local c vars
C_IE_cm_mat = apply(C_IE_cm_mat_ind, c(1,3), function(x) sum(x))

# Make datatable of indirect effects of regional variables
IE_mat = cbind(IE_regS_mat, FH_IE_for_mat, C_IE_FH_mat, C_IE_cm_mat, C_IE_FH_mat_reg, C_IE_FH_mat_loc, C_IE_FH_mat_dir)
indirect_reg = data.frame(predictor=colnames(IE_mat))
indirect_reg$IEvar1 = c(rep('regS',ncol(IE_regS_mat)), rep('PIE.ba.tree', ncol(FH_IE_for_mat)), 
	rep('regS_tree',ncol(C_IE_FH_mat)),
	rep('clim_loc', ncol(C_IE_cm_mat)),
	rep('regS_tree',3*ncol(C_IE_FH_mat)))
indirect_reg$IEvar2 = c(rep(NA, ncol(IE_mat) - 3*ncol(C_IE_FH_mat)),
	rep(c('regS','loc','dir'),each=ncol(C_IE_FH_mat)))
indirect_reg$std.all = apply(IE_mat, 2, mean)
indirect_reg$std.se = apply(IE_mat, 2, function(x) sqrt(var(x)))
indirect_reg$std.ci.lower = apply(IE_mat, 2, function(x) quantile(x, p=0.025))
indirect_reg$std.ci.upper = apply(IE_mat, 2, function(x) quantile(x, p=0.975))


# Calculate total paths of regional climate
TP_C_mat = DE_mat[,Cvars] + IE_regS_mat[,Cvars] + C_IE_FH_mat[,Cvars] + C_IE_cm_mat[,Cvars]

# Calculate total effects of regional climate
TE_C_mat = DE_mat[,Cvars] + IE_regS_mat[,Cvars] + C_IE_FH_mat_reg[,Cvars]

# Make dataframe of total effects
TE_mat = cbind(TE_C_mat, TE_cm_mat, TE_f_mat, TE_FH_mat)
total = data.frame(predictor=colnames(TE_mat))
total$std.all = apply(TE_mat, 2, mean)
total$std.se = apply(TE_mat, 2, function(x) sqrt(var(x)))
total$std.ci.lower = apply(TE_mat, 2, function(x) quantile(x, p=0.025))
total$std.ci.upper = apply(TE_mat, 2, function(x) quantile(x, p=0.975))

# Make dataframe of total paths (include paths through local-regional correlations)
TP_mat = cbind(TP_C_mat, TE_cm_mat, TE_f_mat, TP_FH_mat)
totalp = data.frame(predictor=colnames(TP_mat))
totalp$std.all = apply(TP_mat, 2, mean)
totalp$std.se = apply(TP_mat, 2, function(x) sqrt(var(x)))
totalp$std.ci.lower = apply(TP_mat, 2, function(x) quantile(x, p=0.025))
totalp$std.ci.upper = apply(TP_mat, 2, function(x) quantile(x, p=0.975))

## Save tables and matrices
write.csv(ests, paste(modform,response,'testdata_parameterEstimates.csv', sep='_'), row.names=F)
write.csv(direct_rich, paste(modform,response,'testdata_directeffects_richness.csv', sep='_'), row.names=F)
write.csv(direct_abun, paste(modform,response,'testdata_directeffects_abundance.csv', sep='_'), row.names=F)
write.csv(direct_reg, paste(modform,response,'testdata_directeffects_regS.csv', sep='_'), row.names=F)
write.csv(indirect_abun, paste(modform,response,'testdata_indirecteffects_via_abundance.csv', sep='_'), row.names=F)
write.csv(indirect_for, paste(modform,response,'testdata_indirecteffects_via_forest.csv', sep='_'), row.names=F)
write.csv(indirect_reg, paste(modform,response,'testdata_regionalvars_indirecteffects.csv', sep='_'), row.names=F)
write.csv(total, paste(modform,response, 'testdata_totaleffects.csv', sep='_'), row.names=F)
write.csv(totalp, paste(modform,response, 'testdata_totalpaths.csv', sep='_'), row.names=F)
