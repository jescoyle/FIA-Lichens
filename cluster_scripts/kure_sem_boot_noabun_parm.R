# This script calculates bootstrapped ci for standardized coefficients in my final path model of local lichen richness

library(lavaan, lib.loc="/nas02/home/j/r/jrcoyle/Rlibs/")

options(stringsAsFactors=F)

working_data_test = read.csv('standardized_test_dataset.csv', row.names=1)

# Model
## Define model that will calculate indirect and total effects
noabun_mod = '

	# Latent variables
	lichen_rich =~ sqrt(0.75)*Parm_log

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
	#totalCirc ~ FM4a*wetness + FM4b*rain_lowRH + FM4c*iso + FM4d*pseas + FM4e*mat + FM4f*radiation
	bigTrees ~ FM5a*wetness + FM5b*rain_lowRH + FM5c*iso + FM5d*pseas + FM5e*mat + FM5f*radiation
	light.mean  ~ FM6a*wetness + FM6b*rain_lowRH + FM6c*iso + FM6d*pseas + FM6e*mat + FM6f*radiation
	PC1 ~ FM7a*wetness + FM7b*rain_lowRH + FM7c*iso + FM7d*pseas + FM7e*mat + FM7f*radiation

	# Climate effects on regional richness
	regParm ~ R1*wetness + R2*rain_lowRH + R3*iso + R4*pseas + R5*mat

	# Climate effects on pollution
	totalNS ~ P1*wetness + P2*rain_lowRH

	# Lichen richness regression
	lichen_rich ~ FH1*bark_moist_pct.rao.ba + FH2*wood_SG.rao.ba + FH3*LogSeed.rao.ba +
		FH4*PIE.ba.tree + FH5*propDead + FH6*lightDist.mean + FH7*diamDiversity +
		FM1*bark_moist_pct.ba + FM2*wood_SG.ba + FM3*LogSeed.ba + FM5*bigTrees + FM6*light.mean + FM7*PC1 + 
		C1*wetness + C2*rain_lowRH + C3*iso + C4*pseas + C5*mat + C6*radiation + R*regParm + P*totalNS

	# Indirect effects of climate and local environment on lichen richness
	IE_wetness_FH := FH1a*FH1 + FH2a*FH2 + FH3a*FH3 + FH4a*FH4 + FH5a*FH5 + 
		FH6a*FH6 + FH7a*FH7
	IE_wetness_FM := FM1a*FM1 + FM2a*FM2 + FM3a*FM3 + FM5a*FM5 + 
		FM6a*FM6 + FM7a*FM7
	IE_wetness_R := R1*R
	IE_wetness_P := P1*P
	IE_wetness := IE_wetness_FH + IE_wetness_FM + IE_wetness_R + IE_wetness_P

	IE_rain_lowRH_FH := FH1b*FH1 + FH2b*FH2 +	FH3b*FH3 + FH4b*FH4 + 
		FH5b*FH5 + FH6b*FH6 + FH7b*FH7 
	IE_rain_lowRH_FM := FM1b*FM1 + FM2b*FM2 + FM3b*FM3 + FM5b*FM5 + 
		FM6b*FM6 + FM7b*FM7
	IE_rain_lowRH_R := R2*R
	IE_rain_lowRH_P := P2*P
	IE_rain_lowRH := IE_rain_lowRH_FH + IE_rain_lowRH_FM + IE_rain_lowRH_R + IE_rain_lowRH_P

	IE_iso_FH := FH1c*FH1 + FH2c*FH2 + FH3c*FH3 + FH4c*FH4 + FH5c*FH5 +
		FH6c*FH6 + FH7c*FH7 
	IE_iso_FM := FM1c*FM1 + FM2c*FM2 + FM3c*FM3 + FM5c*FM5 + 
		FM6c*FM6 +  FM7c*FM7
	IE_iso_R := R3*R
	IE_iso := IE_iso_FH + IE_iso_FM + IE_iso_R

	IE_pseas_FH := FH1d*FH1 + FH2d*FH2 + FH3d*FH3 + FH4d*FH4 + FH5d*FH5 + 
		FH6d*FH6 + FH7d*FH7
	IE_pseas_FM := FM1d*FM1 + FM2d*FM2 + FM3d*FM3 + FM5d*FM5 + 
		FM6d*FM6 + FM7d*FM7
	IE_pseas_R := R4*R
	IE_pseas := IE_pseas_FH + IE_pseas_FM + IE_pseas_R
	
	IE_mat_FH := FH1e*FH1 + FH2e*FH2 + FH3e*FH3 + FH4e*FH4 + FH5e*FH5 + 
		FH6e*FH6 + FH7e*FH7
	IE_mat_FM := FM1e*FM1 + FM2e*FM2 + FM3e*FM3 + FM5e*FM5 + 
		FM6e*FM6 + FM7e*FM7
	IE_mat_R := R5*R
	IE_mat := IE_mat_FH + IE_mat_FM + IE_mat_R

	IE_radiation_FH := FH1f*FH1 + FH2f*FH2 + FH3f*FH3 + FH4f*FH4 + FH5f*FH5 + 
		FH6f*FH6 + FH7f*FH7
	IE_radiation_FM := FM1f*FM1 + FM2f*FM2 + FM3f*FM3 + FM5f*FM5 + 
		FM6f*FM6 + FM7f*FM7
	IE_radiation := IE_radiation_FH + IE_radiation_FM

	# Total effects on lichen richness
	TE_wetness := IE_wetness + C1	
	TE_rain_lowRH := IE_rain_lowRH + C2
	TE_iso := IE_iso + C3
	TE_pseas := IE_pseas + C4
	TE_mat := IE_mat + C5
	TE_radiation := IE_radiation + C6

	TE_bark_moist_pct.rao.ba := FH1
	TE_wood_SG.rao.ba := FH2
	TE_LogSeed.rao.ba := FH3
	TE_PIE.ba.tree := FH4
	TE_propDead := FH5
	TE_lightDist.mean := FH6
	TE_diamDiversity := FH7
	TE_bark_moist_pct.ba := FM1
	TE_wood_SG.ba := FM2
	TE_LogSeed.ba := FM3
	#TE_totalCirc := FM4
	TE_bigTrees := FM5
	TE_light.mean := FM6
	TE_PC1 := FM7	

	TE_FH := TE_bark_moist_pct.rao.ba +	TE_wood_SG.rao.ba + TE_LogSeed.rao.ba +
		TE_PIE.ba.tree + TE_propDead + TE_lightDist.mean + TE_diamDiversity
	TE_FM := TE_bark_moist_pct.ba + TE_wood_SG.ba + TE_LogSeed.ba +
		TE_bigTrees + TE_light.mean + TE_PC1

	TE_totalNS := P
'

# Fit model
noabun_fit =  sem(noabun_mod, data=working_data_test, fixed.x=T, estimator='ML', se='robust')


# Bootstrap
noabun_std = bootstrapLavaan(noabun_fit, R=10000, FUN=function(x) c(parameterEstimates(x)$est,standardizedSolution(x)$est.std))

# Save bootstrap output
write.csv(noabun_std, 'boot_parm_noabun.csv', row.names=F)

# Make table of parameter estimates
noabun_ests = parameterEstimates(noabun_fit)[,c('label','lhs','op','rhs')]
nEst = ncol(noabun_std)/2 # number of parameters
noabun_ests$std.all = apply(noabun_std[,(nEst+1):(2*nEst)], 2, mean)
noabun_ests$std.se = apply(noabun_std[,(nEst+1):(2*nEst)], 2, function(x) sqrt(var(x)))
noabun_ests$std.ci.lower = apply(noabun_std[,(nEst+1):(2*nEst)], 2, function(x) quantile(x, p=0.025))
noabun_ests$std.ci.upper = apply(noabun_std[,(nEst+1):(2*nEst)], 2, function(x) quantile(x, p=0.975))

# Make a table of direct effects on lichen richness
direct_rich = subset(noabun_ests, (lhs=='lichen_rich')&(op=='~'))
direct_rich = direct_rich[,c('rhs','std.all','std.se','std.ci.lower','std.ci.upper')] 
names(direct_rich)[1] = 'predictor'
rownames(direct_rich) = direct_rich$predictor
rownames(direct_rich)[rownames(direct_rich)=='regParm'] = 'reg' 

# Make a table of total effects on lichen richness
total = noabun_ests[grep('TE',noabun_ests$label),]
total = total[,c('lhs','std.all','std.se','std.ci.lower','std.ci.upper')]
names(total)[1] = 'predictor'
total$predictor = substring(total$predictor, first=4)
addrows = subset(direct_rich, !(predictor %in% total$predictor))
colnames(addrows) = colnames(total)
total = rbind(total, addrows)
rownames(total) = total$predictor
rownames(total)[rownames(total)=='regParm'] = 'reg' 

# Make a table of indirect effect of climate variable via regional richness
indirect_reg = noabun_ests[grep('IE_[A-za-z._]*_R', noabun_ests$label),]
indirect_reg = indirect_reg[,c('lhs','std.all','std.se','std.ci.lower','std.ci.upper')]
names(indirect_reg)[1] = 'predictor'
indirect_reg$predictor = substring(indirect_reg$predictor, first=4)
indirect_reg$predictor = sapply(indirect_reg$predictor, function (x) substr(x, 1, nchar(x)-2))
rownames(indirect_reg) = indirect_reg$predictor

# Make a table of indirect effect of climate and local environment via forest structure
indirect_for = noabun_ests[grep('IE_[A-za-z._1-9]*_F[MH]', noabun_ests$label),]
indirect_for$Ftype = sapply(indirect_for$label, function(x) substr(x, nchar(x)-1, nchar(x)))
indirect_for = indirect_for[,c('lhs','Ftype','std.all','std.se','std.ci.lower','std.ci.upper')]
names(indirect_for)[1] = 'predictor'
indirect_for$predictor = substring(indirect_for$predictor, first=4)
rownames(indirect_for) = indirect_for$predictor
indirect_for$predictor = sapply(indirect_for$predictor, function (x) substr(x, 1, nchar(x)-3))

# Append variable types column
vartypes = read.csv('var_types.csv', row.names=1)
total$type = vartypes[rownames(total),'type']
direct_rich$type = vartypes[rownames(direct_rich),'type']
indirect_for$type = vartypes[indirect_for$predictor, 'type']
indirect_reg$type = vartypes[indirect_reg$predictor, 'type']

# Save tables
write.csv(noabun_ests, 'Parm_testdata_noabun_parameterEstimates.csv')
write.csv(total, 'Parm_testdata_noabun_totaleffects.csv', row.names=T)
write.csv(direct_rich, 'Parm_testdata_noabun_directeffects_richness.csv', row.names=T)
write.csv(indirect_reg, 'Parm_testdata_noabun_indirecteffects_via_regS.csv', row.names=T)
write.csv(indirect_for, 'Parm_testdata_noabun_indirecteffects_via_forest.csv', row.names=T)

save.image('bootstrap_parm_noabun_model.RData')
