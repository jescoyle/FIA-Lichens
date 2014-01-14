# This script calculates bootstrapped ci for standardized coefficients in my final path model of local lichen richness

library(lavaan, lib.loc="/nas02/home/j/r/jrcoyle/Rlibs/")

options(stringsAsFactors=F)

working_data_test = read.csv('standardized_test_dataset.csv', row.names=1)

working_data_test = subset(working_data_test, !is.na(raoQ))

# Model
path_measerr = '

	# Latent variables- this is retained so that FD model is comparable with richness model
	lichen_rich =~ sqrt(0.75)*raoQ

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

# Fit model
endfit =  sem(path_measerr, data=working_data_test, fixed.x=T, estimator='ML', se='robust')


# Bootstrap
endfit_std = bootstrapLavaan(endfit, R=10000, FUN=function(x) c(parameterEstimates(x)$est,standardizedSolution(x)$est.std))

# Save bootstrap output
write.csv(endfit_std, 'boot_raoQ.csv', row.names=F)

# Make table of parameter estimates
endfit_ests = parameterEstimates(endfit)[,c('label','lhs','op','rhs')]
nEst = ncol(endfit_std)/2 # number of parameters
endfit_ests$std.all = apply(endfit_std[,(nEst+1):(2*nEst)], 2, mean)
endfit_ests$std.se = apply(endfit_std[,(nEst+1):(2*nEst)], 2, function(x) sqrt(var(x)))
endfit_ests$std.ci.lower = apply(endfit_std[,(nEst+1):(2*nEst)], 2, function(x) quantile(x, p=0.025))
endfit_ests$std.ci.upper = apply(endfit_std[,(nEst+1):(2*nEst)], 2, function(x) quantile(x, p=0.975))

# Make a table of direct effects on lichen richness
direct_rich = subset(endfit_ests, (lhs=='lichen_rich')&(op=='~'))
direct_rich = direct_rich[,c('rhs','std.all','std.se','std.ci.lower','std.ci.upper')] 
names(direct_rich)[1] = 'predictor'
rownames(direct_rich) = direct_rich$predictor
rownames(direct_rich)[rownames(direct_rich)=='regS'] = 'reg' 
rownames(direct_rich)[rownames(direct_rich)=='tot_abun_log'] = 'abun_log' 

# Make a table of total effects on lichen richness
total = endfit_ests[grep('TE',endfit_ests$label),]
total = total[,c('lhs','std.all','std.se','std.ci.lower','std.ci.upper')]
names(total)[1] = 'predictor'
total$predictor = substring(total$predictor, first=4)
addrows = subset(direct_rich, !(predictor %in% total$predictor))
colnames(addrows) = colnames(total)
total = rbind(total, addrows)
rownames(total) = total$predictor
rownames(total)[rownames(total)=='regS'] = 'reg' 
rownames(total)[rownames(total)=='tot_abun_log'] = 'abun_log' 


# Make a table of indirect effects on lichen richness via abundance
# For climate and local environment variables that also affect forest structure, this is just the path directly from climate to abundance, not via forest structure
c(paste('IE',c('bark_moist_pct.rao.ba','wood_SG.rao.ba','LogSeed.rao.ba','PIE.ba.tree','propDead','lightDist.mean','diamDiversity',
		'bark_moist_pct.ba','wood_SG.ba','LogSeed.ba','totalCirc','bigTrees','light.mean','PC1','totalNS'), sep='_'),
	paste('IE',c('wetness','rain_lowRH','iso','pseas','mat','radiation'),'A', sep='_'))->indir_vars
indirect = subset(endfit_ests, label %in% indir_vars)
indirect = indirect[,c('lhs','std.all','std.se','std.ci.lower','std.ci.upper')]
names(indirect)[1] = 'predictor'
indirect$predictor = substring(indirect$predictor, first=4)
indirect$predictor = sapply( indirect$predictor, function(x) ifelse(substr(x, nchar(x), nchar(x))=='A', substr(x, 1, nchar(x)-2), x))	
rownames(indirect) = indirect$predictor

# Make a table of direct effects on lichen abundance
direct_abun = subset(endfit_ests, (lhs=='tot_abun_log')&(op=='~'))
direct_abun = direct_abun[,c('rhs','std.all','std.se','std.ci.lower','std.ci.upper')]
names(direct_abun)[1] = 'predictor'
rownames(direct_abun) = direct_abun$predictor

# Make a table of indirect effect of climate variable via regional richness
indirect_reg = endfit_ests[grep('IE_[A-za-z._]*_R', endfit_ests$label),]
indirect_reg = indirect_reg[,c('lhs','std.all','std.se','std.ci.lower','std.ci.upper')]
names(indirect_reg)[1] = 'predictor'
indirect_reg$predictor = substring(indirect_reg$predictor, first=4)
indirect_reg$predictor = sapply(indirect_reg$predictor, function (x) substr(x, 1, nchar(x)-2))
rownames(indirect_reg) = indirect_reg$predictor

# Append variable types column
vartypes = read.csv('var_types.csv', row.names=1)
total$type = vartypes[rownames(total),'type']
direct_rich$type = vartypes[rownames(direct_rich),'type']
direct_abun$type = vartypes[rownames(direct_abun),'type']
indirect$type = vartypes[rownames(indirect),'type']


# Save tables
write.csv(endfit_ests, 'RaoQ_testdata_parameterEstimates.csv')
write.csv(total, 'RaoQ_testdata_totaleffects.csv', row.names=T)
write.csv(direct_rich, 'RaoQ_testdata_directeffects_richness.csv', row.names=T)
write.csv(direct_abun, 'RaoQ_testdata_directeffects_abundance.csv', row.names=T)
write.csv(indirect, 'RaoQ_testdata_indirecteffects_via_abundance.csv', row.names=T)
write.csv(indirect_reg, 'RaoQ_testdata_indirecteffects_via_regS.csv', row.names=T)


save.image('bootstrap_raoQ_model.RData')
