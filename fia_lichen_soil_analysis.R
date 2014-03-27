## This script explores soil data asscoiated with the FIA lichen plots

source('./UNC/Projects/FIA Lichen/GitHub/FIA-Lichens/load_data.R')
source('./GitHub/FIA-Lichens/fia_lichen_analysis_functions.R')

## Read in soil data

statelist = unique(master$state.abbr)

filename = paste('./Data/RawData/FIA_SOIL/',statelist[1],'_SOILS_LAB.CSV', sep='')
soil_lab = read.csv(filename)

# Combine soil lab data into one dataframe
for(s in statelist[2:length(statelist)]){
	filename = paste('./Data/RawData/FIA_SOIL/',s,'_SOILS_LAB.CSV', sep='')
	new_data = read.csv(filename)
	soil_lab = rbind(soil_lab, new_data)
}


## How many plots have matching soil data? And what do they look like?

soil_lab$yrplot.id = make.yrplotid(soil_lab)
sum(rownames(model_data) %in% unique(soil_lab$yrplot.id)) # 703 out of 1923 plots will have soil data

## Where is soil data measured on these plots

# Restrict analysis to plots in focal plots
use_soil = subset(soil_lab, yrplot.id %in% rownames(model_data))

# Restrict analysis to non-litter
use_soil = subset(use_soil, LAYER_TYPE %in% c('MIN_1','MIN_2','ORG_1','ORG_2'))

# Number of plots and locations
length(unique(use_soil$yrplot.id))
soil_plots = subset(master, yrplot.id %in% unique(use_soil$yrplot.id))
table(soil_plots$state.abbr)

# Where were soil samples collected in each plot?
counts = table(use_soil$yrplot.id, use_soil$SMPLNNBR)
hist(rowSums(counts))
sum(counts[,1]==0) # All but 6 plots had at least one sample collected in subplot 2


## Define soil data set to be analyzed

# Only use samples from subplot 2, unless missing, then use other subplot
use_subp = as.numeric(colnames(counts)[apply(counts, 1, function(x) min(grep('1|2', x)))])
names(use_subp) = rownames(counts)

soil_data = data.frame()
for(this_plot in rownames(counts)){
	add_data = subset(use_soil, (yrplot.id==this_plot)&(SMPLNNBR==use_subp[this_plot]))
	soil_data = rbind(soil_data, add_data)
}

# Only use samples in top 0-4 inches
soil_data = soil_data[grep('1', soil_data$LAYER_TYPE),]

# Check to make sure only one observation (line) per plot
table(soil_data$yrplot.id, soil_data$SMPLNNBR)
sum(table(soil_data$yrplot.id)!=1) # 0, all plots sampled once.

write.csv(soil_data, './Data/fia_soil_data.csv', row.names=F)

## Read in parsed soil data
soil_data = read.csv('./Data/fia_soil_data.csv')

# Define soil variables to examine
soil_vars = c('PH_H2O','EXCHNG_K','EXCHNG_CA','EXCHNG_MN','EXCHNG_MG')
soil = soil_data[,c('yrplot.id',soil_vars)]

# Distributions of soil variables
par(mfrow=c(2,3))
for(i in soil_vars){
	hist(soil[,i], main=i)
}

# Log transform ions
for(x in c('EXCHNG_K','EXCHNG_CA','EXCHNG_MN','EXCHNG_MG')){
	soil[,x] = log10(soil[,x]+1)
}

# Scale variables
soil[,soil_vars] = scale(soil[,soil_vars], center=T, scale=T)

## Exploratory analysis of soil factors to determine which to use in analysis
# Use 'fit' data not 'test' data
working_data$yrplot.id = rownames(working_data)

# Format data table
working_data$yrplot.id = rownames(working_data)
working_soil = merge(soil, working_data, all.x=T)
rownames(working_soil) = working_soil$yrplot.id

# Subset to training data set
fitplots = fitplots$yrplot.id
working_soil_fit = working_soil[fitplots[fitplots %in% rownames(working_soil)],]

# Univariate correlations and plots
use_vars = c('lichen.rich_log','PC1','PIE.ba.tree', soil_vars)

cor(working_soil[,use_vars], use='complete.obs')


# Path analysis
library(lavaan)
library(semPlot)

soilmod = '

	# Effect of soil and trees on lichens
	lichen.rich_log ~ T1*PC1 + T2*PIE.ba.tree + S1*PH_H2O + S2*EXCHNG_K + S3*EXCHNG_CA + S4*EXCHNG_MN + S5*EXCHNG_MG

	# Effect of soil on trees
	PC1 ~ T1S1*PH_H2O + T1S2*EXCHNG_K + T1S3*EXCHNG_CA + T1S4*EXCHNG_MN + T1S5*EXCHNG_MG
	PIE.ba.tree ~ T2S1*PH_H2O + T2S2*EXCHNG_K + T2S3*EXCHNG_CA + T2S4*EXCHNG_MN + T2S5*EXCHNG_MG
	
	# Indirect effect of soil
	IE_PH_PC1 := T1*T1S1
	IE_K_PC1 := T1*T1S2
	IE_CA_PC1 := T1*T1S3
	IE_MN_PC1 := T1*T1S4
	IE_MG_PC1 := T1*T1S5
	IE_PH_PIE := T2*T2S1
	IE_K_PIE := T2*T2S2
	IE_CA_PIE := T2*T2S3
	IE_MN_PIE := T2*T2S4
	IE_MG_PIE := T2*T2S5

	IE_PH := IE_PH_PC1 + IE_PH_PIE
	IE_K := IE_K_PC1 + IE_K_PIE
	IE_CA := IE_CA_PC1 + IE_CA_PIE
	IE_MN := IE_MN_PC1 + IE_MN_PIE
	IE_MG := IE_MG_PC1 + IE_MG_PIE

'
soilfit = sem(soilmod, data=working_soil_fit)

summary(soilfit)

semPaths(soilfit, what='est', whatLabel='est', layout='circle2', nCharNodes=0)


## PCA of soil factors
working_soil_noNA = na.omit(working_soil)

soilPCA = prcomp(working_soil_noNA[,soil_vars])
summary(soilPCA)

soilcatPCA = prcomp(working_soil_noNA[,c('EXCHNG_K','EXCHNG_CA','EXCHNG_MN','EXCHNG_MG')])
soilcatPCA
summary(soilcatPCA)

soilpHPCA = prcomp(working_soil_noNA[,c('PH_H2O','EXCHNG_K','EXCHNG_CA','EXCHNG_MG')])
soilpHPCA
summary(soilpHPCA)

# Decided to use PC1 and PC2 from PCA of all soil factors because:
# PC1: high pH + higher K, CA, MG
# PC2: low pH + higher Mn

# Add PC axes to soil data
working_soil$soilPC1 = NA
working_soil$soilPC2 = NA
working_soil[,c('soilPC1','soilPC2')] = predict(soilPCA, working_soil)[,c('PC1','PC2')]

working_soil$pHPC1 = predict(soilpHPCA, working_soil)[,'PC1']
plot(pHPC1~EXCHNG_MN, data=working_soil)

working_soil_fit = working_soil[fitplots[fitplots %in% rownames(working_soil)],]

soilmodPC = '

	# Effect of soil and trees on lichens
	lichen.rich_log ~ T1*PC1 + T2*PIE.ba.tree + S1*soilPC1 + S2*soilPC2

	# Effect of soil on trees
	PC1 ~ T1S1*soilPC1 + T1S2*soilPC2
	PIE.ba.tree ~ T2S1*soilPC1 + T2S2*soilPC2
	
	# Indirect effect of soil
	IE_soil1_PC1 := T1*T1S1
	IE_soil2_PC1 := T1*T1S2
	IE_soil1_PIE := T2*T2S1
	IE_soil2_PIE := T2*T2S2

	IE_soil1 := IE_soil1_PC1 + IE_soil1_PIE
	IE_soil2 := IE_soil2_PC1 + IE_soil2_PIE
'
soilfitPC = sem(soilmodPC, data=working_soil_fit)
summary(soilfitPC, rsq=TRUE, fit.measures=T)
semPaths(soilfitPC, what='std', whatLabel='std', layout='circle2', nCharNodes=0)

soilmodpH1 = '

	# Effect of soil and trees on lichens
	lichen.rich_log ~ T1*PC1 + T2*PIE.ba.tree + S1*pHPC1 + S2*EXCHNG_MN

	# Effect of soil on trees
	PC1 ~ T1S1*pHPC1 + T1S2*EXCHNG_MN
	PIE.ba.tree ~ T2S1*pHPC1 + T2S2*EXCHNG_MN
	
	# Indirect effect of soil
	IE_pH_PC1 := T1*T1S1
	IE_MN_PC1 := T1*T1S2
	IE_pH_PIE := T2*T2S1
	IE_MN_PIE := T2*T2S2

	IE_pH := IE_pH_PC1 + IE_pH_PIE
	IE_MN := IE_MN_PC1 + IE_MN_PIE
'
soilmodpH2 = '

	# Effect of soil and trees on lichens
	lichen.rich_log ~ T1*PC1 + T2*PIE.ba.tree + S1*PH_H2O + S2*EXCHNG_MN

	# Effect of soil on trees
	PC1 ~ T1S1*PH_H2O + T1S2*EXCHNG_MN
	PIE.ba.tree ~ T2S1*PH_H2O + T2S2*EXCHNG_MN
	
	# Indirect effect of soil
	IE_pH_PC1 := T1*T1S1
	IE_MN_PC1 := T1*T1S2
	IE_pH_PIE := T2*T2S1
	IE_MN_PIE := T2*T2S2

	IE_pH := IE_pH_PC1 + IE_pH_PIE
	IE_MN := IE_MN_PC1 + IE_MN_PIE
'
soilfitpH = sem(soilmodpH2, data=working_soil_fit)
summary(soilfitpH, rsq=TRUE, fit.measures=T)
semPaths(soilfitpH, what='std', whatLabel='std', layout='circle2', nCharNodes=0)

##########################################################
### Fit full model with soil variables

# Center and scale soil data
working_soil = working_soil[,colnames(working_soil)!='yrplot.id']
working_soil = data.frame(scale(working_soil, center=T, scale=T))

testplots = testplots$yrplot.id
working_soil_test = working_soil[testplots[testplots %in% rownames(working_soil)],]

# Remove NAs
working_soil_test = na.omit(working_soil_test) # 291 total plots, 11 had missing soil data

# Write out dataset to be used when fitting model
write.csv(working_soil_test, './Soil/working_soil_test.csv', row.names=T)


# See script sem_boot_soilmod_allsp.R for model that was run on Kure

path_soilmod = "

	# Latent variables
	lichen_rich =~ sqrt(0.75)*lichen.rich_log
	
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
		r1cm1*wetness + r1cm2*rain_lowRH + r1cm3*iso + r1cm4*pseas + r1cm5*mat + r1cm6*radiation + r1sm1*soilPC1 + r1sm2*soilPC2 +
		+ r1R1*regS + r1p1*totalNS +	r1a1*tot_abun_log
	
	# Effects on local lichen abundance
	tot_abun_log ~ a1fh1*bark_moist_pct.rao.ba + a1fh2*wood_SG.rao.ba + a1fh3*LogSeed.rao.ba +
		a1fh4*PIE.ba.tree + a1fh5*propDead + a1fh6*lightDist.mean + a1fh7*diamDiversity +
		a1fm1*bark_moist_pct.ba + a1fm2*wood_SG.ba + a1fm3*LogSeed.ba + a1fm4*bigTrees + a1fm5*light.mean + a1fm6*PC1 +
		a1cm1*wetness + a1cm2*rain_lowRH + a1cm3*iso + a1cm4*pseas + a1cm5*mat + a1cm6*radiation + 
		a1sm1*soilPC1 + a1sm2*soilPC2 + a1p1*totalNS 

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

soilmod_fit =  sem(path_soilmod, data=working_soil_fit, fixed.x=T, estimator='ML', se='robust.sem')
use_fit = soilmod_fit

mod_boot = bootstrapLavaan(use_fit, R=4, FUN=function(x) c(parameterEstimates(x)$est,standardizedSolution(x)$est.std))
soilmod_boot = mod_boot

# Save raw bootstrap output and models
response = 'AllSp'
modform = 'soilmod'

summary(soilmod_fit, standardized=T, rsq=TRUE, fit.measures=T)
soilmod_ests = parameterEstimates(endfit, standardized=T, ci=T, level=0.95)
subset(soilmod_ests, lhs=='lichen_rich'&op=='~')

# Read in results from bootstrapped parameter estimates calculated on the Kure computing cluster




############################



# Indirect effect of soil mainly occur through forest mean conditions and lichen abundance
endfit_ests[grep('IE_soil', endfit_ests$label),]

save.image('./Soil/soil_model.Rdata')

# Ran this function on cluster with R=10000
endfit_std = bootstrapLavaan(endfit, R=5, FUN=function(x) c(parameterEstimates(x)$est,standardizedSolution(x)$est.std))

# Save bootstrap output
write.csv(endfit_std, 'boot_soil.csv', row.names=F)

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
		'bark_moist_pct.ba','wood_SG.ba','LogSeed.ba','bigTrees','light.mean','PC1','totalNS'), sep='_'),
	paste('IE',c('wetness','rain_lowRH','iso','pseas','mat','radiation','soilPC1','soilPC2'),'A', sep='_'))->indir_vars
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

# Make a table of indirect effect of climate and local environment via forest structure
indirect_for = endfit_ests[grep('IE_[A-za-z._1-9]*_F[MH]', endfit_ests$label),]
indirect_for$Ftype = sapply(indirect_for$label, function(x) substr(x, nchar(x)-1, nchar(x)))
indirect_for = indirect_for[,c('lhs','Ftype','std.all','std.se','std.ci.lower','std.ci.upper')]
names(indirect_for)[1] = 'predictor'
indirect_for$predictor = substring(indirect_for$predictor, first=4)
rownames(indirect_for) = indirect_for$predictor
indirect_for$predictor = sapply(indirect_for$predictor, function (x) substr(x, 1, nchar(x)-3))

# Append variable types column
vartypes = read.csv('../SEM models/var_types.csv', row.names=1)
total$type = vartypes[rownames(total),'type']
direct_rich$type = vartypes[rownames(direct_rich),'type']
direct_abun$type = vartypes[rownames(direct_abun),'type']
indirect$type = vartypes[rownames(indirect),'type']
indirect_for$type = vartypes[indirect_for$predictor, 'type']

# Save tables
write.csv(endfit_ests, 'Soil_testdata_parameterEstimates.csv')
write.csv(total, 'Soil_testdata_totaleffects.csv', row.names=T)
write.csv(direct_rich, 'Soil_testdata_directeffects_richness.csv', row.names=T)
write.csv(direct_abun, 'Soil_testdata_directeffects_abundance.csv', row.names=T)
write.csv(indirect, 'Soil_testdata_indirecteffects_via_abundance.csv', row.names=T)
write.csv(indirect_reg, 'Soil_testdata_indirecteffects_via_regS.csv', row.names=T)
write.csv(indirect_for, 'Soil_testdata_indirecteffects_via_forest.csv', row.names=T)








