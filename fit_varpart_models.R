## This script performs the variation partitioning analysis on lichen diversity for FIA plots.


source('./UNC/Projects/FIA Lichen/GitHub/FIA-Lichens/load_data.R')
source('./GitHub/FIA-Lichens/fia_lichen_analysis_functions.R')

#save.image('varpart_analysis.Rdata')
load('varpart_analysis.Rdata')

# Color ramp
mycol = read.csv('C:/Users/jrcoyle/Documents/UNC/Projects/blue2red_10colramp.txt')
mycol = apply(mycol,1,function(x) rgb(x[1],x[2],x[3],maxColorValue=256))
mycol = mycol[10:1]

###################################################################################
### Variation Partitioning among climate, forest, regional richness, and pollution

library(hier.part) # combos
library(MuMIn) # r.squaredLR
library(MASS) # glm.nb


# Read in table of predictor variable types
predtypes = read.csv('predictors.csv', row.names=1)

# Define dataset- use trans_data where variables have been log or sqrt transformed but not scaled
use_data = trans_data[,colnames(trans_data) %in% rownames(predtypes)]

# Add other variables used to analyze FD and restricted taxonomic groups
other_data = trans_data[,c('lichen.rich','fric','raoQ','Parmeliaceae','Physciaceae',
	'tot_abun_log','parm_abun_log','phys_abun_log','regS','regParm','regPhys')]

# Define predictor variable and appropriate abundance and regional richness variables
use_response = 'lichen.rich' # This changes based on the analysis
use_reg = 'regS'
use_abun = 'tot_abun_log'
use_data$richness = other_data[,use_response]
use_data$reg = other_data[,use_reg]
use_data$abun_log = other_data[,use_abun]

# Append soil PCA variables
soil = read.csv('./Soil/soil_PCA.csv', row.names=1)
use_data$soilPC1 = soil[rownames(use_data),'PC1']
use_data$soilPC2 = soil[rownames(use_data), 'PC2']

# Testing and training data sets
# Initially use fitplots to determine whether to use quadratic relationships.
# Then, re-run final analysis using only use testing data set
use_data_test = use_data[testplots$yrplot.id,]
use_data_fit = use_data[fitplots$yrplot.id,]

sum(is.na(use_data_fit)) # Checking for NAs- in soil variables

### Test linear, log-linear, Poisson, and Negative Binomial GLMs with different link functions
# define set of variables to consider in models
use_vars = rownames(subset(predtypes, !(label %in% c('','a1','r1'))))
use_vars = use_vars[-grep('soil',use_vars)] # Leave out soil variables

pois_log_mod = glm(richness~., family=poisson(link='log'), data=use_data_fit[,c('richness',use_vars)])
#pois_iden_mod = glm(richness~., family=poisson(link='identity'), start=rep(1,23), data=use_data_fit[,c('richness',use_vars)]) #warnings- may not have fit correctly
nb_log_mod = glm.nb(richness~., link='log', data=use_data_fit[,c('richness',use_vars)])
#nb_iden_mod = glm.nb(richness~., link='identity', init.theta=10, data=use_data_fit[,c('richness',use_vars)])
gaus_log_mod = glm(richness~., family=gaussian(link='log'), data=use_data_fit[,c('richness',use_vars)])
gaus_iden_mod = glm(richness~., family=gaussian(link='identity'), data=use_data_fit[,c('richness',use_vars)])

AIC(gaus_iden_mod, pois_log_mod, nb_log_mod) #gaus_log_mod,
# nb_log mod wins by far for lichen richness, Parmeliaceae, Physcaceae
# gaus_iden wins by far for fric, can't use nb or pois b/c not integer.

### Test linear, log-linear, Poisson, and Negative Binomial GLMs for abundance
abun_pois_log = glm(richness~abun_log, family=poisson(link='log'), data=use_data_fit)
#abun_pois_iden = glm(richness~abun_log, family=poisson(link='identity'), start=rep(1,2), data=use_data_fit) 
abun_nb_log = glm.nb(richness~abun_log, link='log', data=use_data_fit) 
#abun_nb_iden = glm.nb(richness~abun_log, link='identity', init.theta=1, data=use_data_fit)
abun_gaus_log = glm(richness~abun_log, family=gaussian(link='log'), data=use_data_fit)
abun_gaus_iden = glm(richness~abun_log, family=gaussian(link='identity'), data=use_data_fit)

AIC(abun_pois_log, abun_nb_log, abun_gaus_iden) #, abun_gaus_log

# calculate correlation coef for abundance for AllSp and Parm
abun_mod=glm.nb(richness~abun_log, link='log', data=use_data_test)
summary(abun_mod)
r.squaredLR(abun_mod)

# calculate correlation coef for abundance for Phys
abun_mod=glm(richness~abun_log, family=poisson(link='log'), data=use_data_test)
summary(abun_mod)
r.squaredLR(abun_mod)

use_x = seq(min(use_data$abun_log), max(use_data$abun_log), length.out=100)
pois_y = predict(abun_pois_log, list(abun_log=use_x), type='response')
gauslog_y = predict(abun_gaus_log, list(abun_log=use_x), type='response')
gausiden_y = predict(abun_gaus_iden, list(abun_log=use_x), type='response')
nblog_y = predict(abun_nb_log, list(abun_log=use_x), type='response')

plot(richness~abun_log, data=use_data_fit)
lines(use_x, pois_y, col='red', lwd=2)
lines(use_x, gauslog_y, col='blue', lwd=2)
lines(use_x, gausiden_y, col='green', lwd=2)
lines(use_x, nblog_y, col='yellow', lwd=2)

plot(log10(richness)~abun_log, data=use_data_fit, ylim=c(0,2))
lines(use_x, log10(pois_y), col='red', lwd=2)
lines(use_x, log10(gauslog_y), col='blue', lwd=2)
lines(use_x, log10(gausiden_y), col='green', lwd=2)
lines(use_x, log10(nblog_y), col='yellow', lwd=2)

## Test gaussian, Poisson, NB for regS
pois_log = glm(richness~reg, family=poisson(link='log'), data=use_data_fit)
pois_iden = glm(richness~reg, family=poisson(link='identity'), data=use_data_fit)
nb_log = glm.nb(richness~reg, link='log', data=use_data_fit)
nb_iden = glm.nb(richness~reg, link='identity', data=use_data_fit)
gaus_log = glm(richness~reg, family=gaussian(link='log'), data=use_data_fit)
gaus_iden = glm(richness~reg, family=gaussian(link='identity'), data=use_data_fit)

summary(gaus_iden)

AIC(pois_log, pois_iden, nb_log, nb_iden, gaus_iden) # Ended up using nb_log for consistency with other variables
AIC(pois_log, nb_log, gaus_iden)

# Plot nb_log and nb_iden
use_x = seq(0, max(use_data$reg), length.out=200)
log_y = predict(nb_log, list(reg=use_x), type='response')
iden_y = predict(nb_iden, list(reg=use_x), type='response')
plot(richness~reg, data=use_data_fit, xlim=c(0,220))
lines(use_x, log_y, lwd=2, col='blue')
lines(use_x, iden_y, lwd=2, col='red')
r.squaredLR(nb_log, null=glm.nb(richness~1, data=use_data_fit))
r.squaredLR(nb_iden, null=glm.nb(richness~1, data=use_data_fit))


### Univariate models ###

# For lichen richness
use_vars = rownames(subset(predtypes, !(label %in% c('','r1'))))
unimods = sapply(use_vars, function(x){
	
	# Define non-NA observations
	use_obs = rownames(use_data_fit)[!is.na(use_data_fit[,x])]
	
	# By default uses log link function, so coefficients are on log scale
	linear_mod = glm.nb(use_data_fit[use_obs,'richness']~use_data_fit[use_obs,x])
	quad_mod = 	glm.nb(use_data_fit[use_obs,'richness']~use_data_fit[use_obs,x]+I(use_data_fit[use_obs,x]^2))
	null_mod = glm.nb(richness~1, data=use_data_fit[use_obs,])	

	linear_sum = summary(linear_mod)$coefficients
	quad_sum = summary(quad_mod)$coefficients

	c(AIC(linear_mod), AIC(quad_mod), 
		r.squaredLR(linear_mod, null=null_mod), r.squaredLR(quad_mod, null=null_mod),
		coef(linear_mod),linear_mod$theta, linear_mod$SE.theta,
		coef(quad_mod), quad_mod$theta, quad_mod$SE.theta, 
		length(use_obs)
	)
})

unimods = data.frame(t(unimods))
names(unimods) = c('AIC_line','AIC_quad','R2_line','R2_quad',
	'line_int','line_slope','line_theta','line_theta_SE',
	'quad_int','quad_slope','quad_sq','quad_theta','quad_theta_SE','N')

write.csv(unimods, 'univariate_models_Phys.csv', row.names=T)


# For functional diversity
unimods = sapply(rownames(predtypes), function(x){
	
	# Define non-NA observations
	use_obs = rownames(use_data_fit)[!is.na(use_data_fit[,x])]
	
	linear_mod = lm(use_data_fit[use_obs,'richness']~use_data_fit[use_obs,x])
	quad_mod = 	lm(use_data_fit[use_obs,'richness']~use_data_fit[use_obs,x]+I(use_data_fit[use_obs,x]^2))

	linear_sum = summary(linear_mod)$coefficients
	quad_sum = summary(quad_mod)$coefficients

	c(AIC(linear_mod), AIC(quad_mod), 
		summary(linear_mod)$r.squared, summary(quad_mod)$r.squared,
		coef(linear_mod),	coef(quad_mod),  
		length(use_obs)
	)

})

unimods = data.frame(t(unimods))
names(unimods) = c('AIC_line','AIC_quad','R2_line','R2_quad',
	'line_int','line_slope','quad_int','quad_slope','quad_sq','N')

write.csv(unimods, 'univariate_models_fric.csv', row.names=T)

## Examine relationships among regional scale variables and regional richness
use_vars = rownames(subset(predtypes, !(label %in% c('','R1'))&scale=='regional'))

# Plot effects on regional richness
pdf('./Figures/regional vars vs regS.pdf', height=5, width=5)
par(mar=c(4,4,1,1))
for(x in use_vars){
	plot(use_data$reg~use_data[,x], ylab='Regional richness', xlab=x)
}	
dev.off()

unimods_reg = sapply(use_vars, function(x){
	
	# Define non-NA observations
	use_obs = rownames(use_data_fit)[!is.na(use_data_fit[,x])]
	
	# By default uses log link function, so coefficients are on log scale
	linear_mod = lm(use_data_fit[use_obs,'reg']~use_data_fit[use_obs,x])
	quad_mod = 	lm(use_data_fit[use_obs,'reg']~use_data_fit[use_obs,x]+I(use_data_fit[use_obs,x]^2))

	linear_sum = summary(linear_mod)$coefficients
	quad_sum = summary(quad_mod)$coefficients

	c(AIC(linear_mod), AIC(quad_mod), 
		summary(linear_mod)$r.squared, summary(quad_mod)$r.squared,
		coef(linear_mod),	coef(quad_mod),  
		length(use_obs)
	)
})

unimods_reg = data.frame(t(unimods_reg))
names(unimods_reg) = c('AIC_line','AIC_quad','R2_line','R2_quad',
	'line_int','line_slope','quad_int','quad_slope','quad_sq','N')

write.csv(unimods_reg, 'univariate_models_regS_Phys.csv', row.names=T)


## Make a chart of linear vs quadratic and concavity

modcompare = data.frame(deltaAIC = unimods$AIC_line-unimods$AIC_quad) # deltaAIC neg indicates AIC_line < AIC_quad and quadratic term not needed
rownames(modcompare) = rownames(unimods)
modcompare$type = ifelse(modcompare$deltaAIC>2, 'quadratic','linear')
modcompare$coef_sq = unimods$quad_sq
modcompare$concavity = ifelse(unimods$quad_sq>0, 'up','down')
modcompare = cbind(modcompare, predtypes[rownames(modcompare),c('type','scale','mode')])

modcompare = data.frame(predictor=varnames[rownames(modcompare),'midName'], modcompare)
modcompare = modcompare[order(modcompare$scale, modcompare$mode, modcompare$deltaAIC),]

modcompare_reg = data.frame(deltaAIC = unimods_reg$AIC_line-unimods_reg$AIC_quad) # deltaAIC neg indicates AIC_line < AIC_quad and quadratic term not needed
rownames(modcompare_reg) = rownames(unimods_reg)
modcompare_reg$type = ifelse(modcompare_reg$deltaAIC>2, 'quadratic','linear')
modcompare_reg$coef_sq = unimods_reg$quad_sq
modcompare_reg$concavity = ifelse(unimods_reg$quad_sq>0, 'up','down')
modcompare_reg = cbind(modcompare_reg, predtypes[rownames(modcompare_reg),c('type','scale','mode')])

modcompare_reg = data.frame(predictor=varnames[rownames(modcompare_reg),'midName'], modcompare_reg)
modcompare_reg = modcompare_reg[order(modcompare_reg$scale, modcompare_reg$mode, modcompare_reg$deltaAIC),]

# Note: even though coefficients are on log scale, the quadratic models will be 
# concave down whenever the quadratic coefficient is negative- proved using calculus.

write.csv(modcompare, 'Univariate model shapes GLM-NB.csv', row.names=T)
write.csv(modcompare_reg, 'Univariate model shapes of regS.csv', row.names=T)
write.csv(modcompare, 'Univariate model shapes fric.csv', row.names=T)
write.csv(modcompare, 'Univariate model shapes GLM-NB Phys.csv', row.names=T)
write.csv(modcompare_reg, 'Univariate model shapes of regS Phys.csv', row.names=T)

modcompare = read.csv('Univariate model shapes GLM-NB.csv', row.names=1)
modcompare_reg = read.csv('Univariate model shapes of regS.csv', row.names=1)
modcompare = read.csv('Univariate model shapes fric.csv', row.names=1)
modcompare = read.csv('Univariate model shapes GLM-NB Phys.csv', row.names=1)
modcompare_reg = read.csv('Univariate model shapes of regS Phys.csv', row.names=1)


## Calculate r for local ~ regional richness and richness~abundance
cor(use_data_test$richness,use_data_test$reg)
cor(use_data_test$richness,use_data_test$abun_log)


# Plotting some models
mod2 = glm.nb(richness~wetness+I(wetness^2), data=use_data_fit)
coef2 = coef(mod2)
x = seq(-4,3, .1)

y_pred = predict(mod2, list(wetness=x))
y_calc = coef2[1] + coef2[2]*x + coef2[3]*x^2

y_predR = predict(mod2, list(wetness=x), type='response')
y_calcR = exp(coef2[1] + coef2[2]*x + coef2[3]*x^2)

plot(y_predR, y_calcR)

plot(richness~wetness, data=use_data_fit)
lines(y_predR~x, col='red')
lines(y_calcR~x, col='blue')

# Which variables have AIC supported concave-down relationships?
sq_vars = rownames(subset(modcompare, concavity=='down'&type=='quadratic'))
sq_vars_reg = rownames(subset(modcompare_reg, concavity=='down'&type=='quadratic'))

# Use quadratic relationships for concave-down AIC supported models in variation partitioning
sq_df = use_data[,unique(c(sq_vars, sq_vars_reg))]^2
colnames(sq_df) = paste(colnames(sq_df),'2', sep='')

# Add square
use_data = cbind(use_data, sq_df)

## Use test plots for model results:
use_data_test = use_data[testplots$yrplot.id,]

##################################################
### Variation Partitioning by Local / Regional ###

# define set of predictors to be used in models
use_preds = subset(predtypes, !(label %in% c('','r1','a1','p1','P1')))

regionvars = rownames(use_preds)[use_preds$scale=='regional']
localvars = rownames(use_preds)[use_preds$scale=='local']
localvars = localvars[-grep('soil', localvars)] # Leave out soil vars

region_mod = glm.nb(richness~., data=use_data_test[,c('richness', regionvars, paste(sq_vars[sq_vars %in% regionvars],2,sep=''))])
local_mod = glm.nb(richness~., data=use_data_test[,c('richness', localvars, paste(sq_vars[sq_vars %in% localvars],2,sep=''))])

AIC(region_mod, local_mod) # region_mod is best

# Make a list of predictors
predlist = list(Regional=names(region_mod$coefficients[-1]),
	Local=names(local_mod$coefficients[-1])
)

# Define null model
null_mod = glm.nb(richness~1, link='log', data=use_data_test)

# Calculate pseudo R2 for each model containing successive subsets of variables
apply(combos(2)$ragged, 1, function(x){
	use_vars = unlist(predlist[x])
	this_mod = glm.nb(richness~., data=use_data_test[,c('richness',use_vars)], link='log')
	r2 = r.squaredLR(this_mod, null=null_mod)
	attr(r2, 'adj.r.squared')
})->Rs

names(Rs) = names(predlist)

local_regional_partition = partvar2(Rs)

full_mod = glm.nb(richness~., data=use_data_test[,c('richness',unlist(predlist))], link='log')

# For variation partitioning with soil vars
use_plots = !is.na(use_data_test$soilPC1)
region_mod = glm.nb(richness~., data=use_data_test[use_plots,c('richness', regionvars, paste(sq_vars[sq_vars %in% regionvars],2,sep=''))])
local_mod = glm.nb(richness~., data=use_data_test[use_plots,c('richness', localvars, paste(sq_vars[sq_vars %in% localvars],2,sep=''))])
null_mod = glm.nb(richness~1, link='log', data=use_data_test[use_plots,])
predlist = list(Regional=names(region_mod$coefficients[-1]),
	Local=names(local_mod$coefficients[-1])
)
apply(combos(2)$ragged, 1, function(x){
	use_vars = unlist(predlist[x])
	this_mod = glm.nb(richness~., data=use_data_test[use_plots,c('richness',use_vars)], link='log')
	r2 = r.squaredLR(this_mod, null=null_mod)
	attr(r2, 'adj.r.squared')
})->Rs

## Local-regional partitioning without climate variables (which may inflate joint variation explained)

climatevars = rownames(use_preds)[grep('mat|iso|pseas|wetness|rain|totalNS', rownames(use_preds))]
regNCvars = regionvars[!(regionvars %in% climatevars)]
locNCvars = localvars[!(localvars %in% climatevars)]

regNC_mod = glm.nb(richness~., data=use_data_test[,c('richness', regNCvars, paste(sq_vars[sq_vars %in% regNCvars],2,sep=''))])
locNC_mod = glm.nb(richness~., data=use_data_test[,c('richness', locNCvars, paste(sq_vars[sq_vars %in% locNCvars],2,sep=''))])

AIC(region_mod, local_mod, regNC_mod, locNC_mod) # region_mod is best

# Make a list of predictors
predlist = list(Regional=names(regNC_mod$coefficients[-1]),
	'Local'=names(locNC_mod$coefficients[-1])
)

# Define null model
null_mod = glm.nb(richness~1, link='log', data=use_data_test)

# Calculate pseudo R2 for each model containing successive subsets of variables
apply(combos(2)$ragged, 1, function(x){
	use_vars = unlist(predlist[x])
	this_mod = glm.nb(richness~., data=use_data_test[,c('richness',use_vars)], link='log')
	r2 = r.squaredLR(this_mod, null=null_mod)
	attr(r2, 'adj.r.squared')
})->Rs

names(Rs) = names(predlist)

noclimate_partition = partvar2(Rs)

### Variation Partitioning by Heterogeneity / Optimality ###

## At regional scale

# Assess spatial auto-correlation of regional variables
library(gstat); library(sp)
fitdata_sp =  use_data[fitplots$yrplot.id,]
coordinates(fitdata_sp) = master[rownames(fitdata_sp),c('LON','LAT')]
proj4string(fitdata_sp) = CRS("+proj=longlat")
## THIS PART IS NOT DONE
## MAY NEED TO REVISE MODELS TO ACCOUNT FOR SPATIAL STRUCTURE

RHvars = rownames(subset(use_preds, scale=='regional'&mode=='het'))
ROvars = rownames(subset(use_preds, scale=='regional'&mode=='opt'))

RH_mod = lm(reg~., data=use_data_test[,c('reg', RHvars, paste(sq_vars_reg[sq_vars_reg %in% RHvars],2,sep=''))])
RO_mod = lm(reg~., data=use_data_test[,c('reg', ROvars)]) # Don't need to add in sq_vars_reg since no longer including total_NS_r

# Make a list of predictors
predlist = list(Heterogeneity=names(RH_mod$coefficients[-1]),
	Optimality=names(RO_mod$coefficients[-1])
)

# Calculate adjusted R2 for each model containing successive subsets of variables
apply(combos(2)$ragged, 1, function(x){
	use_vars = unlist(predlist[x])
	this_mod = lm(reg~., data=use_data_test[,c('reg',use_vars)])
	summary(this_mod)$adj.r.squared
})->Rs

names(Rs) = names(predlist)

regional_het_opt_partition = partvar2(Rs)

regional_full_mod = lm(reg~., data=use_data_test[,c('reg',unlist(predlist))])


## At local scale
LHvars = rownames(subset(use_preds, scale=='local'&mode=='het'))
LOvars = rownames(subset(use_preds, scale=='local'&mode=='opt'))
LOvars = LOvars[-grep('soil', LOvars)] # Don't include soil

LH_mod = glm.nb(richness~., data=use_data_test[,c('richness', LHvars, paste(sq_vars[sq_vars %in% LHvars],2,sep=''))])
LO_mod = glm.nb(richness~., data=use_data_test[,c('richness', LOvars, paste(sq_vars[sq_vars %in% LOvars],2,sep=''))])

# Make a list of predictors
predlist = list(Heterogeneity=names(LH_mod$coefficients[-1]),
	Optimality=names(LO_mod$coefficients[-1])
)

# Define null model
null_mod = glm.nb(richness~1, link='log', data=use_data_test)

# Calculate pseudo R2 for each model containing successive subsets of variables
apply(combos(2)$ragged, 1, function(x){
	use_vars = unlist(predlist[x])
	this_mod = glm.nb(richness~., data=use_data_test[,c('richness',use_vars)], link='log')
	r2 = r.squaredLR(this_mod, null=null_mod)
	attr(r2, 'adj.r.squared')
})->Rs

names(Rs) = names(predlist)

local_het_opt_partition = partvar2(Rs)

# For variation partitioning with soil vars
LH_mod = glm.nb(richness~., data=use_data_test[use_plots, c('richness', LHvars, paste(sq_vars[sq_vars %in% LHvars],2,sep=''))])
LO_mod = glm.nb(richness~., data=use_data_test[use_plots, c('richness', LOvars, paste(sq_vars[sq_vars %in% LOvars],2,sep=''))])
predlist = list(Heterogeneity=names(LH_mod$coefficients[-1]),
	Optimality=names(LO_mod$coefficients[-1])
)
apply(combos(2)$ragged, 1, function(x){
	use_vars = unlist(predlist[x])
	this_mod = glm.nb(richness~., data=use_data_test[use_plots,c('richness',use_vars)], link='log')
	r2 = r.squaredLR(this_mod, null=null_mod)
	attr(r2, 'adj.r.squared')
})->Rs


### Variation partitioning for Fric ###
niche_mod_fric = lm(richness~., data=use_data_test[,c('richness', nichevars, colnames(sq_df)[sq_vars %in% nichevars])])
res_mod_fric = lm(richness~., data=use_data_test[,c('richness', resvars, colnames(sq_df)[sq_vars %in% resvars])])

predlist = list(Heterogeneity=names(niche_mod$coefficients[-1]),
	Optimality=names(res_mod$coefficients[-1])
)

# Calculate R2 for each model containing successive subsets of variables
apply(combos(2)$ragged, 1, function(x){
	use_vars = unlist(predlist[x])
	this_mod = lm(richness~., data=use_data_test[,c('richness',use_vars)])
	summary(this_mod)$adj.r.squared
})->Rs

names(Rs) = names(predlist)

niche_res_partition_fric = partvar2(Rs)

## Plot variation partioning analysis in three panels
barwide=1.5
use_shade = c('FF','55','99','CC')
use_color = c('#2415B0','#00BF32','#126A71')

svg('./Figures/variation partitioning figure phys.svg', height=5, width=10)
	par(mar=c(0,6,1.5,0))

	# Create plotting window
	plot(1,1, xlim=c(0,5.5), ylim=c(0,1), axes=F, ylab='', xlab='', type='n', cex.lab=2)
	
	## Add background lines
	usr = par('usr')
	abline(h=seq(0,1,.1), lwd=2, col='grey70', lty=3)	

	## Background bars
	rect(0,0,barwide,1,col='white', lwd=3)	
	rect(2,0,2+barwide,1,col='white', lwd=3)
	rect(4,0,4+barwide,1,col='white', lwd=3)

	## Local-Regional Model
	# Add rectangle for first component
	rect(0,0,barwide, sum(local_regional_partition[1]), lwd=3, col=paste(use_color[1],use_shade[4],sep=''))
	# Add rectangle for second component
	rect(0,sum(local_regional_partition[1:2]), barwide,sum(local_regional_partition[1:3]), lwd=3, col=paste(use_color[2],use_shade[4],sep=''))
	# Add rectangle for overlap
	rect(0,local_regional_partition[1],barwide,sum(local_regional_partition[1:2]), lwd=3, col=paste(use_color[3],use_shade[4],sep=''))

	## Regional Heterogeneity-Optimality
	# Add rectangle for first component
	rect(2,0,2+barwide, sum(regional_het_opt_partition[1]), lwd=3, col=paste(use_color[1],use_shade[2],sep=''))
	# Add rectangle for second component
	rect(2,sum(regional_het_opt_partition[1:2]), 2+barwide,sum(regional_het_opt_partition[1:3]), lwd=3, col=paste(use_color[1],use_shade[3],sep=''))
	# Add rectangle for overlap
	rect(2,regional_het_opt_partition[1], 2+barwide,sum(regional_het_opt_partition[1:2]), lwd=3, col=paste(use_color[1],use_shade[4],sep=''))

	## Local Heterogeneity-Optimality
	# Add rectangle for first component
	rect(4,0,4+barwide, sum(local_het_opt_partition[1]), lwd=3, col=paste(use_color[2],use_shade[2],sep=''))
	# Add rectangle for second component
	rect(4,sum(local_het_opt_partition[1:2]), 4+barwide,sum(local_het_opt_partition[1:3]), lwd=3, col=paste(use_color[2],use_shade[3],sep=''))
	# Add rectangle for overlap
	rect(4,local_het_opt_partition[1], 4+barwide,sum(local_het_opt_partition[1:2]), lwd=3, col=paste(use_color[2],use_shade[4],sep=''))

	## Add axis
	axis(2, las=1, cex.axis=2, lwd=3)
	mtext('Variation Explained', 2, 4, cex=2)
	
	# Add partition labels
	lablocs = sapply(1:4, function(x) sum(local_regional_partition[0:(x-1)])+local_regional_partition[x]/2)
	text(barwide/2, lablocs, labels=names(local_regional_partition)[1:4], cex=2)
	lablocs = sapply(1:4, function(x) sum(regional_het_opt_partition[0:(x-1)])+regional_het_opt_partition[x]/2)
	text(2+barwide/2, lablocs, labels=names(regional_het_opt_partition)[1:4], cex=2)
	lablocs = sapply(1:4, function(x) sum(local_het_opt_partition[0:(x-1)])+local_het_opt_partition[x]/2)
	text(4+barwide/2, lablocs, labels=names(local_het_opt_partition)[1:4], cex=2)
	
	# Add panel text A, B, C
	par(xpd=T)
	text(c(0,2,4),1.05, c('A','B','C'), cex=2, adj=c(0,0))
	par(xpd=F)
dev.off()


## Plot variation partioning analysis for models with soil variables
## models fit to 291 plots in testing data set with soil data
barwide=1.5
mycols = matrix(c('#b3b5ffff','#6b6dd7ff','#8dff94ff','#38af4fff'), nrow=2)
colnames(mycols) = c('regional','local')
rownames(mycols) = c('het','opt')


svg('./Figures/variation partitioning figure soil.svg', height=5, width=8)
	par(mar=c(0,6,1.5,0))

	# Create plotting window
	plot(1,1, xlim=c(0,3.5), ylim=c(0,1), axes=F, ylab='', xlab='', type='n', cex.lab=2)
	
	## Add background lines
	usr = par('usr')
	abline(h=seq(0,1,.1), lwd=2, col='grey70', lty=3)	

	## Background bars
	rect(0,0,barwide,1,col='white', lwd=3)	
	rect(2,0,2+barwide,1,col='white', lwd=3)

	## Local-Regional Model
	# Add rectangle for first component
	rect(0,0,barwide, sum(local_regional_partition[1]), lwd=3, col=mycols['opt','regional'])
	# Add rectangle for second component
	rect(0,sum(local_regional_partition[1:2]), barwide,sum(local_regional_partition[1:3]), lwd=3, col=mycols['opt','local'])
	# Add rectangle for overlap
	rect(0,local_regional_partition[1],barwide,sum(local_regional_partition[1:2]), lwd=3, col='#116870cc')

	## Local Heterogeneity-Optimality
	# Add rectangle for first component
	rect(2,0,2+barwide, sum(local_het_opt_partition[1]), lwd=3, col=mycols['het','local'])
	# Add rectangle for second component
	rect(2,sum(local_het_opt_partition[1:2]), 2+barwide,sum(local_het_opt_partition[1:3]), lwd=3, col=mycols['opt','local'])
	# Add rectangle for overlap
	rect(2,local_het_opt_partition[1], 2+barwide,sum(local_het_opt_partition[1:2]), lwd=3, col='#00bd30cc')

	## Add axis
	axis(2, las=1, cex.axis=2, lwd=3)
	mtext('Variation Explained', 2, 4, cex=2)
	
	# Add partition labels
	lablocs = sapply(1:4, function(x) sum(local_regional_partition[0:(x-1)])+local_regional_partition[x]/2)
	text(barwide/2, lablocs, labels=names(local_regional_partition)[1:4], cex=2)
	lablocs = sapply(1:4, function(x) sum(local_het_opt_partition[0:(x-1)])+local_het_opt_partition[x]/2)
	text(2+barwide/2, lablocs, labels=names(local_het_opt_partition)[1:4], cex=2)
	
	# Add panel text A, B
	par(xpd=T)
	text(c(0,2),1.05, c('A','B'), cex=2, adj=c(0,0))
	par(xpd=F)
dev.off()




## Plot local-regional varition partitioning without climate variables.
svg('./Figures/variation partitioning local-regional no climate.svg', height=5, width=4)
	par(mar=c(0,6,1.5,0))

	# Create plotting window
	plot(1,1, xlim=c(0,2), ylim=c(0,1), axes=F, ylab='', xlab='', type='n', cex.lab=2)
	
	## Background bars
	rect(0,0,barwide,1,col='white', lwd=3)

	## Local-Regional Model
	# Add rectangle for first component
	rect(0,0,barwide, sum(noclimate_partition[1]), lwd=3, col=paste(use_color[1],use_shade[4],sep=''))
	# Add rectangle for second component
	rect(0,sum(noclimate_partition[1:2]), barwide,sum(noclimate_partition[1:3]), lwd=3, col=paste(use_color[2],use_shade[4],sep=''))
	# Add rectangle for overlap
	rect(0,noclimate_partition[1],barwide,sum(noclimate_partition[1:2]), lwd=3, col=paste(use_color[3],use_shade[4],sep=''))

	## Add axis
	axis(2, las=1, cex.axis=2, lwd=3)
	mtext('Variation Explained', 2, 4, cex=2)
	
	# Add partition labels
	lablocs = sapply(1:4, function(x) sum(noclimate_partition[0:(x-1)])+noclimate_partition[x]/2)
	text(barwide/2, lablocs, labels=names(noclimate_partition)[1:4], cex=2)

dev.off()

###################################################################
### Interaction between regional environment and local control of local richness

## How does regional heterogeneity affect regional control of local richness?
reghet_vars = rownames(subset(predtypes, scale=='regional'&mode=='het'&label!=''&type=='env'))

fullmod = glm.nb(richness~mat_reg_var*reg + iso_reg_var*reg + 
	pseas_reg_var*reg + wetness_reg_var*reg + rain_lowRH_reg_var*reg +
	regS_tree*reg, data = use_data_test, link='log')

# Using PCA of heterogeneity variables
reghet_pca  = prcomp(use_data[,reghet_vars], center=T, scale=T)
reghet_pc1 = predict(reghet_pca)[,'PC1']
reghet_pc1 = reghet_pc1[rownames(use_data_test)]

pcamod = glm.nb(richness~reg*reghet_pc1, data=use_data_test)

# Plot local-regional richness colored by regional heterogeneity PC1
colorvec = mycol[cut(reghet_pc1, 10, include.lowest=T)]
mycolbw = c('grey80','black')
colorvecbw = colorRampPalette(mycolbw)(100)[cut(reghet_pc1, 100, include.lowest=T)]

par(mar=c(5,5,1,6))
plot(richness~reg, data=use_data_test, pch=16, col=colorvecbw)

usr = par('usr')
plotColorRamp(cols = mycolbw, n = 100, barends = c(usr[2], usr[3], usr[2]+0.05*diff(usr[1:2]), usr[4]),
	labels = seq(-2.5,3.5,.5), uneven.lab=T, labrange=range(reghet_pc1), title='Regional Heterogeneity (PC1)')



## How does regional heterogeneity affect local control of local richness

use_modstring = paste('richness', paste(paste(localvars, '*reghet_pc1', sep=''), collapse='+'), sep='~')

## Opposite interaction for regional heterogeneity on strength of local filters?
reghetloc_mod = glm.nb(use_modstring, data=use_data_test, link='log')
cbind(coef(reghetloc_mod)[c(2,4:21)], coef(reghetloc_mod)[22:40], coef(summary(reghetloc_mod))[22:40,4])

# Probably should only do this for local scale variables that have significant effects.
# Decide using sem?
use_locvars = c('pseas','wetness','bark_moist_pct.rao.ba','PIE.ba.tree','wood_SG.ba','light.mean','PC1')

interactionmods = sapply(use_locvars, function(x){
	
	# Define non-NA observations
	use_obs = rownames(use_data_test)[!is.na(use_data_test[,x])]
	
	# By default uses log link function, so coefficients are on log scale
	linear_mod = glm.nb(use_data_test[use_obs,'richness']~use_data_test[use_obs,x]*reghet_pc1[use_obs])
	#quad_mod = glm.nb(use_data_fit[use_obs,'richness']~use_data_fit[use_obs,x]+I(use_data_fit[use_obs,x]^2))
	null_mod = glm.nb(richness~1, data=use_data_test[use_obs,])	

	linear_sum = summary(linear_mod)$coefficients
	#quad_sum = summary(quad_mod)$coefficients

	c(AIC(linear_mod),
		r.squaredLR(linear_mod, null=null_mod),
		coef(linear_mod),linear_sum[2:4,4],linear_mod$theta, linear_mod$SE.theta,
		length(use_obs)
	)
})

interactionmods = data.frame(t(interactionmods))
names(interactionmods) = c('AIC_line','R2_line',
	'line_int','line_slope','line_heteffect','line_interaction',
	'line_slope_P','line_heteffect_P','line_interaction_P',
	'line_theta','line_theta_SE','N')




reghetloc_mod = glm.nb(richness~wetness*reghet_pc1 + pseas*reghet_pc1 + PIE.ba.tree*reghet_pc1 + 
	bark_moist_pct.rao.ba*reghet_pc1 + wood_SG.ba*reghet_pc1 + light.mean*reghet_pc1 + 
	PC1*reghet_pc1, data=use_data_test, link='log')




## Negative interaction between regional heterogeneity and local abundance?

reghetabun_mod = glm.nb(richness~abun_log*reghet_pc1, data=use_data_test)
reghetregS_mod = glm.nb(richness~abun_log*reg, data=use_data_test)



## Does regional richness env heterogeneity interaction result from pollution differences between east vs west?

plot(totalNS~totalNS_reg, data=use_data_test)

west = master[rownames(use_data_test),'LON'] < -100


## Pollution may alter local-regional richness relationships
svg('./Figures/pollution and reg het effect on local-regional richness.svg', height=6, width=7)
text.cex=1.2

par(mfrow=c(2,2))
par(las=1)
par(lend=1)
par(cex.axis=text.cex)
par(cex.lab=text.cex)
par(mgp=c(2.4,0.7,0))

# color by regional heterogeneity
colorvec = mycol[cut(reghet_pc1, 10, include.lowest=T)]
par(mar=c(2,7,3,2))
plot(richness~reg, data=use_data_test[west,], pch=16, col=colorvec[west], xlim=c(120,210), ylim=c(0,40), 
	axes=F, xlab='', ylab='Local Species Richness')
axis(1, at=seq(125,200,25))
axis(2)
box()
mtext('West',3,font=1, cex=text.cex, line=.5)

par(mar=c(2,1,3,8))
plot(richness~reg, data=use_data_test[!west,], pch=16, col=colorvec[!west], xlim=c(120,210), ylim=c(0,40),
	axes=F, xlab='', ylab='')
axis(1, at=seq(125,200,25))
axis(2)
box()
mtext('East',3,font=1, cex=text.cex, line=.5)

usr = par('usr')
plotColorRamp(cols = mycol, n = 100, barends = c(usr[2], usr[3], usr[2]+0.05*diff(usr[1:2]), usr[4]),
	labels = seq(-2,3,1), uneven.lab=T, labrange=range(reghet_pc1), title='Regional Heterogeneity (PC1)',
	mycex=text.cex)
axis(1, at=seq(125,200,25))
axis(2)
box()

# color by local pollution
colorvec = mycol[cut(master[rownames(use_data_test),'totalNS'], breaks=seq(0,1000,100), include.lowest=T)]
par(mar=c(4,7,1,2))
plot(richness~reg, data=use_data_test[west,], pch=16, col=colorvec[west], xlim=c(120,210), ylim=c(0,40),
	axes=F, xlab='Regional Species Richness', ylab='Local Species Richness')
axis(1, at=seq(125,200,25))
axis(2)
box()

par(mar=c(4,1,1,8))
plot(richness~reg, data=use_data_test[!west,], pch=16, col=colorvec[!west], xlim=c(120,210), ylim=c(0,40),
	axes=F, xlab='Regional Species Richness', ylab='')
usr = par('usr')
plotColorRamp(cols = mycol, n = 100, barends = c(usr[2], usr[3], usr[2]+0.05*diff(usr[1:2]), usr[4]),
	labels = seq(0,1000,200), ndig=0, title='Total N+S Deposition (eq/ha)',
	mycex=text.cex)
axis(1, at=seq(125,200,25))
axis(2)
box()

dev.off()

polmod1 = glm.nb(richness~reg*totalNS, data=use_data_test)
polmod2 = glm.nb(richness~reg*reghet_pc1+totalNS, data=use_data_test)


colorvec = mycol[cut(master[rownames(use_data_test),'totalNS'], breaks=seq(0,1000,100), include.lowest=T)]
plot(reg~reghet_pc1, data=use_data_test, pch=16, col=colorvec, ylim=c(120,210),
	axes=T, xlab='Regional Heterogeneity', ylab='Regional Species Richness')
usr = par('usr')
plotColorRamp(cols = mycol, n = 100, barends = c(usr[2], usr[3], usr[2]+0.05*diff(usr[1:2]), usr[4]),
	labels = seq(0,1000,100), ndig=0, title='Total N+S Deposition (eq/ha)')




##################################################################
### Old Code


barwide=1.5
use_shade = c('FF','55','99','CC')
use_color = c('#2415B0','#00BF32','#126A71')

svg('./Figures/variation partitioning figure v3 lines Phys.svg', height=5, width=10)
	par(mar=c(0,6,0,0))

	# Create plotting window
	plot(1,1, xlim=c(0,5.5), ylim=c(0,1), axes=F, ylab='', xlab='', type='n', cex.lab=2)
	
	## Add background lines
	usr = par('usr')
	abline(h=seq(0,1,.1), lwd=2, col='grey70', lty=3)	

	## Background bars
	rect(2,0,2+barwide,1,col='white', lwd=3)
	rect(0,0,barwide,sum(regional_partition[1:3]),col='white', lwd=3)	
	rect(4,local_regional_partition[1],4+barwide,sum(local_regional_partition[1:3]),col='white', lwd=3)

	## Local-Regional Model
	# Add rectangle for unexplained variation
	rect(2,sum(local_regional_partition[1:3]), 2+barwide, 1, lwd=3, col='white')	
	# Add rectangle for first component
	rect(2,0,2+barwide, sum(local_regional_partition[1]), lwd=3, col=paste(use_color[1],use_shade[4],sep=''))
	# Add rectangle for second component
	rect(2,sum(local_regional_partition[1:2]), 2+barwide,sum(local_regional_partition[1:3]), lwd=3, col=paste(use_color[2],use_shade[4],sep=''))
	# Add rectangle for overlap
	rect(2,local_regional_partition[1],2+barwide,sum(local_regional_partition[1:2]), lwd=3, col=paste(use_color[3],use_shade[4],sep=''))

	## Regional Model (Climate + RegS)
	# Add rectangle for first component
	rect(0,0,barwide, sum(regional_partition[1]), lwd=3, col=paste(use_color[1],use_shade[2],sep=''))
	# Add rectangle for second component
	rect(0,sum(regional_partition[1:2]), barwide,sum(regional_partition[1:3]), lwd=3, col=paste(use_color[1],use_shade[3],sep=''))
	# Add rectangle for overlap
	rect(0,regional_partition[1], barwide,sum(regional_partition[1:2]), lwd=3, col=paste(use_color[1],use_shade[4],sep=''))

	## Regional Model (Climate + RegS)
	# Add rectangle for first component
	rect(4,local_regional_partition[1], 4+barwide, local_regional_partition[1]+niche_res_partition[1], lwd=3, col=paste(use_color[2],use_shade[2],sep=''))
	# Add rectangle for second component
	rect(4,local_regional_partition[1]+sum(niche_res_partition[1:2]), 4+barwide, local_regional_partition[1]+sum(niche_res_partition[1:3]), lwd=3, col=paste(use_color[2],use_shade[3],sep=''))
	# Add rectangle for overlap
	rect(4,local_regional_partition[1]+niche_res_partition[1], 4+barwide, local_regional_partition[1]+sum(niche_res_partition[1:2]), lwd=3, col=paste(use_color[2],use_shade[4],sep=''))

	## Arrows
	jit = 0.05
	arrows(2-jit, c(0,sum(local_regional_partition[1:2])), barwide+jit, c(0,sum(local_regional_partition[1:2])), lwd=3, length=.1, angle=45)
	arrows(2+barwide+jit, c(local_regional_partition[1],sum(local_regional_partition[1:3])), 4-jit, c(local_regional_partition[1],sum(local_regional_partition[1:3])), lwd=3, length=.1, angle=45)

	## Add axis
	axis(2, las=1, cex.axis=2, lwd=3)
	mtext('Variation Explained', 2, 4, cex=2)
	
	# Add partition labels
	lablocs = sapply(1:4, function(x) sum(local_regional_partition[0:(x-1)])+local_regional_partition[x]/2)
	text(2+barwide/2, lablocs, labels=names(local_regional_partition)[1:4], cex=2)
	lablocs = sapply(1:3, function(x) sum(regional_partition[0:(x-1)])+regional_partition[x]/2)
	text(barwide/2, lablocs, labels=names(regional_partition)[1:3], cex=2)
	lablocs = local_regional_partition[1]+sapply(1:3, function(x) sum(niche_res_partition[0:(x-1)])+niche_res_partition[x]/2)
	text(4+barwide/2, lablocs, labels=names(niche_res_partition)[1:3], cex=2)
dev.off()

### Overlap of individual forest predictors and regional richness

regS_mod = glm.nb(richness~reg, data=use_data_test)
null_mod = glm.nb(richness~1, link='log', data=use_data_test)

# Compare each local predictor and regional richness
sapply(localvars, function(x){
	full_mod = glm.nb(richness~., data=use_data_test[,c('richness','reg',x,colnames(sq_df)[sq_vars==x])], link='log')
	loc_mod = update(full_mod, .~.-reg)

	Rs = sapply(list(regS_mod, loc_mod, full_mod), function(m) attr(r.squaredLR(m, null=null_mod),'adj.r.squared') )
	names(Rs) = c('regS','localvar','Full')

	partvar2(Rs)	
})-> local_regS_varpart
local_regS_varpart = data.frame(t(local_regS_varpart))

# Which local variables have the most joint variation with regional richness?
local_regS_varpart[order(local_regS_varpart$Both),]


## Try removing local predictors to see which variable accounts for unique local variance
remove_pred_eff = rep(NA, length(localvars))
names(remove_pred_eff) = localvars

new_local = update(local_mod, .~.-propDead-propDead2)
new_full = update(full_mod, .~.-propDead-propDead2)
Rs_new = c(region=attr(r.squaredLR(region_mod, null=null_mod), 'adj.r.squared'),
	local=attr(r.squaredLR(new_local, null=null_mod), 'adj.r.squared'),
	attr(r.squaredLR(new_full, null=null_mod), 'adj.r.squared'))
partvar2(Rs_new)['local'] -> remove_pred_eff['propDead']



## COLOR FIGURES

# Plot variation partitioning
barwide = .6
use_part = local_regional_partition

svg('./Figures/partition variance among local vs regional include quadratic effects.svg', height=6, width=6*barwide)
	par(mar=c(0,4,0,0))

	# Create plotting window
	plot(1,1, xlim=c(0,barwide), ylim=c(0,1), axes=F, xlab='', ylab='', type='n')
	
	# Add rectangle for unexplained variation
	rect(0,sum(use_part[1:3]),barwide,1, lwd=3, col='white')
	
	# Add rectangle for regional model
	rect(0,0,barwide, sum(use_part[1:2]), lwd=3, col='#AA000088')

	# Add rectangle for local model
	rect(0,use_part[1],barwide,sum(use_part[1:3]), lwd=3, col='#0000AA88')

	# Add axis
	axis(2, las=1, cex.axis=2, lwd=3)
	
	# Add partition labels
	lablocs = sapply(1:4, function(x) sum(use_part[0:(x-1)])+use_part[x]/2)
	text(barwide/2, lablocs, labels=names(use_part), cex=2)
dev.off()


# Plot forest variable partitioning
barwide=.6
use_part = niche_res_partition
svg('./Figures/partition variance among heterogeneity vs optimality include quadratic effects no unexp.svg', height=6, width=6*barwide)
	par(mar=c(0,0,0,5))

	# Create plotting window
	plot(1,1, xlim=c(0,barwide), ylim=c(0,sum(use_part[1:3])), axes=F, xlab='', ylab='', type='n')
	
	# Add rectangle for regional model
	rect(0,0,barwide, sum(use_part[1:2]), lwd=3, col='#00AA0088')

	# Add rectangle for local model
	rect(0,use_part[1],barwide,sum(use_part[1:3]), lwd=3, col='#0000AA88')

	# Add axis
	axis(4, las=1, cex.axis=2, lwd=3)
	
	# Add partition labels
	lablocs = sapply(1:3, function(x) sum(use_part[0:(x-1)])+use_part[x]/2)
	text(barwide/2, lablocs, labels=names(use_part)[1:3], cex=2)
dev.off()








## Variation partitioning by variable type
## Should I be scaling data if I am going to be comparing variation explained among groups?
## Note: not using this any more. hier.part does not give same partitioning as Borcard 1992.
forestvars = rownames(predtypes)[predtypes$type %in% c('FM','FH')]
envvars = rownames(predtypes)[predtypes$type %in% c('C','L')]
envvars = envvars[-grep('soil', envvars)] # not using soil variables here
regvars = rownames(predtypes)[predtypes$type == 'R']
polvars = rownames(predtypes)[predtypes$type == 'P']

forest_mod = glm.nb(richness~. , data=use_data_test[,c('richness', forestvars, colnames(sq_df)[sq_vars %in% forestvars])])
reg_mod = glm.nb(richness~., data=use_data_test[,c('richness', regvars, colnames(sq_df)[sq_vars %in% regvars])])
env_mod = glm.nb(richness~., data=use_data_test[,c('richness', envvars, colnames(sq_df)[sq_vars %in% envvars])])
pol_mod = glm.nb(richness~., data=use_data_test[,c('richness', polvars, colnames(sq_df)[sq_vars %in% polvars])])
AIC(forest_mod, reg_mod, env_mod, pol_mod) # env_mod is best

# Make a list of predictors
predlist = list(environ=names(env_mod$coefficients[-1]),
	forest=names(forest_mod$coefficients[-1]),
	regional=names(reg_mod$coefficients[-1]),
	pollution=names(pol_mod$coefficients[-1])
)

# Define null model
null_mod = glm.nb(richness~1, data=use_data_test)

# Calculate pseudo R2 for each model containing successive subsets of variables
apply(combos(4)$ragged, 1, function(x){
	use_vars = unlist(predlist[x])

	this_mod = glm.nb(richness~., data=use_data_test[,c('richness',use_vars)])
	r.squaredLR(this_mod, null=null_mod)
})->Rs

# Add the null model
nullR = r.squaredLR(null_mod)
Rs = c(nullR,Rs)

# Partition variation among variable sets (from hier.part package)
partition(Rs, 4, var.names = names(predlist))-> varTypes_partition 

png('./Figures/partition variance among predictor types include quadratic effects barplot.png', height=500, width=750)
par(mar=c(4,6,1,1))
barplot(t(as.matrix(varTypes_partition$IJ[,c('I','J')])),
	legend.text = c('Independent','Joint'),ylim=c(0,1),
	las=1, cex.axis=1.5, cex.lab=1.5, cex.names=1.5,col=c('grey30','grey60'),
	args.legend=list(bty='n', cex=1.5, border='transparent'), mgp=c(3,2,1),
	names=c('Environment','Forest\nStructure','Regional\nRichness','Pollution'), border=F
)
mtext('Variation Explained',2,4.5,cex=1.5)
dev.off()












#### Partition among climate, forest heterogeneity and mean conditions

# Define models for heterogeneity (niche) and mean conditions (res)
nichevars = rownames(predtypes)[predtypes$type == 'FH']
resvars = rownames(predtypes)[predtypes$type == 'FM']
niche_mod = glm.nb(richness~. , data=use_data[,c('richness', nichevars, colnames(sq_df)[sq_vars %in% nichevars])])
res_mod = glm.nb(richness~. , data=use_data[,c('richness', resvars, colnames(sq_df)[sq_vars %in% resvars])])

# Make a list of predictors
predlist = list(climate=names(clim_mod$coefficients[-1]),
	niche=names(niche_mod$coefficients[-1]),
	res=names(res_mod$coefficients[-1])
)

# Calculate pseudo R2 for each model containing successive subsets of variables
apply(combos(2)$ragged, 1, function(x){
	use_vars = unlist(predlist[x])

	this_mod = glm.nb(richness~., data=use_data[,c('richness',use_vars)])
	r.squaredLR(this_mod)
})->Rs


# Partition variation
inds = Rs[7]-Rs[6:4]
joint = -1*(Rs[4:6] - Rs[c(1,1,2)] - Rs[c(2,2,3)])
middle = -1*(Rs[7] - sum(inds) -  sum(joint))/2
joint = joint - middle


# Plot Venn Diagram
weights = c(0,inds[1:2],joint[1],inds[3], joint[2:3],middle)
weights = round(weights, 3)

myVenn = Venn(SetNames=c('Climate','Forest\nHeterogeneity','Forest Mean\nConditions'),
	Weight=weights)

# There is no support for changing colors in R so I will do it in Inkscape
pdf('./Figures/climate forest structure venn diagram.pdf', height=8, width=10)
plot(myVenn, doWeights=T, type='ChowRusky')
dev.off()
