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
library(spdep) # spatial autocorrelation
library(sp) # spatial data handling
library(ape) # Moran.I
library(ncf) # correlog - distance-based binning for Moran's I correlogram
library(gstat) # variogram
library(qpcR) # akaike.weights

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

AIC(gaus_log_mod,gaus_iden_mod, pois_log_mod, nb_log_mod) #gaus_log_mod,
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

AIC(pois_log, pois_iden, nb_log, nb_iden, gaus_log, gaus_iden) # Ended up using nb_log for consistency with other variables
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

write.csv(unimods, 'univariate_models_AllSp.csv', row.names=T)


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

write.csv(unimods_reg, 'univariate_models_regS_AllSp.csv', row.names=T)


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
use_data_fit = use_data[fitplots$yrplot.id,]

##################################################
### Variation Partitioning by Local / Regional ###

# Define spatial dataframe for spatial analysis
spdata = master[rownames(use_data_test),c('LAT','LON','lichen.rich','regS')]
coordinates(spdata) = c('LON','LAT'); proj4string(spdata) = CRS("+proj=longlat")

# Create spatial weight matrix for eigenvector decomposition using inverse-distance weighting
invdist_mat = 1/spDists(spdata, longlat=T)
invdist_mat[is.infinite(invdist_mat)]<-0

# Extract eigenvectors from centered spatial weight matrix (from Griffith and Chun 2014)
B = invdist_mat
n = nrow(B)
M = diag(n) - matrix(1,n,n)/n
MBM = M %*% B %*% M
eig = eigen(MBM, symmetric=T)
hist(eig$values/eig$values[1])
EV = as.data.frame(eig$vectors[,eig$values/eig$values[1]>0.25]) # initial selection of candidate set of vectors 
#EV = as.data.frame(eig$vectors) 
colnames(EV) = paste("EV", 1:ncol(EV), sep="")

# Plot candidate set of eigenvectors
EVsp = EV; coordinates(EVsp) = coordinates(spdata); proj4string(EVsp) = CRS("+proj=longlat")
spplot(EVsp, c('EV1', 'EV2','EV3','EV4'), col.regions=mycol, cuts=10)
spplot(EVsp, c('EV42'), col.regions=mycol, cuts=10)

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

# Add model residuals to spatial data 
spdata$fullmod_res = resid(full_mod)

# Calculate Moran's I based on IDW matrix
Moran.I(spdata$fullmod_res, invdist_mat)
plot(spdata$lichen.rich, spdata$fullmod_res)

# Calculate Moran's I for local ve regional models
# greater residual SA in local model
spdata$loc_res = resid(local_mod)
spdata$reg_res = resid(region_mod)
Moran.I(spdata$loc_res, invdist_mat)
Moran.I(spdata$reg_res, invdist_mat)

# Correlogram of residuals (long computation time!)
cg = correlog(spdata$LON, spdata$LAT, spdata$fullmod_res, latlon=T, increment=50)
plot(cg)
abline(h=0)

# Plot residuals
plot_prj = paste("+proj=laea +lat_0=40 +lon_0=-97 +units=km",sep='')
spdata_prj = spTransform(spdata, CRS(plot_prj))
spplot(spdata_prj, 'fullmod_res')
spplot(spdata_prj, 'loc_res')
spplot(spdata_prj, 'reg_res')

# Selection of eigenvectors to use in spatial variation partitioning
# Brute-force determination of set of models with no RSA (residual spatial autocorrelation)
env = use_data_test[,c('richness',unlist(predlist))]
N = ncol(EV)
models = expand.grid(lapply(1:N, function(x) 0:1))

modelset = data.frame(mod=c(), AIC=c(), MC.obs=c(), MC.p = c())
for(i in 1:10){
	this_spat = EV[,as.logical(models[i,])]
	this_mod = glm.nb(richness~., data=cbind(env, this_spat), link='log')
	this_res = resid(this_mod)
	MC = Moran.I(this_res, B)

	modelset = rbind(modelset, data.frame(mod=i, AIC=AIC(this_mod), MC.obs=MC$observed, MC.p=MC$p.value))
	i
}

save(env, EV, B, file='./SpatialAnalysis/fullmod.RData')
# The loop above was run on Kure and the results are below
modelset = read.csv('./SpatialAnalysis/fullmod_spatEV.csv')

# models with no significant SA
noSAset = subset(modelset, MC.p>0.05)

# order by AIC
noSAset = noSAset[order(noSAset$AIC),]
min(modelset$AIC) # note that best model is in this set

# AIC weights within this set
AICwts = akaike.weights(noSAset$AIC)
AICsum = sapply(1:nrow(noSAset), function(x) sum(AICwts$weights[1:x]))
best_mods = noSAset[AICsum <= 0.9,]

colSums(models[best_mods$mod,])
models[noSAset[1:5,'mod'],]

# Choose to use model 82831, has eigenvectors: 2,3,4,8,9,10,15,16
spatmod_full = models[noSAset$mod[1],]

## Variation partitioning with space

# Add eigenvectors to data
use_data_fit = cbind(use_data_fit, EV[,as.logical(spatmod_full)])

# Make a list of predictors
predlist = list(Regional=names(region_mod$coefficients[-1]),
	Local=names(local_mod$coefficients[-1]),
	Spatial=names(EV[,as.logical(spatmod_full)])
)

# Define null model
null_mod = glm.nb(richness~1, link='log', data=use_data_fit)

# Calculate pseudo R2 for each model containing successive subsets of variables
apply(combos(3)$ragged, 1, function(x){
	use_vars = unlist(predlist[x])
	this_mod = glm.nb(richness~., data=use_data_fit[,c('richness',use_vars)], link='log')
	r2 = r.squaredLR(this_mod, null=null_mod)
	attr(r2, 'adj.r.squared')
})->Rs

names(Rs) = apply(combos(3)$binary, 1, function(x) paste(names(predlist)[as.logical(x)], collapse='+'))
partvar3(Rs)

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
RHvars = rownames(subset(use_preds, scale=='regional'&mode=='het'))
ROvars = rownames(subset(use_preds, scale=='regional'&mode=='opt'))

RH_mod = lm(reg~., data=use_data_test[,c('reg', RHvars, paste(sq_vars_reg[sq_vars_reg %in% RHvars],2,sep=''))])
RO_mod = lm(reg~., data=use_data_test[,c('reg', ROvars, paste(sq_vars_reg[sq_vars_reg %in% ROvars],2,sep=''))]) # Don't need to add in sq_vars_reg since no longer including total_NS_r

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

env = use_data_test[,c('reg',unlist(predlist))]
regional_full_mod = lm(reg~., data=use_data_test[,c('reg',unlist(predlist))])


## Look at spatial autocorrelation of residuals in RH vs RO
# SA is high in both models but of comparable magnitude
spdata$RH_res = resid(RH_mod)
spdata$RO_res = resid(RO_mod)
spdata$regmod_res = resid(regional_full_mod)

# Autocorrelation in regional richness
Moran.I(spdata$regmod_res, invdist_mat)
Moran.I(spdata$RH_res, invdist_mat)
Moran.I(spdata$RO_res, invdist_mat)
plot(spdata$regS, spdata$regmod_res)
spplot(spdata, 'regmod_res', col.regions=mycol, cuts=10)

EV_reg = as.data.frame(eig$vectors[,eig$values/eig$values[1]>0.035]) 
EV_reg= data.frame(eig$vectors[,1:150])
colnames(EV_reg) = paste("EV", 1:ncol(EV_reg), sep="")
reg_spat = lm(reg~., data=cbind(env, EV_reg))
Moran.I(resid(reg_spat), B)

plot(regspatmod_res~regS, data=spdata)

spdata$regspatmod_res = resid(reg_spat)
spplot(spdata, c('regspatmod_res','regmod_res'), col.regions=mycol, cuts=10)
vg = variogram(regspatmod_res~1, coordinates(spdata), data=spdata)
plot(vg)

# MCMC exploration of models
this_vars = rep(T, ncol(EV_reg))
this_mod = lm(reg~., data=cbind(env, EV_reg[,this_vars]))
this_MC = Moran.I(resid(this_mod), B)$p.value
this_AIC = AIC(this_mod)
vartrace = data.frame(this_MC, this_AIC)
modtrace = this_vars
for(i in 1:10){
	changevar = sample.int(ncol(EV_reg), 1)
	new_vars = this_vars
	new_vars[changevar] = ifelse(this_vars[changevar], F, T)
	
	this_mod = lm(reg~., data=cbind(env, EV_reg[,this_vars]))
	this_MC = Moran.I(resid(this_mod), B)$p.value
	this_AIC = AIC(this_mod)
	new_mod = lm(reg~., data=cbind(env, EV_reg[,new_vars]))
	new_MC = Moran.I(resid(new_mod), B)$p.value
	new_AIC = AIC(new_mod)

	if((new_MC>0.05)&(this_AIC-new_AIC > -2)) this_vars = new_vars

	vartrace=rbind(vartrace, c(MC_p = this_MC, AIC = this_AIC))
	modtrace= rbind(modtrace, this_vars)
}

save(EV_reg, env, B, file='./SpatialAnalysis/regmod.RData')
# Run this on the cluster multiple times

stats1 = read.csv('./SpatialAnalysis/regmodMCMC_stats-run1.csv')
mods1 = read.csv('./SpatialAnalysis/regmodMCMC_models-run1.csv')
plot(stats2$this_MC, stats2$this_AIC, type='p')
plot(stats2$this_AIC, type='l')
plot(stats2$this_MC, type='l')
hist(colSums(mods3)/nrow(mods3))

# Set of models with lowest AIC - chose based on models with cumulative weight summing to 0.9
stats1 = stats1[order(stats1$this_AIC),]
mods1 = mods1[order(stats1$this_AIC),]
AICwts = akaike.weights(stats1$this_AIC)
AICwts$weights[1:30]
AICsum = sapply(1:nrow(stats1), function(x) sum(AICwts$weights[1:x]))
best_mods = mods1[AICsum <= 0.9,]
best_stats = stats1[AICsum <= 0.9,]
plot(stats1$this_MC, stats1$this_AIC, type='p', col=c(1,2)[as.numeric(AICsum <= 0.9)+1])

# Which variables are in best models?
best_vars = colSums(best_mods)/nrow(best_mods)
hist(best_vars)
use_vars = best_vars>0.5
try_mod = lm(reg~., data=cbind(env, EV_reg[,use_vars]))
AIC(try_mod)
Moran.I(resid(try_mod), B)

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

# Spatial-autocorrelation local models
env = use_data_test[,c('richness',unlist(predlist))]

EV_loc = as.data.frame(eig$vectors[,eig$values/eig$values[1]>0.15]) 
EV_loc= data.frame(eig$vectors[,1:43])
colnames(EV_loc) = paste("EV", 1:ncol(EV_loc), sep="")
loc_spat = glm.nb(richness~., data=cbind(env, EV_loc), link='log')
Moran.I(resid(loc_spat), B)

# Use script testmodels_locmod.R to use MCMC to explore model sets for best model
save(EV_loc, env, B, file='./SpatialAnalysis/locmod.RData')

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

svg('./Figures/variation partitioning figure new regS.svg', height=5, width=10)
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
svg('./Figures/variation partitioning local-regional no climate new regS.svg', height=5, width=4)
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

use_preds = subset(predtypes, !(label %in% c('','r1','a1','p1','P1')))
RHvars = rownames(subset(use_preds, scale=='regional'&mode=='het'))
ROvars = rownames(subset(use_preds, scale=='regional'&mode=='opt'))

## How does regional heterogeneity affect regional control of local richness?

fullmod = glm.nb(richness~mat_reg_var*reg + iso_reg_var*reg + 
	pseas_reg_var*reg + wetness_reg_var*reg + rain_lowRH_reg_var*reg +
	regS_tree*reg, data = use_data_test, link='log')

# Using PCA of heterogeneity variables
reghet_pca  = prcomp(use_data[,RHvars], center=T, scale=T)
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

par(mar=c(5,5,1,6))
plot(richness~reg, data=use_data_test, pch=16, col=colorvec)

usr = par('usr')
plotColorRamp(cols = mycol, n = 100, barends = c(usr[2], usr[3], usr[2]+0.05*diff(usr[1:2]), usr[4]),
	labels = seq(-2.5,3.5,.5), uneven.lab=T, labrange=range(reghet_pc1), title='Regional Heterogeneity (PC1)')

colorvec = mycol[cut(use_data_test$regS, 10, include.lowest=T)]
plot(richness~reghet_pc1, data=use_data_test, pch=16, col=colorvec)

exp(coef(pcamod))

## Partial regression controlling for mean climate

# See p.570 Legendre & Legendre 2012
# y : local richness
# W : regional climate
# X : regional heterogeneity PC1

W = use_data_test[,c(ROvars, paste(sq_vars_reg[sq_vars_reg %in% ROvars],2,sep=''))]

# Normal model
yW_mod = lm(richness~., data=data.frame(richness=use_data_test$richness, W))
yW_res = resid(yW_mod)
XW_mod = lm(reghet_pc1~., data=W)
XW_res = resid(XW_mod)
yspX_mod = lm(use_data_test$richness~XW_res)
yspX_intmod = lm(richness~regS*XW_res, data=use_data_test)
ypX_mod = lm(yW_res~XW_res)
ypX_intmod = lm(yW_res~regS*XW_res, data=use_data_test)

plot(yW_res~regS, data=use_data_test)

# Plot local ~ regional richness colored by reghet residuals
colorvec = colorRampPalette(mycol)(100)[cut(XW_res, 100, include.lowest=T)]
plot(richness~regS, data=use_data_test,
	xlab='Regional Species Richness', ylab='Local Species Richness', 
	las=1, xlim=c(230, 430), ylim=c(0,40),
	col=colorvec, lwd=2, axes=F)
axis(1, at=seq(230,430,50))
axis(2)
box()

## Map reghet and reghet residuals
# Map projection
plot_prj = paste("+proj=laea +lat_0=40 +lon_0=-97 +units=km",sep='')

# N. Am. outline
OUTLINES = readOGR('../../GIS shape files/N Am Outline','na_base_Lambert_Azimuthal')
OUTLINES.laea = spTransform(OUTLINES,CRS(plot_prj))

het_sp = data.frame(reghet_pc1, XW_res)
coordinates(het_sp) = master[rownames(use_data_test), c('LON','LAT')]
proj4string(het_sp) = CRS("+proj=longlat")
het_laea = spTransform(het_sp, CRS(plot_prj))

map('usa')
colorvec = colorRampPalette(mycol)(100)[cut(yW_res, 100, include.lowest=T)]
plot(het_sp, col=colorvec, add=T, pch=16)


## Add regional climate to model to control for east-west
full_mod = glm.nb(richness ~ ., data=cbind(use_data_test[,c('richness','regS')], W))
full_mod = update(full_mod, ~.+regS*reghet_pc1)

# Plot lines
lowhet = quantile(reghet_pc1, 0.25)
highhet = quantile(reghet_pc1, 0.75)
use_coef = coef(full_mod)
RO_term = sum(use_coef[3:8]*colMeans(W))

use_x = seq(min(use_data_test$regS), max(use_data_test$regS), length.out=100)
use_ylow = exp(use_coef[1]+use_coef[2]*use_x+RO_term+use_coef[9]*lowhet+use_coef[10]*lowhet*use_x) 
use_yhigh = exp(use_coef[1]+use_coef[2]*use_x+RO_term+use_coef[9]*highhet+use_coef[10]*highhet*use_x) 

colorvec = colorRampPalette(mycol)(100)[cut(reghet_pc1, 100, include.lowest=T)]
plot(richness~regS, data=use_data_test,
	xlab='Regional Species Richness', ylab='Local Species Richness', 
	las=1, xlim=c(230, 430), ylim=c(0,40),
	col=colorvec, lwd=2, axes=F)
axis(1, at=seq(230,430,50))
axis(2)
box()
lines(use_x, use_ylow, lwd=5, col='white')
lines(use_x, use_yhigh, lwd=5, col=mycol[10])
lines(use_x, use_ylow, lwd=4, col=mycol[1], lty=2)
lines(use_x, use_yhigh, lwd=4, col='white', lty=3)

## Categorical East-West: divide by longitude -100
ew = as.numeric(master[rownames(use_data_test), 'LON'] > -100) # east = 1, west = 0

ew_mod = glm.nb(richness ~ regS*reghet_pc1*ew, data=use_data_test)
summary(ew_mod)


# 3-panel plot showing E vs W models and map of reghet
# Make separately and edit in Inkscape

# Map
pdf('./Figures/map_reghetpc1.pdf', height=4, width=7)
trellis.par.set(axis.line=list(col=NA))
spplot(het_laea, 'reghet_pc1', ylim=c(-1600,1500), main='', panel=function(x,y,subscripts,...){
	sp.polygons(OUTLINES.laea, fill='white')
	panel.pointsplot(x,y,...)
}, cuts=100, cex=1.5, col.regions = colorRampPalette(mycol)(100), auto.key=F)
dev.off()

use_x = seq(min(use_data_test$regS), max(use_data_test$regS), length.out=100)
use_coef = coef(ew_mod)

# 2-panel Models
svg('./Figures/reghet_interaction_models.svg', height=3.5, width=7)
text.cex = 1.2

par(cex.axis=text.cex)
par(cex.lab=text.cex)
par(pch=1)
par(cex=1)
par(las=1)
par(lend=1)
layout(matrix(c(1,2,3), nrow=1), widths=c(0.4, 0.4, 0.2))

colorvec = colorRampPalette(mycol)(100)[cut(reghet_pc1, 100, include.lowest=T)]
medhet = median(reghet_pc1)

par(mar=c(5,4,1.5,0))
plot(richness~regS, data=use_data_test[ew==0,],
	xlab='Regional Species Richness', ylab='Local Species Richness', 
	las=1, xlim=c(230, 430), ylim=c(0,40),
	col=colorvec[ew==0], lwd=2, axes=F, main='West')
axis(1, at=seq(230,430,50))
axis(2)
box()
lowhet = quantile(reghet_pc1[ew==0], 0.25)
highhet = quantile(reghet_pc1[ew==0], 0.75)
use_ylow = exp(cbind(1, use_x, lowhet, 0, use_x*lowhet, 0, 0, 0) %*% use_coef)
use_yhigh = exp(cbind(1, use_x, highhet, 0, use_x*highhet, 0, 0, 0) %*% use_coef) 
lines(use_x, use_ylow, lwd=5, col='white')
lines(use_x, use_yhigh, lwd=5, col=mycol[10])
lines(use_x, use_ylow, lwd=4, col=mycol[1], lty=2)
lines(use_x, use_yhigh, lwd=4, col='white', lty=3)

use_ymed = exp(cbind(1, use_x, medhet, 0, use_x*medhet, 0, 0, 0) %*% use_coef)
lines(use_x, use_ymed, lwd=5, col='black')


par(mar=c(5,3,1.5,1))
plot(richness~regS, data=use_data_test[ew==1,],
	xlab='Regional Species Richness', ylab='', 
	las=1, xlim=c(230, 430), ylim=c(0,40),
	col=colorvec[ew==1], lwd=2, axes=F, main='East')
axis(1, at=seq(230,430,50))
axis(2)
box()
lowhet = quantile(reghet_pc1[ew==1], 0.25)
highhet = quantile(reghet_pc1[ew==1], 0.75)
use_ylow = exp(cbind(1, use_x, lowhet, 1, use_x*lowhet, use_x, lowhet, lowhet*use_x) %*% use_coef)
use_yhigh = exp(cbind(1, use_x, highhet, 1, use_x*highhet, use_x, highhet, highhet*use_x) %*% use_coef) 
lines(use_x, use_ylow, lwd=5, col='white')
lines(use_x, use_yhigh, lwd=5, col=mycol[10])
lines(use_x, use_ylow, lwd=4, col=mycol[1], lty=2)
lines(use_x, use_yhigh, lwd=4, col='white', lty=3)

use_ymed = exp(cbind(1, use_x, medhet, 1, use_x*medhet, use_x, medhet, medhet*use_x) %*% use_coef)
lines(use_x, use_ymed, lwd=5, col='black')

par(mar=c(5,0,1.5,1))
plot.new()
usr = par('usr')
plotColorRamp(cols = mycol, n = 100, barends = c(usr[1], usr[3], usr[1]+0.2*diff(usr[1:2]), usr[4]),
	labels = seq(-3,3,1), uneven.lab=T, labrange=range(reghet_pc1), title='Regional Heterogeneity (PC1)',
	mycex=text.cex, ndig=1)

dev.off()

## 3D plot
library(rgl)
xpoints = seq(min(use_data_test$regS), max(use_data_test$regS), length.out=20)
ypoints = seq(min(reghet_pc1), max(reghet_pc1), length.out=20)


# Create a grid of regS vs reghet
XY = expand.grid(X=xpoints, Y=ypoints)

# Function from fitted model
# r is 0 or 1 for east or west 
Zf = function(X,Y,r){
     exp(cbind(1, X, Y, r, X*Y, X*r, Y*r, X*Y*r) %*% coef(ew_mod)) 
}

# Populate a surface for east and west
west_pred = Zf(XY$X, XY$Y, 0)
east_pred = Zf(XY$X, XY$Y, 1)

# plot

Z = east_pred
zlim = range(Z)
zlen = zlim[2] - zlim[1] + 1

jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
colorzjet <- jet.colors(100)  # 100 separate color 


open3d()
par3d(scale=c(1,15,1))
rgl.surface(x=xpoints, y=matrix(Z,20), 
            coords=c(1,3,2),z=ypoints, 
            color=colorzjet[ findInterval(Z, seq(min(Z), max(Z), length=100))] )

axes3d()

Z = west_pred
zlim = range(Z)
zlen = zlim[2] - zlim[1] + 1

jet.colors = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
colorzjet <- jet.colors(100)  # 100 separate color 

open3d()
par3d(scale=c(5,100,.5))
rgl.surface(x=xpoints, y=matrix(Z,20), 
            coords=c(1,3,2),z=ypoints, 
            color=colorzjet[ findInterval(Z, seq(min(Z), max(Z), length=100))] )

axes3d()

rgl.snapshot("copyMatlabstyle.png")


## Fit two separate models - because my main question isn't 'do they differ', but 'what is the interaction'
w_mod = glm.nb(richness ~ regS*reghet_pc1[ew==0], data=use_data_test[ew==0,])
e_mod = glm.nb(richness ~ regS*reghet_pc1[ew==1], data=use_data_test[ew==1,])

summary(w_mod)
summary(e_mod)

# Plot
svg('./Figures/reghet_interaction_models_separate_regions.svg', height=3.5, width=7)
text.cex = 1.2

par(cex.axis=text.cex)
par(cex.lab=text.cex)
par(pch=1)
par(cex=1)
par(las=1)
par(lend=1)
layout(matrix(c(1,2,3), nrow=1), widths=c(0.4, 0.4, 0.2))

colorvec = colorRampPalette(mycol)(100)[cut(reghet_pc1, 100, include.lowest=T)]
medhet = median(reghet_pc1)

par(mar=c(5,4,1.5,0))
plot(richness~regS, data=use_data_test[ew==0,],
	xlab='Regional Species Richness', ylab='Local Species Richness', 
	las=1, xlim=c(230, 430), ylim=c(0,40),
	col=colorvec[ew==0], lwd=2, axes=F, main='West')
axis(1, at=seq(230,430,50))
axis(2)
box()

use_coef = coef(w_mod)
lowhet = quantile(reghet_pc1[ew==0], 0.25)
highhet = quantile(reghet_pc1[ew==0], 0.75)
use_ylow = exp(cbind(1, use_x, lowhet, use_x*lowhet)%*% use_coef)
use_yhigh = exp(cbind(1, use_x, highhet, use_x*highhet) %*% use_coef) 
lines(use_x, use_ylow, lwd=5, col='white')
lines(use_x, use_yhigh, lwd=5, col=mycol[10])
lines(use_x, use_ylow, lwd=4, col=mycol[1], lty=2)
lines(use_x, use_yhigh, lwd=4, col='white', lty=3)

par(mar=c(5,3,1.5,1))
plot(richness~regS, data=use_data_test[ew==1,],
	xlab='Regional Species Richness', ylab='', 
	las=1, xlim=c(230, 430), ylim=c(0,40),
	col=colorvec[ew==1], lwd=2, axes=F, main='East')
axis(1, at=seq(230,430,50))
axis(2)
box()
use_coef = coef(e_mod)
lowhet = quantile(reghet_pc1[ew==1], 0.25)
highhet = quantile(reghet_pc1[ew==1], 0.75)
use_ylow = exp(cbind(1, use_x, lowhet, use_x*lowhet)%*% use_coef)
use_yhigh = exp(cbind(1, use_x, highhet, use_x*highhet) %*% use_coef) 
lines(use_x, use_ylow, lwd=5, col='white')
lines(use_x, use_yhigh, lwd=5, col=mycol[10])
lines(use_x, use_ylow, lwd=4, col=mycol[1], lty=2)
lines(use_x, use_yhigh, lwd=4, col='white', lty=3)

par(mar=c(5,0,1.5,1))
plot.new()
usr = par('usr')
plotColorRamp(cols = mycol, n = 100, barends = c(usr[1], usr[3], usr[1]+0.2*diff(usr[1:2]), usr[4]),
	labels = seq(-3,3,1), uneven.lab=T, labrange=range(reghet_pc1), title='Regional Heterogeneity (PC1)',
	mycex=text.cex, ndig=1)

dev.off()



## How does regional heterogeneity affect local control of local richness

use_modstring = paste('richness', paste(paste(localvars, '*reghet_pc1', sep=''), collapse='+'), sep='~')

## Opposite interaction for regional heterogeneity on strength of local filters?
reghetloc_mod = glm.nb(use_modstring, data=use_data_test, link='log')
cbind(coef(reghetloc_mod)[c(2,4:21)], coef(reghetloc_mod)[22:40], coef(summary(reghetloc_mod))[22:40,4])


# Probably should only do this for local scale variables that have significant effects.
# Decide using sem?
#use_locvars = c('pseas','wetness','radiation','bark_moist_pct.rao.ba','PIE.ba.tree','wood_SG.ba','light.mean','PC1')

interactionmods = sapply(localvars, function(x){
	
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


sigint_mods = subset(interactionmods, line_interaction_P <0.05)

sigint_mods[,c('line_slope','line_interaction')]

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
