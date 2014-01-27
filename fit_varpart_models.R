## This script performs the variation partitioning analysis on lichen diversity for FIA plots.


source('./UNC/Projects/FIA Lichen/GitHub/FIA-Lichens/load_data.R')
source('./GitHub/FIA-Lichens/fia_lichen_analysis_functions.R')

save.image('varpart_analysis.Rdata')
load('varpart_analysis.Rdata')

###################################################################################
### Variation Partitioning among climate, forest, regional richness, and pollution

library(hier.part)
library(MuMIn)
library(MASS) # glm.nb
library(Vennerable) 


# Read in table of predictor variable types
predtypes = read.csv('./SEM models/var_types.csv', row.names=1)
predtypes = subset(predtypes, !(rownames(predtypes) %in% c('FM','FH')))

# Define dataset- use trans_data where variables have been log or sqrt transformed but not scaled
use_data = trans_data[,colnames(trans_data)%in% rownames(predtypes)]

# Add other variables used to analyze FD and restricted taxonomic groups
other_data = master[rownames(use_data),c('lichen.rich','fric','raoQ','Parmeliaceae','Physciaceae',
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
use_vars = rownames(predtypes)[predtypes$type!='A'] # Leave out abundance
use_vars = use_vars[-grep('soil',use_vars)] # Leave out soil variables

pois_log_mod = glm(richness~., family=poisson(link='log'), data=use_data_fit[,c('richness',use_vars)])
pois_iden_mod = glm(richness~., family=poisson(link='identity'), start=rep(1,23), data=use_data_fit[,c('richness',use_vars)]) #warnings- may not have fit correctly
nb_log_mod = glm.nb(richness~., link='log', data=use_data_fit[,c('richness',use_vars)])
nb_iden_mod = glm.nb(richness~., link='identity', init.theta=10, data=use_data_fit[,c('richness',use_vars)])
gaus_log_mod = glm(richness~., family=gaussian(link='log'), data=use_data_fit[,c('richness',use_vars)])
gaus_iden_mod = glm(richness~., family=gaussian(link='identity'), data=use_data_fit[,c('richness',use_vars)])

AIC(gaus_iden_mod, gaus_log_mod, pois_iden_mod, pois_log_mod, nb_iden_mod, nb_log_mod)
# nb_log mod wins by far for lichen richness
# gaus_iden wins by far for fric, can't use nb or pois b/c not integer.

### Test linear, log-linear, Poisson, and Negative Binomial GLMs for abundance
abun_pois_log = glm(richness~abun_log, family=poisson(link='log'), data=use_data_fit)
#abun_pois_iden = glm(richness~abun_log, family=poisson(link='identity'), start=rep(1,2), data=use_data_fit) 
abun_nb_log = glm.nb(richness~abun_log, link='log', data=use_data_fit) # Doesn't fit b/c data is basically Poisson distributed
#abun_nb_iden = glm.nb(richness~abun_log, link='identity', init.theta=1, data=use_data_fit)
abun_gaus_log = glm(richness~abun_log, family=gaussian(link='log'), data=use_data_fit)
abun_gaus_iden = glm(richness~abun_log, family=gaussian(link='identity'), data=use_data_fit)

AIC(abun_pois_log, abun_gaus_log, abun_gaus_iden) 

pois_y = predict(abun_pois_log, list(abun_log=use_x), type='response')
gauslog_y = predict(abun_gaus_log, list(abun_log=use_x), type='response')
gausiden_y = predict(abun_gaus_iden, list(abun_log=use_x), type='response')

plot(richness~abun_log, data=use_data_test)
lines(use_x, pois_y, col='red')
lines(use_x, gauslog_y, col='blue')
lines(use_x, gausiden_y, col='green')

plot(log10(richness)~abun_log, data=use_data_test, ylim=c(-1,2))
lines(use_x, log10(pois_y), col='red', lwd=2)
lines(use_x, log10(gauslog_y), col='blue', lwd=2)
lines(use_x, log10(gausiden_y), col='green', lwd=2)

## Test gaussian, Poisson, NB for regS
pois_log = glm(richness~reg, family=poisson(link='log'), data=use_data_fit)
pois_iden = glm(richness~reg, family=poisson(link='identity'), data=use_data_fit)
nb_log = glm.nb(richness~reg, link='log', data=use_data_fit)
nb_iden = glm.nb(richness~reg, link='identity', data=use_data_fit)
gaus_log = glm(richness~reg, family=gaussian(link='log'), data=use_data_fit)
gaus_iden = glm(richness~reg, family=gaussian(link='identity'), data=use_data_fit)

summary(gaus_iden)

AIC(pois_log, pois_iden, nb_log, nb_iden, gaus_log, gaus_iden) # Ended up using nb_log for consistency with other variables

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
unimods = sapply(rownames(predtypes), function(x){
	
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

write.csv(unimods, 'univariate_models.csv', row.names=T)

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

# Make a chart of linear vs quadratic and concavity

modcompare = data.frame(deltaAIC = unimods$AIC_line-unimods$AIC_quad) # deltaAIC neg indicates AIC_line < AIC_quad and quadratic term not needed
rownames(modcompare) = rownames(unimods)
modcompare$type = ifelse(modcompare$deltaAIC>2, 'quadratic','linear')
modcompare$coef_sq = unimods$quad_sq
modcompare$concavity = ifelse(unimods$quad_sq>0, 'up','down')
modcompare$varType = predtypes[rownames(modcompare),'type']

modcompare = data.frame(predictor=varnames[rownames(modcompare),'midName'], modcompare)

modcompare = modcompare[order(modcompare$varType, modcompare$deltaAIC),]

# Note: even though coefficients are on log scale, the quadratic models will be 
# concave down whenever the quadratic coefficient is negative- proved using calculus.

write.csv(modcompare, 'Univariate model shapes GLM-NB.csv', row.names=T)
write.csv(modcompare, 'Univariate model shapes fric.csv', row.names=T)

modcompare = read.csv('Univariate model shapes GLM-NB.csv', row.names=1)
modcompare = read.csv('Univariate model shapes fric.csv', row.names=1)

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

# Using quadratic relationships for concave-down AIC supported models in variation partitioning
sq_df = use_data[,sq_vars]^2
colnames(sq_df) = paste(colnames(sq_df),'2', sep='')
use_data = cbind(use_data, sq_df)

## Use test plots for model results:

## Variation Partitioning by Local / Regional

regionvars = rownames(predtypes)[predtypes$scale=='R']
localvars = rownames(predtypes)[predtypes$scale=='L']
localvars = localvars[-grep('soil', localvars)] # Leave out soil vars
localvars = subset(localvars, localvars != 'abun_log') # Leavr out abundance

region_mod = glm.nb(richness~., data=use_data_test[,c('richness', regionvars, colnames(sq_df)[sq_vars %in% regionvars])])
local_mod = glm.nb(richness~., data=use_data_test[,c('richness', localvars, colnames(sq_df)[sq_vars %in% localvars])])

AIC(forest_mod, reg_mod, env_mod, pol_mod, region_mod, local_mod) # region_mod is best

# Make a list of predictors
predlist = list(Regional=names(region_mod$coefficients[-1]),
	Local=names(local_mod$coefficients[-1])
)

# Define null model
null_mod = glm.nb(richness~1, data=use_data_test)

# Calculate pseudo R2 for each model containing successive subsets of variables
apply(combos(2)$ragged, 1, function(x){
	use_vars = unlist(predlist[x])
	this_mod = glm.nb(richness~., data=use_data_test[,c('richness',use_vars)])
	r2 = r.squaredLR(this_mod, null=null_mod)
	attr(r2, 'adj.r.squared')
})->Rs

names(Rs) = names(predlist)

local_regional_partition = partvar2(Rs)

## Variation Partitioning by Forest Heterogeneity / Mean Conditions
nichevars = rownames(predtypes)[predtypes$type=='FH']
resvars = rownames(predtypes)[predtypes$type %in% c('FM','L')]
resvars = resvars[-grep('soil',resvars)] # Don't include soil

niche_mod = glm.nb(richness~., data=use_data_test[,c('richness', nichevars, colnames(sq_df)[sq_vars %in% nichevars])])
res_mod = glm.nb(richness~., data=use_data_test[,c('richness', resvars, colnames(sq_df)[sq_vars %in% resvars])])

AIC(forest_mod, reg_mod, env_mod, pol_mod, region_mod, local_mod, niche_mod, res_mod) # region_mod is best

# Make a list of predictors
predlist = list(Heterogeneity=names(niche_mod$coefficients[-1]),
	Optimality=names(res_mod$coefficients[-1])
)

# Define null model
null_mod = glm.nb(richness~1, data=use_data_test)

# Calculate pseudo R2 for each model containing successive subsets of variables
apply(combos(2)$ragged, 1, function(x){
	use_vars = unlist(predlist[x])
	this_mod = glm.nb(richness~., data=use_data_test[,c('richness',use_vars)])
	r2 = r.squaredLR(this_mod, null=null_mod)
	attr(r2, 'adj.r.squared')
})->Rs

names(Rs) = names(predlist)

niche_res_partition = partvar2(Rs)

## Variation partitioning for Fric
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


# Plot variation vartitioning
barwide = .5
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


# Plot without unexplained variation
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


# Plot variation vartitioning in black and white
barwide = .5
use_part = local_regional_partition

svg('./Figures/partition variance among local vs regional include quadratic effects bw.svg', height=6, width=6*barwide)
	par(mar=c(0,4,0,0))

	# Create plotting window
	plot(1,1, xlim=c(0,barwide), ylim=c(0,1), axes=F, xlab='', ylab='', type='n')
	
	# Add rectangle for unexplained variation
	rect(0,sum(use_part[1:3]),barwide,1, lwd=3, col='white')
	
	# Add rectangle for regional model
	rect(0,0,barwide, sum(use_part[1:2]), lwd=3, density=10, angle=45)

	# Add rectangle for local model
	rect(0,use_part[1],barwide,sum(use_part[1:3]), lwd=3, density=10, angle=-45)

	# Add axis
	axis(2, las=1, cex.axis=2, lwd=3)
	
	# Add partition labels
	lablocs = sapply(1:4, function(x) sum(use_part[0:(x-1)])+use_part[x]/2)
	text(barwide/2, lablocs, labels=names(use_part), cex=2, bg='white')
dev.off()


# Plot without unexplained variation
barwide=.6
use_part = niche_res_partition
svg('./Figures/partition variance among heterogeneity vs optimality include quadratic effects no unexp bw.svg', height=6, width=6*barwide)
	par(mar=c(0,0,0,5))

	# Create plotting window
	plot(1,1, xlim=c(0,barwide), ylim=c(0,sum(use_part[1:3])), axes=F, xlab='', ylab='', type='n')
	
	# Add rectangle for regional model
	rect(0,0,barwide, sum(use_part[1:2]), lwd=3, density=10, angle=20)

	# Add rectangle for local model
	rect(0,use_part[1],barwide,sum(use_part[1:3]), lwd=3, density=10, angle=-20)

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
