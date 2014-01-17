## This script performs the variation partitioning analysis on lichen diversity for FIA plots.

source('./UNC/Projects/FIA Lichen/GitHub/FIA-Lichens/load_data.R')
source('./GitHub/FIA-Lichens/fia_lichen_analysis_functions.R')

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

# Define predictor variable and appropriate abundance and regional richness variables
use_response = 'lichen.rich'
use_reg = 'regS'
use_abun = 'tot_abun_log'
use_data$richness = trans_data[,use_response]
use_data$reg = trans_data[,use_reg]
use_data$abun_log = trans_data[,use_abun]

# Append soil PCA variables
soil = read.csv('./Soil/soil_PCA.csv', row.names=1)
use_data$soilPC1 = soil[rownames(use_data),'PC1']
use_data$soilPC2 = soil[rownames(use_data), 'PC2']

### Test linear, log-linear, Poisson, and Negative Binomial GLMs with different link functions
use_vars = rownames(predtypes)[predtypes$type!='A'] # Leave out abundance
use_vars = use_vars[-grep('soil',use_vars)] # Leav out soil variables

pois_log_mod = glm(richness~., family=poisson(link='log'), data=use_data_fit[,c('richness',use_vars)])
pois_iden_mod = glm(richness~., family=poisson(link='identity'), start=rep(1,23), data=use_data_fit[,c('richness',use_vars)]) #warnings- may not have fit correctly
nb_log_mod = glm.nb(richness~., link='log', data=use_data_fit[,c('richness',use_vars)])
nb_iden_mod = glm.nb(richness~., link='identity', init.theta=10, data=use_data_fit[,c('richness',use_vars)])
gaus_log_mod = glm(richness~., family=gaussian(link='log'), data=use_data_fit[,c('richness',use_vars)])
gaus_iden_mod= glm(richness~., family=gaussian(link='identity'), data=use_data_fit[,c('richness',use_vars)])

AIC(gaus_iden_mod, gaus_log_mod, pois_iden_mod, pois_log_mod, nb_iden_mod, nb_log_mod)
# nb_log mod wins by far.


### Univariate models ###

# Initially use fitplots to determine whether to use quadratic relationships.
# Then, re-run final analysis using only use testing data set
use_data_fit = use_data[fitplots$yrplot.id,]
sum(is.na(use_data_fit)) # Checking for NAs- in soil variables

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

# Make a chart of linear vs quadratic and concavity

modcompare = data.frame(deltaAIC = unimods$AIC_line-unimods$AIC_quad) # deltaAIC neg indicates AIC_line < AIC_quad and quadratic term not needed
rownames(modcompare) = rownames(unimods)
modcompare$type = ifelse(modcompare$deltaAIC>2, 'quadratic','linear')
modcompare$coef_sq = unimods$quad_sq
modcompare$concavity = ifelse(unimods$quad_sq>0, 'up','down')
# Note: even though coefficients are on log scale, the quadratic models will be 
# concave down whenever the quadratic coefficient is negative- proved using calculus.

write.csv(modcompare, 'Univariate model shapes GLM-NB.csv', row.names=T)

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
sq_df = use_data_[,sq_vars]^2
colnames(sq_df) = paste(colnames(sq_df),'2', sep='')
use_data = cbind(use_data, sq_df)





















forestvars = predtypes$predictor[predtypes$type %in% c('FM','FH')]
climvars = predtypes$predictor[predtypes$type == 'C']
regvars = predtypes$predictor[predtypes$type == 'R']
polvars = predtypes$predictor[predtypes$type == 'P']

# Using only linear relationships
forest_mod = glm.nb(richness~., data=use_data[,c('richness', forestvars)] )
reg_mod = glm.nb(richness~., data=use_data[,c('richness', regvars)])
clim_mod = glm.nb(richness~., data=use_data[,c('richness', climvars)])
pol_mod = glm.nb(richness~., data=use_data[,c('richness', polvars)])









forest_mod = glm.nb(richness~. , data=use_data[,c('richness', forestvars, colnames(sq_df)[sq_vars %in% forestvars])])
reg_mod = glm.nb(richness~., data=use_data[,c('richness', regvars, colnames(sq_df)[sq_vars %in% regvars])])
clim_mod = glm.nb(richness~., data=use_data[,c('richness', climvars, colnames(sq_df)[sq_vars %in% climvars])])
pol_mod = glm.nb(richness~., data=use_data[,c('richness', polvars, colnames(sq_df)[sq_vars %in% polvars])])

# Make a list of predictors
predlist = list(climate=names(clim_mod$coefficients[-1]),
	forest=names(forest_mod$coefficients[-1]),
	regional=names(reg_mod$coefficients[-1]),
	pollution=names(pol_mod$coefficients[-1])
)

# Calculate pseudo R2 for each model containing successive subsets of variables
apply(combos(4)$ragged, 1, function(x){
	use_vars = unlist(predlist[x])

	this_mod = glm.nb(richness~., data=use_data[,c('richness',use_vars)])
	r.squaredLR(this_mod)
})->Rs

# Add the null model
nullR = r.squaredLR(glm.nb(richness~1, data=use_data))
Rs = c(nullR,Rs)

# Partition variation among variable sets (from hier.part package)
partition(Rs, 4, var.names = names(predlist))-> var_sets_partition 

pdf('./Figures/partition variance among predictor sets include quadratic effects barplot.pdf', heigh=6, width=7)
par(mar=c(4,6,1,1))
barplot(t(as.matrix(var_sets_partition$IJ[,c('I','J')])),
	legend.text = c('Independent','Joint'),ylim=c(0,.4),
	las=1, cex.axis=1.5, cex.lab=1.5, cex.names=1.5,col=c('grey30','grey60'),
	args.legend=list(bty='n', cex=1.5, border='transparent'), mgp=c(3,2,1),
	names=c('Climate','Forest\nStructure','Regional\nRichness','Pollution'), border=F
)
mtext('Variation Explained',2,4.5,cex=1.5)

dev.off()

### Need to do this for models with quadratic effects


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
apply(combos(3)$ragged, 1, function(x){
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
