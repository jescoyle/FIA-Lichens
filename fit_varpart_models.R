## This script performs the variation partitioning analysis on lichen diversity for FIA plots.

source('load_data.R')
source('fia_lichen_analysis_functions.R')

###################################################################################
### Variation Partitioning among climate, forest, regional richness, and pollution

rownames(model_data) = model_data$yrplot.id

library(hier.part)
library(MuMIn)
library(MASS) # glm.nb
library(Vennerable) 

# Define dataset
use_response = 'lichen.rich'
use_data = trans_data[,colnames(trans_data)!=use_response]

predtypes = allsp[,c('predictor','type')]
predtypes = subset(predtypes, !(predictor %in% c('FM','FH')))
use_data = use_data[, colnames(use_data) %in% predtypes$predictor]
use_data = data.frame(richness=model_data[rownames(trans_data),use_response], use_data)

# Initially use fitplots to determine whether to use quadratic relationships.
# Then, re-run final analysis using only use testing data set
use_data = use_data[fitplots,]
sum(is.na(use_data)) # Checking for NAs

forestvars = predtypes$predictor[predtypes$type %in% c('FM','FH')]
climvars = predtypes$predictor[predtypes$type == 'C']
regvars = predtypes$predictor[predtypes$type == 'R']
polvars = predtypes$predictor[predtypes$type == 'P']

# Using only linear relationships
forest_mod = glm.nb(richness~., data=use_data[,c('richness', forestvars)] )
reg_mod = glm.nb(richness~., data=use_data[,c('richness', regvars)])
clim_mod = glm.nb(richness~., data=use_data[,c('richness', climvars)])
pol_mod = glm.nb(richness~., data=use_data[,c('richness', polvars)])


## Univariate models ##
mymods = sapply(predtypes$predictor, function(x){
	

})




# Using quadratic relationships for concave-down AIC supported models
sq_vars = c('wetness','bark_moist_pct.rao.ba','wood_SG.rao.ba',
	'LogSeed.rao.ba','propDead','lightDist.mean','wood_SG.ba',
	'light.mean','LogSeed.ba','totalNS')
sq_df = use_data[,sq_vars]^2
colnames(sq_df) = paste(colnames(sq_df),'2', sep='')
use_data = cbind(use_data, sq_df)

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
