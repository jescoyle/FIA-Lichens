## This script conducts variation partitioning analyses and variable importance analyses
## use in the second version of the FIA Lichen manuscript


source('./UNC/Projects/FIA Lichen/GitHub/FIA-Lichens/load_data.R')
source('./GitHub/FIA-Lichens/fia_lichen_analysis_functions.R')

#save.image('varpart_analysis.Rdata')
#load('varpart_analysis.Rdata')

# Color ramp
mycol = read.csv('C:/Users/jrcoyle/Documents/UNC/Projects/blue2red_10colramp.txt')
mycol = apply(mycol,1,function(x) rgb(x[1],x[2],x[3],maxColorValue=256))
mycol = mycol[10:1]

###################################################################################
### Read in and organize data

library(MuMIn) # r.squaredLR
library(MASS) # glm.nb
library(spdep) # spatial autocorrelation
library(sp) # spatial data handling

# Read in table of predictor variable types
predtypes = read.csv('predictors.csv', row.names=1)

# Define dataset- use trans_data where variables have been log or sqrt transformed but not scaled
use_data = trans_data[,colnames(trans_data) %in% rownames(predtypes)]

# Add other variables used to analyze FD and restricted taxonomic groups
other_data = trans_data[,c('lichen.rich','fric','raoQ','Parmeliaceae','Physciaceae',
	'tot_abun_log','parm_abun_log','phys_abun_log','regS','regFIA','regParm','regPhys')]

# Define predictor variable and appropriate abundance and regional richness variables
use_response = 'lichen.rich' # This changes based on the analysis
use_reg = 'regFIA' # Main analyses used to be with regS
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

# Define spatial dataframe for spatial analysis at regional scale
spdata = master[rownames(use_data),c('LAT','LON','lichen.rich','regS')]
coordinates(spdata) = c('LON','LAT'); proj4string(spdata) = CRS("+proj=longlat")
spdata_test = spdata[testplots$yrplot.id,]
spdata_fit = spdata[fitplots$yrplot.id,]

##################################################################################
### Univariate models used to decide which variables will have a squared term in subsequent models

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


# Spatial error models
# Models use data scaled to mean=0, var=1 because otherwise models don't fit
use_vars = rownames(subset(predtypes, scale=='regional'&!(label%in%c('','R1'))))
reg_nb = make_regnb(spdata_fit)
reg_listw = make_reglistw(spdata_fit, reg_nb)

sarmods_reg = sapply(use_vars, function(x){
	
	# By default uses log link function, so coefficients are on log scale
	linear_mod = errorsarlm(working_data_fit[,use_reg] ~ working_data_fit[,x], listw=reg_listw)
	quad_mod = 	errorsarlm(working_data_fit[,use_reg]~working_data_fit[,x]+I(working_data_fit[,x]^2), listw = reg_listw)

	c(AIC(linear_mod), AIC(quad_mod), linear_mod$s2, quad_mod$s2, coef(linear_mod),	coef(quad_mod))
})

sarmods_reg = data.frame(t(sarmods_reg))
names(sarmods_reg) = c('AIC_line','AIC_quad','ML_s2_line','ML_s2_quad',
	'line_lambda','line_int','line_slope','quad_lambda', 'quad_int','quad_slope','quad_sq')

write.csv(sarmods_reg, 'univariate_models_regFIA_AllSp_sar.csv', row.names=T)

## Make a chart of linear vs quadratic and concavity

modcompare = data.frame(deltaAIC = unimods$AIC_line-unimods$AIC_quad) # deltaAIC neg indicates AIC_line < AIC_quad and quadratic term not needed
rownames(modcompare) = rownames(unimods)
modcompare$type = ifelse(modcompare$deltaAIC>2, 'quadratic','linear')
modcompare$coef_sq = unimods$quad_sq
modcompare$concavity = ifelse(unimods$quad_sq>0, 'up','down')
modcompare = cbind(modcompare, predtypes[rownames(modcompare),c('type','scale','mode')])

modcompare = data.frame(predictor=varnames[rownames(modcompare),'midName'], modcompare)
modcompare = modcompare[order(modcompare$scale, modcompare$mode, modcompare$deltaAIC),]

modcompare_sar = data.frame(deltaAIC = sarmods_reg$AIC_line-sarmods_reg$AIC_quad) # deltaAIC neg indicates AIC_line < AIC_quad and quadratic term not needed
rownames(modcompare_sar) = rownames(sarmods_reg)
modcompare_sar$type = ifelse(modcompare_sar$deltaAIC>2, 'quadratic','linear')
modcompare_sar$coef_sq = sarmods_reg$quad_sq
modcompare_sar$concavity = ifelse(sarmods_reg$quad_sq>0, 'up','down')
modcompare_sar = cbind(modcompare_sar, predtypes[rownames(modcompare_sar),c('type','scale','mode')])

modcompare_sar = data.frame(predictor=varnames[rownames(modcompare_sar),'midName'], modcompare_sar)
modcompare_sar = modcompare_sar[order(modcompare_sar$scale, modcompare_sar$mode, modcompare_sar$deltaAIC),]

write.csv(modcompare, 'Univariate model shapes GLM-NB.csv', row.names=T)
write.csv(modcompare_sar, 'Univariate model shapes of regFIA SAR.csv', row.names=T)

# Read in previously saved tables
modcompare = read.csv('Univariate model shapes GLM-NB.csv', row.names=1)
modcompare_sar = read.csv('Univariate model shapes of regFIA SAR.csv', row.names=1)

# Which variables have AIC supported concave-down relationships?
sq_vars = rownames(subset(modcompare, concavity=='down'&type=='quadratic'))
sq_vars_sar = rownames(subset(modcompare_sar, concavity=='down'&type=='quadratic'))


##############################################################################
### Models for Variation Partitioning

# Define sets of predictors to be used in models
climvars = c('mat','iso','pseas','wetness','rain_lowRH')
LOvars = c(climvars, 'radiation','bark_moist_pct.ba','wood_SG.ba','light.mean', 'PC1', 'bigTrees')
LHvars = c('bark_moist_pct.rao.ba','wood_SG.rao.ba','lightDist.mean','PIE.ba.tree','propDead', 'diamDiversity')
ROvars = paste(climvars, 'reg_mean', sep='_')
RHvars = c(paste(climvars, 'reg_var', sep='_'), 'regS_tree')
Rvars = c(RHvars, ROvars)
Lvars = c(LHvars, LOvars)
Hvars = c(RHvars, LHvars)
Ovars = c(ROvars, LOvars)
allvars = c(Rvars, Lvars)

### Local vs. regional variance partitioning of Local richness
# Define data for including squared terms
sqdata = use_data_test[,sq_vars]
sqdata = sqdata^2
colnames(sqdata) = paste(colnames(sqdata), 2, sep='')

full_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness','reg',allvars)], sqdata[,paste(sq_vars[sq_vars %in% allvars],2, sep='')]))

RregS_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness','reg',Rvars)], sqdata[,paste(sq_vars[sq_vars %in% Rvars],2, sep='')]))
R_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness',Rvars)], sqdata[,paste(sq_vars[sq_vars %in% Rvars],2, sep='')]))
regS_mod = glm.nb(richness~reg, data=use_data_test)
L_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness', Lvars)], sqdata[,paste(sq_vars[sq_vars %in% Lvars],2,sep='')]))

RH_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness',RHvars)], sqdata[,paste(sq_vars[sq_vars %in% RHvars],2, sep='')]), link='log')
RO_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness',ROvars)], sqdata[,paste(sq_vars[sq_vars %in% ROvars],2, sep='')]), link='log')

LH_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness',LHvars)], sqdata[,paste(sq_vars[sq_vars %in% LHvars],2, sep='')]), link='log')
LO_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness',LOvars)], sqdata[,paste(sq_vars[sq_vars %in% LOvars],2, sep='')]), link='log')

H_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness',Hvars)], sqdata[,paste(sq_vars[sq_vars %in% Hvars],2, sep='')]), link='log')
O_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness',Ovars)], sqdata[,paste(sq_vars[sq_vars %in% Ovars],2, sep='')]), link='log')

# Define null model
null_mod = glm.nb(richness~1, link='log', data=use_data_test)

# Define function for calculating pseduo-R2 from glm
calc_r2 = function(mod, null_mod){
	r2 = r.squaredLR(mod, null=null_mod)
	attr(r2, 'adj.r.squared')
}

# All variables
Rs = sapply(list(RregS_mod, L_mod, full_mod), function(x) calc_r2(x, null_mod)) 
names(Rs) = c('region','local','full')
part_regloc = partvar2(Rs)

# Within heterogeneity variables
Rs = sapply(list(RH_mod, LH_mod, H_mod), function(x) calc_r2(x, null_mod)) 
names(Rs) = c('region','local','full')
part_regloc_het = partvar2(Rs)

# Within optimality variables
Rs = sapply(list(RO_mod, LO_mod, O_mod), function(x) calc_r2(x, null_mod)) 
names(Rs) = c('region','local','full')
part_regloc_opt = partvar2(Rs)

# Print out model stats
modlist = list(RregS_mod, L_mod, full_mod, RO_mod, LO_mod, O_mod, RH_mod, LH_mod, H_mod)
names(modlist) = paste('L', 1:9, sep='')
sapply(modlist, function(x) {
	R2 = attr(r.squaredLR(x, null=null_mod), 'adj.r.squared')
	Deviance = deviance(x)
	aic = AIC(x)
	data.frame(R2, Deviance, aic)
})



### Heterogeneity vs. optimality variance partitioning

## Within local variables
Rs = sapply(list(LH_mod, LO_mod, L_mod), function(x) calc_r2(x, null_mod)) 
names(Rs) = c('het','opt','full')
part_hetopt_loc = partvar2(Rs)

## Within regional variables
# Define data for including squared terms in regional-scale spatial models
sqdata_reg = working_data_test[,sq_vars_sar]
sqdata_reg = sqdata_reg^2
colnames(sqdata_reg) = paste(colnames(sqdata_reg), 2, sep='')

# Make spatial objects for test data set
reg_nb = make_regnb(spdata_test)
reg_listw = make_reglistw(spdata_test, reg_nb)

# Fit spatial error models for regional richness
RH_regmod = errorsarlm(working_data_test[,use_reg] ~ ., data = cbind(working_data_test[,c(RHvars)], sqdata_reg[,paste(sq_vars_sar[sq_vars_sar %in% RHvars],2,sep='')]), listw=reg_listw)
RO_regmod = errorsarlm(working_data_test[,use_reg] ~ ., data = cbind(working_data_test[,c(ROvars)], sqdata_reg[,paste(sq_vars_sar[sq_vars_sar %in% ROvars],2,sep='')]), listw=reg_listw)
R_regmod = errorsarlm(working_data_test[,use_reg] ~ ., data = cbind(working_data_test[,c(RHvars, ROvars)], sqdata_reg[,paste(sq_vars_sar[sq_vars_sar %in% c(ROvars,RHvars)],2,sep='')]), listw=reg_listw)

# Calculate R2
null_mod = errorsarlm(working_data_test[,use_reg] ~ 1, data=working_data_test, listw=reg_listw)
Rs = sapply(list(RH_regmod, RO_regmod, R_regmod), function(x) calc_r2(x, null_mod))
names(Rs) = c('het','opt','full')
part_hetopt_reg = partvar2(Rs)

# Print out K for all models
L_K = lapply(list(RregS_mod, L_mod, full_mod, RO_mod, LO_mod, O_mod, RH_mod, LH_mod, H_mod), logLik)
R_K = lapply(list(RH_regmod, RO_regmod, R_regmod), logLik)

# Print out stats for all models
modlist = list(RH_regmod, RO_regmod, R_regmod)
names(modlist) = paste('R', 1:3, sep='')
sapply(modlist, function(x) {
	R2 = attr(r.squaredLR(x, null=null_mod), 'adj.r.squared')
	Deviance = deviance(x)
	aic = AIC(x)
	data.frame(R2, Deviance, aic)
})

## Figure. Plot local-regional variation partioning analysis in three panels
barwide=1.5
use_shade = c('CC','55','99','FF')
use_color = c('#2415B0','#00BF32','#126A71')

svg('./Figures/New Analysis/variation partitioning loc-reg FIA all obs.svg', height=5, width=10)
	par(mar=c(0,6,1.5,0))

	# Create plotting window
	plot(1,1, xlim=c(0,8), ylim=c(0,1), axes=F, ylab='', xlab='', type='n', cex.lab=2)
	
	## Add background lines
	usr = par('usr')
	segments(-.5, seq(0,1,.1), 5.25, seq(0,1,.1), lwd=2, col='grey70', lty=3)	

	## Background bars
	rect(0,0,barwide,1,col='white', lwd=3)	
	rect(2,0,2+barwide,1,col='white', lwd=3)
	rect(4,0,4+barwide,1,col='white', lwd=3)

	## Full Model
	# Add rectangle for first component
	rect(0,0,barwide, sum(part_regloc[1]), lwd=3, col=paste(use_color[1],use_shade[1],sep=''))
	# Add rectangle for second component
	rect(0,sum(part_regloc[1:2]), barwide,sum(part_regloc[1:3]), lwd=3, col=paste(use_color[2],use_shade[1],sep=''))
	# Add rectangle for overlap
	rect(0,part_regloc[1],barwide,sum(part_regloc[1:2]), lwd=3, col=paste(use_color[3],use_shade[1],sep=''))

	## Among Optimality variables
	# Add rectangle for first component
	rect(2,0,2+barwide, sum(part_regloc_opt[1]), lwd=3, col=paste(use_color[1],use_shade[1],sep=''))
	# Add rectangle for second component
	rect(2,sum(part_regloc_opt[1:2]), 2+barwide,sum(part_regloc_opt[1:3]), lwd=3, col=paste(use_color[2],use_shade[1],sep=''))
	# Add rectangle for overlap
	rect(2,part_regloc_opt[1], 2+barwide,sum(part_regloc_opt[1:2]), lwd=3, col=paste(use_color[3],use_shade[1],sep=''))

	## Among Heterogeneity variables
	# Add rectangle for first component
	rect(4,0,4+barwide, sum(part_regloc_het[1]), lwd=3, col=paste(use_color[1],use_shade[1],sep=''))
	# Add rectangle for second component
	rect(4,sum(part_regloc_het[1:2]), 4+barwide,sum(part_regloc_het[1:3]), lwd=3, col=paste(use_color[2],use_shade[1],sep=''))
	# Add rectangle for overlap
	rect(4,part_regloc_het[1], 4+barwide,sum(part_regloc_het[1:2]), lwd=3, col=paste(use_color[3],use_shade[1],sep=''))

	## Add axis
	axis(2, las=1, cex.axis=2, lwd=3)
	mtext('Variation Explained', 2, 4, cex=2)
	
	# Add partition labels
	#lablocs = sapply(1:4, function(x) sum(part_regloc[0:(x-1)])+part_regloc[x]/2)
	#text(barwide/2, lablocs, labels=names(part_regloc)[1:4], cex=2)
	#lablocs = sapply(1:4, function(x) sum(part_regloc_opt[0:(x-1)])+regional_het_opt_partition[x]/2)
	#text(2+barwide/2, lablocs, labels=names(regional_het_opt_partition)[1:4], cex=2)
	#lablocs = sapply(1:4, function(x) sum(local_het_opt_partition[0:(x-1)])+local_het_opt_partition[x]/2)
	#text(4+barwide/2, lablocs, labels=names(local_het_opt_partition)[1:4], cex=2)
	
	# Add Key
	legend('right', c('Unexplained', 'Local','Shared','Regional'), 
		fill= c('white',paste(use_color[c(2,3,1)], use_shade[1], sep='')),
		cex=2, bty='n')

	# Add panel text A, B, C
	#par(xpd=T)
	#text(c(0,2,4),1.05, c('A','B','C'), cex=2, adj=c(0,0))
	#par(xpd=F)
dev.off()


## Figure. Plot Heterogeneity - Optimality variation partitioning in 2 panels
svg('./Figures/New Analysis/variation partitioning het-opt FIA all obs.svg', height=5, width=8)
	par(mar=c(0,6,1.5,0))

	# Create plotting window
	plot(1,1, xlim=c(0,5.5), ylim=c(0,1), axes=F, ylab='', xlab='', type='n', cex.lab=2)
	
	## Add background lines
	usr = par('usr')
	segments(-.5, seq(0,1,.1), 5.25, seq(0,1,.1), lwd=2, col='grey70', lty=3)	

	## Background bars
	rect(0,0,barwide,1,col='white', lwd=3)	
	rect(2,0,2+barwide,1,col='white', lwd=3)

	## Regional Model
	# Add rectangle for first component
	rect(0,0,barwide, sum(part_hetopt_reg[1]), lwd=3, col=paste(use_color[1],use_shade[2],sep=''))
	# Add rectangle for second component
	rect(0,sum(part_hetopt_reg[1:2]), barwide,sum(part_hetopt_reg[1:3]), lwd=3, col=paste(use_color[1],use_shade[1],sep=''))
	# Add rectangle for overlap
	rect(0,part_hetopt_reg[1],barwide,sum(part_hetopt_reg[1:2]), lwd=3, col=paste(use_color[1],use_shade[3],sep=''))

	## Local Model
	# Add rectangle for first component
	rect(2,0,2+barwide, sum(part_hetopt_loc[1]), lwd=3, col=paste(use_color[2],use_shade[2],sep=''))
	# Add rectangle for second component
	rect(2,sum(part_hetopt_loc[1:2]), 2+barwide,sum(part_hetopt_loc[1:3]), lwd=3, col=paste(use_color[2],use_shade[1],sep=''))
	# Add rectangle for overlap
	rect(2,part_hetopt_loc[1], 2+barwide,sum(part_hetopt_loc[1:2]), lwd=3, col=paste(use_color[2],use_shade[3],sep=''))

	## Add axis
	axis(2, las=1, cex.axis=2, lwd=3)
	mtext('Variation Explained', 2, 4, cex=2)
	
	# Add Text
	#text_height = c(0,.1, .4, .8)
	#text(2+barwide/2, text_height, labels=c('Heterogeneity','Shared','Optimality','Unexplained'), pos=3, offset=0, cex=2)

	# Add panel text A, B, C
	#par(xpd=T)
	#text(c(0,4),1.05, c('A','B'), cex=2, adj=c(0,0))
	#par(xpd=F)
	
	# Add Key
	legend('right', c('Unexplained', 'Heterogeneity','Shared','Optimality'), 
		fill= c('white',paste(use_color[c(2,3,1)], use_shade[1], sep='')),
		cex=2, bty='n')
dev.off()

####################################################################################
### Model averaging and variable relative importance

# Most model averaging was done on the computing cluster using the scripts:
# pdredge_O_models.R
# pdredge_H_models.R
# pdredge_local_rich_models.R


# write out a data set for use on the cluster
#write.csv(use_data_test, './Data/localmods_data.csv', row.names=T)

options(na.action='na.fail') # So that models are not fit to different subsets of data

## Models of local richness

# Scalars use to standardized regression coefficients are based on standard deviations of variables in full data set
sds = apply(use_data[,c('richness', allvars)], 2, function(x) sqrt(var(x, na.rm=T)))
sds = sds/sds['richness']

## Get top variables in each category

# Optimality variables
#Odredge = dredge(O_mod, beta=T, subset=dc('light.mean','light.mean2')&dc('wood_SG.ba','wood_SG.ba2')&dc('radiation','radiation2')&dc('PC1','PC12')&dc('wetness','wetness2')dc('rain_lowRH_reg_mean','rain_lowRH_reg_mean2')&dc('wetness_reg_mean', 'wetness_reg_mean2'))
load('./Data/local_rich_opt_models_best.RData') # Note that there were 1119727 models so best models were subseted on the cluster
O_avg = model.avg(Odredge_best, beta=T)

svg('./Figures/New Analysis/local richness opt model avg coefs.svg', height=13, width=17)
plot_coefs(make_coefTable(O_avg), 'scale')
dev.off()

# Table for publication including variable importance and coefficents
O_coef = format_coefTable(O_avg)
write.table(O_coef, './Figures/New Analysis/local richness opt model avg coefs.txt', sep='\t', row.names=F, quote=F)

# Heterogeneity variables
#Hdredge = dredge(H_mod, beta=T, subset=dc('propDead','propDead2')&dc('wood_SG.rao.ba','wood_SG.rao.ba2')&dc('bark_moist_pct.rao.ba','bark_moist_pct.rao.ba2')&dc('pseas_reg_var', 'pseas_reg_var2')&dc('regS_tree', 'regS_tree2'))
load('./Data/local_rich_het_models_best.RData') # There were 31091 models so best models were subsetted on the cluster
#keep_models = which(Hdredge$weight/Hdredge$weight[1] >= 0.05) # Keep models whose evidence ratio is greater than 0.05

# These are effects of heterogeneity variables on local richness
H_avg = model.avg(Hdredge_best, beta=T)
#H_avg_full = model.avg(Hdredge, beta=T) # average across all models

svg('./Figures/New Analysis/local richness het model avg coefs.svg', height=13, width=17)
plot_coefs(make_coefTable(H_avg), 'scale')
dev.off()

# Table for publication including variable importance and coefficents
H_coef = format_coefTable(H_avg)
write.table(H_coef, './Figures/New Analysis/local richness het model avg coefs.txt', sep='\t', row.names=F, quote=F)

# Dredge full models on cluster using script 'dredge_local_rich_models.R'
#Rdredge = dredge(use_R_mod, beta=T, subset=dc('wetness_reg_mean', 'wetness_reg_mean2')&dc('pseas_reg_var', 'pseas_reg_var2')&dc('regS_tree', 'regS_tree2'), m.min=2)
#Ldredge = dredge(L_mod, beta=T, subset=dc('light.mean','light.mean2')&dc('wood_SG.ba','wood_SG.ba2')&dc('radiation','radiation2')&dc('PC1','PC12')&dc('wetness','wetness2')&
#	dc('propDead','propDead2')&dc('wood_SG.rao.ba','wood_SG.rao.ba2')&dc('bark_moist_pct.rao.ba','bark_moist_pct.rao.ba2'),
#	m.min=2)

load('./Data/local_rich_models_best.RData')

# Average best models (1545 out of 3359214 total)
L_avg = model.avg(Ldredge_best, beta=T)

# Make table of coefficients from model averages
L_coef = format_coefTable(L_avg)
write.table(L_coef, './Figures/New Analysis/local richness model avg coefs.txt', sep='\t', row.names=F, quote=F)

# Calculate model-averaged coefficients for models with 95% cumulative weight
#cum_weight = calc_cumWeight(Rdredge$weight)
#keep_models = 1:(max(which(cum_weight < 0.95))+1)

# These are effects of regional scale variables on local richness
# NOT USING THESE IN THE ANALYSIS
#R_avg = model.avg(Rdredge[keep_models,])
#R_avg_full = model.avg(Rdredge) # average across all models

## Plot comparison of local vs regional climate effects on local richness from model averaging with all optimality variables
keep_vars = c(climvars, ROvars)
O_coefs = make_coefTable(O_avg)
xvar = O_coefs[c(climvars,'wetness2'),]
yvar = O_coefs[c(ROvars,'wetness_reg_mean2'),]


lims = range(rbind(xvar, yvar)[,c('ci.upper','ci.lower')])
lims = c(-1, 1)*max(abs(lims))

pdf('./Figures/New Analysis/Compare climate effects in avg optimality model.pdf', height=4.5, width=4.5)
par(mar=c(4,4.5,1,1))
plot(xvar$est, yvar$est, type='n', xlim=lims, ylim=lims, las=1, xlab='Local variable effect',ylab='')
mtext('Regional variable effect', 2, 3.5)
usr =par('usr')
polygon(c(usr[1],0,usr[1],usr[2],0,usr[2],usr[1]),
	c(usr[3],0,usr[4],usr[4],0,usr[3],usr[3]), col='grey80')
abline(h=0, v=0, lty=2)
segments(xvar$est, yvar$ci.lower, xvar$est, yvar$ci.upper, lwd=2, col='black', lend=1)
segments(xvar$ci.lower, yvar$est, xvar$ci.upper, yvar$est, lwd=2, col='black', lend=1)
points(xvar$est, yvar$est, pch=15)
labels=expression(MAT,ISO,PSEAS,WET,LOWRH,WET^2)
text(xvar$est, yvar$est, labels, pos=c(4,4,2,4,4,2), offset=1)
dev.off()


## Models of regional richness

working_data_test = cbind(working_data_test, sqdata_reg)

# Can change response variable: regS, regFIA
RH_regmod = errorsarlm(regFIA ~ wetness_reg_var + iso_reg_var + rain_lowRH_reg_var + pseas_reg_var + mat_reg_var + regS_tree + regS_tree2, data = working_data_test, listw=reg_listw)
RO_regmod = errorsarlm(regFIA ~ wetness_reg_mean + iso_reg_mean + rain_lowRH_reg_mean + pseas_reg_mean + mat_reg_mean + pseas_reg_mean2, data = working_data_test, listw=reg_listw)
R_regmod = errorsarlm(regFIA ~ wetness_reg_mean + iso_reg_mean + rain_lowRH_reg_mean + pseas_reg_mean + mat_reg_mean + pseas_reg_mean2 + 
	wetness_reg_var + iso_reg_var + rain_lowRH_reg_var + pseas_reg_var + mat_reg_var + regS_tree + regS_tree2, data = working_data_test, listw=reg_listw)


save(R_regmod, RH_regmod, RO_regmod, reg_listw, working_data_test, file='regmods_objects_FIA.RData')

# Do these on cluster using scripts 'dredge_regional_models.R' and 'dredge_regional_model_full.R'
#ROreg_dredge = dredge(RO_regmod, subset=dc('pseas_reg_mean', 'pseas_reg_mean2'))
#RHreg_dredge = dredge(RH_regmod, subset=dc('regS_tree','regS_tree2'))
#R_dredge = dredge(R_regmod, subset=dc('pseas_reg_mean', 'pseas_reg_mean2')&dc('regS_tree','regS_tree2'))

load('./Data/regmod_full_dredged.RData')
load('./Data/regmods_dredged.RData')
load('./Data/regmod_full_dredged_FIA.RData')
load('./Data/regmods_dredged_FIA.RData')


subset(R_dredge, weight >0.01) 
subset(RHreg_dredge, weight >0.01) # Remove regS_tree and its square
subset(ROreg_dredge, weight >0.01) # Possibly remove wetness_reg_mean

# Calculate model-averaged coefficients for models whose evidence ratio is greater than 0.05
keep_models = which(R_dredge$weight/R_dredge$weight[1] >= 0.05) # 11 best models out of 4608 total , 14 for regS with FIA species

Rreg_avg = model.avg(R_dredge[keep_models,])
Rreg_avg_full = model.avg(R_dredge) # average across all models

# Note that regional richness models were fit to standardized variables so coefficients are standardized


## Plot coefficients from model averages
library(lattice)

svg('./Figures/New Analysis/regional richness model avg coefs.svg', height=13, width=17)
plot_coefs(make_coefTable(Rreg_avg), 'mode')
dev.off()

svg('./Figures/New Analysis/regional richness model avg full coefs.svg', height=13, width=17)
plot_coefs(make_coefTable(Rreg_avg_full), 'mode')
dev.off()


svg('./Figures/New Analysis/regional richness model avg full coefs.svg', height=13, width=17)
plot_coefs(make_coefTable(Rreg_avg_full), 'mode')
dev.off()

write.csv(make_impTable(Rreg_avg, length(keep_models)), './Figures/New Analysis/regional richness model var importance FIA.csv', row.names=F)

## Make table of coefficients from model averages
Rreg_coef = format_coefTable(Rreg_avg)
write.table(Rreg_coef, './Figures/New Analysis/regional richness model avg coefs FIA.txt', sep='\t', row.names=F, quote=F)


### Pairwise comparisons of variables
hetopt_pairs = data.frame(het = c('wetness_reg_mean','rain_lowRH_reg_mean','pseas_reg_mean','iso_reg_mean','mat_reg_mean','bark_moist_pct.ba','wood_SG.ba','light.mean','PC1'),
	opt = c('wetness_reg_var','rain_lowRH_reg_var','pseas_reg_var','iso_reg_var','mat_reg_var','bark_moist_pct.rao.ba','wood_SG.rao.ba','lightDist.mean','PIE.ba.tree'),
	scale = c(rep('regional',5), rep('local',4)) )

locreg_pairs = data.frame(loc = c(climvars, 'PIE.ba.tree'), reg=c(ROvars,'regS_tree'))















