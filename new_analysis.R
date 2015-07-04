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
library(qpcR) # akaike.weights

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

## Functions that create neighbor, weights and listw objects given a spatial dataset

# Create a listw object with weights based on overlapping areas of 500km circles
# The weight is the percentage of the area shared with overlapping circle from another observation
A_intersect = function(d, R){ 2*(R^2)*acos(d/(2*R)) - 0.5*d*sqrt(4*(R^2) - (d^2))}
area_dist_func = function(d) A_intersect(d, 500)/(pi*(500^2))
make_regnb = function(sp_data) dnearneigh(sp_data, 0, 1000, row.names=rownames(coordinates(sp_data)))
make_regweights = function(sp_data, reg_nb) lapply(nbdists(reg_nb, sp_data), area_dist_func) 
make_reglistw = function(sp_data, reg_nb) nb2listw(reg_nb, glist=make_regweights(sp_data, reg_nb), style='W')


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
reg_nb = make_regnb(spdata_fit)
reg_listw = make_reglistw(spdata_fit, reg_nb)

sarmods_reg = sapply(use_vars, function(x){
	
	# By default uses log link function, so coefficients are on log scale
	linear_mod = errorsarlm(working_data_fit[, 'regS'] ~ working_data_fit[,x], listw=reg_listw)
	quad_mod = 	errorsarlm(working_data_fit[,'regS']~working_data_fit[,x]+I(working_data_fit[,x]^2), listw = reg_listw)

	c(AIC(linear_mod), AIC(quad_mod), linear_mod$s2, quad_mod$s2, coef(linear_mod),	coef(quad_mod))
})

sarmods_reg = data.frame(t(sarmods_reg))
names(sarmods_reg) = c('AIC_line','AIC_quad','ML_s2_line','ML_s2_quad',
	'line_lambda','line_int','line_slope','quad_lambda', 'quad_int','quad_slope','quad_sq')

write.csv(sarmods_reg, 'univariate_models_regS_AllSp_sar.csv', row.names=T)

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
write.csv(modcompare_sar, 'Univariate model shapes of regS SAR.csv', row.names=T)

# Read in previously saved tables
modcompare = read.csv('Univariate model shapes GLM-NB.csv', row.names=1)
modcompare_sar = read.csv('Univariate model shapes of regS SAR.csv', row.names=1)

# Which variables have AIC supported concave-down relationships?
sq_vars = rownames(subset(modcompare, concavity=='down'&type=='quadratic'))
sq_vars_sar = rownames(subset(modcompare_sar, concavity=='down'&type=='quadratic'))


##############################################################################
### Models for Variation Partitioning

# Define sets of predictors to be used in models
climvars = c('mat','iso','pseas','wetness','rain_lowRH')
LOvars = c(climvars, 'radiation','bark_moist_pct.ba','wood_SG.ba','light.mean', 'PC1')
LHvars = c('bark_moist_pct.rao.ba','wood_SG.rao.ba','lightDist.mean','PIE.ba.tree','propDead')
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
R_mod = glm.nb(richness~., data=cbind(use_data_test[,c('richness',Rvars)], sqdata[,paste(sq_vars[sq_vars %in% Rvars],2, 
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
RH_regmod = errorsarlm(regS ~ ., data = cbind(working_data_test[,c('regS',RHvars)], sqdata_reg[,paste(sq_vars_sar[sq_vars_sar %in% RHvars],2,sep='')]), listw=reg_listw)
RO_regmod = errorsarlm(regS ~ ., data = cbind(working_data_test[,c('regS',ROvars)], sqdata_reg[,paste(sq_vars_sar[sq_vars_sar %in% ROvars],2,sep='')]), listw=reg_listw)
R_regmod = errorsarlm(regS ~ ., data = cbind(working_data_test[,c('regS',RHvars, ROvars)], sqdata_reg[,paste(sq_vars_sar,2,sep='')]), listw=reg_listw)

# Calculate R2
null_mod = errorsarlm(regS ~ 1, data=working_data_test, listw=reg_listw)
Rs = sapply(list(RH_regmod, RO_regmod, R_regmod), function(x) calc_r2(x, null_mod))
names(Rs) = c('het','opt','full')
part_hetopt_reg = partvar2(Rs)

## Figure. Plot local-regional variation partioning analysis in three panels
barwide=1.5
use_shade = c('CC','55','99','FF')
use_color = c('#2415B0','#00BF32','#126A71')

svg('./Figures/New Analysis/variation partitioning loc-reg.svg', height=5, width=10)
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
	par(xpd=T)
	text(c(0,2,4),1.05, c('A','B','C'), cex=2, adj=c(0,0))
	par(xpd=F)
dev.off()


## Figure. Plot Heterogeneity - Optimality variation partitioning in 2 panels
svg('./Figures/New Analysis/variation partitioning het-opt.svg', height=5, width=8)
	par(mar=c(0,6,1.5,0))

	# Create plotting window
	plot(1,1, xlim=c(0,5.5), ylim=c(0,1), axes=F, ylab='', xlab='', type='n', cex.lab=2)
	
	## Add background lines
	usr = par('usr')
	segments(-.5, seq(0,1,.1), 5.25, seq(0,1,.1), lwd=2, col='grey70', lty=3)	

	## Background bars
	rect(0,0,barwide,1,col='white', lwd=3)	
	rect(4,0,4+barwide,1,col='white', lwd=3)

	## Full Model
	# Add rectangle for first component
	rect(0,0,barwide, sum(part_hetopt_reg[1]), lwd=3, col=paste(use_color[1],use_shade[2],sep=''))
	# Add rectangle for second component
	rect(0,sum(part_hetopt_reg[1:2]), barwide,sum(part_hetopt_reg[1:3]), lwd=3, col=paste(use_color[1],use_shade[1],sep=''))
	# Add rectangle for overlap
	rect(0,part_hetopt_reg[1],barwide,sum(part_hetopt_reg[1:2]), lwd=3, col=paste(use_color[1],use_shade[3],sep=''))

	## Among Optimality variables
	# Add rectangle for first component
	rect(4,0,4+barwide, sum(part_hetopt_loc[1]), lwd=3, col=paste(use_color[2],use_shade[2],sep=''))
	# Add rectangle for second component
	rect(4,sum(part_hetopt_loc[1:2]), 4+barwide,sum(part_hetopt_loc[1:3]), lwd=3, col=paste(use_color[2],use_shade[1],sep=''))
	# Add rectangle for overlap
	rect(4,part_hetopt_loc[1], 4+barwide,sum(part_hetopt_loc[1:2]), lwd=3, col=paste(use_color[2],use_shade[3],sep=''))

	## Add axis
	axis(2, las=1, cex.axis=2, lwd=3)
	mtext('Variation Explained', 2, 4, cex=2)
	
	# Add Text
	text_height = c(0,.1, .4, .8)
	text(2+barwide/2, text_height, labels=c('Heterogeneity','Shared','Optimality','Unexplained'), pos=3, offset=0, cex=2)

	# Add panel text A, B, C
	par(xpd=T)
	text(c(0,4),1.05, c('A','B'), cex=2, adj=c(0,0))
	par(xpd=F)
dev.off()










