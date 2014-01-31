## This script conducts supplemental analyses for the FIA Lichen project


# Load data and functions
source('./UNC/Projects/FIA Lichen/GitHub/FIA-Lichens/load_data.R')
source('./GitHub/FIA-Lichens/fia_lichen_analysis_functions.R')


#####################################################################
### Compare models with soil variables


### Variance Partitioning ###
library(MASS) # glm.nb
library(MuMIn) # r.squaredLR
library(hier.part) # combos()

# Read in table of predictor variable types
predtypes = read.csv('./SEM models/var_types.csv', row.names=1)
predtypes = subset(predtypes, !(rownames(predtypes) %in% c('FM','FH')))

# Define dataset- use trans_data where variables have been log or sqrt transformed but not scaled
use_data = trans_data[,colnames(trans_data)%in% rownames(predtypes)]

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

# Variables to include- leave out abundance and total tree circumference
use_vars = rownames(predtypes)[predtypes$type!='A']
use_vars = subset(use_vars, use_vars != 'totalCirc')

# Univariate model shapes
modcompare = read.csv('Univariate model shapes GLM-NB.csv', row.names=1)

# Which variables have AIC supported concave-down relationships?
sq_vars = rownames(subset(modcompare, concavity=='down'&type=='quadratic'))

# Using quadratic relationships for concave-down AIC supported models in variation partitioning
sq_df = use_data[,sq_vars]^2
colnames(sq_df) = paste(colnames(sq_df),'2', sep='')
use_data = cbind(use_data, sq_df)

# Subset data to only plots with soil data
soilRows = which(!is.na(use_data[,'soilPC1']))
use_data = use_data[soilRows,]

# Testing and training data sets
use_data_test = use_data[testplots[testplots$yrplot.id %in% rownames(use_data),'yrplot.id'],]
use_data_fit = use_data[fitplots[fitplots$yrplot.id %in% rownames(use_data),'yrplot.id'],]

## Variation Partitioning by Local / Regional

regionvars = subset(use_vars, use_vars %in% rownames(predtypes)[predtypes$scale=='R'])
localvars = subset(use_vars, use_vars %in% rownames(predtypes)[predtypes$scale=='L'])

region_mod = glm.nb(richness~., data=use_data_test[,c('richness', regionvars, colnames(sq_df)[sq_vars %in% regionvars])])
local_mod = glm.nb(richness~., data=use_data_test[,c('richness', localvars, colnames(sq_df)[sq_vars %in% localvars])])

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

## Variation Partitioning by Forest Heterogeneity / Mean Conditions
nichevars = subset(use_vars, use_vars %in% rownames(predtypes)[predtypes$type=='FH'])
resvars = subset(use_vars, use_vars %in% rownames(predtypes)[predtypes$type %in% c('FM','L')])

niche_mod = glm.nb(richness~., data=use_data_test[,c('richness', nichevars, colnames(sq_df)[sq_vars %in% nichevars])])
res_mod = glm.nb(richness~., data=use_data_test[,c('richness', resvars, colnames(sq_df)[sq_vars %in% resvars])])

AIC(region_mod, local_mod, niche_mod, res_mod) # region_mod is best

# Make a list of predictors
predlist = list(Heterogeneity=names(niche_mod$coefficients[-1]),
	Optimality=names(res_mod$coefficients[-1])
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

niche_res_partition = partvar2(Rs)

# Plot variation vartitioning in black and white
barwide = .6
use_part = local_regional_partition

svg('./Figures/partition variance among local vs regional include soil bw.svg', height=6, width=6*barwide)
	par(mar=c(0,4,0,0))

	# Create plotting window
	plot(1,1, xlim=c(0,barwide), ylim=c(0,1), axes=F, xlab='', ylab='', type='n')
	
	# Add rectangle for unexplained variation
	rect(0,sum(use_part[1:3]),barwide,1, lwd=3, col='white')
	
	# Add rectangle for regional model
	rect(0,0,barwide, sum(use_part[1:2]), lwd=3, col="#00000050")

	# Add rectangle for local model
	rect(0,use_part[1],barwide,sum(use_part[1:3]), lwd=3, col="#00000050")

	# Add axis
	axis(2, las=1, cex.axis=2, lwd=3)
	
	# Add partition labels
	lablocs = sapply(1:4, function(x) sum(use_part[0:(x-1)])+use_part[x]/2)
	text(barwide/2, lablocs, labels=names(use_part), cex=2, bg='white')
dev.off()


# Plot without unexplained variation
barwide=.6
use_part = niche_res_partition
svg('./Figures/partition variance among heterogeneity vs optimality include soil bw.svg', height=6, width=6*barwide)
	par(mar=c(0,0,0,5))

	# Create plotting window
	plot(1,1, xlim=c(0,barwide), ylim=c(0,.5), axes=F, xlab='', ylab='', type='n')
	
	# Add rectangle for regional model
	rect(0,0,barwide, sum(use_part[1:2]), lwd=3, col="#00000050")

	# Add rectangle for local model
	rect(0,use_part[1],barwide,sum(use_part[1:3]), lwd=3, col="#00000050")

	# Add axis
	axis(4, las=1, cex.axis=2, lwd=3, at=seq(0,.5,.1))
	
	# Add partition labels
	lablocs = sapply(1:3, function(x) sum(use_part[0:(x-1)])+use_part[x]/2)
	text(barwide/2, lablocs, labels=names(use_part)[1:3], cex=2)
dev.off()

### SEM ###

allsp = read.csv('./Soil/Soil_testdata_totaleffects.csv', row.names=1)
allsp_d = read.csv('./Soil/Soil_testdata_directeffects_richness.csv', row.names=1)
allsp_da = read.csv('./Soil/Soil_testdata_directeffects_abundance.csv', row.names=1)
allsp_i = read.csv('./Soil/Soil_testdata_indirecteffects_via_abundance.csv', row.names=1)
allsp_ir = read.csv('./Soil/Soil_testdata_indirecteffects_via_regS.csv', row.names=1)
allsp_if = read.csv('./Soil/Soil_testdata_indirecteffects_via_forest.csv', row.names=1)

allsp_ests = read.csv('./Soil/Soil_testdata_parameterEstimates.csv', row.names=1)

# How many predictors have significant total effects?
names(which(apply(allsp[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0))
tot_sig = allsp[which(apply(allsp[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0),]
tot_sig = tot_sig[order(abs(tot_sig$std.all)),]
total[order(abs(total$std.all)),]

allsp_d[order(abs(allsp_d$std.all)),]


## Compare direct effects vs. total effects

# Order variables from lowest to highest total effects
total = allsp # This changes based on what response variable is being analyzed
ordered_vars = rownames(total[order(total$std.all),])
ordered_vars = ordered_vars[!(ordered_vars %in% c('FH','FM'))] # Drop effects of FM and FH categories (they may be non-sensical)

# Put tables in same order
use_direct = allsp_d[ordered_vars,] # This changes based on what response variable is being analyzed
use_total = total[ordered_vars,]

library(lattice)

# Define range limits that will include 95% confidence intervals
myrange = range(c(use_total[,c('std.ci.lower','std.ci.upper')],
	use_direct[,c('std.ci.lower','std.ci.upper')]), na.rm=T)+c(-.04, .04)
myrange[1] = -1

# Color scheme: http://colorschemedesigner.com/#3341SsYrGvyw0
mypch = c(22,23)
mypcols = c('white','grey30')
myadj=.15
mytypes = expression('C'['R'],'C'['L'],'F'['H'],'F'['O'],'C'['R'],'R','') # symbols used in plot to denote variable types
names(mytypes)=c('C','L','FH','FM','P','R','A')

# Total and direct standardized effects on same graph
svg('./Figures/Standardized direct total effects on AllSp richness soil model.svg', height=13, width=19)
dotplot(as.numeric(factor(rownames(use_total), levels = ordered_vars))~std.all, data=use_total, 
	xlab=list('Standardized Effect',cex=3), ylab='',
	main='',cex.lab=3,aspect=1/1, xlim=myrange,
	panel=function(x,y){
	
		# Add horizontal boxes
		panel.rect(-2,1:length(ordered_vars)-.5, 2, 1:length(ordered_vars)+.5,
			col='white', border='grey50')
		
		# Add vertical line at 0
		panel.abline(v=0, col='grey30', lty=2, lwd=2)		
		
		# Add null distribution segments for direct effects
		panel.segments(use_direct$std.ci.lower, y+myadj,
			use_direct$std.ci.upper, y+myadj, 
			col='black', lwd=4.5, lend=1)
		# Add points for direct estimated effects
		panel.points(use_direct$std.all, y+myadj, col='black', fill=mypcols[2], pch=mypch[2], cex=3, lwd=3) 
		
		# Add null distribution segments for total effects
		panel.segments(use_total$std.ci.lower, y,
			use_total$std.ci.upper, y, 
			col='black', lwd=4.5, lend=1)
		# Add points for total estimated effects
		panel.points(x, y, col='black', fill=mypcols[1], pch=mypch[1], cex=3, lwd=3) 
	
		# Add text labeling the variable type
		panel.text(-.95, y, labels=mytypes[use_total$type], cex=2)
		
	},
	scales=list(y=list(labels=varnames[ordered_vars,'midName'], 
		cex=3, col='black'),
		x=list(cex=3, at=seq(-.8,.8,.2))),
	key=list(x=1, y=0, corner=c(1,0), lines=list(type='o', pch=mypch, fill=mypcols, lwd=3, pt.lwd=3),
		text=list(c('Total effect','Direct effect')),
		background='white', cex=3, divide=1, padding.text=5, border='black')
)
dev.off()


## Make table of indirect environmental effects

indirect = allsp_i
indirectR = allsp_ir
indirectF = allsp_if

use_total_sub = subset(use_total, type %in% c('C','L'))
use_direct_sub = subset(use_direct, type %in% c('C','L'))
use_indirect_sub = subset(indirect, type %in% c('C','L'))
order_clim = use_direct_sub[order(use_direct_sub$std.all),'predictor']

use_total_sub = use_total_sub[order_clim,]
use_direct_sub = use_direct_sub[order_clim,]
use_indirect_sub = use_indirect_sub[order_clim,]
use_indirectR_sub = indirectR[order_clim,]
use_indirectF_sub = indirectF[c(paste(order_clim,'FM', sep='_'),paste(order_clim,'FH', sep='_')),]

climEff_tab = data.frame(Predictor=varnames[use_direct_sub$predictor,'displayName'], 
	direct = use_direct_sub$std.all,
	indirectR = use_indirectR_sub$std.all,
	indirectA = use_indirect_sub$std.all,
	indirectFH = subset(use_indirectF_sub, Ftype=='FH')$std.all,
	indirectFM = subset(use_indirectF_sub, Ftype=='FM')$std.all,
	directSig = apply(use_direct_sub[,c('std.ci.lower','std.ci.upper')], 1, prod)>0,
	indirectRSig = apply(use_indirectR_sub[,c('std.ci.lower','std.ci.upper')], 1, prod)>0,
	indirectASig = apply(use_indirect_sub[,c('std.ci.lower','std.ci.upper')], 1, prod)>0,
	indirectFHSig = apply(subset(use_indirectF_sub, Ftype=='FH')[,c('std.ci.lower','std.ci.upper')], 1, prod)>0,
	indirectFMSig = apply(subset(use_indirectF_sub, Ftype=='FM')[,c('std.ci.lower','std.ci.upper')], 1, prod)>0
)

write.csv(climEff_tab, './Soil/Compare effects climate variables soil model.csv', row.names=F)


## Paths for Soil model

# Only look at estimates corresponding to paths
paths = subset(allsp_ests, op=='~')

# Calculate which paths are significant
sigpaths = paths[which(apply(paths[,c('std.ci.lower','std.ci.upper')],1,prod)>0),]

# Histogram of significant paths
hist(sigpaths$std.all)
sigpaths$sigcat = cut(abs(sigpaths$std.all), c(0,.2,.4,.6,.8,1))

# Find all variables in the model
allvars = unique(c(paths$rhs,paths$lhs))

# Find all variables involved in significant relationships
sigvars  = unique(c(sigpaths$lhs, sigpaths$rhs))

## Create a matrix of variable locations
var_locs = matrix(0, nrow=length(allvars), ncol=2, byrow=T)
colnames(var_locs) = c('X','Y')
rownames(var_locs) = allvars

unit = 1

# Start with richness and abundance in center
var_locs['tot_abun_log',] = c(-1,0)
var_locs['lichen_rich',] = c(1,0)

# Add pollution and regS to corners
var_locs['totalNS',] = c(-2.5,2.5)
var_locs['regS',] = c(2.5,2.5)

# Add climate and local environment vars to top and bottom
c_vars = rownames(subset(predtypes[allvars,], type=='C'))
var_locs[c_vars,] = cbind(3/(length(c_vars)-1)*(0:(length(c_vars)-1))-1.5,rep(3,length(c_vars)))
l_vars = rownames(subset(predtypes[allvars,], type=='L'))
var_locs[l_vars,] = cbind(2.5/(length(l_vars)-1)*(0:(length(l_vars)-1))-1.5,rep(-2.5,length(l_vars)))


# Add forest vars to sides
fh_vars = rownames(subset(predtypes[allvars,], type=='FH'))
fm_vars = rownames(subset(predtypes[allvars,], type=='FM'))

var_locs[fh_vars,]=cbind(rep(3,length(fh_vars)), 3.5/(length(fh_vars)-1)*(0:(length(fh_vars)-1))-2)
var_locs[fm_vars,]=cbind(rep(-3,length(fm_vars)), 3.5/(length(fm_vars)-1)*(0:(length(fm_vars)-1))-2)

var_locs=data.frame(var_locs)

# Make plot
mycex=1.5

svg('./Figures/full soil model significant paths.svg', height=5, width=8.5)
par(mar=c(0,0,0,0))
par(lend="butt")
plot(var_locs, xlim=c(-6,6), ylim=c(-5,6), type='n', axes=F, xlab='', ylab='')
for(i in order(abs(sigpaths$std.all))){
	from = var_locs[sigpaths[i,'rhs'],]
	to = var_locs[sigpaths[i,'lhs'],]
	thickness = as.numeric(sigpaths[i,'sigcat'])*mycex
	col=c('red','black')[(sigpaths[i,'std.all']>0)+1]

	arrows(from$X, from$Y, to$X, to$Y, length=0.15, lwd=thickness, col=col)
}
text(var_locs[fm_vars,], labels=varnames[fm_vars,'midName'], pos=2, offset=.25)
text(var_locs[fh_vars,], labels=varnames[fh_vars,'midName'], pos=4, offset=.25)
text(var_locs[c_vars,], labels=varnames[c_vars,'midName'],adj=-.1, srt=45)
text(var_locs[c('lichen_rich','tot_abun_log'),], labels=c('Local\nrichness','Abundance'),pos=1, offset=1)
text(var_locs[l_vars,], labels=varnames[l_vars,'midName'], adj=-.1, srt=-45)
text(var_locs['totalNS',], labels=varnames['totalNS','midName'], pos=2, offset=.25)
text(var_locs['regS',], labels=varnames['regS','midName'], pos=4, offset=.25)

legend(x=-5.5, y=5, xjust=0, yjust=1, c('0.0-0.2','0.2-0.4','0.4-0.6','0.6-0.8'), lwd=mycex*(1:4), 
	bty='n')
legend(x=-6, y=5, xjust=0, yjust=1, rep('', 4), lwd=mycex*(1:4), col='red', bty='n')
text(-5.5, 5, '+', cex=2, adj=-1)
text(-6, 5, '-', cex=2, adj=-2.3)

dev.off()


######################################################################
### Compare models without abundance
noabun_ests = read.csv('./SEM Models/No Abundance/AllSp_testdata_noabun_parameterEstimates.csv', row.names=1)

noabun = read.csv('./SEM Models/No Abundance/AllSp_testdata_noabun_totaleffects.csv', row.names=1)
noabun_d = read.csv('./SEM Models/No Abundance/AllSp_testdata_noabun_directeffects_richness.csv', row.names=1)
noabun_ir = read.csv('./SEM Models/No Abundance/AllSp_testdata_noabun_indirecteffects_via_regS.csv', row.names=1)
noabun_if = read.csv('./SEM Models/No Abundance/AllSp_testdata_noabun_indirecteffects_via_forest.csv', row.names=1)

allsp = read.csv('./SEM models/AllSp_testdata_totaleffects.csv', row.names=1)
allsp_d = read.csv('./SEM models/AllSp_testdata_directeffects_richness.csv', row.names=1)
allsp_ir = read.csv('./SEM models/AllSp_testdata_indirecteffects_via_regS.csv', row.names=1)
allsp_if = read.csv('./SEM models/AllSp_testdata_indirecteffects_via_forest.csv', row.names=1)

ordered_vars = rownames(allsp)[order(allsp$std.all)]

data.frame(predictor=ordered_vars,regular_model=allsp[ordered_vars,'std.all'], noabun_model=noabun[ordered_vars,'std.all'])

## Compare total effect between two models
use_vars = rownames(allsp)
use_vars = subset(use_vars, !(use_vars %in% c('FM','FH','abun_log')))

mypch = c(21:23) # Denote (climate, pollution, regional) or (env, fh, fm)
mycol=c('white','grey30') # Denote regional vs local
myfact = factor(allsp[use_vars,'type'], levels = c('C','P','R','L','FH','FM'))
myfactpch = (as.numeric(myfact)-1)%%3+1
myfactcol = ifelse(allsp[use_vars,'type']%in% c('C','P','R'), 1, 2)

svg('./Figures/compare total effect abun and noabun models.svg', height=6, width=6 )
par(mar=c(4,4,1,1))
plot(allsp[use_vars,'std.all'], noabun[use_vars,'std.all'], type='n', 
	xlim=c(-.6,.6), ylim=c(-.6,.6), las=1, ylab='Model without abundance',
	xlab='Model with abundance', cex.axis=1.2, cex.lab=1.2)
abline(h=0,v=0)
abline(0,1)
arrows(allsp[use_vars,'std.ci.lower'], noabun[use_vars,'std.all'],
	allsp[use_vars,'std.ci.upper'], noabun[use_vars,'std.all'],
	code=3, angle=90, lwd=2, length=0.05)
arrows(allsp[use_vars,'std.all'], noabun[use_vars,'std.ci.lower'],
	allsp[use_vars,'std.all'], noabun[use_vars,'std.ci.upper'],
	code=3, angle=90, lwd=2, length=0.05)
points(allsp[use_vars,'std.all'], noabun[use_vars,'std.all'], 
	pch=mypch[myfactpch], bg=mycol[myfactcol], lwd=2, col='black', cex=2)
legend('topleft',levels(myfact), pch=mypch, pt.bg=rep(mycol, each=3), 
	pt.lwd=2, bg='white', box.lwd=1, ncol=2)
dev.off()

svg('./Figures/compare indirect effect abun and noabun models.svg', height=6, width=12)
par(mar=c(4,4,1,1))
par(mfrow=c(1,2))

# A function that positions text in corners
quadpos = function(x,a,b){
	positions = list(c(-a,-b),c(1+a,-b),c(1+a,1+b),c(-a,1+b))
	positions[[x]]
}

# Indirect effecst via regional richness
noabun_ir = noabun_ir[rownames(allsp_ir),] # Put in same order
plot(allsp_ir[,'std.all'], noabun_ir[,'std.all'], type='n', main='Indirect Effect via Regional Richness', 
	xlim=c(-.05,.2), ylim=c(-.05, .2),las=1, ylab='Model without abundance',
	xlab='Model with abundance', cex.axis=1.2, cex.lab=1.2)
abline(h=0,v=0)
abline(0,1)
arrows(allsp_ir[,'std.ci.lower'], noabun_ir[,'std.all'],
	allsp_ir[,'std.ci.upper'], noabun_ir[,'std.all'],
	code=3, angle=90, lwd=2, length=0.05)
arrows(allsp_ir[,'std.all'], noabun_ir[,'std.ci.lower'],
	allsp_ir[,'std.all'], noabun_ir[,'std.ci.upper'],
	code=3, angle=90, lwd=2, length=0.05)
points(allsp_ir[,'std.all'], noabun_ir[,'std.all'], 
	pch=mypch[2], bg=mycol[1], lwd=2, col='black', cex=2)
textvars=rownames(allsp_ir)
for(i in 1:length(textvars)){
text(allsp_ir[textvars[i],'std.all'],noabun_ir[textvars[i],'std.all'], 
	labels=varnames[textvars[i],'midName'], adj=quadpos(c(3,1,1,1,4)[i], .1, .5)
)
}

noabun_if = noabun_if[rownames(allsp_if),] # Put in same order
fh_vars = rownames(noabun_if)[grep('FH',rownames(noabun_if))]
fm_vars = rownames(noabun_if)[grep('FM',rownames(noabun_if))]

plot(allsp_if[,'std.all'], noabun_if[,'std.all'], type='n', main='Indirect Effect via Forest Structure', 
	xlim=c(-.1,.11), ylim=c(-.1, .11), las=1, ylab='Model without abundance',
	xlab='Model with abundance', cex.axis=1.2, cex.lab=1.2)
abline(h=0,v=0)
abline(0,1)
arrows(allsp_if[,'std.ci.lower'], noabun_if[,'std.all'],
	allsp_if[,'std.ci.upper'], noabun_if[,'std.all'],
	code=3, angle=90, lwd=2, length=0.05)
arrows(allsp_if[,'std.all'], noabun_if[,'std.ci.lower'],
	allsp_if[,'std.all'], noabun_if[,'std.ci.upper'],
	code=3, angle=90, lwd=2, length=0.05)
points(allsp_if[fh_vars,'std.all'], noabun_if[fh_vars,'std.all'], 
	pch=mypch[2], bg=mycol[1], lwd=2, col='black', cex=2)
points(allsp_if[fm_vars,'std.all'], noabun_if[fm_vars,'std.all'], 
	pch=mypch[3], bg=mycol[2], lwd=2, col='black', cex=2)
legend('bottomright',c('Heterogeneity','Optimality'), pt.bg=mycol, pch=mypch[2:3],
	bty='n', pt.cex=2, pt.lwd=2)
dev.off()

### Path Diagram

# Only look at estimates corresponding to paths
paths = subset(noabun_ests, op=='~')

# Calculate which paths are significant
sigpaths = paths[which(apply(paths[,c('std.ci.lower','std.ci.upper')],1,prod)>0),]
sigpaths$sigcat = cut(abs(sigpaths$std.all), c(0,.2,.4,.6,.8,1))

# Find all variables in the model
allvars = unique(c(paths$rhs,paths$lhs))

# Find all variables involved in significant relationships
sigvars  = unique(c(sigpaths$lhs, sigpaths$rhs))

## Create a matrix of variable locations
var_locs = matrix(0, nrow=length(allvars), ncol=2, byrow=T)
colnames(var_locs) = c('X','Y')
rownames(var_locs) = allvars

unit = 1

# Start with richness in center
var_locs['lichen_rich',] = c(0,0)

# Add pollution and regS to corners
var_locs['totalNS',] = c(-2.5,2.5)
var_locs['regS',] = c(2.5,2.5)

# Add climate and local environment vars to top and bottom
c_vars = rownames(subset(predtypes[allvars,], type=='C'))
var_locs[c_vars,] = cbind(3/(length(c_vars)-1)*(0:(length(c_vars)-1))-1.5,rep(3,length(c_vars)))
var_locs['radiation',] = c(0,-2.5)

# Add forest vars to sides
fh_vars = rownames(subset(predtypes[allvars,], type=='FH'))
fm_vars = rownames(subset(predtypes[allvars,], type=='FM'))

var_locs[fh_vars,]=cbind(rep(3,length(fh_vars)), 3.5/(length(fh_vars)-1)*(0:(length(fh_vars)-1))-2)
var_locs[fm_vars,]=cbind(rep(-3,length(fm_vars)), 3.5/(length(fm_vars)-1)*(0:(length(fm_vars)-1))-2)

var_locs=data.frame(var_locs)

# Make plot
mycex=1.5

svg('./Figures/full model no abundance significant paths.svg', height=5, width=8.5)
par(mar=c(0,0,0,0))
par(lend="butt")
plot(var_locs, xlim=c(-6,6), ylim=c(-2.5,6), type='n', axes=F, xlab='', ylab='')
for(i in order(abs(sigpaths$std.all))){
	from = var_locs[sigpaths[i,'rhs'],]
	to = var_locs[sigpaths[i,'lhs'],]
	thickness = as.numeric(sigpaths[i,'sigcat'])*mycex
	col=c('red','black')[(sigpaths[i,'std.all']>0)+1]

	arrows(from$X, from$Y, to$X, to$Y, length=0.15, lwd=thickness, col=col)
}
text(var_locs[fm_vars,], labels=varnames[fm_vars,'midName'], pos=2, offset=.25)
text(var_locs[fh_vars,], labels=varnames[fh_vars,'midName'], pos=4, offset=.25)
text(var_locs[c_vars,], labels=varnames[c_vars,'midName'],adj=-.1, srt=45)
text(var_locs[c('lichen_rich','tot_abun_log'),], labels=c('Local\nrichness','Abundance'),pos=1, offset=1)
text(var_locs['radiation',], labels='Solar radiation', pos=1, offset=.25)
text(var_locs['totalNS',], labels=varnames['totalNS','midName'], pos=2, offset=.25)
text(var_locs['regS',], labels=varnames['regS','midName'], pos=4, offset=.25)

legend(x=-5.5, y=5, xjust=0, yjust=1, c('0.0-0.2','0.2-0.4','0.4-0.6','0.6-0.8'), lwd=mycex*(1:4), 
	bty='n')
legend(x=-6, y=5, xjust=0, yjust=1, rep('', 4), lwd=mycex*(1:4), col='red', bty='n')
text(-5.5, 5, '+', cex=2, adj=-1)
text(-6, 5, '-', cex=2, adj=-2.3)

dev.off()










