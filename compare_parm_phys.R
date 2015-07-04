# This script is used to conduct comparative analyses of FIA lichen data for Parmeliaceae vs Physciaceae richness


#####################################################################
### Read in Data and Functions

source('./UNC/Projects/FIA Lichen/GitHub/FIA-Lichens/load_data.R')
source('./GitHub/FIA-Lichens/fia_lichen_analysis_functions.R')

# Read in table of predictor variable types
predtypes = read.csv('predictors.csv', row.names=1)


####################################################################
### SEM - all analyses conducted on Direct Regional Paths model

library(lavaan)
library(semTools)
library(lattice)

# Table of parameter estimates
allsp_ests = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_AllSp_testdata_parameterEstimates.csv')
parm_ests = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Parm_testdata_parameterEstimates.csv')
phys_ests = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Phys_testdata_parameterEstimates.csv')

# Total effects
allsp = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_AllSp_testdata_totaleffects.csv', row.names=1)
parm = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Parm_testdata_totaleffects.csv', row.names=1)
phys = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Phys_testdata_totaleffects.csv', row.names=1)

# Direct effects on local richness
allsp_d = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_AllSp_testdata_directeffects_richness.csv', row.names=1)
parm_d = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Parm_testdata_directeffects_richness.csv', row.names=1)
phys_d = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Phys_testdata_directeffects_richness.csv', row.names=1)

# Direct effecst on regional richness
allsp_dr = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_AllSp_testdata_directeffects_regS.csv', row.names=1)
parm_dr = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Parm_testdata_directeffects_regS.csv', row.names=1)
phys_dr = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Phys_testdata_directeffects_regS.csv', row.names=1)

# Indirect effects on local richness via local abundance
allsp_i = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_AllSp_testdata_indirecteffects_via_abundance.csv', row.names=1)
parm_i = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Parm_testdata_indirecteffects_via_abundance.csv', row.names=1)
phys_i = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Phys_testdata_indirecteffects_via_abundance.csv', row.names=1)

# Indirect effects of regional scale predictors on local richness
allsp_ir = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_AllSp_testdata_regionalvars_indirecteffects.csv')
parm_ir = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Parm_testdata_regionalvars_indirecteffects.csv')
phys_ir = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Phys_testdata_regionalvars_indirecteffects.csv')

# Direct effects on local abundance
allsp_da = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_AllSp_testdata_directeffects_abundance.csv', row.names=1)
parm_da = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Parm_testdata_directeffects_abundance.csv', row.names=1)
phys_da = read.csv('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Phys_testdata_directeffects_abundance.csv', row.names=1)

# Model objects
load('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Parm_testdata_output.RData')
parm_fit = regTorich_nopol_fit
load('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_Phys_testdata_output.RData')
phys_fit = regTorich_nopol_fit
load('./SEM models/Oct2014 No Measurement Error/regTorich_nopol_AllSp_testdata_output.RData')
allsp_fit = regTorich_nopol_fit

## Assess Model Fit
examine_mod = phys_fit

# Residuals between observed and predicted covariances
res = resid(examine_mod)$cov
res_num = res[lower.tri(res, diag=T)]
resdf = data.frame(residual=res_num, 
	var1=matrix(rownames(res), nrow=nrow(res), ncol=ncol(res))[lower.tri(res, diag=T)],
	var2=matrix(colnames(res), nrow=nrow(res), ncol=ncol(res), byrow=T)[lower.tri(res, diag=T)])
resdf = resdf[order(abs(resdf$residual), decreasing=T),]
big_res = subset(resdf, abs(residual)>0.1)
rich_res = subset(resdf, var1%in%c('Phys_log', 'regPhys')|var2%in%c('Phys_log', 'regPhys')) # This changes based on the model

write.table(big_res, './Parm-Phys/phys regTorich model large residuals.txt', sep='\t', row.names=F)
write.table(rich_res, './Parm-Phys/phys regTorich model richness residuals.txt', sep='\t', row.names=F)



### Compare significance of effects across taxa

# Total effects
tot_sig = data.frame(AllSp_sig = apply(allsp[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0,
	Parm_sig = apply(parm[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0,
	Phys_sig = apply(phys[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0)
tot_eff = data.frame(AllSp = allsp$std.all,
	Parm = parm$std.all,
	Phys = phys$std.all)

tot_sig = cbind(tot_eff, tot_sig)
parm_order = order(abs(tot_sig$Parm), decreasing=T)
phys_order = order(abs(tot_sig$Phys), decreasing=T)
subset(tot_sig[parm_order,], Parm_sig)
tot_sig = cbind(tot_sig, predtypes[rownames(tot_sig),c('mode','scale')])
aggregate(Parm_sig~scale+mode, data=tot_sig, FUN=sum) # counts number of significant effects
aggregate(Phys_sig~scale+mode, data=tot_sig, FUN=sum)

# Direct effects on local richness
dir_sig = data.frame(AllSp = apply(allsp_d[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0,
	Parm = apply(parm_d[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0,
	Phys = apply(phys_d[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0)
dir_sig = cbind(dir_sig, predtypes[rownames(dir_sig),c('mode','scale')])
aggregate(Parm~scale+mode, data=dir_sig, FUN=sum) # counts number of significant effects
aggregate(Phys~scale+mode, data=dir_sig, FUN=sum)

# Direct effects on regional richness
dir_sig = data.frame(AllSp = apply(allsp_dr[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0,
	Parm = apply(parm_dr[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0,
	Phys = apply(phys_dr[,c('std.ci.lower', 'std.ci.upper')], 1, prod)>0)
dir_sig = cbind(dir_sig, predtypes[rownames(dir_sig),c('mode','scale')])
aggregate(Parm~scale+mode, data=dir_sig, FUN=sum) # counts number of significant effects
aggregate(Phys~scale+mode, data=dir_sig, FUN=sum)


### Plot effects

# Total effects
tot_parm = rbind(parm, parm_d[c('regS','tot_abun_log'),])
tot_phys = rbind(phys, phys_d[c('regS','tot_abun_log'),])
ordered_vars = rownames(tot_parm)[order(abs(tot_parm$std.all - tot_phys$std.all))]
use_parm = tot_parm[ordered_vars,]
use_phys = tot_phys[ordered_vars,]

# Direct effects
rownames(parm_d)==rownames(phys_d)
ordered_vars = rownames(parm_d)[order(abs(parm_d$std.all - phys_d$std.all))]
use_parm = parm_d[ordered_vars,]
use_phys = phys_d[ordered_vars,]

# Separate direct local effects
ordered_vars = rownames(parm_d)[order(parm_d$std.all)]
ordered_vars = ordered_vars[(predtypes[ordered_vars,'scale']=='local')|ordered_vars=='regS']
use_parm = parm_d[ordered_vars,]
use_phys = phys_d[ordered_vars,]

# Separate direct regional effects
ordered_vars = rownames(parm_dr)[order(parm_dr$std.all)]
use_parm = parm_dr[ordered_vars,]
use_phys = phys_dr[ordered_vars,]

# Direct effects on abundance
rownames(parm_da)==rownames(phys_da)
ordered_vars = rownames(parm_da)[order(abs(parm_da$std.all - phys_da$std.all))]
use_parm = parm_da[ordered_vars,]
use_phys = phys_da[ordered_vars,]

# Make plots - change the title and x-axis label in the plot
mycols = matrix(c('#b3b5ffff','#6b6dd7ff','#8dff94ff','#38af4fff'), nrow=2)
mycolsbw = c('grey80','white')
colnames(mycols) = c('regional','local')
rownames(mycols) = c('het','opt')
names(mycolsbw) = c('regional','local')
mypch = c(22,23)
mypcols = c('white','grey30')
myadj=.15


# Define range limits that will include 95% confidence intervals
myrange = range(c(use_parm[,c('std.ci.lower','std.ci.upper')],
	use_phys[,c('std.ci.lower','std.ci.upper')]), na.rm=T)+c(-.04, .04)

# Parmeliaceae and Physciaceae effects on same graph
# for all variables use height = 7, for local use height=5, for regional use height=
# for separate regional direct effects comment out colorcombos['regS','mode'] = 'opt'
# for regional direct use aspect 5/8, others use 10/8
svg('./Parm-Phys/Standardized direct effects on Parm Phys local abundance.svg', height=5, width=8)
dotplot(as.numeric(factor(rownames(use_parm), levels = ordered_vars))~std.all, data=use_phys, 
	xlab=list('Standardized Direct Effect',cex=1), ylab='',
	main='',cex.lab=1,aspect=10/8, xlim=myrange,
	panel=function(x,y){
	
		# Add horizontal boxes colored by variable type
		colorcombos = predtypes[ordered_vars, c('mode','scale')]
		#colorcombos['regS','mode'] = 'opt'
		colororder = apply(colorcombos, 1, function(x) mycols[x[1],x[2]])
		panel.rect(myrange[1]-0.01,1:length(ordered_vars)-.5, myrange[2]+0.01, 1:length(ordered_vars)+.5,
			col=colororder, border='grey50')
		
		# Add vertical line at 0
		panel.abline(v=0, col='grey30', lty=2, lwd=2)		
		
		# Add null distribution segments for Parmeliaceae
		panel.segments(use_parm$std.ci.lower, y+myadj,
			use_parm$std.ci.upper, y+myadj, 
			col='black', lwd=2, lend=2)
		# Add points for estimated effects
		panel.points(use_parm$std.all, y+myadj, col='black', fill=mypcols[2], pch=mypch[2], cex=1, lwd=2, ljoin=1) 
		
		# Add null distribution segments for Physciaceae
		panel.segments(use_phys$std.ci.lower, y,
			use_phys$std.ci.upper, y, 
			col='black', lwd=2, lend=2)
		# Add points for estimated effects
		panel.points(x, y, col='black', fill=mypcols[1], pch=mypch[1], cex=1, lwd=2, ljoin=1) 
	
},
	scales=list(y=list(labels=varnames[ordered_vars,'midName'], 
		cex=1, col='black'),
		x=list(cex=1, tick.number=8))#,
	#key=list(x=1, y=0, corner=c(1,0), lines=list(type='o', pch=mypch, fill=mypcols, lwd=2, pt.lwd=2, ljoin=2, lend=2),
	#	text=list(c('Physciaceae','Parmaliaceae')),
	#	background='#FFFFFFCC', cex=1, divide=1, padding.text=5, border='black')
)
dev.off()


### Compare effects by plotting one effect on x-axis and one on y-axis with Parm and Phys side-by-side


# Direct effects on local richness vs indirect effects via abundance
ordered_vars = rownames(parm_i)
use_xparm = parm_d[ordered_vars,]
use_yparm = parm_i[ordered_vars,]
use_xphys = phys_d[ordered_vars,]
use_yphys = phys_i[ordered_vars,]

# Direct effects of local vs regional climate on local richness
ordered_vars = c('mat','iso','pseas','wetness','rain_lowRH')
use_xparm = parm_d[ordered_vars,]
use_yparm = parm_d[paste(ordered_vars, 'reg_mean', sep='_'),]
use_xphys = phys_d[ordered_vars,]
use_yphys = phys_d[paste(ordered_vars, 'reg_mean', sep='_'),]

# Define symbols
mypch = c(22,23)
mycol = c('white','grey30')
myfact = factor(predtypes[ordered_vars,'mode'], levels=c('het','opt'))

# Find variables with significant effecst
sig_xparm = apply(use_xparm[,c('std.ci.lower', 'std.ci.upper')], 1, function(x) prod(x)>0)
sig_yparm = apply(use_yparm[,c('std.ci.lower', 'std.ci.upper')], 1, function(x) prod(x)>0)
sig_xphys = apply(use_xphys[,c('std.ci.lower', 'std.ci.upper')], 1, function(x) prod(x)>0)
sig_yphys = apply(use_yphys[,c('std.ci.lower', 'std.ci.upper')], 1, function(x) prod(x)>0)

# Calculate plot range
myrange = range(c(use_xparm[,c('std.ci.lower','std.ci.upper')], use_yparm[,c('std.ci.lower','std.ci.upper')],
	use_xphys[,c('std.ci.lower','std.ci.upper')], use_yphys[,c('std.ci.lower','std.ci.upper')]))
myrange = max(abs(myrange))*c(-1,1)


# Open graphics file - this name changes depending on what is being compared
svg('./Parm-Phys/compare local regional climate effects parm vs phys.svg', height=5, width=10)
par(mar=c(4,4,2,1))
par(mfrow=c(1,2))

# Make plot - x and y axis labels change based on what is being compared
plot(0,0, type='n', xlim=myrange, ylim=myrange, las=1, ylab='Regional Variable Direct Effect',
	xlab='Local Variable Direct Effect', cex.axis=1, cex.lab=1)
# Draw quadrants
usr=par('usr')
xvals = seq(usr[1],usr[2],length.out=3)
fun1 = function(x) x
fun2 = function(x) -x
polygon(xvals, pmax(fun1(xvals), fun2(xvals)), col='grey80')
polygon(xvals, pmin(fun1(xvals), fun2(xvals)), col='grey80')
abline(h=0,v=0)

# Decide on points to plot
use_pts = rep(T, length(ordered_vars)) #sig_xparm|sig_yparm

# Draw CI
arrows(use_xparm[use_pts,'std.ci.lower'], use_yparm[use_pts,'std.all'],
	use_xparm[use_pts,'std.ci.upper'], use_yparm[use_pts,'std.all'],
	code=3, angle=90, lwd=2, length=0.05)
arrows(use_xparm[use_pts,'std.all'], use_yparm[use_pts,'std.ci.lower'],
	use_xparm[use_pts,'std.all'], use_yparm[use_pts,'std.ci.upper'],
	code=3, angle=90, lwd=2, length=0.05)

# Add points
points(use_xparm[use_pts,'std.all'], use_yparm[use_pts,'std.all'], 
	pch=mypch[myfact[use_pts]], bg=mycol[myfact[use_pts]], lwd=2, col='black', cex=2)

# Add legend
#legend('topright',c('Heterogeneity','Optimality'), pch=mypch, pt.bg=mycol, 
#	pt.lwd=2, bg='white', box.lwd=1, pt.cex=2)

# Add title
mtext('A. Parmeliaceae',3,.5, adj=0)

# Make plot
plot(0,0, type='n', xlim=myrange, ylim=myrange, las=1, ylab='Regional Variable Direct Effect',
	xlab='Local Variable Direct Effect', cex.axis=1, cex.lab=1)
# Draw quadrants
usr=par('usr')
xvals = seq(usr[1],usr[2],length.out=3)
fun1 = function(x) x
fun2 = function(x) -x
polygon(xvals, pmax(fun1(xvals), fun2(xvals)), col='grey80')
polygon(xvals, pmin(fun1(xvals), fun2(xvals)), col='grey80')
abline(h=0,v=0)

# Decide on points to plot
use_pts = rep(T, length(ordered_vars)) #sig_xphys|sig_yphys

# Draw CI
arrows(use_xphys[use_pts,'std.ci.lower'], use_yphys[use_pts,'std.all'],
	use_xphys[use_pts,'std.ci.upper'], use_yphys[use_pts,'std.all'],
	code=3, angle=90, lwd=2, length=0.05)
arrows(use_xphys[use_pts,'std.all'], use_yphys[use_pts,'std.ci.lower'],
	use_xphys[use_pts,'std.all'], use_yphys[use_pts,'std.ci.upper'],
	code=3, angle=90, lwd=2, length=0.05)

# Add points
points(use_xphys[use_pts,'std.all'], use_yphys[use_pts,'std.all'], 
	pch=mypch[myfact[use_pts]], bg=mycol[myfact[use_pts]], lwd=2, col='black', cex=2)

# Add title
mtext('B. Physciaceae',3,.5, adj=0)

# Add legend
#legend('topright',c('Heterogeneity','Optimality'), pch=mypch, pt.bg=mycol, 
#	pt.lwd=2, bg='white', box.lwd=1, pt.cex=2)

dev.off()


### Plot Parm vs Phys effects
predtypes['regS','mode'] = 'other'

# Direct effects on local richness of parm vs phys
ordered_vars = rownames(parm_d)
ordered_vars = ordered_vars[predtypes[ordered_vars, 'scale']=='local']
ordered_vars = c(ordered_vars, 'regS')
use_xdata = parm_d[ordered_vars,]
use_ydata = phys_d[ordered_vars,]
add_text = data.frame(vars = c('tot_abun_log','regS'), pos = c(2,4))

# Direct effects on regional richness of parm vs phys
ordered_vars = rownames(parm_dr)
use_xdata = parm_dr[ordered_vars,]
use_ydata = phys_dr[ordered_vars,]
add_text = data.frame(vars = c('mat_reg_mean','wetness_reg_var','regS_tree','PIE.ba.tree'), pos = c(4,4,1,4))

# Direct effects on local abundance of parm vs phys
ordered_vars = rownames(parm_da)
use_xdata = parm_da[ordered_vars,]
use_ydata = phys_da[ordered_vars,]
add_text = data.frame(vars = c('LogSeed.rao.ba','wetness','mat'), pos = c(4,3,3))



# Define symbols
mycols = matrix(c('#b3b5ffff','#6b6dd7ff','#6b6dd7ff','#8dff94ff','#38af4fff','#38af4fff'), nrow=3)
colnames(mycols) = c('regional','local')
rownames(mycols) = c('het','opt','other')
mypch = c(22,23,21)

use_col = sapply(ordered_vars, function(x) mycols[predtypes[x,'mode'],predtypes[x,'scale']])
use_pch = mypch[factor(predtypes[ordered_vars,'mode'], levels=c('het','opt','other'))]

# Find variables with significant effecst
sig_xdata = apply(use_xdata[,c('std.ci.lower', 'std.ci.upper')], 1, function(x) prod(x)>0)
sig_ydata = apply(use_ydata[,c('std.ci.lower', 'std.ci.upper')], 1, function(x) prod(x)>0)

# Calculate plot range
myrange = range(c(use_xdata[,c('std.ci.lower','std.ci.upper')], use_ydata[,c('std.ci.lower','std.ci.upper')]))
#myrange = max(abs(myrange))*c(-1,1)

# Open graphics file - this name changes depending on what is being compared
svg('./Parm-Phys/compare parm phys local abundance direct effects.svg', height=5, width=5)
par(mar=c(4,4,2,1))

# Make plot - x and y axis labels change based on what is being compared
plot(0,0, type='n', xlim=myrange, ylim=myrange, las=1, ylab='Physciaceae Abundance Direct Effect',
	xlab='Parmeliaceae Abundance Direct Effect', cex.axis=1, cex.lab=1)

# Draw quadrants
abline(h=0,v=0)
abline(0,1)

# Decide on points to plot
use_pts = sig_xdata|sig_ydata #rep(T, length(ordered_vars)) 

# Draw CI
arrows(use_xdata[use_pts,'std.ci.lower'], use_ydata[use_pts,'std.all'],
	use_xdata[use_pts,'std.ci.upper'], use_ydata[use_pts,'std.all'],
	code=3, angle=90, lwd=2, length=0.05)
arrows(use_xdata[use_pts,'std.all'], use_ydata[use_pts,'std.ci.lower'],
	use_xdata[use_pts,'std.all'], use_ydata[use_pts,'std.ci.upper'],
	code=3, angle=90, lwd=2, length=0.05)

# Add points
points(use_xdata[use_pts,'std.all'], use_ydata[use_pts,'std.all'], 
	pch=use_pch[use_pts], bg=use_col[use_pts], lwd=2, col='black', cex=1.5)

# Add labels
text(use_xdata[add_text$vars,'std.all'], use_ydata[add_text$vars,'std.all'],
	labels=varnames[add_text$vars,'shortName'], pos=add_text$pos, offset=1.5)

dev.off()


# Add legend
#legend('topleft',c('Heterogeneity','Optimality'), pch=mypch, pt.bg=mycol, 
#	pt.lwd=2, bg='white', box.lwd=1, pt.cex=2)

### Compare paths that differ between Parm and Phys SEM

physpaths = subset(phys_ests, op=='~') 
parmpaths = subset(parm_ests, op=='~') 

# Which paths differ by more than 0.1?
diff_cutoff = .1
physpaths$diff = abs(parmpaths$std.all - physpaths$std.all)>diff_cutoff
parmpaths$diff = abs(parmpaths$std.all - physpaths$std.all)>diff_cutoff

physpaths$sigcat = cut(abs(physpaths$std.all), c(0,.2,.4,.6,.8,1))
parmpaths$sigcat = cut(abs(parmpaths$std.all), c(0,.2,.4,.6,.8,1))

# Find all variables in the model
allvars = unique(c(parmpaths$rhs,parmpaths$lhs))

# Subset paths to only those that are significant
parmpaths = parmpaths[which(apply(parmpaths[,c('std.ci.lower','std.ci.upper')],1,prod)>0),]
physpaths = physpaths[which(apply(physpaths[,c('std.ci.lower','std.ci.upper')],1,prod)>0),]

subset(parmpaths, diff)
subset(physpaths, diff)

## Make tables of paths by which regional climate variables affect local richness

# Table for regTorich model which includes direct effect
indirect_ir = phys_ir
direct = phys_d

order_clim = c('wetness','rain_lowRH','pseas','mat','iso','radiation')
order_clim_reg = c(paste(order_clim[1:5],'reg_mean', sep='_'), paste(order_clim[1:5], 'reg_var', sep='_'), 'regS_tree')

climReg_IE_tab = data.frame(var = order_clim_reg, predictor = varnames[order_clim_reg,'displayName'])
use_dir = direct[order_clim_reg,]
climReg_IE_tab$dir = use_dir$std.all
use_loc = subset(indirect_ir, IEvar1 %in% c('clim_loc','PIE.ba.tree')); rownames(use_loc) = use_loc$predictor
climReg_IE_tab$loc = use_loc[order_clim_reg, 'std.all']
use_regS = subset(indirect_ir, IEvar1=='regS'); rownames(use_regS) = use_regS$predictor
climReg_IE_tab$regS = use_regS[order_clim_reg,'std.all']
use_FH = subset(indirect_ir, IEvar1=='regS_tree'&is.na(IEvar2)); rownames(use_FH) = use_FH$predictor
climReg_IE_tab$FH = use_FH[order_clim_reg, 'std.all']

climReg_IE_tab$dir_sig = apply(use_dir[,c('std.ci.lower','std.ci.upper')], 1, prod)>0
climReg_IE_tab$loc_sig = apply(use_loc[order_clim_reg,c('std.ci.lower','std.ci.upper')], 1, prod)>0
climReg_IE_tab$regS_sig = apply(use_regS[order_clim_reg,c('std.ci.lower','std.ci.upper')], 1, prod)>0
climReg_IE_tab$FH_sig = apply(use_FH[order_clim_reg,c('std.ci.lower','std.ci.upper')], 1, prod)>0

sub_vars_FH = subset(indirect_ir, (IEvar1=='regS_tree')&!is.na(IEvar2)); colnames(sub_vars_FH)[1] = 'var'; colnames(sub_vars_FH)[colnames(sub_vars_FH)=='std.all'] <- 'subpath_FH'
sub_vars_FH$subFH_sig = apply(sub_vars_FH[,c('std.ci.lower','std.ci.upper')], 1, prod)>0

climReg_IE_tab = merge(climReg_IE_tab, sub_vars_FH[,c('var','IEvar2','subpath_FH','subFH_sig')])
climReg_IE_tab = climReg_IE_tab = climReg_IE_tab[,c('var','IEvar2','predictor','dir','loc','regS','FH','subpath_FH',
	'dir_sig','loc_sig','regS_sig','FH_sig','subFH_sig')]

write.csv(climReg_IE_tab, './Parm-Phys/Indirect regional climate effects phys.csv', row.names=F)

## Compare effects of climate variables when measured locally or regionally

total = phys # This changes based on what response variable is being analyzed
total = rbind(total, phys_d[c('regS','tot_abun_log'),1:4])

cvars = c('wetness','rain_lowRH','pseas','iso','mat')

lvars = c(cvars, 'PIE.ba.tree')
rvars = c(paste(cvars, 'reg_mean', sep='_'), 'regS_tree')

svg('./Parm-Phys/compare local-regional total effects phys.svg', height=5.5, width=5.5)
par(mar=c(5,5,1,1))
par(lwd=2)
par(lend=1)
par(cex.lab=1.2)
plot(total[lvars, 'std.all'], total[rvars, 'std.all'], 	
	xlab='Local Effect', ylab='Regional Effect', type='n', las=1,
	xlim=c(-0.9,0.9), ylim=c(-0.9,0.9))
usr=par('usr')
polygon(usr[c(1,1,2,2)],usr[c(3,4,4,3)], col='#8dff94fb')
polygon(c(usr[1],0,usr[1],usr[2],0,usr[2],usr[1]),
	c(usr[3],0,usr[4],usr[4],0,usr[3],usr[3]), col='#b3b5ffff')
abline(h=0,v=0,lty=2, col='black')
arrows(total[lvars,'std.ci.lower'], total[rvars,'std.all'],
	total[lvars,'std.ci.upper'], total[rvars,'std.all'], col='black', code=3, angle=90, length=0.05)
arrows(total[lvars,'std.all'], total[rvars,'std.ci.lower'],
	total[lvars,'std.all'], total[rvars,'std.ci.upper'], col='black', code=3, angle=90, length=0.05)
points(total[lvars, 'std.all'], total[rvars, 'std.all'], pch=15, col='black')

# Labels
text(total[lvars,'std.all'], total[rvars,'std.all'], varnames[lvars,'midName'],
	pos=c(2,2,4,4,2,4), offset=1.5, col='black')

dev.off()









#############################################################
### Variation Partitioning

library(hier.part) # combos
library(MuMIn) # r.squaredLR
library(MASS) # glm.nb

# Define dataset- use trans_data where variables have been log or sqrt transformed but not scaled
use_data = trans_data[,colnames(trans_data) %in% rownames(predtypes)]

# Define predictor variable and appropriate abundance and regional richness variables
# This changes based on the analysis of Parm vs Phys
use_response = 'Physciaceae' 
use_reg = 'regPhys'
use_abun = 'phys_abun_log'
use_data$richness = trans_data[,use_response]
use_data$reg = trans_data[,use_reg]
use_data$abun_log = trans_data[,use_abun]

# Testing and training data sets
# Initially use fitplots to determine whether to use quadratic relationships.
# Then, re-run final analysis using only use testing data set
use_data_test = use_data[testplots$yrplot.id,]
use_data_fit = use_data[fitplots$yrplot.id,]

## Test linear, log-linear, Poisson, and Negative Binomial GLMs with different link functions

# Define set of variables to consider in models
use_vars = rownames(subset(predtypes, !(label %in% c('','a1','r1','R1'))))
use_vars = use_vars[-grep('soil',use_vars)] # Leave out soil variables

pois_log_mod = glm(richness~., family=poisson(link='log'), data=use_data_fit[,c('richness',use_vars)])
nb_log_mod = glm.nb(richness~., link='log', data=use_data_fit[,c('richness',use_vars)])
#gaus_log_mod = glm(richness~., family=gaussian(link='log'), data=use_data_fit[,c('richness',use_vars)])
gaus_iden_mod = glm(richness~., family=gaussian(link='identity'), data=use_data_fit[,c('richness',use_vars)])

AIC(gaus_iden_mod, pois_log_mod, nb_log_mod) #gaus_log_mod,
# nb_log mod wins by far Parmeliaceae and Physciaceae

# Test models for richness ~ abundance
abun_pois_log = glm(richness~abun_log, family=poisson(link='log'), data=use_data_fit)
abun_nb_log = glm.nb(richness~abun_log, link='log', data=use_data_fit) 
#abun_gaus_log = glm(richness~abun_log, family=gaussian(link='log'), data=use_data_fit)
abun_gaus_iden = glm(richness~abun_log, family=gaussian(link='identity'), data=use_data_fit)

AIC(abun_pois_log, abun_nb_log, abun_gaus_iden) 
# nb_log wins for Parmeliaceae and is equivalent to pois_log for Phys so will use nb for both

# calculate correlation coef for abundance
abun_mod=glm.nb(richness~abun_log, link='log', data=use_data_test)
summary(abun_mod)
r.squaredLR(abun_mod)

## Test gaussian, Poisson, NB for regS
pois_log = glm(richness~reg, family=poisson(link='log'), data=use_data_fit)
pois_iden = glm(richness~reg, family=poisson(link='identity'), data=use_data_fit)
nb_log = glm.nb(richness~reg, link='log', data=use_data_fit)
nb_iden = glm.nb(richness~reg, link='identity', data=use_data_fit)
#gaus_log = glm(richness~reg, family=gaussian(link='log'), data=use_data_fit)
gaus_iden = glm(richness~reg, family=gaussian(link='identity'), data=use_data_fit)

AIC(pois_log, pois_iden, nb_log, nb_iden, gaus_iden)
# nb_log wins for Parmeliaceae and is almost equivalent to nb_iden for Phys so will use nb_log for both


## Univariate models

# For lichen richness - use_vars from above
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

# Save models - names changes based on data set
write.csv(unimods, './Parm-Phys/univariate_models_Phys.csv', row.names=T)

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

# Save models - names changes based on data set
write.csv(unimods_reg, './Parm-Phys/univariate_models_regS_Phys.csv', row.names=T)

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

# Save model comparison - names changes based on data set
write.csv(modcompare, './Parm-Phys/Univariate model shapes GLM-NB Phys.csv', row.names=T)
write.csv(modcompare_reg, './Parm-Phys/Univariate model shapes of regPhys.csv', row.names=T)

# Which variables have AIC supported concave-down relationships?
# Make sure to only run the code that corresponds to the correct data set
sq_vars_parm = rownames(subset(modcompare, concavity=='down'&type=='quadratic'))
sq_vars_reg_parm = rownames(subset(modcompare_reg, concavity=='down'&type=='quadratic'))
sq_vars_phys = rownames(subset(modcompare, concavity=='down'&type=='quadratic'))
sq_vars_reg_phys = rownames(subset(modcompare_reg, concavity=='down'&type=='quadratic'))

save(sq_vars_parm, sq_vars_phys, sq_vars_reg_parm, sq_vars_reg_phys, file='./Parm-Phys/quadratic variables parm-phys.RData')

# Compare quadratic relationships for Parm vs Phys
# The are pretty differect so we will use different sets of quadratic variables in the models for partitioning
sq_vars_parm[!(sq_vars_parm %in% sq_vars_phys)]
sq_vars_phys[!(sq_vars_phys %in% sq_vars_parm)]

sq_vars_reg_parm[!(sq_vars_reg_parm %in% sq_vars_reg_phys)]
sq_vars_reg_phys[!(sq_vars_reg_phys %in% sq_vars_reg_parm)]


### Variation Partitioning

# Get names of variables that should have quadratic relationships in models
load('./Parm-Phys/quadratic variables parm-phys.RData')

# Define dataset- use trans_data where variables have been log or sqrt transformed but not scaled
use_data = trans_data[, colnames(trans_data) %in% rownames(predtypes)]

# Define predictor variable and appropriate abundance and regional richness variables
# This changes based on the analysis of Parm vs Phys
use_response = 'Parmeliaceae' 
use_reg = 'regParm'
use_abun = 'parm_abun_log'
use_data$richness = trans_data[,use_response]
use_data$reg = trans_data[,use_reg]
use_data$abun_log = trans_data[,use_abun]

# Use quadratic relationships for concave-down AIC supported models in variation partitioning
# Use the correct squared terms depending on the data set
sq_vars = sq_vars_parm
sq_vars_reg = sq_vars_reg_parm

sq_df = use_data[,unique(c(sq_vars, sq_vars_reg))]^2
colnames(sq_df) = paste(colnames(sq_df),'2', sep='')

# Add squared terms
use_data_sq = cbind(use_data, sq_df)

# Use test plots for model results:
use_data_test = use_data_sq[testplots$yrplot.id,]

# Predictors to be used in models
use_preds = subset(predtypes, !(label %in% c('','r1','a1','p1','P1')))

## Local-Regional Models
regionvars = rownames(use_preds)[use_preds$scale=='regional']
localvars = rownames(use_preds)[use_preds$scale=='local']
localvars = localvars[-grep('soil', localvars)] # Leave out soil vars

region_mod = glm.nb(richness~., data=use_data_test[,c('richness', regionvars, paste(sq_vars[sq_vars %in% regionvars],2,sep=''))])
local_mod = glm.nb(richness~., data=use_data_test[,c('richness', localvars, paste(sq_vars[sq_vars %in% localvars],2,sep=''))])

AIC(region_mod, local_mod) # local_mod is best for Phys but region_mod is best for Parm

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

# Save partitioning results - this changes based on the data set
local_regional_partition_parm = partvar2(Rs)
local_regional_partition_phys = partvar2(Rs)

## Heterogeneity - Optimality

## At regional scale
RHvars = rownames(subset(use_preds, scale=='regional'&mode=='het'))
ROvars = rownames(subset(use_preds, scale=='regional'&mode=='opt'))

RH_mod = lm(reg~., data=use_data_test[,c('reg', RHvars, paste(sq_vars_reg[sq_vars_reg %in% RHvars],2,sep=''))])
RO_mod = lm(reg~., data=use_data_test[,c('reg', ROvars, paste(sq_vars_reg[sq_vars_reg %in% ROvars],2,sep=''))])
#RO_mod = lm(reg~., data=use_data_test[,c('reg', ROvars)]) # Don't need to add sq_vars_reg for Phys

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

# Save partitioning results - this changes based on the data set
regional_het_opt_partition_parm = partvar2(Rs)
regional_het_opt_partition_phys = partvar2(Rs)

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

# Save partitioning results - this changes based on the data set
local_het_opt_partition_parm = partvar2(Rs)
local_het_opt_partition_phys = partvar2(Rs)

# Save partitioning results for later
save(local_regional_partition_parm, local_regional_partition_phys,
	regional_het_opt_partition_parm, local_het_opt_partition_phys, 
	local_het_opt_partition_parm, local_het_opt_partition_phys,
	file='./Parm-Phys/variation partitioning results parm-phys.RData')



## Plot variation partioning analysis in three panels
load('./Parm-Phys/variation partitioning results parm-phys.RData')

barwide=1.5
mycols = matrix(c('#b3b5ff','#6b6dd7','#8dff94','#38af4f'), nrow=2)
colnames(mycols) = c('regional','local')
rownames(mycols) = c('het','opt')
blend_shade = '99'
blend_col = '#126A71'

use_lr = local_regional_partition_phys
use_ho_loc = local_het_opt_partition_phys
use_ho_reg = regional_het_opt_partition_phys

svg('./Parm-Phys/variation partitioning figure phys.svg', height=5, width=10)
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
	rect(0,0,barwide, sum(use_lr[1]), lwd=3, col=mycols['opt','regional'])
	# Add rectangle for second component
	rect(0,sum(use_lr[1:2]), barwide,sum(use_lr[1:3]), lwd=3, col=mycols['opt','local'])
	# Add rectangle for overlap
	rect(0,use_lr[1],barwide,sum(use_lr[1:2]), lwd=3, col=blend_col)

	## Regional Heterogeneity-Optimality
	# Add rectangle for first component
	rect(2,0,2+barwide, sum(use_ho_reg[1]), lwd=3, col=mycols['het','regional'])
	# Add rectangle for second component
	rect(2,sum(use_ho_reg[1:2]), 2+barwide,sum(use_ho_reg[1:3]), lwd=3, col=paste(mycols['opt','regional'], blend_shade, sep=''))
	# Add rectangle for overlap
	rect(2,use_ho_reg[1], 2+barwide,sum(use_ho_reg[1:2]), lwd=3, col=mycols['opt','regional'])

	## Local Heterogeneity-Optimality
	# Add rectangle for first component
	rect(4,0,4+barwide, sum(use_ho_loc[1]), lwd=3, col=mycols['het','local'])
	# Add rectangle for second component
	rect(4,sum(use_ho_loc[1:2]), 4+barwide,sum(use_ho_loc[1:3]), lwd=3, col=paste(mycols['opt','local'], blend_shade, sep=''))
	# Add rectangle for overlap
	rect(4,use_ho_loc[1], 4+barwide,sum(use_ho_loc[1:2]), lwd=3, col=mycols['opt','local'])

	## Add axis
	axis(2, las=1, cex.axis=2, lwd=3)
	mtext('Variation Explained', 2, 4, cex=2)
	
	# Add partition labels
	lablocs = sapply(1:4, function(x) sum(use_lr[0:(x-1)])+use_lr[x]/2)
	text(barwide/2, lablocs, labels=names(use_lr)[1:4], cex=2)
	lablocs = sapply(1:4, function(x) sum(use_ho_reg[0:(x-1)])+use_ho_reg[x]/2)
	text(2+barwide/2, lablocs, labels=names(use_ho_reg)[1:4], cex=2)
	lablocs = sapply(1:4, function(x) sum(use_ho_loc[0:(x-1)])+use_ho_loc[x]/2)
	text(4+barwide/2, lablocs, labels=names(use_ho_loc)[1:4], cex=2)
	
	# Add panel text A, B, C
	par(xpd=T)
	text(c(0,2,4),1.05, c('A','B','C'), cex=2, adj=c(0,0))
	par(xpd=F)
dev.off()



##########################################################################
### Compare traits
library(stringr) # str_trim
# Read in trait data for FIA species
traits = read.csv('./Data/LichenTraits/FIA_lichen_species_traits_LIAS.csv')
spXtrait = read.csv('./Data/LichenTraits/spXtrait_matrix_foranalysis.csv', row.names=1)

# Read in a list of species that identifies the family they are in
splist = read.csv('lichen_species.csv')
splist$taxon = str_trim(paste(splist$GENUS, splist$SPECIES))

# Make trait matrices for Parm and Physc separately
traitmat_parm = spXtrait[rownames(spXtrait) %in% subset(splist, family=='Parmeliaceae')$taxon,]
traitmat_phys = spXtrait[rownames(spXtrait) %in% subset(splist, family=='Physciaceae')$taxon,]

## Calculate FD within each group




## Calculate trait means within each group

summary(traitmat_parm)







