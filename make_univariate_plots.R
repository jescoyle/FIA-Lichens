## This script makes univariate plots of lichen richness versus environmental variables


source('./UNC/Projects/FIA Lichen/GitHub/FIA-Lichens/load_data.R')
source('./GitHub/FIA-Lichens/fia_lichen_analysis_functions.R')

library(MuMIn) # r.squaredLR
library(MASS) # glm.nb

##############################################################################
### Make Univariate plots

# Read in a list of state codes categorized by regions
regions = read.csv('fia_lichen_regions.csv', row.names=1)

# highlight by region plots
usereg = factor(regions[master[rownames(trans_data),'state.abbr'],'big_region'])
usereg = factor(regions[master[rownames(trans_data),'state.abbr'],'sm_region'])

# Define colors and symbols for different regions
regcol = c('cornflowerblue','darkred','orange','forestgreen','purple','cyan3')
regpch = c(0:5)

# Define variables to be plotted
use_vars = c('regS','wetness','PIE.ba.tree','regS_tree','wetness_reg_mean')

# Define plotting ranges for each of the five variables to be plotted
myranges = data.frame(lower = c(120, -3.5, 0, 0, -3), 
	upper = c(210, 5.5, 1, 170, 2), 
	use.sq = c(F, T, F, T, T), row.names=use_vars)

# Using trans_data for plots b/c has all plots, no outliers, and log/sqrt transformed where necessary
# Should models only be fit to testing dataset even if all plots are shown?
for(i in use_vars){
	
	png(paste('./Figures/New Coordinates/highlight big region ',i,' vs lichen richness.png', sep=''), height=1500, width=1500)
		
		par(mar=c(10,11,4,1))
		
		# Set up plot
		plot(lichen.rich~trans_data[,i], data=trans_data, type='n', las=1, ylim=c(0,40), xlim=as.numeric(myranges[i,c('lower','upper')]),
			xlab='', ylab='',	axes=F
		)
	
		# Add points
		points(lichen.rich~trans_data[,i], data=trans_data, col=regcol[usereg], pch=regpch[usereg], cex=4, lwd=3)
		
		# Add axes and labels
		axis(1, cex.axis=4, las=1, padj=1, lwd=2)
		axis(2, cex.axis=4, las=1, lwd=2, at=seq(0,40,5))
		mtext(varnames[i,'displayName'],1,7, cex=4)
		mtext('Local lichen species richness',2,8, cex=4)

		# Add model
		if(myranges[i,'use.sq']){
			this_data = data.frame(trans_data[,c('lichen.rich',i)], trans_data[,i]^2)
			this_mod = glm.nb(lichen.rich~., data=this_data, link='log')
			use_coef = coef(this_mod)
			use_x = seq(min(trans_data[,i]), max(trans_data[,i]), length.out=100)
			use_y = exp(use_coef[1] + use_coef[2]*use_x + use_coef[3]*(use_x^2))
		} else {
			this_mod = glm.nb(lichen.rich~., data=trans_data[,c('lichen.rich',i)], link='log')
			use_coef = coef(this_mod)
			use_x = seq(min(trans_data[,i]), max(trans_data[,i]), length.out=100)
			use_y = exp(use_coef[1]+use_coef[2]*use_x)
		}
		lines(use_x, use_y, lwd=4, col='black')
		r2 = r.squaredLR(this_mod, null=glm.nb(lichen.rich~1, data=trans_data, link='log'))
		
		use_label = bquote(italic(R)^2 == .(format(r2, nsmall=2, digits=2)))
		mtext(use_label, 3, adj=0, line=0, cex=4)
		
	dev.off()
}

## Abundance vs Richness

svg('./Figures/lichen richness vs abundance.svg', height=8, width=8)
par(mar=c(10,11,4,1))
		
# Set up plot
plot(log10(lichen.rich)~tot_abun_log, data=trans_data, type='n', las=1, xlim=c(0,7),
	xlab='', ylab='',	axes=F#, log='xy'
)

# Add points
points(log10(lichen.rich)~tot_abun_log, data=trans_data, col=regcol[usereg], pch=regpch[usereg], cex=4, lwd=3)
		
# Add axes and labels
axis(1, cex.axis=4, las=1, padj=1, lwd=2)
axis(2, cex.axis=4, las=1, lwd=2)
mtext('Log Abundance',1,7, cex=4)
mtext('Log Richness',2,8, cex=4)
dev.off()

### Black and White figures for 4 or 6-panel figure
### Does not automatically fit best model- need to know ahead of time.
trans_data_test = trans_data[testplots$yrplot.id,]

svg('./Figures/univariate models 4-panel.svg', height=7, width=8)

par(mfrow=c(2,2)) #par(mfrow=c(3,2))
par(mar=c(4,4,1.5,1))
par(mgp=c(2.4,0.7,0))
par(cex.axis=1.2)
par(cex.lab=1.2)
par(pch=1)
par(cex=1)
par(las=1)

# Abundance
# Both richness and abundance are log(x+1) transformed
plot(log(lichen.rich+1)~tot_abun_log, data=trans_data_test,
	xlab='Ln( Abundance + 1 )', ylab='Ln( Species Richness + 1 )', 
	las=1, xlim=c(0,7), col='#00000050', lwd=2)		
use_mod = glm.nb(lichen.rich~tot_abun_log, data=trans_data_test, link='log')
use_coef = coef(use_mod)
use_min = max(min(trans_data_test$tot_abun_log),-use_coef[1]/use_coef[2])
use_x = seq(use_min, max(trans_data_test$tot_abun_log), length.out=100)
use_y = exp(use_coef[1]+use_coef[2]*use_x)
lines(use_x, log(use_y+1), lwd=5, col='white')
lines(use_x, log(use_y+1), lwd=4, col='black')
r2 = r.squaredLR(use_mod, null=glm.nb(lichen.rich~1, data=trans_data_test, link='log'))
r2 = attr(r2, 'adj.r.squared')
use_label = bquote(italic(R)^2 == .(format(r2, nsmall=2, digits=2)))
mtext(use_label, 3, adj=1, line=0)
mtext('A',3,adj=0,line=0, font=2, cex=2)

# Regional richness
plot(lichen.rich~regS, data=trans_data_test,
	xlab='Regional Species Richness', ylab='Local Species Richness', 
	las=1, xlim=c(125,205), ylim=c(0,40),
	col='#00000050', lwd=2, axes=F)
axis(1, at=seq(125,200,25))
axis(2)
box()
use_mod = glm.nb(lichen.rich~regS, data=trans_data_test, link='log')
use_coef = coef(use_mod)
use_x = seq(min(trans_data_test$regS), max(trans_data_test$regS), length.out=100)
use_y = exp(use_coef[1]+use_coef[2]*use_x)
lines(use_x, use_y, lwd=5, col='white')
lines(use_x, use_y, lwd=4, col='black')
r2 = r.squaredLR(use_mod, null=glm.nb(lichen.rich~1, data=trans_data_test, link='log'))
r2 = attr(r2, 'adj.r.squared')
use_label = bquote(italic(R)^2 == .(format(r2, nsmall=2, digits=2)))
mtext(use_label, 3, adj=1, line=0)
mtext('B',3,adj=0,line=0, font=2, cex=2)

# Wetness
plot(lichen.rich~wetness, data=trans_data_test,
	xlab='Wet Climate (local)', ylab='Local Species Richness', 
	las=1, xlim=as.numeric(myranges['wetness',c('lower','upper')]), ylim=c(0,40),
	col='#00000050', lwd=2, axes=F)
axis(1, at=seq(-3,5,2))
axis(1, at=seq(-2,4,2))
axis(2)
box()
use_mod = glm.nb(lichen.rich~wetness+I(wetness^2), data=trans_data_test, link='log')
use_coef = coef(use_mod)
use_x = seq(min(trans_data_test$wetness), max(trans_data_test$wetness), length.out=100)
use_y = exp(use_coef[1]+use_coef[2]*use_x+use_coef[3]*(use_x^2))
lines(use_x, use_y, lwd=5, col='white')
lines(use_x, use_y, lwd=4, col='black')
r2 = r.squaredLR(use_mod, null=glm.nb(lichen.rich~1, data=trans_data_test, link='log'))
r2 = attr(r2, 'adj.r.squared')
use_label = bquote(italic(R)^2 == .(format(r2, nsmall=2, digits=2)))
mtext(use_label, 3, adj=1, line=0)
mtext('C',3,adj=0,line=0, font=2, cex=2)

# Regional Mean Wetness
plot(lichen.rich~wetness_reg_mean, data=trans_data_test,
	xlab='Wet Climate (regional)', ylab='Local Species Richness', 
	las=1, xlim=as.numeric(myranges['wetness_reg_mean',c('lower','upper')]), ylim=c(0,40),
	col='#00000050', lwd=2, axes=F)
axis(1)
axis(2)
box()
use_mod = glm.nb(lichen.rich~wetness_reg_mean+I(wetness_reg_mean^2), data=trans_data_test, link='log')
use_coef = coef(use_mod)
use_x = seq(min(trans_data_test$wetness_reg_mean), max(trans_data_test$wetness_reg_mean), length.out=100)
use_y = exp(use_coef[1]+use_coef[2]*use_x+use_coef[3]*(use_x^2))
lines(use_x, use_y, lwd=5, col='white')
lines(use_x, use_y, lwd=4, col='black')
r2 = r.squaredLR(use_mod, null=glm.nb(lichen.rich~1, data=trans_data_test, link='log'))
use_label = bquote(italic(R)^2 == .(format(r2, nsmall=2, digits=2)))
mtext(use_label, 3, adj=1, line=0)
mtext('D',3,adj=0,line=0, font=2, cex=2)

# Tree diversity
#plot(lichen.rich~PIE.ba.tree, data=trans_data_test,
#	xlab='Tree Diversity (local)', ylab='Local Species Richness', 
#	las=1, xlim=as.numeric(myranges['PIE.ba.tree',c('lower','upper')]), ylim=c(0,40),
#	col='#00000050', lwd=2)
#use_mod = glm.nb(lichen.rich~PIE.ba.tree, data=trans_data_test, link='log')
#use_coef = coef(use_mod)
#use_x = seq(min(trans_data_test$PIE.ba.tree), max(trans_data_test$PIE.ba.tree), length.out=100)
#use_y = exp(use_coef[1]+use_coef[2]*use_x)
#lines(use_x, use_y, lwd=5, col='white')
#lines(use_x, use_y, lwd=4, col='black')
#r2 = r.squaredLR(use_mod, null=glm.nb(lichen.rich~1, data=trans_data_test, link='log'))
#r2 = attr(r2, 'adj.r.squared')
#use_label = bquote(italic(R)^2 == .(format(r2, nsmall=2, digits=2)))
#mtext(use_label, 3, adj=1, line=0)
#mtext('E',3,adj=0,line=0, font=2, cex=2)

# Regional Tree Diversity
#plot(regS~regS_tree, data=trans_data_test,
#	xlab='Tree Richness (regional)', ylab='Regional Species Richness', 
#	las=1, xlim=as.numeric(myranges['regS_tree',c('lower','upper')]), ylim=c(120,220),
#	col='#00000050', lwd=2, axes=F)
#axis(1, at=seq(0,160,40))
#axis(2)
#box()
#use_mod = lm(regS~regS_tree+I(regS_tree^2), data=trans_data_test)
#use_coef = coef(use_mod)
#use_x = seq(min(trans_data_test$regS_tree), max(trans_data_test$regS_tree), length.out=100)
#use_y = use_coef[1]+use_coef[2]*use_x+use_coef[3]*(use_x^2)
#lines(use_x, use_y, lwd=5, col='white')
#lines(use_x, use_y, lwd=4, col='black')
#r2 = summary(use_mod)$adj.r.squared
#use_label = bquote(italic(R)^2 == .(format(r2, nsmall=2, digits=2)))
#mtext(use_label, 3, adj=1, line=0)
#mtext('F',3,adj=0,line=0, font=2, cex=2)

dev.off()

## Two panels- richness (with regional heterogeneity interaction) and abundance
reghet_vars = rownames(subset(predtypes, scale=='regional'&mode=='het'&label!=''&type=='env'))
reghet_pca  = prcomp(use_data[,reghet_vars], center=T, scale=T)
reghet_pc1 = predict(reghet_pca)[,'PC1']
reghet_pc1 = reghet_pc1[rownames(trans_data_test)]
lowhet = quantile(reghet_pc1, 0.25)
highhet = quantile(reghet_pc1, 0.75)
medhet = median(reghet_pc1)


svg('./Figures/univariate models 2-panel.svg', height=7, width=5)

par(mfrow=c(2,1))
par(mar=c(4,4,1.5,6))
par(mgp=c(2.4,0.7,0))
par(cex.axis=1.2)
par(cex.lab=1.2)
par(pch=1)
par(cex=1)
par(las=1)

# Abundance
# Both richness and abundance are log(x+1) transformed
plot(log(lichen.rich+1)~tot_abun_log, data=trans_data_test,
	xlab='Ln( Abundance + 1 )', ylab='Ln( Species Richness + 1 )', 
	las=1, xlim=c(0,7), col='#00000050', lwd=2, pch=1)		
use_mod = glm.nb(lichen.rich~tot_abun_log, data=trans_data_test, link='log')
use_coef = coef(use_mod)
use_min = max(min(trans_data_test$tot_abun_log),-use_coef[1]/use_coef[2])
use_x = seq(use_min, max(trans_data_test$tot_abun_log), length.out=100)
use_y = exp(use_coef[1]+use_coef[2]*use_x)
lines(use_x, log(use_y+1), lwd=5, col='white')
lines(use_x, log(use_y+1), lwd=4, col='black')
r2 = r.squaredLR(use_mod, null=glm.nb(lichen.rich~1, data=trans_data_test, link='log'))
r2 = attr(r2, 'adj.r.squared')
use_label = bquote(italic(R)^2 == .(format(r2, nsmall=2, digits=2)))
mtext(use_label, 3, adj=1, line=0)
mtext('A',3,adj=0,line=0, font=2, cex=2)


# Regional richness with heterogeneity interaction
#par(lend=1)
#plot(lichen.rich~regS, data=trans_data_test,
#	xlab='Regional Species Richness', ylab='Local Species Richness', 
#	las=1, xlim=c(125,205), ylim=c(0,40),
#	col='#00000050', lwd=2, axes=F)
#axis(1, at=seq(125,200,25))
#axis(2)
#box()
#use_mod = glm.nb(lichen.rich~regS*reghet_pc1, data=trans_data_test, link='log')
#use_coef = coef(use_mod)
#use_x = seq(min(trans_data_test$regS), max(trans_data_test$regS), length.out=100)
#use_ylow = exp(use_coef[1]+use_coef[2]*use_x+use_coef[3]*lowhet+use_coef[4]*lowhet*use_x) 
#use_yhigh = exp(use_coef[1]+use_coef[2]*use_x+use_coef[3]*highhet+use_coef[4]*highhet*use_x) 
#use_ymed = exp(use_coef[1]+use_coef[2]*use_x+use_coef[3]*0+use_coef[4]*0*use_x) 
#lines(use_x, use_ylow, lwd=5, col='white')
#lines(use_x, use_ymed, lwd=5, col='white')
#lines(use_x, use_yhigh, lwd=5, col='black')
#lines(use_x, use_ylow, lwd=4, col='black', lty=2)
#lines(use_x, use_ymed, lwd=4, col='black')
#lines(use_x, use_yhigh, lwd=4, col='white', lty=3)

# Regional richness with heterogeneity colored
mycolbw = c('grey80','black')
colorvecbw = colorRampPalette(mycolbw)(100)[cut(reghet_pc1, 100, include.lowest=T)]
plot(lichen.rich~regS, data=trans_data_test,
	xlab='Regional Species Richness', ylab='Local Species Richness', 
	las=1, xlim=c(125,205), ylim=c(0,40),
	col=colorvecbw, lwd=2, axes=F)
axis(1, at=seq(125,200,25))
axis(2)
box()
use_mod = glm.nb(lichen.rich~regS, data=trans_data_test, link='log')
use_coef = coef(use_mod)
use_x = seq(min(trans_data_test$regS), max(trans_data_test$regS), length.out=100)
use_y = exp(use_coef[1]+use_coef[2]*use_x) 
lines(use_x, use_y, lwd=4, col='black')

usr = par('usr')
plotColorRamp(cols = mycolbw, n = 100, barends = c(usr[2], usr[3], usr[2]+0.05*diff(usr[1:2]), usr[4]),
	labels = seq(-2.5,3.5,.5), uneven.lab=T, labrange=range(reghet_pc1), title='Regional Heterogeneity (PC1)',
	mycex=1)


r2 = r.squaredLR(glm.nb(lichen.rich~regS, data=trans_data_test, link='log'), null=glm.nb(lichen.rich~1, data=trans_data_test, link='log'))
r2 = attr(r2, 'adj.r.squared')
use_label = bquote(italic(R)^2 == .(format(r2, nsmall=2, digits=2)))
mtext(use_label, 3, adj=1, line=0)
mtext('B',3,adj=0,line=0, font=2, cex=2)

dev.off()




## With local lichen richness as response for all variables
svg('./Figures/univariate models 6-panel.svg', height=9, width=8)

par(mfrow=c(3,2))
par(mar=c(4,4,1.5,1))
par(mgp=c(2.4,0.7,0))
par(cex.axis=1.2)
par(cex.lab=1.2)
par(pch=1)
par(cex=1)
par(las=1)

# Abundance

plot(log10(lichen.rich)~tot_abun_log, data=trans_data_test,
	xlab='Log Abundance', ylab='Log Species Richness', 
	las=1, xlim=c(0,7), ylim=c(-.1,1.6),
	col='#00000050', lwd=2)		
use_mod = glm.nb(lichen.rich~tot_abun_log, data=trans_data_test, link='log')
use_coef = coef(use_mod)
use_min = max(min(trans_data_test$tot_abun_log),-use_coef[1]/use_coef[2])
use_x = seq(use_min, max(trans_data_test$tot_abun_log), length.out=100)
use_y = exp(use_coef[1]+use_coef[2]*use_x)
lines(use_x, log10(use_y), lwd=5, col='white')
lines(use_x, log10(use_y), lwd=4, col='black')
r2 = r.squaredLR(use_mod, null=glm.nb(lichen.rich~1, data=trans_data_test, link='log'))
r2 = attr(r2, 'adj.r.squared')
use_label = bquote(italic(R)^2 == .(format(r2, nsmall=2, digits=2)))
mtext(use_label, 3, adj=1, line=0)
mtext('A',3,adj=0,line=0, font=2, cex=2)


# Regional richness
plot(lichen.rich~regS, data=trans_data_test,
	xlab='Regional Species Richness', ylab='Local Species Richness', 
	las=1, xlim=c(125,205), ylim=c(0,40),
	col='#00000050', lwd=2, axes=F)
axis(1, at=seq(125,200,25))
axis(2)
box()
use_mod = glm.nb(lichen.rich~regS, data=trans_data_test, link='log')
use_coef = coef(use_mod)
use_x = seq(min(trans_data_test$regS), max(trans_data_test$regS), length.out=100)
use_y = exp(use_coef[1]+use_coef[2]*use_x)
lines(use_x, use_y, lwd=5, col='white')
lines(use_x, use_y, lwd=4, col='black')
r2 = r.squaredLR(use_mod, null=glm.nb(lichen.rich~1, data=trans_data_test, link='log'))
r2 = attr(r2, 'adj.r.squared')
use_label = bquote(italic(R)^2 == .(format(r2, nsmall=2, digits=2)))
mtext(use_label, 3, adj=1, line=0)
mtext('B',3,adj=0,line=0, font=2, cex=2)

# Wetness
plot(lichen.rich~wetness, data=trans_data_test,
	xlab='Wet Climate (local)', ylab='Local Species Richness', 
	las=1, xlim=as.numeric(myranges['wetness',c('lower','upper')]), ylim=c(0,40),
	col='#00000050', lwd=2, axes=F)
axis(1, at=seq(-3,5,2))
axis(1, at=seq(-2,4,2))
axis(2)
box()
use_mod = glm.nb(lichen.rich~wetness+I(wetness^2), data=trans_data_test, link='log')
use_coef = coef(use_mod)
use_x = seq(min(trans_data_test$wetness), max(trans_data_test$wetness), length.out=100)
use_y = exp(use_coef[1]+use_coef[2]*use_x+use_coef[3]*(use_x^2))
lines(use_x, use_y, lwd=5, col='white')
lines(use_x, use_y, lwd=4, col='black')
r2 = r.squaredLR(use_mod, null=glm.nb(lichen.rich~1, data=trans_data_test, link='log'))
r2 = attr(r2, 'adj.r.squared')
use_label = bquote(italic(R)^2 == .(format(r2, nsmall=2, digits=2)))
mtext(use_label, 3, adj=1, line=0)
mtext('C',3,adj=0,line=0, font=2, cex=2)

# Regional Mean Wetness
plot(lichen.rich~wetness_reg_mean, data=trans_data_test,
	xlab='Wet Climate (regional)', ylab='Local Species Richness', 
	las=1, xlim=as.numeric(myranges['wetness_reg_mean',c('lower','upper')]), ylim=c(0,40),
	col='#00000050', lwd=2, axes=F)
axis(1)
axis(2)
box()
use_mod = glm.nb(lichen.rich~wetness_reg_mean+I(wetness_reg_mean^2), data=trans_data_test, link='log')
use_coef = coef(use_mod)
use_x = seq(min(trans_data_test$wetness_reg_mean), max(trans_data_test$wetness_reg_mean), length.out=100)
use_y = exp(use_coef[1]+use_coef[2]*use_x+use_coef[3]*(use_x^2))
lines(use_x, use_y, lwd=5, col='white')
lines(use_x, use_y, lwd=4, col='black')
r2 = r.squaredLR(use_mod, null=glm.nb(lichen.rich~1, data=trans_data_test, link='log'))
r2 = attr(r2, 'adj.r.squared')
use_label = bquote(italic(R)^2 == .(format(r2, nsmall=2, digits=2)))
mtext(use_label, 3, adj=1, line=0)
mtext('D',3,adj=0,line=0, font=2, cex=2)

# Tree diversity
plot(lichen.rich~PIE.ba.tree, data=trans_data_test,
	xlab='Tree Diversity (local)', ylab='Local Species Richness', 
	las=1, xlim=as.numeric(myranges['PIE.ba.tree',c('lower','upper')]), ylim=c(0,40),
	col='#00000050', lwd=2)
use_mod = glm.nb(lichen.rich~PIE.ba.tree, data=trans_data_test, link='log')
use_coef = coef(use_mod)
use_x = seq(min(trans_data_test$PIE.ba.tree), max(trans_data_test$PIE.ba.tree), length.out=100)
use_y = exp(use_coef[1]+use_coef[2]*use_x)
lines(use_x, use_y, lwd=5, col='white')
lines(use_x, use_y, lwd=4, col='black')
r2 = r.squaredLR(use_mod, null=glm.nb(lichen.rich~1, data=trans_data_test, link='log'))
r2 = attr(r2, 'adj.r.squared')
use_label = bquote(italic(R)^2 == .(format(r2, nsmall=2, digits=2)))
mtext(use_label, 3, adj=1, line=0)
mtext('E',3,adj=0,line=0, font=2, cex=2)

# Regional Tree Diversity
plot(lichen.rich~regS_tree, data=trans_data_test,
	xlab='Tree Richness (regional)', ylab='Local Species Richness', 
	las=1, xlim=as.numeric(myranges['regS_tree',c('lower','upper')]), ylim=c(0,40),
	col='#00000050', lwd=2, axes=F)
axis(1, at=seq(0,160,40))
axis(2)
box()
use_mod = glm.nb(lichen.rich~regS_tree+I(regS_tree^2), data=trans_data_test, link='log')
use_coef = coef(use_mod)
use_x = seq(min(trans_data_test$regS_tree), max(trans_data_test$regS_tree), length.out=100)
use_y = exp(use_coef[1]+use_coef[2]*use_x+use_coef[3]*(use_x^2))
lines(use_x, use_y, lwd=5, col='white')
lines(use_x, use_y, lwd=4, col='black')
r2 = r.squaredLR(use_mod, null=glm.nb(lichen.rich~1, data=trans_data_test, link='log'))
r2 = attr(r2, 'adj.r.squared')
use_label = bquote(italic(R)^2 == .(format(r2, nsmall=2, digits=2)))
mtext(use_label, 3, adj=1, line=0)
mtext('F',3,adj=0,line=0, font=2, cex=2)

dev.off()


### Compare local ~ regional richness for AllSp, Parmeliaceae and Physciaceae
# Note that model form (linear or quadratic) was decided from univariate_model_shapes, which was fit to the 'fitting' dataset whereas here we are plotting the 'testing' dataset.

svg('./Figures/compare local-regional richness across taxa.svg', height=4, width=11)

par(mfrow=c(1,3))
par(mar=c(4,4,1.5,1))
par(mgp=c(2.4,0.7,0))
par(cex.axis=1.2)
par(cex.lab=1.2)
par(pch=1)
par(cex=1)
par(las=1)

# All Species
plot(lichen.rich~regS, data=trans_data_test,
	xlab='Regional Species Richness', ylab='Local Species Richness', 
	las=1, xlim=c(125,205), ylim=c(0,40),
	col='#00000050', lwd=2, axes=F)
axis(1, at=seq(125,200,25))
axis(2)
box()
use_mod = glm.nb(lichen.rich~regS, data=trans_data_test, link='log')
use_coef = coef(use_mod)
use_x = seq(min(trans_data_test$regS), max(trans_data_test$regS), length.out=100)
use_y = exp(use_coef[1]+use_coef[2]*use_x)
lines(use_x, use_y, lwd=5, col='white')
lines(use_x, use_y, lwd=4, col='black')
r2 = r.squaredLR(use_mod, null=glm.nb(lichen.rich~1, data=trans_data_test, link='log'))
r2 = attr(r2, 'adj.r.squared')
use_label = bquote(italic(R)^2 == .(format(r2, nsmall=2, digits=2)))
mtext(use_label, 3, adj=1, line=0)
mtext('All Species',3,adj=0,line=0, font=2, cex=1)

# All Species
plot(Parmeliaceae~regParm, data=trans_data_test,
	xlab='Regional Species Richness', ylab='Local Species Richness', 
	las=1, xlim=c(35,75), ylim=c(0,30),
	col='#00000050', lwd=2, axes=F)
axis(1, at=seq(35,75,5))
axis(2)
box()
use_mod = glm.nb(Parmeliaceae~regParm, data=trans_data_test, link='log')
use_coef = coef(use_mod)
use_x = seq(min(trans_data_test$regParm), max(trans_data_test$regParm), length.out=100)
use_y = exp(use_coef[1]+use_coef[2]*use_x)
lines(use_x, use_y, lwd=5, col='white')
lines(use_x, use_y, lwd=4, col='black')
r2 = r.squaredLR(use_mod, null=glm.nb(Parmeliaceae~1, data=trans_data_test, link='log'))
r2 = attr(r2, 'adj.r.squared')
use_label = bquote(italic(R)^2 == .(format(r2, nsmall=2, digits=2)))
mtext(use_label, 3, adj=1, line=0)
mtext('Parmeliaceae',3,adj=0,line=0, font=2, cex=1)

# Physciaceae
plot(Physciaceae~regPhys, data=trans_data_test,
	xlab='Regional Species Richness', ylab='Local Species Richness', 
	las=1, xlim=c(25,50), ylim=c(0,15),
	col='#00000050', lwd=2, axes=F)
axis(1, at=seq(25,50,5))
axis(2)
box()
use_mod = glm.nb(Physciaceae~regPhys+I(regPhys^2), data=trans_data_test, link='log')
use_coef = coef(use_mod)
use_x = seq(min(trans_data_test$regPhys), max(trans_data_test$regPhys), length.out=100)
use_y = exp(use_coef[1]+use_coef[2]*use_x+use_coef[3]*(use_x^2))
lines(use_x, use_y, lwd=5, col='white')
lines(use_x, use_y, lwd=4, col='black')
r2 = r.squaredLR(use_mod, null=glm.nb(Physciaceae~1, data=trans_data_test, link='log'))
r2 = attr(r2, 'adj.r.squared')
use_label = bquote(italic(R)^2 == .(format(r2, nsmall=2, digits=2)))
mtext(use_label, 3, adj=1, line=0)
mtext('Physciaceae',3,adj=0,line=0, font=2, cex=1)

dev.off()



##########################################
### Old Code

png('./Figures/New Coordinates/legend small regions.png', height=600, width=600)
plot.new()
legend('center', levels(usereg), col=regcol[1:nlevels(usereg)], 
	pch=regpch[1:nlevels(usereg)], cex=4, lwd=3, lty=0, bty='n')
dev.off()

# Highlight by fit and test data sets
dataset = factor(rownames(trans_data) %in% testplots)
usecol = c('black','blue')

for(i in c('wetness','regS','PIE.ba.tree','wood_SG.ba','totalNS')){
png(paste('./Figures/New Coordinates/highlight test data ',i,' vs lichen richness.png', sep=''), height=1500, width=1500)
par(mar=c(10,11,1,1))
plot(lichen.rich~trans_data[,i], data=trans_data, type='n', las=1, ylim=c(0,40), xlim=as.numeric(myranges[i,]),
	xlab='', ylab='',	axes=F
)
points(lichen.rich~trans_data[,i], data=trans_data, col=usecol[dataset], pch=1, cex=4, lwd=3)
axis(1, cex.axis=4, las=1, padj=1, lwd=2)
axis(2, cex.axis=4, las=1, lwd=2, at=seq(0,40,5))
mtext(varnames[i,'displayName'],1,7, cex=4)
mtext('Local lichen species richness',2,8, cex=4)
dev.off()
}
