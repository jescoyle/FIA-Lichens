## This script makes univariate plots of lichen richness versus environmental variables


source('./UNC/Projects/FIA Lichen/GitHub/FIA-Lichens/load_data.R')
source('./GitHub/FIA-Lichens/fia_lichen_analysis_functions.R')

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
use_vars = c('regS','wetness','PIE.ba.tree','wood_SG.ba','totalNS')

# Define plotting ranges for each of the six variables to be plotted
myranges = data.frame(lower = c(120, -4, 0, 0.3, 0), 
	upper = c(210, 3, 1, 0.8, 1100), 
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
