## This script makes univariate plots of lichen richness versus environmental variables



##############################################################################
### Make Univariate plots

regions = read.csv('fia_lichen_regions.csv', row.names=1)

# trans_data should have all plots, no outliers, and log/sqrt transformed where necessary

# highlight by region plots
usereg = factor(regions[master[rownames(trans_data),'state.abbr'],'big_region'])
usereg = factor(regions[master[rownames(trans_data),'state.abbr'],'sm_region'])

regcol = c('cornflowerblue','darkred','orange','forestgreen','purple','cyan3')
regpch = c(0:5)

myranges = data.frame(lower = c(-4, 120, 0, 0.3, 0), upper = c(3, 210, 1, 0.8, 1100), row.names=c('wetness','regS','PIE.ba.tree','wood_SG.ba','totalNS'))

for(i in c('wetness','regS','PIE.ba.tree','wood_SG.ba','totalNS')){
png(paste('./Figures/New Coordinates/highlight big region ',i,' vs lichen richness.png', sep=''), height=1500, width=1500)
par(mar=c(10,11,1,1))
plot(lichen.rich~trans_data[,i], data=trans_data, type='n', las=1, ylim=c(0,40), xlim=as.numeric(myranges[i,]),
	xlab='', ylab='',	axes=F
)
points(lichen.rich~trans_data[,i], data=trans_data, col=regcol[usereg], pch=regpch[usereg], cex=4, lwd=3)
axis(1, cex.axis=4, las=1, padj=1, lwd=2)
axis(2, cex.axis=4, las=1, lwd=2, at=seq(0,40,5))
mtext(varnames[i,'displayName'],1,7, cex=4)
mtext('Local lichen species richness',2,8, cex=4)
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
