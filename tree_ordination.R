### This script ordinates tree community composition to be used as a variable in the analysis of FIA lichen species richness.

####################################################
### Load data
source('C:/Users/jrcoyle/Documents/UNC/Projects/FIA Lichen/GitHub/FIA-Lichens/load_data.R')

## Tree site X speces matrices for legacy and current data
trees_current = read.csv('./Data/TreeData/fia_abundance_matrix_ba_legacy_2012-05-16.csv', row.names=1)
trees_legacy = read.csv('./Data/TreeData/fia_abundance_matrix_ba_current_2012-05-13.csv', row.names=1)

## Data frame of forest variables for each plot
master_forest = read.csv('./Data/TreeData/master_data_forest.csv')
years_used = read.csv('match_lichen_tree_current.csv')
spnames = read.csv('C:/Users/jrcoyle/Dropbox/Documentation/REF_SPECIES.CSV')
spgroup = read.csv('C:/Users/jrcoyle/Dropbox/Documentation/REF_SPECIES_GROUP.CSV')

## Write out a list of tree species that are in this data set
treesp = data.frame(SPCD = substring(colnames(siteXsp), first=2))
treesp = merge(treesp, spnames[,c('SPCD','COMMON_NAME','GENUS','SPECIES','VARIETY','SUBSPECIES','JENKINS_SPGRPCD','MAJOR_SPGRPCD')], all.x=T)
treesp$Conifer = treesp$MAJOR_SPGRPCD %in% c(1,2)

write.csv(treesp, './Data/TreeData/lichen_tree_species.csv')

###################################################
### Ordination of trees


## Combine legacy and current matrices

# Vector of all tree species
tree_sp = unique(c(colnames(trees_current), colnames(trees_legacy)))

# Find species missing in each data set
current_miss = !(tree_sp %in% colnames(trees_current))
legacy_miss = !(tree_sp %in% colnames(trees_legacy))

# Make empty matrix to hold data
siteXsp = matrix(0, nrow=nrow(trees_current)+nrow(trees_legacy), ncol=length(tree_sp),
	dimnames=list(c(rownames(trees_current), rownames(trees_legacy)), tree_sp)
)

# Change site X species data frames into matrices
trees_current = as.matrix(trees_current)
trees_legacy = as.matrix(trees_legacy)

# Fill in current and legacy data
siteXsp[rownames(trees_current),!current_miss] = trees_current
siteXsp[rownames(trees_legacy),!legacy_miss] = trees_legacy

# Change rownames in siteXsp matrix to match yrplot.ids of lichen plots: i.e. those used in master_forest
newnames = rownames(siteXsp)
names(newnames) = newnames
newnames[years_used$tree.yrplot.id] = years_used$lich.yrplot.id
rownames(siteXsp) = newnames
missing_plots = master_forest$yrplot.id[which(!(master_forest$yrplot.id %in% rownames(siteXsp)))]

# Only use sites in the data that will be modeled
siteXsp = siteXsp[master_forest$yrplot.id,]
sum(is.na(siteXsp)) # Make sure there is no missing data due to bad matches


## Ordination
library(ape) #pcoa
library(ade4) #is.euclid
library(vegan) # decostand, rda

# Transform abundances using log10 to increase normality
siteXsp_trans = log10(siteXsp+1)

# Use Hellinger transformation: sqrt of relative abundances
siteXsp_hel = decostand(siteXsp_trans, method='hellinger')

# Do PCA on transformed distances
# This is identical to doing PCoA on Hellinger distances, according to vegan package manual
# Hellinger distance is D17 in Legendre & Legendre 2012 and is recommended w/ PCoA by Rao 1995
# This is the Euclidean distance between sqrt transformed relative abunances.
treepca = rda(siteXsp_hel) # since only matrix X is supplied, this will give PCA
treepca2 = prcomp(siteXsp_hel)

# Biplots
biplot(treepca, choices=c(1,2), type=c('text','points'))
biplot(treepca, choices=c(5,6), type=c('text','points'))

# Species scores- loadings
spscore = scores(treepca,1:234,display='species')

PC1species = spscore[order(spscore[,1], decreasing=T),1] # .8 Ulmus rubra, .5 Quercus lyrata, .2 Taxodium ascendens, .2 Prunus serotina
PC2species = spscore[order(spscore[,2], decreasing=T),2] # .9 Ulmus rubra, .6 Quercus lyrata, .5 Acer pensylvanicum, .3 Pinus virginiana
PC3species = spscore[order(spscore[,3], decreasing=T),3] # .5 Pseudotsuga menzeisii, .3 Acer pensylvanicum, .2 Ulmus rubra, .2 Pinus ponderosa, .2 Acer macrophyllum
PC4species = spscore[order(spscore[,4], decreasing=T),4] # .78 Pseudotsuga menzeisii, .4 Pinus contorta, .3 Abies lasiocarpa, .3 Populus tremuloides, .2 Picea engelmanii

# Site scores
sitescore = as.data.frame(scores(treepca, 1:234, display='sites'))

# On which axis does each species load the strongest?
apply(spscore, 1, function(x) which(x==max(x)))

# Which species have the strongest loadings overall and which axes do these correspond to?
maxloads= apply(spscore, 1, max)
maxloads = maxloads[order(maxloads, decreasing=T)

# 1.10: Juniperus osteosperma
# 0.99: Pinus virginiana
# 0.95: Ulmus rubra
# 0.78: Pseudotsuga menziesii
# 0.63: Acer pensylvanicum
# 0.59: Pinus ponderosa
# 0.55: Populus tremuloides
# 0.52: Pinus monophylla
# 0.49: Ulmus americana
# 0.43: Juglans nigra
# 0.39: Prunus serotina
# 0.37: Pinus contorta
# 0.37: Pinus edulis

# Are these the least/most common species?
spabun = colSums(siteXsp) 
plot(1:length(spabun), log10(spabun[order(spabun, decreasing=T)]))
abline(v=c(5,10,15,20), col=2)

# What are the top 20 most abundant species?
top20 = head(spabun[order(spabun, decreasing=T)], n=20)
rownames(spnames) = spnames$SPCD
spnames[substring(names(top20), first=2),c('GENUS','SPECIES')]

# Why is Acer pensylvanicum so abundant- is this just one plot?- no its not.

# How many plots are each of these species in?
colSums(siteXsp[,names(top20)]>0)

# What are the top 20 most prevalent species?
spprev = colSums(siteXsp>0)
top20prev = head(spprev[order(spprev, decreasing=T)], n=20)
cbind(spnames[substring(names(top20prev), first=2),c('GENUS','SPECIES')], top20prev)
plot(1:length(spprev), spprev[order(spprev, decreasing=T)])

# How many species are in in at least 5% of plots? (>111 plots)
sp5pct = names(spprev[spprev>112])

# Are there any plots without any of these species?
sum(rowSums(siteXsp[,sp5pct])==0)

## Map PC1-4

library(lattice)
library(sp)
library(rgdal)

# Define map colors
ncuts=16
mycol = read.csv('C:/Users/jrcoyle/Documents/UNC/Projects/blue2red_10colramp.txt')
mycol = apply(mycol,1,function(x) rgb(x[1],x[2],x[3],maxColorValue=256))
mycol = mycol[10:1]
mycolramp = colorRampPalette(mycol)(ncuts)

# Map projection - Lambert Equal Area
plot_prj = paste("+proj=laea +lat_0=40 +lon_0=-97 +units=km",sep='')

# Read in North America outline and re-project
OUTLINES = readOGR('../../GIS shape files/N Am Outline','na_base_Lambert_Azimuthal')
OUTLINES.laea = spTransform(OUTLINES,CRS(plot_prj))

# Make lichen plot data into spatial data
sp_data = cbind(sitescore, yrplot.id=rownames(sitescore))
sp_data = merge(master[,c('yrplot.id','LAT','LON')], sp_data, all.x=T, all.y=F)
sum(is.na(sp_data[,'PC1'])) # all matched!
sp_data = subset(sp_data, !is.na(LAT)) # remove plots w/o geographic coordinates
coordinates(sp_data) = c('LON','LAT') # define columns that have geographic coordinates
proj4string(sp_data) = CRS("+proj=longlat") # define original projection of data
sp_data = spTransform(sp_data, CRS(plot_prj)) # re-project data


apply(sitescore, 2, range)
colcuts = seq(-.8, .8, .1)

pdf('./Figures/New Coordinates/Tree group PCA Maps.pdf', height=8, width=12)
	trellis.par.set(axis.line=list(col=NA))

	# Define levels for colors
	for(y in paste('PC',1:5, sep='')){
	
	print(
	spplot(sp_data, y, ylim=c(-1600,1500), main='', panel=function(x,y,subscripts,...){
		sp.polygons(OUTLINES.laea, fill='white')
		panel.pointsplot(x,y,...)	
	}, cuts=colcuts, cex=1.5, col.regions = mycolramp, auto.key=F,
	key=list(x=0.05,y=.2, corner=c(0,.5), title=y,
		rectangles=list(col=mycolramp, size=3, border='transparent'),
		text=list(paste(colcuts[1:ncuts],colcuts[2:(ncuts+1)], sep=':'), cex=1.5))
	)
	)
	}
dev.off()

### Do ordination on tree species groups
treesp = read.csv('./Data/TreeData/lichen_tree_species.csv')

## Combine species in siteXsp matrix into species groups
# Put group codes in same order as column names
rownames(treesp) = paste('X',treesp$SPCD, sep='')
groups = treesp[colnames(siteXsp),'JENKINS_SPGRPCD']

# Sum up species in same group for each site
siteXgrp = t(apply(siteXsp, 1, function(x) tapply(x, groups, sum)))

## PCA of Hellinger transformed abundances
# Transform abundances using log10 to increase normality
siteXgrp_trans = log10(siteXgrp+1)

# Use Hellinger transformation: sqrt of relative abundances
siteXgrp_hel = decostand(siteXgrp_trans, method='hellinger')

# PCA
grppca = rda(siteXgrp_hel) 

# Scores
grpscore = scores(grppca,1:10,display='species')
sitescore = as.data.frame(scores(grppca, 1:10, display='sites'))

plot(sitescore[,1], sitescore[,2])

# Write out PC1-3 scores
write.csv(data.frame(yrplot.id = rownames(sitescore), sitescore[,1:3]),'./Data/TreeData/tree_funcgrp_pca1-3.csv', row.names=F)


# Variation explained:
eigenvals(grppca)/sum(eigenvals(grppca))


# Use above code to make maps


# Which species groups are loading on which axes?
jenkins = read.csv('./Data/TreeData/jenkins_codes.csv')
rownames(grpscore) = jenkins[rownames(grpscore),'Description']

grpscore[order(grpscore[,1]),1:3]

biplot(grppca, c(1,2), type=c('text','points'))


# Which axis is lichen richness correlated with?
use_data = cbind(sitescore, yrplot.id=rownames(sitescore))
use_data = merge(use_data, master[,c('yrplot.id','lichen.rich', 'Parmeliaceae','Physciaceae')])

pdf('./Figures/New Coordinates/lichen richness vs tree group pca.pdf', height=6, width=12)
par(mar=c(4,4,2,1))
par(mfrow=c(2,5))
for(y in c('lichen.rich','Parmeliaceae','Physciaceae')){
for(x in paste('PC',1:10, sep='')){
	plot(use_data[,x], use_data[,y], ylab=varnames[y,'displayName'], xlab=x, main='')
	r = format(cor(use_data[,x], use_data[,y]), digits=3)
	mtext(r, 3, -1, adj=1, cex=2, col=2)
}}
dev.off()


###############################################################################
### Calculate conifer area

rownames(treesp) = paste('X',treesp$SPCD, sep='')
coniferTF = treesp[colnames(siteXsp),'Conifer']



unique(spnames[,c('E_SPGRPCD','W_SPGRPCD','MAJOR_SPGRPCD','JENKINS_SPGRPCD')])







