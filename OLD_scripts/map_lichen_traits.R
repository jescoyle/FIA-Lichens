## This script makes maps of traits for 1923 plots used in FIA Lichen analysis

source('./UNC/Projects/FIA Lichen/GitHub/FIA-Lichens/load_data.R')


# Required packages
library(lattice)
library(sp)
library(rgdal)

#################################################################################
## Read in extra data files and define variables##

# File describing infomation on traits: categorical vs numeric
traitnames = read.csv('lias_trait_types.csv', row.names=1)

# Define map colors
ncuts=8
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
sp_data = subset(master, !is.na(LAT)) # remove plots w/o geographic coordinates
coordinates(sp_data) = c('LON','LAT') # define columns that have geographic coordinates
proj4string(sp_data) = CRS("+proj=longlat") # define original projection of data
sp_data = spTransform(sp_data, CRS(plot_prj)) # re-project data

# Read in trait prevalence in plots, summarized as z-scores relative to a null of random assembly
zscores= read.csv('./Data/LichenTraits/trait_z-scores_across_plots.csv', row.names=1)

###################################################################################
## Map binary traits ##

# Define set of traits to map
use_traits = c('pruinose','cilia','isidia','soredia','ascomata',
		'secondary_photo','bullate','macculate','pseudocyphellate','hairy',
		'atranorin','usnic_acid','pulvinic_acid','pulvinic_dilactone','parietin',
		'parietinic_acid')

# Define levels for colors
colcuts = seq(0,1,length.out=ncuts+1)

# Loop through traits and make maps
for(tr in use_traits){
	
	pdf(paste('./Figures/New Coordinates/Trait Maps/map ',tr,'.pdf', sep=''), height=8, width=12)
	trellis.par.set(axis.line=list(col=NA))

	print(
	spplot(sp_data, tr, ylim=c(-1600,1500), main='', panel=function(x,y,subscripts,...){
		sp.polygons(OUTLINES.laea, fill='white')
		panel.pointsplot(x,y,...)	
	}, cuts=colcuts, cex=1.5, col.regions = mycolramp, auto.key=F,
	key=list(x=0.05,y=.2, corner=c(0,.5), title=paste('Proportion of Species\nwith',traitnames[tr,'displayName'],sep=' '),
		rectangles=list(col=mycolramp, size=3, border='transparent'),
		text=list(paste(colcuts[1:ncuts],colcuts[2:(ncuts+1)], sep='-'), cex=1.5)) 
	)
	)
	dev.off()

}


# Merge z-score data with spatial data
zscores$yrplot.id = rownames(zscores)
zscores_sp = merge(master[,c('yrplot.id','LAT','LON')], zscores)
zscores_sp = subset(zscores_sp, !is.na(LAT)) # remove plots w/o geographic coordinates
coordinates(zscores_sp) = c('LON','LAT') # define columns that have geographic coordinates
proj4string(zscores_sp) = CRS("+proj=longlat") # define original projection of data
zscores_sp = spTransform(zscores_sp, CRS(plot_prj)) # re-project data


# Define levels for colors
t(apply(zscores[,1:17], 2, range)) # Examine by hand
colcuts = c(-12,-4,-2,-1,0,1,2,4,50)

# Loop through traits and make maps of z-scores
for(tr in use_traits){
	
	pdf(paste('./Figures/New Coordinates/Trait Maps/z-score map ',tr,'.pdf', sep=''), height=8, width=12)
	trellis.par.set(axis.line=list(col=NA))

	print(
	spplot(zscores_sp, tr, ylim=c(-1600,1500), main=paste('Species\nwith',traitnames[tr,'displayName'],sep=' '), panel=function(x,y,subscripts,...){
		sp.polygons(OUTLINES.laea, fill='white')
		panel.pointsplot(x,y,...)	
	}, cuts=colcuts, cex=1.5, col.regions = mycolramp, auto.key=F,
	key=list(x=0.05,y=.2, corner=c(0,.5), title='z-score',
		rectangles=list(co=mycolramp, size=3, border='transparent'),
		text=list(paste(colcuts[1:ncuts],colcuts[2:(ncuts+1)], sep=' - '), cex=1.5)) 
	)
	)
	dev.off()

}



####################################################################################
## Map numeric traits ##

# Do numerc traits separately:

c('ascospore_length_mid','ascospore_width_mid','conidia_length_mid')











