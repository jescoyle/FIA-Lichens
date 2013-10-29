## This script makes maps of traits for 1923 plots used in FIA Lichen analysis

source('load_data.R')


# Required packages
library(lattice)
library(sp)
library(rgdal)

# Read in extra data files
traitnames = read.csv('lias_trait_types.csv', row.names=1)

# Define map colors
ncuts=8
mycol = read.csv('C:/Users/jrcoyle/Documents/UNC/Projects/blue2red_10colramp.txt')
mycol = apply(mycol,1,function(x) rgb(x[1],x[2],x[3],maxColorValue=256))
mycol = mycol[10:1]
mycolramp = colorRampPalette(mycol)(ncuts)

# Map projection - Lambert Equal Area
plot_prj = paste("+proj=laea +lat_0=40 +lon_0=-97 +units=km",sep='')

# North America outline
OUTLINES = readOGR('../../GIS shape files/N Am Outline','na_base_Lambert_Azimuthal')
OUTLINES.laea = spTransform(OUTLINES,CRS(plot_prj))

# Make data into spatial data
sp_data = subset(master, !is.na(LAT))
coordinates(sp_data) = c('LON','LAT')
proj4string(sp_data) = CRS("+proj=longlat")
sp_data = spTransform(sp_data, CRS(plot_prj))


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

# Do numerc traits separately:

c('ascospore_length_mid','ascospore_width_mid','conidia_length_mid')











