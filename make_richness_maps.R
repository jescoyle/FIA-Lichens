# This script create maps of lichen richness in plots used in the analysis.


source('C:/Users/jrcoyle/Documents/UNC/Projects/FIA Lichen/GitHub/FIA-Lichens/load_data.R')


#######################################################################
# Richness Maps

ncuts=10

mycol = read.csv('C:/Users/jrcoyle/Documents/UNC/Projects/blue2red_10colramp.txt')
mycol = apply(mycol,1,function(x) rgb(x[1],x[2],x[3],maxColorValue=256))
mycol = mycol[10:1]
mycolramp = colorRampPalette(mycol)(ncuts)

## Final lichen richness map
library(lattice)
library(sp)
library(rgdal)

# Map projection
plot_prj = paste("+proj=laea +lat_0=40 +lon_0=-97 +units=km",sep='')

# N. Am. outline
OUTLINES = readOGR('../../GIS shape files/N Am Outline','na_base_Lambert_Azimuthal')
OUTLINES.laea = spTransform(OUTLINES,CRS(plot_prj))

# Forest map: crop and reproject
library(raster)
forest = raster('../GIS Data/forest types/forest_binary.grd')
forest_low48 = crop(forest, extent(matrix(c(-2100000, -2100000,2600000,1000000), nrow=2)))
newext = projectExtent(forest_low48, crs=CRS(plot_prj))
res(newext) = 10
forest_low48_laea = projectRaster(forest_low48, newext, method='ngb')
forest_low48_grid = as(forest_low48_laea, 'SpatialGridDataFrame')

# Make spatial points data frame with only 1923 plots used in models
sp_data = master[rownames(model_data),]
coordinates(sp_data) = c('LON','LAT')
proj4string(sp_data) = CRS("+proj=longlat")
sp_data = spTransform(sp_data, CRS(plot_prj))

pdf('./Figures/Map lichen richness white.pdf', height=8, width=12)
trellis.par.set(axis.line=list(col=NA))
colcuts = c(1,4,8,11,15,19,22,26,29,33,37) #c(1,2,3,4,5,7,8,9,10,11,13) #c(1,3,6,9,11,14,17,19,22,25,28)  
spplot(sp_data, 'lichen.rich', ylim=c(-1600,1500), main='', panel=function(x,y,subscripts,...){
	sp.polygons(OUTLINES.laea, fill='white')
	panel.pointsplot(x,y,...)
}, cuts=colcuts, cex=1.5, col.regions = mycolramp, auto.key=F,
key=list(x=1,y=.3, corner=c(1,.5), title='Species\nRichness',
	rectangles=list(col=mycolramp, size=3, border='transparent'),
	text=list(c(paste(colcuts[1],colcuts[2], sep='-'), paste((colcuts+1)[2:(ncuts)], colcuts[2:ncuts+1], sep='-')))) 
)
dev.off()

pdf('./Figures/Map lichen richness with forest.pdf', height=8, width=12)
trellis.par.set(axis.line=list(col=NA))
spplot(forest_low48_grid, 'layer', ylim=c(-1600,1500), main='',
	panel=function(x,y,z,subscripts,...){
		sp.polygons(OUTLINES.laea, fill='white')
		panel.levelplot(x,y,z,subscripts,...)
		sp.points(sp_data, col=mycolramp[cut(sp_data@data[,'lichen.rich'], colcuts, include.lowest=T)], pch=16, cex=1.5)
	}, col.regions=c(0,'grey60'), colorkey=F,
	key=list(x=1,y=.3, corner=c(1,.5), title='Species\nRichness',
		rectangles=list(col=mycolramp, size=3, border='transparent'),
		text=list(c(paste(colcuts[1],colcuts[2], sep='-'), paste((colcuts+1)[2:(ncuts)], colcuts[2:ncuts+1], sep='-')))) 
)
dev.off()

colcuts =  c(1,3,6,9,11,14,17,19,22,25,28)  #c(1,2,3,4,5,7,8,9,10,11,13) 
pdf('./Figures/Map Parmeliaceae richness with forest.pdf', height=8, width=12)
trellis.par.set(axis.line=list(col=NA))
spplot(forest_low48_grid, 'layer', ylim=c(-1600,1500), main='',
	panel=function(x,y,z,subscripts,...){
		sp.polygons(OUTLINES.laea, fill='white')
		panel.levelplot(x,y,z,subscripts,...)
		sp.points(sp_data, col=mycolramp[cut(sp_data@data[,'Parmeliaceae'], colcuts, include.lowest=T)], pch=16, cex=1.5)
	}, col.regions=c(0,'grey60'), colorkey=F,
	key=list(x=1,y=.3, corner=c(1,.5), title='Species\nRichness',
		rectangles=list(col=mycolramp, size=3, border='transparent'),
		text=list(c(paste(colcuts[1],colcuts[2], sep='-'), paste((colcuts+1)[2:(ncuts)], colcuts[2:ncuts+1], sep='-')))) 
)
dev.off()

## regional richness
mycolramp100 = colorRampPalette(mycol)(100)

pdf('./Figures/Map regional lichen richness white.pdf', height=8, width=12)
trellis.par.set(axis.line=list(col=NA))
colcuts = seq(120,220,10)  
spplot(sp_data, 'regS', ylim=c(-1600,1500), main='', panel=function(x,y,subscripts,...){
	sp.polygons(OUTLINES.laea, fill='white')
	panel.pointsplot(x,y,...)
}, cuts=100, cex=1.5, col.regions = mycolramp100, auto.key=F)

dev.off()




