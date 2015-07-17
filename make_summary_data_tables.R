## This script makes summary tables of lichen data

source('load_data.R')


################################################################################
### Information about plots


yrXstate = table(master$state.abbr, master$YEAR) # may need to change to MEASYEAR
write.csv(yrXstate[rowSums(yrXstate)>0,], './Data/All lichen plots by state and year.csv')

master.recent = subset(master, YEAR >=1997)

coordsXstate = table(master.recent$state.abbr,is.na(master.recent$LAT))

write.csv(coordsXstate[rowSums(coordsXstate)>0,], './Data/recent lichen plots by state and coordinates.csv')

################################################################################
### Site x Species abundance matrix

# Get list of species actually in data
lichen = read.csv('./Data/FIA_lichens_parsed_2012-05-13.csv') ## Current Data ##
lichen.current = subset(lichen, yrplot.id %in% master$yrplot.id)[,c('yrplot.id','ABUNDANCE_CLASS','LICH_SPPCD')]
load('./Data/LegacyData/legacy_lichen_data_parsed_2012-05-15.Rdata') ## Legacy Data ##
lichen.legacy = subset(lichen, yrplot.id %in% master$yrplot.id)[,c('yrplot.id','ABUNDANCE_CLASS','LICH_SPPCD')]

lichen = rbind(lichen.current, lichen.legacy)
siteXsp = as.matrix(xtabs(ABUNDANCE_CLASS~yrplot.id+LICH_SPPCD, data=lichen))

# Remove species that never occur
siteXsp = siteXsp[,colSums(siteXsp)!=0]

# Remove unknown species
siteXsp = siteXsp[,-which(colnames(siteXsp)=='9999')]


### Where are the most species rich plots?
ordered_data = model_data[order(model_data$lichen.rich, decreasing=T),] # aslo examined fric and genus.rich
ordered_data30 = ordered_data[ordered_data$lichen.rich>=30,c('yrplot.id','STATE_NAME','COUNTYNM','LAT','LON','lichen.rich','genus.rich','regS','fric')]

write.csv(ordered_data30, './Paper/Plots with at least 30 species.csv', row.names=F)

write.csv(siteXsp, 'siteXsp_matrix_current+legacy.csv', row.names=T)

siteXsp = read.csv('./Data/siteXsp_matrix_current+legacy.csv', row.names=1, check.names=F)

### How many taxa are just IDed to genus
lichsp = read.csv('./Data/LICHEN_SPECIES_SUMMARY.CSV')
model_data = read.csv('./Data/fia_lichen_model_data.csv', row.names=1)
gen_cds = subset(lichsp, SPECIES=='')$LICH_SPPCD

siteXsp_gen = siteXsp[rownames(model_data),colnames(siteXsp) %in% as.character(gen_cds)]
model_data$justgenus = rowSums(siteXsp_gen)
table(model_data$justgenus)

master = read.csv('./Data/fia_lichen_master_data.csv', row.names=1)
spdata = data.frame(master[rownames(model_data), c('LON','LAT')], justgenus = model_data$justgenus)

library(sp)
library(rgdal)

coordinates(spdata) = c('LON','LAT')
proj4string(spdata) = CRS('+proj=longlat')

mycol = read.csv('C:/Users/jrcoyle/Documents/UNC/Projects/blue2red_10colramp.txt')
mycol = apply(mycol,1,function(x) rgb(x[1],x[2],x[3],maxColorValue=256))
mycol = mycol[10:1]

# Map projection
plot_prj = paste("+proj=laea +lat_0=40 +lon_0=-97 +units=km",sep='')

# N. Am. outline
OUTLINES = readOGR('../../GIS shape files/N Am Outline','na_base_Lambert_Azimuthal')
OUTLINES.laea = spTransform(OUTLINES,CRS(plot_prj))

spdata = spTransform(spdata, CRS(plot_prj))

colcuts = seq(1, 15, 2)
ncuts=7
mycolramp = colorRampPalette(mycol)(ncuts)

pdf('./Figures/Maps/Just genus ID richness map.pdf', height=8, width=12)
spplot(spdata, 'justgenus', ylim=c(-1600,1500), main='', panel=function(x,y,subscripts,...){
	sp.polygons(OUTLINES.laea, fill='white')
	panel.pointsplot(x,y,...)
}, cuts=colcuts, cex=1.5, col.regions = mycolramp, auto.key=F,
key=list(x=1,y=.3, corner=c(1,.5), title='# Taxa Ided\nto Genus only',
	rectangles=list(col=mycolramp, size=3, border='transparent'),
	text=list(c(paste(colcuts[1],colcuts[2], sep='-'), paste((colcuts+1)[2:(ncuts)], colcuts[2:ncuts+1], sep='-')))) 
)
dev.off()

png('./Figures/richness vs just genus taxa.png', height=300, width=300)
par(mar=c(4,4,1,1))
plot(justgenus~lichen.rich, data=model_data, las=1, xlab='Species richness', ylab='# Taxa IDed to Genus only', 
	pch=16, col='#00000040')
dev.off()

###################################################################
### How many occurences are there of species missing from LIAS? ###

# Species list
myspecies = read.table('C:/Users/jrcoyle/Documents/Projects/Lichen Traits/fia_lichen_master_list.txt', header=T, sep='\t')
noLias = subset(myspecies, is.na(LIAS_ID))

## Current Data ##
lichen = read.csv('./Data/FIA_lichens_parsed_2012-05-13.csv')

lichen.noLias = subset(lichen, LICH_SPPCD %in% noLias$LICH_SPPCD)
noLias$nplots.current = sapply(noLias$LICH_SPPCD, function(x) nrow(subset(lichen.noLias, LICH_SPPCD==x)))

## Legacy Data ##
load('./Data/LegacyData/legacy_lichen_data_parsed_2012-05-15.Rdata')

lichen.noLias = subset(lichen, LICH_SPPCD %in% noLias$LICH_SPPCD)
noLias$nplots.legacy = sapply(noLias$LICH_SPPCD, function(x) nrow(subset(lichen.noLias, LICH_SPPCD==x)))

lichen.noLias = subset(lichen, (LICH_SPPCD %in% noLias$LICH_SPPCD)&(year>=1997))
noLias$nplots.legacyPost97 = sapply(noLias$LICH_SPPCD, function(x) nrow(subset(lichen.noLias, LICH_SPPCD==x)))

# Number of plots in data I will be analyzing (post 1997)
noLias$nplots = noLias$nplots.current + noLias$nplots.legacyPost97

# Ordered

noLias[order(noLias$nplots, decreasing=T), ]

write.csv(noLias[order(noLias$nplots, decreasing=T), ], 'frequency of missing LIAS species in FIA.csv', row.names=F)
