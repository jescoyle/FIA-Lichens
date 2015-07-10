# This script diversity of lichens from occurence records downloaded from
# Consortium of North American Lichen Herbaria in the fall of 2014.

library(sp)
library(rgdal)
library(vegan) # rarefy
library(stringr) #str_trim
#library(rgeos) # gBuffer

source('C:/Users/jrcoyle/Documents/UNC/Projects/FIA Lichen/GitHub/FIA-Lichens/load_data.R')
# Read in locations of master data of FIA plots with plot data and environmental data
master_locs = read.csv('./Data/fia_lichen_plot_locations.csv')

# Move to CNALH Directory
setwd('C:/Users/jrcoyle/Documents/UNC/Projects/CNALH Diversity/')

# Read in CNALH records
records_raw = read.csv('./Lower48_2014-10-08_specieslist/checklist.csv') #361502 records

# Read in list of lichen species from FIA
fia_splist = read.csv('../Lichen Traits/fia_lichen_master_list.csv')

# Read in list of accepted species names
splist = read.csv('./Lower48_2014-10-08_specieslist/checklist.csv')

target_fams = c('Parmeliaceae','Physciaceae','Cladoniaceae','Coccocarpiaceae','Collemataceae',
	'Lobariaceae','Caliciaceae','Pannariaceae','Lecideaceae','Nephromataceae','Peltigeraceae',
	'Placynthiaceae','Ramalinaceae','Sphaerophoraceae','Teloschistaceae','Umbilicariaceae',
	'Candelariaceae','Lecanoraceae')

target_genera = unique(c(fia_splist$Genus, fia_splist$syn.genus))
target_genera = target_genera[target_genera!='']

# Only look at Canada, USA, Mexico
records_raw$country = factor(records_raw$country)
levels(records_raw$country)
mycountries = c('Canada',"CANADA",'canada','Meixco','Mexcio','mexico',
	'Mexico','MEXICO',"México",'Meixco', 'Mexcio','Mwxico', 'U.S.A.',
	'U.S.A..','USA - HAWAII','United States','United States of America',
	'usa','USA','Usa','US','united States','United states','United STates')
records = subset(records_raw, country %in% mycountries) #360373 records

# Standardize capitalization of genera names
sapply(records$genus, function(x){
	cap.first(gsub(' ','',x))
}) -> records$genus

# Examine families
rec_fams = unique(records$family)
rec_fams[order(rec_fams)]

# Fix some misspelled families
records[records$family=="Brigantiaeaceae",'family']<- 'Brigantiaceae'
records[records$family %in% c("Chrysothrichaceae","Chrysotrichaceae"),'family']<- 'Chrysothricaceae'
records[records$family=="Haematommaceae",'family']<- "Haematommataceae" 
records[records$family=="Sphinctrinacaea",'family']<- "Sphinctrinaceae" 
records[records$family=="Stictidacaea",'family']<- "Stictidaceae" 

# Remove records not in targeted families
records = subset(records, family %in% target_fams) # 260853 records

# Examine genera
unique(records$genus)[order(unique(records$genus))]

# Fix mis-spelled genera
records[records$genus=="Tuckermanopsis",'genus']<- 'Tuckermannopsis'

# Records that do not match central taxonomy in CNALH
good_names = str_trim(unique(paste(str_trim(splist$genus), str_trim(splist$specificEpithet))))
record_names = get.latinbinom(records$scientificName)
records_bad = subset(records, !(record_names %in% good_names))
unique(records_bad$scientificName)

# Add species to good_names that clearly are correct in CNALH
add_names = c('Gyalolechia xanthostigmoidea', 'Fuscopannaria coralloidea', 'Umbilicaria proboscidea',
	'Buellia lacteoidea', 'Lecidea','Bacidina egenuloidea')
good_names = c(good_names, add_names)

# Fix mis-spelled or old species names in records
records[records$scientificName=='Buellia lepidastroidea','scientificName'] <- 'Buellia sequax'
records[records$scientificName=='Bacidia egenuloidea','scientificName'] <- 'Bacidina egenuloidea'
records[records$scientificName=='Caloplaca xanthostigmoidea', 'scientificName'] <- 'Gyalolechia xanthostigmoidea'
records[records$scientificName=='Lecidea insularis', 'scientificName'] <- 'Rimularia insularis'

## Remove records whose scientificNames are not in good_names (derived from genus, species in species list)
record_names = get.latinbinom(records$scientificName)
sum(record_names %in% good_names)
records = subset(records, record_names %in% good_names) # 257295

# Remove species not in targeted genera - may want to do this for targeted species instead
records = subset(records, genus %in% target_genera) # 181904 records

# Divide records into a set with coordinates and a set located to county level
records$decimalLatitude = as.numeric(records$decimalLatitude)
records$decimalLongitude = as.numeric(records$decimalLongitude)

# Add a column flagging whether coordinates changed
records$locEdit = F

# Assign coordinates to redacted localities using centroids of counties
records_county = subset(records, is.na(decimalLatitude)) #169 records
records_county[,c('country','stateProvince','county','locality')]
find_locs = unique(records_county[,c('country','stateProvince','county')])
find_locs$Lat = NA; find_locs$Lon = NA
find_locs[find_locs$county=='Edmonson',c('Lat','Lon')] <- c(37.202580, -86.218226)
find_locs[find_locs$county=='Avery',c('Lat','Lon')] <- c(36.083778, -81.916961)
find_locs[find_locs$county=='Haywood',c('Lat','Lon')] <- c(35.553605, -82.966312)
find_locs[find_locs$county=='Bladen',c('Lat','Lon')] <- c(34.602780, -78.539721)
find_locs[find_locs$county=='Carteret',c('Lat','Lon')] <- c(34.793599, -76.587769)
find_locs[find_locs$county=='Florence',c('Lat','Lon')] <- c(34.024970, -79.690412)
find_locs[find_locs$county=='Darlington',c('Lat','Lon')] <- c(34.329154, -79.952160)
find_locs[find_locs$county=='Wake',c('Lat','Lon')] <- c(35.793767, -78.666188)
find_locs[find_locs$county=='Teton',c('Lat','Lon')] <- c(43.948921, -110.580647)
find_locs[find_locs$county=='Park',c('Lat','Lon')] <- c(44.518208, -109.454910)
find_locs[find_locs$county=='Orange',c('Lat','Lon')] <- c(36.066422, -79.126231)
find_locs[find_locs$county=='Pierce',c('Lat','Lon')] <- c(47.033251, -122.159689)
find_locs[find_locs$county=='Chatham',c('Lat','Lon')] <- c(35.701793, -79.291376)
find_locs[find_locs$county=='Montgomery',c('Lat','Lon')] <- c(35.349764, -79.887426)

# Replace NAs in records data
# and indicate that record locality info was changed
for(i in 1:nrow(find_locs)){
	use_info = find_locs[i,]

	records[is.na(records$decimalLatitude)&(records$stateProvince==use_info$stateProvince)&(records$county==use_info$county), 'decimalLatitude'] <-  use_info$Lat
	records[is.na(records$decimalLongitude)&(records$stateProvince==use_info$stateProvince)&(records$county==use_info$county), 'decimalLongitude'] <-  use_info$Lon
	records[is.na(records$decimalLongitude)&(records$stateProvince==use_info$stateProvince)&(records$county==use_info$county), 'locEdit'] <- T
}

# Records to spatial data
records_geo = subset(records, !is.na(decimalLatitude)) #181900
records_ll = records_geo[,c('id','family','scientificName','genus','specificEpithet','taxonRank',
	'infraspecificEpithet','country','stateProvince','county','locality','decimalLatitude','decimalLongitude',
	'geodeticDatum')]
coordinates(records_ll) = c('decimalLongitude','decimalLatitude')
proj4string(records_ll) = CRS("+proj=longlat")

# Examine and remove records that fall in the water
NA_outline = readOGR('C:/Users/jrcoyle/Documents/UNC/GIS shape files/N Am Outline','na_base_Lambert_Azimuthal')
NA_buff = gBuffer(NA_outline, width=2000) # buffer extending approx. 2km from coast
NA_buff10 = gBuffer(NA_buff, width=8000)
save(NA_buff, NA_buff10, file='north_america_outline_buffers.RData')
plot(NA_buff)

# Remove records not in N. America

# Read in file with new coordinates (from previous analysis of Parm/Phys data)
bad_coords = read.csv('reset_coordinates.csv')
reset_coords = subset(bad_coords, !is.na(lat))

for(i in 1:nrow(reset_coords)){
	x = reset_coords[i,]
	
	use_recs = which(records_geo$country==x$country & records_geo$stateProvince==x$stateProvince & records_geo$county==x$county & records_geo$locality==x$locality)
	
	changed = (records_geo[use_recs,'decimalLongitude']==x$lon) * (records_geo[use_recs,'decimalLatitude']==x$lat)

	records_geo[use_recs,'decimalLongitude'] <- x$lon
	records_geo[use_recs,'decimalLatitude'] <- x$lat
	records_geo[use_recs,'locEdit'] <- changed==0
}

# Make into spatial data again
records_ll = records_geo[,c('id','family','scientificName','genus','specificEpithet','taxonRank',
	'infraspecificEpithet','country','stateProvince','county','locality','decimalLatitude','decimalLongitude',
	'geodeticDatum')]
coordinates(records_ll) = c('decimalLongitude','decimalLatitude')
proj4string(records_ll) = CRS("+proj=longlat")
records_sp = spTransform(records_ll, CRS(proj4string(NA_outline)))

par(mar=c(0,0,0,0))
plot(NA_buff10)
plot(records_sp, pch=15, cex=.3, add=T)

rec_poly = over(records_sp, NA_buff10)
water_recs = records_sp[is.na(rec_poly),]

par(mar=c(0,0,0,0))
plot(NA_buff10)
plot(records_sp, pch=15, cex=.3, add=T)
plot(water_recs, pch=16, cex=.5, col=2, add=T)

water_places = unique(water_recs[,c('country','stateProvince','county')]@data)
water_places[order(water_places$country, water_places$stateProvince),]
water_places$country = as.character(water_places$country)

# Find water_places that have not already been scanned in reset_coords (from Parm/Phys)
reset_add = unique(water_recs[,c("country", "stateProvince", "county", "locality")]@data)
reset = rbind(reset_coords[,1:4], reset_add)
new_places = duplicated(reset)[(nrow(reset_coords)+1):nrow(reset)]==F
reset_add = reset_add[new_places,]
reset_add$lat = NA
reset_add$lon = NA

#write.csv(reset_add, 'reset_coordinates_fiagenera.csv', row.names=F)

pdf('water_records_maps_fiagenera.pdf', height=5, width=5)
for(i in 1:nrow(water_places)){
	par(mar=c(0,0,0,0))
	map('world',c('canada','usa','mexico','hawaii'), xlim=c(-180,-50), ylim=c(0,90))
	abline(v=seq(-180,-50,5), col='grey')
	abline(h=seq(0,90,5), col='grey')
	axis(1)
	axis(2)	
	use_place = water_places[i,]
	use_recs = subset(water_recs, country==use_place$country & stateProvince==use_place$stateProvince & county==use_place$county)$id
	plot(subset(records_ll, id %in% use_recs), pch=16, col=2, add=T)
	text(-175,5,labels=paste(use_place, collapse=' '), pos=4)
	text(-50, 85, pos=2, cex=0.7, labels=paste(subset(records_ll, id %in% use_recs)$locality, collapse='\n'))
}
dev.off()

# Manually check locations against Google Maps
i=86
use_place = water_places[i,]
use_recs = subset(water_recs, country==use_place$country & stateProvince==use_place$stateProvince & county==use_place$county)$id
subset(records_ll, id %in% use_recs)
unique(coordinates(subset(records_ll, id %in% use_recs)[,c('stateProvince','locality')]))

# Read back in edited coordinates
reset_add = read.csv('reset_coordinates_fiagenera.csv')
reset_coords2 = subset(reset_add, !is.na(lat))

for(i in 1:nrow(reset_coords2)){
	x = reset_coords2[i,]
	
	use_recs = which(records_geo$country==x$country & records_geo$stateProvince==x$stateProvince & records_geo$county==x$county & records_geo$locality==x$locality)
	
	changed = (records_geo[use_recs,'decimalLongitude']==x$lon) * (records_geo[use_recs,'decimalLatitude']==x$lat)

	records_geo[use_recs,'decimalLongitude'] <- x$lon
	records_geo[use_recs,'decimalLatitude'] <- x$lat
	records_geo[use_recs,'locEdit'] <- changed==0
}

# Write out record data
write.csv(records_geo, '../FIA Lichen/Data/CNALH_records_fia_genera_NAm_2014-10-08.csv', row.names=F)

#############################################################
### Calculate Regional Richness, All Species and Parmeliaceae/Physciaceae

parmphys = read.csv('Parm_Phys_records_2014-09-20.csv')
allsp = read.csv('../FIA Lichen/Data/CNALH_records_fia_genera_NAm_2014-10-08.csv')
parmphys_fia = subset(allsp, family %in% c('Parmeliaceae','Physciaceae'))

# Calculate number of species in records
species = unique(allsp[,'scientificName'])
species.binom = unique(get.latinbinom(species)) # 1757 species (including just genus)
genera = unique(sapply(strsplit(species.binom, ' '), function(x) x[1])) #80 genera

# Calculate how many species in each family in records
species_f = unique(allsp[,c('family','scientificName','genus')])
genera_f = unique(species_f[,c('family','genus')])
binom_f = unique(data.frame(family=allsp$family, latin_binom = get.latinbinom(allsp$scientificName)))

table(genera_f$family) # 44 Parmeliaceae, 9 Physciaceae
table(binom_f$family) # 868 Parmeliaceae, 227 Physciaceae

## Calculate regional richness for FIA plots using rarefaction

# Subset to the family of interest
parm = subset(parmphys_fia, family=='Parmeliaceae')
phys = subset(parmphys_fia, family=='Physciaceae')

# Make records into spatial data
allsp_sp = allsp
parm_sp = parm
phys_sp = phys
coordinates(allsp_sp) = c('decimalLongitude','decimalLatitude')
proj4string(allsp_sp) = CRS("+proj=longlat")
coordinates(parm_sp) = c('decimalLongitude','decimalLatitude')
proj4string(parm_sp) = CRS("+proj=longlat")
coordinates(phys_sp) = c('decimalLongitude','decimalLatitude')
proj4string(phys_sp) = CRS("+proj=longlat")


# Only calculate regional richness for plots used in models
# This will exclude AK plots where I did not download records for
fia_geo = subset(master_locs, yrplot.id %in% rownames(model_data))

# Calculate for FIA plots with spatial coordinates
coordinates(fia_geo) = c('LON','LAT')
proj4string(fia_geo) = CRS("+proj=longlat")

# Determine which records will be used to calculate regional richness
recs = allsp_sp

# Calculate the number of records within 500 km buffer (Could easily check other buffers)
dist500 = c()
for(i in 1:nrow(fia_geo)){
	dist500 = c(dist500, nrow(find_recs(fia_geo[i,], recs, 500)))
} # min=2622 all sp, 882 parm, 374 phys
names(dist500) = fia_geo$yrplot.id
min(dist500) # Min num records for all FIA plots 

recordnum = data.frame(parm = dist500)
recordnum$phys = dist500
write.csv(recordnum, '../FIA Lichen/Data/Regional Richness/num_records_CNALH_2014-10-08_fiagenera.csv', row.names=F)

## Calculate Regional Richness using rarefaction (vegan) only for plots used in analysis
regS = data.frame(yrplot.id=rownames(model_data))

# Define size of subsample (nsamps)
nsamp = 350 # use 350 for regParm and regPhys, use 2500 for allsp

recs = phys_sp

richness = c()

for(i in 1:nrow(model_data)){
	usepoint = fia_geo[fia_geo$yrplot.id==rownames(model_data)[i],]
	userecs = find_recs(usepoint, recs, 500)
	comm = tab_species(userecs$scientificName)
	rich = rarefy(comm, nsamp)
	richness = c(richness, as.numeric(rich))
}


# Compare richness with sampling effort (# records)
png('../FIA Lichen/Figures/ regS (500 record rarefaction) vs sampling effort.png', height=600, width=600)
par(mar=c(5,6,1,1.5))
plot(richness~dist500[rownames(model_data)], las=1, xlab='Number of records w/in 500km radius', ylab='Rarefied richness for 350 records', cex=2, cex.lab=1.5, cex.axis=1.5)
dev.off()


# Save data
regS$regS = richness
regS$regPhys = richness
regS$regParm = richness
write.csv(regS, '../FIA Lichen/Data/Regional Richness/fia_lichen_reg_richness_CNALH-2014-09-20.csv', row.names=F)


#
fia_geo2 = merge(fia_geo, regS)

# Quick Plot 
spplot(fia_geo2, 'regPhys', col.regions=colorRampPalette(c('dark blue','blue','green','yellow','orange','red'))(10), cuts=10)



#######################################################################3
## Assess how sampling radius (km) affects richness estimates
## see script 'calc_regional_richnee_multi-radii.R' for actual calculation run on cluster

# Load R objects with richness and number of observations computed at different scales and number of sampled records
load('./data/Regional Richness/regS_across_scales.Rdata')

# Objects are regS and nobs, scales
# Define scales over which richness was calculated
scales = c(50, 100, 250, 400, 500)
nsamps = c(25, 50, 100, 200, 250, 500, 750, 1000, 1500, 2000, 2500, 3000)

# Calculate % of plots with nsamp at each scale
pcts = array(NA, dim=c(length(nsamps), length(scales)), dimnames=list(nsamp=nsamps, scale=paste('R',scales, sep='')))
for(i in nsamps){
for(j in colnames(pcts)){
	pcts[as.character(i),j] = sum(nobs[,j] >= i)/dim(regS)[1]	
}}

# Format and save table
pcts_tab = apply(format(pcts*100, digits=1), c(1,2), function(x) paste(x, '%', sep=''))
write.table(pcts_tab, './Figures/pct plots with different num records across scales.txt', sep='\t', row.names=T, quote=F)

# Plot example rarefaction curves from two different regions
most_samps = rownames(nobs)[order(nobs[,'R400'], decreasing=T)]
master[most_samps,c('yrplot.id','state.abbr')]
# use 1998_55_3_9514 (from WI) and 2007_4_1_89307 (form AZ) and 1999_53_37_10205 (from WA) and 1998_50_21_7453 (from VT)

ex_plots = c('1998_55_3_9514','2007_4_1_89307','1999_53_37_10205','1998_50_21_7453')
names(ex_plots) = master[ex_plots,'state.abbr']
nobs[ex_plots,]

use_regS = regS[,'R400',]

pdf('./Figures/regional richness rarefaction.pdf', height=5, width=5)
par(mar=c(4,4,1,1))
plot(nsamps, use_regS[ex_plots[1],], type='n', las=1, 
	xlab='# Records Sampled', ylab='Species Richness', ylim=c(0,max(use_regS[ex_plots,])))
for(i in 1:4) lines(nsamps, use_regS[ex_plots[i],], lty=i, lwd=3)
rect(2400, 0, 2600, 430, col='white', border='white')
text(2500, use_regS[ex_plots,'2500'], labels=names(ex_plots))
dev.off()

nobs[ex_plots, ]

# Calculate correlations between richness calculated at different scales
# Compare 100 samples at 250km and 2500 samples at 500km

cor(regS[,'R250','500'], regS[,'R500','2500'], method='spearman', use='complete.obs')
cor(regS[,'R250','500'], regS[,'R500','2500'], method='pearson', use='complete.obs')

cor(regS[,'R250','100'], regS[,'R500','2500'], method='spearman', use='complete.obs')
cor(regS[,'R250','100'], regS[,'R500','2500'], method='pearson', use='complete.obs')

## Map regional richness using 250km sampling radius and 500 subsamples

# Colors
ncuts=10
mycol = read.csv('C:/Users/jrcoyle/Documents/UNC/Projects/blue2red_10colramp.txt')
mycol = apply(mycol,1,function(x) rgb(x[1],x[2],x[3],maxColorValue=256))
mycol = mycol[10:1]
mycolramp = colorRampPalette(mycol)(ncuts)
mycolrampbw = colorRampPalette(c('grey80','black'))(ncuts)

# Map projection
plot_prj = paste("+proj=laea +lat_0=40 +lon_0=-97 +units=km",sep='')

# N. Am. outline
OUTLINES = readOGR('../../GIS shape files/N Am Outline','na_base_Lambert_Azimuthal')
OUTLINES.laea = spTransform(OUTLINES,CRS(plot_prj))

# Make spatial points data frame with only 1923 plots used in models
reg500 = regS[,'R500','2500']
reg250 = regS[,'R250','500']
sp_data = data.frame(master[dimnames(regS)$yrplot.id,c('LON','LAT')], reg500, reg250)
coordinates(sp_data) = c('LON','LAT')
proj4string(sp_data) = CRS("+proj=longlat")
sp_data = spTransform(sp_data, CRS(plot_prj))

range(sp_data$reg250, na.rm=T)

pdf('./Figures/Maps/Map regional lichen richness R250-S500.pdf', height=8, width=12)
trellis.par.set(axis.line=list(col=NA))
colcuts = seq(90, 210, 12) #seq(230,430,20)  
spplot(sp_data, 'reg250', ylim=c(-1600,1500), main='', panel=function(x,y,subscripts,...){
		sp.polygons(OUTLINES.laea, fill='white')
		panel.pointsplot(x,y,...)
	}, cuts=colcuts, cex=1.5, col.regions = mycolramp, auto.key=F,
	key=list(x=1,y=.3, corner=c(1,.5), title='Species\nRichness',
		rectangles=list(col=mycolramp, size=3, border='transparent'),
		text=list(c(paste(colcuts[1],colcuts[2], sep='-'), paste((colcuts+1)[2:(ncuts)], colcuts[2:ncuts+1], sep='-')))) 
)

dev.off()

#################################################
## Number of species in CNALH data set

cnalh_sp = get_species(records$scientificName)

unique(sapply(strsplit(cnalh_sp, ' '), function(x) x[1]))


### Plot map of county-level and overlayed plot-level regional richness
county.sh$uniqueID = paste(county.sh$STATE_NAME, county.sh$NAME)
rownames(county_loc) = paste(county_loc$STATE_NAME, county_loc$COUNTYNM)

# Make county polygons w/ regS
fia_poly = county.sh[county.sh$uniqueID %in% rownames(county_loc),]
fia_poly$regS = county_loc[fia_poly$uniqueID,'regS']

# Change to laea projection
mybbox = bbox(fia_geo)
plot_prj = paste("+proj=laea +lat_0=",mean(mybbox[2,])," +lon_0=",mean(mybbox[1,])," +units=km",sep='')
plot_prj = paste("+proj=laea +lat_0=40 +lon_0=-97 +units=km",sep='')

fia_poly = spTransform(fia_poly, CRS(plot_prj))
fia_geo = spTransform(fia_geo, CRS(plot_prj))

# N. Am. outline
OUTLINES = readOGR('../../GIS shape files/N Am Outline','na_base_Lambert_Azimuthal')
OUTLINES.laea = spTransform(OUTLINES,CRS(plot_prj))

# Define graphical parameters
mycols = c('dark blue','blue','green','yellow','orange','red')
ncuts = 10
mycolramp = colorRampPalette(mycols)(ncuts)
colcuts = seq(min(master$regS), max(master$regS), length.out=ncuts)

png('./Figures/lichen regional richness for fia plots and counties.png', height=1000, width=1200)
spplot(fia_geo, 'regS', ylim=c(-1600,3400), panel=function(x,y,...){
	sp.polygons(OUTLINES.laea, fill='grey80')
	sp.polygons(fia_poly, col='black', fill=mycolramp[cut(fia_poly$regS, colcuts, include.lowest=T)])
	panel.pointsplot(x,y,...)
},key.space='right', col.regions = paste(mycolramp, '98', sep=''), cuts=colcuts, cex=1.5,
colorkey=list(text=list(cex=3), points=list(cex=2))
)
dev.off()

##########################################
### Calculate binned species richness map

setwd('C:/Users/jrcoyle/Documents/UNC/Projects/CNALH Diversity')
#bindata = read.csv('Richness all FIA genera 5 degree bins.csv')

# records does not contain 58 records w/o coordinates.  May want to change this.
records = read.csv('Parm_Phys_records_2014-09-20.csv')

# Make Lat/lon bins
binwidth = 2
mylons = seq(-167, -50, binwidth)
mylats = seq(15, 80, binwidth)

records$lonbin = cut(records$decimalLongitude, breaks=mylons)
records$latbin = cut(records$decimalLatitude, breaks=mylats)

# Decide which group to use
recs = subset(records, family=='Physciaceae')

# Each cell is a list of the community occuring there
species_distribution = tapply(recs[,c('scientificName')], list(recs$latbin, recs$lonbin), FUN=tab_species)

# Calculate richness
richness_raw = apply(species_distribution, c(1,2), function(x) length(unlist(x)))
class(richness_raw) = "table"

bindata = as.data.frame(richness_raw)
names(bindata) = c('latbin','lonbin','richness.raw')
bindata$lat = rep(mylats[1:(length(mylats)-1)],length(mylons)-1)
bindata$lon = rep(mylons[1:(length(mylons)-1)],each=length(mylats)-1)


# Plot on map

#mycolors = colorRampPalette(c("white","yellow","red", "darkred"))(7)

mycolors = c('#FFFFFF',colorRampPalette(c('#2c7bb6','#abd9e9','#fdae61','#ca373b'))(7))
mycolorsT = paste(mycolors, "80", sep='')

colorvar = cut(bindata$richness.raw, breaks=c(0,seq(1,max(bindata$richness.raw), length.out=8)), right=F, include.lowest=T)

tiff('physciaceae_species_richness_uncorrected_map_2-deg_CNALH-2014-09-20.tif', width=800, height=800)
map('world', c('canada','usa','mexico'),col='grey50', xlim=c(-180,-45))
points(decimalLatitude~decimalLongitude, data=records, pch='.', cex=1)
rect(bindata$lon,bindata$lat,bindata$lon+binwidth,bindata$lat+binwidth,
	col=mycolorsT[colorvar], border='transparent')

legend('bottomleft',levels(colorvar)[2:8], pch=15, col=mycolorsT[2:8], bty='n', cex=2)
title(main='Physciaceae Species Richness\n(Uncorrected)', cex.main=2)
dev.off()


# Tally number of samples in each bin
samples = xtabs(~latbin+lonbin, data=recs)
samples = as.data.frame(samples)
names(samples) = c('latbin','lonbin','n.obs')

hist(samples$n.obs[samples$n.obs>0], breaks=20)

bindata = merge(bindata, samples)


# Calculate richness using rarefaction only for samples with at least 100 records
nsamp = c(25, 50, 100, 200)

richness_rarefy = sapply(nsamp, function(n){

apply(bindata, 1, function(x){
	latbin = x['latbin']
	lonbin = x['lonbin']
	
	comm = unlist(species_distribution[latbin,lonbin])

	if(length(comm)>0){
		as.numeric(rarefy(comm, n))
	} else 0
})

})

colnames(richness_rarefy) = paste('rich', nsamp, sep='')

bindata = cbind(bindata, richness_rarefy)

head(bindata)

write.csv(phys2_bin, 'Phys_richness_2-deg_CNALH-2014-09-20.csv', row.names=F)
phys5_bin = bindata
phys2_bin = bindata
parm5_bin = bindata
parm2_bin = bindata


# Map subsampled richness
bindata = parm2_bin

# Make sure bins are correct size

tiff('parmeliaceae_species_richness_samp50_map.tif', width=800, height=800)

bordervar = bindata$n.obs>=50
colorvar = cut(bindata$rich50, breaks=c(0,seq(min(bindata$rich50[bordervar]),max(bindata$rich50[bordervar]), length.out=8)), right=F, include.lowest=T)
colorvar[!bordervar]<-NA

map('world',c('canada','usa','mexico'), col='grey50', xlim=c(-180,-50))
points(decimalLatitude~decimalLongitude, data=records, pch='.', cex=3, col='grey50')
rect(bindata$lon,bindata$lat,bindata$lon+binwidth,bindata$lat+binwidth,
	col=mycolorsT[colorvar], border="transparent")
legend('bottomleft',levels(colorvar)[2:8], pch=15, col=mycolorsT[2:8], bty='n', cex=2)
title(main='Parmeliaceae Species Richness\n(50 records)', cex.main=2)
dev.off()


## NOT DONE YET ##
# Map Number of observations
tiff('fia_genera_nsamps_map.tif', width=800, height=800)
png('parmeliaceae_nsamps_map.png', height=1200, width=1200, bg='transparent')
colorvar = cut(bindata$n.obs, breaks=c(0,10,100,500,1000,2000,5000,10000), right=F, include.lowest=T)

map('world',c('canada','usa','mexico'), col='grey50', xlim=c(-180,-50))
points(decimalLatitude~decimalLongitude, data=records, pch='.', cex=3)
rect(bindata$lon,bindata$lat,bindata$lon+binwidth,bindata$lat+binwidth,
	col=mycolorsT[colorvar], border="transparent")
legend('bottomleft',c("<10","[10,100)","[100,500)","[500,1000)","[1000,2000)","[2000,5000)",">5000"), pch=15, col=mycolorsT, bty='n', cex=2)
title(main='Number of Records (All FIA Genera)', cex.main=2)
dev.off()

# Map plain observation locations
parmelia.sp = parmelia
coordinates(parmelia.sp) = c("decimalLongitude","decimalLatitude")
proj4string(parmelia.sp) = CRS("+proj=longlat +ellps=WGS84")
mybbox = bbox(parmelia.sp)

midpoint = apply(mybbox, 1, mean)
myprojstring = paste("+proj=laea +lat_=",midpoint[2]," +lon_0=",midpoint[1],"+x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=km +no_defs", sep='')
parmelia.laea = spTransform(parmelia.sp, CRS(myprojstring))

na.outlines = readOGR('C:/Users/jrcoyle/Documents/UNC/GIS shape files/N Am Outline','na_base_Lambert_Azimuthal')
na.outlines = spTransform(na.outlines,CRS(myprojstring))
mybbox.laea = bbox(parmelia.laea)

png('parmeliaceae_obslocs_map.png', width=1200, height=1200, bg="transparent")
plot(na.outlines,xlim=mybbox.laea[1,], ylim=mybbox.laea[2,], col='white', bg='transparent', border='grey50', lwd=2)
points(parmelia.laea,pch='.', col=rgb(0,0,.4), cex=2)
dev.off()


######################################################
### Continuous richness maps based on records within 500km radii
library(maps)
library(maptools)

# Map projection - equal area
plot_prj = paste("+proj=laea +lat_0=40 +lon_0=-97 +units=km",sep='')

# Read records back in as a dataframe
records_sp = read.csv('Parm_Phys_records_2014-09-20.csv') # this is for parmeliaceae and physciaceae
records_sp = read.csv() # this is for all fia genera

# Make records into spatial data
coordinates(records_sp) = c('decimalLongitude','decimalLatitude')
proj4string(records_sp) = CRS("+proj=longlat")

# Make map of USA 
usa = map('usa')
usa_sp = map2SpatialPolygons(usa, IDs=usa$names, proj4string=CRS("+proj=longlat +datum=WGS84"))
usa_laea = spTransform(usa_sp, CRS(plot_prj))

# Map for north america
nam_outline = readOGR('C:/Users/jrcoyle/Documents/UNC/GIS shape files/N Am Outline','na_base_Lambert_Azimuthal')
nam_outline = subset(nam_outline, COUNTRY %in% c('USA','MEXICO','CANADA'))
nam_laea = spTransform(nam_outline,CRS(plot_prj))

# Generate grid of points over which to sample
mybox = bbox(usa_laea)
mybox = bbox(nam_laea)

x = seq(mybox[1,1]+50, mybox[1,2], 100)
y = seq(mybox[2,1]+50, mybox[2,2], 100)

mygrid = expand.grid(x,y)
plot(nam_laea)
points(mygrid)



# Exclude grid not in North America
nam_poly = SpatialPolygons(nam_laea@polygons, proj4string=CRS(plot_prj))
nam_poly1 = unionSpatialPolygons(nam_poly, IDs = rep(1,length(nam_poly)))

mygrid_sp = SpatialPoints(mygrid, proj4string=CRS(plot_prj))
usegrid = over(mygrid_sp, nam_poly1)

mygrid_sp = mygrid_sp[!is.na(usegrid),]

plot(nam_poly1)
plot(mygrid_sp, add=T)


# Covert back to lat lon so that find_recs() works properly
mygrid_ll = spTransform(mygrid_sp, CRS("+proj=longlat +datum=WGS84"))

# Define data frame where going to store data
mygrid_data = data.frame(coordinates(mygrid_ll))
names(mygrid_data) = c('lon','lat')

# Define group for which to calculate richness
recs = subset(records_sp, family=='Physciaceae') # Parm=89487 Phys=45918
#recs = subset(records_sp, family=='Parmeliaceae' & genus %in% target_genera)

# Calculate number of records with given radius of each grid point
# will use this for determining number of records for rarefaction when calculating richness
rad = 250 # 500 km radius]

nobs = sapply(1:length(mygrid_ll), function(i){
	userecs = find_recs(mygrid_ll[i,], recs, rad)
	nrow(userecs)
})

mygrid_data$nobsParm = nobs
mygrid_data$nobsPhys = nobs

# Examine maps where pixels are colored by whether number of records is below a certain threshold
nsamp = 50
plot(mygrid_sp, col=as.numeric(nobs<nsamp)+1, pch=15, cex=.5)

## Calculate rarefied regional species richness across grid

# Define number of records to use in rarefaction 
nsamp = c(25, 50, 100, 150, 200, 300, 500)

richness = data.frame()
for(i in 1:length(mygrid_ll)){
	userecs = find_recs(mygrid_ll[i,], recs, 500) # 500km radius

	sapply(nsamp, function(n){
		if(nrow(userecs)>=n){
			comm = tab_species(userecs$scientificName)
			as.numeric(rarefy(comm, n))
		} else {
			rich = NA
		}
	}) ->  rich

	richness = rbind(richness, rich)
}

names(richness) = paste('regParm',nsamp, sep='') # this changes depending on the taxonomic group
names(richness) = paste('regPhys',nsamp, sep='')

mygrid_data = cbind(mygrid_data, richness)

#write.csv(mygrid_data, 'grid_map_Parm_Phys_500km_regS.csv', row.names=F)
write.csv(mygrid_data, 'grid_map_Parm_Phys_250km_regS.csv', row.names=F)


# Plot
mybreaks = seq(40, 130, by=10)
colorvar = cut(mygrid_data$regParm300, breaks=mybreaks)
usecolors = colorRampPalette(mycolors[2:8])(length(mybreaks-1))

png('parmeliaceae_richness_subsample300_250km.png', height=700, width=700)
par(mar=c(0,4,0,7))
plot(nam_poly1)
plot(mygrid_sp, col=usecolors[colorvar], add=T, pch=15, cex=.8)
mtext('Parmeliaceae Species Richness\n(from 300 records w/in 250km radius)', 3, -6, cex=1.5)
usr = par('usr')
par(xpd=F)
xwide = (usr[2]-usr[1])/30
ypad = (usr[4]-usr[3])/3
plotColorRamp(usecolors, length(mybreaks)-1, c(usr[2], usr[3]+ypad, usr[2]+xwide, usr[4]-ypad), 
	labels=mybreaks, title='Num. Species', mycex=1, uneven.lab = F, labrange = NA, ndig=1)

dev.off()

rect()


windows()
par(mar=c(0,0,0,0))

plot(0,0, xlim=c(0,1), ylim=c(-.2,1.2), type='n')
rect(rep(0,8), seq(0, .7, .1), rep(.1, 8), seq(.1, .8, .1), col=mycolors)
text(rep(0.1, 9), seq(0,.9,.1), labels=seq(40, 140, length.out=9), pos=4, cex=2)

par('usr')



### CODE BELOW NOT YET CORRECT FOR THIS SRCIPT



 Calc  regS for each grid cell 

# Define size of subsample (nsamps) and number of times to resample (reps)
nsamp = 350


regS = c()

for(i in 1:length(mygrid_ll)){
	userecs = find_recs(mygrid_ll[i,], records_sp, 500)

	
	regS = c(regS, calc_rich_samp(userecs, nsamp, reps))
}


regS_sp = SpatialPointsDataFrame(mygrid_sp, data.frame(regS))



#### Map


# Color Ramp
ncuts=10

mycol = read.csv('C:/Users/jrcoyle/Documents/UNC/Projects/blue2red_10colramp.txt')
mycol = apply(mycol,1,function(x) rgb(x[1],x[2],x[3],maxColorValue=256))
mycol = mycol[10:1]
mycolramp = colorRampPalette(mycol)(ncuts)


usecol = mycolramp[cut(regS, c(100,110,120,130,140,150,160,170,180,190,204))]

png('./Data/Regional Richness/regS_lichen_grid.png', height=800, width=1200)
par(mar=c(0,0,0,0))
plot(regS_sp, col=usecol, pch=15, cex=4) 
dev.off()

png('./Data/Regional Richness/regS_lichen_grid_key.png', height=400, width=200)
par(mar=c(0,0,0,0))
plot(0,0, xlim=c(0,1), ylim=c(-.2,1.2), type='n')
rect(rep(0,10), seq(0, .9, .1), rep(.1, 10), seq(.1, 1, .1), col=mycolramp)
text(rep(0.1, 11), seq(0,1,.1), labels=seq(100,200,10), pos=4, cex=2)
dev.off()


save.image('./Data/Regional Richness/lichen_regS_grid.RData')






######################################################
### Species Ranges

# Map range of one species

psulcata  = subset(parmelia, scientificName=="Parmelia sulcata")
apallidula = subset(parmelia, scientificName=="Ahtiana pallidula")
aaurescens = subset(parmelia, scientificName=="Ahtiana aurescens")
asphaerosporella = subset(parmelia, scientificName=="Ahtiana sphaerosporella")

ahtiana = rbind(apallidula, aaurescens, asphaerosporella)
coordinates(ahtiana) = c("decimalLongitude","decimalLatitude")
proj4string(ahtiana) = CRS("+proj=longlat +ellps=WGS84")
mybbox = bbox(ahtiana)

midpoint = apply(mybbox, 1, mean)
myprojstring = paste("+proj=laea +lat_=",midpoint[2]," +lon_0=",midpoint[1],"+x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=km +no_defs", sep='')
ahtiana.laea = spTransform(ahtiana, CRS(myprojstring))

na.outlines = readOGR('C:/Users/jrcoyle/Documents/UNC/GIS shape files/N Am Outline','na_base_Lambert_Azimuthal')
na.outlines = spTransform(na.outlines,CRS(myprojstring))
mybbox.laea = bbox(ahtiana.laea)

load('ahtiana map.png', height=)
ap = bio.stack[[12]]
test.ext = projectExtent(ap, CRS(myprojstring))
ap.laea = projectRaster(from=ap, to=test.ext, crs=CRS(myprojstring))

mycols = c('black','red','blue')
mypch = c(1,2,0)

png('ahtiana_map.png', width=800, height=400, bg=rgb(1,1,1,.8))
par(mar=c(1,1,1,1))
par(oma=c(0,0,0,4))
plot(ap.laea, axes=F, col=colorRampPalette(c('darkorange3','yellow','green','cyan2','dodgerblue','deepskyblue3','darkblue'))(50),
	ext=mybbox.laea+c(-200,-200,200,200), legend.lab="Annual Precipitation (mm)")
plot(na.outlines,xlim=mybbox.laea[1,], ylim=mybbox.laea[2,], col='transparent', bg='transparent', border='grey50', lwd=2, add=T)
points(ahtiana.laea,lwd=2,pch = mypch[factor(ahtiana.laea$scientificName)], col=mycols[factor(ahtiana.laea$scientificName)], cex=3)

legend('topright', levels(factor(ahtiana.laea$scientificName)), col=mycols, pch=mypch, pt.cex=3, cex=2, pt.lwd=2,box.col='transparent', bg=rgb(1,1,1,.8))
dev.off()
sps = unique(parmelia$scientificName)

sps[order(sps)]





############################################################################
### Functions

# A functions that returns Genus species for a scientificName that may contain subspecific taxonomic info
get.latinbinom = function(x){
	namelist = strsplit(str_trim(x), ' ')
	genus = str_trim(sapply(namelist, function(y) y[1]))
	species = str_trim(sapply(namelist, function(y) ifelse(is.na(y[2]), '',y[2])))
	str_trim(paste(genus, species, sep=' '))
}


# A function that capitalizes the first letter of a word
cap.first = function(x){
	lastpart = substr(x,2,nchar(x))
	firstletter = substr(x,1,1)
	x=paste(toupper(firstletter),tolower(lastpart), sep='')

	return(x)
}



# A function that finds all records within a given distance buffer from a points
find_recs = function(point, recs, buff){
	dists = spDistsN1(recs, point, longlat=T) # In km
	userecs = recs[dists<buff,]

	return(userecs)
}

# A function that returns a species list from a vector of species names
# It removes genus spp. that are already represented by a species
get_species = function(x){
	new_x = get.latinbinom(x)	
	splitnames = strsplit(new_x, " ")
	justgenus = unique(new_x[as.numeric(lapply(splitnames, length))==1])
	allspecies = new_x[as.numeric(lapply(splitnames, length))>1]
	allspecies = unique(allspecies)
	spgenera = unique(sapply(strsplit(allspecies, " "), function(y) y[1]))
	keepgenera = justgenus[!(justgenus %in% spgenera)]

	c(allspecies,keepgenera)
}

# A function that returns a vector counting the number of records for each species binomial
tab_species = function(x){
	if(length(x)>0){

	new_x = get.latinbinom(x)
	splitnames = strsplit(new_x, " ")
	justgenus = new_x[as.numeric(lapply(splitnames, length))==1]
	allspecies = new_x[as.numeric(lapply(splitnames, length))>1]

	if(length(allspecies)>0){
		spgenera = unique(sapply(strsplit(allspecies, " "), function(y) y[1]))
		keepgenera = justgenus[!(justgenus %in% spgenera)]
		allspecies = c(allspecies, keepgenera)
	} else {
		allspecies = justgenus
	}

	as.matrix(table(allspecies))
	} else NULL	
}


# A function that calculates average species richness of a sample of records
calc_rich_samp = function(recs, nsamp, reps){

	if(nrow(recs)>nsamp){

	mysamples = sapply(1:reps, function(x){
		species_list = get_species(recs[sample(1:nrow(recs), nsamp, replace=F),]$scientificName)
		length(species_list)
	})

	return(median(mysamples))

	} else {

	return(length(get_species(recs$scientificName)))	

	}
}

# A function that plots a vertical color ramp on the side of a plot
# cols    : the colors to use
# n       : number of divisions
# barends : location of whole bar c(xleft, ybottom, xright, ytop)
# labels    : vector of labels for bar, assumes 1st and last numbers correspond to 1st and last colors
# title   : title to print above bar
# mycex   : size of label and title text
# uneven.lab : TRUE when numeric labels are not evenly spaced along length of color bar
# labrange : use when uneven.lab=T to specify minimum and maximum values of color range c(min, max)
# ndig    : number of digits beyond decimal placeto print for numeric labels 
plotColorRamp = function(cols, n, barends, labels=NA, title=NA, mycex=1.5, uneven.lab = F, labrange = NA, ndig=1){
	dX = barends[3] - barends[1]
	dY = barends[4] - barends[2]
	dy = dY/n
	
	xpd.old = par('xpd')
	par(xpd=T)

	lend.old = par('lend')
	par(lend=1)

	usecols = colorRampPalette(cols)(n)

	for(i in 1:n){
		rect(barends[1], barends[2]+dy*(i-1), barends[3], barends[2]+dy*i, col=usecols[i], border=usecols[i])
	}

	if(!is.na(labels)){
		if(is.numeric(labels)){
			labels.round = format(round(labels, ndig), nsmall=ndig, trim=F)

			if(uneven.lab){
				dz = dY/diff(labrange)
				Yposition = barends[2] + dz*(labels-labrange[1])
			} else {
				dZ = labels[length(labels)]-labels[1]
				dz = dY/dZ
				Yposition = barends[2] + dz*(labels-labels[1])
			}

			text(barends[3]+dX*0.5, Yposition, labels.round, pos=4, cex=mycex)

		} else {
			labels.round = labels
			dz = dY/length(labels)
			Yposition = barends[2] + dz*(0:(length(labels)-1))

			text(barends[3]+dX*0.5, Yposition, labels, pos=4, cex=mycex)
		}

		segments(barends[3], Yposition, barends[3]+dX*0.5, Yposition)	
	}
	if(!is.na(title)){
		
		
		## Determine how many characters away to place title
		digits = max(nchar(labels.round)) # Maximum number of digits in a label
		largest = labels.round[which(nchar(labels.round)==digits)] # Which labels are longest
		
		small.chars = grep('[-.]', largest) # Does the largest label have a small character?
			if(length(small.chars)==length(largest)) digits = digits-0.6 # Discount the size of the largest label by 0.6 a character
		
		text(barends[3]+dX*0.5+par('cxy')[1]*mycex*(digits+.5), barends[2]+0.5*dY, labels=title, srt=-90, cex=mycex)
	}
	par(xpd=xpd.old)
	par(lend=lend.old)
}
