# This script diversity of lichens from occurence records downloaded from
# Consortium of North American Lichen Herbaria on 2-19-2012.

library(sp)
library(rgdal)
library(vegan) # rarefy

setwd('C:/Users/jrcoyle/Documents/UNC/Projects/CNALH Diversity/')
#setwd('C:/Users/jrcoyle/Documents/Projects/CNALH Diversity/')
options(stringsAsFactors=F)

# Read in CNALH records
records_raw = read.csv('./Lower48_2014-10-08/occurrences.csv') #355119 records

# Read in master data of FIA plots with plot data and environmental data
master = read.csv('../FIA Lichen/Data/master_data.csv')
master_locs = read.csv('../FIA Lichen/Data/fia_lichen_plot_locations.csv')

## NEED TO LOAD MODEL_DATA 


# Read in list of lichen species from FIA
fia_splist = read.csv('../Lichen Traits/fia_lichen_master_list.csv')

# Read in list of accepted species names
cnalh_splist = read.csv('CNALH_allrecords_NAm_2013-02-19_checklist.csv')

target_fams = c('Parmeliaceae','Physciaceae','Cladoniaceae','Coccocarpiaceae','Collemataceae',
	'Lobariaceae','Caliciaceae','Pannariaceae','Lecideaceae','Nephromataceae','Peltigeraceae',
	'Placynthiaceae','Ramalinaceae','Sphaerophoraceae','Teloschistaceae','Umbilicariaceae',
	'Candelariaceae','Lecanoraceae')

target_genera = unique(c(fia_splist$Genus, fia_splist$syn.genus))
target_genera = target_genera[target_genera!='']

# Only look at Canada, USA, Mexico
records_raw$country = factor(records_raw$country)
levels(records_raw$country)
mycountries = c('Canada',"CANADA",'canada','Meixco','Mexcio','mexico','Mexico','MEXICO',"México",'U.S.A.','U.S.A..',
	'USA - HAWAII','United States','United States of America','usa','USA','Usa','US')
records = subset(records_raw, country %in% mycountries) #331244 records

# Standardize capitalization of genera names
sapply(records$genus, function(x){
	cap.first(gsub(' ','',x))
}) -> records$genus

# Examine families
rec_fams = unique(records$family)
rec_fams[order(rec_fams)]

# Fix misspelled families
records[records$family=="Brigantiaeaceae",'family']<- 'Brigantiaceae'
records[records$family %in% c("Chrysothrichaceae","Chrysotrichaceae"),'family']<- 'Chrysothricaceae'
records[records$family=="Haematommaceae",'family']<- "Haematommataceae" 
records[records$family=="Sphinctrinacaea",'family']<- "Sphinctrinaceae" 
records[records$family=="Stictidacaea",'family']<- "Stictidaceae" 

# Remove records not in targeted families
records = subset(records, family %in% target_fams) # 239781 records

# Examine genera
unique(records$genus)[order(unique(records$genus))]

# Fix misspelled genera
records[records$genus=="Diplotemma",'genus']<- "Diplotomma" 
records[records$genus %in% c("Hypotrachnya","Hypotrachna"),'genus']<-  "Hypotrachyna"
records[records$genus=="Lasillia",'genus']<- "Lasallia" 
records[records$genus=="Tuckermanopsis",'genus']<- 'Tuckermannopsis'
records[records$genus=="Vulpicidia",'genus']<- 'Vulpicida'

# Remove species not in targeted genera - may want to do this for targeted species instead
records = subset(records, genus %in% target_genera) # 158448 records

# Remove records that do not match standardized naming conventions according to CNALH
# May want to try to deal with matching names later
records = subset(records, scientificName %in% cnalh_splist$ScientificName) # 148285 records

# Divide records into a set with coordinates and a set located to county level
records$decimalLatitude = as.numeric(records$decimalLatitude)
records$decimalLongitude = as.numeric(records$decimalLongitude)

records_geo = subset(records, !is.na(decimalLatitude)) #148230
records_county = subset(records, is.na(decimalLatitude)) #55 records from South Carolina

# Records to spatial data
records_sp = records_geo[,c('family','scientificName','genus','specificEpithet','taxonRank',
	'infraspecificEpithet','country','stateProvince','county','decimalLatitude','decimalLongitude',
	'geodeticDatum')]


# Probably don't need to worry about converting datums
datums = unique(records_sp$geodeticDatum)

coordinates(records_sp) = c('decimalLongitude','decimalLatitude')
proj4string(records_sp) = CRS("+proj=longlat")

plot(records_sp)

# Examine and remove records that fall in the water
NA_outline = readOGR('C:/Users/jrcoyle/Documents/GIS shape files/N Am Outline','na_base_Lambert_Azimuthal')

records_sp = spTransform(records_sp, CRS(proj4string(NA_outline)))
rec_poly = overlay(records_sp, NA_outline)

plot(NA_outline)
plot(subset(records_sp, is.na(rec_poly)), pch=1, add=T)

water_recs = select.spatial(records_sp)

records_geo[water_recs,c('decimalLatitude','decimalLongitude','verbatimCoordinates','country','stateProvince','county')]

records

# Fix bad coordinates
records[records$verbatimCoordinates=="44°29'02'N, 9353'15'W",'decimalLongitude'] = -93.8875 
records[records$verbatimCoordinates=="44°41'46'N, 8610'57'W",'decimalLongitude'] = -86.1825
records[records$verbatimCoordinates=="38°27'N, 172°36'02'W",'decimalLongitude'] = -122.6006
records[records$verbatimCoordinates=="7°25'N, 147°10'W",c('decimalLatitude','decimalLongitude')]= c(NA,NA) # South Carolina records- can't figure out what lat lon was supposed to be
records[(records$locality=="Mission Falls area, Flathead Indian Reservation. E of St. Ignatius, MT.")&(records$decimalLongitude<=-120),'decimalLongitude']<- -114
records[(records$locality=="JUNEAU, MENDENHALL LAKE, 5 MI N OF AUKE BAY"),'decimalLatitude'] <- 58.42
records[(records$locality=="Eight miles S of Furlong Bay Provincial Campground along Highway 25, N of Kitimat"),'decimalLatitude'] <- 54.26667
records = records[(records$locality!="Australia; Langi Ghiran State Park, 15km SE of Ararat."),]
records[(records$locality=="Arctic Alaska; East side Jago River, crest of Jag Mt."),'decimalLatitude'] <- 69 # From Google maps
records[(records$locality=="Queen Charlotte Islands, Moresby Island, Port Alliford near Sandspit"),'decimalLongitude'] <- -131.9667 
records[(records$locality=="Island of Maui, Olinda Koolan Forest Reserve, North Haleakala"),'decimalLongitude'] <- -156.2361
records[(records$locality=="NORTH HALEAKALA METROSIDEROS; OLINDA FOREST RESERVE"),'decimalLongitude'] <- -156.2361
records[(records$locality=="NORTH HALEAKALA, OLINDA KOOLAN F; HALEAKALA NP"),'decimalLongitude'] <- -156.2361
records = records[(records$locality!="Mt. Wilhem, Eastern Highlands, Kombugomambuno, Keglsugl-Pindaunde trail"),]
records[records$verbatimCoordinates=="15°07'S, 167°01'E",c('decimalLatitude','decimalLongitude')]<- c(NA,NA) # Carson pass, Alpine CA had wrong coordinates
records[(records$locality=="Island of Kaua'i, Kaaweiki Ridge, near Methodist camp"),'decimalLongitude'] <- -159.6833
records = records[(records$locality!="Hungary Creek Watershed. Raven Lake Mountain Caribou Study. Prince George Region."),] # Bad coordinates and county info
records[(records$locality=="Ocala National Forest, along Co. Rd. 316, ca. 0.8 km E of Oklawaha River bridge at Eureka, 0.32 km W of Forest Service Rd. 67"),'decimalLatitude'] <- 39.47
records[records$verbatimCoordinates=="4113'56'N, 81°35'51'W",'decimalLatitude'] = 41.232222
records[(records$locality=="CHELAN, S OF"),'decimalLatitude'] <- 47.75 
records[(records$locality=="Algoma District:  Lake Superior Provincial Park; cliffs and spruce forest overlooking Old Woman Bay,"),'decimalLatitude'] <- 47.75
records[(records$locality=="Chadron State Park, 19 km S of Chadron off US Route 385, near primative campground"),'decimalLongitude'] <-  -103.0217
records[(records$locality=="0.2 km up dirt road, S side of Hwy 16,  60 km E of Tensleep"),'decimalLongitude'] <- -106.9333
records[(records$locality=="Camp Meriwether Boy Scout Camp"),'decimalLongitude'] <- -124
records[(records$locality=='HIGHLANDS')&(records$verbatimCoordinates=="69°20'N, 145°00'W"),c('decimalLongitude','decimalLatitude')] <- c(-83.2, -83.2, 35,35)

# Used for checking records
bad_recs = records_geo[water_recs,c('decimalLatitude','decimalLongitude','verbatimCoordinates','country','stateProvince','county','locality')]
subset(bad_recs, stateProvince %in% c('Colorado', 'NORTH CAROLINA','Nunavut'))

# Write out record data
write.csv(records, '../FIA Lichen/Data/CNALH_records_fia_genera_NAm.csv', row.names=F)


### Calculate regional species richness for each FIA plot
records = read.csv('../FIA Lichen/Data/CNALH_records_fia_genera_NAm.csv')

library('spatstat') #disc
library('rgeos')

# Define spatial record data
records_geo = subset(records, !is.na(decimalLatitude)) #148224
records_county = subset(records, is.na(decimalLatitude)) #58 records, 55 from South Carolina

# Records to spatial data
records_sp = records_geo[,c('family','scientificName','genus','specificEpithet','taxonRank',
	'infraspecificEpithet','country','stateProvince','county','decimalLatitude','decimalLongitude',
	'geodeticDatum')]
coordinates(records_sp) = c('decimalLongitude','decimalLatitude')
proj4string(records_sp) = CRS("+proj=longlat")

## Make county records into spatial data

# Read in county data for U.S>
county.sh = readOGR('../FIA Lichen/GIS','counties')

records_county$county = cap.first(records_county$county)
records_county$stateProvince = sapply(strsplit(records_county$stateProvince, ' '), function(x){
	if(length(x)==1) return(cap.first(x[1]))
	else return(paste(cap.first(x[1]), cap.first(x[2])))
})

records_poly = records_county[,c('family','scientificName','genus','specificEpithet','taxonRank',
	'infraspecificEpithet','country','stateProvince','county','decimalLatitude','decimalLongitude',
	'geodeticDatum')]

# Assign lat/lon based on the center of the county
coords = data.frame()
for(i in 1:nrow(records_poly)){
	x = records_poly[i,]

	use_poly = county.sh[(county.sh$STATE_NAME==x$stateProvince)&(county.sh$NAME==x$county),]
	these_coords= gCentroid(use_poly)
	coords = rbind(coords, these_coords@coords)
}

records_poly[,c('decimalLongitude','decimalLatitude')] = coords

# Make into spatial data an concatenate with other records
coordinates(records_poly) = c('decimalLongitude','decimalLatitude')
proj4string(records_poly) = CRS("+proj=longlat")

records_sp = rbind(records_sp, records_poly) # 148282 records

# Write out record data
write.csv(records_sp, '../FIA Lichen/Data/CNALH_records_fia_genera_NAm.csv', row.names=F)



#############################################################
### Calculate Regional Richness, All Species and Parmeliaceae/Physciaceae

parmphys = read.csv('Parm_Phys_records_2014-09-20.csv')
allsp = read.csv()

# Subset to genera in FIA data
parmphys_fia = subset(parmphys, genus %in% target_genera)

# Calculate how many species in Physciaceae and Parmeliaceae in records
species = unique(parmphys_fia[,c('family','scientificName','genus')])

species = subset(species, family %in% c('Parmeliaceae','Physciaceae'))
genera = unique(species[,c('family','genus')])
species = unique(species[,c('family','scientificName')])

table(genera$family) # 44 Parmeliaceae, 9 Physciaceae
table(species$family) # 937 Parmeliaceae, 243 Physciaceae

## Calculate regional richness for FIA plots using rarefaction

# Subset to the family of interest
parm = subset(parmphys_fia, family=='Parmeliaceae')
phys = subset(parmphys_fia, family=='Physciaceae')

# Make records into spatial data
parm_sp = parm
phys_sp = phys
coordinates(parm_sp) = c('decimalLongitude','decimalLatitude')
proj4string(parm_sp) = CRS("+proj=longlat")
coordinates(phys_sp) = c('decimalLongitude','decimalLatitude')
proj4string(phys_sp) = CRS("+proj=longlat")

plot(phys_sp)

# Divide FIA data into plots with and without coordinates
fia_geo = subset(master_locs, !is.na(LAT))

# Calculate for FIA plots with spatial coordinates
coordinates(fia_geo) = c('LON','LAT')
proj4string(fia_geo) = CRS("+proj=longlat")

# Determine which records will be used to calculate regional richness
recs = parm_sp

# Calculate the number of records within 500 km buffer (Could easily check other buffers)
dist500 = c()
for(i in 1:nrow(fia_geo)){
	dist500 = c(dist500, nrow(find_recs(fia_geo[i,], recs, 500)))
} # min=499 all sp, min=144 Parmeliaceae, min=8 Physciaceae, but used 144 since only 26 plots with fewer than 144.
names(dist500) = fia_geo$yrplot.id

min(dist500) # Min num records for all FIA plots 
# Parm = 598 Phys = 60
min(dist500[fia_geo$yrplot.id %in% rownames(model_data)]) # Min num records for plots used in analysis
# Parm = 886 Phys = 374


## Calculate Regional Richness using rarefaction (vegan) only for plots used in analysis
regS = data.frame(yrplot.id=rownames(model_data))

# Define size of subsample (nsamps)
nsamp = 350

recs = parm_sp

richness = c()

for(i in 1:nrow(model_data)){
	usepoint = fia_geo[fia_geo$yrplot.id==rownames(model_data)[i],]
	userecs = find_recs(usepoint, recs, 500)
	comm = tab_species(userecs$scientificName)
	rich = rarefy(comm, nsamp)
	richness = c(richness, as.numeric(rich))
}


# Compare richness with sampling effort (# records)
png('../FIA Lichen/Figures/ regS Parm (500 record rarefaction) vs sampling effort.png', height=600, width=600)
par(mar=c(5,6,1,1.5))
plot(richness~dist500[rownames(model_data)], las=1, xlab='Number of records w/in 500km radius', ylab='Rarefied richness for 350 records', cex=2, cex.lab=1.5, cex.axis=1.5)
dev.off()


# Save data
regS$regPhys = richness
regS$regParm = richness
write.csv(regS, '../FIA Lichen/Data/Regional Richness/ParmPhys_regS_CNALH-2014-09-20.csv', row.names=F)


#

# Quick Plot 
spplot(fia_geo, 'regS', col.regions=colorRampPalette(c('dark blue','blue','green','yellow','orange','red'))(10), cuts=10)

# Still biased toward higher richness in more frequently sampled areas
png('../FIA Lichen/Figures/regional richness (500 records subsample) vs sampling effort.png', height=600, width=600)
par(mar=c(5,6,1,1))
plot(regS~dist500, las=1, xlab='Number of records w/in 500km radius', ylab='Avg. richness of 499 records', cex=2, cex.lab=1.5, cex.axis=1.5)
dev.off()



## Assess how sampling radius (km) affects richness estimates
#Note: subsampled different number of records in each case b/c more records in larger radii

rich1000 = read.csv('./Data/Regional Richness/fia_lichen_reg_richness_1000km.csv')
rich700 = read.csv('./Data/Regional Richness/fia_lichen_reg_richness_700km.csv')
rich600 = read.csv('./Data/Regional Richness/fia_lichen_reg_richness_600km.csv')
rich500 = read.csv('./Data/Regional Richness/fia_lichen_reg_richness.csv')
rich400 = read.csv('./Data/Regional Richness/fia_lichen_reg_richness_400km.csv')
rich200 = read.csv('./Data/Regional Richness/fia_lichen_reg_richness_200km.csv')


plot(rich1000$regS~rich500$regS)

ranks = data.frame(R1000 = rank(rich1000$regS),
	R700 = rank(rich700$regS),
	R600 = rank(rich600$regS),
	R500 = rank(rich500$regS),
	R400 = rank(rich400$regS),
	R200 = rank(rich200$regS)
)

rich = data.frame(R1000 = rich1000$regS,
	R700 = rich700$regS,
	R600 = rich600$regS,
	R500 = rich500$regS,
	R400 = rich400$regS, 
	R200 = rich200$regS
)

plot(ranks$R500,ranks$R400)
plot(ranks$R500,ranks$R600)
plot(ranks$R500,ranks$R700)

cor(ranks)











## Calculate for FIA plots without spatial coordinates

# Assign plot location to centroid of counties
coords = data.frame()
for(i in 1:nrow(fia_county)){
	x = fia_county[i,]

	use_poly = county.sh[(county.sh$STATE_NAME==x$STATE_NAME)&(county.sh$NAME==x$COUNTYNM),]
	these_coords= gCentroid(use_poly)
	coords = rbind(coords, these_coords@coords)
}

# Calculate richness for each county 
county_locs = unique(cbind(fia_county[,c('STATE_NAME','COUNTYNM')], coords))
coordinates(county_locs) = c('x','y')
proj4string(county_locs) = CRS("+proj=longlat")

# Remove counties that have already been calculated
prev_county = read.csv('../FIA Lichen/Data/regional richness fia county level.csv')
prev_county[,c('STATE_NAME','COUNTYNM')]
merge(county_locs, prev_county) # Looks like all counties have already bee calculated

# Calculate the number of records within various buffers
dist500 = c()
for(i in 1:nrow(county_locs)){
	dist500 = c(dist500, nrow(find_recs(county_locs[i,], records_sp, 500)))
} # min=1131, min=530 for Parm, min=221 for Phys

nsamp = 499 # minimum from fia plots with coordinates 
reps = 500

regS = c()
for(i in 1:nrow(county_locs)){
	userecs = find_recs(county_locs[i,], records_sp, 500)
	regS = c(regS, calc_rich_samp(userecs, nsamp, reps))
}
county_locs$regS = regS
county_locs$regParm = regS #STOPPed HERE
county_locs$regPhys = regS

# Save county-level regional richness
write.csv(county_locs, '../FIA Lichen/Data/regional Physciaceae richness fia county level.csv', row.names=F)
county_loc = read.csv('../FIA Lichen/Data/regional Physciaceae richness fia county level.csv')

fia_county = fia_county[,-which(names(fia_county)=='regS')]

fia_county = merge(fia_county, county_loc, all.x=T)

# Combine county and geo data and put into master dataframe
regS_all = rbind(fia_county[,c('yrplot.id','regS')],fia_geo@data[,c('yrplot.id','regS')])

# Merge data when only a few regS have been calculated for plots with missing values
master[regS_all$yrplot.id, 'regS'] = regS_all$regS

# Merge data when all regS are calculated at once across all plots
rownames(master) = master$yrplot.id)
master = merge(master, regS_all)

## Merge data when only Parmeliaceae or Physciaceae regional richness have been calculated

# merge county-level data with plot ids
smdata = merge(fia_county[,c('yrplot.id','STATE_NAME','COUNTYNM')], county_loc[,c('STATE_NAME','COUNTYNM','regPhys')], all.x=T)

# create small data frame of just Parm or Phys regional richness to append to main data frame
newdata = rbind(smdata[,c('yrplot.id','regPhys')], fia_geo@data[,c('yrplot.id','regPhys')])

# This was done on cluster, while Parm was done on laptop
write.csv(newdata, 'phys reg rich.csv', row.names=F)

# Then I merged Phys and Parm data together
phys_reg = read.csv('./Data/phys reg rich.csv')
newdata = merge(newdata, phys_reg, all.x=T)

# The I merged both with the master data
master = merge(master, newdata, all.x=T)
model_data = merge(model_data, newdata, all.x=T)

write.csv(master,'../FIA Lichen/Data/master_data.csv', row.names=F)
write.csv(master[,c('yrplot.id','regS','regParm','regPhys')], '../FIA Lichen/Data/regional_lichen_rich.csv', row.names=F)

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

write.csv(bindata, 'Phys_richness_2-deg_CNALH-2014-09-20.csv', row.names=F)
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


# plot

colorvar = cut(mygrid_data$regParm300, breaks=seq(40, 130, length.out=9))

plot(nam_poly1)
plot(mygrid_sp, col=mycolors[colorvar], add=T, pch=15, cex=.4)


windows()
par(mar=c(0,0,0,0))
plot(0,0, xlim=c(0,1), ylim=c(-.2,1.2), type='n')
rect(rep(0,8), seq(0, .7, .1), rep(.1, 8), seq(.1, .8, .1), col=mycolors)
text(rep(0.1, 9), seq(0,.9,.1), labels=seq(40, 140, length.out=9), pos=4, cex=2)

par('usr')



### CODE BELOW NOT YET CORRECT FOR THIS SRCIPT






# Calc  regS for each grid cell 

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
	splitnames = strsplit(x, " ")
	justgenus = unique(x[as.numeric(lapply(splitnames, length))==1])
	allspecies = splitnames[as.numeric(lapply(splitnames, length))>1]
	allspecies = unique(sapply(allspecies, function(s) paste(s[1],s[2])))
	spgenera = unique(sapply(strsplit(allspecies, " "), function(y) y[1]))
	keepgenera = justgenus[!(justgenus %in% spgenera)]

	c(allspecies,keepgenera)
}

# A function that returns a vector counting the number of records for each species binomial
tab_species = function(x){
	if(length(x)>0){

	splitnames = strsplit(x, " ")
	justgenus = x[as.numeric(lapply(splitnames, length))==1]
	allspecies = splitnames[as.numeric(lapply(splitnames, length))>1]
	allspecies = sapply(allspecies, function(s) paste(s[1],s[2]))
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









