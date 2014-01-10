# This script diversity of lichens from occurence records downloaded from
# Consortium of North American Lichen Herbaria on 2-19-2012.

library(sp)
library(rgdal)

#setwd('C:/Users/jrcoyle/Documents/UNC/Projects/CNALH Diversity/')
setwd('C:/Users/jrcoyle/Documents/Projects/CNALH Diversity/')
options(stringsAsFactors=F)

# Read in CNALH records
records_raw = read.csv('CNALH_allrecords_NAm_2013-02-19.csv') #355119 records

# Read in master data of FIA plots with plot data and environmental data
master = read.csv('../FIA Lichen/Data/master_data.csv')
master_locs = read.csv('../FIA Lichen/Data/fia_lichen_plot_locations.csv')

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

### Calculate how many species in Physciaceae and Parmeliaceae in records
species = unique(records_sp@data[,c('family','scientificName','genus')])
species = species[order(species$scientificName),]

species = subset(species, family %in% c('Parmeliaceae','Physciaceae'))
genera = unique(species[,c('family','genus')])
species = unique(species[,c('family','scientificName')])

table(genera$family) # 45 Parmeliaceae, 9 Physciaceae
table(species$family) # 698 Parmeliaceae, 193 Physciaceae

### Calculate regional richness

# Read records back in as a dataframe
records_sp = read.csv('../FIA Lichen/Data/CNALH_records_fia_genera_NAm.csv')


# When calculating just richness within a family, do this subsetting
records_sp = subset(records_sp, family=='Parmeliaceae')

# Make records into spatial data
coordinates(records_sp) = c('decimalLongitude','decimalLatitude')
proj4string(records_sp) = CRS("+proj=longlat")

# Divide FIA data into plots with and without coordinates
# I am no longer using the code for county-level data
fia_geo = subset(master_locs, !is.na(LAT))
fia_county = subset(master_locs, is.na(LAT))

## Calculate for FIA plots with spatial coordinates
coordinates(fia_geo) = c('LON','LAT')
proj4string(fia_geo) = CRS("+proj=longlat")

# Calculate the number of records within various buffers
dist500 = c()
for(i in 1:nrow(fia_geo)){
	dist500 = c(dist500, nrow(find_recs(fia_geo[i,], records_sp, 500)))
} # min=499 all sp, min=144 Parmeliaceae, min=8 Physciaceae, but used 144 since only 26 plots with fewer than 144.

dist1000 = c()
for(i in 1:nrow(fia_geo)){
	dist1000 = c(dist1000, nrow(find_recs(fia_geo[i,], records_sp, 1000)))
} # min=1865

dist700 = c()
for(i in 1:nrow(fia_geo)){
	dist700 = c(dist700, nrow(find_recs(fia_geo[i,], records_sp, 700)))
} # min=772

dist600 = c()
for(i in 1:nrow(fia_geo)){
	dist600 = c(dist600, nrow(find_recs(fia_geo[i,], records_sp, 600)))
} # min=655

dist400 = c()
for(i in 1:nrow(fia_geo)){
	dist400 = c(dist400, nrow(find_recs(fia_geo[i,], records_sp, 400)))
} # min=422

dist200 = c()
for(i in 1:nrow(fia_geo)){
	dist200 = c(dist200, nrow(find_recs(fia_geo[i,], records_sp, 200)))
} # min=12, but used 100

dist100 = c()
for(i in 1:nrow(fia_geo)){
	dist100 = c(dist100, nrow(find_recs(fia_geo[i,], records_sp, 100)))
} # 12 plots have no nearby records

dist50 = c()
for(i in 1:nrow(fia_geo)){
	dist50 = c(dist50, nrow(find_recs(fia_geo[i,], records_sp, 50)))
} # Many plots have no nearby records.

# Only calculate for plots that are missing regional richness data
# This was done after the initial calculation
#fia_geo = subset(fia_geo, is.na(regS))
#fia_county = subset(fia_county, is.na(regS))

# Define size of subsample (nsamps) and number of times to resample (reps)
nsamp = min(dist500)
reps = 500

regS = c()

for(i in 1:nrow(fia_geo)){
	userecs = find_recs(fia_geo[i,], records_sp, 500)
	regS = c(regS, calc_rich_samp(userecs, nsamp, reps))
}

fia_geo$regS = regS
fia_geo$regParm = regS
fia_geo$regPhys = regS

## Plot 
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
bindata = read.csv('Richness all FIA genera 5 degree bins.csv')

# records does not contain 58 records w/o coordinates.  May want to change this.
records = read.csv('../FIA Lichen/Data/CNALH_records_fia_genera_NAm.csv')

# Make Lat/lon bins
binwidth = 5
mylons = seq(-167, -50, binwidth)
mylats = seq(15, 80, binwidth)

records$lonbin = cut(records$decimalLongitude, breaks=mylons)
records$latbin = cut(records$decimalLatitude, breaks=mylats)


# Each cell is a list of the species occuring there
species.distribution = tapply(records[,c('scientificName')], list(records$latbin, records$lonbin), FUN=get_species)

# Calculate richness
richness.raw = apply(species.distribution,c(1,2), function(x) length(unlist(x)))
class(richness.raw) = "table"

bindata = as.data.frame(richness.raw)
names(bindata) = c('latbin','lonbin','richness.raw')
bindata$lat = rep(mylats[1:(length(mylats)-1)],length(mylons)-1)
bindata$lon = rep(mylons[1:(length(mylons)-1)],each=length(mylats)-1)


# Plot on map

#mycolors = colorRampPalette(c("white","yellow","red", "darkred"))(7)

mycolors = c('#FFFFFF',colorRampPalette(c('#2c7bb6','#abd9e9','#fdae61','#ca373b'))(7))
mycolorsT = paste(mycolors, "80", sep='')

colorvar = cut(bindata$richness.raw, breaks=c(0,seq(1,max(bindata$richness.raw), length.out=8)), right=F, include.lowest=T)

tiff('fia_genera_species_richness_uncorrected_map.tif', width=800, height=800)
map('world', c('canada','usa','mexico'),col='grey50', xlim=c(-180,-45))
points(decimalLatitude~decimalLongitude, data=records, pch='.', cex=3)
rect(bindata$lon,bindata$lat,bindata$lon+binwidth,bindata$lat+binwidth,
	col=mycolorsT[colorvar], border='transparent')

legend('bottomleft',levels(colorvar)[2:8], pch=15, col=mycolorsT[2:8], bty='n', cex=2)
title(main='FIA Genera Species Richness\n(Uncorrected)', cex.main=2)
dev.off()


# Calculate richness after subsampling
ITERATIONS=100

tapply(records[,c('scientificName')], list(records$latbin, records$lonbin), FUN=function(x){
	if(length(x)>=10){
	r.dist = c()
	for(i in 1:ITERATIONS){
		use.x = sample(x, 10, replace=F) 
		r.dist = c(r.dist, length(get_species(use.x)))
	}
	r.dist
	} else {
	r.dist=rep(length(get_species(x)),ITERATIONS)
	}
})->subsample

# Calculate the median of the samples
richness.samp = apply(subsample, c(1,2), function(x) ifelse(is.null(unlist(x)),0,median(unlist(x))))
class(richness.samp) = 'table'
richness.samp = as.data.frame(richness.samp)
names(richness.samp) = c('latbin','lonbin','richness.samp10')
bindata = merge(bindata, richness.samp)

# Tally number of samples in each bin
samples = xtabs(~latbin+lonbin, data=records)
samples = as.data.frame(samples)
names(samples) = c('latbin','lonbin','n.obs')

bindata = merge(bindata, samples)

which(bindata$n.obs==max(bindata$n.obs))
bindata[order(bindata$n.obs, decreasing=T),]
plot(richness.samp200~log10(n.obs), data=bindata)

#write.csv(bindata, 'Richness all FIA genera 5 degree bins.csv')
#write.csv(bindata, 'Richness all FIA genera 2 degree bins.csv')

# Map subsampled richness
tiff('fia_genera_species_richness_samp200_map.tif', width=800, height=800)
#png('parmeliaceae_species_richness_samp200_map.png', height=1200, width=1200, bg="transparent")

bordervar = bindata$n.obs>=200
colorvar = cut(bindata$richness.samp200, breaks=c(0,seq(min(bindata$richness.samp200[bordervar]),max(bindata$richness.samp200[bordervar]), length.out=8)), right=F, include.lowest=T)
colorvar[!bordervar]<-NA

map('world',c('canada','usa','mexico'), col='grey50', xlim=c(-180,-50))
points(decimalLatitude~decimalLongitude, data=records, pch='.', cex=3)
rect(bindata$lon,bindata$lat,bindata$lon+binwidth,bindata$lat+binwidth,
	col=mycolorsT[colorvar], border="transparent")
legend('bottomleft',levels(colorvar)[2:8], pch=15, col=mycolorsT[2:8], bty='n', cex=2)
title(main='FIA Genera Species Richness\n(subsample 200 observations)', cex.main=2)
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
######################################################################
### Calculate a few collectors curves

subset(bindata, richness.raw>100)

i=147

i=160

REPS=100

this.bin = subset(parmelia, (latbin==bindata$latbin[i])&(lonbin==bindata$lonbin[i]))
nrow(this.bin)

sapply(seq(10, nrow(this.bin),50), function(x){
		
	sapply(1:REPS, function(t) length(get.species(sample(this.bin$scientificName, x))))	

})->subsamples


png(paste('parmeliaceae',bindata$lat[i],bindata$lon[i],'.png',sep='_'), 
	height=500, width=500, bg='transparent') 
par(mar=c(5,5.5,1,1))
plot(rep(seq(10, nrow(this.bin),50), each=REPS), subsamples, las=1, xlab='Number of Samples',
	ylab='', cex.lab=2,cex.axis=1.8)
title(ylab='Number of Species', line=4, cex.lab=2)
text(nrow(this.bin)/2,10,paste('Lat = ',bindata$latbin[i],'  Lon = ',bindata$lonbin[i],sep=''), cex=1.5)

dev.off()




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









