### This script extract environmental data for FIA phase 3 plots used in my lichen species richness analysis
library(sp)
library(rgdal)
library(raster)

##########################################################################
### Variables ###

#mydir = 'C:/Users/jrcoyle/Documents/UNC/Projects/FIA Lichen'
mydir = 'C:/Users/jrcoyle/Documents/Projects/FIA Lichen'

my.bio = c(1,3,12,15) # Use mean annual temp, isothermality, annual precip, and precip seasonality from WorldClim

##########################################################################
### Data ###

setwd(mydir)

# Plot locations
master_loc = read.csv('./Data/fia_lichen_plot_locations.csv')

# Lichen plot tree data
master_forest = read.csv('./Data/master_data_forest.csv')

# Previously calculated data
master = read.csv('./Data/master_data.csv')
#master = merge(master, master_forest)

# WorldClim  data
#load('Z:/GIS/WorldClim data/bioclim.Rdata')
#load('//BioArk/HurlbertLab/GIS/WorldClim data/bioclim.Rdata')
load('../../DBDGS/bioclim.Rdata')

# PET - probably won't use
#setwd('C:/Users/jrcoyle/Documents/UNC/Projects/GIS Data/Global PET - Annual/PET_he_annual')
#pet = raster('pet_he_yr')


# PRISM Modeled humidity data
setwd('C:/Users/jrcoyle/Documents/Projects/GIS Data/PRISM')
#setwd('C:/Users/jrcoyle/Documents/UNC/Projects/GIS Data/PRISM')
vp = raster('mean_vapor_pressure_Jan-Dec_1997-2008')
rh = raster('mean_relative_humidity_Jan-Dec_1997-2008.grd')
hum = stack(vp, rh)
names(hum) = c('vp','rh')

# GSOD Alaska humidity weather station data
setwd('C:/Users/jrcoyle/Documents/Projects/GIS Data/GSOD/GSOD_SE_AK_1997-2008')
#setwd('C:/Users/jrcoyle/Documents/UNC/Projects/GIS Data/GSOD/GSOD_SE_AK_1997-2008')
humidityAK = read.csv('avg_humidity_SE_AK_1997-2008.csv')
humidityAK = subset(humidityAK, n_years>0)

# Modeled pollution data
nitsul = raster('C:/Users/jrcoyle/Documents/Projects/GIS Data/Air Pollution/NTN/ntn_NplusS_1997-2008_avg.grd')
ntn = read.csv('C:/Users/jrcoyle/Documents/Projects/GIS Data/Air Pollution/NTN/ntn_site_data_1997-2008.csv', skip=3)

# Shapefiles of US counties
setwd(mydir)
county.sh = readOGR('./GIS','counties')

# FIA county name data
setwd(mydir)
county.names = read.csv('./Data/FIA_county_data_match_shapefiles.csv') # File has been altered from what was originally downloaded from FIA to match county shapefiles.

# FIA county codes and state codes
#setwd('C:/Users/jrcoyle/Documents/UNC/Projects/FIA_Fall2011/Data')
setwd('C:/Users/jrcoyle/Documents/Projects/FIA_Fall2011/Data')
fia.stcd = read.csv('FIA_state_codes.csv')


#########################################################################
### Code ###

### Focus data layers only to extent of US data
bio.stack = crop(bio.stack, extent(county.sh))
#pet = crop(pet, extent(county.sh))

# Use only WorldClim variables of interest
bio.stack = stack(bio.stack)
bio.stack = subset(bio.stack, my.bio)


### Divide dataset into plots with points and plots in counties
## There are now only 31 plots without coordinates so we will not be using the calculations for county level data.
master.point = subset(master_loc, !is.na(LAT))
master.county = subset(master_loc, is.na(LAT))

### Extract values for plots with lat/lon

# Convert dataframe to spatial data
coordinates(master.point) = c('LON','LAT')
projection(master.point) = CRS("+proj=longlat +datum=NAD83")

## Extract from WorldClim
wc.pts = extract(bio.stack, master.point)

# Assign NA values to value of nearest cell - Only need to do this for mat (bio.stack[[1]]) #9 plots in ID(5), WY(3), MT(1)
nas = which(is.na(wc.pts))
napts = master.point[nas,]
napts = spTransform(napts, CRS(projection(bio.stack)))

for(i in 1:nrow(napts)){
	useval = findNearestVals(bio.stack[[1]], napts[i,])
	wc.pts[nas[i],1]<-useval
}

## Extract from PET
#pet.pts = extract(pet, master.point)

## Extract from humidity data
hum.pts = extract(hum, master.point)

# Assign NA values to value of nearest # 7 plots in ME(4), NC(1), WA(1), and OR(1)
nas = unique(which(is.na(hum.pts), arr.ind=T)[,'row'])
napts = master.point[nas,] # All AK plots- will be dealt with later.

## Extract from pollution data
nitsul.pts = extract(nitsul, master.point)

## Assign AK points a value based on Juneau NTN site 'AK02'
ntnAK = subset(ntn, SiteID=='AK02')

# Convert NH4, NO3, SO4 to total N+S
AK_ns = ((55.44*ntnAK$NH4) + (16.13*ntnAK$NO3) + (20.83*ntnAK$SO4))*(ntnAK$Ppt.)*.001

nitsul.pts[which(master.point$state.abbr=='AK')]<-mean(AK_ns)

# Assign NA values to value of nearest # 7 plots in ME(4), NC(1), WA(1), and OR(1)
nas = which(is.na(nitsul.pts))
napts = master.point[nas,]
napts = spTransform(napts, CRS(projection(nitsul)))

for(i in 1:nrow(napts)){
	useval = findNearestVals(nitsul, napts[i,])
	nitsul.pts[nas[i]]<-useval
}

# Combine environmental columns into a dataframe
env.pts = data.frame(yrplot.id = master.point$yrplot.id, wc.pts, hum.pts, totalNS = nitsul.pts) # temp seems low to me!

# For AK data, match with closest weather station and use that data
matchplotids = subset(env.pts, is.na(env.pts$vp))$yrplot.id
matchplots = subset(master.point, yrplot.id %in% matchplotids)

coordinates(humidityAK) = c('lon','lat')
projection(humidityAK) = CRS("+proj=longlat +datum=NAD83")

dists = spDists(matchplots, humidityAK, longlat=T)
apply(dists, 1, function(x){
	c(min(x),which(x==min(x)))
})->hummatch

humdata = data.frame(matchplots$yrplot.id, t(hummatch))
names(humdata) = c('yrplot.id','match_dist','match')
humdata = cbind(humdata, humidityAK[humdata$match,c('vp','rh','n_years','name','elev','station')])
rownames(humdata) = as.character(humdata$yrplot.id)

# This file pairs shows which station each plot was paired with including the elevation of the plot (ELEVATION) and station (elev)
humdata = merge(humdata, master[,c('yrplot.id','ELEVATION')], all.x=T)
write.csv(humdata, 'AK_matched_humidity_data.csv', row.names=F)

env.pts[is.na(env.pts$vp),c('vp','rh')]<-humdata[as.character(matchplotids),c('vp','rh')]

setwd(mydir)
write.csv(env.pts, './Data/fia_lichen_env_data_points.csv', row.names=F)


#################################################################
### Extract values for plots only known to county

## This code is no longer neaded

## Link FIA county names and shape file county names

# Add state names to county.names
county.names = merge(county.names, fia.stcd)

# Only worry about counties in the lichen data set
lichen.counties = merge(unique(master[,c('STATE','COUNTY')]), county.names, by.x=c('STATE','COUNTY'), by.y=c('STATECD','COUNTYCD'), all.x=T)

# Check for mis-matches by state: Re-run after adding in additional Alaska counties
no.match = data.frame()
county.df = data.frame()
for(s in unique(lichen.counties$state.name)){
	fia.ctys = subset(lichen.counties, state.name==s)
	sh.ctys = subset(county.sh@data, STATE_NAME==s)

	matches = merge(fia.ctys, sh.ctys, all.x=T, by.x='COUNTYNM',by.y='NAME')
	
	no.match = rbind(no.match,subset(matches, is.na(STATE_NAME)))
	county.df = rbind(county.df, matches)
}

# Write out county data once new AK counties have been added
write.csv(county.df, 'fia_lichen_county_data.csv', row.names=F)


# There are only sixteen- matched by hand in new fia county file: 'FIA_county_data_match_shapefiles.csv'
# La Salle -> LaSalle in Illinois
# La Porte -> LaPorte in Indiana
# King And Quees -> King and Queen in Virginia
# Bedford -> Bedford City in Virginia
# Suffolk city -> Suffolk in Virginia
# Okanagon -> Okanogan in Washington
# Haines Borough -> Haines in Alaska
# Juneau Borough -> Juneau in Alaska
# Kenai Peninsula Borough -> Kenai Peninsula in Alaska
# Ketchikan Gateway Borough -> Ketchikan Gateway in Alaska
# Prince of Wales-Outer Ketchikan Census Area -> Prince of Wales-Hyder in Alaska
# Sitka Borough -> Sitka in Alaska
# Valdez-Cordova Census Area -> Valdez-Cordova in Alaska
# Yakutat Borough -> Yakutat in Alaska

# Add two combined counties - # CANT GET THIS WORKING
# Skagway-Hoonah-Angoon Census Area  = Hoonah-Angoon + Skagway in Alaska

HA = county.sh[county.sh$NAME=='Hoonah-Angoon',]
SKWY = county.sh[county.sh$NAME=='Skagway',]
HASK = SpatialPolygons(list(Polygons(c(SKWY@polygons[[1]]@Polygons,HA@polygons[[1]]@Polygons), "Skagway-Hoonah-Angoon Census Area")), proj4string=CRS(proj4string(county.sh)))

HASKdata = HA@data
HASKdata$NAME="Skagway-Hoonah-Angoon Census Area"
HASKdata[,c('POP2000','POP2010','SQMI')] = HA@data[,c('POP2000','POP2010','SQMI')]+SKWY@data[,c('POP2000','POP2010','SQMI')]
HASKdata[,c('POP00_SQMI','POP10_SQMI')] = HASKdata[,c('POP2000','POP2010')]/HASKdata[,'SQMI']
rownames(HASKdata) = HASKdata$NAME
HASKdf = SpatialPolygonsDataFrame(HASK, HASKdata)

county.sh = rbind(county.sh, HASKdf)

# Wrangell-Petersburg Census Area = Wrangell + Petersburg in Alaska

WR = county.sh[county.sh$NAME=='Wrangell',]
PB = county.sh[(county.sh$NAME=='Petersburg')&(county.sh$STATE_NAME=='Alaska'),]
WRPB= SpatialPolygons(list(Polygons(c(WR@polygons[[1]]@Polygons,PB@polygons[[1]]@Polygons), "Wrangell-Petersburg Census Area")), proj4string=CRS(proj4string(county.sh)))

WRPBdata = WR@data
WRPBdata$NAME="Wrangell-Petersburg Census Area"
WRPBdata[,c('POP2000','POP2010','SQMI')] = WR@data[,c('POP2000','POP2010','SQMI')]+PB@data[,c('POP2000','POP2010','SQMI')]
WRPBdata[,c('POP00_SQMI','POP10_SQMI')] = WRPBdata[,c('POP2000','POP2010')]/WRPBdata[,'SQMI']
rownames(WRPBdata) = WRPBdata$NAME
WRPBdf = SpatialPolygonsDataFrame(WRPB, WRPBdata)

county.sh = rbind(county.sh, WRPBdf)

# Re-format county names to characters rather than factors
county.df$COUNTYNM = as.character(county.df$COUNTYNM)
county.df$STATE_NAME = as.character(county.df$STATE_NAME)
county.df$state.name = as.character(county.df$state.name)


## Extract environmental variables as means across all counties

# Read in fixed county data
county.df = read.csv('./Data/fia_lichen_county_data.csv', stringsAsFactors=F)

county.envmeans = data.frame()

for(i in 1:nrow(county.df)){
	x = county.df[i,c('STATE_NAME','COUNTYNM')]
	x = as.character(x)

	# Find county in shapefile
	use.poly = county.sh[(county.sh$STATE_NAME==x[1])&(county.sh$NAME==x[2]),]

	# Extract data
	wclim = extract(bio.stack, use.poly, fun=mean, na.rm=T, small=T)
	if(x[1]=='Alaska'){
		this_hum = data.frame(vp=NA, rh=NA)
		totalNS = NA
	} else {
		this_hum = extract(hum, use.poly, fun=mean, na.rm=T, small=T)
		totalNS = extract(nitsul, use.poly, fun=mean, na.rm=T, small=T)
	}
	
	county.envmeans = rbind(county.envmeans, data.frame(wclim, this_hum, totalNS))
	
}

county.df = cbind(county.df[,c('COUNTYNM','STATE_NAME','STATE','COUNTY','state.abbr','SQMI','POP10_SQMI')], county.envmeans)

# Re-name env vars to indicate that they are county-level means
names(county.df)[8:ncol(county.df)] = paste(names(county.df)[8:ncol(county.df)], 'county', sep='_')

write.csv(county.df, './Data/fia_lichen_env_data_county.csv', row.names=F)

## Merge point and county level data with plots data

env.pts = read.csv('./Data/fia_lichen_env_data_points.csv')
county.df = read.csv('./Data/fia_lichen_env_data_county.csv')
regS = read.csv('./Data/regional_lichen_rich.csv')

# Merge data
master = merge(master_forest, county.df, all.x=T)
master = merge(master, env.pts, all.x=T)
master = merge(master, regS, all.x=T)

mycols = c('darkblue','blue','lightblue','yellow','orange','red')

## Compare point vs. county-level mean environment
png('./Figures/Spring 2013 Analyses/Env vars county vs point estimates color by Lat.png', height=600, width=1000)
par(mfrow=c(2,2))
plot(mat~mat_county, data=master, las=1, xlab='County-level average', 
	ylab='At plot coordinates', main='Mean Annual Temperature (C)', 
	col=colorRampPalette(mycols)(20)[cut(subset(master, !is.na(mat))[,'LAT'], breaks=20, include.lowest=T)], cex=2)
abline(0,1, lwd=2)
plot(ap~ap_county, data=master, las=1, xlab='County-level average', 
	ylab='At plot coordinates', main='Annual Precipitation (mm)',
	col=colorRampPalette(mycols)(20)[cut(subset(master, !is.na(ap))[,'LAT'], breaks=20, include.lowest=T)], cex=2)
abline(0,1, lwd=2)
plot(rh~rh_county, data=master, las=1, xlab='County-level average', 
	ylab='At plot coordinates', main='Relative Humidity',
	col=colorRampPalette(mycols)(20)[cut(subset(master, !is.na(ap))[,'LAT'], breaks=20, include.lowest=T)], cex=2)
abline(0,1, lwd=2)
plot(vp~vp_county, data=master, las=1, xlab='County-level average', 
	ylab='At plot coordinates', main='Vapor Pressure (mbar)',
	col=colorRampPalette(mycols)(20)[cut(subset(master, !is.na(ap))[,'LAT'], breaks=20, include.lowest=T)], cex=2)
abline(0,1, lwd=2)

dev.off()


# Write out master data

setwd(mydir)
write.csv(master, './Data/master_data.csv', row.names=F)



#############################################################
### Maps
mymapcolors = c(rgb(69,117,180,maxColorValue=255),rgb(145,191,219,maxColorValue=255),
rgb(224,243,248,maxColorValue=255),rgb(254,224,144,maxColorValue=255),rgb(252,141,89,maxColorValue=255),
rgb(215,48,39,maxColorValue=255))




## Relative humidity map

png('Avg annual relative humidity 1997-2008.png', height=1000, width=1000)
spplot(trim(rh), axes=F, colorkey=list(labels=list(cex=3), title='Avg. annual relative humidity'),
	col.regions=colorRampPalette(mymapcolors[length(mymapcolors):1])(20), border=F)
dev.off()






#####################################################################
### Functions

# A function that finds the closest non-NA raster cells to a point and calculates their average values
# Usually this is just on cell
# May not be the best way to do this

findNearestVals = function(rast, point){
	ptdist = distanceFromPoints(rast, point)
	ptdist = mask(ptdist, rast)
	ptdist_na = ptdist==minValue(ptdist)
	ptdist_na[ptdist_na==0]<-NA
	useval = mask(rast, ptdist_na)
	useval = cellStats(useval, 'mean', na.rm=T)
}



