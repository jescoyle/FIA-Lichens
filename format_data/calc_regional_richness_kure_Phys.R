
library(sp)
library(rgdal)
library('spatstat') #disc
library('rgeos')
options(stringsAsFactors=F)

### Functions ###

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

### CNALH Records ###
# Read in CNALH records
records_sp = read.csv('../CNALH_records_fia_genera_NAm.csv')

# Subset by family
records_sp = subset(records_sp, family=='Physciaceae')

# Make records into spatial data
coordinates(records_sp) = c('decimalLongitude','decimalLatitude')
proj4string(records_sp) = CRS("+proj=longlat")


### FIA Data ###
master_locs = read.csv('../fia_lichen_plot_locations.csv')
fia_geo = subset(master_locs, !is.na(LAT))
coordinates(fia_geo) = c('LON','LAT')
proj4string(fia_geo) = CRS("+proj=longlat")

## Calculate richness by subsampling
nsamp = 144
reps = 500

regS = c()
for(i in 1:nrow(fia_geo)){
	userecs = find_recs(fia_geo[i,], records_sp, 500)
	regS = c(regS, calc_rich_samp(userecs, nsamp, reps))
	print(i)
}

save.image('calc_regPhys.Rdata')
write.csv(data.frame(yrplot.id=fia_geo$yrplot.id, regPhys=regS), 'fia_lichen_reg_Phys_richness.csv', row.names=F)