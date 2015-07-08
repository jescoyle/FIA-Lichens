## This script calculates regional species richness across several different radius buffers
## It uses records downloaded from CNALH in the Fall of 2014 which have previously been cleaned up to only contain genera in the FIA data set.
## See calc_regional_richness.R script.

library(sp)
library(rgdal)
library(vegan) # rarefy
library(stringr) # str_trim

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

##########################################################################
### Data & Analysis

# Read in CNALH records
allsp = read.csv('CNALH_records_fia_genera_NAm_2014-10-08.csv')

# Read in plot locations
master_locs = read.csv('../fia_lichen_plot_locations.csv')

# Read in data frame used for modeling
model_data = read.csv('fia_lichen_model_data.csv', row.names=1)

# Make records into spatial data
allsp_sp = allsp
coordinates(allsp_sp) = c('decimalLongitude','decimalLatitude')
proj4string(allsp_sp) = CRS("+proj=longlat")

# Only calculate regional richness for plots used in models
# This will exclude AK plots where I did not download records for
fia_geo = subset(master_locs, yrplot.id %in% rownames(model_data))
coordinates(fia_geo) = c('LON','LAT')
proj4string(fia_geo) = CRS("+proj=longlat")

# Define scales over which richness should be calculated
scales = c(50, 100, 250, 400, 500)
nsamps = c(25, 50, 100, 200, 250, 500, 750, 1000, 1500, 2000, 2500, 3000)

# Define number of 

# Make array to hold results
regS = array(NA, dim=c(nrow(fia_geo), length(scales), length(nsamps)), 
	dimnames=list(yrplot.id=fia_geo$yrplot.id, scale=paste('R',scales, sep=''),nsamp = nsamps))

nobs = array(NA, dim=c(nrow(fia_geo), length(scales)), 
	dimnames=list(yrplot.id=fia_geo$yrplot.id, scale=paste('R',scales, sep='')))


# Define records to use
recs = allsp_sp

# Loop through scales
for(r in scales){
		
	# Loop through FIA plots
	for(i in 1:nrow(model_data)){
		usepoint = fia_geo[fia_geo$yrplot.id==rownames(model_data)[i],]
		userecs = find_recs(usepoint, recs, r)
		N = nrow(userecs)
		if(N>0){
			comm = tab_species(userecs$scientificName)
			rich = rarefy(comm, nsamps)
		} else {
			rich = rep(0, length(nsamps))
		}
		rich[nsamps>N] = NA
		this_r = paste('R',r,sep='')
		regS[i,this_r,] = rich
		nobs[i,this_r] = N
	}
	save(regS, nobs, file='regS_across_scales.RData')
	print(paste('Done with',r))
}

save(regS, nobs, file='regS_across_scales.RData')








