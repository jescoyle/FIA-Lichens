## This script calculates functional diversity of lichens on FIA plots use LIAS trait data.
options(stringsAsFactors=F)

setwd('./UNC/Projects/FIA Lichen')

###########################################
### Read in Data

# Site X species matrix of lichen species abundances (by abundance categories)
siteXsp = read.csv('./Data/siteXsp_matrix_current+legacy.csv', row.names=1, check.names=F)

# Lichen species list from FIA that I provided to LIAS
lichenSp = read.csv('../Lichen Traits/fia_lichen_master_list.csv')

# Combine Cladonia squamosa (1236) with Cladonia squamosa var. subsquamosa (1238)
siteXsp[,'1236'] = apply(siteXsp[,c('1236','1238')],1,combine.abundance)
siteXsp = siteXsp[, colnames(siteXsp)!='1238']

# Combine Cladonia gracilis (1260) and Cladonia gracilis ssp. turbinata (1261) because LIAS didn't give me back different traits for these taxa.
siteXsp[,'1260'] = apply(siteXsp[,c('1260','1261')],1,combine.abundance)
siteXsp = siteXsp[, colnames(siteXsp)!='1261']

# Combine Bryoria mystery oline sp (627) and Bryoria (600)
siteXsp[,'600'] = apply(siteXsp[,c('627','600')],1,combine.abundance)
siteXsp = siteXsp[, colnames(siteXsp)!='627']

# Fix one species whose name changed in transit:
lichenSp[(lichenSp$Genus=='Cetraria')&(lichenSp$Species=='americana'),'Genus'] = 'Tuckermannopsis'

# Remove unidentified lichens (9999) and Cladina-form (1180)
siteXsp = siteXsp[, !(colnames(siteXsp) %in% c('1180','9999'))]

# Only use species in siteXsp matrix
lichenSp = subset(lichenSp, LICH_SPPCD %in% colnames(siteXsp))

# Trait data from LIAS
traits_raw = read.delim('../Lichen Traits/LIAS Traits/LIAS_traits.txt', check.names=F)

# Merge species lists to update my species list with new LIAS species codes
old = paste(lichenSp$Genus, lichenSp$Species)
new = paste(traits_raw$Genus, traits_raw$Species)
old[which(!(old %in% new))] # Check to make sure all names match

traits = merge(lichenSp[,c('LICH_SPPCD','Genus','Species')], traits_raw, all.x=T, all.y=F)
nrow(traits)==ncol(siteXsp)

### Clean up trait data

# Reassign no data value
traits = replace(traits, traits=='n.d.', NA)

colnames(traits)
traits[,11:18] = replace(traits[,11:18], traits[,11:18]==0, NA)

traits[,32:34] = replace(traits[,32:34], traits[,32:34]=='', NA)


colnames(traits)[11:34] = c('pruinose','cilia','isidia','soredia','ascomata','primary_photo','secondary_photo',
	'metabolites','atranorin','usnic_acid','pulvinic_acid','pulvinic_dilactone','pulvinic_acid_lactone',
	'unk_pulvinic_acid','pulvinic_acid_derivative','parietin','parietinic_acid','bullate','macculate',
	'pseudocyphellate','hairy','ascospore_length','ascospore_width','conidia_length')

write.csv(traits, './Data/LichenTraits/FIA_lichen_species_traits_LIAS.csv', row.names=F)
write.csv(siteXsp, './Data/siteXsp_matrix_current+legacy_cleaned.csv', row.names=T)

######################################################
### Recode traits

traits = read.csv('./Data/LichenTraits/FIA_lichen_species_traits_LIAS.csv')

## Calculate genus-level traits
genera = subset(traits, Species=='')

for(i in c('pruinose','cilia','isidia','soredia','ascomata','primary_photo','secondary_photo','metabolites')){
	sapply(genera$Genus, function(x){
		usesp = subset(traits, Genus==x)
		avg.3state(usesp[,i], 1:3)
	}) -> genus_avg

	genera[,i] = genus_avg
}


for(i in c('atranorin','usnic_acid','pulvinic_acid','pulvinic_dilactone','pulvinic_acid_lactone',
	'unk_pulvinic_acid','pulvinic_acid_derivative','parietin','parietinic_acid','bullate','macculate',
	'pseudocyphellate','hairy','ascospore_length','ascospore_width','conidia_length')){
	
	sapply(genera$Genus, function(x){
		usesp = subset(traits, Genus==x)
		avg.3state(usesp[,i], 0:2)
	}) -> genus_avg

	genera[,i] = genus_avg
}

for(i in c('ascospore_length','ascospore_width','conidia_length')){
	sapply(genera$Genus, function(x){
		usesp = subset(traits, Genus==x)
		tot.range(usesp[,i])
	}) -> genus_avg

	genera = cbind(genera, t(genus_avg))
}

sapply(c('ascospore_length','ascospore_width','conidia_length'), function(x){
	c(paste(x, 'min', sep='_'),paste(x, 'max', sep='_'))
})->newnames
colnames(genera)
colnames(genera)[35:40] =  as.character(newnames)

## Unstack numeric trait ranges from single column into min and max columns

for(i in c('ascospore_length','ascospore_width','conidia_length')){
	
	newranges = split.range(traits[,i])

	traits = cbind(traits, newranges)
}
sapply(c('ascospore_length','ascospore_width','conidia_length'), function(x){
	c(paste(x, 'min', sep='_'),paste(x, 'max', sep='_'))
})->newnames
colnames(traits)
colnames(traits)[35:40] =  as.character(newnames)


## Remove old genera w/o traits and replace with genus level traits
traits = subset(traits, Species!='')
traits = rbind(traits, genera)

### Recode traits that are 1:3 as 0:2
traits[,c('pruinose','cilia','isidia','soredia','ascomata',
	'primary_photo','secondary_photo','metabolites')] <-
traits[,c('pruinose','cilia','isidia','soredia','ascomata',
	'primary_photo','secondary_photo','metabolites')] - 1

summary(traits)


write.csv(traits, './Data/LichenTraits/FIA_lichen_species_traits_LIAS.csv', row.names=F)


###############################################
### Calculate functional diversity
library(stringr) # str_trim white space remover
library(cluster) # daisy
library(ade4) # dist.ktab Gower distance extension
library(abind)
library(FD) # dbFD calculates functional diversity

load('./Data/LichenFD.Rdata') # Imports lichen FD calculated across all traits


siteXsp = read.csv('./Data/siteXsp_matrix_current+legacy_cleaned.csv', row.names=1, check.names=F)
traits = read.csv('./Data/LichenTraits/FIA_lichen_species_traits_LIAS.csv')

## Rename siteXsp columns with species names and
## Check which species are not included in trait distance matrix
rownames(traits) = as.character(traits$LICH_SPPCD)
missingsp = colnames(siteXsp)[!(colnames(siteXsp) %in% traits$LICH_SPPCD)] #none missing, safe to rename

# Rename siteXsp columns
colnames(siteXsp) = str_trim(paste(traits[colnames(siteXsp),'Genus'], traits[colnames(siteXsp),'Species']))

# For now, use midpoint of numeric traits
traits$ascospore_length_mid = apply(traits[,c('ascospore_length_min','ascospore_length_max')], 1, mean)
traits$ascospore_width_mid = apply(traits[,c('ascospore_width_min','ascospore_width_max')], 1, mean)
traits$conidia_length_mid = apply(traits[,c('conidia_length_min','conidia_length_max')], 1, mean)


# Define list of traits and whether they are binary or numeric
# Removed pulvinic_acid_lactone, unk_pulvinic_acid, and pulvinic_acid_derivative b/c no variation among these species
traitdf = data.frame(
	trait=c('pruinose','cilia','isidia','soredia','ascomata',
		'secondary_photo','bullate','macculate','pseudocyphellate',
		'hairy','metabolites','atranorin','usnic_acid',
		'pulvinic_acid','pulvinic_dilactone','parietin',
		'parietinic_acid', 'ascospore_length_mid', 'ascospore_width_mid',
		'conidia_length_mid'),
	type = c('binary','binary','binary','binary','binary',
		'binary','binary','binary','binary',
		'binary','binary','binary','binary',
		'binary','binary','binary',
		'binary', 'numeric', 'numeric',
		'numeric'),
	category = c('morph','morph','reproduce','reproduce','reproduce',
		'photo','morph','morph','morph',
		'morph','chem','chem','chem',
		'chem','chem','chem',
		'chem', 'reproduce', 'reproduce',
		'reproduce')
)

write.csv('lias_trait_types.csv', row.names=F)

summary(traits[,traitdf$trait])


## Make species X trait matrix
spXtrait = traits[,traitdf$trait]
rownames(spXtrait) = str_trim(paste(traits$Genus, traits$Species))

## Remove species that have no trait information:
badsp = c('Cetraria viridis','Cladonia brevis','Physcia leptalea','Melanelia','Vulpicida')
spXtrait[badsp,] # No information on any of these species
for(i in badsp){
	print(i)
	print(table(siteXsp[,i]))
}
# Cetraria viridis, Cladonia brevis, Vulpicida occur in one plot at abundances of 1,2,1
# Physcia leptalea occurs in 4 plots, 3 of which at high abundance (3)
# Melanelia occurs in 30 plots, 18 of which are high abundance (3), but it doesn't have any daugther species to get traits from.

spXtrait = spXtrait[!(rownames(spXtrait) %in% badsp),]
siteXsp = siteXsp[,!(colnames(siteXsp) %in% badsp)]

## Save matrices
write.csv(siteXsp,'./Data/siteXsp_matrix_foranalysis.csv', row.names=T)
write.csv(spXtrait,'./Data/spXtrait_matrix_foranalysis.csv', row.names=T)

spXtrait = read.csv('./Data/spXtrait_matrix_foranalysis.csv', row.names=1)
siteXsp = read.csv('./Data/siteXsp_matrix_foranalysis.csv', row.names=1)

### Format trait data from calculating distance matrix
# Define dataframes of binary and numeric data
bindat = spXtrait[,subset(traitdf, type=='binary')$trait]
numdat = spXtrait[,subset(traitdf, type=='numeric')$trait]

## Make binary traits into multistate nominal variables (since some are coded as 'both')
sapply(colnames(bindat), simplify = "array", FUN=function(x){
	usedata = bindat[,x]
	newdata = sapply(usedata, function(y){
		if(is.na(y)){
			absent = NA
			present = NA
		} else {
			absent = ifelse(y %in% c(0,2), 1, 0)
			present = ifelse(y %in% 1:2, 1, 0)
		}
		c(absent, present)
	})
	
	newdata = t(newdata)
	colnames(newdata) = c('absent','present')
	rownames(newdata) = rownames(bindat)
	
	newdata
}) -> binary_trait_arr

# Unstack array into 2-D data frames of different categories of binary traits
bintraits = subset(traitdf, type=='binary')

morph_traits_arr = binary_trait_arr[,,subset(bintraits, category=='morph')$trait]
morph_traits_df = morph_traits_arr[,,1]
for( i in 2:dim(morph_traits_arr)[3]){
	morph_traits_df = cbind(morph_traits_df, morph_traits_arr[,,i])
}
morph_traits_df = prep.binary(data.frame(morph_traits_df), 
	col.blocks = rep(2, dim(morph_traits_arr)[[3]]),
	label = dimnames(morph_traits_arr)[[3]])

rep_traits_arr = binary_trait_arr[,,subset(bintraits, category=='reproduce')$trait]
rep_traits_df = rep_traits_arr[,,1]
for( i in 2:dim(rep_traits_arr)[3]){
	rep_traits_df = cbind(rep_traits_df, rep_traits_arr[,,i])
}
rep_traits_df = prep.binary(data.frame(rep_traits_df), 
	col.blocks = rep(2, dim(rep_traits_arr)[[3]]),
	label = dimnames(rep_traits_arr)[[3]])

chem_traits_arr = binary_trait_arr[,,subset(bintraits, category=='chem')$trait]
chem_traits_df = chem_traits_arr[,,1]
for( i in 2:dim(chem_traits_arr)[3]){
	chem_traits_df = cbind(chem_traits_df, chem_traits_arr[,,i])
}
chem_traits_df = prep.binary(data.frame(chem_traits_df), 
	col.blocks = rep(2, dim(chem_traits_arr)[[3]]),
	label = dimnames(chem_traits_arr)[[3]])

photo_traits_df = binary_trait_arr[,,subset(bintraits, category=='photo')$trait]
photo_traits_df = prep.binary(data.frame(photo_traits_df), 
	col.blocks = 2,
	label=subset(bintraits, category=='photo')$trait)

## Define quantitative traits
quant_traits = numdat

## Create ktable and calculate Gower distance
# Note that this weights trait categories equally, but not the traits themselves
ktab = ktab.list.df(list(quant_traits, morph_traits_df, rep_traits_df, chem_traits_df, photo_traits_df),
	tabnames=c('spores','morphology','reproduction','chemistry','photobiont'))

ktab

# Calculate Gower distance
traitdist = dist.ktab(ktab, type=c('Q','B','B','B','B'), option="scaledBYrange")

hist(traitdist)


# Reorder into same as trait distance matrix
colnames(siteXsp) %in% labels(traitdist)
siteXsp = siteXsp[,labels(traitdist)]

### Calculate Functional Diversity

# Across all traits
lichenFD = dbFD(traitdist, siteXsp, w.abun=F)

save(lichenFD, file='./Data/lichenFD.Rdata')

lichenFD_df = data.frame(yrplot.id = names(lichenFD$RaoQ), raoQ = lichenFD$RaoQ,
	fric = lichenFD$FRic, feve = lichenFD$FEve, fdiv = lichenFD$FDiv, fdis = lichenFD$FDis)

write.csv(lichenFD_df, './Data/lichen_FD.csv', row.names=F)


##########################################
### Explore lichen FD

master = merge(master, lichenFD_df, all.x=T)

plot(raoQ~rh, data=master)



##############################################
### Calculate trait means across plots

spXtrait = read.csv('./Data/spXtrait_matrix_foranalysis.csv', row.names=1)
siteXsp = read.csv('./Data/siteXsp_matrix_foranalysis.csv', row.names=1)
lichenFD_df = read.csv('./Data/lichen_FD.csv')

## Just based on presence absence
siteXsp_pres = siteXsp>0

## Caluculate the proportion of species on each site that have each binary trait.
## NA is not counted in total
## Both is added to both

sapply(subset(traitdf, type=='binary')$trait, function(x){
	calc.bintraitmean(spXtrait[,x], siteXsp_pres)
}) -> bintrait_means

sapply(subset(traitdf, type=='numeric')$trait, function(x){
	calc.traitmean(spXtrait[,x], siteXsp_pres)
}) -> numtrait_means

lichen_trait_means = cbind(bintrait_means, numtrait_means)
rownames(lichen_trait_means) = rownames(siteXsp)

rownames(lichen_trait_means) == lichenFD_df$yrplot.id
lichenLIAS = cbind(lichenFD_df, lichen_trait_means)

write.csv(lichen_trait_means, './Data/fia_lichen_trait_means.csv', row.names=T)
write.csv(lichenLIAS, './Data/fia_lichen_LIAS_means_diversity.csv')



##############################################
### Functions

combine.abundance = function(x){
	
	final.abun = NA
	
	x = x[order(x)]
	
	if( sum(x >= 4) > 0 ){ final.abun = 4 } else {
		if( sum(x >= 3) > 0 ){ final.abun = 3 } else {
			if( length(x) < 3 ){
				final.abun = sum(x)
			}
		}
	}
	
	final.abun

}

# Compute a genus-level average trait from:
# x : a vector of species level traits
# states : a vector of length 3 with the characters representing FALSE, TRUE, BOTH
avg.3state = function(x, states){
	x = factor(x, levels=states)
	
	counts = table(x)
	counts = counts>0
	
	if(sum(counts)==0){
		output = NA
	} else {
		output = ifelse(sum(counts)>1,states[3],states[which(counts)])
	}

	output
}

# Compute a genus-level total trait range from:
# x : a vector of trait ranges num-num or num
tot.range = function(x){
	
	splitlist = strsplit(x, '-')
	
	low = NA
	hi = NA
	
	for(n in 1:length(splitlist)){
		thisrange = splitlist[[n]]
		
		low = c(low, ifelse(length(thisrange)==2,thisrange[1],thisrange))
		hi = c(hi, ifelse(length(thisrange)==2,thisrange[2],thisrange))
	}

	if(sum(!is.na(low))==0){ output = c(NA,NA) } else {
		output = c(min(as.numeric(low), na.rm=T), max(as.numeric(hi), na.rm=T))
	}

	output
}

# Split trait ranges into min and max
split.range = function(x){
	usesplit = strsplit(x, '-')
	
	sapply(usesplit, function(y){
		low = as.numeric(ifelse(length(y)==2, y[1], y))
		hi = as.numeric(ifelse(length(y)==2, y[2],y))

		c(low,hi)
	}) -> useranges

	t(useranges)	
}


# Calculates distance between two instances of a binary trait when a trait can have both states
# x : vector of length 2 with the trait values
# states : vector of length 3 with the characters representing FALSE, TRUE, BOTH

dist.3state = function(x, states){
	
	sum(x==states[3])

}

# Calculates the proportion of species in a site that have a binary trait
# x is vector of trait values
# abun is the siteXsp abundance matrix
calc.bintraitmean = function(x, abun){
	x_pres = x %in% 1:2
	x_abs = x %in% c(0,2)
	x_notNA = !is.na(x)

	abun = as.matrix(abun)

	(abun %*% x_pres) / ((abun %*% x_pres) + (abun %*% x_abs))
#	(abun %*% x_pres) / (abun %*% x_notNA)
}

calc.traitmean = function(x, abun){
	x_noNA = x
	x_noNA[is.na(x)] <- 0

	x_notNA = !is.na(x)

	abun = as.matrix(abun)

	(abun %*% x_noNA) / (abun %*% x_notNA)

}

