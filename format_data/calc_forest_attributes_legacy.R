# This script calculates various forest structure attricbutes for the lichen legacy plots


#########################################################
### Read in Data
setwd('./UNC/Projects/FIA Lichen')

load('./Data/LegacyData/legacy_lichen_data_parsed_2012-05-15.Rdata')

setwd('C://Users/jrcoyle/Documents/UNC/Projects/FIA Lichen')
master.data.tree = read.csv('./Data/TreeData/master_data_tree_legacy.csv')

treetraits = read.csv('./Data/TreeData/tree_traits_all.csv')
treetraits = subset(treetraits, needed)
treetraits$LogSeed = log10(treetraits$Seed)

#######################################################################################################
### Calculate Forest Characteristics ###
mydirectory = 'C:/Users/jrcoyle/Documents/UNC/Projects/FIA Lichen'
setwd(mydirectory)

options(stringsAsFactors=F)

# Read in tree trait data
#treetraits = read.csv('./Data/TreeData/REF_SPECIES.CSV')
#treetraits = treetraits[,1:4]
#woodbark = read.csv('./Data/TreeData/tree_wood_bark.csv')
#woodbark = woodbark[,c(1,seq(2,ncol(woodbark),2))]
#treetraits = merge(treetraits, woodbark, by='SPCD')

# Read in DBDGS trait data

#dbdgstraits = read.csv('C://Users/jrcoyle/Documents/UNC/DBDGS/SVN/Trait Data/derived_data/FIATraitData_EastSubset_LooseyGoosey.csv', row.names=1)
#dbdgstraits = dbdgstraits[,c('Seed','Nitrogen')]
#dbdgstraits$SPCD = rownames(dbdgstraits)

#traits = merge(treetraits, dbdgstraits, all=T)

# Mark which species are needed for legacy and for regular data
#legabun = read.csv('./Data/TreeData/fia_abundance_matrix_stems_legacy_2012-05-10.csv', row.names=1, check.names=F)
#regabun = read.csv('./Data/TreeData/fia_abundance_matrix_2012-01-03.csv', row.names=1, check.names=F)
#splist = unique(c(names(legabun), names(regabun)))
#traits$needed = traits$SPCD %in% splist

#write.csv(traits,'tree_traits_all.csv')

######## Need to find out whether legacy FIA plots were measured using same layout as regular FIA plots

# Plots sampled 1994 - 1999
plotXyr = table(plot$plot.id, plot$YEAR)

table(plot$STATE_ABBREV, plot$YEAR)




#################
### Functions ###

library('FD')
library('picante')

# A function that calculates the abundance of each species for a particular plot.
#  Seedlings are not counted.

calc.plot.abun = function(use.trees){
	
	# scaling constant to convert counts/area on microplot to counts/area on subplot
	AREA.SCALING.FACTOR = (24^2)/(6.8^2)
	
	return.data = c()
	
	# generate species list for subplot
	tree.species = unique(use.trees$FHM_SPECIES)
		
	# catches the exception when there are no trees on plot
	if(length(tree.species)>0){
	
	# calculate stats for each species
	for(s in tree.species){
	
		woodland=F	
		diams = ifelse(is.na(use.trees$DBH)[1], 'PREVIOUS_DBH', 'DBH')

		if(is.na(use.trees[1,diams])){

			diams = ifelse(is.na(use.trees$DRC)[1], 'PREVIOUS_DRC', 'DRC')
	
			woodland=T
		}


			
		# big trees measures on 24 ft^2 radius subplot, small trees and seedlings measured on 6.8 ft^2 microplot
		small.trees = subset(use.trees,(use.trees$FHM_SPECIES == s)&(use.trees[,diams] < 5))
		big.trees = subset(use.trees, (use.trees$FHM_SPECIES == s)&(use.trees[,diams] >= 5))
		
		# small trees are scaled up as if they were measured on the whole subplot
		n.stems = ifelse(woodland,
				sum(small.trees$NBR_STEMS)*AREA.SCALING.FACTOR + sum(big.trees$NBR_STEMS),
				nrow(small.trees)*AREA.SCALING.FACTOR + nrow(big.trees))
		
		
		# in cm^2
		basal.area = (2.54^2)*(ifelse(nrow(big.trees)>0,sum(sapply(big.trees[,diams], function(x) pi*((x/2)^2))),0) + 
				ifelse(nrow(small.trees)>0,sum(sapply(small.trees[,diams], function(x) pi*((x/2)^2))),0)*AREA.SCALING.FACTOR)

		return.data = rbind(return.data, data.frame(FHM_SPECIES = s, n.stems=n.stems, basal.area=basal.area))

	} # closes species loop
	} # closes if statement for that excludes instances where there are no trees nothe subplot
	
	return(return.data)

}# closes function


# A function that calculates the Simpson index for each plot
#	use.data : a dataframe containing abun.data for one plot
#	abun.column : the column in the data to be used as a proxy for abundance

calc.simpson = function(use.data, abun.column){
	
	p = sapply(use.data[,abun.column], function(y) y/sum(use.data[,abun.column]))
	
	D = sum(p^2)

	D
}

# A function that calculates Rao's quadratic entropy for a set of communities.
#	abundances :  matrix of species abundances in communities.  colnames must equal species names as in rownames of traits.
#	traits : matrix of species mean trait values S X T.  rownames must equal species names as in colnames of abundances.
#	standard : whether to standardize traits to 0 and unit variance

calc.Rao = function(traits, abundances, standard){

require(picante)
require(FD)

## Put trait matrix and abudance matrix in the same order
if(is.null(dim(traits))){

	if(length(traits)!=ncol(abundances)){ stop("Error: traits and abundances have different numbers of species.") }

	if(is.null(names(traits))|is.null(colnames(abundances))){ stop("Error: Species names not in matrices.  Matrices cannot be compared.") }

	if(sum(!(colnames(abundances)%in%names(traits)))>0){ stop("Error: Species names don't match. Species missing from trait matrix.") }

	traits = traits[colnames(abundances)]

} else {

	if(nrow(traits)!=ncol(abundances)){ stop("Error: traits and abundances have different numbers of species.") }

	if(is.null(rownames(traits))|is.null(names(abundances))){ stop("Error: Species names not in matrices.  Matrices cannot be compared.") }

	if(sum(!(names(abundances)%in%rownames(traits)))>0){ stop("Error: Species names don't match. Species missing from trait matrix.") }

	traits = traits[names(abundances),]
}

temp = dbFD(traits, abundances, w.abun=T, stand.x=standard, calc.FRic=F, scale.RaoQ=F, calc.FGR=F, calc.CWM=F, calc.FDiv=F)

temp$RaoQ

} # closes calc.Rao function

make.plotid = function(x){
	paste(x$STATECD, x$COUNTYCD, x$PLOT, sep='_')
}

# note that in FIA, I have used MEASYEAR
make.yrplotid = function(x){
	paste(x$INVYR, make.plotid(x), sep='_')
}
#################################

##########################################
### Code - Calculate Forest Attributes ###

options(stringsAsFactors=F)

## Calculate tree abundance data on plots (plots without trees will not be included)
abun.data = data.frame()

plots = plot; trees = tree # Use for legacy data

for(p in plots$yrplot.id){
	use.trees = subset(trees, yrplot.id == p)
	
	if(nrow(use.trees)>0){
		this.abun.data = calc.plot.abun(use.trees)
	
		this.abun.data$yrplot.id = p

		abun.data = rbind(abun.data, this.abun.data)
	} else {
		#abun.data = rbind(abun.data, c(rep(NA,ncol(abun.data)-1), yrplot.id = p))
	}
	
}

#write.csv(abun.data, './Data/TreeData/tree_abundances.csv', row.names=F)
write.csv(abun.data, './Data/TreeData/tree_abundances_legacy.csv', row.names=F)

abun.data = read.csv('./Data/TreeData/tree_abundances_legacy.csv')

## Calculate tree diversity on each plot

# Remove non-species
unknown.sp = c(998,999,299)
abun.data.noUnk = subset(abun.data, !(FHM_SPECIES %in% unknown.sp))

# Calculate diversity (trees not IDed to species are assumed to be distinct species)

plot.list = unique(abun.data.noUnk$yrplot.id)

tree.diversity = data.frame()

# This calculation assumes plot are laid out in the same way as the recent FIA plots
for(p in plot.list){
	
	use.data = subset(abun.data.noUnk, yrplot.id == p)
	
	S = length(unique(use.data$FHM_SPECIES))
	D.abun = calc.simpson(use.data, 'n.stems')
	D.area = calc.simpson(use.data, 'basal.area')

	tree.diversity = rbind(tree.diversity, data.frame(yrplot.id = p, S = S, D.abun = D.abun, D.area = D.area))
}

master.data.tree = tree.diversity

names(master.data.tree)[2:4] = sapply(names(master.data.tree)[2:4], function(x) paste(x,'tree',sep='.'))

# Get max diameter for each plot
master.data.tree$maxDiam = sapply(master.data.tree$yrplot.id, function(x){
	max(trees[trees$yrplot.id==x,c('DBH','PREVIOUS_DBH','DRC')], na.rm=T)
})

# Get the number of trees > 5in in each plot
trees$DIA = ifelse(is.na(trees$DRC),ifelse(is.na(trees$DBH), trees$PREVIOUS_DBH, trees$DBH),trees$DRC)

tree.diamat = table(trees$yrplot.id, trees$TREE_TYPE) # type 1 = DBH>=5in, type 2 = DBH < 5in
tree.diamat = tree.diamat[master.data.tree$yrplot.id,] # put in same order as master data
master.data.tree$numTreesBig = tree.diamat[,1]
master.data.tree$numTreesSm = tree.diamat[,2]

# Get the total circumference of all trees in each plot
master.data.tree$totalCirc = sapply(master.data.tree$yrplot.id, function(x){
	use.dias = subset(trees, yrplot.id==x)$DIA

	sum(use.dias*pi, na.rm=T)
})


# Get the proportion of trees that are dead: excludes cut trees
tree.statmat = table(trees$yrplot.id, trees$CURRENT_TREE_HISTORY)
tree.statmat = tree.statmat[master.data.tree$yrplot.id,]
master.data.tree$propDead = apply(tree.statmat, 1, function(x) sum(x[c('5','6')])/sum(x[as.character(c(1:6,8,10,12))]))
master.data.tree$numDead = apply(tree.statmat, 1, function(x) sum(x[c('5','6')]))
master.data.tree$numCut = tree.statmat[,'7']

master.data.tree$PIE.stems.tree = 1-master.data.tree$D.abun.tree
master.data.tree$PIE.ba.tree = 1-master.data.tree$D.area.tree


## Calculate tree trait diversity

# make matrices
traits<-treetraits[,c('wood_moist_pct','bark_moist_pct','wood_SG','bark_SG','LogSeed')]
rownames(traits) = treetraits$SPCD
traits = as.matrix(traits)
abundances = reshape(abun.data.noUnk[,c('yrplot.id','FHM_SPECIES','n.stems')], 
	v.names='n.stems', idvar='yrplot.id', timevar='FHM_SPECIES', direction='wide')
row.names(abundances)<-abundances$yrplot.id
abundances = abundances[,-1]
names(abundances) = unlist(lapply(strsplit(names(abundances), '.', fixed=T), function(x) x[3]))
abundances[is.na(abundances)]<-0

# For basal area
abundances = reshape(abun.data.noUnk[,c('yrplot.id','FHM_SPECIES','basal.area')], 
	v.names='basal.area', idvar='yrplot.id', timevar='FHM_SPECIES', direction='wide')
row.names(abundances)<-abundances$yrplot.id
abundances = abundances[,-1]
names(abundances) = unlist(lapply(strsplit(names(abundances), '.', fixed=T), function(x) x[3]))
abundances[is.na(abundances)]<-0

# Checking that all species are in the trait table
names(abundances)[which(!(names(abundances) %in% rownames(traits)))]-> missing.trees
## There were 3 species not in trait table.  This is what I did with them:
subset(trees, FHM_SPECIES %in% missing.trees)[,]

# 376: western paper birch, 4 trees in two plots.  Is a subspecies of paper birch (375).  Changed to 275.
# 477: hairy mountain mahogany, 1 tree. May be subspecies of curlleaf mountain-mahogany (475). Changed to 475.
# 476: alder-leaf mountain mahogany, 3 trees in 3 plots. Kept name the same, but gave same traits as curlleaf mountain mahogany (475).

# Checking that codes are the same
unique(trees[,c('FHM_SPECIES','SPECIES_COMMON_NAME')])-> ref_fhm
ref_fhm[order(ref_fhm$FHM_SPECIES),]

traits = traits[names(abundances),]

# Checking that there are no trees with zero abundance
which(colSums(abundances) == 0)
subset(trees, FHM_SPECIES==299)

write.csv(abundances, paste('./Data/TreeData/fia_abundance_matrix_stems_legacy_',Sys.Date(),'.csv',sep=''), row.names=T)
write.csv(abundances, paste('./Data/TreeData/fia_abundance_matrix_ba_legacy_',Sys.Date(),'.csv',sep=''), row.names=T)
write.csv(traits, paste('./Data/TreeData/fia_trait_matrix_legacy_',Sys.Date(),'.csv',sep=''), row.names=T)

abundances = read.csv('./Data/TreeData/fia_abundance_matrix_stems_legacy_2012-05-10.csv', row.names=1, check.names=F)
abundances = read.csv('./Data/TreeData/fia_abundance_matrix_ba_legacy_2012-05-10.csv', row.names=1, check.names=F)


## Calculate diversity indices - did for both abundances based on n.stems and on basal.area

# For wood/bark traits
# Standardize moisture percentages by 200% and specific gravities by 1
traits.std = traits

# Standardize to fall between 0 and 1
for(i in colnames(traits)){
	traits.std[,i]=(traits[,i]-min(treetraits[,i], na.rm=T))/diff(range(treetraits[,i], na.rm=T))
}

sapply(colnames(traits.std), function(n) calc.Rao(traits.std[,n], abundances, standard=F))->rao.abun
colnames(rao.abun) = paste(colnames(rao.abun),'rao.abun',sep='.')

presences = ifelse(abundances>0,1,0)

sapply(colnames(traits.std), function(n) calc.Rao(traits.std[,n], presences, standard=F))->rao.pres
colnames(rao.pres) = paste(colnames(rao.pres),'rao.pres',sep='.')

tree.traitDiv = data.frame(yrplot.id = rownames(rao.abun),rao.abun,rao.pres)

# Calculate mean trait values
mean.traits = (as.matrix(abundances) %*% as.matrix(traits))/rowSums(abundances)
tree.traitDiv = cbind(tree.traitDiv, mean.traits)

write.csv(tree.traitDiv, './Data/TreeData/tree_FD_stems_legacy.csv', row.names=F)
write.csv(tree.traitDiv, './Data/TreeData/tree_FD_ba_legacy.csv', row.names=F)

treeDiv.stems = read.csv('./Data/TreeData/tree_FD_stems_legacy.csv')
treeDiv.ba = tree.traitDiv

plot(treeDiv.stems$bark_moist_pct~treeDiv.ba$bark_moist_pct); abline(0,1,col=2)

names(treeDiv.stems) = gsub('abun', 'stems', names(treeDiv.stems))
names(treeDiv.stems)[names(treeDiv.stems) %in% colnames(traits)] = paste(names(treeDiv.stems)[names(treeDiv.stems) %in% colnames(traits)], 'stems', sep='.')

names(treeDiv.ba) = gsub('abun', 'ba', names(treeDiv.ba))
names(treeDiv.ba)[names(treeDiv.ba) %in% colnames(traits)] = paste(names(treeDiv.ba)[names(treeDiv.ba) %in% colnames(traits)], 'ba', sep='.')

treeDiv = merge(treeDiv.stems, treeDiv.ba, all=T)
write.csv(treeDiv, './Data/TreeData/tree_FD_legacy.csv', row.names=F)


## Merge with master.data
master.data.tree = merge(master.data.tree, treeDiv, by='yrplot.id', all=T)


### Calculate tree size diversity (only use trees > 5 in)

sapply(master.data.tree$yrplot.id, function(x){
	
	these.trees = subset(trees, yrplot.id==x)[,'DIA']
	big.trees = these.trees[these.trees>=5]

	if(length(big.trees)>0){ mean(as.numeric(dist(big.trees)), na.rm=T) } else {NA}

}) -> master.data.tree$diamDist.mean


## Calculate total basal area
cbind(master.data.tree, t(sapply(master.data.tree$yrplot.id, function(x){
	these.trees = subset(trees, yrplot.id==x)

	colSums(calc.plot.abun(these.trees)[,c('n.stems','basal.area')])
}))) -> master.data.tree


# function that calculates the surface area of a tree based on its diameter and height
coneArea = function(d,h){
	r = d/2
	s = sqrt((r^2) + (h^2))
	pi*r*s
}

## Calculate variation in light penetration and average light penetration
sapply(master.data.tree$yrplot.id, function(x){
	these.trees = subset(trees, yrplot.id==x)
	big.trees = subset(these.trees, DIA >= 5)

	if(nrow(big.trees)>0){
		
	# If 50% of tree >5in have CROWN_DENSITY data, then do calculations
	if(sum(!is.na(big.trees$CROWN_DENSITY))/nrow(big.trees) > .5){
		light.mean = mean(big.trees$CROWN_DENSITY, na.rm=T)
		light.var = var(big.trees$CROWN_DENSITY, na.rm=T)
		lightDist.mean = mean(as.numeric(dist(big.trees$CROWN_DENSITY)), na.rm=T)		
	} else {
		light.mean = NA
		light.var = NA
		lightDist.mean = NA
	}
	} else {
		light.mean = NA
		light.var = NA
		lightDist.mean = NA
	}

	c(light.mean, light.var, lightDist.mean)
	
}) -> tree.light

tree.light = data.frame(t(tree.light))
names(tree.light)= c('light.mean','light.var','lightDist.mean')

#master.data.tree = cbind(master.data.tree, tree.light) # NaN results from the distance function not working ononly one big tree
master.data.tree[,c('light.mean','light.var','lightDist.mean')] = tree.light # Used for updating master data if re-calculating lighht values
master.data.tree$lightDist.mean[is.nan(master.data.tree$lightDist.mean)]<-NA

summary(trees$FOLIAGE_TRANSPARENCY)

hist(subset(trees, is.na(trees$FOLIAGE_TRANSPARENCY))[,'DIA'])
hist(subset(trees, !is.na(trees$FOLIAGE_TRANSPARENCY))[,'DIA'])


## Add forest type code
cond$yrplot.id = paste(cond$year, cond$plot.id, sep='_')

# find plots without conditions
missing_cond = !(master.data.tree$yrplot.id %in% cond$yrplot.id)
master.data.tree$yrplot.id[missing_cond] # all from 1994, just ignore

FORTYPCD = sapply(master.data.tree$yrplot.id, function(x){
	use_conds = subset(cond, yrplot.id==x)
	
	FORTYPCD=NA
	
	if(nrow(use_conds)==0){ } else {

	if(nrow(use_conds)==1){
		FORTYPCD = use_conds$FOREST_TYPE
	} else {
		FORTYPCD = use_conds[use_conds$CONDITION_CLASS_NBR==1, 'FOREST_TYPE']
	}

	}

	FORTYPCD
})

master.data.tree$FORTYPCD = as.numeric(FORTYPCD)



## Write out master forest data
write.csv(master.data.tree, './Data/TreeData/master_data_tree_legacy.csv', row.names=F)

############################################################
### Calculate Stand Age - I DID NOT DO THIS YET

## Get stand age from raw legacy data
filelist_raw = list.files('./Data/LegacyData')
stateXyrs = t(sapply(strsplit(filelist_raw, '_'), function(x) c(x[1], substr(x[length(x)],1,4))))

standage = data.frame()
for(i in 1:nrow(stateXyrs)){
	filename = paste(stateXyrs[i,1],'_PLOT_COND_CLASS_',stateXyrs[i,2],'.CSV', sep='')
	
	usefile = read.csv(paste('./Data/LegacyData/',filename, sep=''))
	
	make.yrplotid(usefile) # STOPPED HERE because cond tables don't have the necessary columns to make plot ids

}





