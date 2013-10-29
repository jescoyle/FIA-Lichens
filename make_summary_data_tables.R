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
