### This script adds geographic coordinates to legacy lichen data from a file
### that Sarah Jovan sent me on 5/6/2013


options(stringsAsFactors=F)

setwd('C:/Users/jrcoyle/Documents/UNC/Projects/FIA Lichen')


master = read.csv('./Data/master_data.csv')
jovan_plots = read.delim('./Data/From Sarah Jovan/AllLichenPlots_forJES.txt')
legacy = read.csv('./Data/LegacyData/legacy_lichen_plots.csv')

legacy_plots = unique(legacy[,c('yrplot.id','P3ID')])
nolatlon = subset(master, is.na(LAT))

nolatlon = merge(nolatlon, legacy[,c('yrplot.id','P3ID')], all.x=T)

nolatlon[is.na(nolatlon$P3ID),'plot.id']->notinleg # Just 2 plots to check by hand, but these are excluded later anyways


jovan_plots$yrplot.id = paste(jovan_plots$MEASYEAR, jovan_plots$STATECD, 
	jovan_plots$COUNTYCD, jovan_plots$P3ID, sep='_')

jovan_plots$plot.id = paste(jovan_plots$STATECD, 
	jovan_plots$COUNTYCD, jovan_plots$P3ID, sep='_')

new_coords = unique(subset(jovan_plots, plot.id %in% nolatlon$plot.id)[,c('plot.id','FZ_LAT','FZ_LON')])

length(unique(nolatlon$plot.id))==nrow(nolatlon) # Check that plot.id is a unique indetifier
rownames(nolatlon) = nolatlon$plot.id

nolatlon[new_coords$plot.id,c('LAT','LON')] = new_coords[,2:3]

plot(nolatlon$LON, nolatlon$LAT)
sum(is.na(nolatlon$LON)) # Only 384 plots without coordinates

# Add coordinates to master data
length(unique(master$plot.id))==nrow(master) # Check that plot.id is a unique indetifier
rownames(master) = master$plot.id
master[new_coords$plot.id,c('LAT','LON')] = new_coords[,2:3]


sum(is.na(master$LON))
plot(master$LON, master$LAT)

write.csv(master, './Data/master_data.csv', row.names=F)


### Find all plots that are in the correct STATE COUNTY and YEAR but do not have P3IDs
nolatlon = subset(master, is.na(LAT))
nolatlon = merge(nolatlon, legacy[,c('yrplot.id','P3ID')], all.x=T)

# A list of plots that may not have matched correctly due to NA in P3ID or COUNTYCD
badID = jovan_plots[is.na(jovan_plots$P3ID)|is.na(jovan_plots$COUNTYCD),]

# A list of  STATE COUNTY YEAR combos to search through
searchareas = unique(nolatlon[,c('STATE','COUNTY','YEAR','P3ID')])
rownames(searchareas) = paste(searchareas$STATE, searchareas$COUNTY, searchareas$P3ID, sep='_')
searchareas$LAT = NA
searchareas$LON = NA
searchareas$foundplot = NA

for(i in 1:nrow(searchareas)){
	parm = as.numeric(searchareas[i,])

	potentials = subset(badID, (STATECD==parm[1])&((MEASYEAR==parm[3])|(is.na(MEASYEAR)))&((COUNTYCD==parm[2])|(is.na(COUNTYCD))))
	
	if(sum(potentials$P3ID==parm[4], na.rm=T)==1){
		searchareas[i,'foundplot'] = T
		searchareas[i,c('LAT','LON')] = potentials[which(potentials$P3ID==parm[4]),c('FZ_LAT','FZ_LON')]
	} else {
		searchareas[i,'foundplot'] = F
	}
}

# These are the plots that matched. Only 31 now without coordinates.
searchareas[searchareas$foundplot,]

# Make a new dataframe of plot locations that will be used to extract environmental data.
master_coords = master[,c('yrplot.id','plot.id','LAT','LON','state.abbr','COUNTYNM')]
master_coords[rownames(subset(searchareas, foundplot)), c('LAT','LON')] = subset(searchareas, foundplot)[,c('LAT','LON')]

write.csv(master_coords, './Data/fia_lichen_plot_locations.csv', row.names=F)

subset(master_coords, is.na(LAT)) # Plots without coordinates.

#####################################################################
## Old Code prior to getting data file with P3IDs on 7-2-2013

sum(nolatlon$yrplot.id %in% jovan_plots$yrplot.id)

sum(nolatlon$plot.id %in% jovan_plots$plot.id)

sum(master$plot.id %in% jovan_plots$plot.id)

sum(jovan_plots$plot.id %in% master$plot.id)

sum(master$yrplot.id %in% jovan_plots$yrplot.id)


head(master[master$yrplot.id %in% jovan_plots$yrplot.id, ])

nolatlon[nolatlon$plot.id %in% jovan_plots$plot.id,]
subset(jovan_plots, plot.id=='16_49_80175')



new_coords = unique(subset(jovan_plots, plot.id %in% nolatlon$plot.id)[,c('plot.id','FZ_LAT','FZ_LON')])


# Define data set to send to Sarah showing the plots I am trying to get
nolatlon = subset(nolatlon, !(plot.id %in% new_coords$plot.id))
nolatlon = subset(nolatlon, YEAR>=1997)
nolatlon = subset(nolatlon, numTreesBig>1)

nolatlon = nolatlon[, c('yrplot.id','plot.id','STATE','COUNTY','P3ID','YEAR','LAT','LON','STATE_NAME','COUNTYNM')]


write.csv(nolatlon, './Data/plots_missing_latlon_for_Sarah.csv', row.names=F)





