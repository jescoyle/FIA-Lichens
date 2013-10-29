
## This script uses slope, aspect and basal area data from Sarah Jovan to calculate site-scale environmental variables for lichen plots.
options(stringsAsFactors=F)
setwd('C:/Users/jrcoyle/Documents/UNC/Projects/FIA Lichen')

master = read.csv('./Data/fia_lichen_master_data.csv')
jovan_plots = read.table('./Data/From Sarah Jovan/AllLichenPlots_forJES.txt', header=T, sep='\t')
legacy = read.csv('./Data/LegacyData/legacy_lichen_plots.csv')
current = read.csv('./Data/FIA_lichen_summary_parsed.csv')

# Substitute Plot numbers for P3IDs when missing
jovan_plots[is.na(jovan_plots$P3ID),'P3ID'] <- jovan_plots[is.na(jovan_plots$P3ID),'PLOT']

jovan_plots$yrplot.id = paste(jovan_plots$MEASYEAR, jovan_plots$STATECD, 
	jovan_plots$COUNTYCD, jovan_plots$P3ID, sep='_')
jovan_plots$plot.id = paste(jovan_plots$STATECD, 
	jovan_plots$COUNTYCD, jovan_plots$P3ID, sep='_')


### Make a table matching master data yrplot.id to jovan data IDs

## Match plots that do not have NA for P3ID or county or state
useplots = subset(jovan_plots, (!is.na(P3ID))&(!is.na(COUNTYCD))&(!is.na(STATECD))) 

matchplots = subset(useplots, yrplot.id %in% master$yrplot.id)[,c('yrplot.id','ID')]

## Deal with plots that have NA in P3ID or county

# A list of plots that have NA in P3ID or COUNTYCD

# A list of  STATE COUNTY YEAR combos to search through
useplots = subset(current, (yrplot.id %in% master$yrplot.id)&!(yrplot.id %in% matchplots$yrplot.id))
searchareas_c = unique(useplots[,c('STATECD','COUNTYCD','MEASYEAR','PLOT')])
useplots = subset(legacy, (yrplot.id %in% master$yrplot.id)&!(yrplot.id %in% matchplots$yrplot.id))
searchareas_l = unique(useplots[,c('STATE','COUNTY','YEAR','P3ID')])
names(searchareas_c) = names(searchareas_l)
searchareas = rbind(searchareas_c, searchareas_l)

rownames(searchareas) = paste(searchareas$YEAR, searchareas$STATE, searchareas$COUNTY, searchareas$P3ID, sep='_')
searchareas$ID = NA
searchareas$foundplot = NA

for(i in 1:nrow(searchareas)){
	parm = as.numeric(searchareas[i,])

	potentials = subset(jovan_plots, (STATECD==parm[1])&((MEASYEAR==parm[3])|(is.na(MEASYEAR)))&((COUNTYCD==parm[2])|(is.na(COUNTYCD))))
	
	if(sum(potentials$P3ID==parm[4], na.rm=T)==1){
		searchareas[i,'foundplot'] = T
		searchareas[i,'ID'] = potentials[which(potentials$P3ID==parm[4]),'ID']
	} else {
		searchareas[i,'foundplot'] = F
	}
}

# These are the plots that matched. 
searchareas[!searchareas$foundplot,] #27 without match
table(searchareas$foundplot)

# Add onto matchplots
matchplots_more = data.frame(yrplot.id = rownames(searchareas), ID = searchareas$ID)
matchplots = rbind(matchplots, matchplots_more) # 2228 rows

# Save matching
write.csv(matchplots, './Data/From Sarah Jovan/match plot IDs.csv', row.names=F)

## Add Aspect, Slope, Elevation and tree composition to master data
matchplots = read.csv('./Data/From Sarah Jovan/match plot IDs.csv')

rownames(jovan_plots) = jovan_plots$ID
rownames(matchplots) = sapply(strsplit(matchplots$yrplot.id, '_'), function(x) paste(x[2],x[3],x[4], sep='_'))

newvars = merge(matchplots, jovan_plots[,c('ID','MEASYEAR','Slope....','Aspect..degrees.','Basal.area.from.conifers....','FZ_ELEV')], all.x=T, all.y=F)
newvars$plot.id = sapply(strsplit(newvars$yrplot.id, '_'), function(x) paste(x[2],x[3],x[4], sep='_'))

rownames(newvars) = newvars$plot.id # B/C yrplot.id in newvars is different from in master
colnames(newvars)[3:7] = c('MEASYEAR','Slope','Aspect','Conifer_area','Elevation')

# Convert elevation to meters
newvars$Elevation = newvars$Elevation*0.3048

# Extract Elevation data from DEM for plots missing data
library('raster')

dem = raster('Z:/GIS/DEM/alt_2-5m_bil/alt.bil')

# Only extract data for plots already missing data- can do this b/c master and newvars are in same order
newvars = newvars[master$plot.id,] # Puts matches in same order as master data
useplots = subset(master, !is.na(LAT)&is.na(newvars$Elevation))

# Extract data
newelev = extract(dem, useplots[,c('LON','LAT')])

# Put into data frame
rownames(newvars) = newvars$yrplot.id
newvars[useplots$yrplot.id, 'Elevation'] = newelev

# Save data
write.csv(newvars[,c('plot.id','ID','MEASYEAR','Slope','Aspect','Elevation','Conifer_area','radiation')], './Data/fia_lichen_env_data_plots.csv', row.names=F)

##################################################
### Calculate Solar Radiation

library('insol')
options(digits=12)


# Set times over which to integrate
myJD = JD(seq(ISOdate(2010,1,1), ISOdate(2010,12,31), by='30 min'))

# Set air visibility in km and relative humidity
vis = 50
rh = 20

# Set up vector in which to save radiation values
radiation = rep(NA, nrow(master))

for(i in 1:nrow(master)){
	useplot = master[i,]

	lat = useplot$LAT
	lon = useplot$LON
	slope = useplot$Slope
	aspect = useplot$Aspect
	elev = useplot$Elevation

	# If Slope is NA, use Slope = 0
	if(is.na(slope)){
		slope=0
		aspect=0
	}
	if(is.na(aspect)){
		aspect=0
	}

	# Calculate sun positions over which to integrate	
	mysun = sunpos(sunvector(jd=myJD, latitude=lat, longitude=lon, timezone=0))

	# Calculate insolation on flat surface
	myinsol = data.frame(insolation(zenith = mysun[,'zenith'], jd = myJD, height = elev, visibility=vis, RH=rh, tempK=278.15, O3=.03, alphag=0))

	# Calculate vector normal to slope 
	surfacevec = array(normalvector(slope, aspect), dim=c(1,1,3))

	# Calculate the illumination of the slope for each time point (0 indicates self-shading)
	sunvec = sunvector(jd=myJD, latitude=lat, longitude=lon, timezone=0)
	shade = apply(sunvec, 1, function(x)hillshading(surfacevec, x)) 

	# Calculate daily direct incident radiation
	myrad = shade*myinsol$In
	myradtot = sum(myrad)

	# Save total radiation
	radiation[i] = myradtot
}

newvars$radiation = radiation

# Save data
write.csv(newvars[,c('plot.id','ID','MEASYEAR','Slope','Aspect','Elevation','Conifer_area','radiation')], './Data/fia_lichen_env_data_plots.csv', row.names=F)

plot(master$LON, master$LAT, cex=.7, 
	col=colorRampPalette(c('dark blue','blue','yellow','orange','red'))(20)[cut(newvars$radiation, 20)]
)










