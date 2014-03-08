# This script calculated regional tree species richness from FIA plots


rawtrees = read.csv('fia_all_trees_1997-2008.csv')

# Make unique plot identifier

rawtrees$yrplot.id = paste(rawtrees$INVYR, rawtrees$STATECD, rawtrees$COUNTYCD, rawtrees$PLOT, sep='_')


# Read in Plot locations data
statelist = c('AL','AK','AZ','AR','CA','CO','CT','DE','FL','GA','ID','IL','IN','IA','KS','KY','LA','ME','MD','MA','MI','MN','MS','MO','MT','NE','NV','NH','NJ','NM','NY','NC','ND','OH','OK','OR','PA','RI','SC','SD','TN','TX','UT','VT','VA','WA','WV','WI','WY')
use_cols = c('CN','PREV_PLT_CN','INVYR','STATECD','COUNTYCD','PLOT','PLOT_STATUS_CD','DESIGNCD','LAT','LON','ELEV','MACRO_BREAKPOINT_DIA','QA_STATUS')

plots = matrix(ncol=13, dimnames=list(NA,use_cols))
for(s in statelist){
	these_plots = read.csv(paste('./PLOT/',s,'_PLOT.CSV', sep=''))

	these_plots = these_plots[,use_cols]	
	
	plots = rbind(plots, these_plots)

	print(s)
}

plots = plots[-1,]

write.csv(plots, 'fia_plots.csv', row.names=F)

# Make unique plot identifier

plots$yrplot.id = paste(plots$INVYR, plots$STATECD, plots$COUNTYCD, plots$PLOT, sep='_')

# Subset plots to list of those in tree data set

plotsWtrees = unique(rawtrees$yrplot.id)
plots = subset(plots, yrplot.id %in% plotsWtrees)


# Make site X sp matrix and order it to be the same as the plots

siteXsp = table(rawtrees$yrplot.id, rawtrees$SPCD)
siteXsp = siteXsp[plots$yrplot.id,]

write.csv(siteXsp, 'fia_tree_siteXsp.csv', row.names=T)

dim(siteXsp) #199381 plots with 387 'species'


# Read in locations of lichen plots
master = read.csv('fia_lichen_master_data.csv')
lichen_locs = subset(master, state.abbr != 'AK')[,c('yrplot.id','LAT','LON')]
lichen_locs = subset(lichen_locs, !is.na(LAT))

# Convert to spatial data
library(sp)
library(rgdal)

coordinates(lichen_locs) = c('LON','LAT')
proj4string(lichen_locs) = CRS("+proj=longlat")
coordinates(plots) = c('LON','LAT')
proj4string(plots) = CRS("+proj=longlat")

### Calculate minimum number of FIA plots around each lichen plot
sapply(1:nrow(lichen_locs), function(i){
	p = lichen_locs[i,]
	these_plots = find_recs(p, plots, 500)

	nrow(these_plots)
}) -> nplots

# min: 1723, max 54080

# Use 1000 plots

# Remove unidentified trees from the tree list
use_trees = subset(rawtrees, !(SPCD %in% c(999,998,299)))
REPS = 501
NSAMP = 1000

regS_tree = c()
for(i in 1:nrow(lichen_locs)){
	p = lichen_locs[i,]
	these_plots = find_recs(p, plots, 500)
	
	mysamples = sapply(1:REPS, function(x){
		plot_subsamp = these_plots[sample(1:nrow(these_plots), NSAMP, replace=F),]$yrplot.id
		tree_subsamp = subset(use_trees, yrplot.id %in% plot_subsamp)
		
		length(unique(tree_subsamp$SPCD))
	})

	return(median(mysamples))

}





####################################################################
### Functions ###

# A function that finds all records within a given distance buffer from a points
find_recs = function(point, recs, buff){
	dists = spDistsN1(recs, point, longlat=T) # In km
	userecs = recs[dists<buff,]

	return(userecs)
}
