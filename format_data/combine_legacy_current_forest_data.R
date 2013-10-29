### This script combines forest structure and lichen diversity measurements from 
### legacy and current data sets.


options(stringsAsFactors=F)

setwd('C:/Users/jrcoyle/Documents/UNC/Projects/FIA Lichen')

years.used= read.csv('match_lichen_tree_current.csv')

# Legacy Data
lichDiv = read.csv('./Data/lichen_richness_legacy.csv')
plot = read.csv('./Data/lichen_plotData_legacy.csv')
treeDiv = read.csv('./Data/TreeData/master_data_tree_legacy.csv')

legacy = merge(plot[,c('yrplot.id','plot.id','YEAR','STATE','COUNTY','ELEVATION')], lichDiv, all.y=T)
legacy = merge(legacy, treeDiv, all.x=T)


# Current Data
state.matcher = read.csv('states.csv')
lichDiv = read.csv('./Data/lichen_richness_current.csv')
lichen.plots = read.csv('lichen_plot_list.csv')
treeDiv = read.csv('./Data/TreeData/master_data_tree_current.csv')

current = merge(years.used, lichen.plots, by.x='lich.yrplot.id', by.y='yrplot.id')
current = merge(current, lichDiv, by.x=c('lich.yrplot.id','plot.id'),by.y=c('yrplot.id','plot.id'))
current = merge(current, treeDiv, by.x='tree.yrplot.id',by.y='yrplot.id', all.x=T)

# 7 plots without tree data
subset(current, is.na(S.tree))$plot.id->missing.plots
subset(trees, plot.id %in% missing.plots)


## Create a dataframe to be used in models (no duplicate plots)
length(unique(legacy$plot.id)); nrow(legacy)
length(unique(current$plot.id)); nrow(current)

current.new = data.frame()
for(p in unique(current$plot.id)){
	these.plots = subset(current, plot.id == p)	
	
	if(nrow(these.plots) > 1){

		# Choose the year with less missing data
		nas = apply(these.plots, 1, function(x) sum(is.na(x)))
		use.plot = these.plots[which(nas==min(nas)),]

		# If multiple years have few nas, use the most recent year
		if(nrow(use.plot) > 1){
			years = use.plot$MEASYEAR
			use.plot = subset(use.plot, MEASYEAR==years[years == max(years)])
		}
		
	} else {
		use.plot = these.plots
	}

	current.new = rbind(current.new, use.plot) 
}

legacy.new = data.frame()
for(p in unique(legacy$plot.id)){
	these.plots = subset(legacy, plot.id == p)	
	
	if(nrow(these.plots) > 1){

		# Choose the year with less missing data
		nas = apply(these.plots, 1, function(x) sum(is.na(x)))
		use.plot = these.plots[which(nas==min(nas)),]

		# If multiple years have few nas, use the most recent year
		if(nrow(use.plot) > 1){
			years = use.plot$YEAR
			use.plot = subset(use.plot, YEAR==years[years == max(years)])
		}
		
	} else {
		use.plot = these.plots
	}

	legacy.new = rbind(legacy.new, use.plot) 
}


# Combine data
unique(c(names(legacy), names(current)))

current.new = current.new[,c(2,3,7,8,10:14,18:79)]
names(current.new)[c(1,3,4,5,8)] = c('yrplot.id','STATE','COUNTY','YEAR','ELEVATION')

legacy.new = data.frame(legacy.new[,c(1,2,4,5,3)], LAT=NA, LON=NA, legacy.new[,c(6:69)])
legacy.new = legacy.new[,names(current.new)]

master.data = rbind(current.new, legacy.new)

write.csv(master.data, 'master_data_forest.csv',row.names=F)


