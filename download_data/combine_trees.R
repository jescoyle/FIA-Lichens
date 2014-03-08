


statelist = c('AL','AK','AZ','AR','CA','CO','CT','DE','FL','GA','ID','IL','IN','IA','KS','KY','LA','ME','MD','MA','MI','MN','MS','MO','MT','NE','NV','NH','NJ','NM','NY','NC','ND','OH','OK','OR','PA','RI','SC','SD','TN','TX','UT','VT','VA','WA','WV','WI','WY')

use_cols = c('CN','PLT_CN','PREV_TRE_CN','INVYR','STATECD','COUNTYCD','PLOT','SUBP','TREE','CONDID','AZIMUTH','DIST','STATUSCD','SPCD','SPGRPCD','DIA') 

trees = matrix(ncol=16, dimnames=list(NA,use_cols))

for(s in statelist){
	these_trees = read.csv(paste('./TREE/',s,'_TREE.ZIP/',s,'_TREE.CSV', sep=''))

	these_trees = these_trees[,use_cols]	
	
	these_trees = subset(these_trees, SUBP %in% 1:4)
	these_trees = subset(these_trees, INVYR<=2008&INVYR>=1997)

	trees = rbind(trees, these_trees)

	print(s)
}


trees = trees[-1,]

write.csv(trees, 'fia_all_trees_1997-2008.csv', row.names=F)