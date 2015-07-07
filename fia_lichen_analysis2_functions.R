## This script contains functions used in the analyses for FIA Lichen manuscript version 2




## Functions that create neighbor, weights and listw objects given a spatial dataset

# Create a listw object with weights based on overlapping areas of 500km circles
# The weight is the percentage of the area shared with overlapping circle from another observation
A_intersect = function(d, R){ 2*(R^2)*acos(d/(2*R)) - 0.5*d*sqrt(4*(R^2) - (d^2))}
area_dist_func = function(d) A_intersect(d, 500)/(pi*(500^2))
make_regnb = function(sp_data) dnearneigh(sp_data, 0, 1000, row.names=rownames(coordinates(sp_data)))
make_regweights = function(sp_data, reg_nb) lapply(nbdists(reg_nb, sp_data), area_dist_func) 
make_reglistw = function(sp_data, reg_nb) nb2listw(reg_nb, glist=make_regweights(sp_data, reg_nb), style='W')


# A functions for calculating cumulative weight of models from a model dredging table
calc_cumWeight = function(x) sapply(1:length(x), function(i) sum(x[1:i]))



## Functions for creating coefficient plots
library(lattice)

# A function that makes the coefficient table from a model averaging object
make_coefTable = function(x){

	var = names(coef(x))
	tab = as.data.frame(coefTable(x))
	tab$var = rownames(tab)
	names(tab)[1:2] = c('est','se')
	tab = tab[,c(3,1,2)]
	tab$ci.lower = qnorm(0.025)*tab$se + tab$est
	tab$ci.upper = qnorm(0.975)*tab$se + tab$est

	tab	
}	

plot_coefs = function(tab, shadeby){

	myrange = range(tab[,c('ci.upper','ci.lower')]) + c(-0.05, 0.05)
	
	# Split coef table into linear and quadratic terms
	sq = F
	sqterms = tab$var[grep('2', tab$var)]		
	if(length(sqterms)>0){
		sq = T
		tab2 = tab[sqterms,]
		tab2$var = substr(tab2$var, 1, nchar(tab2$var)-1)
		rownames(tab2) = tab2$var 
	}

	# Subset coef table to only linear effects
	tab = subset(tab, var %in% rownames(predtypes))
	ordered_vars = tab$var[order(tab$est)]
	
	
	# Color scheme
	#mycols = matrix(c('#b3b5ffff','#6b6dd7ff','#8dff94ff','#38af4fff'), nrow=2)
	#colnames(mycols) = c('regional','local')
	#rownames(mycols) = c('het','opt')
	mycolsbw = c('grey80','white')
	names(mycolsbw) = unique(predtypes[tab$var,shadeby])

	# Offset for squared terms
	nudge = 0.06
	
	dotplot(as.numeric(factor(var, levels = ordered_vars))~est, data=tab, 
		xlab=list('Standardized Effect',cex=3), ylab='',
		main='',cex.lab=3,aspect=5/3, xlim=myrange,
		panel=function(x,y){
	
		# Add horizontal boxes
		colorcombos = predtypes[ordered_vars, c('mode','scale')]
		if('regS' %in% tab$var) colorcombos['regS','mode'] = 'opt'
		#colororder = apply(colorcombos, 1, function(x) mycols[x[1],x[2]])
		colororder = mycolsbw[colorcombos[ordered_vars,shadeby]]

		panel.rect(myrange[1]-0.01,1:length(ordered_vars)-.5, myrange[2]+0.01, 1:length(ordered_vars)+.5,
			col=colororder, border='grey50')
	
		# Add vertical line at 0
		panel.abline(v=c(-.2,0,.2), col='grey30', lty=2, lwd=2)		
		
		# Add null distribution segments for total effects
		panel.segments(tab[ordered_vars,'ci.lower'], y,
			tab[ordered_vars,'ci.upper'], y, 
			col='black', lwd=4.5, lend=1)
		# Add points for total estimated effects
		panel.points(x, y, col='black', fill='white', pch=21, cex=3, lwd=3) 
	
		# Add points for squared terms if they exist
		if(sq){
			panel.segments(tab2[ordered_vars,'ci.lower'], y+nudge,
				tab2[ordered_vars,'ci.upper'], y+nudge, 
				col='black', lwd=4.5, lend=1)
			panel.points(tab2$est, nudge + as.numeric(factor(tab2$var, levels=ordered_vars)), col='black', fill='grey20', pch=21, cex=2, lwd=3)
		}
		},
		scales=list(y=list(labels=varnames[ordered_vars,'midName'], 
			cex=3, col='black'),
			x=list(cex=3, tick.number=8))
	)
}

# A function that formats a table of variable importances from a model averaging object
make_impTable = function(x, Nmod){
	imp = importance(x)
	tab = data.frame(imp)
	tab$nmods = paste(attr(imp, 'n.models'), Nmod, sep='/')
	tab$var = rownames(tab)
	
	sqterms = grep(2, tab$var)
	if(length(sqterms)>0) tab$var[sqterms] =  substr(tab$var[sqterms] , 1, nchar(tab$var[sqterms])-1)

	tab$term = 'linear'
	tab$term[sqterms] = 'quadratic'

	tab$scale = predtypes[tab$var,'scale']
	tab$mode = predtypes[tab$var, 'mode']
	
	tab$displayName = varnames[tab$var, 'midName'] 

	tab = tab[,c('var','displayName','term','imp','nmods','scale','mode')]
}





