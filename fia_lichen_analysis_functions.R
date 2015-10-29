## This file contains functions used to analyze FIA lichen data.


######################################################################################################
### Functions


# Two functions that create unique plot identifiers by pasting together state, county plot, and year information.
make.plotid = function(x){
	paste(x$STATECD, x$COUNTYCD, x$PLOT, sep='_')
}
make.yrplotid = function(x){
	paste(x$INVYR, make.plotid(x), sep='_')
}


# Two functions that calculate Spearman rank correlation coefficients between two sets of variables
# and returns a dataframe with correlations or an array of permutation-based confidence intervals
calc.corrs = function(predictors, responses, data){

	richCors = sapply(predictors, function(x){
		noNA = data[!is.na(data[,x]),]
		Xvar = noNA[,x]
		cor(Xvar, noNA[,responses], method="spearman")
	})

	richCors = data.frame(t(richCors))
	names(richCors) = responses

	richCors$predictor = rownames(richCors)
	richCors$nplots = sapply(rownames(richCors), function(x) sum(!is.na(data[,x])))

	return(richCors)
}

calc.corrs.null = function(predictors, responses, data, ITERATIONS=1000, intervals = c(0.025, 0.975)){
	# Calculate null distribution quantiles of rho based on permutation
	richCors.null = array(NA, dim=c(ITERATIONS, length(responses), length(predictors)),
		dimnames=list(iteration=1:ITERATIONS, response=responses, predictor=predictors))
	for(i in 1:ITERATIONS){
		sapply(predictors, function(x){
			noNA = data[!is.na(data[,x]),]
			Xvar = noNA[,x]

			Xvar = sample(Xvar)
		
			cor(Xvar, noNA[,responses], method="spearman")
		}) -> richCors.null[i,,]
	}
	quant = apply(richCors.null, c(2,3), function(x) quantile(x, probs=intervals))

	return(quant)
}


## A function that partitions variation between two models
## Must be supplied with R2 for [A, B, AB]
partvar2 = function(R2s){
	
	a = R2s[3]-R2s[2]
	c = R2s[3]-R2s[1]
	b = R2s[1] - a
	d = 1-R2s[3]

	this_partition = c(a,b,c,d)
	names(this_partition) = c(names(R2s)[1], 'Both', names(R2s)[2], 'Unexplained')

	this_partition
}

## A function that partitions variation between three models
## Must be supplied with R2 for [A, B, C, AB, AC, BC, ABC]
partvar3 = function(R2s){

	setmat = rbind(c(1,0,0,1,1,0,1),
		c(0,1,0,1,0,1,1),
		c(0,0,1,0,1,1,1),
		c(1,1,0,1,1,1,1),
		c(1,0,1,1,1,1,1),
		c(0,1,1,1,1,1,1),
		c(1,1,1,1,1,1,1)
	)
	solve(setmat, R2s) #using matrix algebra
}



# A function which plots any type of polynomial given the coefficients.
#  a is a vector of coefficients
polyfunc = function(x,a){
	series = c()
	
	for(p in 1:length(a)){
		series = cbind(series, a[p]*(x^(p-1)))
	}

	rowSums(series)
}



# A function that plots a vertical color ramp on the side of a plot
# cols    : the colors to use
# n       : number of divisions
# barends : location of whole bar c(xleft, ybottom, xright, ytop)
# labels    : vector of labels for bar, assumes 1st and last numbers correspond to 1st and last colors
# title   : title to print above bar
# mycex   : size of label and title text
# uneven.lab : TRUE when numeric labels are not evenly spaced along length of color bar
# labrange : use when uneven.lab=T to specify minimum and maximum values of color range c(min, max)
# ndig    : number of digits beyond decimal placeto print for numeric labels 
plotColorRamp = function(cols, n, barends, labels=NA, title=NA, mycex=1.5, uneven.lab = F, labrange = NA, ndig=1){
	dX = barends[3] - barends[1]
	dY = barends[4] - barends[2]
	dy = dY/n
	
	xpd.old = par('xpd')
	par(xpd=T)

	lend.old = par('lend')
	par(lend=1)

	usecols = colorRampPalette(cols)(n)

	for(i in 1:n){
		rect(barends[1], barends[2]+dy*(i-1), barends[3], barends[2]+dy*i, col=usecols[i], border=usecols[i])
	}

	if(!is.na(labels)){
		if(is.numeric(labels)){
			labels.round = format(round(labels, ndig), nsmall=ndig, trim=F)

			if(uneven.lab){
				dz = dY/diff(labrange)
				Yposition = barends[2] + dz*(labels-labrange[1])
			} else {
				dZ = labels[length(labels)]-labels[1]
				dz = dY/dZ
				Yposition = barends[2] + dz*(labels-labels[1])
			}

			text(barends[3]+dX*0.5, Yposition, labels.round, pos=4, cex=mycex)

		} else {
			labels.round = labels
			dz = dY/length(labels)
			Yposition = barends[2] + dz*(0:(length(labels)-1))

			text(barends[3]+dX*0.5, Yposition, labels, pos=4, cex=mycex)
		}

		segments(barends[3], Yposition, barends[3]+dX*0.5, Yposition)	
	}
	if(!is.na(title)){
		
		
		## Determine how many characters away to place title
		digits = max(nchar(labels.round)) # Maximum number of digits in a label
		largest = labels.round[which(nchar(labels.round)==digits)] # Which labels are longest
		
		small.chars = grep('[-.]', largest) # Does the largest label have a small character?
			if(length(small.chars)==length(largest)) digits = digits-0.6 # Discount the size of the largest label by 0.6 a character
		
		text(barends[3]+dX*0.5+par('cxy')[1]*mycex*(digits+.5), barends[2]+0.5*dY, labels=title, srt=-90, cex=mycex)
	}
	par(xpd=xpd.old)
	par(lend=lend.old)
}


### The following functions were used in the revised version of the manuscript:


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

# A function that formats a coefficient table for publication from a model average object x
format_coefTable = function(x){

	# Get coefficient table
	tab = make_coefTable(x)

	# Add variable importance and number of models
	imp = importance(x)
	imp_df = as.data.frame(imp)
	imp_df$nmods = attr(imp, 'n.models')
	tab = cbind(tab, imp_df[rownames(tab),])
	
	# Split coef table into linear and quadratic terms
	sq = F
	sqterms = tab$var[grep('2', tab$var)]		
	if(length(sqterms)>0){
		sq = T
		tab2 = tab[sqterms,]
		tab2$var = substr(tab2$var, 1, nchar(tab2$var)-1)
		rownames(tab2) = tab2$var
		tab2$Term = 'quadratic'
	}

	# Subset coef table to only linear effects of predictors
	tab = subset(tab, var %in% rownames(predtypes))
	tab$Term = 'linear'
	ordered_vars = tab$var[order(abs(tab$est), decreasing=T)]

	# Combine tables
	newtab = rbind(tab[ordered_vars,], tab2)
	
	# Add variable labels
	newtab$Predictor = varnames[newtab$var,'midName']

	# Add variable mode and scale
	newtab$Scale = toupper(substr(predtypes[newtab$var,'scale'],1,1))
	newtab$Mode = toupper(substr(predtypes[newtab$var,'mode'],1,1))

	# Add variable importance and number of models
	newtab$Importance = format(newtab$imp, digits=2)
	newtab$'Num. Models' = paste(newtab$nmods, nrow(x$msTable), sep='/')

	# Format estimate and SE
	newtab$Estimate = format(newtab$est, digits=2)
	newtab$SE = format(newtab$se, digits=2)

	newtab[,c('Predictor','Term','Scale','Mode','Importance','Num. Models','Estimate','SE')]
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

