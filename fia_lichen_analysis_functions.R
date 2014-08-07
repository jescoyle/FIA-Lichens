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


