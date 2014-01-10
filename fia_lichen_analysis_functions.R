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


## A function that partitions variation between three models
## Assumes same set of observations and same underlying model structure
## For some reason this function isn't working. It returns the error that it can't find dataB.

partvar3 = function(modlist){
	#require(MuMIn)
	modA = modlist[[1]]
	modB = modlist[[2]]
	modC = modlist[[3]]

	dataA = as.matrix(modA$model)
	response = dataA[,1]
	dataA = dataA[,2:ncol(dataA)]
	dataB = as.matrix(modB$model)
	dataB = dataB[,2:ncol(dataB)]	
	dataC = as.matrix(modC$model)
	dataC = dataC[,2:ncol(dataC)]

	modAB = update(modA, ~.+dataB)
	modAC = update(modA, ~.+dataC)
	modBC = update(modB, ~.+dataC)
	modABC = update(modAB, ~.+dataC)

	mods = list(modA, modB, modC, modAB, modAC, modBC, modABC)
	r2s = sapply(mods, r.squaredLR)

	components = r2s[7]-r2s[6:4]
	components = c(components, -(r2s[4]-r2s[1]-r2s[2]),
		-(r2s[6]-r2s[2]-r2s[3]), -(r2s[5]-r2s[1]-r2s[3]))
	components = c(components, -(r2s[7]-sum(components[1:6]))/2)
	components[4:6] = components[4:6]-components[7]
	
	names(components) = c('A','B','C','AplusB','BplusC','AplusC','all')
	
	return(components)
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
