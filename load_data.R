## This script loads data needed for models of lichen diversity in FIA plots.


setwd('C://Users/jrcoyle/Documents/UNC/Projects/FIA Lichen')
options(stringsAsFactors=F)
varnames=read.csv('varnames.csv', row.names=1)

records_sp = read.csv('./Data/Regional Richness/CNALH_records_fia_genera_NAm.csv')
save.image('FIA_lichen_analysis_final.Rdata')

master = read.csv('./Data/fia_lichen_master_data.csv')
rownames(master) = master$yrplot.id
model_data = read.csv('./Data/fia_lichen_model_data.csv', row.names=1)
trans_data = read.csv('./Data/fia_lichen_trans_data.csv', row.names=1)
working_data = read.csv('./Data/fia_lichen_working_data.csv', row.names=1)

# Read in lists of test and fit plots
testplots = read.csv('./Data/model test plots.csv')
fitplots = read.csv('./Data/model fit plots.csv')

# Subset working data to be in testing or fitting data sets
working_data_test = working_data[testplots$yrplot.id,]
working_data_fit = working_data[fitplots$yrplot.id,]



## Generic Functions ##

# A function that plots a vertical color ramp on the side of a plot
# cols    : the colors to use
# n       : number of divisions
# barends : location of whole bar c(xleft, ybottom, xright, ytop)
# labels    : vector of labels for bar, assumes 1st and last numbers correspond to 1st and last colors
# title   : title to print above bar
plotColorRamp = function(cols, n, barends, labels=NA, title=NA, mycex=1.5, uneven.lab = F, labrange = NA){
	dX = barends[3] - barends[1]
	dY = barends[4] - barends[2]
	dy = dY/n
	
	xpd.old = par('xpd')
	par(xpd=T)

	usecols = colorRampPalette(cols)(n)

	for(i in 1:n){
		rect(barends[1], barends[2]+dy*(i-1), barends[3], barends[2]+dy*i, col=usecols[i], border=NA)
	}

	if(!is.na(labels)){
		if(is.numeric(labels)){

			if(uneven.lab=T){
				dz = dY/diff(labrange)
				Yposition = barends[2] + dz*(labels-labrange[1])
			} else {
				dZ = labels[length(labels)]-labels[1]
				dz = dY/dZ
				Yposition = barends[2] + dz*(labels-labels[1])
			}

			text(barends[3]+dX*0.5, Yposition, round(labels,2), pos=4, cex=mycex)

		} else {
			dz = dY/length(labels)
			Yposition = barends[2] + dz*(0:(length(labels)-1))

			text(barends[3]+dX*0.5, Yposition, labels, pos=4, cex=mycex)
		}

		segments(barends[3], Yposition, barends[3]+dX*0.5, Yposition)	
	}
	if(!is.na(title)){
		labels.round = round(labels, 2)
		
		## Determine how many characters away to place title
		digits = max(nchar(round(labels, 2))) # Maximum number of digits in a label
		largest = labels.round[which(nchar(labels.round)==digits)] # Which labels are longest
		no.decimal = sum(largest == floor(largest))>0 # Does one largest label lack a decimal?
			if(!no.decimal) digits = digits-0.6 # Discount the size of the largest label by 0.6 a character
		no.negative = sum(largest >= 0)>0 # Does one largest label lack a negative sign?
			if(!no.negative) digits = digits-0.6 # Discount the size of the largest label by 0.6 a character
		
		text(barends[3]+dX*0.5+par('cxy')[1]*mycex*(digits+.5), barends[2]+0.5*dY, labels=title, srt=-90, cex=mycex)
	}
	par(xpd=xpd.old)
}
