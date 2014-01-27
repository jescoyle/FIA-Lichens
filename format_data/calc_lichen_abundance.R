# This script calculates lichen abundance on plots from abundance codes in site x species matrix

setwd('./UNC/Projects/FIA Lichen')
source('./GitHub/FIA-Lichens/load_data.R')

################################################################
### Plot abundance vs. richness

siteXsp = read.csv('./Data/siteXsp_matrix_current+legacy.csv', row.names=1, check.names=F)

# Abundance codes:
# 1 = 1-3 ind
# 2 = 4-10 ind
# 3 = > 10 ind, but < 50% of boles
# 4 = > 50% boles

numtrees = rowSums(master[rownames(siteXsp),c('numTreesSm','numTreesBig')])
sum(!(rownames(siteXsp) %in% rownames(master))) #93 plots won't be counted b/c not included in master table, but that is ok b/c all are in model_data

## Calculate abundance based on occupancy on trees
# 1 = min(2/num trees, .3)
# 2 = min(7/num trees, .4)
# 3 = ifelse(num trees > 20, mean(10/num trees, .5), .5)
# 4 = .75

# New calculation based on occupancy of trees
tot_abun = c()
for(i in rownames(siteXsp)){
	ntrees = numtrees[i]
	newabun = c(min(2/ntrees, .3), min(7/ntrees, .4), ifelse(ntrees > 20, 5/ntrees + .25, .5), .75)

	x = siteXsp[i,]
	spabun = newabun[x[x!=0]]
	tot_abun = c(tot_abun, sum(spabun)*ntrees)
}
names(tot_abun) = rownames(siteXsp)


# NOTE: 1/24/2014 Actually this is the abundance that has been used in all of the analyses.
# Old calculation based on exponentiating abundance classes
tot_abun_exp = apply(siteXsp, 1, function(x){
	sum(exp(x[x!=0]))
})

avg_abun = apply(siteXsp, 1, function(x){
	mean(exp(x[x!=0]))
})

## Calculate abundance of lichen families
splist = read.csv('lichen_species.csv') # List of lichen species that identifies the family they are in
sum(!(colnames(siteXsp) %in% splist$LICH_SPPCD)) # All species accounted for

parmsp = subset(splist, family=='Parmeliaceae')$LICH_SPPCD
physsp = subset(splist, family=='Physciaceae')$LICH_SPPCD

parmsp = parmsp[parmsp %in% colnames(siteXsp)]
physsp = physsp[physsp %in% colnames(siteXsp)]

parm_abun = c()
for(i in rownames(siteXsp)){
	ntrees = numtrees[i]
	newabun = c(min(2/ntrees, .3), min(7/ntrees, .4), ifelse(ntrees > 20, 5/ntrees + .25, .5), .75)

	x = siteXsp[i,as.character(parmsp)]
	spabun = newabun[x[x!=0]]
	parm_abun = c(parm_abun, sum(spabun)*ntrees)
}
names(parm_abun) = rownames(siteXsp)

phys_abun = c()
for(i in rownames(siteXsp)){
	ntrees = numtrees[i]
	newabun = c(min(2/ntrees, .3), min(7/ntrees, .4), ifelse(ntrees > 20, 5/ntrees + .25, .5), .75)

	x = siteXsp[i,as.character(physsp)]
	spabun = newabun[x[x!=0]]
	phys_abun = c(phys_abun, sum(spabun)*ntrees)
}
names(phys_abun) = rownames(siteXsp)


# Add abundances to data frame
newabun = data.frame(tot_abun,parm_abun, phys_abun)

# Log transform for use in models
newabun$tot_abun_log = log(newabun$tot_abun + 1)
newabun$parm_abun_log = log(newabun$parm_abun + 1)
newabun$phys_abun_log = log(newabun$phys_abun + 1)
newabun$yrplot.id = rownames(newabun)

# Write out abundance data
write.csv(newabun, './Data/lichen abundance based on tree occupancy.csv', row.names=F)


### Plot the tree-occupance -> abundance function ###
func1 = function(N) min(2/N, .3)
func2 = function(N) min(7/N, .4)
func3 = function(N) ifelse(N > 20, 5/N + .25, .5)
func4 = function(N) .75

ntrees = 2:100

pdf('./Paper/supp fig abundance calc.pdf', height=4, width=7)
par(mfrow=c(1,2))
par(mar=c(5,5,1,1))
plot(ntrees, sapply(ntrees, func4), pch=16, las=1, xlab='Number of Trees', 
	ylab='Proportion of Trees Occupied', ylim=c(0,1))
points(ntrees, sapply(ntrees, func3), pch=16, col='dodgerblue')
points(ntrees, sapply(ntrees, func2), pch=16, col='forestgreen')
points(ntrees, sapply(ntrees, func1), pch=16, col='darkorange')
text(1, c(func1(2), func2(2), func3(2), func4(2)), paste('Class', 1:4), adj=c(0,-.5))

plot(ntrees, sapply(ntrees, func4)*ntrees, pch=16, las=1, xlab='Number of Trees', 
	ylab='Number of Trees Occupied', ylim=c(0,100))
points(ntrees, sapply(ntrees, func3)*ntrees, pch=16, col='dodgerblue')
points(ntrees, sapply(ntrees, func2)*ntrees, pch=16, col='forestgreen')
points(ntrees, sapply(ntrees, func1)*ntrees, pch=16, col='darkorange')
text(100, c(func1(100), func2(100), func3(100), func4(100))*100, paste('Class', 1:4), adj=c(1,-.5))
dev.off()


###########################################################
### Expected correlation between richness and abundance

library(picante)
library(MASS)
library(MuMIn)

## Based on summing the exponentiation of the abundance class.

freq = table(siteXsp[siteXsp!=0])
freq = freq/sum(freq)

richness = rowSums(siteXsp>0)
abun_mat = matrix(NA, nrow=nrow(siteXsp), ncol=1000)
rownames(abun_mat) = rownames(siteXsp)

for(i in 1:1000){

# Sample abundance categories according to their global frequency
# then calculate abundance
sapply(richness, function(n){
	abuns = sample(1:4, n, replace=T, prob=freq)
	sum(exp(abuns))
}) -> randAbun

# Save in table
abun_mat[,i] = randAbun
}

# Only look at plots in our data set and transform them as they were transformed in our dataset
use_abun = log(abun_mat[testplots$yrplot.id,]+1)]

trans_data_test = trans_data[testplots$yrplot.id,]

# Calculate correlation between abundance and richness
apply(use_abun, 2, function(x){
	usemod = glm(trans_data_test$lichen.rich~x, family=poisson(link='log'))
	attr(r.squaredLR(usemod, null=glm(lichen.rich~1, family=poisson(link='log'), data=trans_data_test)), 'adj.r.squared')
}) -> R2s
obs_mod = glm(lichen.rich~tot_abun_log, data=trans_data_test, family=poisson(link='log'))
obs.R2 = r.squaredLR(obs_mod, null = glm(lichen.rich~1, data=trans_data_test, family=poisson(link='log')))

apply(use_abun, 2, function(x){
	cor(trans_data_test$lichen.rich,x)
}) -> r2s
obs.r2 = cor(trans_data_test$lichen.rich,trans_data_test$tot_abun_log)

hist(R2s, xlim=c(.95,1))
abline(v=obs.R2)

hist(r2s, xlim=c(.8,1))
abline(v=obs.r2)

# What about when abundances are drawn randomly?
abun_mat2 = abun_mat
for(i in 1:1000){

sapply(richness, function(n){
	abuns = sample(1:4, n, replace=T)
	sum(exp(abuns))
}) -> randAbun

# Save in table
abun_mat[,i] = randAbun
}

# I DON'T UNDERSTAND WHY THE CORRELATION IS SO STRONG AND THE OBSERVED IS SO WEAK

# abundance is calculated in the same way!
cbind(trans_data[testplots$yrplot.id,'tot_abun_log'], log(tot_abun_exp[testplots$yrplot.id]+1))


### Based on tree occupancy


# New calculation based on occupancy of trees
tot_abun_occu = c()
for(i in rownames(siteXsp)){
	ntrees = numtrees[i]
	newabun = c(min(2/ntrees, .3), min(7/ntrees, .4), ifelse(ntrees > 20, 5/ntrees + .25, .5), .75)

	x = siteXsp[i,]
	spabun = newabun[x[x!=0]]
	tot_abun_occu = c(tot_abun_occu, sum(spabun)*ntrees)
}
names(tot_abun_occu) = rownames(siteXsbp)

# Correlation between richness and number of trees
plot(richness~numtrees)
abline(lm(richness~numtrees))
cor(numtrees,richness, use='complete.obs')


# Relationship between abundance and richness
plot(log(richness)~log(tot_abun_occu))
hist(log(tot_abun_occu))

# Decide on model
fit_rich = richness[fitplots$yrplot.id] 
fit_abun = log(tot_abun_occu[fitplots$yrplot.id]+1)

pois_log = glm(fit_rich~fit_abun, family=poisson(link='log'))
#pois_iden = glm(fit_rich~fit_abun, family=poisson(link='identity'), start=rep(1,2))
gaus_log = glm(fit_rich~fit_abun, family=gaussian(link='log'))
gaus_iden = glm(fit_rich~fit_abun, family=gaussian(link='identity'))
nb_log = glm.nb(fit_rich~fit_abun, link='log')
#nb_iden = glm.nb(fit_rich~fit_abun, link='identity',init.theta=1)

AIC(pois_log, gaus_log, gaus_iden, nb_log)

use_x = seq(min(fit_abun), max(fit_abun), length.out=100)

pois_y = predict(pois_log, list(fit_abun=use_x), type='response')
gauslog_y = predict(gaus_log, list(fit_abun=use_x), type='response')
gausiden_y = predict(gaus_iden, list(fit_abun=use_x), type='response')
nblog_y = predict(nb_log, list(fit_abun=use_x), type='response')

plot(log(fit_rich)~fit_abun)
lines(use_x, log(pois_y), col='red', lwd=3)
lines(use_x, log(gauslog_y), col='blue', lwd=3)
lines(use_x, log(gausiden_y), col='green', lwd=3)
lines(use_x, log(nblog_y), col='orange', lwd=3)

# Make randomized abundance matrices

abun_mat_occu = matrix(NA, nrow=nrow(siteXsp), ncol=1000)
rownames(abun_mat_occu) = rownames(siteXsp)
for(j in 1:1000){

# Sample abundance categories according to their global frequency
# then calculate abundance
sapply(names(richness), function(i){
	n = richness[i]
	x = sample(1:4, n, replace=T, prob=freq)

	ntrees = numtrees[i]
	newabun = c(min(2/ntrees, .3), min(7/ntrees, .4), ifelse(ntrees > 20, 5/ntrees + .25, .5), .75)

	spabun = newabun[x]
	sum(spabun)*ntrees

}) -> randAbun

# Save in table
abun_mat_occu[,j] = randAbun
}

## Calculate correlation for randomized abundance matrix
use_abun_mat = log(abun_mat_occu2[fitplots$yrplot.id,]+1)

apply(use_abun_mat, 2, function(x){
	usemod = glm(fit_rich~x, family=poisson(link='log'))
	attr(r.squaredLR(usemod, null=glm(fit_rich~1, family=poisson(link='log'))), 'adj.r.squared')
}) -> R2s
obs_mod = glm(fit_rich~fit_abun, family=poisson(link='log'))
obs.R2 = r.squaredLR(obs_mod, null = glm(lichen.rich~1, data=trans_data_test, family=poisson(link='log')))

hist(R2s, xlim=c(.92,1))
abline(v=obs.R2, lwd=2, col='red')

# correlation
apply(use_abun_mat, 2, function(x){
	cor(log(fit_rich+1),x)
}) -> cors
obs.cor= cor(log(fit_rich+1),fit_abun)

hist(cors, xlim=c(.75,.9))
abline(v=obs.cor, lwd=2, col=2)

# Make randomized abundance matrices when abundance classes are drawn from uniform distribution

abun_mat_occu2 = matrix(NA, nrow=nrow(siteXsp), ncol=1000)
rownames(abun_mat_occu2) = rownames(siteXsp)
for(j in 1:1000){

# Sample abundance categories according to their global frequency
# then calculate abundance
sapply(names(richness), function(i){
	n = richness[i]
	x = sample(1:4, n, replace=T)

	ntrees = numtrees[i]
	newabun = c(min(2/ntrees, .3), min(7/ntrees, .4), ifelse(ntrees > 20, 5/ntrees + .25, .5), .75)

	spabun = newabun[x]
	sum(spabun)*ntrees

}) -> randAbun

# Save in table
abun_mat_occu2[,j] = randAbun
}
























