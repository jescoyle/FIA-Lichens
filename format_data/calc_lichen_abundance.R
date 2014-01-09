# This script calculates lichen abundance on plots from abundance codes in site x species matrix



################################################################
### Plot abundance vs. richness

siteXsp = read.csv('./Data/siteXsp_matrix_current+legacy.csv', row.names=1, check.names=F)

# Abundance codes:
# 1 = 1-3 ind
# 2 = 4-10 ind
# 3 = > 10 ind, but < 50% of boles
# 4 = > 50% boles

numtrees = rowSums(master[rownames(siteXsp),c('numTreesSm','numTreesBig')])
sum(!(rownames(siteXsp) %in% rownames(master))) #93 plots won;t be counted b/c not included in master table, but that is ok b/c all are in model_data

## Calculate abundance based on occupancy on trees
# 1 = min(2/num trees, .5)
# 2 = min(7/num trees, .5)
# 3 = ifelse(num trees > 20, mean(10/num trees, .5), .5)
# 4 = .75

# New calculation based on occupancy of trees
tot_abun = c()
for(i in rownames(siteXsp)){
	ntrees = numtrees[i]
	newabun = c(min(2/ntrees, .3), min(7/ntrees, .4), ifelse(x > 20, 5/x + .25, .5), .75)

	x = siteXsp[i,]
	spabun = newabun[x[x!=0]]
	tot_abun = c(tot_abun, sum(spabun)*ntrees)
}
names(tot_abun) = rownames(siteXsp)


# Old calculation based on exponentiating abundance classes
tot_abun = apply(siteXsp, 1, function(x){
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
	newabun = c(min(2/ntrees, .3), min(7/ntrees, .4), ifelse(x > 20, 5/x + .25, .5), .75)

	x = siteXsp[i,as.character(parmsp)]
	spabun = newabun[x[x!=0]]
	parm_abun = c(parm_abun, sum(spabun)*ntrees)
}
names(parm_abun) = rownames(siteXsp)

phys_abun = c()
for(i in rownames(siteXsp)){
	ntrees = numtrees[i]
	newabun = c(min(2/ntrees, .3), min(7/ntrees, .4), ifelse(x > 20, 5/x + .25, .5), .75)

	x = siteXsp[i,as.character(physsp)]
	spabun = newabun[x[x!=0]]
	phys_abun = c(phys_abun, sum(spabun)*ntrees)
}
names(phys_abun) = rownames(siteXsp)


# Add abundances to data frame
model_data$tot_abun = tot_abun[rownames(model_data)]
model_data$parm_abun = parm_abun[rownames(model_data)]
model_data$phys_abun = phys_abun[rownames(model_data)]
#model_data$avg_abun = avg_abun[rownames(model_data)]

# Log transform for use in models
model_data$tot_abun_log = log(model_data$tot_abun + 1)
model_data$parm_abun_log = log(model_data$parm_abun + 1)
model_data$phys_abun_log = log(model_data$phys_abun + 1)

# Write out abundance data
abun_data = model_data[,c('tot_abun','parm_abun','phys_abun','tot_abun_log','parm_abun_log','phys_abun_log')]
write.csv(abun_data, './Data/lichen abundance on plots.csv', row.names=T)


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

