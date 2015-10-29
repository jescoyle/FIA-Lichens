# This script explores patterns of trait variation geographically among different lichen families.
# Analyses are performed both at local (FIA data) and regional (CNALH data) scales


source('C://Users/jrcoyle/Documents/UNC/Projects/FIA Lichen/GitHub/FIA-Lichens/load_data.R')
source('./GitHub/FIA-Lichens/fia_lichen_analysis_functions.R')

library(stringr) # str_trim white space remover

###########################################################
### Data

# File describing infomation on traits: categorical vs numeric
traitnames = read.csv('lias_trait_types.csv', row.names=1)

# File with species abundances at each FIA plot
siteXsp = read.csv('./Data/siteXsp_matrix_current+legacy_cleaned.csv', row.names=1, check.names=F)

# File describing species listed in FIA data set
lichensp_raw = read.csv('lichen_species.csv')
lichensp = subset(lichensp_raw, LICH_SPPCD %in% as.numeric(colnames(siteXsp))) 
lichensp = subset(lichensp, is.na(YEAREND))
lichensp_nogen = subset(lichensp, SPECIES!='')

# Include genera not represented by any species
genera = subset(lichensp, SPECIES=='')
gen_nosp = genera[which(!genera$GENUS %in% lichensp_nogen$GENUS),]
lichensp_nogen = rbind(lichensp_nogen, gen_nosp)

# Num. species per family names in FIA data
fam_count = table(lichensp_nogen$family)
fams = rownames(fam_count)

# List of species in each family
fam_list = lapply(fams, function(x){
	use_sp = subset(lichensp_nogen, family==x)
	str_trim(paste(use_sp$GENUS, use_sp$SPECIES))
})
names(fam_list) = fams

# File with trait data for FIA species
traits = read.csv('./Data/LichenTraits/FIA_lichen_species_traits_LIAS.csv')


## For mapping:
# Define map colors
ncuts=8
mycol = read.csv('C:/Users/jrcoyle/Documents/UNC/Projects/blue2red_10colramp.txt')
mycol = apply(mycol,1,function(x) rgb(x[1],x[2],x[3],maxColorValue=256))
mycol = mycol[10:1]
mycolramp = colorRampPalette(mycol)(ncuts)

# Map projection - Lambert Equal Area
plot_prj = paste("+proj=laea +lat_0=40 +lon_0=-97 +units=km",sep='')

# Read in North America outline and re-project
OUTLINES = readOGR('../../GIS shape files/N Am Outline','na_base_Lambert_Azimuthal')
OUTLINES.laea = spTransform(OUTLINES,CRS(plot_prj))

##########################################################
### Trait Summaries Across Families

use_fams = fams[fam_count>3]

# Rename siteXsp columns with species names and
# First: Check which species are not included in trait distance matrix
rownames(traits) = as.character(traits$LICH_SPPCD)
missingsp = colnames(siteXsp)[!(colnames(siteXsp) %in% traits$LICH_SPPCD)] #none missing, safe to rename
colnames(siteXsp) = str_trim(paste(traits[colnames(siteXsp),'Genus'], traits[colnames(siteXsp),'Species']))

# Calculate midpoint of numeric traits
traits$ascospore_length_mid = apply(traits[,c('ascospore_length_min','ascospore_length_max')], 1, mean)
traits$ascospore_width_mid = apply(traits[,c('ascospore_width_min','ascospore_width_max')], 1, mean)
traits$conidia_length_mid = apply(traits[,c('conidia_length_min','conidia_length_max')], 1, mean)

# Make species X trait matrix
spXtrait = traits[,rownames(traitnames)]
rownames(spXtrait) = str_trim(paste(traits$Genus, traits$Species))


## Trait summary boxplots by family
bintraits = rownames(subset(traitnames, type=='binary'))

bintrait_props = sapply(bintraits, function(tr){
	sapply(use_fams, function(f){
		use_trvec = na.omit(spXtrait[fam_list[[f]],tr])
		sum(use_trvec==1)/length(use_trvec)
	})
})

pdf('./Figures/trait prevalence in families.pdf', height=6, width=6)
par(mar=c(8,3,2,1))
for(tr in colnames(bintrait_props)){
	bp = barplot(bintrait_props[,tr], las=2, ylim=c(0,1), main=traitnames[tr, 'displayName'], col=1)
	text(bp, bintrait_props[use_fams, tr]+.05, paste(fam_count[use_fams],' (',missing_traits_byfam[tr,use_fams],')', sep=''), cex=0.6)
}
dev.off()

## Missing data by family
missing_traits_byfam = sapply(fams, function(x){
	use_tr = spXtrait[fam_list[[x]],]
	apply(use_tr, 2, function(y) sum(is.na(y)))
})

rowSums(missing_traits_byfam)

