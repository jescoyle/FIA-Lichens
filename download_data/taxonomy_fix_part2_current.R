## This script was used to make further taxonomic changes in order to use both
## current and legacy FIA lichen data. Data will range from 1994 - 2009
setwd('C:/Users/jrcoyle/Documents/UNC/Projects/FIA Lichen')

load('lichen_analysis_basic.Rdata')

state.matcher = read.csv("C:/Users/jrcoyle/Documents/UNC/Projects/FIA_Fall2011/Data/FIA_state_codes.csv")
unique(lichenPlot$STATECD)

NorthReg = c("IL","IN",'IA','KS','MI','MN','MO','NE','ND','SD','WI','DE','CT','MA','MD','ME','NH','NJ','NY','OH','PA','RI','VT','WV')
SouthReg = c('AL','AR','FL','GA','KY','LA','MS','NC','OK','SC','TN','TX','VA')
PNWReg = c('AK','CA','HI','OR','WA')
IntWestReg = c('AZ','CO','ID','MT','NM','NV','WY','UT')

weststates = state.matcher[state.matcher$state.name %in% c(PNWReg,IntWestReg), 'STATECD']
eaststates = state.matcher[state.matcher$state.name %in% c(NorthReg,SouthReg), 'STATECD']
southstates = state.matcher[state.matcher$state.name %in% SouthReg, 'STATECD']



### Combine Species

# a comment follows each line of code indicating whether the code was run on the data or how many observations were affected by the change.

data = lichen

### Take care of complicated Action Code 2 (combos)

# Melanohalea multispora.  4010 Melanohalea multispora should be combined into 4017 M. subolivacea. 4010 Melanelia multispora = 4010 Melanohalea multispora. The latter name was used starting with 2005 data. Data from 1994-98 normally include 4010 M. multispora with 4017 M. subolivacea. Starting with 1999, if this species is recorded spores have been observed in thin section and this represents a correct name. 4010 Melanohalea multispora should be combined into 4017 M. subolivacea. 4017 Melanelia subolivacea = 4017 Melanohalea subolivacea. The latter name was used starting with 2005 data. Western US: In data from 1994-8 this taxon normally represents a group including 4017 M. subolivacea and 4010 M. multispora, using the name Melanelia subolivacea group.  Also see notes under 4010. In arid habitats small specimens of 4017 M. subolivacea [group] may be difficult to distinguish from small specimens of 1008 Tuckermanella fendleri. Narrower lobes, pale lower surface, uniform color for the latter are the most reliable characters in absence of elevated pycnidia on margins. Follow directions in “Lichen Specialist Identification Procedures/Isomorphs.”  

data[(data$LICH_SPPCD==4010)&(data$MEASYEAR<1999),'LICH_SPPCD']<-4017 #None

### Take care of Action Code 3 (analyses crossing multiple years)

## Pseudocyphellaria sp. For analyzing data from multiple years crossing 2003, 6408 Pseudosyphellaria perpetua should be combined into 6404 P. crocata. This name (Miadlikowska et al. 2002) was first applied in the 2004 data.  Prior to this, it was most likely identified as 6404 P. crocata. 6408 P. perpetua has marginal soralia, a yellow medulla and tends to have a grayer upper cortex, while 6404 P. crocata has mostly laminal soredia, a white medulla and it tends to have a browner upper cortex.

# If the analysis will be comparing plots before and after 2003, then run this line

data[data$LICH_SPPCD==6404, ]
data[data$LICH_SPPCD==6408, 'LICH_SPPCD']<-6404 #None

# If the analysis will not be comparing plots before and after 2003, then run this line

data[(data$LICH_SPPCD==6408)&(data$MEASYEAR<2004), 'LICH_SPPCD']<-6404 #not run



## Sticta sp. For analyzing data from multiple years crossing 2004, 7507 Sticta carolinensis,and 7508 S. fragilinata should be combined into 7501 S. beauvoisii. Segregated in 2004 from 7501 S. beauvoisii, by having phyllidiate isidia (McDonald et al. 2003).  It differs from 7508 S. fragilinata in having smaller cyphellae with a white basal membrane and in having no secondary metabolites.   

# If the analysis compares plots before and after 2004, then run this line

data[data$LICH_SPPCD %in% c(7507,7508),'LICH_SPPCD']<-7501 #None


##Usnea sp. For analyzing data from multiple years crossing   2000, 8023 U. diplotypus should be combined into 8000 Usnea. ACTION 4: For analyzing pre-2000 data, 8023 Usnea diplotypus should be combined into 8000 Usnea. This species was used widely in PNW and CA 1998-99 as per concepts in McCune and Geiser, 1997 (Macrolichens of the Pacific Northwest).  In 2000, based on considerable progress in Usnea, this name was considered to be incorrectly applied in previous years. The name and code 8023 U. diplotypus should only be used for specimens identified with TLC.

# If the analysis uses plots from before 2000, then run this line

data[data$LICH_SPPCD==8023,'LICH_SPPCD']<-8000 #not run

# If the analysis compares plot from before and after 2000, then run this line

data[(data$LICH_SPPCD==8014) & (data$STATECD %in% weststates), 'LICH_SPPCD']<-8000 #None


## Usnea chaetophora. For analyzing data from multiple years crossing 1998, 8087 Usnea chaetophora should be combined into 8029 U. filipendula. In the western US, this taxon represents part of the U. filipendula group, which includes U. plicata and several others.  This name began to be applied widely in the western US in 1998 data. This species group replaces the old (pre-1998) taxon U. plicata.

# If the analysis compres plots before and after 1998, then run this line:

data[data$LICH_SPPCD==8087, 'LICH_SPPCD'] <- 8029 # 1


## Usnea esperantiana. For analyzing data from multiple years crossing 2000, 8088 Usnea esperantiana should be combined into 8036 U. glabrata. This taxon was first applied in 2000 in the west.  It is a recognizable syndrome, but possibly intergrades with 8036 U. glabrata

# If the analysis compares plots before and after 2000, then run this line:

data[data$LICH_SPPCD==8088,'LICH_SPPCD'] <- 8026 #None


## Xanthoria falax pre-1997. For data analysis for multiple years crossing 1997, 8210 Xanthomendoza fulva, 8219 X. galericulata, 8213 X. mendozae, 8215 X. oregana, and 8218 X. ulophyllodes should be combined into 8203 X. fallax. 8203 Xanthoria fallax = 8203 Xanthomendoza fallax. The latter name was used starting with 2004 data. Xanthoria fallax was applied broadly prior to Lindblom (1997). Starting with 1997 data, Lindblom's much narrower concept of X. fallax was applied. Western US: Before 1997, 8210 X. fulva, 8213 X. mendozae, 8215 X. oregana, and 8218 X. ulophyllodes were mostly included in the concept of 8203 X. fallax. In most western states 8210 X. fulva is common, so pre-1995 8203 X. fallax probably includes many specimens of 8210 X. fulva. Colorado 1992-1996 8203 X. fallax were re-examined; only a few were found to be 8210 X. fulva and no pre-1997 data entries were changed. Eastern US: 8210  X. fulva and 8218 X. ulophyllodes are moderately common, so pre-1997 X. fallax probably includes many specimens of 8210 X. fulva and 8218 X. ulophyllodes.  

# If the analysis compares data from before and after 1997, then run this line:

data[data$LICH_SPPCD %in% c(8210,8219,8213,8215,8218),'LICH_SPPCD'] <- 8203 # 153


## Xanthoria fulva. In pre-2004 data this taxon may have included 8219 X. galericulata (also see notes under 8215 Xanthomendoza oregana). Before 1997, 8210 X. fulva, if present, was mostly included with 8203 X. fallax. In 1997 and after, Lindblom's much narrower concept of 8203 X. fallax was applied and the name 8210 X. fulva was used more frequently. ACTION 3: For data analysis for multiple years crossing 2004, 8219 X. galericulata should be combined into 8210 X. fulva. 
## Xanthomendoza galericulata. For data analysis for multiple years crossing 2004, 8219 X. galericulata should be combined into 8210 X. fulva. This is a distinct species included in Nash et al. (2004...Sonoran Desert, v2) and Lindblom 2006. It differs from 8210 X. fulva and 8215 X. oregana in having soredia produced from the underside of hood-like or helmet-shaped lobe tips.  At this time we have no recommendation for how to compare with pre-2004 data.

# If the analysis compares data from before and after 2004, then run this line:

data[data$LICH_SPPCD == 8210,]
data[data$LICH_SPPCD == 8219, 'LICH_SPPCD'] <- 8210 #None


## Candelaria pacifica. For analyzing data from multiple years crossing 2002, 8303 Candelaria pacifica should be combined into 8301 C. concolor. This is a distinct species segregated from 8301 C. concolor  not officially named yet; in Nash et al. (2002...Sonoran Desert, v1). It differs from 8301 C. concolor in having greenish soredia from the underside, subfruticose habit, and an underside which mostly lacks cortex.

# If the analysis compares data before and after 2002, then run this line:

data[data$LICH_SPPCD == 8303,]
data[data$LICH_SPPCD == 8303, 'LICH_SPPCD'] <- 8301 #None


### Take care of Action Code 5 species (differences between east and western naming)


## Cladonia sp.  West - 1211 Cladonia coniocraea should be combined into 1228 C. ochrochlora. ACTION 2: East - 1228 C. ochrochlora should be combined into 1211 C. coniocraea. Western US:  In Colorado in 1992-1996, then 2000 data on, and in the West Coast region starting with 1998 data, the name 1228 C. ochrochlora was applied to include both 1211 C. coniocraea and 1228 C. ochrochlora. Eastern US: in 1993-1998 1211 C. coniocraea was applied as the default name; in the East it is the more common taxon. Only very characteristic specimens were labeled 1228 C. ochrochlora. Starting with 1999 data, a broader concept of 1228 C. ochrochlora has been applied; if any podetium of a specimen meets this species' criteria, the specimen is assigned 1228.

data[(data$LICH_SPPCD==1211) & (data$STATECD %in% weststates),'LICH_SPPCD']<-1228 #None

data[(data$LICH_SPPCD==1228) & (data$STATECD %in% eaststates),'LICH_SPPCD']<-1211 #None


## Phaeophyscia endococcina & endococcinodes. West - for analyzing data from multiple years crossing 2005, 5616 Phaeophyscia endococcina should be combined into 5618 Phaeophyscia endococcinodes. ACTION 4: West - for data analysis pre-2005, 5616 Phaeophyscia endococcina should be combined into 5618 Phaeophyscia endococcinodes. ACTION 0:  East. This species is considered exceedingly rare in North America. Esslinger considers this a different species from 5616 Phaeophyscia endococcinodes in Nash et al. (2004. Sonoran Desert, v2), though Moberg considers it a synonym. At this point the two species should be treated separately. 

# If the analysis will not be comparing plots before and after 2005, then run this line

data[(data$LICH_SPPCD==5616) & (data$STATECD %in% weststates) & (data$MEASYEAR<2005),'LICH_SPPCD']<-5618 #not run

# If the analysis will be comparing plots before and after 2005, then run this line

data[(data$LICH_SPPCD==5616) & (data$STATECD %in% weststates),'LICH_SPPCD']<-5618 #None


## Physconia sp. West - for analyzing data from multiple years crossing 1998, 5907 Physconia isidiigera, 5906 P. perisidiosa, and 5911 P. leucoleiptes should be combined into 5901 P. detersa. ACTION 0: East. Western US: Where 5901 P. detersa occurs in data prior to 1998 data, it represents a composite including one or all of these taxa: 5907 P. isidiigera, 5911 P. leucoleiptes and 5906 P. perisidiosa.  From 1998 onward, these taxa were distinguished in all data sets. Eastern US: P. detersa, P. leucoleiptes, and P. perisidiosa are distinguished in all data sets.

# If the analysis will be comparing plots before and after 1998, then run this line

data[(data$LICH_SPPCD %in% c(5907,5906,5911)) & (data$STATECD %in% weststates), 'LICH_SPPCD']<-5901 # 38

# If the analysis will not be comparing plots before and after 1998, then run this line

data[(data$LICH_SPPCD %in% c(5907,5906,5911)) & (data$STATECD %in% weststates) & (data$MEASYEAR < 1998), 'LICH_SPPCD']<-5901 #not run


## Xanthoria polycarpa & hasseana. West - for data analysis for multiple years crossing 1997, 8204 Xanthomendoza hasseana and 8214 X. montana should be combined into 8207 Xanthoria polycarpa group. ACTION 3: East - for data analysis for multiple years crossing 1997, 8207 X. polycarpa should be combined into 8204 X. hasseana. 8207 X. polycarpa was applied broadly prior to Lindblom (1997). Starting with 1997 data, 8207 X. polycarpa has been applied much more narrowly. Western US:  Prior to 1997, this taxon may have included 8207 X. polycarpa, 8214 X. montana, 8204 X. hasseana, and 8213 X. mendozae.  Colorado 1992-1996 specimens were examined and reassigned to 8214 X. montana. Also see notes for 8214 X. montana and 8204 X. hasseana. Eastern US: 8207 X. polycarpa for pre-1997 data is almost all the same taxon as 8204 X. hasseana, the name used for 1997 and later data.

# If the analysis will be comparing data from before and after 1997, then run these lines:

data[(data$LICH_SPPCD %in% c(8204,8214))&(data$STATECD %in% weststates),'LICH_SPPCD'] <- 8207 # 326
data[(data$LICH_SPPCD == 8207)&(data$STATECD %in% eaststates),'LICH_SPPCD'] <- 8204 # 3


### Taking care of Action Code 7 (complicated problems)

# Pyxine sp. For any analysis where Southern (S) FIA Region is included - for analyzing data from multiple years crossing 1999, 6801 Pyxine albovirens, 6803 P. caesiopruinosa, 6809 P. subcinerea should be combined into 6809 P. subcinerea. ACTION 4: For any analysis when Southern (S) FIA Region is included - for analyzing pre-2000 data, 6801 Pyxine albovirens, 6803 P. caesiopruinosa, 6809 P. subcinerea should be combined into 6809 P. subcinerea. Eastern US: The names 6801 P. albovirens, 6803  P. caesiopruinosa, and 6809 P. subcinerea may have been misapplied in early years in the SE. See Amtoft 2002 for a key to UV+ Pyxine in the East. 6803 P. caesiopruinosa should have a K+ purple medulla that is orange and 'marginal dactyls and coarse soredia'; found in the SE coastal plain. 

# If the analysis uses plots from southern states and compares plots before and after 1999, then run this line.  Note that non southern states may be biased toward higher richness.

data[(data$LICH_SPPCD %in% c(6801,6803)) & (data$STATECD %in% southstates), 'LICH_SPPCD']<-6809 #None

# If the analysis only uses data prior to 2000 in southern states, then run this line. Note that non southern states may be biased toward higher richness. 

data[(data$LICH_SPPCD %in% c(6801,6803)) & (data$STATECD %in% southstates) & (data$MEASYEAR < 2000), 'LICH_SPPCD']<-6809  #None


### Simple Species Changes
data2 = data

recode = read.csv('./Parsing Data/lichen_sppcd_recode.csv')
recode = recode[,-4]
recode$n.changes = NA

for(i in 1:nrow(recode)){
	if(sum(data$LICH_SPPCD == recode$LICH_SPPCD_IN[i]) > 0){
		recode$n.changes[i] <- sum(data$LICH_SPPCD == recode$LICH_SPPCD_IN[i])
		data2[data$LICH_SPPCD == recode$LICH_SPPCD_IN[i],'LICH_SPPCD'] <- recode[i,'LICH_SPPCD_OUT']
	} else { recode$n.changes[i] = 0 }
}

# remove all entries with NA as LICH_SPPCD
data2 = data2[!is.na(data2$LICH_SPPCD),]

write.csv(recode, paste('lichen_sppcd_recode_',Sys.Date(),'.csv',sep=''), row.names=F)

lichen2 = data2


### Take care of changes for lichen table where combinations have resulted in two rows for one species

# Find which plots have multiple rows for the same species
table(lichen2$yrplot.id, lichen2$LICH_SPPCD) -> plotXsp
doubles = which(plotXsp > 1, arr.ind=T)
get.plots = data.frame(yrplot.id = rownames(plotXsp)[doubles[,1]], LICH_SPPCD = colnames(plotXsp)[doubles[,2]])

lichen3 = lichen2
get.plots$rows.combined = NA

for(i in 1:nrow(get.plots)){

	# Subset out only lichen data that has multiple rows per species
	data.to.combine = subset(lichen3, (yrplot.id == get.plots[i,'yrplot.id'])&(LICH_SPPCD==get.plots[i,'LICH_SPPCD']))	
	get.plots$rows.combined[i] = nrow(data.to.combine)


	# Calculate the new abundance by combining the rows
	new.abun = combine.abundance(data.to.combine$ABUNDANCE_CLASS)

	# Delete all but one row and subsitute the new abundance
	remove.row = which((lichen3$yrplot.id == get.plots[i,'yrplot.id'])&(lichen3$LICH_SPPCD==get.plots[i,'LICH_SPPCD']))
	
	lichen3[remove.row[1],'ABUNDANCE_CLASS']<-new.abun
	lichen3 = lichen3[-remove.row[-1],]
}

# Checking that it worked:
sum(get.plots$rows.combined)-nrow(get.plots) + nrow(lichen3) == nrow(lichen2)



### Record what data looks like after changes:

table(lichen3$LICH_SPPCD)-> lichen.abun.post
lichen.abun.post = data.frame(LICH_SPPCD = names(lichen.abun.post), n.plots = as.vector(lichen.abun.post))
lichen.abun.post = merge(lichen.abun.post, spref[,c('LICH_SPPCD','YEARSTART','YEAREND','GENUS','SPECIES')], all.x=T)
lichen.abun.post = lichen.abun.post[order(lichen.abun.post$n.plots, decreasing=T),]
write.csv(lichen.abun.post, paste('lichen_abundance_combined_',Sys.Date(),'.csv',sep=''), row.names=F)


### Save changes
save(list=ls(all=T), file=paste('lichen_data_parsed_',Sys.Date(),'.Rdata',sep=''))

write.csv(lichen3,paste('FIA_lichens_parsed_',Sys.Date(),'.csv',sep=''), row.names=F)

