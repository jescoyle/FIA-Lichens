## This script loads data needed for models of lichen diversity in FIA plots.


setwd('C://Users/jrcoyle/Documents/UNC/Projects/FIA Lichen')
options(stringsAsFactors=F)
varnames=read.csv('varnames.csv', row.names=1)

records_sp = read.csv('./Data/Regional Richness/CNALH_records_fia_genera_NAm.csv')

master = read.csv('./Data/fia_lichen_master_data_2015-09-19.csv')
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

