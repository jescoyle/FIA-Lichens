# FIA-Lichens
Code for manipulating and analyzing lichen data from the FIA.

-------------------------------------------------------------------
Directories:

Data
  Derived data tables for conducting final analyses.

OLD_scripts
  Scripts used in prior versions of the project, including SEMs and analysis of soil attributes.

cluster_scripts
  Scipts that were run on the computing cluster, including model dredging, sems, calculation of regional species richness
  
download_data
  Scripts for downloading data from FIA and CNALH websites and conducting quality control (e.g. taxonomic revisions).

format_data
  Scripts for calculating environmental covariates and richness estimates and combing FIA data sets.
  Includes: forest structure attributes, regional species richness, local richness, solar insolation, tree species composition, and climate variables.
  
-----------------------------------------------------------------
Scripts:

analysis.R
  Contains variance partitioning and model averagaing analysis reported in Coyle & Hurlbert "Environmental optimality, not heterogeneity, drives regional and local species richness in lichen epiphytes"
  
fia_lichen_analysis_function.R
  Contains functions used in analyses.

load_data.R
  Loads data needed to conduct analyses into an R workspace.

make_data_tables_for_analysis.R
  Combines lichen data and environmental data into tables for analysis (see Data directory).

make_richness_maps.R
  Plots local and regional macrolichen species richness as maps.
