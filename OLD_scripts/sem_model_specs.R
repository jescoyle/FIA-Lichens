## This script contains model specifications for FIA Lichen SEMs (All taxa and functional richness)

#####################################################################################
## Models without pollution variables

# BaseMod, AllSp
path_nopol = "

	# Local environment/climate effects on forest structure
	bark_moist_pct.rao.ba ~ fh1cm1*wetness + fh1cm2*rain_lowRH + fh1cm3*iso + fh1cm4*pseas + fh1cm5*mat + fh1cm6*radiation
	wood_SG.rao.ba ~ fh2cm1*wetness + fh2cm2*rain_lowRH + fh2cm3*iso + fh2cm4*pseas + fh2cm5*mat + fh2cm6*radiation
	LogSeed.rao.ba ~ fh3cm1*wetness + fh3cm2*rain_lowRH + fh3cm3*iso + fh3cm4*pseas + fh3cm5*mat + fh3cm6*radiation
	PIE.ba.tree ~ fh4cm1*wetness + fh4cm2*rain_lowRH + fh4cm3*iso + fh4cm4*pseas + fh4cm5*mat + fh4cm6*radiation
	propDead ~ fh5cm1*wetness + fh5cm2*rain_lowRH + fh5cm3*iso + fh5cm4*pseas + fh5cm5*mat + fh5cm6*radiation
	lightDist.mean ~ fh6cm1*wetness + fh6cm2*rain_lowRH + fh6cm3*iso + fh6cm4*pseas + fh6cm5*mat + fh6cm6*radiation
	diamDiversity ~ fh7cm1*wetness + fh7cm2*rain_lowRH + fh7cm3*iso + fh7cm4*pseas + fh7cm5*mat + fh7cm6*radiation
	bark_moist_pct.ba ~ fm1cm1*wetness + fm1cm2*rain_lowRH + fm1cm3*iso + fm1cm4*pseas + fm1cm5*mat + fm1cm6*radiation
	wood_SG.ba ~ fm2cm1*wetness + fm2cm2*rain_lowRH + fm2cm3*iso + fm2cm4*pseas + fm2cm5*mat + fm2cm6*radiation
	LogSeed.ba ~ fm3cm1*wetness + fm3cm2*rain_lowRH + fm3cm3*iso + fm3cm4*pseas + fm3cm5*mat + fm3cm6*radiation
	bigTrees ~ fm4cm1*wetness + fm4cm2*rain_lowRH + fm4cm3*iso + fm4cm4*pseas + fm4cm5*mat + fm4cm6*radiation
	light.mean  ~ fm5cm1*wetness + fm5cm2*rain_lowRH + fm5cm3*iso + fm5cm4*pseas + fm5cm5*mat + fm5cm6*radiation
	PC1 ~ fm6cm1*wetness + fm6cm2*rain_lowRH + fm6cm3*iso + fm6cm4*pseas + fm6cm5*mat + fm6cm6*radiation

	# Regional climate, pollution, and regional forest heterogeneity effects on regional richness
	regS ~ R1CM1*wetness_reg_mean + R1CM2*rain_lowRH_reg_mean + R1CM3*iso_reg_mean + R1CM4*pseas_reg_mean + R1CM5*mat_reg_mean +
		R1CH1*wetness_reg_var + R1CH2*rain_lowRH_reg_var + R1CH3*iso_reg_var + R1CH4*pseas_reg_var + R1CH5*mat_reg_var +
		R1FH1*regS_tree 

	# Regional climate effects on regional forest heterogeneity
	regS_tree ~ FH1CM1*wetness_reg_mean + FH1CM2*rain_lowRH_reg_mean + FH1CM3*iso_reg_mean + FH1CM4*pseas_reg_mean + FH1CM5*mat_reg_mean +
		FH1CH1*wetness_reg_var + FH1CH2*rain_lowRH_reg_var + FH1CH3*iso_reg_var + FH1CH4*pseas_reg_var + FH1CH5*mat_reg_var

	# Effects on local lichen richness
	lichen.rich_log ~ r1fh1*bark_moist_pct.rao.ba + r1fh2*wood_SG.rao.ba + r1fh3*LogSeed.rao.ba +
		r1fh4*PIE.ba.tree + r1fh5*propDead + r1fh6*lightDist.mean + r1fh7*diamDiversity +
		r1fm1*bark_moist_pct.ba + r1fm2*wood_SG.ba + r1fm3*LogSeed.ba + r1fm4*bigTrees + r1fm5*light.mean + r1fm6*PC1 +
		r1cm1*wetness + r1cm2*rain_lowRH + r1cm3*iso + r1cm4*pseas + r1cm5*mat + r1cm6*radiation + r1R1*regS + r1a1*tot_abun_log
	
	# Effects on local lichen abundance
	tot_abun_log ~ a1fh1*bark_moist_pct.rao.ba + a1fh2*wood_SG.rao.ba + a1fh3*LogSeed.rao.ba +
		a1fh4*PIE.ba.tree + a1fh5*propDead + a1fh6*lightDist.mean + a1fh7*diamDiversity +
		a1fm1*bark_moist_pct.ba + a1fm2*wood_SG.ba + a1fm3*LogSeed.ba + a1fm4*bigTrees + a1fm5*light.mean + a1fm6*PC1 +
		a1cm1*wetness + a1cm2*rain_lowRH + a1cm3*iso + a1cm4*pseas + a1cm5*mat + a1cm6*radiation

	## Covariances between exogenous and endogenous variables
	# Dont need to specificy for Climate L-R b/c these are exogenous and are calculated automatically in lavaan

	# Local-regional forest structure
	regS_tree ~~ FH1fh4*PIE.ba.tree

	## Covariances among endogenous predictors in the same group
	# Dont need to specificy for Climate LO b/c these are exogenous and are calculated automatically in lavaan
	
	# Forest LH
	PIE.ba.tree ~~ bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba + propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba ~~ propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba ~~ wood_SG.rao.ba + LogSeed.rao.ba
	wood_SG.rao.ba ~~ LogSeed.rao.ba
	propDead ~~ lightDist.mean + diamDiversity
	lightDist.mean ~~ diamDiversity

	# Forest LO
	bark_moist_pct.ba + wood_SG.ba + LogSeed.ba ~~ bigTrees + light.mean + PC1
	bark_moist_pct.ba ~~ wood_SG.ba + LogSeed.ba 
	wood_SG.ba ~~ LogSeed.ba 
	bigTrees ~~ light.mean + PC1
	light.mean ~~ PC1
"

# BaseMod, Fric
path_nopol_fric = "
	
	# Local environment/climate effects on forest structure
	bark_moist_pct.rao.ba ~ fh1cm1*wetness + fh1cm2*rain_lowRH + fh1cm3*iso + fh1cm4*pseas + fh1cm5*mat + fh1cm6*radiation
	wood_SG.rao.ba ~ fh2cm1*wetness + fh2cm2*rain_lowRH + fh2cm3*iso + fh2cm4*pseas + fh2cm5*mat + fh2cm6*radiation
	LogSeed.rao.ba ~ fh3cm1*wetness + fh3cm2*rain_lowRH + fh3cm3*iso + fh3cm4*pseas + fh3cm5*mat + fh3cm6*radiation
	PIE.ba.tree ~ fh4cm1*wetness + fh4cm2*rain_lowRH + fh4cm3*iso + fh4cm4*pseas + fh4cm5*mat + fh4cm6*radiation
	propDead ~ fh5cm1*wetness + fh5cm2*rain_lowRH + fh5cm3*iso + fh5cm4*pseas + fh5cm5*mat + fh5cm6*radiation
	lightDist.mean ~ fh6cm1*wetness + fh6cm2*rain_lowRH + fh6cm3*iso + fh6cm4*pseas + fh6cm5*mat + fh6cm6*radiation
	diamDiversity ~ fh7cm1*wetness + fh7cm2*rain_lowRH + fh7cm3*iso + fh7cm4*pseas + fh7cm5*mat + fh7cm6*radiation
	bark_moist_pct.ba ~ fm1cm1*wetness + fm1cm2*rain_lowRH + fm1cm3*iso + fm1cm4*pseas + fm1cm5*mat + fm1cm6*radiation
	wood_SG.ba ~ fm2cm1*wetness + fm2cm2*rain_lowRH + fm2cm3*iso + fm2cm4*pseas + fm2cm5*mat + fm2cm6*radiation
	LogSeed.ba ~ fm3cm1*wetness + fm3cm2*rain_lowRH + fm3cm3*iso + fm3cm4*pseas + fm3cm5*mat + fm3cm6*radiation
	bigTrees ~ fm4cm1*wetness + fm4cm2*rain_lowRH + fm4cm3*iso + fm4cm4*pseas + fm4cm5*mat + fm4cm6*radiation
	light.mean  ~ fm5cm1*wetness + fm5cm2*rain_lowRH + fm5cm3*iso + fm5cm4*pseas + fm5cm5*mat + fm5cm6*radiation
	PC1 ~ fm6cm1*wetness + fm6cm2*rain_lowRH + fm6cm3*iso + fm6cm4*pseas + fm6cm5*mat + fm6cm6*radiation

	# Regional climate, pollution, and regional forest heterogeneity effects on regional richness
	regS ~ R1CM1*wetness_reg_mean + R1CM2*rain_lowRH_reg_mean + R1CM3*iso_reg_mean + R1CM4*pseas_reg_mean + R1CM5*mat_reg_mean +
		R1CH1*wetness_reg_var + R1CH2*rain_lowRH_reg_var + R1CH3*iso_reg_var + R1CH4*pseas_reg_var + R1CH5*mat_reg_var +
		R1FH1*regS_tree 

	# Regional climate effects on regional forest heterogeneity
	regS_tree ~ FH1CM1*wetness_reg_mean + FH1CM2*rain_lowRH_reg_mean + FH1CM3*iso_reg_mean + FH1CM4*pseas_reg_mean + FH1CM5*mat_reg_mean +
		FH1CH1*wetness_reg_var + FH1CH2*rain_lowRH_reg_var + FH1CH3*iso_reg_var + FH1CH4*pseas_reg_var + FH1CH5*mat_reg_var

	# Effects on local lichen richness
	fric ~ r1fh1*bark_moist_pct.rao.ba + r1fh2*wood_SG.rao.ba + r1fh3*LogSeed.rao.ba +
		r1fh4*PIE.ba.tree + r1fh5*propDead + r1fh6*lightDist.mean + r1fh7*diamDiversity +
		r1fm1*bark_moist_pct.ba + r1fm2*wood_SG.ba + r1fm3*LogSeed.ba + r1fm4*bigTrees + r1fm5*light.mean + r1fm6*PC1 +
		r1cm1*wetness + r1cm2*rain_lowRH + r1cm3*iso + r1cm4*pseas + r1cm5*mat + r1cm6*radiation + r1R1*regS + r1a1*tot_abun_log
	
	# Effects on local lichen abundance
	tot_abun_log ~ a1fh1*bark_moist_pct.rao.ba + a1fh2*wood_SG.rao.ba + a1fh3*LogSeed.rao.ba +
		a1fh4*PIE.ba.tree + a1fh5*propDead + a1fh6*lightDist.mean + a1fh7*diamDiversity +
		a1fm1*bark_moist_pct.ba + a1fm2*wood_SG.ba + a1fm3*LogSeed.ba + a1fm4*bigTrees + a1fm5*light.mean + a1fm6*PC1 +
		a1cm1*wetness + a1cm2*rain_lowRH + a1cm3*iso + a1cm4*pseas + a1cm5*mat + a1cm6*radiation

	## Covariances between exogenous and endogenous variables
	# Dont need to specificy for Climate L-R b/c these are exogenous and are calculated automatically in lavaan

	# Local-regional forest structure
	regS_tree ~~ FH1fh4*PIE.ba.tree

	## Covariances among endogenous predictors in the same group
	# Dont need to specificy for Climate LO b/c these are exogenous and are calculated automatically in lavaan
	
	# Forest LH
	PIE.ba.tree ~~ bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba + propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba ~~ propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba ~~ wood_SG.rao.ba + LogSeed.rao.ba
	wood_SG.rao.ba ~~ LogSeed.rao.ba
	propDead ~~ lightDist.mean + diamDiversity
	lightDist.mean ~~ diamDiversity

	# Forest LO
	bark_moist_pct.ba + wood_SG.ba + LogSeed.ba ~~ bigTrees + light.mean + PC1
	bark_moist_pct.ba ~~ wood_SG.ba + LogSeed.ba 
	wood_SG.ba ~~ LogSeed.ba 
	bigTrees ~~ light.mean + PC1
	light.mean ~~ PC1
"

# RegToRich, AllSp
path_regTorich_nopol = "

	# Local environment/climate effects on forest structure
	bark_moist_pct.rao.ba ~ fh1cm1*wetness + fh1cm2*rain_lowRH + fh1cm3*iso + fh1cm4*pseas + fh1cm5*mat + fh1cm6*radiation
	wood_SG.rao.ba ~ fh2cm1*wetness + fh2cm2*rain_lowRH + fh2cm3*iso + fh2cm4*pseas + fh2cm5*mat + fh2cm6*radiation
	LogSeed.rao.ba ~ fh3cm1*wetness + fh3cm2*rain_lowRH + fh3cm3*iso + fh3cm4*pseas + fh3cm5*mat + fh3cm6*radiation
	PIE.ba.tree ~ fh4cm1*wetness + fh4cm2*rain_lowRH + fh4cm3*iso + fh4cm4*pseas + fh4cm5*mat + fh4cm6*radiation
	propDead ~ fh5cm1*wetness + fh5cm2*rain_lowRH + fh5cm3*iso + fh5cm4*pseas + fh5cm5*mat + fh5cm6*radiation
	lightDist.mean ~ fh6cm1*wetness + fh6cm2*rain_lowRH + fh6cm3*iso + fh6cm4*pseas + fh6cm5*mat + fh6cm6*radiation
	diamDiversity ~ fh7cm1*wetness + fh7cm2*rain_lowRH + fh7cm3*iso + fh7cm4*pseas + fh7cm5*mat + fh7cm6*radiation
	bark_moist_pct.ba ~ fm1cm1*wetness + fm1cm2*rain_lowRH + fm1cm3*iso + fm1cm4*pseas + fm1cm5*mat + fm1cm6*radiation
	wood_SG.ba ~ fm2cm1*wetness + fm2cm2*rain_lowRH + fm2cm3*iso + fm2cm4*pseas + fm2cm5*mat + fm2cm6*radiation
	LogSeed.ba ~ fm3cm1*wetness + fm3cm2*rain_lowRH + fm3cm3*iso + fm3cm4*pseas + fm3cm5*mat + fm3cm6*radiation
	bigTrees ~ fm4cm1*wetness + fm4cm2*rain_lowRH + fm4cm3*iso + fm4cm4*pseas + fm4cm5*mat + fm4cm6*radiation
	light.mean  ~ fm5cm1*wetness + fm5cm2*rain_lowRH + fm5cm3*iso + fm5cm4*pseas + fm5cm5*mat + fm5cm6*radiation
	PC1 ~ fm6cm1*wetness + fm6cm2*rain_lowRH + fm6cm3*iso + fm6cm4*pseas + fm6cm5*mat + fm6cm6*radiation

	# Regional climate, pollution, and regional forest heterogeneity effects on regional richness
	regS ~ R1CM1*wetness_reg_mean + R1CM2*rain_lowRH_reg_mean + R1CM3*iso_reg_mean + R1CM4*pseas_reg_mean + R1CM5*mat_reg_mean +
		R1CH1*wetness_reg_var + R1CH2*rain_lowRH_reg_var + R1CH3*iso_reg_var + R1CH4*pseas_reg_var + R1CH5*mat_reg_var +
		R1FH1*regS_tree

	# Regional climate effects on regional forest heterogeneity
	regS_tree ~ FH1CM1*wetness_reg_mean + FH1CM2*rain_lowRH_reg_mean + FH1CM3*iso_reg_mean + FH1CM4*pseas_reg_mean + FH1CM5*mat_reg_mean +
		FH1CH1*wetness_reg_var + FH1CH2*rain_lowRH_reg_var + FH1CH3*iso_reg_var + FH1CH4*pseas_reg_var + FH1CH5*mat_reg_var

	# Effects on local lichen richness
	lichen.rich_log ~ r1fh1*bark_moist_pct.rao.ba + r1fh2*wood_SG.rao.ba + r1fh3*LogSeed.rao.ba +
		r1fh4*PIE.ba.tree + r1fh5*propDead + r1fh6*lightDist.mean + r1fh7*diamDiversity +
		r1fm1*bark_moist_pct.ba + r1fm2*wood_SG.ba + r1fm3*LogSeed.ba + r1fm4*bigTrees + r1fm5*light.mean + r1fm6*PC1 +
		r1cm1*wetness + r1cm2*rain_lowRH + r1cm3*iso + r1cm4*pseas + r1cm5*mat + r1cm6*radiation + r1R1*regS +
		r1CM1*wetness_reg_mean + r1CM3*rain_lowRH_reg_mean + r1CM3*iso_reg_mean + r1CM4*pseas_reg_mean + r1CM5*mat_reg_mean +
		r1CH1*wetness_reg_var + r1CH2*rain_lowRH_reg_var + r1CH3*iso_reg_var + r1CH4*pseas_reg_var + r1CH5*mat_reg_var +
		r1FH1*regS_tree + r1a1*tot_abun_log
	
	# Effects on local lichen abundance
	tot_abun_log ~ a1fh1*bark_moist_pct.rao.ba + a1fh2*wood_SG.rao.ba + a1fh3*LogSeed.rao.ba +
		a1fh4*PIE.ba.tree + a1fh5*propDead + a1fh6*lightDist.mean + a1fh7*diamDiversity +
		a1fm1*bark_moist_pct.ba + a1fm2*wood_SG.ba + a1fm3*LogSeed.ba + a1fm4*bigTrees + a1fm5*light.mean + a1fm6*PC1 +
		a1cm1*wetness + a1cm2*rain_lowRH + a1cm3*iso + a1cm4*pseas + a1cm5*mat + a1cm6*radiation

	## Covariances between exogenous and endogenous variables
	# Dont need to specificy for Climate L-R b/c these are exogenous and are calculated automatically in lavaan

	# Local-regional forest structure
	regS_tree ~~ FH1fh4*PIE.ba.tree

	## Covariances among endogenous predictors in the same group
	# Dont need to specificy for Climate LO b/c these are exogenous and are calculated automatically in lavaan
	
	# Forest LH
	PIE.ba.tree ~~ bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba + propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba ~~ propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba ~~ wood_SG.rao.ba + LogSeed.rao.ba
	wood_SG.rao.ba ~~ LogSeed.rao.ba
	propDead ~~ lightDist.mean + diamDiversity
	lightDist.mean ~~ diamDiversity

	# Forest LO
	bark_moist_pct.ba + wood_SG.ba + LogSeed.ba ~~ bigTrees + light.mean + PC1
	bark_moist_pct.ba ~~ wood_SG.ba + LogSeed.ba 
	wood_SG.ba ~~ LogSeed.ba 
	bigTrees ~~ light.mean + PC1
	light.mean ~~ PC1
"


# RegToRich, Fric
path_regTorich_nopol_fric = "

	# Local environment/climate effects on forest structure
	bark_moist_pct.rao.ba ~ fh1cm1*wetness + fh1cm2*rain_lowRH + fh1cm3*iso + fh1cm4*pseas + fh1cm5*mat + fh1cm6*radiation
	wood_SG.rao.ba ~ fh2cm1*wetness + fh2cm2*rain_lowRH + fh2cm3*iso + fh2cm4*pseas + fh2cm5*mat + fh2cm6*radiation
	LogSeed.rao.ba ~ fh3cm1*wetness + fh3cm2*rain_lowRH + fh3cm3*iso + fh3cm4*pseas + fh3cm5*mat + fh3cm6*radiation
	PIE.ba.tree ~ fh4cm1*wetness + fh4cm2*rain_lowRH + fh4cm3*iso + fh4cm4*pseas + fh4cm5*mat + fh4cm6*radiation
	propDead ~ fh5cm1*wetness + fh5cm2*rain_lowRH + fh5cm3*iso + fh5cm4*pseas + fh5cm5*mat + fh5cm6*radiation
	lightDist.mean ~ fh6cm1*wetness + fh6cm2*rain_lowRH + fh6cm3*iso + fh6cm4*pseas + fh6cm5*mat + fh6cm6*radiation
	diamDiversity ~ fh7cm1*wetness + fh7cm2*rain_lowRH + fh7cm3*iso + fh7cm4*pseas + fh7cm5*mat + fh7cm6*radiation
	bark_moist_pct.ba ~ fm1cm1*wetness + fm1cm2*rain_lowRH + fm1cm3*iso + fm1cm4*pseas + fm1cm5*mat + fm1cm6*radiation
	wood_SG.ba ~ fm2cm1*wetness + fm2cm2*rain_lowRH + fm2cm3*iso + fm2cm4*pseas + fm2cm5*mat + fm2cm6*radiation
	LogSeed.ba ~ fm3cm1*wetness + fm3cm2*rain_lowRH + fm3cm3*iso + fm3cm4*pseas + fm3cm5*mat + fm3cm6*radiation
	bigTrees ~ fm4cm1*wetness + fm4cm2*rain_lowRH + fm4cm3*iso + fm4cm4*pseas + fm4cm5*mat + fm4cm6*radiation
	light.mean  ~ fm5cm1*wetness + fm5cm2*rain_lowRH + fm5cm3*iso + fm5cm4*pseas + fm5cm5*mat + fm5cm6*radiation
	PC1 ~ fm6cm1*wetness + fm6cm2*rain_lowRH + fm6cm3*iso + fm6cm4*pseas + fm6cm5*mat + fm6cm6*radiation

	# Regional climate, pollution, and regional forest heterogeneity effects on regional richness
	regS ~ R1CM1*wetness_reg_mean + R1CM2*rain_lowRH_reg_mean + R1CM3*iso_reg_mean + R1CM4*pseas_reg_mean + R1CM5*mat_reg_mean +
		R1CH1*wetness_reg_var + R1CH2*rain_lowRH_reg_var + R1CH3*iso_reg_var + R1CH4*pseas_reg_var + R1CH5*mat_reg_var +
		R1FH1*regS_tree

	# Regional climate effects on regional forest heterogeneity
	regS_tree ~ FH1CM1*wetness_reg_mean + FH1CM2*rain_lowRH_reg_mean + FH1CM3*iso_reg_mean + FH1CM4*pseas_reg_mean + FH1CM5*mat_reg_mean +
		FH1CH1*wetness_reg_var + FH1CH2*rain_lowRH_reg_var + FH1CH3*iso_reg_var + FH1CH4*pseas_reg_var + FH1CH5*mat_reg_var

	# Effects on local lichen richness
	fric ~ r1fh1*bark_moist_pct.rao.ba + r1fh2*wood_SG.rao.ba + r1fh3*LogSeed.rao.ba +
		r1fh4*PIE.ba.tree + r1fh5*propDead + r1fh6*lightDist.mean + r1fh7*diamDiversity +
		r1fm1*bark_moist_pct.ba + r1fm2*wood_SG.ba + r1fm3*LogSeed.ba + r1fm4*bigTrees + r1fm5*light.mean + r1fm6*PC1 +
		r1cm1*wetness + r1cm2*rain_lowRH + r1cm3*iso + r1cm4*pseas + r1cm5*mat + r1cm6*radiation + r1R1*regS +
		r1CM1*wetness_reg_mean + r1CM3*rain_lowRH_reg_mean + r1CM3*iso_reg_mean + r1CM4*pseas_reg_mean + r1CM5*mat_reg_mean +
		r1CH1*wetness_reg_var + r1CH2*rain_lowRH_reg_var + r1CH3*iso_reg_var + r1CH4*pseas_reg_var + r1CH5*mat_reg_var +
		r1FH1*regS_tree + r1a1*tot_abun_log
	
	# Effects on local lichen abundance
	tot_abun_log ~ a1fh1*bark_moist_pct.rao.ba + a1fh2*wood_SG.rao.ba + a1fh3*LogSeed.rao.ba +
		a1fh4*PIE.ba.tree + a1fh5*propDead + a1fh6*lightDist.mean + a1fh7*diamDiversity +
		a1fm1*bark_moist_pct.ba + a1fm2*wood_SG.ba + a1fm3*LogSeed.ba + a1fm4*bigTrees + a1fm5*light.mean + a1fm6*PC1 +
		a1cm1*wetness + a1cm2*rain_lowRH + a1cm3*iso + a1cm4*pseas + a1cm5*mat + a1cm6*radiation

	## Covariances between exogenous and endogenous variables
	# Dont need to specificy for Climate L-R b/c these are exogenous and are calculated automatically in lavaan

	# Local-regional forest structure
	regS_tree ~~ FH1fh4*PIE.ba.tree

	## Covariances among endogenous predictors in the same group
	# Dont need to specificy for Climate LO b/c these are exogenous and are calculated automatically in lavaan
	
	# Forest LH
	PIE.ba.tree ~~ bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba + propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba ~~ propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba ~~ wood_SG.rao.ba + LogSeed.rao.ba
	wood_SG.rao.ba ~~ LogSeed.rao.ba
	propDead ~~ lightDist.mean + diamDiversity
	lightDist.mean ~~ diamDiversity

	# Forest LO
	bark_moist_pct.ba + wood_SG.ba + LogSeed.ba ~~ bigTrees + light.mean + PC1
	bark_moist_pct.ba ~~ wood_SG.ba + LogSeed.ba 
	wood_SG.ba ~~ LogSeed.ba 
	bigTrees ~~ light.mean + PC1
	light.mean ~~ PC1
"

# Recip, AllSp
path_nopol_recip = "

	# Local environment/climate effects on forest structure
	bark_moist_pct.rao.ba ~ fh1cm1*wetness + fh1cm2*rain_lowRH + fh1cm3*iso + fh1cm4*pseas + fh1cm5*mat + fh1cm6*radiation
	wood_SG.rao.ba ~ fh2cm1*wetness + fh2cm2*rain_lowRH + fh2cm3*iso + fh2cm4*pseas + fh2cm5*mat + fh2cm6*radiation
	LogSeed.rao.ba ~ fh3cm1*wetness + fh3cm2*rain_lowRH + fh3cm3*iso + fh3cm4*pseas + fh3cm5*mat + fh3cm6*radiation
	PIE.ba.tree ~ fh4cm1*wetness + fh4cm2*rain_lowRH + fh4cm3*iso + fh4cm4*pseas + fh4cm5*mat + fh4cm6*radiation
	propDead ~ fh5cm1*wetness + fh5cm2*rain_lowRH + fh5cm3*iso + fh5cm4*pseas + fh5cm5*mat + fh5cm6*radiation
	lightDist.mean ~ fh6cm1*wetness + fh6cm2*rain_lowRH + fh6cm3*iso + fh6cm4*pseas + fh6cm5*mat + fh6cm6*radiation
	diamDiversity ~ fh7cm1*wetness + fh7cm2*rain_lowRH + fh7cm3*iso + fh7cm4*pseas + fh7cm5*mat + fh7cm6*radiation
	bark_moist_pct.ba ~ fm1cm1*wetness + fm1cm2*rain_lowRH + fm1cm3*iso + fm1cm4*pseas + fm1cm5*mat + fm1cm6*radiation
	wood_SG.ba ~ fm2cm1*wetness + fm2cm2*rain_lowRH + fm2cm3*iso + fm2cm4*pseas + fm2cm5*mat + fm2cm6*radiation
	LogSeed.ba ~ fm3cm1*wetness + fm3cm2*rain_lowRH + fm3cm3*iso + fm3cm4*pseas + fm3cm5*mat + fm3cm6*radiation
	bigTrees ~ fm4cm1*wetness + fm4cm2*rain_lowRH + fm4cm3*iso + fm4cm4*pseas + fm4cm5*mat + fm4cm6*radiation
	light.mean  ~ fm5cm1*wetness + fm5cm2*rain_lowRH + fm5cm3*iso + fm5cm4*pseas + fm5cm5*mat + fm5cm6*radiation
	PC1 ~ fm6cm1*wetness + fm6cm2*rain_lowRH + fm6cm3*iso + fm6cm4*pseas + fm6cm5*mat + fm6cm6*radiation

	# Regional climate, pollution, and regional forest heterogeneity effects on regional richness
	regS ~ R1CM1*wetness_reg_mean + R1CM2*rain_lowRH_reg_mean + R1CM3*iso_reg_mean + R1CM4*pseas_reg_mean + R1CM5*mat_reg_mean +
		R1CH1*wetness_reg_var + R1CH2*rain_lowRH_reg_var + R1CH3*iso_reg_var + R1CH4*pseas_reg_var + R1CH5*mat_reg_var +
		R1FH1*regS_tree + R1r1*lichen.rich_log

	# Regional climate effects on regional forest heterogeneity
	regS_tree ~ FH1CM1*wetness_reg_mean + FH1CM2*rain_lowRH_reg_mean + FH1CM3*iso_reg_mean + FH1CM4*pseas_reg_mean + FH1CM5*mat_reg_mean +
		FH1CH1*wetness_reg_var + FH1CH2*rain_lowRH_reg_var + FH1CH3*iso_reg_var + FH1CH4*pseas_reg_var + FH1CH5*mat_reg_var

	# Effects on local lichen richness
	lichen.rich_log ~ r1fh1*bark_moist_pct.rao.ba + r1fh2*wood_SG.rao.ba + r1fh3*LogSeed.rao.ba +
		r1fh4*PIE.ba.tree + r1fh5*propDead + r1fh6*lightDist.mean + r1fh7*diamDiversity +
		r1fm1*bark_moist_pct.ba + r1fm2*wood_SG.ba + r1fm3*LogSeed.ba + r1fm4*bigTrees + r1fm5*light.mean + r1fm6*PC1 +
		r1cm1*wetness + r1cm2*rain_lowRH + r1cm3*iso + r1cm4*pseas + r1cm5*mat + r1cm6*radiation + r1R1*regS + r1a1*tot_abun_log
	
	# Effects on local lichen abundance
	tot_abun_log ~ a1fh1*bark_moist_pct.rao.ba + a1fh2*wood_SG.rao.ba + a1fh3*LogSeed.rao.ba +
		a1fh4*PIE.ba.tree + a1fh5*propDead + a1fh6*lightDist.mean + a1fh7*diamDiversity +
		a1fm1*bark_moist_pct.ba + a1fm2*wood_SG.ba + a1fm3*LogSeed.ba + a1fm4*bigTrees + a1fm5*light.mean + a1fm6*PC1 +
		a1cm1*wetness + a1cm2*rain_lowRH + a1cm3*iso + a1cm4*pseas + a1cm5*mat + a1cm6*radiation

	## Covariances between exogenous and endogenous variables
	# Dont need to specificy for Climate L-R b/c these are exogenous and are calculated automatically in lavaan

	# Local-regional forest structure
	regS_tree ~~ FH1fh4*PIE.ba.tree

	## Covariances among endogenous predictors in the same group
	# Dont need to specificy for Climate LO b/c these are exogenous and are calculated automatically in lavaan
	
	# Forest LH
	PIE.ba.tree ~~ bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba + propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba ~~ propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba ~~ wood_SG.rao.ba + LogSeed.rao.ba
	wood_SG.rao.ba ~~ LogSeed.rao.ba
	propDead ~~ lightDist.mean + diamDiversity
	lightDist.mean ~~ diamDiversity

	# Forest LO
	bark_moist_pct.ba + wood_SG.ba + LogSeed.ba ~~ bigTrees + light.mean + PC1
	bark_moist_pct.ba ~~ wood_SG.ba + LogSeed.ba 
	wood_SG.ba ~~ LogSeed.ba 
	bigTrees ~~ light.mean + PC1
	light.mean ~~ PC1
"

# RegToRich, Recip, AllSp
path_regTorich_nopol_recip = "

	# Local environment/climate effects on forest structure
	bark_moist_pct.rao.ba ~ fh1cm1*wetness + fh1cm2*rain_lowRH + fh1cm3*iso + fh1cm4*pseas + fh1cm5*mat + fh1cm6*radiation
	wood_SG.rao.ba ~ fh2cm1*wetness + fh2cm2*rain_lowRH + fh2cm3*iso + fh2cm4*pseas + fh2cm5*mat + fh2cm6*radiation
	LogSeed.rao.ba ~ fh3cm1*wetness + fh3cm2*rain_lowRH + fh3cm3*iso + fh3cm4*pseas + fh3cm5*mat + fh3cm6*radiation
	PIE.ba.tree ~ fh4cm1*wetness + fh4cm2*rain_lowRH + fh4cm3*iso + fh4cm4*pseas + fh4cm5*mat + fh4cm6*radiation
	propDead ~ fh5cm1*wetness + fh5cm2*rain_lowRH + fh5cm3*iso + fh5cm4*pseas + fh5cm5*mat + fh5cm6*radiation
	lightDist.mean ~ fh6cm1*wetness + fh6cm2*rain_lowRH + fh6cm3*iso + fh6cm4*pseas + fh6cm5*mat + fh6cm6*radiation
	diamDiversity ~ fh7cm1*wetness + fh7cm2*rain_lowRH + fh7cm3*iso + fh7cm4*pseas + fh7cm5*mat + fh7cm6*radiation
	bark_moist_pct.ba ~ fm1cm1*wetness + fm1cm2*rain_lowRH + fm1cm3*iso + fm1cm4*pseas + fm1cm5*mat + fm1cm6*radiation
	wood_SG.ba ~ fm2cm1*wetness + fm2cm2*rain_lowRH + fm2cm3*iso + fm2cm4*pseas + fm2cm5*mat + fm2cm6*radiation
	LogSeed.ba ~ fm3cm1*wetness + fm3cm2*rain_lowRH + fm3cm3*iso + fm3cm4*pseas + fm3cm5*mat + fm3cm6*radiation
	bigTrees ~ fm4cm1*wetness + fm4cm2*rain_lowRH + fm4cm3*iso + fm4cm4*pseas + fm4cm5*mat + fm4cm6*radiation
	light.mean  ~ fm5cm1*wetness + fm5cm2*rain_lowRH + fm5cm3*iso + fm5cm4*pseas + fm5cm5*mat + fm5cm6*radiation
	PC1 ~ fm6cm1*wetness + fm6cm2*rain_lowRH + fm6cm3*iso + fm6cm4*pseas + fm6cm5*mat + fm6cm6*radiation

	# Regional climate, pollution, and regional forest heterogeneity effects on regional richness
	regS ~ R1CM1*wetness_reg_mean + R1CM2*rain_lowRH_reg_mean + R1CM3*iso_reg_mean + R1CM4*pseas_reg_mean + R1CM5*mat_reg_mean +
		R1CH1*wetness_reg_var + R1CH2*rain_lowRH_reg_var + R1CH3*iso_reg_var + R1CH4*pseas_reg_var + R1CH5*mat_reg_var +
		R1FH1*regS_tree + R1r1*lichen.rich_log

	# Regional climate effects on regional forest heterogeneity
	regS_tree ~ FH1CM1*wetness_reg_mean + FH1CM2*rain_lowRH_reg_mean + FH1CM3*iso_reg_mean + FH1CM4*pseas_reg_mean + FH1CM5*mat_reg_mean +
		FH1CH1*wetness_reg_var + FH1CH2*rain_lowRH_reg_var + FH1CH3*iso_reg_var + FH1CH4*pseas_reg_var + FH1CH5*mat_reg_var

	# Effects on local lichen richness
	lichen.rich_log ~ r1fh1*bark_moist_pct.rao.ba + r1fh2*wood_SG.rao.ba + r1fh3*LogSeed.rao.ba +
		r1fh4*PIE.ba.tree + r1fh5*propDead + r1fh6*lightDist.mean + r1fh7*diamDiversity +
		r1fm1*bark_moist_pct.ba + r1fm2*wood_SG.ba + r1fm3*LogSeed.ba + r1fm4*bigTrees + r1fm5*light.mean + r1fm6*PC1 +
		r1cm1*wetness + r1cm2*rain_lowRH + r1cm3*iso + r1cm4*pseas + r1cm5*mat + r1cm6*radiation + r1R1*regS +
		r1CM1*wetness_reg_mean + r1CM3*rain_lowRH_reg_mean + r1CM3*iso_reg_mean + r1CM4*pseas_reg_mean + r1CM5*mat_reg_mean +
		r1CH1*wetness_reg_var + r1CH2*rain_lowRH_reg_var + r1CH3*iso_reg_var + r1CH4*pseas_reg_var + r1CH5*mat_reg_var +
		r1FH1*regS_tree + r1a1*tot_abun_log
	
	# Effects on local lichen abundance
	tot_abun_log ~ a1fh1*bark_moist_pct.rao.ba + a1fh2*wood_SG.rao.ba + a1fh3*LogSeed.rao.ba +
		a1fh4*PIE.ba.tree + a1fh5*propDead + a1fh6*lightDist.mean + a1fh7*diamDiversity +
		a1fm1*bark_moist_pct.ba + a1fm2*wood_SG.ba + a1fm3*LogSeed.ba + a1fm4*bigTrees + a1fm5*light.mean + a1fm6*PC1 +
		a1cm1*wetness + a1cm2*rain_lowRH + a1cm3*iso + a1cm4*pseas + a1cm5*mat + a1cm6*radiation

	## Covariances between exogenous and endogenous variables
	# Dont need to specificy for Climate L-R b/c these are exogenous and are calculated automatically in lavaan

	# Local-regional forest structure
	regS_tree ~~ FH1fh4*PIE.ba.tree

	## Covariances among endogenous predictors in the same group
	# Dont need to specificy for Climate LO b/c these are exogenous and are calculated automatically in lavaan
	
	# Forest LH
	PIE.ba.tree ~~ bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba + propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba ~~ propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba ~~ wood_SG.rao.ba + LogSeed.rao.ba
	wood_SG.rao.ba ~~ LogSeed.rao.ba
	propDead ~~ lightDist.mean + diamDiversity
	lightDist.mean ~~ diamDiversity

	# Forest LO
	bark_moist_pct.ba + wood_SG.ba + LogSeed.ba ~~ bigTrees + light.mean + PC1
	bark_moist_pct.ba ~~ wood_SG.ba + LogSeed.ba 
	wood_SG.ba ~~ LogSeed.ba 
	bigTrees ~~ light.mean + PC1
	light.mean ~~ PC1
"


# No Abundance, AllSp
path_noabun_nopol = "
	
	# Local environment/climate effects on forest structure
	bark_moist_pct.rao.ba ~ fh1cm1*wetness + fh1cm2*rain_lowRH + fh1cm3*iso + fh1cm4*pseas + fh1cm5*mat + fh1cm6*radiation
	wood_SG.rao.ba ~ fh2cm1*wetness + fh2cm2*rain_lowRH + fh2cm3*iso + fh2cm4*pseas + fh2cm5*mat + fh2cm6*radiation
	LogSeed.rao.ba ~ fh3cm1*wetness + fh3cm2*rain_lowRH + fh3cm3*iso + fh3cm4*pseas + fh3cm5*mat + fh3cm6*radiation
	PIE.ba.tree ~ fh4cm1*wetness + fh4cm2*rain_lowRH + fh4cm3*iso + fh4cm4*pseas + fh4cm5*mat + fh4cm6*radiation
	propDead ~ fh5cm1*wetness + fh5cm2*rain_lowRH + fh5cm3*iso + fh5cm4*pseas + fh5cm5*mat + fh5cm6*radiation
	lightDist.mean ~ fh6cm1*wetness + fh6cm2*rain_lowRH + fh6cm3*iso + fh6cm4*pseas + fh6cm5*mat + fh6cm6*radiation
	diamDiversity ~ fh7cm1*wetness + fh7cm2*rain_lowRH + fh7cm3*iso + fh7cm4*pseas + fh7cm5*mat + fh7cm6*radiation
	bark_moist_pct.ba ~ fm1cm1*wetness + fm1cm2*rain_lowRH + fm1cm3*iso + fm1cm4*pseas + fm1cm5*mat + fm1cm6*radiation
	wood_SG.ba ~ fm2cm1*wetness + fm2cm2*rain_lowRH + fm2cm3*iso + fm2cm4*pseas + fm2cm5*mat + fm2cm6*radiation
	LogSeed.ba ~ fm3cm1*wetness + fm3cm2*rain_lowRH + fm3cm3*iso + fm3cm4*pseas + fm3cm5*mat + fm3cm6*radiation
	bigTrees ~ fm4cm1*wetness + fm4cm2*rain_lowRH + fm4cm3*iso + fm4cm4*pseas + fm4cm5*mat + fm4cm6*radiation
	light.mean  ~ fm5cm1*wetness + fm5cm2*rain_lowRH + fm5cm3*iso + fm5cm4*pseas + fm5cm5*mat + fm5cm6*radiation
	PC1 ~ fm6cm1*wetness + fm6cm2*rain_lowRH + fm6cm3*iso + fm6cm4*pseas + fm6cm5*mat + fm6cm6*radiation

	# Regional climate, pollution, and regional forest heterogeneity effects on regional richness
	regS ~ R1CM1*wetness_reg_mean + R1CM2*rain_lowRH_reg_mean + R1CM3*iso_reg_mean + R1CM4*pseas_reg_mean + R1CM5*mat_reg_mean +
		R1CH1*wetness_reg_var + R1CH2*rain_lowRH_reg_var + R1CH3*iso_reg_var + R1CH4*pseas_reg_var + R1CH5*mat_reg_var +
		R1FH1*regS_tree 

	# Regional climate effects on regional forest heterogeneity
	regS_tree ~ FH1CM1*wetness_reg_mean + FH1CM2*rain_lowRH_reg_mean + FH1CM3*iso_reg_mean + FH1CM4*pseas_reg_mean + FH1CM5*mat_reg_mean +
		FH1CH1*wetness_reg_var + FH1CH2*rain_lowRH_reg_var + FH1CH3*iso_reg_var + FH1CH4*pseas_reg_var + FH1CH5*mat_reg_var

	# Effects on local lichen richness
	lichen.rich_log ~ r1fh1*bark_moist_pct.rao.ba + r1fh2*wood_SG.rao.ba + r1fh3*LogSeed.rao.ba +
		r1fh4*PIE.ba.tree + r1fh5*propDead + r1fh6*lightDist.mean + r1fh7*diamDiversity +
		r1fm1*bark_moist_pct.ba + r1fm2*wood_SG.ba + r1fm3*LogSeed.ba + r1fm4*bigTrees + r1fm5*light.mean + r1fm6*PC1 +
		r1cm1*wetness + r1cm2*rain_lowRH + r1cm3*iso + r1cm4*pseas + r1cm5*mat + r1cm6*radiation + r1R1*regS	
	
	## Covariances between exogenous and endogenous variables
	# Dont need to specificy for Climate L-R b/c these are exogenous and are calculated automatically in lavaan

	# Local-regional forest structure
	regS_tree ~~ FH1fh4*PIE.ba.tree

	## Covariances among endogenous predictors in the same group
	# Dont need to specificy for Climate LO b/c these are exogenous and are calculated automatically in lavaan
	
	# Forest LH
	PIE.ba.tree ~~ bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba + propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba ~~ propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba ~~ wood_SG.rao.ba + LogSeed.rao.ba
	wood_SG.rao.ba ~~ LogSeed.rao.ba
	propDead ~~ lightDist.mean + diamDiversity
	lightDist.mean ~~ diamDiversity

	# Forest LO
	bark_moist_pct.ba + wood_SG.ba + LogSeed.ba ~~ bigTrees + light.mean + PC1
	bark_moist_pct.ba ~~ wood_SG.ba + LogSeed.ba 
	wood_SG.ba ~~ LogSeed.ba 
	bigTrees ~~ light.mean + PC1
	light.mean ~~ PC1
"

# SoilMod, AllSp
path_soilmod_nopol = "

	# Local environment/climate effects on forest structure
	bark_moist_pct.rao.ba ~ fh1cm1*wetness + fh1cm2*rain_lowRH + fh1cm3*iso + fh1cm4*pseas + fh1cm5*mat + fh1cm6*radiation + fh1sm1*soilPC1 + fh1sm2*soilPC2
	wood_SG.rao.ba ~ fh2cm1*wetness + fh2cm2*rain_lowRH + fh2cm3*iso + fh2cm4*pseas + fh2cm5*mat + fh2cm6*radiation + fh2sm1*soilPC1 + fh2sm2*soilPC2
	LogSeed.rao.ba ~ fh3cm1*wetness + fh3cm2*rain_lowRH + fh3cm3*iso + fh3cm4*pseas + fh3cm5*mat + fh3cm6*radiation + fh3sm1*soilPC1 + fh3sm2*soilPC2
	PIE.ba.tree ~ fh4cm1*wetness + fh4cm2*rain_lowRH + fh4cm3*iso + fh4cm4*pseas + fh4cm5*mat + fh4cm6*radiation + fh4sm1*soilPC1 + fh4sm2*soilPC2
	propDead ~ fh5cm1*wetness + fh5cm2*rain_lowRH + fh5cm3*iso + fh5cm4*pseas + fh5cm5*mat + fh5cm6*radiation + fh5sm1*soilPC1 + fh5sm2*soilPC2
	lightDist.mean ~ fh6cm1*wetness + fh6cm2*rain_lowRH + fh6cm3*iso + fh6cm4*pseas + fh6cm5*mat + fh6cm6*radiation + fh6sm1*soilPC1 + fh6sm2*soilPC2
	diamDiversity ~ fh7cm1*wetness + fh7cm2*rain_lowRH + fh7cm3*iso + fh7cm4*pseas + fh7cm5*mat + fh7cm6*radiation + fh7sm1*soilPC1 + fh7sm2*soilPC2
	bark_moist_pct.ba ~ fm1cm1*wetness + fm1cm2*rain_lowRH + fm1cm3*iso + fm1cm4*pseas + fm1cm5*mat + fm1cm6*radiation + fm1sm1*soilPC1 + fm1sm2*soilPC2
	wood_SG.ba ~ fm2cm1*wetness + fm2cm2*rain_lowRH + fm2cm3*iso + fm2cm4*pseas + fm2cm5*mat + fm2cm6*radiation + fm2sm1*soilPC1 + fm2sm2*soilPC2
	LogSeed.ba ~ fm3cm1*wetness + fm3cm2*rain_lowRH + fm3cm3*iso + fm3cm4*pseas + fm3cm5*mat + fm3cm6*radiation + fm3sm1*soilPC1 + fm3sm2*soilPC2
	bigTrees ~ fm4cm1*wetness + fm4cm2*rain_lowRH + fm4cm3*iso + fm4cm4*pseas + fm4cm5*mat + fm4cm6*radiation + fm4sm1*soilPC1 + fm4sm2*soilPC2
	light.mean  ~ fm5cm1*wetness + fm5cm2*rain_lowRH + fm5cm3*iso + fm5cm4*pseas + fm5cm5*mat + fm5cm6*radiation + fm5sm1*soilPC1 + fm5sm2*soilPC2
	PC1 ~ fm6cm1*wetness + fm6cm2*rain_lowRH + fm6cm3*iso + fm6cm4*pseas + fm6cm5*mat + fm6cm6*radiation + fm6sm1*soilPC1 + fm6sm2*soilPC2

	# Regional climate, pollution, and regional forest heterogeneity effects on regional richness
	regS ~ R1CM1*wetness_reg_mean + R1CM2*rain_lowRH_reg_mean + R1CM3*iso_reg_mean + R1CM4*pseas_reg_mean + R1CM5*mat_reg_mean +
		R1CH1*wetness_reg_var + R1CH2*rain_lowRH_reg_var + R1CH3*iso_reg_var + R1CH4*pseas_reg_var + R1CH5*mat_reg_var +
		R1FH1*regS_tree 

	# Regional climate effects on regional forest heterogeneity
	regS_tree ~ FH1CM1*wetness_reg_mean + FH1CM2*rain_lowRH_reg_mean + FH1CM3*iso_reg_mean + FH1CM4*pseas_reg_mean + FH1CM5*mat_reg_mean +
		FH1CH1*wetness_reg_var + FH1CH2*rain_lowRH_reg_var + FH1CH3*iso_reg_var + FH1CH4*pseas_reg_var + FH1CH5*mat_reg_var

	# Effects on local lichen richness
	lichen.rich_log ~ r1fh1*bark_moist_pct.rao.ba + r1fh2*wood_SG.rao.ba + r1fh3*LogSeed.rao.ba +
		r1fh4*PIE.ba.tree + r1fh5*propDead + r1fh6*lightDist.mean + r1fh7*diamDiversity +
		r1fm1*bark_moist_pct.ba + r1fm2*wood_SG.ba + r1fm3*LogSeed.ba + r1fm4*bigTrees + r1fm5*light.mean + r1fm6*PC1 +
		r1cm1*wetness + r1cm2*rain_lowRH + r1cm3*iso + r1cm4*pseas + r1cm5*mat + r1cm6*radiation + r1sm1*soilPC1 + r1sm2*soilPC2 +
		+ r1R1*regS + r1a1*tot_abun_log
	
	# Effects on local lichen abundance
	tot_abun_log ~ a1fh1*bark_moist_pct.rao.ba + a1fh2*wood_SG.rao.ba + a1fh3*LogSeed.rao.ba +
		a1fh4*PIE.ba.tree + a1fh5*propDead + a1fh6*lightDist.mean + a1fh7*diamDiversity +
		a1fm1*bark_moist_pct.ba + a1fm2*wood_SG.ba + a1fm3*LogSeed.ba + a1fm4*bigTrees + a1fm5*light.mean + a1fm6*PC1 +
		a1cm1*wetness + a1cm2*rain_lowRH + a1cm3*iso + a1cm4*pseas + a1cm5*mat + a1cm6*radiation + 
		a1sm1*soilPC1 + a1sm2*soilPC2 

	## Covariances between exogenous and endogenous variables
	# Dont need to specificy for Climate L-R b/c these are exogenous and are calculated automatically in lavaan

	# Local-regional forest structure
	regS_tree ~~ FH1fh4*PIE.ba.tree

	## Covariances among endogenous predictors in the same group
	# Dont need to specificy for Climate LO b/c these are exogenous and are calculated automatically in lavaan
	
	# Forest LH
	PIE.ba.tree ~~ bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba + propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba ~~ propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba ~~ wood_SG.rao.ba + LogSeed.rao.ba
	wood_SG.rao.ba ~~ LogSeed.rao.ba
	propDead ~~ lightDist.mean + diamDiversity
	lightDist.mean ~~ diamDiversity

	# Forest LO
	bark_moist_pct.ba + wood_SG.ba + LogSeed.ba ~~ bigTrees + light.mean + PC1
	bark_moist_pct.ba ~~ wood_SG.ba + LogSeed.ba 
	wood_SG.ba ~~ LogSeed.ba 
	bigTrees ~~ light.mean + PC1
	light.mean ~~ PC1	
"

# SoilMod, RegToRich, AllSp
path_soilmod_regTorich_nopol = "


	# Local environment/climate effects on forest structure
	bark_moist_pct.rao.ba ~ fh1cm1*wetness + fh1cm2*rain_lowRH + fh1cm3*iso + fh1cm4*pseas + fh1cm5*mat + fh1cm6*radiation + fh1sm1*soilPC1 + fh1sm2*soilPC2
	wood_SG.rao.ba ~ fh2cm1*wetness + fh2cm2*rain_lowRH + fh2cm3*iso + fh2cm4*pseas + fh2cm5*mat + fh2cm6*radiation + fh2sm1*soilPC1 + fh2sm2*soilPC2
	LogSeed.rao.ba ~ fh3cm1*wetness + fh3cm2*rain_lowRH + fh3cm3*iso + fh3cm4*pseas + fh3cm5*mat + fh3cm6*radiation + fh3sm1*soilPC1 + fh3sm2*soilPC2
	PIE.ba.tree ~ fh4cm1*wetness + fh4cm2*rain_lowRH + fh4cm3*iso + fh4cm4*pseas + fh4cm5*mat + fh4cm6*radiation + fh4sm1*soilPC1 + fh4sm2*soilPC2
	propDead ~ fh5cm1*wetness + fh5cm2*rain_lowRH + fh5cm3*iso + fh5cm4*pseas + fh5cm5*mat + fh5cm6*radiation + fh5sm1*soilPC1 + fh5sm2*soilPC2
	lightDist.mean ~ fh6cm1*wetness + fh6cm2*rain_lowRH + fh6cm3*iso + fh6cm4*pseas + fh6cm5*mat + fh6cm6*radiation + fh6sm1*soilPC1 + fh6sm2*soilPC2
	diamDiversity ~ fh7cm1*wetness + fh7cm2*rain_lowRH + fh7cm3*iso + fh7cm4*pseas + fh7cm5*mat + fh7cm6*radiation + fh7sm1*soilPC1 + fh7sm2*soilPC2
	bark_moist_pct.ba ~ fm1cm1*wetness + fm1cm2*rain_lowRH + fm1cm3*iso + fm1cm4*pseas + fm1cm5*mat + fm1cm6*radiation + fm1sm1*soilPC1 + fm1sm2*soilPC2
	wood_SG.ba ~ fm2cm1*wetness + fm2cm2*rain_lowRH + fm2cm3*iso + fm2cm4*pseas + fm2cm5*mat + fm2cm6*radiation + fm2sm1*soilPC1 + fm2sm2*soilPC2
	LogSeed.ba ~ fm3cm1*wetness + fm3cm2*rain_lowRH + fm3cm3*iso + fm3cm4*pseas + fm3cm5*mat + fm3cm6*radiation + fm3sm1*soilPC1 + fm3sm2*soilPC2
	bigTrees ~ fm4cm1*wetness + fm4cm2*rain_lowRH + fm4cm3*iso + fm4cm4*pseas + fm4cm5*mat + fm4cm6*radiation + fm4sm1*soilPC1 + fm4sm2*soilPC2
	light.mean  ~ fm5cm1*wetness + fm5cm2*rain_lowRH + fm5cm3*iso + fm5cm4*pseas + fm5cm5*mat + fm5cm6*radiation + fm5sm1*soilPC1 + fm5sm2*soilPC2
	PC1 ~ fm6cm1*wetness + fm6cm2*rain_lowRH + fm6cm3*iso + fm6cm4*pseas + fm6cm5*mat + fm6cm6*radiation + fm6sm1*soilPC1 + fm6sm2*soilPC2

	# Regional climate, pollution, and regional forest heterogeneity effects on regional richness
	regS ~ R1CM1*wetness_reg_mean + R1CM2*rain_lowRH_reg_mean + R1CM3*iso_reg_mean + R1CM4*pseas_reg_mean + R1CM5*mat_reg_mean +
		R1CH1*wetness_reg_var + R1CH2*rain_lowRH_reg_var + R1CH3*iso_reg_var + R1CH4*pseas_reg_var + R1CH5*mat_reg_var +
		R1FH1*regS_tree 

	# Regional climate effects on regional forest heterogeneity
	regS_tree ~ FH1CM1*wetness_reg_mean + FH1CM2*rain_lowRH_reg_mean + FH1CM3*iso_reg_mean + FH1CM4*pseas_reg_mean + FH1CM5*mat_reg_mean +
		FH1CH1*wetness_reg_var + FH1CH2*rain_lowRH_reg_var + FH1CH3*iso_reg_var + FH1CH4*pseas_reg_var + FH1CH5*mat_reg_var

	# Effects on local lichen richness
	lichen.rich_log ~ r1fh1*bark_moist_pct.rao.ba + r1fh2*wood_SG.rao.ba + r1fh3*LogSeed.rao.ba +
		r1fh4*PIE.ba.tree + r1fh5*propDead + r1fh6*lightDist.mean + r1fh7*diamDiversity +
		r1fm1*bark_moist_pct.ba + r1fm2*wood_SG.ba + r1fm3*LogSeed.ba + r1fm4*bigTrees + r1fm5*light.mean + r1fm6*PC1 +
		r1cm1*wetness + r1cm2*rain_lowRH + r1cm3*iso + r1cm4*pseas + r1cm5*mat + r1cm6*radiation + r1sm1*soilPC1 + r1sm2*soilPC2 +
		r1R1*regS + r1a1*tot_abun_log + r1CM1*wetness_reg_mean + r1CM3*rain_lowRH_reg_mean + r1CM3*iso_reg_mean + r1CM4*pseas_reg_mean + r1CM5*mat_reg_mean +
		r1CH1*wetness_reg_var + r1CH2*rain_lowRH_reg_var + r1CH3*iso_reg_var + r1CH4*pseas_reg_var + r1CH5*mat_reg_var +
		r1FH1*regS_tree
	
	# Effects on local lichen abundance
	tot_abun_log ~ a1fh1*bark_moist_pct.rao.ba + a1fh2*wood_SG.rao.ba + a1fh3*LogSeed.rao.ba +
		a1fh4*PIE.ba.tree + a1fh5*propDead + a1fh6*lightDist.mean + a1fh7*diamDiversity +
		a1fm1*bark_moist_pct.ba + a1fm2*wood_SG.ba + a1fm3*LogSeed.ba + a1fm4*bigTrees + a1fm5*light.mean + a1fm6*PC1 +
		a1cm1*wetness + a1cm2*rain_lowRH + a1cm3*iso + a1cm4*pseas + a1cm5*mat + a1cm6*radiation + 
		a1sm1*soilPC1 + a1sm2*soilPC2 

	## Covariances between exogenous and endogenous variables
	# Dont need to specificy for Climate L-R b/c these are exogenous and are calculated automatically in lavaan

	# Local-regional forest structure
	regS_tree ~~ FH1fh4*PIE.ba.tree

	## Covariances among endogenous predictors in the same group
	# Dont need to specificy for Climate LO b/c these are exogenous and are calculated automatically in lavaan
	
	# Forest LH
	PIE.ba.tree ~~ bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba + propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba ~~ propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba ~~ wood_SG.rao.ba + LogSeed.rao.ba
	wood_SG.rao.ba ~~ LogSeed.rao.ba
	propDead ~~ lightDist.mean + diamDiversity
	lightDist.mean ~~ diamDiversity

	# Forest LO
	bark_moist_pct.ba + wood_SG.ba + LogSeed.ba ~~ bigTrees + light.mean + PC1
	bark_moist_pct.ba ~~ wood_SG.ba + LogSeed.ba 
	wood_SG.ba ~~ LogSeed.ba 
	bigTrees ~~ light.mean + PC1
	light.mean ~~ PC1	
"


#############################################################################
## Models with pollution

# BaseMod, AllSp
path_pol = "

	# Local environment/climate effects on forest structure
	bark_moist_pct.rao.ba ~ fh1cm1*wetness + fh1cm2*rain_lowRH + fh1cm3*iso + fh1cm4*pseas + fh1cm5*mat + fh1cm6*radiation
	wood_SG.rao.ba ~ fh2cm1*wetness + fh2cm2*rain_lowRH + fh2cm3*iso + fh2cm4*pseas + fh2cm5*mat + fh2cm6*radiation
	LogSeed.rao.ba ~ fh3cm1*wetness + fh3cm2*rain_lowRH + fh3cm3*iso + fh3cm4*pseas + fh3cm5*mat + fh3cm6*radiation
	PIE.ba.tree ~ fh4cm1*wetness + fh4cm2*rain_lowRH + fh4cm3*iso + fh4cm4*pseas + fh4cm5*mat + fh4cm6*radiation
	propDead ~ fh5cm1*wetness + fh5cm2*rain_lowRH + fh5cm3*iso + fh5cm4*pseas + fh5cm5*mat + fh5cm6*radiation
	lightDist.mean ~ fh6cm1*wetness + fh6cm2*rain_lowRH + fh6cm3*iso + fh6cm4*pseas + fh6cm5*mat + fh6cm6*radiation
	diamDiversity ~ fh7cm1*wetness + fh7cm2*rain_lowRH + fh7cm3*iso + fh7cm4*pseas + fh7cm5*mat + fh7cm6*radiation
	bark_moist_pct.ba ~ fm1cm1*wetness + fm1cm2*rain_lowRH + fm1cm3*iso + fm1cm4*pseas + fm1cm5*mat + fm1cm6*radiation
	wood_SG.ba ~ fm2cm1*wetness + fm2cm2*rain_lowRH + fm2cm3*iso + fm2cm4*pseas + fm2cm5*mat + fm2cm6*radiation
	LogSeed.ba ~ fm3cm1*wetness + fm3cm2*rain_lowRH + fm3cm3*iso + fm3cm4*pseas + fm3cm5*mat + fm3cm6*radiation
	bigTrees ~ fm4cm1*wetness + fm4cm2*rain_lowRH + fm4cm3*iso + fm4cm4*pseas + fm4cm5*mat + fm4cm6*radiation
	light.mean  ~ fm5cm1*wetness + fm5cm2*rain_lowRH + fm5cm3*iso + fm5cm4*pseas + fm5cm5*mat + fm5cm6*radiation
	PC1 ~ fm6cm1*wetness + fm6cm2*rain_lowRH + fm6cm3*iso + fm6cm4*pseas + fm6cm5*mat + fm6cm6*radiation

	# Local climate effects on pollution
	totalNS ~ p1cm1*wetness + p1cm2*rain_lowRH

	# Regional climate effects on pollution
	totalNS_reg ~ P1CM1*wetness_reg_mean + P1CM2*rain_lowRH_reg_mean

	# Regional climate, pollution, and regional forest heterogeneity effects on regional richness
	regS ~ R1CM1*wetness_reg_mean + R1CM2*rain_lowRH_reg_mean + R1CM3*iso_reg_mean + R1CM4*pseas_reg_mean + R1CM5*mat_reg_mean +
		R1CH1*wetness_reg_var + R1CH2*rain_lowRH_reg_var + R1CH3*iso_reg_var + R1CH4*pseas_reg_var + R1CH5*mat_reg_var +
		R1FH1*regS_tree + R1P1*totalNS_reg

	# Regional climate effects on regional forest heterogeneity
	regS_tree ~ FH1CM1*wetness_reg_mean + FH1CM2*rain_lowRH_reg_mean + FH1CM3*iso_reg_mean + FH1CM4*pseas_reg_mean + FH1CM5*mat_reg_mean +
		FH1CH1*wetness_reg_var + FH1CH2*rain_lowRH_reg_var + FH1CH3*iso_reg_var + FH1CH4*pseas_reg_var + FH1CH5*mat_reg_var

	# Effects on local lichen richness
	lichen.rich_log ~ r1fh1*bark_moist_pct.rao.ba + r1fh2*wood_SG.rao.ba + r1fh3*LogSeed.rao.ba +
		r1fh4*PIE.ba.tree + r1fh5*propDead + r1fh6*lightDist.mean + r1fh7*diamDiversity +
		r1fm1*bark_moist_pct.ba + r1fm2*wood_SG.ba + r1fm3*LogSeed.ba + r1fm4*bigTrees + r1fm5*light.mean + r1fm6*PC1 +
		r1cm1*wetness + r1cm2*rain_lowRH + r1cm3*iso + r1cm4*pseas + r1cm5*mat + r1cm6*radiation + r1R1*regS + r1p1*totalNS+ r1a1*tot_abun_log
	
	# Effects on local lichen abundance
	tot_abun_log ~ a1fh1*bark_moist_pct.rao.ba + a1fh2*wood_SG.rao.ba + a1fh3*LogSeed.rao.ba +
		a1fh4*PIE.ba.tree + a1fh5*propDead + a1fh6*lightDist.mean + a1fh7*diamDiversity +
		a1fm1*bark_moist_pct.ba + a1fm2*wood_SG.ba + a1fm3*LogSeed.ba + a1fm4*bigTrees + a1fm5*light.mean + a1fm6*PC1 +
		a1cm1*wetness + a1cm2*rain_lowRH + a1cm3*iso + a1cm4*pseas + a1cm5*mat + a1cm6*radiation + a1p1*totalNS

	## Covariances between exogenous and endogenous variables
	# Dont need to specificy for Climate L-R b/c these are exogenous and are calculated automatically in lavaan
	
	# Local-regional pollution
	totalNS ~~ p1P1*totalNS_reg

	# Local-regional forest structure
	regS_tree ~~ FH1fh4*PIE.ba.tree

	## Covariances among endogenous predictors in the same group
	# Dont need to specificy for Climate LO b/c these are exogenous and are calculated automatically in lavaan
	
	# Forest LH
	PIE.ba.tree ~~ bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba + propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba + wood_SG.rao.ba + LogSeed.rao.ba ~~ propDead + lightDist.mean + diamDiversity
	bark_moist_pct.rao.ba ~~ wood_SG.rao.ba + LogSeed.rao.ba
	wood_SG.rao.ba ~~ LogSeed.rao.ba
	propDead ~~ lightDist.mean + diamDiversity
	lightDist.mean ~~ diamDiversity

	# Forest LO
	bark_moist_pct.ba + wood_SG.ba + LogSeed.ba ~~ bigTrees + light.mean + PC1
	bark_moist_pct.ba ~~ wood_SG.ba + LogSeed.ba 
	wood_SG.ba ~~ LogSeed.ba 
	bigTrees ~~ light.mean + PC1
	light.mean ~~ PC1
"


