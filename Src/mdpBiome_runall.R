# Run MDPBiome for all datasets 

source("initMDPBiome.R")

# Running Ballou2016
source("../Data/Ballou2016/config_Ballou2016.R")
#source("../Data/Ballou2016_Dominant0.5perc/config_Ballou2016_dominant0.5perc.R")
#source("../Data/Ballou2016_Dominant1perc/config_Ballou2016_dominant1perc.R")
#source("../Data/Ballou2016_NonDominant0.5perc/config_Ballou2016_nondominant0.5perc.R")
#source("../Data/Ballou2016_NonDominant1perc/config_Ballou2016_nondominant1perc.R")

mdpBiomeBase()
mdpBiomeLoocv("preferGood",goalVar = "diversity")

#============================================================================================

# Running LaRosa 2014 ::: Comment out one config file to execute 
source("../Data/LaRosa2014_S/config_LaRosa2014_S.R")
#source("../Data/LaRosa2014_Dominant0.5perc/config_LaRosa2014_dominant0.5perc.R")
#source("../Data/LaRosa2014_Dominant1perc/config_LaRosa2014_dominant1perc.R")
#source("../Data/LaRosa2014_NonDominant0.5perc/config_LaRosa2014_nondominant0.5perc.R")
#source("../Data/LaRosa2014_NonDominant1perc/config_LaRosa2014_nondominant1perc.R")


mdpBiomePreAnalysis(data.norm,data.raw)
mdpBiomeBase()
mdpBiomeLoocv("preferGood",goalVar = "diversity")

#============================================================================================

# Run Gajer2012
source("../Data/Gajer2012/config_Gajer2012.R")
#source("../Data/Gajer2012_AllTaxa/config_Gajer2012_AllTaxa.R")
#source("../Data/Gajer2012_Dominant0.5perc/config_Gajer2012_dominant0.5perc.R")
#source("../Data/Gajer2012_Dominant1perc/config_Gajer2012_dominant1perc.R")
# No raw data to compute diversity, niether taxa at phylum level
# mdpBiomePreAnalysis(data.norm)
mdpBiomeBase(goalDiversity = FALSE)
mdpBiomeLoocv(goal_preference,"avoidBad",goalVar = "Inv.Nugent")

# TO RUN OVERNIGHT
# HEREEEEEEEEE
source("../Data/Gajer2012_AllTaxa/config_Gajer2012_AllTaxa.R")
mdpBiomeBase(goalDiversity = FALSE)
mdpBiomeLoocv(goal_preference,"avoidBad",goalVar = "Inv.Nugent")

source("../Data/Gajer2012_Dominant0.5perc/config_Gajer2012_dominant0.5perc.R")
mdpBiomeBase(goalDiversity = FALSE)
mdpBiomeLoocv(goal_preference,"avoidBad",goalVar = "Inv.Nugent")

source("../Data/Gajer2012_Dominant1perc/config_Gajer2012_dominant1perc.R")
mdpBiomeBase(goalDiversity = FALSE)
mdpBiomeLoocv(goal_preference,"avoidBad",goalVar = "Inv.Nugent")

