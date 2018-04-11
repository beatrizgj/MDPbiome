
###############################################################
# MAIN PROGRAM (local calls to several datasets at once)
###############################################################
  
library(ggplot2,quietly=TRUE)
library(phyloseq)
# Warnings with ggplot2 before phyloseq (I hope it is not a problem later)
#Warning messages:
#  1: replacing previous import ‘BiocGenerics::Position’ by ‘ggplot2::Position’ when loading ‘phyloseq’ 
#  2: replacing previous import ‘S4Vectors::Position’ by ‘ggplot2::Position’ when loading ‘DESeq2’
library(lattice,quietly=TRUE)
library(cluster,quietly=TRUE)
library(fpc,quietly=TRUE)
library(reshape2,quietly=TRUE) # Load the reshape2 package (for the melt() function)

source('robust.clustering.metagenomics.functions.r')


#######################################
# OPTION 1: Only one call per datasets
#######################################
robust.clustering.all.steps("~/RobustClustering/David2014","~/RobustClustering/David2014/data.norm_David2014.RData",'David2014',"COLLECTION_DAY")

robust.clustering.all.steps("~/RobustClustering/Ballou2016","~/RobustClustering/Ballou2016/data.norm_Ballou2016.RData",'Ballou2016',"age")

robust.clustering.all.steps("~/RobustClustering/Gajer2012","~/RobustClustering/Gajer2012/data.norm_Gajer2012.RData",'Gajer2012',"Time.in.study")

robust.clustering.all.steps("~/RobustClustering/LaRosa2014","~/RobustClustering/LaRosa2014/data.norm_LaRosa2014.RData",'LaRosa2014',"age")


# Calls for taxa subset, with different dominant percentage, with Ballou2016 dataset:
robust.clustering.all.steps("Dominant0.5perc","data.norm_Ballou2016.RData",'Ballou2016',"age","dominant",0.5)
robust.clustering.all.steps("Dominant0.5perc","data.norm_Ballou2016.RData",'Ballou2016',"age","nonDominant",0.5)
robust.clustering.all.steps("Dominant0.5perc","data.norm_Gajer2012.RData",'Gajer2012',"Time.in.study","dominant",0.5)
robust.clustering.all.steps("Dominant0.5perc","data.norm_Gajer2012.RData",'Gajer2012',"Time.in.study","nonDominant",0.5)
robust.clustering.all.steps("Dominant0.5perc","data.norm_LaRosa2014.RData",'LaRosa2014',"age","dominant",0.5)
robust.clustering.all.steps("Dominant0.5perc","data.norm_LaRosa2014.RData",'LaRosa2014',"age","nonDominant",0.5)
robust.clustering.all.steps("Dominant0.5perc","data.norm_David2014.RData",'David2014',"COLLECTION_DAY","dominant",0.5)
robust.clustering.all.steps("Dominant0.5perc","data.norm_David2014.RData",'David2014',"COLLECTION_DAY","nonDominant",0.5)


# Calls for taxa subset, with David2014 dataset:
# Subset dominant taxa:
robust.clustering.all.steps("~/RobustClustering/David2014","~/RobustClustering/David2014/data.norm_David2014.RData",'David2014',"COLLECTION_DAY",'dominant',1)
# Subset nonDominant taxa:
robust.clustering.all.steps("~/RobustClustering/David2014","~/RobustClustering/David2014/data.norm_David2014.RData",'David2014',"COLLECTION_DAY",'nonDominant')
# Subset all taxa at genes level:
robust.clustering.all.steps("~/RobustClustering/David2014","~/RobustClustering/David2014/data.norm_David2014.RData",'David2014',"COLLECTION_DAY",'all','genus')
# Subset dominant taxa at genes level:
robust.clustering.all.steps("~/RobustClustering/David2014","~/RobustClustering/David2014/data.norm_David2014.RData",'David2014',"COLLECTION_DAY",'dominant','genus')
# Subset nonDominant taxa at genes level:
robust.clustering.all.steps("~/RobustClustering/David2014","~/RobustClustering/David2014/data.norm_David2014.RData",'David2014',"COLLECTION_DAY",'nonDominant','genus')


###########################################################
# OPTION 2: Call to the independent functions, step by step
###########################################################
# A) ROBUST CLUSTERING COMPUTATION
setwd("RobustClustering/David2014")
load('data.norm_David2014.RData')
list.results <- robust.clustering(data.norm,"COLLECTION_DAY","David2014")
rm(data.norm)

setwd("RobustClustering/Ballou2016")
load('data.norm_Ballou2016.RData')
list.results <- robust.clustering(data.norm,"age","Ballou20106")
rm(data.norm)

setwd("RobustClustering/Gajer2012")
load('data.norm_Gajer2012.RData')
list.results <- robust.clustering(data.norm,"Time.in.study","Gajer2012")
rm(data.norm)

setwd("RobustClustering/LaRosa2014")
load('data.norm_LaRosa2014.RData')
list.results <- robust.clustering(data.norm,"age","LaRosa2014")
rm(data.norm)

########################################
# B) PLOT GRAPHS
# Using next loop requires you are in a directory where there are a folder named as the dataset (ej. 'RobustClustering/' with subfolders 'David2014', 'Ballou2016', etc). Each subfolder should have at least the RData with the phyloseq object called data.norm (in these examples <data.norm_XXdatasetXX.RData>).
for (dataset in c('David2014','Ballou2016','Gajer2012','LaRosa2014')){
  setwd(dataset)
  load(paste('eval.arrays_',dataset,'.RData',sep=''))
  plot.robust.clustering.all.together.formatted(eval.array4d, eval.SIpoints, eval.Kvalues)
  setwd('..')
} # end-for


########################################
# C) AUTOMATIC CLUSTERING DECISION
for (dataset in c('David2014','Ballou2016','Gajer2012','LaRosa2014')){
  setwd(dataset)
  load(paste('data.norm_',dataset,'.RData',sep=''))
  load(paste('eval.arrays_',dataset,'.RData',sep=''))
  switch(dataset,
         David2014={ var.graph="COLLECTION_DAY" },
         Chicks={ var.graph="age" },
         Gajer2012={ var.graph="Time.in.study" },
         LaRosa2014={ var.graph="age" },
         { var.graph=NULL;
           print('ERROR: dataset without var.color defined!!')}
  ) #end-switch-alg
  data.norm <- robust.clustering.decision(data.norm,eval.array4d,var.graph)
    rm(data.norm)
  setwd('..')
} # end-for
