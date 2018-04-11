# MDPBiome init file.  This load the required library and function files
# Please set the wworking directory to the Src path

library(phyloseq)
library(Hmisc)
library(MDPtoolbox)
library(reshape2)
library(ggplot2)


sourceFiles <- c("support.R",
                 "mdpPreprocess.R",
                 "mdpcreation.R",
                 "mdpOutput.R",
                 "mdpEvaluation.R",
                 "mdpBiomBase.R",
                 "mdpgraphs.R")
sapply(sourceFiles,source,.GlobalEnv)