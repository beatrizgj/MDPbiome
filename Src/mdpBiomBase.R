# Main Call to MDPBiome

# General call for MDPBiome
#   goalDiversity: TRUE if reward table is based on computed alpha diversity. When FALSE the user should
#                  provide cluster preferences
#   utilityVar   : Attribute containing the "goodness" of samples. The per-cluster average will be used as
#                  utility function

mdpBiomeBase <- function(goalDiversity=TRUE, utilityVar=NULL){ 

  mdpstates <- levels(sample_data(data.norm)$cluster)
  num.states <- length(mdpstates)
  
  genResult <- data.frame(Perturbation=character(),Reward=character(),Stability=numeric(),
                          stringsAsFactors = FALSE)
  
  if (goalDiversity){
    alphaDiversity <- computeDiversity(data.raw, data.norm)
    writeStateTable(alphaDiversity, mdpstates, "alpha_diversity", diroutput)
    cluster_preference <- alphaDiversity
  }
  
  for (iPerturbation in Perturbations){
    print(paste ("Processing perturbation ", iPerturbation, "..."))
    sample_data(data.norm)$action <- as.factor(get_variable(data.norm,iPerturbation))
    
    transObserved <- transitionTriplets(data.norm, iPerturbation, 
                                        timeStepId = stepVar, subjectId=subjectVar)
    
    mdpactions <- levels(transObserved$Action)
    
    num.actions <- length(mdpactions)
    mdpdim.names <- list(mdpstates, mdpstates, mdpactions)
    
    transCounts <- transitionCounts(transObserved, mdpdim.names)
    writeMDPTable(transCounts,"table_counts","all",iPerturbation,diroutput, saveObject = TRUE)
    
    QendPrior <- endStatePrior(transCounts)
    print("---Prior---")
    print(QendPrior)
    transFrequencies <- transitionFrequencies(transCounts, mdpdim.names)
    writeMDPTable(transFrequencies,"table_basefrequencies","all",iPerturbation,diroutput)
    
    #fill in all the transitions with prob=0 for the same initial state, and whatever action.
    mdpTransitions <- setEmptyTransitions(transFrequencies, transCounts)
    writeMDPTable(mdpTransitions,"table_transitions","all",iPerturbation,diroutput, saveObject = TRUE)
    
    mdpGraph(mdpTransitions, iPerturbation, diroutput)
    mdpChrodGraph(mdpTransitions, iPerturbation, diroutput)
    
    # Iterate for the different reward function classes
    for(jreward in c("preferGood","avoidBad","proportional")){
      rewardsFunction <-rewardTable(mdpstates, mdpactions, cluster_preference, jreward)  
      writeMDPTable(rewardsFunction,"table_reward",jreward,iPerturbation,diroutput)
      
      mdpbiome.policy <- mdp_value_iteration(mdpTransitions, rewardsFunction, 0.9)
      
      print("****OPTIMAL POLICY****")
      print(paste("  Actions:", paste(mdpactions, collapse = " ")))
      print(paste("  Policy:", paste(mdpbiome.policy$policy, collapse = " ")))
      
      writeMDPsolution(mdpbiome.policy,iPerturbation,mdpstates,mdpactions,jreward,diroutput)
      
      res <- evalPerturbation(transCounts, QendPrior, mdpdim.names, iPerturbation, rewardsFunction, jreward,
                              diroutput, nsamples = 1000)
      print("")
      print("============================================")
      print(paste("Schema:", jreward, " -- Results for Perturbation:",iPerturbation))
      print(res)
      print("")
      # the second elements is the stability ratio for actions
      genResult[(nrow(genResult)+1),] <- c(iPerturbation,jreward,res[2])
      
    }# end-for ireward
  }# end-for iPerturbation
  
  # General result as table
  sink(paste(diroutput,"genResult.tsv",sep = ""))
  print(dcast(genResult, Perturbation ~ Reward))
  sink()
  
  return(genResult)
}


# Leave one out evaluation for MDPBBiome
# One subject is left out to verify if the policy is consistent to what the out of sample subject is doing.

#   goalDiversity: TRUE if reward table is based on computed alpha diversity. When FALSE the user should
#                  provide cluster preferences

mdpBiomeLoocv <- function(rewardType, goalVar=NULL, goal_preference=NULL, loadPolicies=FALSE){ 

  mdpstates <- levels(sample_data(data.norm)$cluster)
  num.states <- length(mdpstates)
  
  #genResult <- data.frame(Perturbation=character(),Reward=character(),Stability=numeric(),
  #                        stringsAsFactors = FALSE)
  
  if (goalVar == "diversity"){
    alphaDiversity <- computeDiversity(data.raw, data.norm)
    writeStateTable(alphaDiversity, mdpstates, "alpha_diversity", diroutput_loocv)
    loocv_preference <-  normalize(alphaDiversity)
    cluster_preference <- loocv_preference
    sample_data(data.norm)$diversity <- estimate_richness(data.raw,measures = "Shannon")$Shannon
  }else{
    loocv_preference <- goal_preference
  }
  
  list.subjects <- unique(sort(get_variable(data.norm, subjectVar)))
  num.subjects <- length(list.subjects)
  
  followvars <- c("follow-better","follow-equal","follow-worse", 
                  "notfollow-better","notfollow-equal","notfollow-worse")
  overallFollowFreq <- array(dim=c(length(Perturbations),6), dimnames = list(Perturbations, followvars)) 
  overallFollowCounts <- data.frame(Perturbation=character(),
                                   Follow=character(),
                                   Better=numeric(), Equal=numeric(),Worse=numeric(),
                                   stringsAsFactors = FALSE)  
  
  for (iPerturbation in Perturbations){
    print(paste ("Processing perturbation ", iPerturbation, "..."))
  
    loocvResults <- array(dim= c(num.subjects,4), 
                          dimnames = list(list.subjects, c("policy","Avg.Utility",
                                                           "Follow.Utility","Not-follow.Utility"))) 
    loocvFollowCounts <- array(dim=c(num.subjects,6), dimnames = list(list.subjects, followvars)) 
                                
    # This takes too long. We can reload from file if nothing else has changed
    if (loadPolicies){
      filePolicies <- paste(diroutput_loocv,iPerturbation,"/policies_loocv-all.txt", sep="")
      policies <- as.matrix(read.table(filePolicies))
    }else{
      policies <- loocvPolicies(data.norm, rewardType, cluster_preference, iPerturbation, diroutput_loocv)
    }
    
    for (jsubject in list.subjects){
      
      percFollowed <- subjectFollowedPolicy(data.norm, jsubject, policies[jsubject,], iPerturbation)
      percInGoal <- subjectInGoal(data.norm, jsubject, loocv_preference ,iPerturbation)
      deltaGoal <- subjectDeltaGoal(data.norm, jsubject, policies[jsubject,], goalVar, iPerturbation)
      followCounts <- subjectFollowProgress(data.norm, jsubject, policies[jsubject,], 
                                         loocv_preference, iPerturbation)
      
      loocvResults[jsubject,] <- c(percFollowed, deltaGoal)
      loocvFollowCounts[jsubject,] <- followCounts
    }
    
    print(paste(" LOOCV RESULTS FOR ",iPerturbation))
    print(loocvFollowCounts)
    writeMDPTable(loocvFollowCounts, "loocv_follow_result",rewardType,iPerturbation, diroutput_loocv)
    print("--------------------------------------------")
    totalFollow <- colSums(loocvFollowCounts)
    print(totalFollow)
    followpercent <- totalFollow[1:3]/sum(totalFollow[1:3])
    notfollowpercent <- totalFollow[4:6]/sum(totalFollow[4:6])
    
    overallFollowFreq[iPerturbation,] <- c(followpercent, notfollowpercent)
    
    n <- nrow(overallFollowCounts)
    overallFollowCounts[n+1,1:2] <- c(iPerturbation,"F") 
    overallFollowCounts[n+1,3:5] <- totalFollow[1:3]
    
    overallFollowCounts[n+2,1:2] <- c(iPerturbation,"nF") 
    overallFollowCounts[n+2,3:5] <- totalFollow[4:6]

    #overallFollowCounts[n+2,] <- c(iPerturbation,"no",totalFollow[4:6])
    
    
    #followPolicyGraph(as.data.frame(loocvResults),iPerturbation,diroutput_loocv)
  }# end-for iPerturbation
  
  print(" OVERALL LOOCV RESULTS")
  writeMDPTable(overallFollowFreq, "loocv_follow_overall_freq_result",rewardType,"", diroutput_loocv)
  writeMDPTable(overallFollowCounts, "loocv_follow_overall_count_result",rewardType,"", diroutput_loocv)
  followProb <- summarizeFollow(overallFollowCounts,eqType = 0)

  followBarsGraph(followProb, diroutput_loocv)
  
  return(overallFollowCounts)
}


mdpBiomePreAnalysis <- function(normData, rawData=NULL){
  if (!missing(rawData)){
    
    sample_data(rawData)$cluster <- get_variable(data.norm,"cluster")
    
    # Alpha diversity by State 
    pdf(paste(diranalysis,"alphaDiversityByState.pdf",sep=""),width=10)
    print(plot_richness(rawData, x="cluster"))
    dev.off()
    png(paste(diranalysis,"alphaDiversityByState.png",sep=""))
    print(plot_richness(rawData, x="cluster"))
    dev.off()
    
  }
  
  pdf(paste(diranalysis,"barplot_cluster_phylum_merge.pdf",sep=""),width=10)
  pg2<-plot_bar(data.norm, "cluster", fill="Phylum")
  print(pg2 + geom_bar(aes(color=Phylum, fill=Phylum, order=factor(Phylum)), 
                       stat="identity", position="stack"))
  dev.off()
  
  
}


