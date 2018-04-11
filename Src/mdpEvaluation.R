# Functions for MDP Evaluation
library(DirichletReg)
library(MDPtoolbox)

# This function sample a transition table following
# a Dirichlet distribution for the p(Q.end | Qini,Action)
# The parameters are the observed counts from the original data
#   basecounts:  Array [Ini, End, Action] with the observed occurrrences
#   priorConcentration: vector with the frequency of reaching end state.
#                       this will be used as the initial concentration parameters to avoid 0 counts
sampleTransitionFun <- function(baseCounts, priorConcentration){
  
  states <- dimnames(baseCounts)[[1]]
  actions <- dimnames(baseCounts)[[3]]
  
  sampledTransitions <- array(data = 0, dim = c(length(states), length(states), length(actions)),
                                                dimnames = dimnames(baseCounts))
  
  for(istate in states){
    for(jaction in actions){
      ijcount <- baseCounts[istate, ,jaction] + priorConcentration 
      sampledTransitions[istate, ,jaction] <- rdirichlet(1,ijcount)
    }
  }
  return(sampledTransitions)
}

# Compute the difference between the optimal policy and the policy
# sampled from the evaluation schema, in terms of the evaluated metrics
#    varNames: data.frame colnames to return the result in the same format 
comparePolicies <- function(bestsol, samplesol, varNames){
  
  areEqual <- TRUE
  actionsMatch <- 0
  nstates <- length(bestsol$policy)
  individualMatch <- rep(0,nstates) 
  
  for (i in 1:nstates){
    if (bestsol$policy[i] == samplesol$policy[i]){
      actionsMatch <- actionsMatch + 1
      individualMatch[i] <- individualMatch[i] + 1
    }else{
      areEqual <- FALSE
    }
  }
  
  sampleres <- data.frame(as.list(c(areEqual, actionsMatch, individualMatch)))
  colnames(sampleres) <- varNames
  
  #print("Comparing...")
  #print(bestsol$policy)
  #print(samplesol$policy)
  #print(sampleres)
  
  return(sampleres)
}


evalPerturbation <- function(transCounts, priorConcentration, mdpdim.names, actionVar, rewardfn, rewardType,
                             diroutput, nsamples = 10){

  states <- mdpdim.names[[1]]
  actions <- mdpdim.names[[3]]
  num.actions <- length(actions)
  num.states <- length(states)
  
  transFrequencies <- transitionFrequencies(transCounts, mdpdim.names)
  
  #fill in all the transitions with prob=0 for the same initial state, and whatever action.
  transitionfn <- setEmptyTransitions(transFrequencies, transCounts)
 
  optimalPolicy <- mdp_value_iteration(transitionfn, rewardfn, 0.9)
  write.table(optimalPolicy$policy, file=paste(diroutput,actionVar,"/policy_raw_",rewardType,".tsv",sep=""))
  
  print("****OPTIMAL EVALUATION POLICY****")
  print(paste("  Actions:", paste(actions, collapse = " ")))
  print(paste("  Policy:", paste(optimalPolicy$policy, collapse = " ")))
  
  
  sampleResults <- data.frame(policymatch=logical(),
                              nactionsmatch=integer(),
                              lapply(states, function(x){numeric()}))
  colnames(sampleResults) <- c("policymatch","nactionsmatch",
                               sapply(states,function(s){paste("state.",s,sep="")}))
  
  policyCounts <- array(data=0, dim=c(num.states, num.actions), dimnames = list(states, actions))
  
  for (i in 1:nsamples){
    itransfn <- sampleTransitionFun(transCounts, priorConcentration)

    sink("/dev/null")
    ipolicy <- mdp_value_iteration(itransfn, rewardfn, 0.9)
    sink()
  
    for(j in 1:num.states){
      policyCounts[states[j],ipolicy$policy[j]] <- policyCounts[states[j],ipolicy$policy[j]]+1   
    } 
    
    diffpolicy <- comparePolicies(optimalPolicy, ipolicy, colnames(sampleResults))
    sampleResults <- rbind(sampleResults, diffpolicy)
  }
 
  print("---*-SAMPLE RESULTS-*----")
  #print(sampleResults)
  policyCounts <- policyCounts / nsamples
  print(policyCounts)

  freqMatched <- sum(sampleResults$policymatch)/nsamples
  freqActions <- sum(sampleResults$nactionsmatch)/(num.states*nsamples)
  
  otext <- paste("Perturbation  :", actionVar)
  otext <- paste("RewardType    :", rewardType)
  otext <- c(otext, paste("Samples       :",nsamples))
  otext <- c(otext, paste("% Matched     :", freqMatched))
  otext <- c(otext, paste("% Actions_equal:", freqActions))
  otext <- c(otext, "---------------------")
  
  freqByAction <- numeric()
  for (i in 1:length(states)){
   #the per-action matching start at third columns 
    iMatchedFreq <- sum(sampleResults[[i+2]])/nsamples
    otext <- c(otext, paste("% Matched",states[i],":", iMatchedFreq))
    freqByAction <- c(freqByAction, iMatchedFreq)
  }
  
  writeMDPTable(policyCounts,"table_policy_actionfreq",rewardType,actionVar, diroutput)
  actionEvalGraph(policyCounts, optimalPolicy$policy, rewardType, actionVar, diroutput)
    
  actiondir <- paste(diroutput,actionVar,"/",sep="")
  setdirectory(actiondir,paste(actionVar,"output"))
  resultFileName <- paste(actiondir, "sampling_evaluation",rewardType,".txt", sep="")
  ostream <- file(resultFileName, open = "w")
  writeLines(otext, ostream)
  close(ostream)
  
  return(c(freqMatched, freqActions, freqByAction))
  
} 


# Computes optimal policy for each subject, excluding samples tha belongs to the subject
loocvPolicies <- function(phyloObject, rewardType, clusterPref, iPerturbation, dirouteval){
  mdpstates <- levels(sample_data(phyloObject)$cluster)
  num.states <- length(mdpstates)

  print(paste ("Processing perturbation ", iPerturbation, "..."))
  sample_data(phyloObject)$action <- as.factor(get_variable(phyloObject,iPerturbation))
  mdpactions <- levels(sample_data(phyloObject)$action)
    
  num.actions <- length(mdpactions)
  mdpdim.names <- list(mdpstates, mdpstates, mdpactions)
    
  rewardsFn <-rewardTable(mdpstates, mdpactions, clusterPref, rewardType)
  writeMDPTable(rewardsFn,"table_reward",rewardType,iPerturbation,dirouteval)
    
  list.subjects <- unique(sort(get_variable(phyloObject, subjectVar)))
  num.subjects <- length(list.subjects)
    
  policy.matrix <- array(dim = c(num.subjects, num.states),  dimnames = list(list.subjects, mdpstates))
    
  for (isubject in list.subjects){
    print(paste("Processing MDP Leave-1-out evaluation: Excluding subject",isubject))
    trainData <- phyloFilterOut(phyloObject,subjectVar,isubject)
      
    transObserved <- transitionTriplets(trainData, iPerturbation, 
                                        timeStepId = stepVar, subjectId=subjectVar, 
                                        verbose = FALSE, saveResult = FALSE)
  
    transCounts <- transitionCounts(transObserved, mdpdim.names)
    #writeMDPTable(transCounts,"table_counts",isubject,iPerturbation,dirouteval)
      
    transFrequencies <- transitionFrequencies(transCounts, mdpdim.names)
    #writeMDPTable(transFrequencies,"table_basefrequencies","all",iPerturbation,diroutput)
      
    mdpTransitions <- setEmptyTransitions(transFrequencies, transCounts)
    #writeMDPTable(mdpTransitions,"table_transitions",isubject,iPerturbation,dirouteval)
      
    ipolicy.best <- mdp_value_iteration(mdpTransitions, rewardsFn, 0.9)
    policy.matrix[isubject,] <- ipolicy.best$policy
  }# end-for isubject
  
  writeMDPTable(policy.matrix,"policies_loocv","all",iPerturbation,dirouteval)
  
  return(policy.matrix)
  
}

# Computes the frequency a subject is following a given policy
subjectFollowedPolicy <- function(phyloObject, subject, policy, perturbation){
  
  print(paste("Testing subject following policy",subject))
  
  mdp.names <- getMDPnames(phyloObject, perturbation)
  mdp.states <- mdp.names[[1]]
  
  testData <- phyloSubset(phyloObject,subjectVar,subject)
  transObserved <- transitionTriplets(testData, perturbation, 
                                      timeStepId = stepVar, subjectId=subjectVar, verbose=FALSE)
  
  transCounts <- transitionCounts(transObserved, mdp.names)
  followCounts <- 0
  
  for(i in 1:length(mdp.states)){
    followCounts <- followCounts + sum(transCounts[i, ,policy[i]])
  }
  
  followed <- followCounts/sum(transCounts)
  
  return(followed)
}


# Compute the percentage of times a subject is in a goal state
subjectInGoal <- function(phyloObject, subject, statePref, perturbation){
  print(paste("Testing subject being in Goal state:",subject))
  
  mdp.names <- getMDPnames(phyloObject, perturbation)
  mdp.states <- mdp.names[[1]]
  
  testData <- phyloSubset(phyloObject,subjectVar,subject)
  transObserved <- transitionTriplets(testData, perturbation, 
                                      timeStepId = stepVar, subjectId=subjectVar, verbose=FALSE)
  
  transCounts <- transitionCounts(transObserved, mdp.names)
  
  inGoalCounts <- 0
  
  for(i in 1:length(mdp.states)){
    inGoalCounts <- inGoalCounts + statePref[i]* sum(transCounts[ ,i,])
  }
  
  freqInGoal <- inGoalCounts/sum(transCounts)
  
  return(freqInGoal)
  
}

# Some vectors with zero observation (pollicy not followed) should return NA value
quietMean <- function(vec){
  if (length(vec) > 0){
    return(mean(vec))
  }else{
    return(NA)
  }
}

subjectDeltaGoal <- function(phyloObject, subject, policy, goalVar, perturbation){
  print(paste("Testing subject for Delta",goalVar, ":",subject))
  mdp.names <- getMDPnames(phyloObject, perturbation)
  mdp.states <- mdp.names[[1]]
  mdp.actions <- mdp.names[[3]]
  
  labelPolicy <- mdp.actions[policy]
  names(labelPolicy) <- mdp.states
  
  testData <- phyloSubset(phyloObject,subjectVar,subject)
  transObserved <- transitionTriplets(testData, perturbation, 
                                      timeStepId = stepVar, subjectId=subjectVar, goalVar=goalVar,
                                      verbose=FALSE, saveResult = FALSE)
  
  utility_follow <- numeric()
  utility_notfollow <- numeric()
  
  for (i in 1:nrow(transObserved)){
    itran <- transObserved[i,]
    if (labelPolicy[as.character(itran$Q.ini)] == itran$Action){
      utility_follow <- c(utility_follow, itran$utility)
    }else{
      utility_notfollow <- c(utility_notfollow, itran$utility)
    }
  }
  
  avg_follow <- quietMean(utility_follow)
  avg_notfollow <- quietMean(utility_notfollow)
  avg_utility <- mean(transObserved$utility)
  
  return(c(avg_utility, avg_follow, avg_notfollow))
  
}

# Return a vector [better, equal, worse] for S->S' transition valuation
transitionProgress <- function(qini, qend, cluster_pref){
  if (qini == qend){
    value <- c(0,1,0)  
  }else if (cluster_pref[qini] > cluster_pref[qend]){
    value <- c(0,0,1)
  }else{
    value <- c(1,0,0)
  }
  return(value)
}


subjectFollowProgress <- function(phyloObject, subject, policy, cluster_pref, perturbation){
  print(paste("Testing subject for following [yes/no] Transitions:",subject))
  mdp.names <- getMDPnames(phyloObject, perturbation)
  mdp.states <- mdp.names[[1]]
  mdp.actions <- mdp.names[[3]]
  
  labelPolicy <- mdp.actions[policy]
  names(labelPolicy) <- mdp.states
  
  testData <- phyloSubset(phyloObject,subjectVar,subject)
  transObserved <- transitionTriplets(testData, perturbation, 
                                      timeStepId = stepVar, subjectId=subjectVar, goalVar=goalVar,
                                      verbose=FALSE, saveResult = FALSE)
  
  follow<- c(0,0,0) 
  not_follow <- c(0,0,0)
  
  for (i in 1:nrow(transObserved)){
    itran <- transObserved[i,]
    if (labelPolicy[as.character(itran$Q.ini)] == itran$Action){
      follow <- follow + transitionProgress(as.character(itran$Q.ini),
                                            as.character(itran$Q.end), cluster_pref)
      #print(itran)
      #print(follow)
    }else{
      not_follow <- not_follow + transitionProgress(as.character(itran$Q.ini),
                                                    as.character(itran$Q.end), cluster_pref)
      #print(itran)
      #print(not_follow)
    }
  }
  
  return(c(follow, not_follow))
}
