

# list: vector of integers of time slots
# curr: current time slot
next.time.slot <- function(list, curr){
  i=1
  while(curr != list[i]){
    i=i+1
  } 
  return(list[i+1])
} 

# Simple normalization for a vector
normalize <- function(vec){
  normalvec <- (vec - min(vec))/(max(vec) - min(vec))
  return(normalvec)
}


phyloSubset <- function(phyloOb, Var, Value){
  remove_idx <- as.character(get_variable(phyloOb, Var)) == as.character(Value)
  filteredOb <- prune_samples(remove_idx, phyloOb)
  return(filteredOb)
}


phyloFilterOut <- function(phyloOb, Var, Value){
  remove_idx <- as.character(get_variable(phyloOb, Var)) != as.character(Value)
  filteredOb <- prune_samples(remove_idx, phyloOb)
  return(filteredOb)
}


# create the [state,state,action] dimension names
getMDPnames <- function(phyloObject, perturbation){
  mdpstates <- levels(sample_data(phyloObject)$cluster)
  
  actions <- as.factor(get_variable(phyloObject,perturbation))
  mdpactions <- levels(actions)
  
  mdpdim <- list(mdpstates, mdpstates, mdpactions)
  return(mdpdim)
}


# Write transition triples from phyloseq object with cluster annotation
#   actionVar: name of the actions.
#   timeStepId: name of the variable to use as the timestep identifier
transitionTriplets <- function(phyloData, actionVar, timeStepId, subjectId="Subject", 
                               goalVar=NULL, verbose=TRUE, saveResult=TRUE, dirout=diroutput){

  sample_data(phyloData)$timeStep <- as.factor(get_variable(phyloData,timeStepId))
  sample_data(phyloData)$action <- as.factor(get_variable(phyloData,actionVar))
  
  if (!missing(goalVar)){
    sample_data(phyloData)$goalvalue <- get_variable(phyloData,goalVar)
  }
  
  transitions <- data.frame(Q.ini=integer(), Action=character(), Q.end=integer(), stringsAsFactors=FALSE)
  utility <- numeric()
  
  list.subjects <- unique(sort(get_variable(phyloData, subjectId)))
  
  i <- 0
  for (isubject in list.subjects){
    if(i %% 10 == 0 && verbose){
      print(paste(' Parsing examples from subject:',isubject, "... please wait"))
    }
    i <- i + 1
    
    phyloData.1subject <- NULL
    
    phyloData.1subject <- phyloSubset(phyloData,subjectId,isubject)
    
    timeSlots <- sort(get_variable(phyloData.1subject,timeStepId))
    max.time <- timeSlots[length(timeSlots)] 
    for (day in timeSlots){
      newrow <- NULL
      sample.ini <- sample_data(phyloSubset(phyloData.1subject,"timeStep",day))
    
      # We assume that perturbation is related to the observed sample, as opposite of 
      # traditional MDP where the action is stored as it is being applied in the current state
        
      if(day != max.time){ # There are possible new transitions
        next.day <- next.time.slot(timeSlots,day)
        next.sample.subset <- phyloSubset(phyloData.1subject,"timeStep",next.day)
        next.sample <- sample_data(next.sample.subset)
        
        # NA comes as a valid value for factors, but this means the action meta-data has not been registered
        if (next.sample$action[1] != 'NA'){
          newrow <- data.frame(Q.ini=sample.ini$cluster, Action=next.sample$action, Q.end=next.sample$cluster)
          #print(paste(sample.ini$cluster,next.sample$action,next.sample$cluster))
          transitions <- rbind(transitions, newrow)
          if (!missing(goalVar)){
            deltaGoalValue <- next.sample$goalvalue - sample.ini$goalvalue
            utility <- c(utility, deltaGoalValue)
          }
        } #end-if NA check
      } #end-if new transition	
    } #end-for days
  } #end-for subjects
  
  if(!missing(goalVar)){
    transitions$utility <- utility
  }
  
  if(saveResult){
    actiondir <- paste(dirout,actionVar,"/",sep="")
    setdirectory(actiondir,paste(actionVar,"output"))
    transition_file <- paste(actiondir,'transitions_',tolower(actionVar),'.txt',sep='') 
    write.table(transitions,transition_file,row.names=FALSE,quote=FALSE)
  }
  
  return(transitions)
}


# Counts the number of occurrence from the observed transitions.
# For some computation we need the absolute value rather than the freq.
transitionCounts <- function(transitions, dimnames){
  
  num.states <- length(dimnames[[1]])
  num.actions <- length(dimnames[[3]])
  counts <- array(data = 0, dim = c(num.states,num.states,num.actions), dimnames = dimnames)
  
  # Counts the number of transitions: Fill in the transition probability table with the absolute number of transitions
  # between each pair of states, per action.
  for (i in 1:nrow(transitions)) {
    # as.character is important, to take the name and not the index of the value!!!!!!
    q.ini <- as.character(transitions[i, "Q.ini"]) 
    q.end <- as.character(transitions[i, "Q.end"])
    action <- as.character(transitions[i, "Action"])
    counts[q.ini,q.end,action] <- counts[q.ini,q.end,action] + 1
  }
  return(counts)
}

# Compute the probability from the observed frequency of each transition
#    transtitions: data.frame with triples (Q.ini, Action, Q.end)
# dimnames = list(sort(levels(sample_data(data.norm)$cluster)),
#                 sort(levels(sample_data(data.norm)$cluster)),
#                  sort(levels(transitions$Action))) 

transitionFrequencies <- function(counts, dimnames, correction=""){

  num.states <- length(dimnames[[1]])
  num.actions <- length(dimnames[[3]])
  
  # Correction for avoiding 0-probabilities of unobserved events
  if(correction=="Laplace"){
    adj.counts <- counts + 1
  }else{
    adj.counts <- counts
  }
  
  # Frequencies: Dividing each row in each table by the total number of transitions per row, to get the frecuencies per state and action.
  P <- array(data = 0, dim = c(num.states,num.states,num.actions), dimnames = dimnames)
  
  for (action in 1:num.actions){
    P[,,action] <- t(apply(adj.counts[,,action], 1, function(x){x/sum(x)}))
  } 
  P[is.na(P)] <- 0
 
  return(P) 
}



setEmptyTransitions <- function(transTable, countsTable){
  states <- dimnames(transTable)[[1]]
  actions <- dimnames(transTable)[[3]]
  
  for (istate in states){
    for (jaction in actions){
      if(sum(transTable[istate,,jaction])==0){
        print(paste("Updating Void Transitions:",istate,jaction))
        for(kstate in states){
          # The probability of P(S'|S) which is independent of the applied action
          probNextState <- sum(countsTable[istate,kstate,])/sum(countsTable[istate,,])
          transTable[istate,kstate,jaction]  <- probNextState
        }
        print(transTable[istate,,jaction])
        #P[q.ini,,action] <- (1/num.states)
      } 
    } 
  } 
  return(transTable)
}
  


computeDiversity <- function(prenormData, phyloData){

  sample_data(prenormData)$cluster <- sample_data(phyloData)$cluster   
  # Compute alpha diversity by cluster
  alphaDiversity <- NULL
  #sink('output_alphaDiversityValues.txt')
  print('Alpha diversity per cluster:')
  for (i in levels(sample_data(phyloData)$cluster)){
    #print(estimate_richness(prenormData,split=TRUE,measure='Shannon')) # All samples individually
    data.clus.i <- phyloSubset(prenormData,"cluster",i)
    
    idiversity <- mean(estimate_richness(data.clus.i,split=TRUE,measure='Shannon')$Shannon)
    alphaDiversity <- c(alphaDiversity, idiversity)
    
    print(paste('cluster',i,': ',idiversity))
  } 
  #sink()
  names(alphaDiversity) <- levels(sample_data(phyloData)$cluster)
  return(alphaDiversity)
}

# Creates the reward function for the MDP
#   states : Mdp states labels
#   actopms: MDP action labels
#   clusterDiversity: alpha diversity pre-computed for each cluster
#   method: rule for setting reward function as follow:
#           "avoidBad" to escape the worst state. Reward 1 to the rest
#           "preferGood" to go to the best state. Reward 1 only to the most diverse cluster
#           "proportinal" to set reward proportionally to the diversity
rewardTable <- function(states, actions, clusterDiversity, method = "preferGood"){
  if (!is.element(method, c("proportional","avoidBad", "preferGood"))){
    stop(paste("Reward function method not recognized:",method))
  }
  
  num.states <- length(states)
  num.actions <- length(actions)
  bestDiversity <- max(clusterDiversity)
  worstDiversity <- min(clusterDiversity)
  
  rewards <- array(data = 0, dim = c(num.states, num.states, num.actions), 
                   dimnames = list(states, states, actions))
  
  if (method =="proportional"){
    normDiversity <- normalize(clusterDiversity)
  
    for (i in 1:length(states)){
      rewards[,states[i],] <- normDiversity[i]
    }
  }
  else{
    
    if (method=="avoidBad"){
      goalStates <- states[clusterDiversity > worstDiversity]
    }else if(method=="preferGood"){
      goalStates <- states[clusterDiversity == bestDiversity]    
    }
    
    for (qf in goalStates){
      rewards[,qf,] <- 1
    }
  }
  
  return(rewards)
}


# It Computes the probability of reaching and end state regardless the evidence
# this will be used as the based concentration parameters for Dirichlet sampling
endStatePrior <- function(countsTable){
  
  mdpstates <- dimnames(countsTable)[[1]]
  #negligile prop above zero
  
  freq <- rep(0.0001, length(mdpstates))
  names(freq) <- mdpstates
  
  for (istate in mdpstates){
    estimatedFreq <- sum(countsTable[,istate,])/sum(countsTable)  
    if(estimatedFreq > 0){
      freq[istate] <- estimatedFreq
    }
  }
  return(freq)
}



# Computes the time serie of cluster for a subject
stateSubjectSerie <- function(phyloData, subject,  timeStepId, subjectId="Subject"){
  sample_data(phyloData)$timeStep <- as.factor(get_variable(phyloData,timeStepId))
  
  subject.data <- phyloSubset(phyloData,subjectId,subject)
  timeSlotsIdx <- order(get_variable(subject.data,timeStepId))
  
  clusSerie <- sample_data(subject.data)[,c(timeStepId,"cluster")]
  sortedSerie <- clusSerie[timeSlotsIdx,]
  print(sortedSerie)
  return(as.character(sortedSerie$cluster))
}


# Computes the table Time X Subjects representing the state serie progression of subjects

stateSerieTable <- function(phyloData, timeStepId, subjectId="Subject", dirout=diroutput){
  
  sample_data(phyloData)$timeStep <- as.factor(get_variable(phyloData,timeStepId))
  
  #transitions <- data.frame(Q.ini=integer(), Action=character(), Q.end=integer(), stringsAsFactors=FALSE)
  #utility <- numeric()
  
  list.subjects <- unique(sort(get_variable(phyloData, subjectId)))
  list.series <- lapply(list.subjects, stateSubjectSerie, phyloData = phyloData, 
                        timeStepId = stepVar, subjectId= subjectVar)
  maxlen <- max(sapply(list.series,length))
  tableserie <- sapply(list.series, function(x){c(x)[1:maxlen]})
  colnames(tableserie) <- list.subjects
  
  return(tableserie)
}


