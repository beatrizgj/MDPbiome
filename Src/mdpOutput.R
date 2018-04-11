# MDPBiome 
# Functions for Output and Evaluate MDPs

# Writes basic MDP solutions to output files 
writeMDPsolution <- function(mdp.sol, actionVar, states, actions, rewardType, diroutput = ""){
  actiondir <- paste(diroutput,actionVar,"/",sep="")
  
  policyFileName <- paste(actiondir, "policy-",rewardType,".txt", sep="")
  
  otext <- paste("Action:",actionVar)
  otext <- c(otext,paste("Rewards:", rewardType))
  otext <- c(otext,"-----------------------------")
  otext <- c(otext,"Optimal Policy (State : Action)")
  for(i in 1:length(states)){
    bestActIdx <- mdp.sol$policy[i]
    otext <- c(otext, paste(states[i], actions[bestActIdx], sep=" : "))
  }
  
  # Print to console up to this point
  for(istr in otext){
    print(istr)
  }
  
  otext <- c(otext,"-----------------------------")
  otext <- c(otext,"Value Function")
  for(i in 1:length(states)){
    otext <- c(otext, paste(states[i], mdp.sol$V[i], sep=" : "))
  }
  
  ostream <- file(policyFileName, open = "w")
  writeLines(otext, ostream)
  close(ostream)
}

# Function to write any kind of MDP information in the form [state, state, action]
#  tag:  string to differentiate tables. eg. "rewards" or "counts"    
writeMDPTable <- function(tableFn, tableTag, rewardType, actionVar, diroutput, saveObject=FALSE){
  actiondir <- paste(diroutput,actionVar,"/",sep="")
  setdirectory(actiondir,paste(actionVar,"output"))
  
  fileName <- paste(actiondir, tableTag,"-",rewardType,".txt", sep="")
  sink(fileName)
  print(tableFn)
  sink()
  
  if(saveObject){
    fileObj <- paste(actiondir, tableTag,"-",rewardType,".RData", sep="")
    save(tableFn, file = fileObj)
  }
}


writeStateTable <- function(tableVal, states, tableTag, diroutput){
  otext <- paste("Info", tableTag)
  
  for(i in 1:length(states)){
    otext <- c(otext, paste(states[i],tableVal[i],sep= " : "))
  }
  
  fileName <- paste(diroutput, tableTag,".txt", sep="")
  ostream <- file(fileName, open = "w")
  writeLines(otext, ostream)
  close(ostream)
}

timeSerieHeatMap <- function(tableSerie, dirouput, fname){
  
  if(missing(fname)){
    fname <- "StateSequence_heatmap.pdf"
  }
  
  pdf(paste(dirouput,fname, sep=""))
  
  graph.data <- melt(tableSerie)
  colnames(graph.data) <- c("Time","Subject","State")
  graph.data$State <- as.factor(graph.data$State)
  
  heatmg <- ggplot(graph.data, aes(x=Time, y=Subject, fill=State)) + geom_tile()
  print(heatmg)
  
  dev.off()
}

