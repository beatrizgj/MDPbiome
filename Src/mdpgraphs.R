# MDPbIOME
# Function for generating output graphs
library(ggplot2)
library(circlize)

theme_update(plot.title = element_text(hjust = 0.5)) # center title in ggplot

# Writes the graph file for the MDP and creates a PDF it calling dot
#   transitionFreqs: array [state, state, action]-> Prob, with MDP transition table
#   actionVar      : name of the action group
#   fileName       : (optional) for changing default path
mdpGraph <- function(transitionFreqs, actionVar, diroutput="", fileName=NULL){
  
  actiondir <- paste(diroutput,actionVar,"/",sep="")
  setdirectory(actiondir,paste(actionVar,"output"))
  transition_file <- paste(actiondir,'transitions_',tolower(actionVar),'.txt',sep='') 
  
  if(missing(fileName)){
    graphName = paste("mdp_diagram_",tolower(actionVar), sep = "")
    fileName <- paste(actiondir, graphName,".gv",sep='')
    fileNamePdf <- paste(actiondir,graphName,".pdf",sep='')
  }
  
  states <- dimnames(transitionFreqs)[[1]]
  actions <- dimnames(transitionFreqs)[[3]]
  
  unlink(fileName)
  fout<-file(fileName,open="a")
  write("digraph{", fout,append=TRUE)
  
  for (qi in states){  
    for (qf in states){ 
      for (act in actions){
        freq <- round(transitionFreqs[qi,qf,act], 4)
        if(freq > 0){			
          graphinfo <- paste('\t',gsub("-","_",qi),' -> ',gsub("-","_",qf),' [label="',
                             act,':',freq,';"]',sep='')          
          write(graphinfo,fout,append=TRUE)
        } 
      } 
    } 
  } 
  write("// title",fout,append=TRUE)
  write("labelloc=\"t\";",fout,append=TRUE)
  write(paste("label=\"",actionVar,"\";",sep=""),fout,append=TRUE)        
  write("}",fout,append=TRUE)
  close(fout)
  system(paste("dot -Tpdf '",fileName,"' -o '",fileNamePdf,"'",sep=''))	

}


#### Plot transition diagram with CHROD:

mdpChrodGraph <- function(mdpTrans, actionVar, diroutput="", online=FALSE){
  # From Brewer. We set 10 colors and then select depending on the number of states
  colorPalette <- c('#1f78b4','#33a02c','#e31a1c','#ff7f00','#6a3d9a',
                    '#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#cab2d6')
  
  states <- dimnames(mdpTrans)[[1]]
  actions <- dimnames(mdpTrans)[[3]]

  # https://cran.r-project.org/web/packages/circlize/vignettes/visualize_relations_by_chord_diagram.pdf
  
  actiondir <- paste(diroutput,actionVar,"/",sep="")

  for (iact in actions){
    mat=matrix(mdpTrans[,,iact], length(states))
    rownames(mat)=rownames(mdpTrans)
    colnames(mat)=colnames(mdpTrans)
    
    if (!online){
      pdf(paste(actiondir, 'chrod_',tolower(actionVar),'_',iact,'.pdf',sep=''))
    }
    
    par(cex=2, mar = c(0.5, 1, 1, 1))
    circos.par(start.degree=135)
    # Colors: http://i0.wp.com/bc.bojanorama.pl/wp-content/uploads/2013/04/rcolorsheet-0.png
    #grid.col = c("darkblue", "darkgreen", "orangered", "greenyellow", "red")
    #grid.col = c("darkblue", "red")
    #grid.col = c("#1b9e77","#d95f02","#7570b3","#e7298a")
    grid.col <- colorPalette[1:length(states)]
    
    chordDiagram(mat, grid.col=grid.col, directional = TRUE, direction.type=c("diffHeight+arrows"), 
                 link.arr.type="big.arrow", annotationTrack=c('grid','name'), 
                 order=sort(states,decreasing=TRUE),
                 transparency=0.2)
    title(main=paste(actionVar,": ",iact,sep=''),cex.main=1)
    circos.clear()
    if(!online){
      dev.off()
    }
  } 
}

combinedTransitionGraph <- function(actionVar, outputdir){
  load(paste(diroutput, actionVar, "/table_transitions-all.RData", sep=""))
  png(paste(outputdir,"combinedTransitions_",actionVar,".png",sep=""))
  par(mfrow=c(2,2))
  mdpChrodGraph(tableFn, actionVar, online = TRUE)
  dev.off()
}



# Creates barplots for showing the percentage of time an action belongs 
# to a policy in the evaluation
#   freqs: State x Action matrix containing the frequency each action belongs to a optimal policy
#   policy    : Optimal policy to mark the "nominal" action 
#   actionVar : Perturbation variables
actionEvalGraph <- function(freqs, policy, reward, actionVar, diroutput, graphname=""){
  actiondir <- paste(diroutput,actionVar,"/",sep="")

  states <- dimnames(freqs)[[1]]
  actions <- dimnames(freqs)[[2]]
  npoints <- length(states) * length(actions)

  actionFreq <- as.data.frame(freqs)
  
  actionFreq$state <- rownames(actionFreq)
  dataFreq <- melt(actionFreq, id.vars = "state")
  colnames(dataFreq) <- c("state","action","stability")
  dataFreq$mark <- rep(FALSE,npoints)
  for (i in 1:length(states)){
    dataFreq$mark[dataFreq$state==states[i] & dataFreq$action==actions[policy[i]]] <- TRUE
  }
  
  pdf(paste(actiondir, graphname, tolower(actionVar),'_', reward, '.pdf',sep=''))
  gbar <- ggplot(dataFreq, aes(x=state,y=stability, fill=action, alpha=as.factor(mark))) + 
    geom_bar(stat="identity") +
    scale_alpha_manual(values=c("FALSE"=0.4, "TRUE"=1.0), guide='none') + 
    theme(panel.background = element_rect(fill = 'white'),
          text = element_text(size=24)) + 
    ggtitle(actionVar) 

  print(gbar)
  dev.off()
}

# This shows the percentage 
followPolicyGraph <- function(followdf,actionVar, actiondir){
  pdf(paste(actiondir, 'corrfollowedpolicy_graph_',tolower(actionVar),'.pdf',sep=''))
  #colnames(followdf) <- c("% Followed Policy","% Leading to GoalState")
  gcor <- ggplot(followdf, aes(x=policy,y=Avg.Utility)) + geom_point(size=2)   
  print(gcor)
  dev.off()
  
  # Double points for following or not following policies.
  double.df <- followdf[,c(1,3,4)]
  graphdf <- melt(double.df, id.vars = "policy")
  names(graphdf) <- c("Freq.Policy","Criteria","Utility")
  
  pdf(paste(actiondir, 'corrsplit_followedpolicy_graph_',tolower(actionVar),'.pdf',sep=''))
  gdouble <- ggplot(graphdf, aes(x=Freq.Policy,y=Utility, colour=Criteria)) + geom_point(size=2)
  print(gdouble)
  dev.off()
  
}

followBarsGraph <- function(followResults, dirout, fname=NULL, glegend=TRUE){
  gdata <- melt(followResults, id.vars = c("Perturbation","Follow"))
  colnames(gdata) <- c("Perturbation","Follow","Transition","Frequency")
  if (missing(fname)){
    fname <- 'followBars_graph.pdf'
  }
  pdf(paste(dirout, fname,sep=''))
  gobj <- ggplot(gdata, aes(x=Follow, y=Frequency, fill=Transition)) + 
          #geom_bar(stat = 'identity', position = 'stack') + 
          geom_bar(stat = 'identity') +
          theme(text = element_text(size=20)) +
          scale_fill_manual(values=c('#4daf4a','#377eb8','#e41a1c')) +
          facet_grid(~ Perturbation) + 
          ggtitle(titledata)
  
  if(!glegend){
    gobj + guides(fill=FALSE)
  }
  #panel.background = element_rect(fill = 'white'),
  print(gobj)
  dev.off()
}


summarizeFollow <- function(followdf, eqType=0){
  if (eqType==1){
    resdf <- data.frame(Perturbation = followdf$Perturbation,
                Follow = followdf$Follow,
                Better = followdf$Better / (followdf$Better + followdf$Equal + followdf$Worse),
                EqWorse = (followdf$Equal + followdf$Worse) / (followdf$Better + followdf$Equal + followdf$Worse))
                        
    
  }else if (eqType == 2){
    resdf <- data.frame(Perturbation = followdf$Perturbation,
                        Follow = followdf$Follow,
                        EqBetter = (followdf$Better + followdf$Equal) / (followdf$Better + followdf$Equal + followdf$Worse),
                        Worse = followdf$Worse / (followdf$Better + followdf$Equal + followdf$Worse))
  }else{
    resdf <- data.frame(Perturbation = followdf$Perturbation,
                        Follow = followdf$Follow,
                        Better = followdf$Better / (followdf$Better + followdf$Equal + followdf$Worse),
                        Equal = followdf$Equal / (followdf$Better + followdf$Equal + followdf$Worse),
                        Worse = followdf$Worse / (followdf$Better + followdf$Equal + followdf$Worse))
    
  }
  return(resdf)
}

