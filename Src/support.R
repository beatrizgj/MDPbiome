
diranalysis <- ""
diroutput <- ""
diroutput_loocv <- ""

createTreeDir <- function(basedir, Perturbation){
  diranalysis <<- paste(basedir,"output_analysis/",sep = "")
  diroutput <<-paste(basedir,"output/",sep = "")
  diroutput_loocv <<- paste(basedir,"output_loocv/",sep = "")
  setdirectory(diranalysis, "analysis graphs", quiet = FALSE)
  setdirectory(diroutput, "output data", quiet = FALSE)
  setdirectory(diroutput_loocv, "output loocv data", quiet = FALSE)
  
  for (idir in c(diroutput, diroutput_loocv)){
    for (jsubdir in Perturbation){
      jpath <- paste(idir,jsubdir, sep="")
      setdirectory(jpath,"",quiet = TRUE)
    }
  }
}
    
setdirectory <- function(newpath, tagmsg, quiet =TRUE){
  if (dir.exists(newpath)){
    if (!quiet){
      print(paste("Local dir for",tagmsg,"already exists"))
    }
  }else{
    dir.create(newpath)
    print(paste("New local dir for",tagmsg))
  }
}

datapath <- function(filename, raw=FALSE){
  if (raw){
    dir <- paste(dirdata,filename, sep = "")
  }else{
    dir <- paste(dircleandata,filename, sep = "")
  }
  return (dir)
}


