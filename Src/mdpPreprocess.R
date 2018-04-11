# Pre-processing steps

MDPpreprocess <- function(phyloObject, processDir=dirdata){
 
  otutablePrenormfile <- paste(processDir,"otuprenorm.tsv",sep="")
  otutablePostnormfile <- paste(processDir,"otunormalized.tsv",sep="")
  
  phyloObject <- prune_taxa(taxa_sums(phyloObject) > 0, phyloObject)

  otuTablePrenorm <- as.matrix(t(otu_table(phyloObject)))
  write.table(otuTablePrenorm,otutablePrenormfile,quote=FALSE,sep='\t',col.names=FALSE,row.names=FALSE)

  # Call to python script to apply the normalization computed by [David,2014]:
  print ("Applying normalization (David et al, 2014)...")
  pycmd <- paste("python normalize_David2014.py -i", otutablePrenormfile,"-o", otutablePostnormfile)
  system(pycmd)

  dfPostnorm <- read.table(otutablePostnormfile,sep='\t')
  # To put the names from otuTablePrenorm
  rownames(dfPostnorm) <- colnames(otu_table(phyloObject))
  colnames(dfPostnorm) <- rownames(otu_table(phyloObject))
  otuTablePostnorm <- t(data.matrix(dfPostnorm[1:(nrow(dfPostnorm)),1:(ncol(dfPostnorm))]))
  OTU <- otu_table(otuTablePostnorm,taxa_are_rows=TRUE)
  data.norm <- phyloseq(OTU, tax_table(phyloObject), sample_data(phyloObject))
 
  return(data.norm) 
}
