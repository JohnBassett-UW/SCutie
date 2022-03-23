data.path <- "C:/Users/jbassett/Desktop/github/Data/UW36/filtered_feature_bc_matrix.h5"
pipeline <- function(data.path){
  ptm = proc.time()
  counts <- import10x(data.path)
  sce <- newSC_obj(counts)
  sce <- Attach_QC(sce)
  sce <- detect_anomalies(sce)
  sce <- rm_anomalies(sce)
  
  cat("Total \n time ")
  print((proc.time() - ptm)[3]) #print time elapsed
  return(sce)
}

sce <- pipeline(data.path)
