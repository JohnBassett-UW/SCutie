#add credit to solitude
isolate <- function(sc_obj, max_depth = NULL, verbose = F){
  if(is.null(max_depth)){
    max_depth = (ceiling(log2(nrow(sc_obj[[]]))))
  }
  anomaly_rate = 1
  while(anomaly_rate >= doublet_rate(sc_obj) | max_depth > 500){
    max_depth = 2*max_depth
    x = sc_obj[[]]
    iforest <- solitude::isolationForest$new(sample_size = ceiling(nrow(x)/10), 
                                             seed = 101,
                                             num_trees = 100,
                                             max_depth = max_depth
    )
    iforest$fit(x)
    scores = iforest$predict(x)
    Amin <- min(scores[,3])
    Amax <- max(scores[,3])
    if(verbose){cat("Anomaly Score Range: ", Amin, "to" , Amax, "\n")}
    if(Amin > 0.5){
      warning("Min anomaly score > 0.5")
      stop("Outlier Detection failed")
    }
    index <- (which(scores[,3]>0.5))
    anomalies <- rep(F, nrow(sc_obj[[]]) )
    anomalies[index] <- T
    anomaly_rate = length(which(anomalies > 0))/nrow(sc_obj[[]])
  }
  sc_obj <- addMetaData(sc_obj, anomalies, col.name = "Anomaly")
  if(verbose){cat("Predicted anomaly rate: ", anomaly_rate, "\n")}
  return(sc_obj)
}