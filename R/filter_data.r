# violin <- function(data, y, split_by, quants = NULL){
#   vplot <- ggplot2::ggplot(data, ggplot2::aes_string(y = y, x = split_by)) +
#     ggplot2::geom_violin(ggplot2::aes(fill = "green"), show.legend = F) +
#     ggplot2::geom_jitter(size = 0.5)
#     if(!is.null(quants)){vplot <- vplot + ggplot2::geom_hline(yintercept = quants)}
#   return(vplot)
# }

scatter<- function(data, x, y, color = NULL){
  Splot <- ggplot2::ggplot(data, ggplot2::aes_string(x = x, y = y, color = color)) + ggplot2::geom_point() #attach plot to obj
  return(Splot)
}

Attach_QC <- function(sc_obj, pattern = "default"){
  if(pattern == "default"){
    pattern = c("^MT-", "^RP")
  }

  glist <- list()

  for(p in pattern){ #iterate through patterns and compare grouped transcripts as fractions of cell nCount totals
    fgenes <- grep(dimnames(sc_obj)[[1]], pattern = p)
    fcounts <- Matrix::colSums(sc_obj[fgenes,])
    fpercent <- fcounts/sc_obj[[,"nCount_RNA"]] * 100
    col.name <- paste("percent", gsub(pattern = "[^[:alnum:]]", "", p), sep = "")
    sc_obj <- addMetaData(sc_obj, fpercent, col.name = col.name)

    fplot <- scatter(data = sc_obj[[]], y = col.name, x = "nCount_RNA")
    glist <- c(glist, list(fplot))
  }

  nCount_qc <- scatter(data = sc_obj[[]], x = "nCount_RNA", y = "nFeatures_RNA" )

  glist<- c(glist, list(nCount_qc))

  complete_qc <- gridExtra::grid.arrange(grobs = glist)

  sc_obj <- addGraph(x = sc_obj, value = complete_qc)

  return(sc_obj)
}

detect_anomalies <- function(sc_obj, method = "isolation"){
  sc_obj <- switch(method,

                   "isolation" = {isolate(sc_obj)},

                   "doublet_scoring"= {anomalies <- est_doublets(sc_obj)
                                       cat("Predicted doublet rate: ", anomalies[[1]], "\n")
                                       addMetaData(sc_obj, anomalies[[2]], col.name = "Anomaly")
                                      }
                   )

  glist <- list()
  QC_cols <- colnames(sc_obj[[,-1]])
  QC_cols <- QC_cols[-length(QC_cols)]
  for(value in QC_cols){
    fplot <- scatter(sc_obj[[]], x = "nCount_RNA", y = value, color = "Anomaly")
    glist <- c(glist, list(fplot))
  }
  complete_qc <- gridExtra::grid.arrange(grobs = glist)
  sc_obj <- addGraph(x = sc_obj, value = complete_qc)

  return(sc_obj)
}

rm_anomalies <- function(sc_obj){
  if(!"Anomaly" %in% colnames(sc_obj[[]])){
    stop("Object 'Anomaly' not found in metadata, have you tried running detect_anomalies()?")
  }
  data.rm_anom <- sc_obj[][,which(!sc_obj[[,"Anomaly"]])] #TODO check access method
  sc_obj@assay <- list(anomalies_removed = data.rm_anom) #TODO change access method
  cat(ncol(sc_obj@raw.data)-ncol(sc_obj@assay$anomalies_removed), "anomalies not carried over")
  return(sc_obj)
}
