violin <- function(data, y, split_by, quants = NULL){
  vplot <- ggplot2::ggplot(data, ggplot2::aes_string(y = y, x = split_by)) +
    ggplot2::geom_violin(ggplot2::aes(fill = "green"), show.legend = F) +
    ggplot2::geom_jitter(size = 0.5)
    if(!is.null(quants)){vplot <- vplot + ggplot2::geom_hline(yintercept = quants)}
  return(vplot)
}

scatter<- function(data, x, y){
  Splot <- ggplot2::ggplot(data, ggplot2::aes_string(x = x, y = y)) + ggplot2::geom_point() #attach plot to obj
  return(Splot)
}

Perform_QC <- function(sc_obj, pattern = "default", dbl_params = NULL){
  if(pattern == "default"){
    pattern = c("^MT-", "^RP")
  }
  if(!is.null(dbl_params)){
    sc_obj <- addMetaData(sc_obj, dbl_params[[2]], col.name = "sim_doublets")
  }

  glist <- list()

  for(p in pattern){ #iterate through patterns and compare grouped transcripts as fractions of cell nCount totals
    fgenes <- grep(dimnames(sc_obj)[[1]], pattern = p)
    fcounts <- Matrix::colSums(sc_obj[fgenes,])
    fpercent <- fcounts/sc_obj[["nCount_RNA"]] * 100
    col.name <- paste("percent", gsub(pattern = "[^[:alnum:]]", "", p), sep = "")
    sc_obj <- addMetaData(sc_obj, fpercent, col.name = col.name)

    fplot <- scatter(data = sc_obj[[]], y = col.name, x = "nCount_RNA")
    glist <- c(glist, list(fplot))
  }

  nCount_qc <- violin(data = sc_obj[[]], y = "nCount_RNA", split_by = "sim_doublets")
  nfeature_qc <- violin(data = sc_obj[[]], y = "nFeatures_RNA", split_by = "sim_doublets")

  glist<- c(glist, list(nCount_qc, nfeature_qc))

  complete_qc <- gridExtra::grid.arrange(grobs = glist)

  sc_obj <- addGraph(x = sc_obj, value = complete_qc)

  return(sc_obj)
} #TODO remove plots from prefom qc
#Add new function for plotting graphs
#use grid arrange to compile all initial qc plots in a single instance
#Add all initial graphs to graphs category in one instance


QC_filter <- function(sc_obj){ #TODO filter data based on doublet estimation
  est_doublets <- sc_obj[[sc_obj[[,"sim_doublets"]] == T, ]]
  qfeatures <- quantile(est_doublets[,"nFeatures_RNA"])
  qcounts <- quantile(est_doublets[,"nCount_RNA"])

  #mitochondrial outliers
  qMT <- quantile(sc_obj[[,"percentMT"]], probs = seq(0, 1, 0.01))
  outsMT <- min(qMT[qMT >= sd(qMT)])

  #ribosomal outliers
  qRP <- quantile(sc_obj[[,"percentRP"]], probs = seq(0.5, 1, 0.01))
  IQRrp <- IQR(sc_obj[[,"percentRP"]])
  outsRP <- qRP[[4]]+(1.7* IQRrp)

  violin(sc_obj[[]], y = "nFeatures_RNA", split_by = "sim_doublets", quants = qfeatures)
  violin(sc_obj[[]], y = "nCount_RNA", split_by = "sim_doublets", quants = qfeatures)
  violin(sc_obj[[]], y = "percentMT", split_by = "sim_doublets", quants = outsMT)
  violin(sc_obj[[]], y = "percentRP", split_by = "sim_doublets", quants = outsRP)
  hist(sc_obj[[,"percentMT"]], breaks = 100)
  hist(sc_obj[[,"percentRP"]], breaks = 500)


  #set filter thresholds
  feature_thresh  = qfeatures[[4]]
  count_thresh  = qcounts[[4]]
  MT_thresh = 0.20
  RP_thresh = 0.20


}
