#' scatter
#'
#' Covenience function for generating scatter plots with ggplot2
#'
#' @param data
#' @param x
#' @param y
#' @param color
#'
#' @return ggplot
#' @importFrom ggplot2 ggplot aes_string geom_point
#'
scatter<- function(data, x, y, color = NULL){
  Splot <- ggplot2::ggplot(data, ggplot2::aes_string(x = x, y = y, color = color)) + ggplot2::geom_point(size = 0.25) #attach plot to obj
  return(Splot)
}

#' Attach_QC
#'
#' Automates quality control (QC) steps and attaches QC metrics to sc_obj
#'
#' @param sc_obj an sc_obj that containing a single cell RNA sequencing count matrix
#' @param pattern a regular expression or vector of regular expressions to select
#' the features for which outliers should be removed. The default pattern is
#'   pattern = ("^MT-", "^RP").
#'
#' @importFrom Matrix colSums
#' @importFrom gridExtra arrangeGrob
#'
#' @return sc_obj
#' @export
#'
#' @examples
#' data.path <- "path/to/file/"
#' counts <- import10x(data.path)
#' sce <- newSC_obj(counts)
#' sce <- Attach_QC(sce, pattern = "^MT")
#'
#'
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

  cat("generating raw data QC... \n")
  complete_qc <- gridExtra::arrangeGrob(grobs = glist)

  sc_obj <- addGraph(x = sc_obj, value = complete_qc)

  return(sc_obj)
}

#' detect_anomalies
#'
#' Uses isolation forest to perform multivariate anomaly detection across all columns
#' present in sc_obj metadata
#'
#' @param sc_obj an sc_obj with attached QC
#' @param method By default method is set to "isolation" for anomaly detection. "doublet_scoring" may be used
#' to perform troubleshooting tasks but will not detect anomalies.
#' @param verbose
#'
#' @importFrom gridExtra arrangeGrob
#'
#' @return an sc_obj with identities of outlier events included in the metadata.
#' Visualization for anomaly detection is included in the sc_obj graphs.
#' @export
#'
#' @examples
#'
#'   sce <- detect_anomalies(sce)
#'
detect_anomalies <- function(sc_obj, method = "isolation", verbose = T){
  cat("simulating doublets... \n")
  sc_obj <- switch(method,

                   "isolation" = {
                      doublets = est_doublets(sc_obj, method = 'cxds')
                      sc_obj <- set_doublet_rate(sc_obj, doublets[[1]])
                      cat("detecting anomalies... \n")
                      isolate(sc_obj)
                     },

                   "doublet_scoring"= {
                      doublets <- est_doublets(sc_obj)
                      cat("Predicted doublet rate: ", doublets[[1]], "\n")
                      sc_obj <- addMetaData(sc_obj, doublets[[2]], col.name = "Anomaly")
                      set_doublet_rate(sc_obj, doublets[[1]])
                    }
                   )

  glist <- list()
  QC_cols <- colnames(sc_obj[[,-1]])
  QC_cols <- QC_cols[-length(QC_cols)]
  for(value in QC_cols){
    fplot <- scatter(sc_obj[[]], x = "nCount_RNA", y = value, color = "Anomaly")
    glist <- c(glist, list(fplot))
  }
  cat("generating anomaly detection QC...  \n")
  complete_qc <- gridExtra::arrangeGrob(grobs = glist)
  sc_obj <- addGraph(x = sc_obj, value = complete_qc)

  return(sc_obj)
}

#' rm_anomalies
#'
#' Convenience function which generate a new sc_obj assay containing the dataset
#' after anomalies have been removed. Raw data remains intact.
#'
#' @param sc_obj
#'
#' @return sc_obj with anomalies removed from the active dataset
#' @export
#'
#' @examples
#'   sce <- rm_anomalies(sce)
rm_anomalies <- function(sc_obj){
  if(!"Anomaly" %in% colnames(sc_obj[[]])){
    stop("Object 'Anomaly' not found in metadata, have you tried running detect_anomalies()?")
  }
  data.rm_anom <- sc_obj[][,which(!sc_obj[[,"Anomaly"]])] #TODO check access method
  sc_obj@assay <- list(anomalies_removed = data.rm_anom) #TODO change access method
  cat("[*_*] done. \n")
  cat(ncol(sc_obj@raw.data)-ncol(sc_obj@assay$anomalies_removed), "anomalies not carried over \n")
  return(sc_obj)
}
