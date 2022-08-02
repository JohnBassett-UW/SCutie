#' scatter
#'
#' Covenience function for generating scatter plots with ggplot2
#'
#' @param data dataset to generate plot from
#' @param x data to map to x axis
#' @param y data to map to y axis
#' @param color color
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
    metaData(sc_obj, col.name) <-  fpercent

    fplot <- scatter(data = sc_obj[[]], y = col.name, x = "nCount_RNA")
    glist <- c(glist, list(fplot))
  }

  nCount_qc <- scatter(data = sc_obj[[]], x = "nCount_RNA", y = "nFeatures_RNA" )

  glist<- c(glist, list(nCount_qc))

  cat("generating raw data QC... \n")
  complete_qc <- list(gridExtra::arrangeGrob(grobs = glist))
  names(complete_qc) <- "Quality Check"
  graphs(sc_obj) <- complete_qc

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
#' @param verbose By default is set to True
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
detect_anomalies <- function(sc_obj, method = "none", verbose = T){
  if(!method %in% c("hybrid", "bcds", "cxds")){
    cat("setting defualt method \n")
    method = "cxds"
    }

  cat("simulating doublets using method ", method, "... \n", sep = "")
  doublets = est_doublets(sc_obj, method = method)
  est_doublet_rate(sc_obj) <- doublets[[1]]
  cat("detecting anomalies... \n")
  sc_obj <- isolate(sc_obj)

  cat("generating anomaly detection QC...  \n")
  glist <- list()
  QC_cols <- colnames(sc_obj[[,2:4]])
  for(value in QC_cols){
    fplot <- scatter(sc_obj[[]], x = "nCount_RNA", y = value, color = "Anomaly")
    glist <- c(glist, list(fplot))
  }
  complete_qc <- list(gridExtra::arrangeGrob(grobs = glist))
  names(complete_qc) <- c("Anomalies")
  graphs(sc_obj) <- complete_qc

  cat("Attaching doublet estimates...  \n")
  metaData(sc_obj, col.name = "doublet_prediction") <- doublets[[2]]

  cat("generating doublet detection QC... \n")
  glist <- list()
  for(value in QC_cols){
    fplot <- scatter(sc_obj[[]], x = "nCount_RNA", y = value, color = "doublet_prediction")
    glist <- c(glist, list(fplot))
  }
  complete_qc <- list(gridExtra::arrangeGrob(grobs = glist))
  names(complete_qc) <- c("Anomalies")
  graphs(sc_obj) <- complete_qc


  return(sc_obj)
}

#' rm_anomalies
#'
#' Convenience function which generate a new sc_obj assay containing the dataset
#' after anomalies have been removed. Raw data remains intact.
#'
#' @param sc_obj an sc_obj with attached QC.
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

#' colnames.unique
#' Takes a vector of meta data column names that contain duplicates
#'  and returns a vector of unique column names
#'
#'  This function relies on coercion of character values -> NA to detect
#'  names that have already been made unique and insure that they are not further
#'  modified.
#'
#' @param char_names
#'
#' @return
#' @export
#'
#' @examples
colnames.unique <- function(char_names){ #accepts character vector of column names
  char_list <- strsplit(char_names, "\\.")
  name_terminus <- sapply(char_list, FUN = function(x){
    x[length(x)]
  })
  suppressWarnings(remove_terminus <- !is.na(as.numeric(name_terminus)))
  base_names <- c()
  for(i in seq_along(remove_terminus)){
    if(remove_terminus[i]){
      base_names <- c(base_names,
                      paste(char_list[[i]][1:length(char_list[[i]])-1],
                            collapse = "."))
    }else{
      base_names <- c(base_names, char_names[i])
    }
  }
  return(make.unique(base_names))
}
