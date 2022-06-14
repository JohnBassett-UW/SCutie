#' est_doublets
#' dispatches count matrix to scds to estimate doublet rate
#' refer to https://github.com/kostkalab/scds
#' doi: 10.1093/bioinformatics/btz698
#'
#' @param counts single cell RNAseq count matrix where columns are cell ID's and rows are genes
#' @param counts alternately accepts sc_obj
#' @param counts alternately accepts SingleCellExperiment Object
#' @param method defaults to "hybrid". Alternately "cxds" or "bcds" can be used
#' @importFrom scds cxds bcds cxds_bcds_hybrid
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @return estimate of the proportion of doublets
#' @export
#'
#' @examples
#' est_doublets allows you to chose any of the three methods used to score doublets
#' provided by Bais et al.
#'
#' est_doublets(counts, method = "hybrid")
#'
est_doublets <- function(counts, method = "hybrid"){
  if(inherits(counts, "SC_obj")){
    counts <- counts[]
  }
  if(!inherits(counts, "SingleCellExperiment")){
    sce = SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  }
  else{
    sce = counts
  }

  switch(method,
         "cxds" = {
           sce = scds::cxds(sce, estNdbl =TRUE)
           fDoublets = sce@metadata$cxds$ndbl["balanced",][3]
           calls = sce$cxds_call
         },
         "bcds" = {
           sce = scds::cxds(sce, estNdbl =TRUE)
           fDoublets = sce@metadata$bcds$ndbl["balanced",][3]
           calls = sce$bcds_call
         },
         "hybrid" = {
           sce = scds::cxds_bcds_hybrid(sce, estNdbl = TRUE)
           fDoublets = sce@metadata$hybrid$ndbl["balanced",][3]
           calls = sce$hybrid_call
         }
        )
  cat("approximate doublet rate: ", fDoublets, "\n")
  return(list(as.numeric(fDoublets), calls))
}
