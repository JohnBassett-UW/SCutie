#' HTO_obj
#'
#' @slot HTO.matrix ANY.
#' @slot HTO.CLR_norm ANY.
#' @slot UMAP matrix.
#' @slot meta.data data.frame.
#' @slot graphs list.
#'
#' @export
#'
HTO_obj <- setClass(Class = "HTO_obj",
                   slots = c(
                     HTO.matrix = 'ANY', #contains HTOcounts after filtering out outliers from sc_data
                     HTO.CLR_norm = 'ANY', #clr normalized filtered HTO counts
                     UMAP = 'matrix',
                     classifications = 'data.frame',
                     statistics = 'data.frame',
                     graphs = 'list'
                   )
)

#' HTO.format
#'
#' Format function equipped to handle the outputs of CITE-seq count or the MULTIseq Alignment Suite
#'
#' HTO.format() performs the following manipulations:
#' casts data input to a matrix using as.matrix()
#' transposes the input data if necessary
#' removes any uninformative columns attached during the alignment
#' @param UMI_count_matrix raw count output
#'
#' @return formatted count matrix where each column corresponds to a barcode and each row to a cell UMI
#' @export
#'
#' @examples
#' f.bar.table <- data.format(bar.table)
HTO.format <- function(UMI_count_matrix){
  #Check if input is type matrix array
  if(!is(UMI_count_matrix, "matrix")){
    cat("Attempting to cast object type ", class(UMI_count_matrix), " to matrix \n")
    #Attempt to cast to sparse matrix
    UMI_count_matrix <- tryCatch(as.matrix(UMI_count_matrix),

                                 error = function(e){
                                   message("Object ", class(UMI_count_matrix)," could not be cast to matrix: \n", e)
                                   stop()
                                 },

                                 warning = function(w){
                                   message("Warnings detected during casting to matrix\n", w)
                                 },

                                 finally = {
                                   cat("Succeeded \n")
                                 })
  }
  else{message("input class ", class(UMI_count_matrix))}
  #
  #assumes there are fewer hashes than UMIs and transposes matrix so that UMIs are expressed as rows
  if( ncol(UMI_count_matrix) > nrow(UMI_count_matrix) ){
    UMI_count_matrix <- t(UMI_count_matrix)
    cat("input matrix transposed \n")
  }
  #
  #removes uninformative columns
  columns.remove <- c("nUMI", "nUMI_total", "unmapped") #list of known uninformative column names
  columns.retain <- setdiff(colnames(UMI_count_matrix), columns.remove)
  columns.remove <- intersect(colnames(UMI_count_matrix), columns.remove) #columns.remove now reflects columns which will be removed
  UMI_count_matrix <- UMI_count_matrix[,columns.retain]

  if(length(columns.remove) == 0 ){
    cat("No bad columns detected \n")}
  else if(length(columns.remove) == 1){
    cat("one uninformative column removed: ", paste(columns.remove, collapse = " "), "\n")}
  else{
    cat(length(columns.remove), " uninformative columns removed: ", paste(columns.remove, collapse = " "), "\n")
  }

  cat("rows: ", nrow(UMI_count_matrix), ", columns: ", ncol(UMI_count_matrix), " \n")
  return(UMI_count_matrix)
}

#'Perform center logratio transform on the input data
#'
#'The CLR function provided is nearly identical to the one used in the Seurat package
#'except log1p() has been replaced with log(). The rational for this change is
#'that valid barcode counts are always positive and substantially greater than 0.
#'
#' @param bar.table formatted sample barcode table
#'
#' @return center logratio transformed bar.table as class Matrixarray
#' @export
#'
#' @examples
#' n.bar.table <- CLR(bar.table)
CLR <- function(bar.table) {
  bar.table <- apply(bar.table, MARGIN = 2, FUN = function(x){
    log(x=x/exp(x=sum(log(x=x[x>0]), na.rm = TRUE)/length(x=x)))
  })
  bar.table[is.infinite(bar.table[])==TRUE] = 0
  return(bar.table)
}

#' generate_UMAP
#' Performs dimension reduction on a formatted sample barcode table
#'
#' @param bar.table formatted sample barcode table
#' @param normalize logical, should the return table contain CLR normalized values. defaults to False.
#' @param n_neighbors n neighbors parameter for the umap function from package uwot
#' @param min_dist minimum distance parameter for the umap function from package uwot
#'
#' @importFrom uwot umap
#'
#' @return sample barcode table with UMAP dimensions UMAP1 & UMAP2 inserted as first two columns.
#' @export
#'
#' @examples
#' bar.UMAP <- generate_UMAP(bar.table, n_neighbors =50, min_dist = 0.1)
#' bar.UMAP <- generate_umap(bar.table)
generate_UMAP <- function(n.bar.table, n_neighbors = 50, min_dist = 0.1){

  #remove infinite values
  n.bar.table[is.infinite(n.bar.table[])==TRUE] = 0
  #Generate UMAP via uwot package using default parameters
  cat("performing dimension reduction \n")
  cat("this may take a minute... \n")
  UMAP.res <- uwot::umap(n.bar.table,
                         n_neighbors = n_neighbors, #Size of Local neighborhood to constrain manifold learning
                         n_components = 2, #dimensions to embed to
                         n_epochs = NULL, #default to 500 for datasets containing >= 10k vertices. 200 otherwise.
                         min_dist = min_dist, #minimum distance apart points are allowed to be
                         n_trees = 50) #number of trees to build when constructing nearest neighbors. sensible values are 10-100. larger -> better

  colnames(UMAP.res) <- c("UMAP1", "UMAP2")

  return(UMAP.res)
}


#' Covenience function for visualizing output of barUMAP()
#'
#' @param bar.UMAP output matrix from barUMAP() function
#' @param plt logical, should plots be printed to std out. defaults to True.
#' Otherwise ggplot gobs are returned as a list.
#'
#' @importFrom ggplot2 geom_point theme_classic scale_color_gradient theme ggtitle
#' @importFrom gridExtra grid.arrange
#'
#' @return If plt is FALSE returns list of ggplot grobs otherwise no return value
#' @export
#'
#' @examples
#' plotHashes(bar.UMAP)
plotHashes <- function(n.bar.table, UMAP.dims){

  #ggplot requires data.frame
  UMAP.dims <- tryCatch(as.data.frame(UMAP.dims),
                       error = function(e){
                         stop("ggplot requires type data.frame")
                       })
  #plotHashes function requires normalized UMAP
  bars <- tryCatch(as.data.frame(n.bar.table),
                   error = function(e){
                     stop("Encountered error evaluating Normalization status")
                   },
                   warning = function(w){
                     message("Normalized data detected")
                     bars
                   }#,
                   # finally = function(f){
                   #   message("Good status")
                   # }
  )

  #Remove all values below the column geometric means
  bars[bars<0] <- NaN

  #bars <- as.data.frame(bars)
  cat(ncol(bars)," barcodes \n")

  ##populate plist with barcode columns from bar.UMAP
  plist <- as.list(bars)
  #iterate over barcode columns to create list of UMAPs for each barcode with ggplot
  counter <- 0
  plist <- lapply(plist, function(x){
    counter <<- counter + 1
    ggplot2::ggplot(UMAP.dims, ggplot2::aes(x=UMAP1, y=UMAP2, col=x)) +
      ggplot2::geom_point(size=0.5) +
      ggplot2::theme_classic() +
      ggplot2::scale_color_gradient(low="blue", high="red") +
      ggplot2::theme(legend.position = "none") +
      ggplot2::ggtitle(names(bars[counter]))
  })

  cat("generating, plots... \n")
  rowsArranged = max(1, ncol(bars)%/%8)
  complete_UMAPS <- gridExtra::grid.arrange(grobs = plist , nrow=rowsArranged , ncol=min(c(ncol(bars),8))) #arrange UMAPs side by side
  complete_UMAPS <- list(complete_UMAPS)
  names(complete_UMAPS) <- "HTO_QC"
  cat("Done. \n [*_*] \n")
  return(complete_UMAPS)
}

