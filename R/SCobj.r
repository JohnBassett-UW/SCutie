#' SCobj
#'
#' @slot raw.data ANY.
#' @slot HTO.counts ANY.
#' @slot assay list.
#' @slot meta.data data.frame.
#' @slot graphs list.
#' @slot ident factor.
#' @slot est_doublet_rate numeric.
#'
#' @return S4 object of class SC_obj
#' @export
SC_obj <- setClass(Class = "SC_obj",
                   contains = 'HTO_obj',
                   slots = c(
                     raw.data = 'ANY',
                     HTO.counts = 'ANY',
                     HTO.dmplex = 'HTO_obj',
                     assay = 'list',
                     meta.data = 'data.frame',
                     graphs = 'list',
                     ident = 'factor',
                     est_doublet_rate = 'numeric'
                   )
)

setGeneric("newSC_obj", function(GeX_counts, ...){
  standardGeneric("newSC_obj")})

setGeneric("metaData", function(x){
  standardGeneric("metaData")
})

setGeneric("metaData<-", function(x, ...){
  standardGeneric("metaData<-")
})

setGeneric('graphs', function(x){
  standardGeneric("graphs")
})

setGeneric('graphs<-', function(x, ...){
  standardGeneric("graphs<-")
})

setGeneric("est_doublet_rate", function(x){
  standardGeneric("est_doublet_rate")
})

setGeneric("est_doublet_rate<-", function(x, value){
  standardGeneric("est_doublet_rate<-")
})

setGeneric("Graph", function(x, ...){
  standardGeneric("Graph")
})

setGeneric("Assay", function(x, pos){
  standardGeneric("Assay")
})

setGeneric("Assay<-", function(x, assay.name, value){
  standardGeneric("Assay<-")
})

setGeneric("HTO_raw", function(x){
  standardGeneric("HTO_raw")
})

setGeneric("HTO.dmplex<-", function(x, value){
  standardGeneric("HTO.dmplex<-")
})

setGeneric("HTO.build", function(x){
  standardGeneric("HTO.build")
})

setGeneric("HTOs", function(x){
  standardGeneric("HTOs")
})

setMethod(f = "HTOs",
          signature = ("SC_obj"),
          definition = function(x){
            return(x@HTO.dmplex@HTO.matrix)
          })

setMethod(f = "HTO.dmplex<-",
          signature = signature(x = "SC_obj"),
          definition = function(x, value){
            slot(x, "HTO.dmplex") <- value
            return(x)
})

setMethod(f = "HTO.build",
          signature = signature("SC_obj"),
          definition = function(x){
            if(is.null(Assay(x, "anomalies_removed"))){
              stop("Failed to build HTO object: assay does not exist")
            }
            cells.whiteList <- colnames(Assay(x, "anomalies_removed"))
            filt.HTO <- HTO_raw(x)[,cells.whiteList]
            filt.HTO <- HTO.format(filt.HTO)
            CLR.HTO <- CLR(filt.HTO)
            UMAP.dims <- generate_UMAP(CLR.HTO)
            HTO.dmplex(x) <- new(Class = "HTO_obj",
                                 HTO.matrix = filt.HTO,
                                 HTO.CLR_norm = CLR.HTO,
                                 UMAP = UMAP.dims,
                                 graphs = plotHashes(CLR.HTO, UMAP.dims)
                                 )
            return(x)
          })


setMethod(f = "HTO_raw",
          signature = signature(x = "SC_obj"),
          definition = function(x){
            return(slot(x, "HTO.counts"))
          })

setMethod(f = "Assay",
          signature = signature(x = "SC_obj"),
          definition = function(x, pos){
            return(slot(x, "assay")[[pos]])
          })

setMethod(f = "Assay<-",
          signature = signature(x = "SC_obj"),
          definition = function(x, assay.name, value){
            slot(x, "assay")[[assay.name]] <- value
            return(x)
          })

#' Graph
#'
#' @param x SC_obj.
#'
#' @importFrom gridExtra grid.arrange
#'
#' @return gridExtra output for stored graphs
#'
#' @examples
#' Graph(sc_obj, "Quality Check")
setMethod(f = "Graph",
          signature = signature(x = "SC_obj"),
          definition = function(x, index){
            gridExtra::grid.arrange(slot(x, "graphs")[[index]])
          })

setMethod("initialize",
          signature= signature(.Object = "SC_obj"),
          function(.Object,
                   raw.data,
                   HTO.counts,
                   meta.data,
                   ...){
            if(!missing(raw.data)){
              .Object@raw.data <- raw.data
            }else(
              stop("SC_obj class requires a gene expression count matrix")
            )
            if(!missing(HTO.counts)){
              .Object@HTO.counts <- HTO.counts
            }

            .Object@meta.data <- meta.data

            return(.Object)
          })

setMethod(f = "[", signature = signature(x= "SC_obj"),
           definition = function(x, i, j, ...){
             if(missing(i)){
               if(missing(j)){
                 return(slot(x, "raw.data"))
               }
               return(slot(x, "raw.data")[,j])
             }
             else if(missing(j)){
               return(slot(x, "raw.data")[i,])
             }
             else{
               return(slot(x, "raw.data")[i,j])
             }
           })

setMethod(f = "[[",
           signature = signature(x = "SC_obj"),
           definition = function(x, i, j, ...){
             if(missing(i)){
               if(missing(j)){
                 return(x@meta.data)
               }
               return(x@meta.data[,j])
             }else if(missing(j)){
               if(is.numeric(i)){
                 return(x@meta.data[i,])
               }
               return(x@meta.data[,i])
             }
             else{
               return(x@meta.data[i,j])
             }
           })

setMethod(f = "[[<-",
           signature = signature(x = "SC_obj"),
           definition = function(x, i, j, value) {
            if(missing(i)){
              if(missing(j)){
                slot(x, "meta.data") <- value
              }
              slot(x, "meta.data")[,j] <- value
            }else if(missing(j)){
              if(is.numeric(i)){
                slot(x, "meta.data")[i,] <- value
              }
              slot(x, "meta.data")[,i] <- value
            }
            else{
              slot(x, "meta.data")[i,j] <- value
            }
            return(x)
          })

setMethod(f = "dimnames",
          signature = signature(x = "SC_obj"),
          definition = function(x){
            return(dimnames(x@raw.data))
          })


setMethod(f = "est_doublet_rate",
          signature = signature(x = "SC_obj"),
          definition = function(x){
            return(slot(x, "est_doublet_rate"))
          })

setMethod(f = "est_doublet_rate<-",
          signature = signature(x = "SC_obj", value = "numeric"),
          definition = function(x, value){
            slot(x, "est_doublet_rate") <- value
            return(x)
          })

#' newSC_obj
#' generates a single cell object.
#'
#' @param GeX_counts Takes a matrix of GeX counts in the form of a dgCMatrix.
#' @param HTO_counts Optionally accepts a matrix of HTO counts in the form of a dgCMatrix
#'
#' @return a single cell object with generated meta data.
#' @export
#'
#' @examples
#' newSC_obj(GeX_counts)
#' ##optionally##
#' newSC_obj(Gex_counts, HTO_counts)
setMethod(f = "newSC_obj",
          signature = signature(GeX_counts = "dgCMatrix"),
          definition = function(GeX_counts, HTO_counts){
            #SUBROUTINE#
            #function to remove trailing "-1" from cell ranger data
            remove1 <- function(GeX_counts){
              dimnames(GeX_counts)[[2]] <- sub("-1", "", dimnames(GeX_counts)[[2]])
              return(GeX_counts)
            }
            #SUBROUTINE#
            #function to initialize nCount_RNA and nFeaturs_RNA as meta data
            initializeMetaData <- function(GeX_counts){ ##function for initializing metadata from count matrix for sc_obj
              meta.data<- data.frame(Matrix::colSums(GeX_counts))
              meta.data[,2] <- diff(GeX_counts@p)
              names(meta.data) <- c("nCount_RNA", "nFeatures_RNA")
              return(meta.data)
            }

            #MAIN#
            #detection and removal of trailing "-1" from Gex_counts
            GeXnames <- unlist(
              strsplit(colnames(GeX_counts), "-"),
              use.names = T)
            if(all(GeXnames[seq.int(2, length(GeXnames),2)] == "1")){
              GeX_counts <- remove1(GeX_counts)
            }
            #detection and removal of trailing "-1" from HTO.counts
            HTOnames <- unlist(
              strsplit(colnames(HTO_counts), "-"),
              use.names = T)
            if(all(HTOnames[seq.int(2, length(HTOnames),2)] == "1")){
              HTO_counts <- remove1(HTO_counts)
            }
            #Check that GeX and HTO cell IDS Match
            if(!missing(HTO_counts)){
              ID.diffs <- match(colnames(GeX_counts), colnames(HTO_counts))
              if(!identical(colnames(GeX_counts),
                            colnames(HTO_counts)[ID.diffs])){
                message("HTO and GeX Cell ID's do not match")
                outMessage <- all.equal(colnames(GeX_counts), colnames(HTO_counts)[ID.diffs])
                outMessage <- gsub("current", "HTO_cell_IDs",
                                   gsub("target", "GeX_cell_IDs", outMessage))
                warning(outMessage)
                if(dim(GeX_counts)[2] > dim(HTO_counts)[2]){
                  message(dim(GeX_counts)[2]-dim(HTO_counts)[2], " cell Ids removed")
                  GeX_counts <- GeX_counts[,colnames(HTO_counts)]
                }
              }
            }#create new SC_obj
            new(Class = "SC_obj",
                raw.data = GeX_counts,
                meta.data = initializeMetaData(GeX_counts),
                HTO.counts = HTO_counts
            )
          })

setMethod(f = "metaData",
          signature = signature(x = "SC_obj"),
          definition = function(x){
            return(slot(x, "meta.data"))
          })

setMethod(f = "metaData<-",
          signature = signature(x = "SC_obj"),
          definition = function(x, col.name, value){
            slot(x, "meta.data")[,col.name] <- value
            return(x)
          })

setMethod(f= "graphs",
          signature = signature(x = "SC_obj"),
          definition = function(x){
            if(is.null(names(slot(x, "graphs")))){
              return(slot(x, "graphs"))
            }else{
              return(names(slot(x, "graphs")))
            }
          })

setMethod(f= "graphs<-",
          signature = signature(x ="SC_obj"),
          definition = function(x, value){
            slot(x, "graphs") <- append(slot(x, "graphs"), value)
            return(x)
          })
