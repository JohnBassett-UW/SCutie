#' SC_obj
#'
#' @slot raw.data ANY.
#' @slot assay list.
#' @slot meta.data data.frame.
#' @slot graphs list.
#' @slot ident factor.
#' @slot est_doublet_rate numeric.
#'
#' @return
#' @export
#'
#' @examples
SC_obj <- setClass(Class = "SC_obj",
                  slots = c(
                    raw.data = 'ANY',
                    assay = 'list',
                    meta.data = 'data.frame',
                    graphs = 'list',
                    ident = 'factor',
                    est_doublet_rate = 'numeric'
                  )
)

setMethod( f = "[",
           signature = c("SC_obj"),
           definition = function(x, i, j, ...){
             if(missing(i)){
               if(missing(j)){
                 return(x@raw.data)
               }
               x@raw.data[,j]
             }
             else if(missing(j)){
               x@raw.data[i,]
             }
             else{
               return(x@raw.data[i,j])
             }
           })

setMethod( f = "[[",
           signature = c("SC_obj"),
           definition = function(x, i, j, ...){
             if(missing(i)){
               if(missing(j)){
                 return(x@meta.data)
               }
               x@meta.data[,j]
             }
             else if(missing(j)){
               x@meta.data[i,]
             }
             else{
               return(x@meta.data[i,j])
             }
           })

setMethod( f = "[[<-",
  signature = c("SC_obj"),
  definition = function(x,i,j,value) {
    if(missing(i)){
      if(missing(j)){
        x@meta.data <- value
        print("whoops")
      }
      print("here")
      x@meta.data[,j] <- value
    }
    else if(missing(j)){
      x@meta.data[i,] <- value
    }
    else{
      x@meta.data[i,j] <- value
    }
  })

setMethod(f = "dimnames",
          signature = c("SC_obj"),
          definition = function(x){
            return(dimnames(x@raw.data))
          })

###############################################################################
#' select
#'
#' @param x
#' @param meta.logical
#'
#' @return
#' @export
#'
#' @examples
select <- function(x, meta.logical){ #this is wonky af
  y <- x@meta.data[x@meta.data[,meta.logical],]
  y <- rownames(y) #TODO make this less error prone
  return(x@raw.data[,y])
}

#' addGraph
#'
#' @param x
#' @param value
#'
#' @return
#' @export
#'
#' @examples
addGraph <-function(x, value){
  x@graphs <- c(x@graphs, list(value))
  return(x)
}

#' addMetaData
#'
#' @param x
#' @param value
#' @param col.name
#'
#' @return
#' @export
#'
#' @examples
addMetaData <- function(x, value, col.name){
  x@meta.data[col.name] <- value
  return(x)
}

#' remove1
#'
#' @param data.10x
#'
#' @return
#' @export
#'
#' @examples
remove1 <- function(data.10x){
  dimnames(data.10x)[[2]] <- sub("-1", "", dimnames(data.10x)[[2]])
  return(data.10x)
}

#' initializeMetaData
#'
#' @param data.10x
#'
#' @return
#'
#'
#' @examples
initializeMetaData <- function(data.10x){ ##function for initializing metadata from count matrix for sc_obj
  meta.data<- data.frame(Matrix::colSums(data.10x))
  meta.data[,2] <- diff(data.10x@p)
  names(meta.data) <- c("nCount_RNA", "nFeatures_RNA")
  return(meta.data)
}

#' newSC_obj
#'
#' @param data.10x
#' @param remove1
#'
#' @return
#' @export
#'
#' @examples
newSC_obj <- function(data.10x, remove1 = T){
  if(remove1 == T){
    data.10x <- remove1(data.10x)
  }
  new(Class = "SC_obj",
      raw.data = data.10x,
      meta.data = initializeMetaData(data.10x))
}

#' set_doublet_rate
#'
#' @param sc_obj
#' @param doublet_rate
#'
#' @return
#' @export
#'
#' @examples
set_doublet_rate <- function(sc_obj, doublet_rate){
  sc_obj@est_doublet_rate = doublet_rate
  return(sc_obj)
}

#' doublet_rate
#'
#' @param sc_obj
#'
#' @return
#' @export
#'
#' @examples
doublet_rate <- function(sc_obj){
  return(sc_obj@est_doublet_rate)
}
