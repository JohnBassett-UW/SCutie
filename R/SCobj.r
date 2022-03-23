SC_obj <- setClass(Class = "SC_obj",
                  slots = c(
                    raw.data = 'ANY',
                    assay = 'list',
                    meta.data = 'data.frame',
                    graphs = 'list',
                    ident = 'factor'
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
select <- function(x, meta.logical){ #this is wonky af
  y <- x@meta.data[x@meta.data[,meta.logical],]
  y <- rownames(y) #TODO make this less error prone
  return(x@raw.data[,y])
}

addGraph <-function(x, value){
  x@graphs <- c(x@graphs, list(value))
  return(x)
}

addMetaData <- function(x, value, col.name){
  x@meta.data[col.name] <- value
  return(x)
}

remove1 <- function(data.10x){
  dimnames(data.10x)[[2]] <- sub("-1", "", dimnames(data.10x)[[2]])
  return(data.10x)
}

initializeMetaData <- function(data.10x){ ##function for initializing metadata from count matrix for sc_obj
  meta.data<- data.frame(Matrix::colSums(data.10x))
  meta.data[,2] <- diff(data.10x@p)
  names(meta.data) <- c("nCount_RNA", "nFeatures_RNA")
  return(meta.data)
}

newSC_obj <- function(data.10x, remove1 = T){
  if(remove1 == T){
    data.10x <- remove1(data.10x)
  }
  new(Class = "SC_obj",
      raw.data = data.10x,
      meta.data = initializeMetaData(data.10x))
}
