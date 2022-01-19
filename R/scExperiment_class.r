scExperiment <- setClass(Class = "scExperiment",
                         slots = c(
                           raw.data = 'ANY',
                           assay = 'list',
                           meta.data = 'data.frame',
                           reductions = 'list',
                           ident = 'factor'
                          )
                         )

setGeneric("newDemultiplex", function(x) standardGeneric("newDemultiplex"))

setGeneric("meta.data", function(x) standardGeneric("meta.data"))
setGeneric("meta.data<-", function(x, value, col.name) standardGeneric("meta.data<-"))


setMethod("meta.data", "scExperiment", function(x) x@meta.data)
setMethod("meta.data<-", "scExperiment", function(x, value){
  cbind(x@meta.data, value)
  x@meta.data
})
