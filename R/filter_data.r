violin <- function(data, y){
  vplot <- ggplot2::ggplot(data, ggplot2::aes_string(y = y, x = 0)) +
    ggplot2::geom_violin(ggplot2::aes(fill = "darkred"), show.legend = F) +
    ggplot2::geom_jitter(size = 0.5)
  return(vplot)
}

scatter<- function(data, x, y){
  Splot <- ggplot2::ggplot(data, ggplot2::aes_string(x = x, y = y)) + ggplot2::geom_point() #attach plot to obj
  return(Splot)
}

Perform_QC <- function(sc_obj, pattern = "default"){
  if(pattern == "default"){
    pattern = c("^MT-", "^RP")
  }

  glist <- list()

  for(p in pattern){ #iterate through patterns and compare grouped transcripts as fractions of cell nCount totals
    fgenes <- grep(dimnames(sc_obj)[[1]], pattern = p)
    fcounts <- colSums(sc_obj[fgenes,])
    fpercent <- fcounts/sc_obj[["nCount_RNA"]] * 100
    col.name <- paste("percent", gsub(pattern = "[^[:alnum:]]", "", p), sep = "")
    sc_obj <- addMetaData(sc_obj, fpercent, col.name = col.name)

    fplot <- ggplot2::ggplot(sc_obj[[]], ggplot2::aes_string(x = "nCount_RNA", y = col.name)) + ggplot2::geom_point() #TODO move to own function
    glist <- c(glist, list(fplot))
  }

  nCount_qc <- violin(data = sc_obj[[]], y = "nCount_RNA" )
  nfeature_qc <- violin(data = sc_obj[[]], y = "nFeatures_RNA" )

  glist<- c(glist, list(nCount_qc, nfeature_qc))

  complete_qc <- gridExtra::grid.arrange(grobs = glist )

  sc_obj <- addGraph(x = sc_obj, value = complete_qc)

  return(sc_obj)
} #TODO remove plots from perfom qc
#Add new function for plotting graphs
#use grid arrange to compile all initial qc plots in a single instance
#Add all initial graphs to graphs category in one instance


QC_filter <- function(SC_obj){ #TODO filter data based on user input parameters

}

sc_obj <- newSC_obj(data.10x)
sc_obj <- Perform_QC(sc_obj)
View(sc_obj)



