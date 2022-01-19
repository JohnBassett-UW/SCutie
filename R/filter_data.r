newDemultiplex <- function(data_10x){
  meta.data["nCount_RNA"] <- colSums(data_10x)
  meta.data["nfeature_RNA"] <- # number of non-zero values in cols of data_10x
  new(Class = "scExperiment", raw.data = data_10x, meta.data = )
}


remove1 <- function(char_cell_barcodes){
  sub("-1", "", char_cell_barcodes)
}

fraction_MT <- function(raw.data){
  #calculates the fraction of mito genes in cell
}
