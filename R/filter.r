library(Seurat)
Seurat.data <- Read10X_h5("../Data/UW36/filtered_feature_bc_matrix.h5")
data <- import10x("../Data/UW36/filtered_feature_bc_matrix.h5")
data <- import10x("")

check_matrices(data, Seurat.data)

data@Dimnames[[1]][duplicated(data@Dimnames[[1]])]

which(data@Dimnames[[1]] == "TBCE")

percentMT <- function(data_10x){

}




dupes <- duplicated(dimnames(data_10x)[[1]])
char_duplicates <- dimnames(data_10x)[[1]][dupes]

gen_names_unique <- function(char_duplicates, count = 0){
  count = count + 1
  if(count > 1){

  }
  for(i in 1:length(char_duplicates)){
  char_duplicates[i] <- paste(char_duplicates[i], count, sep = ".")
  }

  if(any(duplicated(char_duplicates))){
    char_duplicates[duplicated(char_duplicates)] <- gen_names_unique(char_duplicates, count)
  }

  return(char_duplicates)
}
