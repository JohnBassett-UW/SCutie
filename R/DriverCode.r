data.path <- "C:/Users/jbassett/Desktop/github/Data/UW36/filtered_feature_bc_matrix.h5"
counts <- import10x(data.path)
sce <- newSC_obj(counts)
dbl_params <- est_doublets(sce)
sc_obj <- sce
sce <- Perform_QC(sce, dbl_params = dbl_params)
