data.path <- "C:/Users/jbassett/Desktop/github/Data/UW36/filtered_feature_bc_matrix.h5"
counts <- import10x(data.path)
sce <- newSC_obj(counts)
sce <- Attach_QC(sce)
sce <- detect_anomalies(sce)
sce <- rm_anomalies(sce)

