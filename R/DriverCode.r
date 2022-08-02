HTO_path <- "Z:/paddison_lab/John/scData/UW44scx/HTO_lib_Trt(UW44scx02)/U44scx02_HTO/umi_count"
GeX_path <- "C:/Users/basse/Documents/UW/scData_local/UW44/UW44scx_02/outs/filtered_feature_bc_matrix.h5"

HTO_counts <- importX(HTO_path)
GeX_counts <- importX(GeX_path)
sce <- newSC_obj(GeX_counts, HTO_counts)
sce <- Attach_QC(sce)
Graph(sce, "Quality Check")
sce <- detect_anomalies(sce, method = 'cxds')
Graph(sce, "Anomalies")
sce <- rm_anomalies(sce)


#Import and qc filter
QC_pipeline <- function(data.path){
  ptm = proc.time()
  counts <- import10x(data.path)
  sce <- newSC_obj(GeX_counts, HTO_counts)
  sce <- Attach_QC(sce)
  sce <- detect_anomalies(sce)
  sce <- rm_anomalies(sce)

  cat("Total \n time ")
  print((proc.time() - ptm)[3]) #print time elapsed
  return(sce)
}

#sc_obj to Seurat object
sc_to_Seurat <- function(sc_obj){
  Seurat_Object <- CreateSeuratObject(counts = sc_obj@assay[[1]])
  return(Seurat_Object)
}

#ccSeurat
ccSeurat <- function(){
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes

  obj_seurat <- NormalizeData(obj_seurat)
  obj_seurat <- FindVariableFeatures(obj_seurat, selection.method = "vst")
  obj_seurat <- ScaleData(obj_seurat, features = rownames(obj_seurat))
  obj_seurat <- RunPCA(obj_seurat, features = VariableFeatures(obj_seurat))

  DimHeatmap(obj_seurat)

  obj_seurat <- CellCycleScoring(obj_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  # view cell cycle scores and phase assignments
  head(obj_seurat[[]])
  # Visualize the distribution of cell cycle markers across
  RidgePlot(obj_seurat, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

  obj_seurat <- RunPCA(obj_seurat, features = c(s.genes, g2m.genes))
  DimPlot(obj_seurat)
}


#sce <- pipeline(data.path)
#obj_seurat <- sc_to_Seurat(sce)
gridExtra::grid.arrange(sce@graphs[[1]])


#counts <- import10x(data.path)
#sce2 <- newSC_obj(dat)
#sce2 <- Attach_QC(sce2)
#sce2 <- detect_anomalies(sce2)
#sce2 <- rm_anomalies(sce2)

