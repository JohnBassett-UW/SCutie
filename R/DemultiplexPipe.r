cells.whiteList <- colnames(Assay(sce, "anomalies_removed"))
filt.HTO <- HTO_counts(sce)[,cells.whiteList]
Assay(sce, "Demultiplex_HTOs") <- data.format(filt.HTO)
bar.UMAP <- barUMAP(Assay(sce, "Demultiplex_HTOs"))
plotHashes(bar.UMAP)
