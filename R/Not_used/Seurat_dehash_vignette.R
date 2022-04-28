library(Seurat)

# Load in the UMI matrix
pbmc.umis <- readRDS("../Data/HTODemuxFiles/pbmc_umi_mtx.rds")

# For generating a hashtag count matrix from FASTQ files, please refer to
# https://github.com/Hoohm/CITE-seq-Count.  Load in the HTO count matrix
pbmc.htos <- readRDS("../Data/HTODemuxFiles/pbmc_hto_mtx.rds")

# Select cell barcodes detected by both RNA and HTO In the example datasets we have already
# filtered the cells for you, but perform this step for clarity.
joint.bcs <- intersect(colnames(pbmc.umis), colnames(pbmc.htos))

# Subset RNA and HTO counts by joint cell barcodes
pbmc.umis <- pbmc.umis[, joint.bcs]
pbmc.htos <- as.matrix(pbmc.htos[, joint.bcs])

# Confirm that the HTO have the correct names
rownames(pbmc.htos)

# Setup Seurat object
pbmc.hashtag <- CreateSeuratObject(counts = pbmc.umis)

# Normalize RNA data with log normalization
pbmc.hashtag <- NormalizeData(pbmc.hashtag)
# Find and scale variable features
pbmc.hashtag <- FindVariableFeatures(pbmc.hashtag, selection.method = "mean.var.plot")
pbmc.hashtag <- ScaleData(pbmc.hashtag, features = VariableFeatures(pbmc.hashtag))

# Add HTO data as a new assay independent from RNA
pbmc.hashtag[["HTO"]] <- CreateAssayObject(counts = pbmc.htos)
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO", normalization.method = "CLR")

# If you have a very large dataset we suggest using k_function = 'clara'. This is a k-medoid
# clustering function for large applications You can also play with additional parameters (see
# documentation for HTODemux()) to adjust the threshold for classification Here we are using
# the default settings
pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO", positive.quantile = 0.99)

# Global classification results
table(pbmc.hashtag$HTO_classification.global)

# Group cells based on the max HTO signal
Idents(pbmc.hashtag) <- "HTO_maxID"
RidgePlot(pbmc.hashtag, assay = "HTO", features = rownames(pbmc.hashtag[["HTO"]])[1:5], ncol = 2)

FeatureScatter(pbmc.hashtag, feature1 = "hto_HTO-A", feature2 = "hto_HTO-C")

Idents(pbmc.hashtag) <- "HTO_classification.global"
VlnPlot(pbmc.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# First, we will remove negative cells from the object
pbmc.hashtag.subset <- subset(pbmc.hashtag, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(pbmc.hashtag.subset) <- "HTO"
pbmc.hashtag.subset <- ScaleData(pbmc.hashtag.subset, features = rownames(pbmc.hashtag.subset),
                                 verbose = FALSE)
pbmc.hashtag.subset <- RunPCA(pbmc.hashtag.subset, features = rownames(pbmc.hashtag.subset), approx = FALSE)
pbmc.hashtag.subset <- RunTSNE(pbmc.hashtag.subset, dims = 1:8, perplexity = 100)
DimPlot(pbmc.hashtag.subset)

# To increase the efficiency of plotting, you can subsample cells using the num.cells argument
HTOHeatmap(pbmc.hashtag, assay = "HTO", ncells = 5000)

# Extract the singlets
pbmc.singlet <- subset(pbmc.hashtag, idents = "Singlet")

# Select the top 1000 most variable features
pbmc.singlet <- FindVariableFeatures(pbmc.singlet, selection.method = "mean.var.plot")

# Scaling RNA data, we only scale the variable features here for efficiency
pbmc.singlet <- ScaleData(pbmc.singlet, features = VariableFeatures(pbmc.singlet))

# Run PCA
pbmc.singlet <- RunPCA(pbmc.singlet, features = VariableFeatures(pbmc.singlet))

# We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
pbmc.singlet <- FindNeighbors(pbmc.singlet, reduction = "pca", dims = 1:10)
pbmc.singlet <- FindClusters(pbmc.singlet, resolution = 0.6, verbose = FALSE)
pbmc.singlet <- RunTSNE(pbmc.singlet, reduction = "pca", dims = 1:10)

# Projecting singlet identities on TSNE visualization
DimPlot(pbmc.singlet, group.by = "HTO_classification")
