
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)
library(glmGamPoi)

# s03
s03.data <- Read10X(data.dir = ".../s03/04.Matrix", gene.column=1)
s03 <- CreateSeuratObject(counts = s03.data, project = "s03", min.cells = 3, min.features = 200)
dim(s03)
# [1] 23724  2852
s03[["percent.mt"]] <- PercentageFeatureSet(s03, pattern = "^MT-")
p1 <- VlnPlot(s03, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.05)
s03_v1 <- subset(s03, subset = nFeature_RNA < 9000 & percent.mt < 10)
dim(s03_v1)
# [1] 23724  2654
p2 <- VlnPlot(s03_v1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.05)

# s08
s08.data <- Read10X(data.dir = ".../s08/04.Matrix", gene.column=1)
s08 <- CreateSeuratObject(counts = s08.data, project = "s08", min.cells = 3, min.features = 200)
dim(s08)
# [1] 26481  3872
s08[["percent.mt"]] <- PercentageFeatureSet(s08, pattern = "^MT-")
p3 <- VlnPlot(s08, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.05)
s08_v1 <- subset(s08, subset = nFeature_RNA < 9000 & percent.mt < 10)
dim(s08_v1)
# [1] 26481  3397
p4 <- VlnPlot(s08_v1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0.05)


# Perform normalization and dimensionality reduction
s03_sct <- SCTransform(s03_v1, vst.flavor = "v2", verbose = FALSE)
s08_sct <- SCTransform(s08_v1, vst.flavor = "v2", verbose = FALSE)

sct.list <- list(s03 = s03_sct, s08 = s08_sct)
features <- SelectIntegrationFeatures(object.list = sct.list , nfeatures = 3000)
sct.list  <- PrepSCTIntegration(object.list = sct.list, anchor.features = features)

sct.anchors <- FindIntegrationAnchors(object.list = sct.list, normalization.method = "SCT",
                                      anchor.features = features)
s38.combined.sct <- IntegrateData(anchorset = sct.anchors, normalization.method = "SCT")

# Perform an integrated analysis
DefaultAssay(s38.combined.sct)
# [1] "integrated"
s38.combined.sct <- RunPCA(s38.combined.sct)

# dims=40
s38.combined.sct <- RunUMAP(s38.combined.sct, reduction = "pca", dims = 1:40)
s38.combined.sct <- FindNeighbors(s38.combined.sct, reduction = "pca", dims = 1:40)
s38.combined.sct <- FindClusters(s38.combined.sct, resolution = seq(0.5,1.2,by = 0.1))
s38.combined.sct$clusteronce <- s38.combined.sct$integrated_snn_res.0.8

# 
DefaultAssay(s38.combined.sct) <- "SCT"
Idents(s38.combined.sct) <- s38.combined.sct$clusteronce

s38.combined.sct <- PrepSCTFindMarkers(s38.combined.sct)
All.marker <- FindAllMarkers(s38.combined.sct, assay = "SCT", slot = "data")
Idents(s38.combined.sct) <- s38.combined.sct$clusteronce

s38.combined.sct <- RenameIdents(s38.combined.sct,
                                 "0" = "α cells",
                                 "1" = "α cells",
                                 "2" = "β cells",
                                 "3" = "β cells",
                                 "4" = "β cells",
                                 "5" = "EC cells",
                                 "6" = "Pancreatic progenitor",
                                 "7" = "δ cells",
                                 "8" = "α cells",#ARX+IRX2+
                                 "9" = "Polyhormonal cells",#ARX+IRX2+INS+
                                 "10" = "Proliferation cells",
                                 "11" = "ε cells")
s38.combined.sct$cell_type <- Idents(s38.combined.sct)

saveRDS(s38.combined.sct, file = ".../s38.combined.sct_afterRename.rds")