#### Load packages ----
library(tidyverse)
library(Seurat)
library(patchwork)
library(qs)

#### Load data ----
# WT
wt_counts <- Read10X(data.dir = "./data/WT/")
colnames(wt_counts) <- paste("WT", colnames(wt_counts), sep = "_")

wt_seurat_obj <-
  CreateSeuratObject(
    counts = counts,
    project = "WT",
    min.cells = 3,
    min.features = 200
  )

# KO
ko_counts <- Read10X(data.dir = "./data/KO/")
colnames(ko_counts) <- paste("KO", colnames(ko_counts), sep = "_")

ko_seurat_obj <-
  CreateSeuratObject(
    counts = ko_counts,
    project = "KO",
    min.cells = 3,
    min.features = 200
  )

merged <- merge(wt_seurat_obj, ko_seurat_obj)

#### QC and selecting cells for further analysis ----
merged[["percent.mt"]] <-
  PercentageFeatureSet(merged, pattern = "^mt-")

VlnPlot(
  merged,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

plot1 <-
  FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <-
  FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

merged <-
  subset(merged, subset = nFeature_RNA > 200 &
    nFeature_RNA < 5000 & percent.mt < 7.5)

#### Normalizing the data ----
merged <-
  NormalizeData(merged,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )

#### Identification of highly variable features (feature selection) ----
merged <-
  FindVariableFeatures(merged,
    selection.method = "vst",
    nfeatures = 2500
  )

#### Scaling the data ----
merged <- ScaleData(merged, vars.to.regress = "percent.mt")

#### Perform linear dimensional reduction ----
merged <-
  RunPCA(merged,
    features = VariableFeatures(object = merged)
  )

ElbowPlot(merged, ndims = 30)

#### Cluster the cells ----
merged <- FindNeighbors(merged, dims = 1:30)
merged <- FindClusters(merged, resolution = 0.7)
head(Idents(merged), 5)

#### Run non-linear dimensional reduction (UMAP/tSNE) ----
merged <- RunTSNE(merged, dims = 1:30)

DimPlot(merged,
  reduction = "tsne",
  label = TRUE,
  pt.size = 1.5
)

plot1 <-
  DimPlot(
    subset(merged, subset = orig.ident == "WT"),
    reduction = "tsne",
    label = TRUE,
    pt.size = 1.5
  )
plot2 <-
  DimPlot(
    subset(merged, subset = orig.ident == "KO"),
    reduction = "tsne",
    label = TRUE,
    pt.size = 1.5
  )
plot1 + plot2

DimPlot(
  merged,
  reduction = "tsne",
  label = TRUE,
  group.by = "orig.ident",
  pt.size = 1.5
)

qsave(merged, "./output/seurat/merged.qs")

#### Finding differentially expressed features (cluster biomarkers) ----
markers <-
  FindAllMarkers(
    merged,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
write_csv(markers, "./output/seurat/marker.csv")

#### Assigning cell type identity to clusters ----
Idents(merged) <- factor(merged@meta.data$seurat_clusters)
new.cluster.ids <- c(
  "naïve T cells",
  "M1 macrophage",
  "NK cells",
  "M2 macrophage",
  "CD8 effector",
  "M1 macrophage",
  "M2 macrophage",
  "cDCs",
  "B cells",
  "Treg",
  "Monocytes",
  "pDCs",
  "NK cells",
  "Ki67+ CD8 T cells",
  "Neutrophils",
  "gdT cells",
  "M1 macrophage",
  "Migratory DCs",
  "Mast cells"
)

names(new.cluster.ids) <- levels(merged)
merged <- RenameIdents(merged, new.cluster.ids)
DimPlot(merged, label = T) + NoLegend()
qsave(merged, "./output/seurat/merged_anno.qs")
