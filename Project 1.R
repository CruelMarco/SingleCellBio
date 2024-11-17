## Week 1

#Installing libraries
library(Seurat)
library(dplyr)
library(spatstat.core)
library(Seurat)
library(patchwork)
library(DoubletFinder)
library(SingleR)
library(enrichR)
library(CellChat)
library(SingleCellExperiment)
library(SeuratWrappers)
library(tidyverse)
library(celldex)
library(ggplot2)
# library(monocle3)
set.seed(42)

#setting directory - need to write code to install dataset from website if not in directory to generalise it

setwd("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1")

#create object and add metadata
BMMC_D1T1 <- readRDS("~/Desktop/WiSe24/singlecell/Assignment1/scbi_ds1/GSM4138872_scRNA_BMMC_D1T1.rds")
BMMC_D1T2 <- readRDS("~/Desktop/WiSe24/singlecell/Assignment1/scbi_ds1/GSM4138873_scRNA_BMMC_D1T2.rds")
CD34_D2T1 <- readRDS("~/Desktop/WiSe24/singlecell/Assignment1/scbi_ds1/GSM4138874_scRNA_CD34_D2T1.rds")
CD34_D3T1 <- readRDS("~/Desktop/WiSe24/singlecell/Assignment1/scbi_ds1/GSM4138875_scRNA_CD34_D3T1.rds")

BMMC_D1T1 <- CreateSeuratObject(counts = BMMC_D1T1, project = "BMMC_D1T1")
BMMC_D1T2 <- CreateSeuratObject(counts = BMMC_D1T2, project = "BMMC_D1T2")
CD34_D2T1 <- CreateSeuratObject(counts = CD34_D2T1, project = "CD34_D2T1")
CD34_D3T1 <- CreateSeuratObject(counts = CD34_D3T1, project = "CD34_D3T1")

metadata <- data.frame(
  Sample = c("BMMC_D1T1", "BMMC_D1T2", "CD34_D2T1", "CD34_D3T1"),
  Donor = c("D1", "D1", "D2", "D3"),
  Replicate = c("T1", "T2", "T1", "T1"),
  Sex = c("F", "F", "M", "F")
)

BMMC_D1T1$Sample <- "BMMC_D1T1"
BMMC_D1T1$Donor <- "D1"
BMMC_D1T1$Replicate <- "T1"
BMMC_D1T1$Sex <- "F"

BMMC_D1T2$Sample <- "BMMC_D1T2"
BMMC_D1T2$Donor <- "D1"
BMMC_D1T2$Replicate <- "T2"
BMMC_D1T2$Sex <- "F"

CD34_D2T1$Sample <- "CD34_D2T1"
CD34_D2T1$Donor <- "D2"
CD34_D2T1$Replicate <- "T1"
CD34_D2T1$Sex <- "M"

CD34_D3T1$Sample <- "CD34_D3T1"
CD34_D3T1$Donor <- "D3"
CD34_D3T1$Replicate <- "T1"
CD34_D3T1$Sex <- "F"

#Summary Information
cell_num_bmmc1 <- ncol(BMMC_D1T1)
cell_num_bmmc2 <- ncol(BMMC_D1T2)
cell_num_cd341 <- ncol(CD34_D2T1)
cell_num_cd342 <- ncol(CD34_D3T1)

print(cell_num_bmmc1)
print(cell_num_bmmc2)
print(cell_num_cd341)
print(cell_num_cd342)

cell_gene_bmmc1 <- nrow(BMMC_D1T1)
cell_gene_bmmc2 <- nrow(BMMC_D1T2)
cell_gene_cd341 <- nrow(CD34_D2T1)
cell_gene_cd342 <- nrow(CD34_D3T1)

print(cell_gene_bmmc1)
print(cell_gene_bmmc2)
print(cell_gene_cd341)
print(cell_gene_cd342)

print(colnames(BMMC_D1T1@meta.data))
print(colnames(BMMC_D1T2@meta.data))
print(colnames(CD34_D2T1@meta.data))
print(colnames(CD34_D3T1@meta.data))

saveRDS(BMMC_D1T1, file="/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/bmmc_d1t1.rds")
saveRDS(BMMC_D1T2, file="/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/bmmc_d1t2.rds")
saveRDS(CD34_D2T1, file="/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/cd34_d2t1.rds")
saveRDS(CD34_D3T1, file="/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/cd34_d3t1.rds")

## Week 2

#quality control
BMMC_D1T1[["percent.mt"]] <- PercentageFeatureSet(BMMC_D1T1, pattern = "^MT-")
BMMC_D1T2[["percent.mt"]] <- PercentageFeatureSet(BMMC_D1T2, pattern = "^MT-")
CD34_D2T1[["percent.mt"]] <- PercentageFeatureSet(CD34_D2T1, pattern = "^MT-")
CD34_D3T1[["percent.mt"]] <- PercentageFeatureSet(CD34_D3T1, pattern = "^MT-")

QC_violin_plot <- VlnPlot(
  BMMC_D1T1,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
)

ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/QC_violin_plot_BMMC1.png", plot = QC_violin_plot)

QC_violin_plot <- VlnPlot(
  BMMC_D1T2,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
)

ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/QC_violin_plot_BMMC2.png", plot = QC_violin_plot)

QC_violin_plot <- VlnPlot(
  CD34_D2T1,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
)

ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/QC_violin_plot_CD341.png", plot = QC_violin_plot)

QC_violin_plot <- VlnPlot(
  CD34_D3T1,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
)

ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/QC_violin_plot_CD342.png", plot = QC_violin_plot)

#Filtering

BMMC_D1T1 <- subset(BMMC_D1T1, 
                       subset = nFeature_RNA >= 200 & nFeature_RNA <= 4000 &
                         nCount_RNA < 10000)
BMMC_D1T2 <- subset(BMMC_D1T2, 
                    subset = nFeature_RNA >= 200 & nFeature_RNA <= 4000 &
                      nCount_RNA < 10000)
CD34_D2T1 <- subset(CD34_D2T1, 
                    subset = nFeature_RNA <= 5000 &
                      nCount_RNA < 12000)
CD34_D3T1 <- subset(CD34_D3T1, 
                    subset = nFeature_RNA <= 4000 &
                      nCount_RNA < 10000)

BMMC_D1T1 <- subset(BMMC_D1T1, 
                    subset = nFeature_RNA >= 200 & nFeature_RNA <= 4000 &
                      nCount_RNA < 10000)
BMMC_D1T2 <- subset(BMMC_D1T2, 
                    subset = nFeature_RNA >= 200 & nFeature_RNA <= 4000 &
                      nCount_RNA < 10000)
CD34_D2T1 <- subset(CD34_D2T1, 
                    subset = nFeature_RNA <= 5000 &
                      nCount_RNA < 12000)
CD34_D3T1 <- subset(CD34_D3T1, 
                    subset = nFeature_RNA <= 4000 &
                      nCount_RNA < 10000)

#pipeline for preprocessing steps
BMMC_D1T1 <- BMMC_D1T1 %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(features = VariableFeatures(object = BMMC_D1T1))

BMMC_D1T2 <- BMMC_D1T2 %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(features = VariableFeatures(object = BMMC_D1T2))

CD34_D2T1 <- CD34_D3T1 %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(features = VariableFeatures(object = CD34_D2T1))

CD34_D3T1 <- CD34_D3T1 %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>%
  RunPCA(features = VariableFeatures(object = CD34_D3T1))

bmmc1_elbow <- ElbowPlot(BMMC_D1T1, ndims = 50) +
                         theme_minimal() +
                         geom_vline(xintercept = 20, linetype = "dashed", color = "black", size = 0.5)
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/bmmc1_elbow.png", plot = bmmc1_elbow, width = 6, height = 4)

bmmc2_elbow <- ElbowPlot(BMMC_D1T2, ndims = 50) +
                         theme_minimal() +
                         geom_vline(xintercept = 20, linetype = "dashed", color = "black", size = 0.5)
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/bmmc2_elbow.png", plot = bmmc2_elbow, width = 6, height = 4)

cd341_elbow <- ElbowPlot(CD34_D2T1, ndims = 50) +
                         theme_minimal() +
                         geom_vline(xintercept = 20, linetype = "dashed", color = "black", size = 0.5)
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/cd341_elbow.png", plot = cd341_elbow, width = 6, height = 4)

cd342_elbow <- ElbowPlot(CD34_D3T1, ndims = 50) +
                         theme_minimal() +
                         geom_vline(xintercept = 20, linetype = "dashed", color = "black", size = 0.5)
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/cd342_elbow.png", plot = cd342_elbow, width = 6, height = 4)

bmmc1_dim <- DimPlot(BMMC_D1T1, reduction = "pca", group.by = "orig.ident") +
  theme_minimal()
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/bmmc1_dim.png", plot = bmmc1_dim, width = 6, height = 4)

bmmc2_dim <- DimPlot(BMMC_D1T2, reduction = "pca", group.by = "orig.ident") +
  theme_minimal()
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/bmmc2_dim.png", plot = bmmc2_dim, width = 6, height = 4)

cd341_dim <- DimPlot(CD34_D2T1, reduction = "pca", group.by = "orig.ident") +
  theme_minimal()
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/cd341_dim.png", plot = cd341_dim, width = 6, height = 4)

cd342_dim <- DimPlot(CD34_D3T1, reduction = "pca", group.by = "orig.ident") +
  theme_minimal()
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/cd342_dim.png", plot = cd342_dim, width = 6, height = 4)

#we select 20 features for each because that's the x-value after which the dots become stable overall. it helps compare easily too.

#DoubletFinder

sweep.res.list_bmmc1 <- paramSweep(BMMC_D1T1, PCs = 1:20, sct = FALSE)
sweep.stats_bmmc1 <- summarizeSweep(sweep.res.list_bmmc1, GT = FALSE)
bcmvn_bmmc1 <- find.pK(sweep.stats_bmmc1)
optimal_pK_bmmc1 <- as.numeric(as.character(bcmvn_bmmc1$pK[which.max(bcmvn_bmmc1$BCmetric)]))
pk_bmmc1 <- ggplot(bcmvn_bmmc1, aes(x = pK, y = BCmetric)) +
  geom_line(group = 1, color = "lightblue", size = 1) +
  geom_point(size = 2, color = "blue") +
  annotate("text", x = optimal_pK_bmmc1, y = max(bcmvn_bmmc1$BCmetric),
           label = paste("Optimal pK:", optimal_pK_bmmc1), hjust = -0.1, vjust = -0.5) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white", color = NA)) +
  labs(title = "Optimal pK Selection for BMMC_D1T1", x = "pK", y = "BC Metric")
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/pK_select_bmmc1.png", plot = pk_bmmc1, width = 12, height = 6)

doublet_rate_for_all <- 0.075 #7.5% assumed
nExp_bmmc1 <- round(ncol(BMMC_D1T1) * doublet_rate_for_all)
cat("Expected doublets number:", nExp_bmmc1)

BMMC_D1T1 <- doubletFinder(
      BMMC_D1T1,
      PCs = 1:20,
      pN = 0.25,
      pK = optimal_pK_bmmc1,
      nExp = nExp_bmmc1
    )

doublet_dim_bmmc1 <- DimPlot(BMMC_D1T1, group.by = "DF.classifications_0.25_0.27_470") +
  labs(title = "Doublet Classification") +
  theme_minimal()
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/doublet_dim_bmmc1.png", plot = doublet_dim_bmmc1, width = 12, height = 6)

sweep.res.list_bmmc2 <- paramSweep(BMMC_D1T2, PCs = 1:20, sct = FALSE)
sweep.stats_bmmc2 <- summarizeSweep(sweep.res.list_bmmc2, GT = FALSE)
bcmvn_bmmc2 <- find.pK(sweep.stats_bmmc2)
optimal_pK_bmmc2 <- as.numeric(as.character(bcmvn_bmmc2$pK[which.max(bcmvn_bmmc2$BCmetric)]))
pk_bmmc2 <- ggplot(bcmvn_bmmc2, aes(x = pK, y = BCmetric)) +
  geom_line(group = 1, color = "lightblue", size = 1) +
  geom_point(size = 2, color = "blue") +
  annotate("text", x = optimal_pK_bmmc2, y = max(bcmvn_bmmc2$BCmetric),
           label = paste("Optimal pK:", optimal_pK_bmmc2), hjust = -0.1, vjust = -0.5) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white", color = NA)) +
  labs(title = "Optimal pK Selection for BMMC_D1T2", x = "pK", y = "BC Metric")
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/pK_select_bmmc2.png", plot = pk_bmmc2, width = 12, height = 6)

doublet_rate_for_all <- 0.075 #7.5% assumed
nExp_bmmc2 <- round(ncol(BMMC_D1T2) * doublet_rate_for_all)
cat("Expected doublets number:", nExp_bmmc2)

BMMC_D1T2 <- doubletFinder(
  BMMC_D1T2,
  PCs = 1:20,
  pN = 0.25,
  pK = optimal_pK_bmmc2,
  nExp = nExp_bmmc2
)

doublet_dim_bmmc2 <- DimPlot(BMMC_D1T2, group.by = "DF.classifications_0.25_0.04_475") +
  labs(title = "Doublet Classification") +
  theme_minimal()
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/doublet_dim_bmmc2.png", plot = doublet_dim_bmmc2, width = 12, height = 6)

sweep.res.list_cd341 <- paramSweep(CD34_D2T1, PCs = 1:20, sct = FALSE)
sweep.stats_cd341 <- summarizeSweep(sweep.res.list_cd341, GT = FALSE)
bcmvn_cd341 <- find.pK(sweep.stats_cd341)
optimal_pK_cd341 <- as.numeric(as.character(bcmvn_cd341$pK[which.max(bcmvn_cd341$BCmetric)]))
pk_cd341 <- ggplot(bcmvn_cd341, aes(x = pK, y = BCmetric)) +
  geom_line(group = 1, color = "lightblue", size = 1) +
  geom_point(size = 2, color = "blue") +
  annotate("text", x = optimal_pK_cd341, y = max(bcmvn_cd341$BCmetric),
           label = paste("Optimal pK:", optimal_pK_cd341), hjust = -0.1, vjust = -0.5) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white", color = NA)) +
  labs(title = "Optimal pK Selection for CD34_D2T1", x = "pK", y = "BC Metric")
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/pK_select_cd341.png", plot = pk_cd341, width = 12, height = 6)

doublet_rate_for_all <- 0.075 #7.5% assumed
nExp_cd341 <- round(ncol(CD34_D2T1) * doublet_rate_for_all)
cat("Expected doublets number:", nExp_cd341)

CD34_D2T1 <- doubletFinder(
  CD34_D2T1,
  PCs = 1:20,
  pN = 0.25,
  pK = optimal_pK_cd341,
  nExp = nExp_cd341
)

doublet_dim_cd341 <- DimPlot(CD34_D2T1, group.by = "DF.classifications_0.25_0.005_431") +
  labs(title = "Doublet Classification") +
  theme_minimal()
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/doublet_dim_cd341.png", plot = doublet_dim_cd341, width = 12, height = 6)

sweep.res.list_cd342 <- paramSweep(CD34_D3T1, PCs = 1:20, sct = FALSE)
sweep.stats_cd342 <- summarizeSweep(sweep.res.list_cd342, GT = FALSE)
bcmvn_cd342 <- find.pK(sweep.stats_cd342)
optimal_pK_cd342 <- as.numeric(as.character(bcmvn_cd342$pK[which.max(bcmvn_cd342$BCmetric)]))
pk_cd342 <- ggplot(bcmvn_cd342, aes(x = pK, y = BCmetric)) +
  geom_line(group = 1, color = "lightblue", size = 1) +
  geom_point(size = 2, color = "blue") +
  annotate("text", x = optimal_pK_cd342, y = max(bcmvn_cd342$BCmetric),
           label = paste("Optimal pK:", optimal_pK_cd342), hjust = -0.1, vjust = -0.5) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white", color = NA)) +
  labs(title = "Optimal pK Selection for CD34_D3T1", x = "pK", y = "BC Metric")
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/pK_select_cd342.png", plot = pk_cd342, width = 12, height = 6)

doublet_rate_for_all <- 0.075 #7.5% assumed
nExp_cd342 <- round(ncol(CD34_D3T1) * doublet_rate_for_all)
cat("Expected doublets number:", nExp_cd342)

CD34_D3T1 <- doubletFinder(
  CD34_D3T1,
  PCs = 1:20,
  pN = 0.25,
  pK = optimal_pK_cd342,
  nExp = nExp_cd342
)

doublet_dim_cd342 <- DimPlot(CD34_D3T1, group.by = "DF.classifications_0.25_0.005_431") +
  labs(title = "Doublet Classification") +
  theme_minimal()
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/doublet_dim_cd342.png", plot = doublet_dim_cd342, width = 12, height = 6)

#INTERESTING OBSERVATION - THE PCA AND DOUBLETS FOR CD34_D2T1 AND CD34_D3T1 ARE THE SAME OR SIMILAR

#Refiltering data to remove doublets

classification_col_bmmc1 <- "DF.classifications_0.25_0.27_470"
table(BMMC_D1T1@meta.data[[classification_col_bmmc1]])
BMMC_D1T1 <- subset(BMMC_D1T1, subset = DF.classifications_0.25_0.27_470 == "Singlet")

classification_col_bmmc2 <- "DF.classifications_0.25_0.04_475"
table(BMMC_D1T2@meta.data[[classification_col_bmmc2]])
BMMC_D1T2 <- subset(BMMC_D1T2, subset = DF.classifications_0.25_0.04_475 == "Singlet")

classification_col_cd341 <- "DF.classifications_0.25_0.005_431"
table(CD34_D2T1@meta.data[[classification_col_cd341]])
CD34_D2T1 <- subset(CD34_D2T1, subset = DF.classifications_0.25_0.005_431 == "Singlet")

classification_col_cd342 <- "DF.classifications_0.25_0.005_431"
table(CD34_D3T1@meta.data[[classification_col_cd342]])
CD34_D3T1 <- subset(CD34_D3T1, subset = DF.classifications_0.25_0.005_431 == "Singlet")

#Merging the data and Batch Correction
#Merge without batch correction
merged_data_without <- merge(BMMC_D1T1, y = c(BMMC_D1T2, CD34_D2T1, CD34_D3T1), add.cell.ids = c("BMMC_D1T1", "BMMC_D1T2", "CD34_D2T1","CD34_D3T1"), project = "Project1")

merged_data_without <- NormalizeData(merged_data_without) %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
  #RunUMAP(dims = 1:50)

merged_without_dim <- DimPlot(merged_data_without, group.by = "orig.ident") + ggtitle("Merged Without Batch Correction")
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/merged_without_dim.png", plot = merged_without_dim, width = 12, height = 6)

#Merge with seurat

dataset_list <- list(BMMC_D1T1, BMMC_D1T2, CD34_D2T1, CD34_D3T1)
dataset_list <- lapply(dataset_list, NormalizeData)
anchors <- FindIntegrationAnchors(object.list = dataset_list, dims = 1:50)
merged_data_batch <- IntegrateData(anchorset = anchors)

#Preprocess and PCA again
merged_data_batch <- ScaleData(merged_data_batch) %>%
  FindVariableFeatures() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:50)
merged_with_dim <- DimPlot(merged_data_batch, group.by = "orig.ident") + ggtitle("Merged With Batch Correction")
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/merged_with_dim.png", plot = merged_with_dim, width = 12, height = 6)
merged_with_umap <- DimPlot(merged_data_batch, reduction = "umap")
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/merged_with_umap.png", plot = merged_with_umap, width = 12, height = 6)
merged_elbow <- ElbowPlot(merged_data_batch, ndims = 50) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = "black"),
        plot.background = element_rect(fill = "white", color = NA)) +
  geom_vline(xintercept = 20, linetype = "dashed", color = "black", size = 0.5)
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/merged_elbow.png", plot = merged_elbow, width = 6, height = 4)

#The elbow plot stable after 20 

#Side-by-side comparison of merging without batch correction and with batch correction
combined_dim_plot <- merged_without_dim | merged_with_dim
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/comparison_batch_plot.png", plot = combined_dim_plot, width = 12, height = 7)

#Dimensionality reduction and clustering

## repeating the pca and umap here (CAN REMOVE THIS AND DIRECTLY START WITH CLUSTERING AND FIND NEIGHBOURS)
merged_data_batch <- RunPCA(merged_data_batch)
merged_pca_dim <- DimPlot(merged_data_batch, reduction = "pca") + ggtitle("Merged Data Final PCA")
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/merged_pca_dim.png", plot = merged_pca_dim, width = 12, height = 6)
merged_umap_dim <- DimPlot(merged_data_batch, reduction = "umap") + ggtitle("Merged Data Final UMAP")
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/merged_umap_dim.png", plot = merged_umap_dim, width = 12, height = 6)

#Clustering
merged_data_batch <- FindNeighbors(merged_data_batch, dims = 1:20)
merged_data_batch <- FindClusters(merged_data_batch, resolution = 0.2)
merged_data_batch <- RunUMAP(merged_data_batch, dims = 1:20)
merged_cluster_umap <- DimPlot(merged_data_batch, reduction = "umap", label = TRUE) + ggtitle("Merged data clustering")
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/merged_cluster_umap.png", plot = merged_cluster_umap, width = 12, height = 6)

## Week 3
#Automatic annotation

reference_hpcad <- celldex::HumanPrimaryCellAtlasData()
exp_matrix <- GetAssayData(merged_data_batch, assay = "integrated")
singleR_result <- SingleR(test = exp_matrix,
                          ref = reference_hpcad,
                          labels = reference_hpcad$label.main)
merged_data_batch[["SingleR.labels"]] <- singleR_result$labels
labelled_umap <- DimPlot(merged_data_batch, group.by = "SingleR.labels", label = TRUE, repel = TRUE) + 
  ggtitle("UMAP with labels") +
  theme_minimal()
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/labelled_umap.png", plot = labelled_umap, width = 12, height = 6)

#Manual Annotation
c_markers <- FindAllMarkers(
  merged_data_batch,
  only.pos = TRUE,        
  min.pct = 0.25,           
  logfc.threshold = 0.25    # Minimum log fold change threshold
)


write.csv(c_markers, "cluster_markers.csv")

markers_type <- list(
    HSC = c("CD34", "CD38", "Sca1", "Kit"),
    LMPP = c("CD38", "CD52", "CSF3R", "ca1", "Kit", "CD34", "Flk2"),
    CLP = c("IL7R"),
    GMP_Neutrophils = c("ELANE"),
    CMP = c("IL3", "GM-CSF", "M-CSF"),
    B_Cells = c("CD19"),
    Pre_B = c("CD19", "CD34"),
    Plasma = c("SDC1", "IGHA1", "IGLC1", "MZB1", "JCHAIN"),
    T_CD8 = c("CD3D", "CD3E", "CD8A", "CD8B"),
    T_CD4 = c("CD3D", "CD3E", "CD4"),
    NKC = c("FCGR3A", "NCAM1", "NKG7", "KLRB1"),
    Eryth = c("GATA1", "HBB", "HBA1", "HBA2"),
    pDC = c("IRF8", "IRF4", "IRF7"),
    cDC = c("CD1C", "CD207", "ITGAM", "NOTCH2", "SIRPA"),
    CD14_M = c("CD14", "CCL3", "CCL4", "IL1B"),
    CD16_M = c("FCGR3A", "CD68", "S100A12"),
    Baso = c("GATA2"))

match_cluster_to_cell_type <- function(cluster_markers, markers_type) {
  cell_type_scores <- sapply(names(markers_type), function(cell_type) {
    markers <- markers_type[[cell_type]]
    sum(cluster_markers %in% markers) 
  })
  best_match <- names(cell_type_scores)[which.max(cell_type_scores)]
  return(best_match)
}
clusters <- levels(merged_data_batch)

#manual plotting
cluster_annotations <- vector("character", length = length(clusters))
names(cluster_annotations) <- clusters

for (cluster in clusters) {
  cluster_markers <- c_markers %>%
    filter(cluster == !!cluster) %>%   
    #arrange(p_val_adj) %>%              
    pull(gene) %>%                     
    unique()
  cell_type <- match_cluster_to_cell_type(cluster_markers, markers_type)
  cluster_annotations[cluster] <- paste0(cell_type, " (", cluster, ")")
}

cluster_annotations

str(Idents(merged_data_batch))
merged_data_batch$ManualAnnotations <- factor(Idents(merged_data_batch), 
                                              levels = names(cluster_annotations), 
                                              labels = cluster_annotations)
table(merged_data_batch$ManualAnnotations)

manual_map <- DimPlot(merged_data_batch, group.by = "ManualAnnotations", label = TRUE) +
  ggtitle("UMAP with manual Cell type annotations")
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/manual_map.png", plot = manual_map, width = 12, height = 6)

meta_clusters <- factor(merged_data_batch$ManualAnnotations) 
meta_cluster_ids <- as.integer(meta_clusters) 

merged_data_batch$MetaClusters <- meta_clusters

automated_ann <- DimPlot(merged_data_batch, reduction = "umap", group.by = "SingleR.labels", label = FALSE, repel = TRUE) + ggtitle("Automatic annotations") + theme_minimal()
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/automated_compare.png", plot = automated_ann, width = 12, height = 6)

manual_ann <- DimPlot(merged_data_batch, reduction = "umap", group.by = "ManualAnnotations", label = FALSE, repel = TRUE) + ggtitle("Manual annotations") + theme_minimal()
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/manual_compare.png", plot = manual_ann, width = 12, height = 6)

automated_ann + manual_ann #comparison of automated and manual annotations

#We selected CD38 because HSC and LMPP are most occurring clusters and both have CD38 as a common marker. CD19 is a marker for B cells and Progenitor B cells. CD14 is from CD14+ monocytes  and also  occur twice in manual annotation clusters.

marker_plot <- VlnPlot(merged_data_batch,
                       features = c("CD38", "CD19", "CD14"),
                       group.by = "ManualAnnotations",
                       pt.size = 0) + 
  theme_minimal()
ggsave("/Users/kathanpandya/Desktop/WiSe24/singlecell/Assignment1/fotos/marker_plot.png", plot = marker_plot, width = 12, height = 6)

#Differential Expression Analysis
b_vs_t <- FindMarkers(merged_data_batch, 
                      ident.1 = "B cells", 
                      ident.2 = "T cells", 
                      assay = "RNA",
                      group.by = "ManualAnnotations",
                      logfc.threshold = 0.25, 
                      min.pct = 0.1)

t_vs_mono <- FindMarkers(merged_data_batch, 
                         ident.1 = "T cells", 
                         ident.2 = "Monocytes", 
                         assay = "RNA", 
                         logfc.threshold = 0.25, 
                         min.pct = 0.1)

