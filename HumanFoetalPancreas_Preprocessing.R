### written by Maria Cristina Garrone (Dr Rocio Sancho Lab)
# Analysis of the Human Foetal Pancreas Dataset from Yu et al 2021

# Load libraries
library(DESeq2)
library(Seurat)
library(clustree)
library(dplyr)
library(tidyverse)
library(ggplot2)

#Human_10x.UMI.txt.gz

# Load data 
Xu_paper <- read.table("~/Rotation3/scRNAseq_other_datasets/Xu_2021/RawData/Human_10x.UMI.txt", sep = "\t", header = TRUE)
Xu_paper <- as.data.frame(Xu_paper)

rownames(Xu_paper) <- Xu_paper$Symbol
Xu_paper2 <- Xu_paper[,-1]
head(Xu_paper2[ , 1:3])

# Save original matrix
save(Xu_paper2, file = "~/Rotation3/scRNAseq_other_datasets/Xu_2021/PreProcessing/Xu_UMI_matrix.Rda")

# to load
load(file = "~/Rotation3/scRNAseq_other_datasets/Xu_2021/PreProcessing/obj/Xu_UMI_matrix.Rda")

# subset each sample individually
Xu_paper2_W8 <- Xu_paper2[ , grepl( "^W8" , names( Xu_paper2 ) ) ]
Xu_paper2_W10 <- Xu_paper2[ , grepl( "^W10" , names( Xu_paper2 ) ) ]
Xu_paper2_W12 <- Xu_paper2[ , grepl( "^W12_" , names( Xu_paper2 ) ) ]
Xu_paper2_W14 <- Xu_paper2[ , grepl( "^W14" , names( Xu_paper2 ) ) ]
Xu_paper2_W16 <- Xu_paper2[ , grepl( "^W16" , names( Xu_paper2 ) ) ]
Xu_paper2_W18 <- Xu_paper2[ , grepl( "^W18" , names( Xu_paper2 ) ) ]
Xu_paper2_W19 <- Xu_paper2[ , grepl( "^W19" , names( Xu_paper2 ) ) ]

# Create Seurat Object for each sample, select a minimum of 3 cells and 100 features
Xu_W8_seurat <- CreateSeuratObject(counts = Xu_paper2_W8, 
                                   min.cells = 3,
                                   min.features = 100, 
                                   project = "Xu_2021")

Xu_W10_seurat <- CreateSeuratObject(counts = Xu_paper2_W10, 
                                    min.cells = 3,
                                    min.features = 100, 
                                    project = "Xu_2021")

Xu_W12_seurat <- CreateSeuratObject(counts = Xu_paper2_W12, 
                                    min.cells = 3,
                                    min.features = 100, 
                                    project = "Xu_2021")

Xu_W14_seurat <- CreateSeuratObject(counts = Xu_paper2_W14, 
                                   min.cells = 3,
                                   min.features = 100, 
                                   project = "Xu_2021")

Xu_W16_seurat <- CreateSeuratObject(counts = Xu_paper2_W16, 
                                   min.cells = 3,
                                   min.features = 100, 
                                   project = "Xu_2021")

Xu_W18_seurat <- CreateSeuratObject(counts = Xu_paper2_W18, 
                                   min.cells = 3,
                                   min.features = 100, 
                                   project = "Xu_2021")

Xu_W19_seurat <- CreateSeuratObject(counts = Xu_paper2_W19, 
                                   min.cells = 3,
                                   min.features = 100, 
                                   project = "Xu_2021")

# add sample timepoint to the metadata
Xu_W8_seurat$sample <- "W8"
Xu_W10_seurat$sample <- "W10"
Xu_W12_seurat$sample <- "W12"
Xu_W14_seurat$sample <- "W14"
Xu_W16_seurat$sample <- "W16"
Xu_W18_seurat$sample <- "W18"
Xu_W19_seurat$sample <- "W19"

# Merge all samples into one Seurat Object 
all.weeks <- merge(Xu_W8_seurat, c(Xu_W10_seurat, Xu_W12_seurat,Xu_W14_seurat, Xu_W16_seurat, Xu_W18_seurat, Xu_W19_seurat),
                  add.cell.ids = c("W8","W10", "W12", "W14", "W16", "W18", "W19"),
                  project = "merged_Xu")

# Reorder the timepoint
all.weeks$sample <- factor(x = all.weeks$sample, levels = c("W8", "W10", "W12", "W14", "W16", "W18", "W19"))
levels(x = all.weeks$sample)

# Perform quality Control
all.weeks <- PercentageFeatureSet(all.weeks, pattern = "^MT-", col.name = "percent.mt")
VlnPlot(all.weeks, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, group.by = "sample")
all.weeks <- subset(all.weeks, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & nCount_RNA > 3000 & nCount_RNA < 20000 & percent.mt < 15)

# Normalize and scale the data
all.weeks <- SCTransform(all.weeks, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = T)

# Dimensionality reduction and clustering
all.weeks <- RunPCA(all.weeks, features = VariableFeatures(all.weeks), verbose = T)
all.weeks <- RunUMAP(all.weeks, reduction = "pca", dims = 1:30, verbose = F)
all.weeks <- FindNeighbors(all.weeks, reduction = "pca", dims = 1:30)
all.weeks <- FindClusters(all.weeks, resolution = 0.2)
table(all.weeks$orig.ident)

# Relabel clusters
DimPlot(all.weeks)
cell.type <- c( "THY1", "COL3A1", "PDX1", "PROX1", "EPCAM","PECAM1", "ASCL1", "PTPRC")
DotPlot(all.weeks, features = cell.type, cols= c("gray", "red"), dot.scale = 8) + RotatedAxis() 

all.weeks <- RenameIdents(object = all.weeks, 
                         "0" = "Mesenchyme",
                         "1" = "Epithelium",
                         "2" = "Mesenchyme",
                         "3" = "Mesenchyme",
                         "4" = "Endothelium",
                         "5" = "Epithelium",
                         "6" = "Epithelium",
                         "7" = "Epithelium",
                         "8" = "Endothelium",
                         "9" = "Neurons",
                         "10" = "Immune",
                         "11" = "Immune",
                         "12" = "Immune",
                         "13" = "Endothelium",
                         "14" = "Mesenchyme",
                         "15" = "Mesenchyme", 
                         "16" = "Mesenchyme",
                         "17" = 'Endothelium')

DimPlot(all.weeks, group.by = "seurat_clusters")

saveRDS(all.weeks, "~/Rotation3/scRNAseq_other_datasets/Xu_2021/all.weeks//Xu_all.rds")
all.weeks <- readRDS("~/Rotation3/scRNAseq_other_datasets/Xu_2021/all.weeks//Xu_all.rds")

# Subset epithelial cells only
all.weeks.epi <- subset(x = all.weeks, idents = "Epithelium")
all.weeks.epi <- RunPCA(all.weeks.epi, features = VariableFeatures(all.weeks.epi), verbose = T)
all.weeks.epi <- RunUMAP(all.weeks.epi, reduction = "pca", dims = 1:30, verbose = F)
all.weeks.epi <- FindNeighbors(all.weeks.epi, reduction = "pca", dims = 1:30)

# Select resolution with clustree 
clustree(all.weeks.epi)
all.weeks.epi <- FindClusters(all.weeks.epi, resolution = seq(0,1,0.1))
DimPlot(all.weeks.epi)

all.weeks.epi <- FindClusters(all.weeks.epi, resolution = 0.5)
table(all.weeks.epi$sample)
DimPlot(all.weeks.epi, group.by = "sample")
DimPlot(all.weeks.epi)

# Perform DGE 
all <- FindAllMarkers(all.weeks.epi, only.pos = T)
all$pct.diff <- (all$pct.1 - all$pct.2)
all <- dplyr::arrange(all, cluster, desc(avg_logFC), desc(pct.diff))
all$gene <- rownames(all)
write_csv(all, "~/Rotation3/FinalRNAseq/Xu_all_tp/all.clusters.csv")

# Relabel clusters
cells <- c("HES1", "SAT1","RBPJL","GP2","NEUROG3", "GCG", "INS","SST", "PPY")
DotPlot(all.week.epi, features = cells, col.min = 0, cols= c("gray", "red"), dot.scale = 8, group.by = "seurat_clusters") + RotatedAxis()  
DimPlot(all.weeks.epi) + DimPlot(all.weeks.epi, group.by = "sample") 

all.weeks.epi <- RenameIdents(object = all.weeks.epi, 
                          "0" = "Acinar",
                          "1" = "Beta",
                          "2" = "Acinar",
                          "3" = "Trunk/Progenitor",
                          "4" = "Tip",
                          "5" = "Acinar",
                          "6" = "Tip",
                          "7" = "Alpha",
                          "8" = "Beta",
                          "9" = "Ductal",
                          "10" = "Delta",
                          "11" = "Beta",
                          "12" = "PP",
                          "13" = "Beta",
                          "14" = "Beta",
                          "15" = "EP", 
                          "16" = "Trunk/Progenitor")

all.weeks.epi@active.ident <- factor(all.weeks.epi.2@active.ident,
                                levels= c("Trunk/Progenitor", "Ductal","Tip", "Acinar","EP", "Alpha", "Beta","Delta", "PP"))
DimPlot(all.weeks.epi)

saveRDS(all.weeks.epi, "~/Rotation3/FinalRNAseq/Xu_all_tp/all.weeks.epi.RDS")
all.weeks.epi <- readRDS("~/Documents/Rotation3/Results/RNAseq/Xu/Objects/all.weeks.epi.RDS")

# Subset endocrine cells only
endocrine <- subset(all.week.epi, idents = c("EP","Alpha", "Beta", "Delta","PP") )
DimPlot(endocrine)
endocrine <- RunPCA(endocrine, features = VariableFeatures(endocrine), verbose = T)
endocrine <- RunUMAP(endocrine, reduction = "pca", dims = 1:30, verbose = F)
endocrine <- FindNeighbors(endocrine, reduction = "pca", dims = 1:30)
DimPlot(endocrine)
DimPlot(endocrine, group.by = "sample")


