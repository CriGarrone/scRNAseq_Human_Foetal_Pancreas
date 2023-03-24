### written by Maria Cristina Garrone (Dr Rocio Sancho Lab)
# Figure for Manea et al 2023

########################################################################################################################
########################################################################################################################
# Figure 5 - USP7 expression during human embryonic pancreas development. 
########################################################################################################################
########################################################################################################################

# load libraries and data 

library(Seurat)
library(ggplot2)
library(dplyr)
library(scCustomize)

all.weeks.epi <- readRDS("~/Documents/Projects/USP7/HumanFoetalDatasets/Yu_et_al_2021/Objects/all.weeks.epi.RDS")

########################################################################################################################
# A
########################################################################################################################

pdf("~/Documents/Projects/USP7/HumanFoetalDatasets/Yu_et_al_2021/Results/5A.pdf")
DimPlot(all.weeks.epi)
dev.off()

########################################################################################################################
# B
########################################################################################################################

pdf("~/Documents/Projects/USP7/HumanFoetalDatasets/Yu_et_al_2021/Results/5B.pdf")
DimPlot(all.weeks.epi, split.by = "sample")
dev.off()

########################################################################################################################
# C
########################################################################################################################

pdf("~/Documents/Projects/USP7/HumanFoetalDatasets/Yu_et_al_2021/Results/5C.pdf")
FeaturePlot(all.weeks.epi, features = c("USP7","NEUROG3"), blend=T, blend.threshold = 0, cols = c("gray", "blue", "red"), order = T, combine = T)
dev.off()

########################################################################################################################
# D
########################################################################################################################

DefaultAssay(endocrine) <- "RNA"
endocrine <- NormalizeData(endocrine)
all.genes <- rownames(endocrine)
endocrine <- ScaleData(endocrine, features = all.genes)

VlnPlot(endocrine2, features = c("USP7","PDX1", "NEUROG3"), slot = "data")
DotPlot(endocrine2, features = c("PDX1", "USP7","NEUROG3"),  cols = c("gray", "red"),col.min = 0,dot.scale = 8) + RotatedAxis() + coord_flip() 
DotPlot(endocrine2, features = c("PDX1", "USP7","NEUROG3" ), group.by = "sample", cols = c("gray", "red"),col.min = 0, dot.scale = 8) + RotatedAxis() + coord_flip() 

########################################################################################################################
# E
########################################################################################################################

# Subset endocrine cells only
endocrine <- subset(all.weeks.epi, idents = c("EP","Alpha", "Beta", "Delta","PP") )

DefaultAssay(endocrine) <- "RNA"
endocrine <- NormalizeData(endocrine)
all.genes <- rownames(endocrine)
endocrine <- ScaleData(endocrine, features = all.genes)

DefaultAssay(endocrine) <- "SCT"
endocrine <- RunPCA(endocrine, features = VariableFeatures(endocrine), verbose = T)
endocrine <- RunUMAP(endocrine, reduction = "pca", dims = 1:30, verbose = F)

pdf("~/Documents/Projects/USP7/HumanFoetalDatasets/Yu_et_al_2021/Results/5E.pdf")
DimPlot(endocrine)
Dev.off()

########################################################################################################################
# F
########################################################################################################################

pdf("~/Documents/Projects/USP7/HumanFoetalDatasets/Yu_et_al_2021/Results/5F.pdf")
FeaturePlot(endocrine, features = c("USP7","NEUROG3"), blend=T, blend.threshold = 0, cols = c("gray", "blue", "red"), order = T, combine = T) 
dev.off()

########################################################################################################################
# G
########################################################################################################################

DefaultAssay(endocrine) <- "RNA"

pdf("~/Documents/Projects/USP7/HumanFoetalDatasets/Yu_et_al_2021/Results/5G.pdf")
DotPlot(endocrine, features = c("PDX1", "USP7","NEUROG3"),  cols = c("gray", "red"),col.min = 0,dot.scale = 8) + RotatedAxis() + coord_flip() 
dev.off()

########################################################################################################################
# H
########################################################################################################################

# calculate percentage of cells expressing NGN3 and USP7

endocrine.ngn3 <- WhichCells(endocrine, expression = NEUROG3 >0)
endocrine.ngn3.neg <- WhichCells(endocrine, expression = NEUROG3 ==0)
pos_cells <- subset(endocrine,cells=endocrine.ngn3)
neg_cells <- subset(endocrine, cells = endocrine.ngn3.neg)

percent_express <- Percent_Expressing(seurat_object = pos_cells, features = c("USP7"),group_by = "sample")
percent_express_entire  <- Percent_Expressing(seurat_object = pos_cells, features = c("USP7"),entire_object = T)
percent_express_entire_neg  <- Percent_Expressing(seurat_object = neg_cells, features = c("USP7"),entire_object = T)
percent_express_neg  <- Percent_Expressing(seurat_object = neg_cells, features = c("USP7"), group_by = "sample")

percent_express <- as.data.frame(t(percent_express))
percent_express$sample <- rownames(percent_express)
percent_express$Percentage_USP7Pos <- percent_express$USP7

percent_express_neg <- as.data.frame(t(percent_express_neg))
percent_express_neg$sample <- rownames(percent_express_neg)
percent_express_neg$Percentage_USP7Pos <- percent_express_neg$USP7

# percentage of ngn3- cells that are also USP7+
levels <- c("W8","W10","W12", "W14", "W16", "W18", "W19")
percent_express_neg$sample <- factor(x= percent_express_neg$sample, levels = levels)

pdf("~/Documents/Projects/USP7/HumanFoetalDatasets/Yu_et_al_2021/Results/5H.pdf")
ggplot(percent_express_neg, aes(x = sample, y = Percentage_USP7Pos, fill = sample))+
  geom_bar(stat = "identity")
dev.off()

########################################################################################################################
# I
########################################################################################################################

# percentage of ngn3+ cells that are also USP7+

pdf("~/Documents/Projects/USP7/HumanFoetalDatasets/Yu_et_al_2021/Results/5I.pdf")
ggplot(percent_express, aes(x = sample, y = Percentage_USP7Pos, fill = sample))+
  geom_bar(stat = "identity")
dev.off()

########################################################################################################################
########################################################################################################################
# Figure S6 - Quality control and clustering of scRNA-Seq dataset. 
########################################################################################################################
########################################################################################################################

# load libraries and data 

library(Seurat)
library(ggplot2)
library(dplyr)

all.weeks <- readRDS("~/Documents/Projects/USP7/HumanFoetalDatasets/Yu_et_al_2021/Objects/all.weeks.RDS")


########################################################################################################################
# A
########################################################################################################################

pdf("~/Documents/Projects/USP7/HumanFoetalDatasets/Yu_et_al_2021/Results/S1A.pdf")
VlnPlot(all.weeks, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0, group.by = "sample")
dev.off()

########################################################################################################################
# B
########################################################################################################################

pdf("~/Documents/Projects/USP7/HumanFoetalDatasets/Yu_et_al_2021/Results/S1B.pdf")
DimPlot(all.weeks)
dev.off()

########################################################################################################################
# C
########################################################################################################################

pdf("~/Documents/Projects/USP7/HumanFoetalDatasets/Yu_et_al_2021/Results/S1C.pdf")
cell.type <- c( "THY1", "COL3A1", "PDX1", "PROX1", "EPCAM","PECAM1", "ASCL1", "PTPRC")
DotPlot(all.weeks, features = cell.type, cols= c("gray", "red"), dot.scale = 8) + RotatedAxis() 
dev.off()

########################################################################################################################
# D
########################################################################################################################

pdf("~/Documents/Projects/USP7/HumanFoetalDatasets/Yu_et_al_2021/Results/S1D.pdf")
DimPlot(all.weeks.epi, group.by = "seurat_clusters")
dev.off()

########################################################################################################################
# E
########################################################################################################################

pdf("~/Documents/Projects/USP7/HumanFoetalDatasets/Yu_et_al_2021/Results/S1E.pdf")
DefaultAssay(all.weeks.epi) <- "SCT"
cells <- c("HES1", "SAT1","RBPJL","GP2","NEUROG3", "GCG", "INS","SST", "PPY")
DotPlot(all.weeks.epi, features = cells, col.min = 0, cols= c("gray", "red"), dot.scale = 8, group.by = "seurat_clusters") + RotatedAxis()  
dev.off()

########################################################################################################################
# F
########################################################################################################################

pdf("~/Documents/Projects/USP7/HumanFoetalDatasets/Yu_et_al_2021/Results/S1F.pdf")
FeaturePlot(all.weeks.epi, features = "NEUROG3", order=T, cols = c("gray","red"))
dev.off()

########################################################################################################################
# G
########################################################################################################################

pdf("~/Documents/Projects/USP7/HumanFoetalDatasets/Yu_et_al_2021/Results/S1G.pdf")
FeaturePlot(all.weeks.epi, features = "USP7", order=T, cols = c("gray","blue"))
dev.off()
