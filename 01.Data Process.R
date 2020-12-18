#Final script for Toshie's data
##Authro: Xin Zhou
#Created: 2019-07-31

#last updated 2019-02-04

#load required package
library("Seurat")
library("dplyr")
library("Matrix")
library("ggpubr")
library("DoubletFinder")
library("ggsci")
library("DESeq2")
library("ggrepel")
library("scales")
library("cowplot")

#input file
pbmc.data <- Read10X(data.dir = "~/Box/XinZhouFiles/Projects/ZXP1_PAH/scRNA_DC_Toshie (Rabinovitch lab)/hg19/") 

#Start of the pipeline
pbmc <- CreateSeuratObject(pbmc.data, min.cells = 3, min.features=700, project = "10X_PBMC_DC")
pbmc@meta.data$subject <- sub("^([^-]+-)","S",row.names(pbmc@meta.data))
pbmc@meta.data$condition <-  pbmc@meta.data$subject
pbmc@meta.data$condition[pbmc@meta.data$subject %in% c("S1", "S2", "S3", "S4")]<- "CONT"
pbmc@meta.data$condition[pbmc@meta.data$subject %in% c("S5", "S6", "S7")]<- "PAH"

#https://science.sciencemag.org/content/352/6282/189 
hkgene <- read.table("~/Box/XinZhouFiles/Projects/ZXP1_PAH/Single Cell Related/House_keeping_Gene.txt")
hkgenes <- as.vector(hkgene$V1)
hkgenes.found <- which(toupper(rownames(pbmc@assays$RNA)) %in% hkgenes)
n.expressed.hkgenes <- Matrix::colSums(pbmc@assays$RNA[hkgenes.found, ] > 0)
pbmc <- AddMetaData(object = pbmc, metadata = n.expressed.hkgenes, col.name = "n.exp.hkgenes")

colnames(pbmc)

mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@assays$RNA), value = TRUE)

percent.mito <- Matrix::colSums(pbmc@assays$RNA@counts[mito.genes, ]) / Matrix::colSums(pbmc@assays$RNA@counts)
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")

pbmc <- subset(pbmc, subset= nFeature_RNA > 700 & nFeature_RNA <2500 &  percent.mito < 0.05 & nCount_RNA >2000 & n.exp.hkgenes > 60)
pbmc

#create a dataset for reference
pbmc4 <- subset(pbmc, subject=="S1"|subject=="S2"|subject=="S3"|subject=="S4")
pbmc4 <- subset(pbmc4, subset= nFeature_RNA > 1000 & nFeature_RNA <1400 & percent.mito < 0.035 & percent.mito > 0.02 & nCount_RNA >3000)
pbmc4 

pbmc4@meta.data$subject <- "S8"
pbmc4@meta.data$condition <- "control"

#merge file together
DClist <- SplitObject(pbmc, split.by = "subject")
DClist <- c(DClist, pbmc4)
names(DClist) <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8")
names(DClist) 
DClist

#use sct transform before integret
for (i in names(DClist)) {
  DClist[[i]] <- SCTransform(DClist[[i]], verbose = T)
}

DC.features <- SelectIntegrationFeatures(object.list = DClist, nfeatures = 1000)
DClist <- PrepSCTIntegration(object.list = DClist, anchor.features = DC.features)

reference_dataset <- which(names(DClist) == c("S8"))

DC.anchors <- FindIntegrationAnchors(object.list = DClist, normalization.method = "SCT", dims = 1:30,
                                     anchor.features = DC.features, reference = reference_dataset)

DC.integrated <- IntegrateData(anchorset = DC.anchors, normalization.method = "SCT")

#remove the reference data
DC.integrated <- subset(DC.integrated, subject != "S8")
DC.integrated <- FindVariableFeatures(DC.integrated, selection.method = "vst", nfeatures = 2000, assay = "SCT")

all.genes <- rownames(DC.integrated)
DC.integrated <- ScaleData(DC.integrated, features = all.genes)
DC.integrated <- RunPCA(DC.integrated, features = VariableFeatures(object = DC.integrated))

DC.integrated <- FindNeighbors(DC.integrated, dims = 1:20)
DC.integrated <- FindClusters(DC.integrated, resolution = 0.7)
head(Idents(DC.integrated), 5)
DimPlot(DC.integrated)

DC.integrated <- RunUMAP(DC.integrated, dims = 1:15)
DimPlot(DC.integrated, reduction = "umap", label = T)
DimPlot(DC.integrated, reduction = "umap", label = T, split.by = "subject", ncol=4)
VlnPlot(DC.integrated, features = c("IFITM3", "CLEC9A", "CD1C", "THBD","CD14", "GAPDH", "ACTB"), assay = "integrated", pt.size=0.001, ncol=3, group.by = "subject")
FeaturePlot(DC.integrated, features = "SIRPA")

saveRDS(DC.integrated, file="~/Box/XinZhouFiles/Projects/ZXP1_PAH/scRNA_DC_Toshie (Rabinovitch lab)/Figure result/Final analysis/01.Data Process/DCsuratObect.Final.RData")





