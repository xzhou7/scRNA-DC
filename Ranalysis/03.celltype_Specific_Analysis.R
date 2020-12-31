#Final script for Toshie's data
##Authro: Xin Zhou
#Created: 2019-07-31

#last updated 2019-02-10

#load necessary package
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
library("ggsci")
library("monocle")
library("MAST")

setwd("~/Box/XinZhouFiles/Projects/ZXP1_PAH/scRNA_DC_Toshie (Rabinovitch lab)/Figure result/Final analysis/")
DC.integrated <- readRDS(file="~/Box/XinZhouFiles/Projects/ZXP1_PAH/scRNA_DC_Toshie (Rabinovitch lab)/Figure result/Final analysis/02.UMAP and Cell Identity/DCsuratObect.Final02.RData")
DC.integrated

Idents(DC.integrated) <- "celltype"
Idents(DC.integrated)
#Levels: FCER1A_DC CD1C_A S100A9_DC CEAC9_DC_R ANXA1_DC CEAC9_DC_A IFIT_DC LYZ_DC Intermediate_Cell
DC.LYZ     <- subset(DC.integrated, idents="LYZ_DC")
DC.FCER1A  <- subset(DC.integrated, idents="FCER1A_DC")
DC.CD1C    <- subset(DC.integrated, idents="CD1C_A")
DC.S100A9  <- subset(DC.integrated, idents="S100A9_DC")
DC.ANXA1   <- subset(DC.integrated, idents="ANXA1_DC")
DC.IFIT    <- subset(DC.integrated, idents="IFIT_DC")
DC.CEAC9_R <- subset(DC.integrated, idents="CEAC9_DC_R")
DC.CEAC9_A <- subset(DC.integrated, idents="CEAC9_DC_A")
DC1        <- subset(DC.integrated, idents=c("LYZ_DC", "FCER1A_DC","CD1C_A","ANXA1_DC","IFIT_DC"))
DC.CEAC9   <- subset(DC.integrated, idents=c("CEAC9_DC_R","CEAC9_DC_A"))

#analysis by cluster
DC.LYZ <- subset(DC.integrated, idents="LYZ_DC")
DefaultAssay(DC.LYZ) <- "RNA"
Idents(DC.LYZ) <- "condition"
DefaultAssay(DC.LYZ)
Idents(DC.LYZ)
table(DC.LYZ$condition)

LYZ.markers.RNA <- FindMarkers(DC.LYZ, ident.1 = "CONT", ident.2 = "PAH", min.pct = 0.01,
                                                test.use="DESeq2", assay="RNA",logfc.threshold=0.3, max.cells.per.ident =180)
LYZ.markers.RNA
LYZ.markers.RNA$genes <- row.names(LYZ.markers.RNA)
cluster3.markers <- LYZ.markers.RNA

p <- ggplot(cluster3.markers, aes(x=avg_logFC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("Deseq Result on LYZ DC Control/PAH DEGs")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC < -0.17 & p_val_adj < 0.001), color="#8C1515")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC > 0.17 & p_val_adj < 0.001), color="#09425A")
p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC < -0.3  & p_val_adj < 0.001), aes(label=genes), color="#8C1515")
p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC >0.3 & p_val_adj < 0.001), aes(label=genes), color="#09425A")
p <- p + xlim(-1,1)
p

write.csv(file = "./03.Celltypeanalysis/LYZ_DEG.csv", LYZ.markers.RNA)
ggsave2(filename = "./03.Celltypeanalysis/LYZ_DEG.pdf", p)

#FCER1A
DefaultAssay(DC.FCER1A) <- "RNA"
Idents(DC.FCER1A) <- "condition"
DefaultAssay(DC.FCER1A)
Idents(DC.FCER1A)
table(DC.FCER1A$condition)

FCER1A.markers.RNA <- FindMarkers(DC.FCER1A, ident.1 = "CONT", ident.2 = "PAH",
                               min.pct = 0.01,test.use="DESeq2", assay="RNA",logfc.threshold=0.3, max.cells.per.ident =1300)

FCER1A.markers.RNA
FCER1A.markers.RNA$genes <- row.names(FCER1A.markers.RNA)
cluster3.markers <- FCER1A.markers.RNA

p <- ggplot(cluster3.markers, aes(x=avg_logFC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("Deseq Result on FCER1A DC Control/PAH DEGs")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC < -0.17 & p_val_adj < 0.001), color="#8C1515")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC > 0.17 & p_val_adj < 0.001), color="#09425A")
p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC < -0.3  & p_val_adj < 0.001), aes(label=genes), color="#8C1515")
p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC >0.3 & p_val_adj < 0.001), aes(label=genes), color="#09425A")
p <- p + xlim(-1,1)
p

write.csv(file = "./03.Celltypeanalysis/FCER1A_DEG.csv", FCER1A.markers.RNA)
ggsave2(filename = "./03.Celltypeanalysis/FCER1A_DEG.pdf", p)

#CD1C
DefaultAssay(DC.CD1C) <- "RNA"
Idents(DC.CD1C) <- "condition"
DefaultAssay(DC.CD1C)
Idents(DC.CD1C)
table(DC.CD1C$condition)

CD1C.markers.RNA <- FindMarkers(DC.CD1C, ident.1 = "CONT", ident.2 = "PAH",
                                  min.pct = 0.01,test.use="DESeq2", assay="RNA",logfc.threshold=0.3, max.cells.per.ident =900)

CD1C.markers.RNA
CD1C.markers.RNA$genes <- row.names(CD1C.markers.RNA)
cluster3.markers <- CD1C.markers.RNA

p <- ggplot(cluster3.markers, aes(x=avg_logFC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("Deseq Result on CD1C DC Control/PAH DEGs")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC < -0.17 & p_val_adj < 0.001), color="#8C1515")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC > 0.17 & p_val_adj < 0.001), color="#09425A")
p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC < -0.3  & p_val_adj < 0.001), aes(label=genes), color="#8C1515")
p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC >0.3 & p_val_adj < 0.001), aes(label=genes), color="#09425A")
p

write.csv(file = "./03.Celltypeanalysis/CD1C_DEG.csv", CD1C.markers.RNA)
ggsave2(filename = "./03.Celltypeanalysis/CD1C_DEG.pdf", p)


#S100A9
DefaultAssay(DC.S100A9) <- "RNA"
Idents(DC.S100A9) <- "condition"
DefaultAssay(DC.S100A9)
Idents(DC.S100A9)
table(DC.S100A9$condition)

S100A9.markers.RNA <- FindMarkers(DC.S100A9, ident.1 = "CONT", ident.2 = "PAH",
                                min.pct = 0.01,test.use="DESeq2", assay="RNA",logfc.threshold=0.3, max.cells.per.ident =600)

S100A9.markers.RNA
S100A9.markers.RNA$genes <- row.names(S100A9.markers.RNA)
cluster3.markers <- S100A9.markers.RNA

p <- ggplot(cluster3.markers, aes(x=avg_logFC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("Deseq Result on S100A9 DC Control/PAH DEGs")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC < -0.17 & p_val_adj < 0.001), color="#8C1515")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC > 0.17 & p_val_adj < 0.001), color="#09425A")
p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC < -0.3  & p_val_adj < 0.001), aes(label=genes), color="#8C1515")
p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC >0.3 & p_val_adj < 0.001), aes(label=genes), color="#09425A")
p

write.csv(file = "./03.Celltypeanalysis/S100A9_DEG.csv", S100A9.markers.RNA)
ggsave2(filename = "./03.Celltypeanalysis/S100A9_DEG.pdf", p)

#DC.ANXA1
DefaultAssay(DC.ANXA1) <- "RNA"
Idents(DC.ANXA1) <- "condition"
DefaultAssay(DC.ANXA1)
Idents(DC.ANXA1)
table(DC.ANXA1$condition)

ANXA1.markers.RNA <- FindMarkers(DC.ANXA1, ident.1 = "CONT", ident.2 = "PAH",
                                  min.pct = 0.01,test.use="DESeq2", assay="RNA",logfc.threshold=0.3, max.cells.per.ident =300)

ANXA1.markers.RNA
ANXA1.markers.RNA$genes <- row.names(ANXA1.markers.RNA)
cluster3.markers <- ANXA1.markers.RNA

p <- ggplot(cluster3.markers, aes(x=avg_logFC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("Deseq Result on ANXA1 DC Control/PAH DEGs")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC < -0.17 & p_val_adj < 0.001), color="#8C1515")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC > 0.17 & p_val_adj < 0.001), color="#09425A")
p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC < -0.3  & p_val_adj < 0.001), aes(label=genes), color="#8C1515")
p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC >0.3 & p_val_adj < 0.001), aes(label=genes), color="#09425A")
p

write.csv(file = "./03.Celltypeanalysis/ANXA1_DEG.csv", ANXA1.markers.RNA)
ggsave2(filename = "./03.Celltypeanalysis/ANXA1_DEG.pdf", p)

#DC.IFIT
DefaultAssay(DC.IFIT) <- "RNA"
Idents(DC.IFIT) <- "condition"
DefaultAssay(DC.IFIT)
Idents(DC.IFIT)
table(DC.IFIT$condition)

IFIT.markers.RNA <- FindMarkers(DC.IFIT, ident.1 = "CONT", ident.2 = "PAH",
                                 min.pct = 0.01,test.use="DESeq2", assay="RNA",logfc.threshold=0.3, max.cells.per.ident =230)

IFIT.markers.RNA
IFIT.markers.RNA$genes <- row.names(IFIT.markers.RNA)
cluster3.markers <- IFIT.markers.RNA

p <- ggplot(cluster3.markers, aes(x=avg_logFC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("Deseq Result on IFIT DC Control/PAH DEGs")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC < -0.17 & p_val_adj < 0.001), color="#8C1515")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC > 0.17 & p_val_adj < 0.001), color="#09425A")
p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC < -0.3  & p_val_adj < 0.001), aes(label=genes), color="#8C1515")
p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC >0.3 & p_val_adj < 0.001), aes(label=genes), color="#09425A")
p

write.csv(file = "./03.Celltypeanalysis/IFIT_DEG.csv", IFIT.markers.RNA)
ggsave2(filename = "./03.Celltypeanalysis/IFIT_DEG.pdf", p)

#DC.CEAC9_A
DefaultAssay(DC.CEAC9_A) <- "RNA"
Idents(DC.CEAC9_A) <- "condition"
DefaultAssay(DC.CEAC9_A)
Idents(DC.CEAC9_A)
table(DC.CEAC9_A$condition)

CEAC9_A.markers.RNA <- FindMarkers(DC.CEAC9_A, ident.1 = "CONT", ident.2 = "PAH",
                                min.pct = 0.01,test.use="DESeq2", assay="RNA",logfc.threshold=0.3, max.cells.per.ident =170)

CEAC9_A.markers.RNA
CEAC9_A.markers.RNA$genes <- row.names(CEAC9_A.markers.RNA)
cluster3.markers <- CEAC9_A.markers.RNA

p <- ggplot(cluster3.markers, aes(x=avg_logFC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("Deseq Result on CEAC9_A DC Control/PAH DEGs")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC < -0.17 & p_val_adj < 0.001), color="#8C1515")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC > 0.17 & p_val_adj < 0.001), color="#09425A")
p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC < -0.3  & p_val_adj < 0.001), aes(label=genes), color="#8C1515")
p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC >0.3 & p_val_adj < 0.001), aes(label=genes), color="#09425A")
p

write.csv(file = "./03.Celltypeanalysis/CEAC9_A_DEG.csv", CEAC9_A.markers.RNA)
ggsave2(filename = "./03.Celltypeanalysis/CEAC9_A_DEG.pdf", p)

#DC.CEAC9_R
DefaultAssay(DC.CEAC9_R) <- "RNA"
Idents(DC.CEAC9_R) <- "condition"
DefaultAssay(DC.CEAC9_R)
Idents(DC.CEAC9_R)
table(DC.CEAC9_R$condition)

CEAC9_R.markers.RNA <- FindMarkers(DC.CEAC9_R, ident.1 = "CONT", ident.2 = "PAH",
                          min.pct = 0.01,test.use="DESeq2", assay="RNA",logfc.threshold=0.3, max.cells.per.ident =250)

CEAC9_R.markers.RNA
CEAC9_R.markers.RNA$genes <- row.names(CEAC9_R.markers.RNA)
cluster3.markers <- CEAC9_R.markers.RNA

p <- ggplot(cluster3.markers, aes(x=avg_logFC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("Deseq Result on CEAC9_R DC Control/PAH DEGs")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC < -0.17 & p_val_adj < 0.001), color="#8C1515")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC > 0.17 & p_val_adj < 0.001), color="#09425A")
p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC < -0.3  & p_val_adj < 0.001), aes(label=genes), color="#8C1515")
p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC >0.3 & p_val_adj < 0.001), aes(label=genes), color="#09425A")
p

write.csv(file = "./03.Celltypeanalysis/CEAC9_A_DEG.csv", CEAC9_R.markers.RNA)
ggsave2(filename = "./03.Celltypeanalysis/CEAC9_A_DEG.pdf", p)


#DC.2
DefaultAssay(DC.CEAC9) <- "RNA"
Idents(DC.CEAC9) <- "condition"
DefaultAssay(DC.CEAC9)
Idents(DC.CEAC9)
table(DC.CEAC9$condition)

CEAC9.markers.RNA <- FindMarkers(DC.CEAC9, ident.1 = "CONT", ident.2 = "PAH", min.pct = 0.01,
                                  test.use="DESeq2",logfc.threshold= 0.3, assay="RNA", max.cells.per.ident = 400)
CEAC9.markers.RNA
CEAC9.markers.RNA$genes <- row.names(CEAC9.markers.RNA)
cluster3.markers <- CEAC9.markers.RNA

p <- ggplot(cluster3.markers, aes(x=avg_logFC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("Deseq Result on DC2 Control/PAH DEGs")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC < -0.17 & p_val_adj < 0.001), color="#8C1515")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC > 0.17 & p_val_adj < 0.001), color="#09425A")
p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC < -0.3  & p_val_adj < 0.001), aes(label=genes), color="#8C1515")
p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC >0.3 & p_val_adj < 0.001), aes(label=genes), color="#09425A")
p

write.csv(file = "./03.Celltypeanalysis/DC2_DEG.csv", CEAC9.markers.RNA)
ggsave2(filename = "./03.Celltypeanalysis/DC2_DEG.pdf", p)


#DC.1
DefaultAssay(DC1) <- "RNA"
Idents(DC1) <- "condition"
DefaultAssay(DC1)
Idents(DC1)
table(DC1$condition)

DC1.markers.RNA <- FindMarkers(DC1, ident.1 = "CONT", ident.2 = "PAH", min.pct = 0.01,
                                 test.use="DESeq2",logfc.threshold= 0.3, assay="RNA", max.cells.per.ident = 3000)
DC1.markers.RNA
DC1.markers.RNA$genes <- row.names(DC1.markers.RNA)
cluster3.markers <- DC1.markers.RNA

p <- ggplot(cluster3.markers, aes(x=avg_logFC, y=-log(p_val_adj))) + geom_point() + theme_set(theme_cowplot(12))
p <- p + xlab("Average Log Fold Change") + ylab("Negative Log Adjusted P Value") + ggtitle("Deseq Result on DC1 Control/PAH DEGs")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC < -0.17 & p_val_adj < 0.001), color="#8C1515")
p <- p + geom_point(data=subset(cluster3.markers, avg_logFC > 0.17 & p_val_adj < 0.001), color="#09425A")
p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC < -0.3  & p_val_adj < 0.001), aes(label=genes), color="#8C1515")
p <- p + geom_label_repel(data=subset(cluster3.markers, avg_logFC >0.3 & p_val_adj < 0.001), aes(label=genes), color="#09425A")
p <- p + xlim(-1,1)
p

write.csv(file = "./03.Celltypeanalysis/DC1_DEG.csv", DC1.markers.RNA)
ggsave2(filename = "./03.Celltypeanalysis/DC1_DEG.pdf", p)


#UMAP for single gene
DefaultAssay(DC.integrated) <- "RNA"
p=FeaturePlot(object = DC.integrated, features = "GPR183",pt.size=0.4, sort.cell = T, combine = T, split.by = "condition")
p

#RidgePlot for single gene
p2 <- RidgePlot(object = DC.integrated, features = "IFITM3", pt.size=0, group.by = "subject")
p2


#VlnPlot for single gene
p31 <- VlnPlot(object = DC.integrated, features = c("PF4","ANXA1", "ACTB", "GPR183", "CXCR4", "CD1C", "GAPDH", "RPS26", "IFITM3"), ncol=2, combine = T, pt.size=0.01, group.by = "subject", assay="SCT")
p31
p32 <- VlnPlot(object = DC.integrated, features = c("GPR183", "SMDT1", "CD52", "CCR7", "CCL22","CCR5"), ncol=2, combine = T, pt.size=0.01, group.by = "subject", assay="SCT")
p32

p32 <- VlnPlot(object = DC.integrated, features = markers1, ncol=4, combine = T, pt.size=0.01, group.by = "subject", assay="SCT")
p32

#Plot for Centain Marker groups
VlnPlot(object = DC.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mito", "n.exp.hkgenes"), ncol=2, combine = T, pt.size=0.01, group.by = "subject")

markers1 <- grep(pattern = "^CXCL", x = rownames(x = DC.integrated@assays$RNA), value = TRUE)
markers1 <- grep(pattern = "^SDF", x = rownames(x = DC.integrated@assays$RNA), value = TRUE)
markers1 <- grep(pattern = "^CXCR", x = rownames(x = DC.integrated@assays$RNA), value = TRUE)

markers1 <- grep(pattern = "^CCL", x = rownames(x = DC.integrated@assays$RNA), value = TRUE)
markers1 <- grep(pattern = "^CCR", x = rownames(x = DC.integrated@assays$RNA), value = TRUE)

markers1 <- grep(pattern = "^TLR", x = rownames(x = DC.integrated@assays$RNA), value = TRUE)

markers1 <- grep(pattern = "^IL", x = rownames(x = DC.integrated@assays$RNA), value = TRUE)
markers1



DC1.markers.RNA[DC1.markers.RNA$genes == "IL17RA",]


#####added Sep.2020
DefaultAssay(DC.integrated) <- "SCT"
FeaturePlot(object = DC.integrated, features = "LYZ", pt.size = 0.01, split.by = "condition", order=T) 

DC.markers <-  c("IDO1", "HAVCR2", "CD163", "SAMHD1", "STAT1", "STAT3")

for(gene in DC.markers){
  print(gene)
  p=FeaturePlot(object = DC.integrated, features = c(gene), pt.size = 0.01, split.by = "condition", order=T) 
  ggsave(p,filename = paste0("./DCmarkers/Marker.",gene,".pdf"),height=3,width=6, dpi=300)
}


p <- VlnPlot(DC.integrated, features = c("IDO1", "HAVCR2", "CD163", "SAMHD1", "STAT1", "STAT3"), assay = "SCT", ncol=3, pt.size = 0.01, combine = T)
p


markers.to.plot <- c("CLEC9A", "C1orf54", "CADM1", "CAMK2D", "XCR1", "S100A9", "S100A8", "VCAN",
                     "LYZ", "ANXA1", "CD1C", "FCER1A","ISG15", "CLEC10A", "FTL", "SERPINA1", "LST1", "HLA-DQA1", 
                     "AXL", "PPP1R14A", "IFITM3", "IFITM2", "IFIT1", "IFI6", "DUSP1", "DUSP2", "CD14", "GAPDH", 
                     "IDO1", "HAVCR2", "CD163", "SAMHD1", "STAT1", "STAT3")

pDotPlot <- DotPlot(DC.integrated, features = rev(markers.to.plot), cols = c("grey", "red"), dot.scale = 8, assay = "integrated", col.min = 0.01, scale.min=10) + RotatedAxis()
pDotPlot 

pDotPlot2 <- DotPlot(DC.integrated, features = rev(markers.to.plot), cols = c("grey", "red"), dot.scale = 8, assay = "SCT", col.min = 0.01, scale.min=10) + RotatedAxis()
pDotPlot2 


