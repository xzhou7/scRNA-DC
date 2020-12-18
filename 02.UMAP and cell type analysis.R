#Final script for Toshie's data
##Authro: Xin Zhou
#Created: 2019-07-31

#last updated 2019-02-04

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

setwd("~/Box/XinZhouFiles/Projects/ZXP1_PAH/scRNA_DC_Toshie (Rabinovitch lab)/Figure result/Final analysis/")

DC.integrated_all <- readRDS(file="~/Box/XinZhouFiles/Projects/ZXP1_PAH/scRNA_DC_Toshie (Rabinovitch lab)/Figure result/Final analysis/01.Data Process/DCsuratObect.Final.RData")
DC.integrated_all <- DC.integrated
#p <- VlnPlot(DC.integrated_all, features = c("CD1C", "HLA-DQA1", "LYZ","ITGAX", "IFITM2", "CLEC9A", "CD14","TLR4", "CD3D"), assay = "integrated", ncol=3, pt.size = 0.01, combine = T)
#p

DC.integrated <- subset(DC.integrated, ident = c(0:7, 11))
DC.integrated
DimPlot(DC.integrated, label = T)
DC.clean.list <- CellSelector(DimPlot(DC.integrated))
DC.integrated <- subset(DC.integrated, cells = DC.clean.list)


markers.to.plot <- c("CLEC9A", "C1orf54", "CADM1", "CAMK2D", "BATF3", "XCR1", "S100A9", "S100A8", "VCAN", "EEF1A1", "TPT1", "CD163","GPR183",
                     "LYZ", "ANXA1", "CD1C", "FCER1A","S100A4", "SELL","CLEC10A", "FTL", "SERPINA1", "LST1", "HLA-DQA1","HLA-DRB1",
                     "AXL", "PPP1R14A", "IFITM3", "IFITM2", "IFIT1", "DUSP1", "DUSP2", "CD14", "GAPDH")

pDotPlot <- DotPlot(DC.integrated, features = rev(markers.to.plot), cols = c("grey", "red"), dot.scale = 8, assay = "integrated", col.min = 0.01, scale.min=10) + RotatedAxis()
pDotPlot 

pUMAP <- DimPlot(DC.integrated, reduction = "umap", label = T, ncol=4) + scale_color_aaas()
pUMAP

all.genes <- rownames(DC.integrated)
DC.integrated <- ScaleData(DC.integrated, features = all.genes)
DC.integrated <- FindNeighbors(DC.integrated, dims = 1:20)
DC.integrated <- FindClusters(DC.integrated, resolution = 0.7)
DimPlot(DC.integrated)

DefaultAssay(DC.integrated) <- "integrated"
DC.cluster.Marker <- FindAllMarkers(DC.integrated, logfc.threshold = 0.25,test.use = "LR",only.pos = T)
write.csv(file = "./DC.Cluster.Marker.csv", DC.cluster.Marker)

DC.integrated <- RenameIdents(DC.integrated,  
                              "0"="EEF1A1_DC", 
                              "1"="CD1C_A", 
                              "2"="S100A9_DC", 
                              "3"= "FCER1A_DC",
                              "4"="CEAC9_DC_R", 
                              "5"= "ANXA1_DC", 
                              "6"="CEAC9_DC_A", 
                              "7"="IFIT_DC", 
                              "8"="LYZ_DC", 
                              "9"="Intermediate_Cell")

DC.integrated$celltype <- factor(DC.integrated$celltype, levels = c("S100A9_DC","ANXA1_DC","EEF1A1_DC","FCER1A_DC","IFIT_DC",
                                                                    "LYZ_DC","CD1C_A","Intermediate_Cell","CEAC9_DC_R","CEAC9_DC_A"))

Idents(DC.integrated) <- DC.integrated$celltype 
DC.integrated$celltype <- Idents(DC.integrated)

pUMAP <- DimPlot(DC.integrated, reduction = "umap", label = T, ncol=4)# + scale_color_aaas(alpha = 0.7)
pUMAP
ggsave2(filename = "./F1.UMAP_new.pdf", pUMAP, scale=0.5)

pDotPlot <- DotPlot(DC.integrated, features = rev(markers.to.plot), cols = c("grey", "red"), dot.scale = 8, assay = "integrated", col.min = 0.01, scale.min=10) + RotatedAxis()
pDotPlot
ggsave2(filename = "./F2.markersdot.pdf", pDotPlot, scale=0.77)

count2 <- group_by(DC.integrated@meta.data, subject, celltype) %>%
  dplyr::summarise(count = n()) %>% 
  mutate(freq = count / sum(count))
colnames(count2)[2] <- "group"
count2$condition[count2$subject %in% c("S1", "S2", "S3", "S4")]<- "CONT"
count2$condition[count2$subject %in% c("S5", "S6", "S7")]<- "PAH"
write.csv(count2, file="./percentagetable.csv")

count2$group
p2 <- ggplot(count2, aes(x=condition, y=freq, color=condition)) +
  geom_point(position = position_dodge(width=0.75)) +
  geom_boxplot(alpha = 0.1, width=0.75)
p2 <- p2 + scale_y_continuous(labels = percent_format()) + scale_color_lancet()
p2 <- p2 + facet_wrap(.~group, scales="free_y", ncol=5)
p2 <- p2 + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2 <- p2 + stat_compare_means(label =  "p.format", method = "t.test") + theme_cowplot()
p2

p3 <- ggplot(count2, aes(x=condition, y=count, fill=group)) + geom_bar(stat="identity", position = "fill") + scale_y_continuous(labels=scales::percent)
p3

ggsave2(filename = "./F3.percentcell.pdf", p2, scale=0.77)
ggsave2(filename = "./F3.2.percentcell.pdf", p3, width = 4, height = 6, scale=0.77)

cDC1count <- subset(count2,group=="CEAC9_DC_R")
cDC1count$Ac <- subset(count2,group=="CEAC9_DC_A")$freq
cDC1count$cDC1 <- cDC1count$freq + cDC1count$Ac 
cDC1count$cDC2 <- 1 - cDC1count$cDC1
cDC1count$cDC2_cDC1_Ratio <- cDC1count$cDC2/cDC1count$cDC1

t.test(cDC1count$cDC1[1:4], cDC1count$cDC1[5:7])
t.test(cDC1count$cDC2_cDC1_Ratio[1:4], cDC1count$cDC2_cDC1_Ratio[5:7])

DC.markers <- FindAllMarkers(DC.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

marker.table <- DC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(marker.table, file="./marker.table10.csv")

#antigen presenting ability of DCs
MHC2gene <-  c("HLA-DRA","HLA-DRB5","HLA-DRB1","HLA-DQA1", "HLA-DQB1","HLA-DQA2",
               "HLA-DQB2","HLA-DOB","HLA-DMB","HLA-DMA", "HLA-DPA1","HLA-DPB1")
p3 <- VlnPlot(DC.integrated, features = MHC2gene, assay = "integrated", ncol=4, pt.size = 0.00, split.by = "condition")
p3
ggsave2(filename = "./F4.HLA.pdf", p3, scale=1)

#featureplot
DefaultAssay(DC.integrated) <- "RNA"
p=FeaturePlot(object = DC.integrated, features = c("BMPR2"),cols=c("lightgrey", "#931F21"),pt.size=0.4, sort.cell = T, combine = T, split.by = "condition")
p

DC.BMPR2 <- subset(DC.integrated, BMPR2 > 0)
count3 <- group_by(DC.BMPR2@meta.data, subject, celltype) %>%
  dplyr::summarise(count = n()) %>% 
  mutate(freq = count / sum(count))
colnames(count3)[2] <- "group"
count3$condition[count3$subject %in% c("S1", "S2", "S3", "S4")]<- "CONT"
count3$condition[count3$subject %in% c("S5", "S6", "S7")]<- "PAH"

p2 <- ggplot(count3, aes(x=condition, y=count,color=condition)) +
  geom_point(position = position_dodge(width=0.75)) +
  geom_boxplot(alpha = 0.1, width=0.75)
p2 <- p2 + facet_wrap(.~group, scales="free_y", ncol=5)
p2

#check interferon signal
gmarkersIFN <- grep(pattern = "^IFN", x = rownames(x = DC.integrated@assays$RNA), value = TRUE)

vplot4 <- VlnPlot(object = DC.integrated, features = gmarkersIFN, ncol=3, combine = T, pt.size=0.00, group.by = "condition", cols =c("#1F77B4FF","#FF7F0EFF"), assay = "SCT" ) 
vplot4

ggsave2(filename = "./F5.IFN.pdf", vplot4, scale=1)


#IL17 related
marksgrepIL17 <- append(grep(pattern = "^IL17", x = rownames(x = DC.integrated@assays$RNA), value = TRUE), 
                        grep(pattern = "^ROR", x = rownames(x = DC.integrated@assays$RNA), value = TRUE))
vplot3 <- VlnPlot(object = DC.integrated, features = marksgrepIL17, ncol=3, combine = T, pt.size=0.00, group.by = "condition", cols =c("#1F77B4FF","#FF7F0EFF"), assay = "SCT" ) 
vplot3

ggsave2(filename = "./F6.IL17.pdf", vplot3, scale=1)

DC.IL17RA <- subset(DC.integrated, IL17RA > 0)
DC.IL17RA

count4.1 <- group_by(DC.IL17RA@meta.data, subject) %>%
  dplyr::summarise(count = n()) %>% 
  mutate(freq = count / sum(count))
count4.2 <- group_by(DC.integrated@meta.data, subject) %>%
  dplyr::summarise(count = n()) %>% 
  mutate(freq = count / sum(count))

count4 <- cbind(count4.1, count4.2)
colnames(count4) <- c("subject","count1","freq1","subject2","count2","freq2" )
count4$Perc <- count4$count1/count4$count2

count4$condition[count4$subject %in% c("S1", "S2", "S3", "S4")]<- "CONT"
count4$condition[count4$subject %in% c("S5", "S6", "S7")]<- "PAH"

p5 <- ggplot(count4, aes(x=condition, y=Perc)) +
  geom_point(position = position_dodge(width=0.75)) +
  geom_boxplot(alpha = 0.1, width=0.75) + ggtitle("Percent of IL17RA positive DC")
p5 <- p5 + theme_cowplot() +stat_compare_means(label =  "p.format", method = "t.test")
p5
ggsave2(filename = "./F7.IL17_perfect.pdf", p5, scale=0.3)


linmarkers <- c("ISG15", "CD83", "AXL", "LYZ", "CLEC9A", "C1orf54", "CADM1", "CAMK2D", "XCR1", "S100A9", "S100A8", "VCAN", 
                "LYZ", "ANXA1", "CD1C", "FCER1A", "CLEC10A", "FTL", "SERPINA1", "LST1", "HLA-DQA1", 
                "AXL", "PPP1R14A", "IFITM3", "IFITM2", "IFIT1", "DUSP1", "DUSP2", "CD14", "GAPDH")

DefaultAssay(DC.integrated) <- "SCT"
for(gene in linmarkers){
  print(gene)
  p=FeaturePlot(object = DC.integrated, features = c(gene), split.by="condition",cols=c("lightgrey", "#931F21"), pt.size=1.0, sort.cell=T)+
    theme(axis.title.x = element_text( vjust=0.5, size=24,face="bold"),axis.text.x  = element_text( vjust=0.5, size=24,face="bold"),
          axis.title.y = element_text( vjust=0.5, size=24,face="bold"),axis.text.y  = element_text( vjust=0.5, size=24,face="bold"),
          legend.text = element_text( size = 20, face = "bold"), plot.title = element_text(lineheight=.8, size=40, face="bold.italic"))
  ggsave2(p,filename = paste0("./features/",gene,".pdf"),height=4,width=12)
}


saveRDS(DC.integrated, file="~/Box/XinZhouFiles/Projects/ZXP1_PAH/scRNA_DC_Toshie (Rabinovitch lab)/Figure result/Final analysis/02.UMAP and Cell Identity/DCsuratObect.Final02.RData")





