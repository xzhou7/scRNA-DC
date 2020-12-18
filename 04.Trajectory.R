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
library("monocle3")
library("MAST")
options(stringsAsFactors = FALSE)

setwd("~/Box/XinZhouFiles/Projects/ZXP1_PAH/scRNA_DC_Toshie (Rabinovitch lab)/Figure result/Final analysis/")

seurat  <- readRDS(file="./02.UMAP and Cell Identity/DCsuratObect.Final02.RData")

nPC <- 15
cluster.res <- 0.6

gene_annotation <- as.data.frame(rownames(seurat@reductions[["pca"]]@feature.loadings), row.names = rownames(seurat@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"

cell_metadata <- as.data.frame(seurat@assays[["RNA"]]@counts@Dimnames[[2]], row.names = seurat@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"

New_matrix <- seurat@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(seurat@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix

cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)

recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition


list_cluster <- seurat@meta.data[[sprintf("ClusterNames_%s_%sPC", cluster.res, nPC)]]
names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

#cds_from_seurat@reduce_dim_aux@listData[["UMAP"]] <-seurat@reductions[["umap"]]@cell.embeddings

cds_from_seurat@preprocess_aux$gene_loadings <- seurat@reductions[["pca"]]@feature.loadings

cds_from_seurat <- preprocess_cds(cds_from_seurat, num_dim = 100)
cds_from_seurat <- reduce_dimension(cds_from_seurat,preprocess_method="PCA",reduction_method = "UMAP")
cds_from_seurat = cluster_cells(cds_from_seurat, resolution=0.0001)

cds_from_seurat <- learn_graph(cds_from_seurat)
cds_from_seurat <- order_cells(cds_from_seurat)

plot_cells(cds_from_seurat, color_cells_by="cluster")

cluster <- plot_cells(cds_from_seurat, 
                         color_cells_by = 'cluster')
cluster

pseudotime <- plot_cells(cds_from_seurat, 
                            color_cells_by = 'pseudotime')

ciliated_genes <- c("XCR1",
                    "S100A9",
                    "IFITM3",
                    "CD1C",
                    "FCER1A",
                    "LYZ")

plot_cells(cds_from_seurat,
           genes=ciliated_genes,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

#sudo time
cds_from_seurat <- order_cells(cds_from_seurat)
plot_cells(cds_from_seurat,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)

plot_cells(cds_from_seurat, label_groups_by_cluster=FALSE,  color_cells_by = "pseudotime")
