# Differential Gene Expression between MeV+/MeV- cells 
# Author: Andy Chan
# 13 June 2023

# Set working directory ---------------------------------------

setwd("~/Johns Hopkins Immunology/Diane Lab/BCR_clonal_analysis_pipelines_Mar_2022/cluster_files_1Mar23/")

# Open libraries -------------------------
# Can always run sessionInfo() to look at and save the version of the softward used

library("dplyr")
library("tidyr")
library("ggpubr")
library("Seurat")
library("tibble")
library("patchwork")
library("cowplot")
library("ggplot2")
library("fields")
library('escape')
library('dittoSeq')
library('DESeq2')
library('pheatmap')
library("EnhancedVolcano")
library("stringr")
library("stringi")
library("ggrepel")

# Load data for input ---------------------------------------
# vgm_mev is the B cell only Seurat object after filtering for cells with both GEX and VDJ information on each single cell

N_barcode_raw <- read.csv("~/Johns Hopkins Immunology/Diane Lab/BCR_clonal_analysis_pipelines_Mar_2022/cluster_files_1Mar23/N_barcodes_new.txt")

vgm_mev <- readRDS("~/Johns Hopkins Immunology/Diane Lab/BCR_clonal_analysis_pipelines_Mar_2022/cluster_files_1Mar23/vgm_mev_GEX_only.RDS")

# Data transformation ---------------------------------------

# Transform N_barcode_raw
N_barcode <- N_barcode_raw %>%
  mutate (n_positive_b_cells = "n_positive_b_cell") %>%
  mutate(orig_barcode = str_sub(Barcode, end=-3)) #here I also transformed the dataset for N positive barcodes

# To compare DGE between N+ cells in D11 (43F) vs uninfected B cells in all samples
vgm_mev@meta.data <- vgm_mev@meta.data %>%
  mutate(n_positive_or_negative = ifelse(orig_barcode %in% N_barcode$orig_barcode & sample_id == "s6", "n_positive_b_cell", "n_negative_b_cell"))

# RUN THIS ONLY IF YOU compare DGE between N+ cells in D11 (43F) vs uninfected B cells in D11 (43F) 
vgm_mev <- subset (vgm_mev, sample_id == "s6")
vgm_mev@meta.data <- vgm_mev@meta.data %>%
  mutate(n_positive_or_negative = ifelse(orig_barcode %in% N_barcode$orig_barcode & sample_id == "s6", "n_positive_b_cell", "n_negative_b_cell"))

table(vgm_mev$n_positive_or_negative)
table(vgm_mev$n_positive_or_negative, vgm_mev$seurat_clusters)

# my_volcano_plot function to look at DGE between groups --------------------------
# Can modify parameters of your liking
my_volcano_plot <- function(clus) {
  clus$diff_expressed <- "NO"
  clus$diff_expressed[clus$avg_log2FC > 0.4 & clus$p_val_adj < 10e-2] <- "UP"
  clus$diff_expressed[clus$avg_log2FC < -0.4 & clus$p_val_adj < 10e-2] <- "DOWN"
  clus$delabel <- NA
  clus$gene_symbol <- rownames(clus)
  clus$delabel[clus$diff_expressed != "NO"] <- clus$gene_symbol[clus$diff_expressed != "NO"]
  
  return (ggplot(data=clus, aes(x=avg_log2FC, y=-log10(p_val_adj), col = diff_expressed, label= delabel)) +
            geom_point(alpha=0.5) + 
            theme_minimal() +
            theme(axis.line = element_line(colour = "black"))+
            theme(axis.text=element_text(size=12),
                  axis.title=element_text(size=14,face="bold"))+
            theme(legend.position="top")+
            geom_text_repel(max.overlaps = Inf) +
            # scale_color_manual(values=c("blue", "black", "grey" , "red")) +
            scale_color_manual(values=c("blue", "grey" , "red")) +
            geom_vline(xintercept=c(-0.4, 0.4), col="black", linetype = "dashed") +
            geom_hline(yintercept=-log10(10e-2), col="red", linetype = "dashed"))
}

Idents(vgm_mev) <- "n_positive_or_negative"

n_positive_or_negative_deg <- FindMarkers(vgm_mev, ident.1 = "n_positive_b_cell", ident.2 = "n_negative_b_cell", test.use="wilcox", min.pct = 0.25, logfc.threshold = log(0.4))

my_volcano_plot(n_positive_or_negative_deg)

# heat map, dot plot and violin plots to compare DGE identified from the volcano plot ------------------------

alldata <- ScaleData(vgm_mev, assay = "RNA")
Seurat :: DoHeatmap(subset(alldata, downsample = 32),assay = "RNA", features = c("HBB", "IFI27", "ENSMMUG00000041453"), label = T, size = 5)+ scale_fill_gradientn(colors = c("yellow", "white", "purple")) 
DotPlot(subset(alldata, downsample = 32), features = c("HBB", "IFI27", "ENSMMUG00000041453"), group.by = "n_positive_or_negative", assay = "RNA") + coord_flip()
VlnPlot(alldata, features = c("HBB", "IFI27", "ENSMMUG00000041453"), ncol = 3, group.by = "n_positive_or_negative", assay = "RNA", pt.size = 0) 

# heat map, dot plot and violin plots to compare top DGE
# need to first create an ordered gene list 

Idents(vgm_mev) <- "n_positive_or_negative"
find_all_markers <- FindAllMarkers(object = vgm_mev, 
                                   only.pos = TRUE,
                                   logfc.threshold = 0.25)
top_genes_lfc_0.4 <- find_all_markers %>%
  arrange(-avg_log2FC) %>%
  group_by(cluster) %>%
  dplyr :: filter (avg_log2FC > 0.40) %>%  # can change to any numbers
  dplyr :: slice(1:15) # can change to any numbers

alldata <- ScaleData(vgm_mev, features = as.character(unique(top_genes_lfc_0.4$gene)), assay = "RNA")
Seurat :: DoHeatmap(subset(alldata, downsample = 32),assay = "RNA", features = as.character(unique(top_genes_lfc_0.4$gene)), label = T, size = 5)+ scale_fill_gradientn(colors = c("yellow", "white", "purple")) 
DotPlot(subset(alldata, downsample = 100), features = rev(as.character(unique(top_genes_lfc_0.4$gene))), group.by = "n_positive_or_negative", assay = "RNA") + coord_flip()
VlnPlot(alldata, features = as.character(unique(top_genes_lfc_0.4$gene)), ncol = 10, group.by = "n_positive_or_negative", assay = "RNA", pt.size = 0)

