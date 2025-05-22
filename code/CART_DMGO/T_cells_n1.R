# Set seed for reproducibility
set.seed(123)

# Define file paths for input data and output results
data_loc <- "D:/SURF_2/Shared/Dream3DLab (Groupfolder)/1.Projects/BRO_BrainOrganoids/3.Analysis/Github repo/Data/"
res_path <- "D:/SURF_2/Shared/Dream3DLab (Groupfolder)/1.Projects/BRO_BrainOrganoids/3.Analysis/Github repo/Repeat_output/"

# Create output directory if it doesn't exist
if (!dir.exists(res_path)) {
  dir.create(res_path, recursive = TRUE)
}

# Load required libraries
#library(SeuratDisk)
library(Seurat)
library(ggplot2)
library(pheatmap)
library(VISION)
library(readxl)
library(dplyr)

# Load pre-processed T cell Seurat object (cycling cells removed)
seurat_T_cells_BRO <- readRDS(paste0(data_loc, "All_Tcells_Cycling_Removed.rds"))

# Load signature gene set (Serial killer cells from Supplementary Table 5)
Supplementary_Table_5 <- read_excel(paste0(data_loc, "41587_2022_1397_MOESM7_ESM.xlsx"))
Supplementary_Table_5$value <- 1  # Set all genes to weight 1
sig_killer_temp <- setNames(as.list(Supplementary_Table_5$value), Supplementary_Table_5$Gene)
sig_killer <- createGeneSignature(name = "TEG serial killer", sigData = sig_killer_temp)

# Run VISION on Seurat object with the killer signature
mySignatures <- c(sig_killer)
vision.obj <- Vision(seurat_T_cells_BRO, assay = "RNA", signatures = mySignatures)
vision.obj <- calcSignatureScores(vision.obj)

# Extract signature scores and UMAP coordinates
score_sig <- vision.obj@SigScores
umap <- getProjections(vision.obj)[["Seurat_umap"]]

# Cap outliers in signature scores at 1st and 99th percentiles
sigScore2 <- score_sig
for (i in 1:ncol(sigScore2)) {
  a <- quantile(sigScore2[, i], 0.99)
  b <- quantile(sigScore2[, i], 0.05)
  sigScore2[sigScore2[, i] > a, i] <- a
  sigScore2[sigScore2[, i] < b, i] <- b
}

# Generate UMAP plots colored by signature expression
signatures_to_plot <- colnames(sigScore2)
pdf(paste0(res_path, "/TEG_Own_signatures_", Sys.Date(), ".pdf"))
for (sig in signatures_to_plot) {
  df2 <- data.frame(vision.obj@metaData, umap1 = umap[,1], umap2 = umap[,2], signature = sigScore2[,sig])
  df2 <- df2 %>% mutate(signature = ifelse(signature < quantile(signature, 0.05), NA, signature))
  df2 <- subset(df2, !is.na(signature))
  
  print(
    ggplot(df2) +
      geom_point(aes(x = umap1, y = umap2, color = signature), size = 1.5, alpha = 0.3) +
      theme_bw() +
      ggtitle(sig) +
      scale_color_stepsn(
        colors = c("#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#fddbc7", "#f4a582", "#d6604d", "#b2182b"),
        n.breaks = 80,
        na.value = "white"
      ) +
      theme(
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank()
      )
  )
}
dev.off()

# ========================
# Single Cell Clustering
# ========================

# Find neighbors and clusters
seurat_T_cells_BRO <- FindNeighbors(seurat_T_cells_BRO, reduction = "pca", dims = 1:10)
seurat_T_cells_BRO <- FindClusters(seurat_T_cells_BRO, resolution = 0.45, algorithm = 1)

# Plot cluster UMAP
DimPlot(seurat_T_cells_BRO) + theme(aspect.ratio = 1)

# Rename cluster identities
new_names <- c(
  "0" = "Tund",
  "1" = "Tgzk",
  "2" = "Til2",
  "3" = "Tmi",
  "4" = "Tms",
  "5" = "Tcyt",
  "6" = "Tisg",
  "7" = "Tpr",
  "8" = "Tex"
)
seurat_T_cells_BRO <- RenameIdents(seurat_T_cells_BRO, new_names)

# Reorder clusters
levels(seurat_T_cells_BRO)  # Check current order of clusters

# Define the new order you want, e.g.:
new_order <- c("Tund", "Tmi", "Tms","Tisg" , "Tpr", "Til2","Tex" ,  "Tgzk",  "Tcyt"  )  # specify cluster names/IDs in desired order

# Reorder the Idents factor levels
seurat_T_cells_BRO <- SetIdent(seurat_T_cells_BRO, value = factor(Idents(seurat_T_cells_BRO), levels = new_order))

# Save cluster UMAP
pdf(paste0(res_path, "/n1_UMAP_", Sys.Date(), ".pdf"))
DimPlot(seurat_T_cells_BRO) + theme(aspect.ratio = 1)
dev.off()

# ============================
# Differential Gene Expression
# ============================

# Find DEGs between clusters
all.markers <- FindAllMarkers(seurat_T_cells_BRO)
all.markers_up <- subset(all.markers, p_val_adj < 0.05 & avg_log2FC > 0)
all.markers_dn <- subset(all.markers, p_val_adj < 0.05 & avg_log2FC < 0)

# Save DEG results
write.csv(all.markers, paste0(res_path, "/New_clustering_DEG_renamed_", Sys.Date(), ".csv"))
write.csv(all.markers_up, paste0(res_path, "/New_clustering_DEG_up_renamed_", Sys.Date(), ".csv"))
write.csv(all.markers_dn, paste0(res_path, "/New_clustering_DEG_dn_renamed_", Sys.Date(), ".csv"))

# Save updated Seurat object
saveRDS(seurat_T_cells_BRO, paste0(data_loc, "/T_Cells_BRO_n1_", Sys.Date(), ".rds"))

# ========================
# DotPlot: Signature Genes
# ========================

# Load signature gene list
Signature_genes <- read_excel(paste0(data_loc, "Signature genes T cell clusters_MA.xlsx"))

# Dot plot of signature genes across clusters
pdf(paste0(res_path, "/Dotplot_sig_genes_per_cl_", Sys.Date(), ".pdf"), width = 12, height = 12)
DotPlot(seurat_T_cells_BRO, features = unique(Signature_genes$gene), dot.scale = 5, cols = c("deepskyblue", "coral2")) +
  scale_size_continuous(range = c(1.5, 6)) +
  theme(aspect.ratio = 1, axis.text.x = element_text(size = 7, angle = 90, vjust = 1, hjust = 1))
dev.off()

# ====================
# Killer Gene Dot Plot
# ====================
killer_genes <- c("GZMK", "GZMB", "PRF1", "IFNG")
pdf(paste0(res_path, "/Dotplot_Killer_genes_per_cl_", Sys.Date(), ".pdf"), width = 12, height = 12)
DotPlot(seurat_T_cells_BRO, features = killer_genes, dot.scale = 5, cols = c("deepskyblue", "coral2")) +
  scale_size_continuous(range = c(1.5, 6)) +
  theme(aspect.ratio = 1, axis.text.x = element_text(size = 7, angle = 90, vjust = 1, hjust = 1))
dev.off()

# ========================
# Exhaustion Gene Dot Plot
# ========================
exhaustion_genes <- c("LAG3", "HAVCR2", "TIGIT", "SELPLG", "PRDM1")
pdf(paste0(res_path, "/Dotplot_exhaus_genes_per_cl_", Sys.Date(), ".pdf"), width = 12, height = 12)
DotPlot(seurat_T_cells_BRO, features = exhaustion_genes, dot.scale = 5, cols = c("deepskyblue", "coral2")) +
  scale_size_continuous(range = c(1.5, 6)) +
  theme(aspect.ratio = 1, axis.text.x = element_text(size = 7, angle = 90, vjust = 1, hjust = 1))
dev.off()


# ========================
# Project the T cells signatures that we find on  a Pan-cancet T cell atlas including brain tumors https://pubmed.ncbi.nlm.nih.gov/37248301/
# ========================
# 1. Load the Wang Pan-cancer T cell atlas Seurat object and set identity to cell types
seurat_Wang <- readRDS(paste0(data_loc, "Chu_CD8.rds"))
seurat_Wang <- SetIdent(object = seurat_Wang, value = "cell.type")

# Visualize UMAP of Wang dataset
DimPlot(seurat_Wang, reduction = "umap") + theme(aspect.ratio = 1)

# 2. Import our DEG (differentially expressed genes) signatures from clustering results
New_clustering_DEG <- read.csv(paste0(res_path, "N1_DEG.csv"))
# Filter to only highly significant DEGs
New_clustering_DEG <- subset(New_clustering_DEG, p_val_adj < 1e-5)

# 3. Function to create gene signatures per cluster for VISION
get_sig <- function(cell_type) {
  cell_pop_df <- subset(New_clustering_DEG, cluster == cell_type)
  # Assign +1 if upregulated, -1 if downregulated
  cell_pop_df$value <- ifelse(cell_pop_df$avg_log2FC > 0, 1, -1)
  sig_temp <- setNames(as.list(cell_pop_df$value), cell_pop_df$gene)
  sig <- createGeneSignature(name = cell_type, sigData = sig_temp)
  return(sig)
}

# Create a list of gene signatures, one per cluster
mySignatures <- lapply(unique(New_clustering_DEG$cluster), get_sig)
names(mySignatures) <- unique(New_clustering_DEG$cluster)

# 4. Run VISION analysis on Wang dataset using the created signatures
vision.obj <- Vision(seurat_Wang, assay = "RNA", dimRed = "umap", signatures = mySignatures)
vision.obj <- calcSignatureScores(vision.obj)

# 5. Process signature scores to remove outliers (cap at 5th and 99th percentiles)
score_sig <- vision.obj@SigScores
umap <- getProjections(vision.obj)[["Seurat_umap"]]

sigScore_BRO <- score_sig
for (i in seq_len(ncol(sigScore_BRO))) {
  upper_bound <- quantile(sigScore_BRO[, i], 0.99)
  lower_bound <- quantile(sigScore_BRO[, i], 0.05)
  sigScore_BRO[sigScore_BRO[, i] > upper_bound, i] <- upper_bound
  sigScore_BRO[sigScore_BRO[, i] < lower_bound, i] <- lower_bound
}

# Normalize each signature by subtracting the minimum value
sigScore_BRO <- apply(sigScore_BRO, 2, function(col) col - min(col))

# 6. Define which signatures to plot (first 9 in this case)
signatures_to_plot_BRO <- colnames(sigScore_BRO)[1:9]

# 7. Add signature scores as metadata columns in Seurat object
seurat_Wang$unresponsive <- sigScore_BRO[, "unresponsive"]
seurat_Wang$midterm_activated_HLA_expressing <- sigScore_BRO[, "mid-term activated_HLA-expressing"]
seurat_Wang$IL_2_responsive <- sigScore_BRO[, "IL-2 responsive"]
seurat_Wang$migrating_interacting <- sigScore_BRO[, "migrating_interacting"]
seurat_Wang$metabolically_stressed <- sigScore_BRO[, "metabolically stressed"]
seurat_Wang$sustained_cytotoxic <- sigScore_BRO[, "sustained cytotoxic"]
seurat_Wang$IFN_I_responsive_ISG_expressing <- sigScore_BRO[, "IFN-I responsive_ISG-expressing"]
seurat_Wang$proliferating <- sigScore_BRO[, "proliferating"]
seurat_Wang$cytotoxic_exhausted <- sigScore_BRO[, "cytotoxic-exhausted"]

# 8. Save plots to PDF
pdf(paste0(res_path, "/BRO_signatures_on_Wang_", Sys.Date(), ".pdf"))

# Plot UMAP with clusters
print(DimPlot(seurat_Wang) + theme(aspect.ratio = 1))

# Prepare data for violin plots
sigScore_BRO_df <- reshape2::melt(sigScore_BRO)

# Violin plots of signature score distributions
print(
  ggplot(sigScore_BRO_df, aes(x = Var2, y = value, color = Var2)) +
    geom_hline(yintercept = 0) +
    geom_violin(trim = FALSE) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
)

# Dot plot of signature scores on Seurat object
print(
  DotPlot(seurat_Wang, features = signatures_to_plot_BRO, scale = FALSE) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
)

# 9. Plot individual signatures on UMAP (filtered for low values)
for (sig in signatures_to_plot_BRO) {
  df2 <- data.frame(vision.obj@metaData, umap1 = umap[, 1], umap2 = umap[, 2], signature = sigScore_BRO[, sig])
  df2 <- df2 %>%
    mutate(signature = ifelse(signature < quantile(signature, 0.05), NA, signature)) %>%
    filter(!is.na(signature))
  
  print(
    ggplot(df2) +
      geom_point(aes(x = umap1, y = umap2, color = signature), size = 1.5, alpha = 0.3) +
      theme_bw() +
      ggtitle(sig) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank()
      ) +
      scale_color_gradientn(
        limits = c(0, 2),
        breaks = seq(0, 2, by = 0.2),
        labels = seq(0, 2, by = 0.2),
        colors = c("#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#fddbc7", "#f4a582", "#d6604d", "#b2182b"),
        na.value = "white"
      )
  )
  print("done")
}

dev.off()

