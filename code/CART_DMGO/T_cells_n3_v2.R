# Set seed for reproducibility
set.seed(123)

# Define data and result paths
data_loc <- "D:/SURF_2/Shared/Dream3DLab (Groupfolder)/1.Projects/BRO_BrainOrganoids/3.Analysis/Github repo/Data/"
res_path <- "D:/SURF_2/Shared/Dream3DLab (Groupfolder)/1.Projects/BRO_BrainOrganoids/3.Analysis/Github repo/Repeat_output/"

# Create result directory if it does not exist
if (!dir.exists(res_path)) {
  dir.create(res_path, recursive = TRUE)
}

# Load required libraries
library(Seurat)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(ggtrace)
library(viridis)
library(openxlsx)

# -------------------------------
# Define heatmap generation function
generate_heatmap <- function(seurat_object, idents, features, color_palette = c("blue2", "coral1", "coral3")) {
  # Extract cells belonging to specified clusters
  cells <- WhichCells(seurat_object, idents = idents)
  
  # Get cluster annotations for selected cells
  cluster_annotations <- as.data.frame(Idents(seurat_object)[cells])
  colnames(cluster_annotations) <- "Cluster"
  cluster_annotations <- cluster_annotations[order(cluster_annotations$Cluster), , drop = FALSE]
  
  # Sort cells by cluster order
  sorted_cells <- rownames(cluster_annotations)
  
  # Extract normalized expression data for selected features and cells
  expression_data <- GetAssayData(
    object = seurat_object,
    slot = "data",
    assay = "RNA"
  )[features, sorted_cells]
  expression_matrix <- as.matrix(expression_data)
  
  # Log normalization (log(1 + x))
  log_normalized_data <- log1p(expression_matrix)
  
  # Scale log-normalized data by rows (genes)
  scaled_data <- t(apply(log_normalized_data, 1, function(x) (x - mean(x)) / sd(x)))
  
  # Prepare annotations for heatmap columns (cells)
  annotation_col <- cluster_annotations
  
  # Generate heatmap
  pheatmap(
    scaled_data,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    treeheight_row = 0,
    cellwidth = 0.05,
    cellheight = 10,
    annotation_col = annotation_col,
    color = colorRampPalette(color_palette)(100)
  )
}

# -------------------------------
# Import pre-processed datasets
seurat_T_cellBROn3 <- readRDS(paste0(data_loc, "/T_cell_BRO_unprocessed_n3.rds"))
seurat_T_cells_BRO_n1 <- readRDS(paste0(data_loc, "/T_Cells_BRO_n1_processed.rds"))


# -------------------------------
# Prepare metadata for n1 and n3 datasets for integration

# Extract old cluster identities from n1
desired_order <- unique(cell_type_distribution$`seurat_T_cells_BRO_n1@active.ident`)
cell_identn1 <- as.data.frame(seurat_T_cells_BRO_n1@active.ident)
colnames(cell_identn1) <- "old_clusters"
cell_identn1$cellname <- rownames(cell_identn1)
cell_identn1$Sample <- seurat_T_cells_BRO_n1@meta.data$Sample

# Extract and modify n3 metadata
cell_identn3 <- as.data.frame(seurat_T_cellBROn3@meta.data)
cell_identn3_subn1 <- subset(cell_identn3, sequencing_batch == 1)
cell_identn3_subn1$cellname <- rownames(cell_identn3_subn1)
cell_identn3_subn1$cellname <- sub("^[^_]*_[^_]*_[^_]*_", "APC_", cell_identn3_subn1$cellname)
cell_identn3_subn1$cellname <- gsub("_CAR", "", cell_identn3_subn1$cellname)

# Join with n1 identities
library(dplyr)
temp1 <- left_join(cell_identn3_subn1, cell_identn1, by = "cellname")
rownames(temp1) <- rownames(cell_identn3_subn1)

cell_identn3 <- merge(cell_identn3, temp1[c("old_clusters", "Sample")], by = "row.names", all.x = TRUE)
cell_identn3$old_clusters <- as.character(cell_identn3$old_clusters)
cell_identn3$old_clusters[is.na(cell_identn3$old_clusters)] <- "unknown"
cell_identn3$Sample[is.na(cell_identn3$Sample)] <- "unknown"
rownames(cell_identn3) <- cell_identn3$Row.names
cell_identn3$old_clusters <- factor(cell_identn3$old_clusters)

seurat_T_cellBROn3 <- AddMetaData(seurat_T_cellBROn3, metadata = cell_identn3[c("old_clusters", "Sample")])

# -------------------------------
# Add metadata from N2 experiment
seurat_Ncam1_Tcell <- readRDS(paste0(data_loc, "/NCAM_combined_clustering_seurat.rds"))

# Cluster with resolution 0.10 using algorithm 1
seurat_Ncam1_Tcell <- FindClusters(seurat_Ncam1_Tcell, resolution = 0.10, algorithm = 1)

# Plot UMAP by sample and cluster
DimPlot(seurat_Ncam1_Tcell, group.by = "sample")
DimPlot(seurat_Ncam1_Tcell)

# Define treatment groups
seurat_Ncam1_Tcell$treatment <- ifelse(
  Idents(seurat_Ncam1_Tcell) %in% c("1", "2"),
  "no exposure",
  "NA"
)

seurat_Ncam1_Tcell$treatment <- ifelse(
  seurat_Ncam1_Tcell$sample == "BRO_NCAM_1" & seurat_Ncam1_Tcell$treatment == "NA",
  "NCAM +",
  ifelse(
    seurat_Ncam1_Tcell$sample == "BRO_NCAM_2" & seurat_Ncam1_Tcell$treatment == "NA",
    "NCAM -",
    seurat_Ncam1_Tcell$treatment
  )
)

# UMAP plot colored by treatment
DimPlot(seurat_Ncam1_Tcell, reduction = "umap", group.by = "treatment") +
  theme(aspect.ratio = 1)

# Prepare treatment metadata for n3 dataset
cell_identn2 <- as.data.frame(seurat_Ncam1_Tcell@meta.data)
cell_identn2$cellname <- rownames(cell_identn2)

cell_identn3$treatment <- ifelse(cell_identn3$sequencing_batch == "2", NA, cell_identn3$treatment)
cell_identn3 <- cell_identn3 %>%
  mutate(cellname = rownames(cell_identn3)) %>%
  left_join(cell_identn2 %>% select(cellname, treatment), by = "cellname", suffix = c(".old", ".new")) %>%
  mutate(treatment = ifelse(is.na(treatment.old), treatment.new, treatment.old)) %>%
  select(-treatment.old, -treatment.new)

rownames(cell_identn3) <- cell_identn3$cellname
cell_identn3$cellname <- NULL

seurat_T_cellBROn3 <- AddMetaData(seurat_T_cellBROn3, metadata = cell_identn3[c("treatment")])

# -------------------------------
# CCA integration and downstream analysis

seurat_T_cellBROn3.2 <- seurat_T_cellBROn3
DefaultAssay(seurat_T_cellBROn3.2) <- "RNA"

# Normalize and scale data
seurat_T_cellBROn3.2 <- seurat_T_cellBROn3.2 %>%
  NormalizeData() %>%
  ScaleData()

# Identify differentially expressed markers by treatment
test_vf <- SetIdent(seurat_T_cellBROn3.2, value = "treatment")
markers_dif <- FindAllMarkers(test_vf, min.pct = 0.3, logfc.threshold = 1.3)

# Identify variable features per sequencing batch and combine
VariableFeatures(seurat_T_cellBROn3.2) <- split(row.names(seurat_T_cellBROn3.2@meta.data), seurat_T_cellBROn3.2@meta.data$sequencing_batch) %>%
  lapply(function(cells_use) {
    seurat_T_cellBROn3.2[, cells_use] %>%
      FindVariableFeatures(selection.method = "vst", nfeatures = 500) %>%
      VariableFeatures()
  }) %>%
  unlist() %>%
  unique()

# Combine variable features and differential markers
features_var <- c(VariableFeatures(seurat_T_cellBROn3.2), markers_dif$gene) %>% unique()

# Run PCA on selected features
seurat_T_cellBROn3.2 <- RunPCA(seurat_T_cellBROn3.2, features = features_var, npcs = 50, verbose = FALSE)

# Perform CCA integration
seurat_T_cellBROn3.2[["RNA"]] <- split(seurat_T_cellBROn3.2[["RNA"]], f = seurat_T_cellBROn3.2$sequencing_batch)
seurat_T_cellBROn3.2 <- IntegrateLayers(
  seurat_T_cellBROn3.2,
  features = features_var,
  method = CCAIntegration,
  orig.reduction = 'pca',
  scale.layer = "scale.data",
  new.reduction = "integrated.cca"
)
seurat_T_cellBROn3.2 <- JoinLayers(seurat_T_cellBROn3.2, assay = "RNA")

# Find neighbors, clusters, and run UMAP
seurat_T_cellBROn3.2 <- FindNeighbors(seurat_T_cellBROn3.2, reduction = "integrated.cca", dims = 1:30)
seurat_T_cellBROn3.2 <- FindClusters(seurat_T_cellBROn3.2, resolution = 0.75)
seurat_T_cellBROn3.2 <- RunUMAP(seurat_T_cellBROn3.2, dims = 1:30, reduction = "integrated.cca", reduction.name = "umap")

# -------------------------------
# Visualization

DimPlot(seurat_T_cellBROn3.2) +
  theme(aspect.ratio = 1, strip.text = element_text(size = 5))

DimPlot(seurat_T_cellBROn3.2, reduction = "umap", split.by = "sequencing_batch")

FeaturePlot(seurat_T_cellBROn3.2, reduction = "umap", features = c("MKI67", "TOP2A")) + coord_fixed()

DimPlot(seurat_T_cellBROn3.2, reduction = "umap", group.by = "treatment")

DimPlot(seurat_T_cellBROn3.2, reduction = "umap", split.by = "old_clusters", ncol = 3) +
  theme(aspect.ratio = 1, strip.text = element_text(size = 5))

DimPlot(seurat_T_cellBROn3.2, split.by = "sequencing_batch", ncol = 3) +
  theme(aspect.ratio = 1, strip.text = element_text(size = 5))

DimPlot(seurat_T_cellBROn3.2, split.by = "sequencing_batch", group.by = "treatment", ncol = 3) +
  theme(aspect.ratio = 1, strip.text = element_text(size = 5))

# -------------------------------
# Cluster cell type distribution plots

seurat_T_cellBROn3.2$seurat_clusters <- seurat_T_cellBROn3.2@active.ident

# Cell type distribution by cluster and treatment
set.seed(123)
cell_type_distribution <- seurat_T_cellBROn3.2@meta.data %>%
  group_by(treatment) %>%
  ungroup() %>%
  group_by(seurat_clusters, old_clusters, treatment) %>%
  summarise(cell_count = n()) %>%
  filter(old_clusters != "unknown")

ggplot(cell_type_distribution, aes(x = seurat_clusters, y = cell_count, fill = old_clusters)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cluster", y = "Proportion", fill = "Cell Type") +
  theme_minimal() +
  ggtitle("Proportion of Cell Types in Each Cluster") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))

# Treatment proportions by cluster
set.seed(321)
cell_type_distribution <- seurat_T_cellBROn3.2@meta.data %>%
  group_by(treatment) %>%
  sample_n(1000) %>%
  group_by(seurat_clusters, treatment, sequencing_batch) %>%
  summarise(cell_count = n())

ggplot(cell_type_distribution, aes(x = seurat_clusters, y = cell_count, fill = treatment)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cluster", y = "Proportion", fill = "Treatment") +
  theme_minimal() +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)
  ) +
  ggtitle("Proportion of Treatment in Each Cluster")

ggplot(cell_type_distribution, aes(x = treatment, y = cell_count, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Treatment", y = "Proportion", fill = "Cluster") +
  theme_minimal() +
  theme(aspect.ratio = 1,
    axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))
    

    ### Rename clusters based on contingency tables and assign meaningful names
    
    # Define new cluster names mapping
    new_cluster_names <- c(
      "0" = "Tund",
      "1" = "Tex",
      "2" = "Tmi",
      "3" = "Tgzk",
      "4" = "Til2",
      "5" = "Tisg",
      "6" = "Tcyt",
      "7" = "Tpr",
      "8" = "Tms",
      "9" = "Tpr",
      "10" = "non-exposed",
      "11" = "Tpr",
      "12" = "Tpr",
      "13" = "Tpr",
      "14" = "Ths"
    )
    
    # Rename cluster identities in Seurat object
    seurat_T_cellBROn3.2 <- RenameIdents(seurat_T_cellBROn3.2, new_cluster_names)
    #save first order clustering
    seurat_T_cellBROn3.2$first_order_cl <- seurat_T_cellBROn3.2@active.ident

    # Update meta data to reflect new cluster names
    seurat_T_cellBROn3.2$seurat_clusters <- seurat_T_cellBROn3.2@active.ident
    
    # Set seed for reproducibility
    set.seed(123)
    
    # Create contingency table of clusters by treatment groups
    cell_type_distribution <- seurat_T_cellBROn3.2@meta.data %>%
      group_by(seurat_clusters, treatment) %>%
      summarise(cell_count = n())
    
    # Plot proportion of treatments within each cluster
    ggplot(cell_type_distribution, aes(x = seurat_clusters, y = cell_count, fill = treatment)) +
      geom_bar(stat = "identity", position = "fill") +  # stacked proportional bars
      labs(x = "Cluster", y = "Proportion", fill = "Treatment") +
      theme_minimal() +
      theme(aspect.ratio = 1, axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
      ggtitle("Proportion of Treatment in Each Cluster")
    
    # Alternate bar plot showing proportion of clusters within each treatment group
    ggplot(cell_type_distribution, aes(fill = seurat_clusters, y = cell_count, x = treatment)) +
      geom_bar(stat = "identity", position = "fill") +
      labs(x = "Treatment", y = "Proportion", fill = "Cluster") +
      theme_minimal() +
      theme(aspect.ratio = 1, axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
      ggtitle("Proportion of Clusters in Each Treatment")
    
    ### Highlight exhausted and non-exposed cells
    
    # Copy the Seurat object to work on highlighting
    seurat_T_cellBROn3.2_play <- seurat_T_cellBROn3.2
    
    # Create a new metadata column 'new_clusters' from the renamed clusters
    seurat_T_cellBROn3.2_play$second_order_cl <- as.character(seurat_T_cellBROn3.2_play$first_order_cl)
    
    # Reclassify 'exhausted_exp_nexp' cluster based on treatment status
    seurat_T_cellBROn3.2_play$second_order_cl <- ifelse(
      seurat_T_cellBROn3.2_play$second_order_cl == "Tex" & seurat_T_cellBROn3.2_play$treatment == "no exposure",
      "Tex-non-exposed",
      ifelse(
        seurat_T_cellBROn3.2_play$second_order_cl == "Tex" & seurat_T_cellBROn3.2_play$treatment != "no exposure",
        "Tex-DMGO-exposed",
        seurat_T_cellBROn3.2_play$second_order_cl
      )
    )
    
    # Plot UMAP colored by new cluster labels
    DimPlot(seurat_T_cellBROn3.2_play, reduction = "umap", group.by = "second_order_cl", pt.size = 1)
    
    # Define clusters to highlight
    clusters_to_highlight <- c("Tex-non-exposed", "Tex-DMGO-exposed")
    
    # Create a column for highlighted clusters, label others as "Other"
    seurat_T_cellBROn3.2_play$new_clusters_highlighted <- ifelse(
      seurat_T_cellBROn3.2_play$second_order_cl %in% clusters_to_highlight,
      seurat_T_cellBROn3.2_play$second_order_cl,
      "Other"
    )
    
    # Save a PDF with a UMAP plot highlighting the exhausted clusters
    pdf(paste0(res_path, "/Exhausted_cluster_zoom_out_points_new_colors.pdf"))
    
    DimPlot(
      seurat_T_cellBROn3.2_play,
      reduction = "umap",
      group.by = "new_clusters_highlighted",
      pt.size = 1,
      alpha = 0.5,
      cols = c(
        "Tex-non-exposed" = "seagreen",
        "Tex-DMGO-exposed" = "hotpink3",
        "Other" = "grey"
      )
    ) +
      theme(aspect.ratio = 1) +
      ggtitle("UMAP Plot Highlighting Exhausted Clusters")
    
    dev.off()
    
    # Update original object to include the new clusters and set identities accordingly
    seurat_T_cellBROn3.2 <- seurat_T_cellBROn3.2_play
    seurat_T_cellBROn3.2 <- SetIdent(seurat_T_cellBROn3.2, value = "second_order_cl")
     # Save the processed Seurat object for downstream analyses
    saveRDS(seurat_T_cellBROn3.2, paste0(data_loc, "T_Cells_BRO_n3_processed.rds"))
    
    ### Additional plots and data export
    
    library(readxl)
    library(ggplot2)
    library(writexl)
    
    # Update seurat_clusters in metadata for plotting
    seurat_T_cellBROn3.2$seurat_clusters <- seurat_T_cellBROn3.2$first_order_cl
    
    set.seed(123)
    pdf(paste0(res_path, "/Sankey_diagram.pdf"))
    
    # Create a contingency table grouped by clusters, original clusters and treatment (excluding unknown old clusters)
    cell_type_distribution <- seurat_T_cellBROn3.2@meta.data %>%
      ungroup() %>%
      group_by(seurat_clusters, old_clusters, treatment) %>%
      summarise(cell_count = n()) %>%
      filter(old_clusters != "unknown")
    
    # Plot stacked barplot of old clusters in each new cluster
    ggplot(cell_type_distribution, aes(x = seurat_clusters, y = cell_count, fill = old_clusters)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(x = "Cluster", y = "Cell Count", fill = "Old Cluster") +
      theme_minimal() +
      ggtitle("Proportion of Old Clusters in Each New Cluster") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))
    

    
    # Recalculate cell type distribution for treatments
    cell_type_distribution <- seurat_T_cellBROn3.2@meta.data %>%
      group_by(seurat_clusters, treatment) %>%
      summarise(cell_count = n())
    
    # Save proportions data to Excel
    cell_type_distribution2 <- cell_type_distribution %>%
      group_by(treatment) %>%
      mutate(total_n = sum(cell_count)) %>%
      ungroup() %>%
      mutate(proportion = cell_count / total_n)%>%filter(treatment%in%c("Bulk", "no exposure"))
    
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(ggsankey)  # for Sankey plots
    library(tidyverse) # for general data manipulation
    
    # Prepare full combinations of clusters and treatments to fill missing data with zeros
    full_combinations <- expand.grid(
      seurat_clusters = unique(cell_type_distribution2$seurat_clusters),
      treatment = unique(cell_type_distribution2$treatment)
    )
    
    # Join full combinations with the original dataset and replace missing proportions with 0
    cell_type_distribution2_full <- full_combinations %>%
      left_join(cell_type_distribution2, by = c("seurat_clusters", "treatment")) %>%
      mutate(proportion = ifelse(is.na(proportion), 0, proportion))
    
    # Ensure 'treatment' factor has consistent ordering
    cell_type_distribution2_full$treatment <- factor(cell_type_distribution2_full$treatment, levels = c("no exposure", "Bulk"))
    
    # Sankey plot showing the proportion flow between treatments and clusters
    ggplot(cell_type_distribution2_full, aes(
      x = treatment,
      node = seurat_clusters,
      fill = seurat_clusters,
      value = proportion)) +
      geom_sankey_bump(space = 0, type = "alluvial", color = "transparent", smooth = 15) +
      theme_sankey_bump(base_size = 16) +
      labs(x = NULL, y = "Proportion", fill = "Cluster") +
      theme(aspect.ratio = 1, legend.position = "right") +
      ggtitle("Sankey Plot of Cluster Proportions by Treatment")
    
    # Visualize clustering and treatment with UMAP plots
    DimPlot(seurat_T_cellBROn3.2, reduction = "umap", group.by = "sequencing_batch", pt.size = 0.5) + theme(aspect.ratio = 1)
    DimPlot(seurat_T_cellBROn3.2, reduction = "umap", group.by = "treatment", pt.size = 0.5) + theme(aspect.ratio = 1)
    DimPlot(seurat_T_cellBROn3.2, reduction = "umap") + theme(aspect.ratio = 1)
    
    dev.off()  # Close any open graphic device
    
    # Identify marker genes for each cluster using Seurat's FindAllMarkers
    markers <- FindAllMarkers(seurat_T_cellBROn3.2, min.pct = 0.25, logfc.threshold = 1)
    
    library(writexl)
    # Save marker genes to Excel
    write_xlsx(markers, path = paste0(res_path, "N3_T_cell_DEG", Sys.Date(), ".xlsx"))
    

    # Define gene sets for exhaustion and other signatures to visualize
    ex_markers <- c("TNFRSF9", "TNFRSF18", "LAYN", "CD83", "CTLA4", "DUSP4", "GNLY", "ENTPD1", "IGFLR1", "NR4A2")
    ex_tf <- c("SOX4", "RGS16", "RGS1", "SNX9")
    
    # Create DotPlots or heatmaps for these markers in specified clusters
    pdf(paste0(res_path, "/Dotplot exhausted cluster markers", Sys.Date(), ".pdf"))
    
    # Set factor levels to order clusters in plots
    seurat_T_cellBROn3.2@active.ident <- factor(
      seurat_T_cellBROn3.2@active.ident,
      levels = c(
        "mid_term_activated", "metabolically_stressed", "unresponsive",
        "exhaust_non_exp", "exhaust_exp", "IL2_responsive", "IFN-I responsive",
        "non-exposed", "prolif", "sustained_cytotoxic", "migrat_interact", "NCAM-"
      )
    )
    
    generate_heatmap(seurat_T_cellBROn3.2, idents = c("exhaust_non_exp", "exhaust_exp"), features = ex_markers)
    generate_heatmap(seurat_T_cellBROn3.2, idents = c("exhaust_non_exp", "exhaust_exp"), features = ex_tf)
    
    dev.off()
    
    # Load external dataset from Wang et al. for comparison
    seurat_Wang <- readRDS(paste0(data_loc, "Chu_CD8.rds"))
    seurat_Wang <- SetIdent(seurat_Wang, value = "cell.type")
    
    DimPlot(seurat_Wang)
    
    # Load DEG results for Ths cluster to project onto Wang dataset
    NCAM1_neg_sig <- subset(markers, cluster == "Ths" & avg_log2FC > 0 & p_val_adj < 0.05)
    
    # Add module scores for NCAM1 negative signature
    seurat_Wang <- AddModuleScore(seurat_Wang, features = list(NCAM1_neg_sig$gene), name = "NCAM1_neg_n3")
    
    pdf(paste0(res_path, "/NCAM1 neg projection on Wang.pdf"), width = 10, height = 15)
    DimPlot(seurat_Wang) + theme(aspect.ratio = 1)
    FeaturePlot(seurat_Wang, features = "NCAM1_neg_n31", cols = c("deepskyblue", "coral3", "coral")) + theme(aspect.ratio = 1)
    dev.off()
    
    # Reload the saved Seurat object for further analysis
    seurat_T_cellBROn3.2 <- readRDS("D:/SURF_2/Shared/Dream3DLab (Groupfolder)/1.Projects/BRO_BrainOrganoids/3.Analysis/Results_Maria/revision20241204_all_1000features/seurat_T_cellBROn3_1000_features.rds")
    
    # Subset to exposed populations with Bulk treatment only
    seurat_T_cellBROn3.2_sub <- subset(
      seurat_T_cellBROn3.2,
      idents = c(
        "IFN-I responsive", "prolif", "metabolically_stressed", "unresponsive",
        "non-exposed", "exhaust_exp", "IL2_responsive", "sustained_cytotoxic",
        "mid_term_activated", "migrat_interact"
      )
    )
    seurat_T_cellBROn3.2_sub <- subset(seurat_T_cellBROn3.2_sub, treatment %in% c("Bulk"))
    
    # Load Chu et al signatures from Excel
    Chu_sig <- read_excel(paste0(data_loc, "41591_2023_2371_MOESM3_ESM_Chu.xlsx"), sheet = "Table S4", skip = 1)
    
    # Add module scores dynamically for each gene list column in Chu_sig
    for (col_name in colnames(Chu_sig)) {
      gene_list <- list(Chu_sig[[col_name]])
      seurat_T_cellBROn3.2_sub <- AddModuleScore(
        object = seurat_T_cellBROn3.2_sub,
        features = gene_list,
        name = col_name
      )
    }
    
    # Collect module score names generated by AddModuleScore (AddModuleScore appends '1' by default)
    module_names <- paste0(colnames(Chu_sig), "1")
    
    # Extract DotPlot data for module scores
    dot_data <- DotPlot(seurat_T_cellBROn3.2_sub, features = module_names)$data
    
    # Order clusters for plotting
    dot_data$id <- factor(dot_data$id, levels = c(
      "non-exposed", "exhaust_exp", "IL2_responsive", "mid_term_activated", "sustained_cytotoxic",
      "migrat_interact", "IFN-I responsive", "prolif", "metabolically_stressed", "unresponsive"
    ))
    
    # Plot heatmap-style tile plot for Chu signatures across clusters
    pdf(paste0(res_path, "Additional Chu signatures projection", Sys.Date(), ".pdf"), height = 15, width = 10)
    ggplot(dot_data, aes(y = id, x = features.plot, fill = avg.exp.scaled)) +
      geom_tile() +
      scale_fill_gradient2(low = "deepskyblue", mid = "white", high = "coral2", midpoint = 0.5) +
      coord_flip() +
      ggtitle("n3 Chu sig projection Normalized by Column") +
      theme(
        aspect.ratio = 1,
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 5)
      )
    dev.off()
    