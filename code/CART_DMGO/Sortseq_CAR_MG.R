set.seed(123)

data_loc <- "D:/SURF_2/Shared/Dream3DLab (Groupfolder)/1.Projects/BRO_BrainOrganoids/3.Analysis/Github repo/Data/"
res_path <- "D:/SURF_2/Shared/Dream3DLab (Groupfolder)/1.Projects/BRO_BrainOrganoids/3.Analysis/Github repo/Repeat_output/"
# Check if the directory exists, if not, create it
if (!dir.exists(res_path)) {
  dir.create(res_path, recursive = TRUE)  # Use recursive = TRUE to create all necessary parent directories
}

## function for enirchment analysis


library(enrichR)
library(openxlsx)
library(ggplot2)  # For plotting enrichment results
# Subset the data for the identities of interest
perform_ttest <- function(seurat_object, ident_values, features, test_type = "anova") {
  # Extract relevant data from the Seurat object
  data <- FetchData(seurat_object, vars = c("ident", features))
  
  # Subset the data for the selected 'ident' values
  data <- data[data$ident %in% ident_values, ]
  
  if (test_type == "anova") {
    # Perform ANOVA and Tukey's HSD test
    results <- lapply(features, function(feature) {
      # Perform ANOVA for the given feature
      anova_model <- aov(data[[feature]] ~ data$ident)
      
      # Perform Tukey's HSD test on the ANOVA model
      tukey_result <- TukeyHSD(anova_model)
      
      # Extract group comparisons, p-values, and test statistics
      tukey_data <- data.frame(tukey_result$`data$ident`)  # Adjust variable name here
      tukey_data$Comparison <- rownames(tukey_data)
      
      # Prepare a data frame with relevant results for the feature
      cbind(Feature = feature, tukey_data)
    })
    
    # Combine the results from all features into a single data frame
    test_results <- do.call(rbind, results)
  } else {
    # Perform t-tests
    results <- lapply(features, function(feature) {
      # Perform t-test for the given feature
      t_test <- t.test(data[[feature]] ~ data$ident)
      
      # Return a data frame with the feature name, p-value, and test statistic
      data.frame(
        Feature = feature,
        Statistic = t_test$statistic,
        PValue = t_test$p.value
      )
    })
    
    # Combine the results from all features into a single data frame
    test_results <- do.call(rbind, results)
    
    # Adjust p-values for multiple testing using the Benjamini-Hochberg method
    test_results$AdjustedPValue <- p.adjust(test_results$PValue, method = "BH")
  }
  
  # Return the test results
  return(test_results)
}


generate_heatmap <- function(seurat_object, idents, features, num_cells_per_ident = 164,to_cl=FALSE, color_palette = c("dodgerblue4","dodgerblue3","dodgerblue", "white", "red","red3","red4")) {
  library(dplyr)
  library(pheatmap)
  library(tibble)
  
  # Extract cells belonging to the specified clusters
  cells <- WhichCells(seurat_object, idents = idents)
  
  # Retrieve the cluster identities for the selected cells
  cluster_annotations <- as.data.frame(Idents(seurat_object)[cells])
  colnames(cluster_annotations) <- "Cluster"  # Rename column for clarity
  features<-features[features%in%rownames(seurat_object)]
  # If num_cells_per_ident is specified, sample cells
  if (!is.null(num_cells_per_ident)) {
    cluster_annotations <- cluster_annotations %>%
      rownames_to_column("Cell") %>%
      group_by(Cluster) %>%
      sample_n(num_cells_per_ident, replace = FALSE) %>%  # Sample cells without replacement
      column_to_rownames("Cell")
  }
  
  # Order the data frame by the Cluster column
  cluster_annotations <- cluster_annotations[order(cluster_annotations$Cluster), , drop = FALSE]
  
  # Sort cells by the new order of clusters
  sorted_cells <- rownames(cluster_annotations)  # Use rownames, which are cell names, now ordered
  
  # Extract normalized expression data for sorted cells
  expression_data <- GetAssayData(
    object = seurat_object, 
    slot = "data", # Use normalized data
    assay = "RNA"
  )[features, sorted_cells]  # Retrieve expression data for the selected markers and sorted cells
  expression_matrix <- as.matrix(expression_data)

  # Step 2: Scaling each row (gene) independently to [-1, 1]
  # Scale each gene (row) to match DotPlot behavior
  scaled_data <- t(apply(expression_matrix, 1, function(x) {
    if (sd(x) == 0) {
      return(rep(0, length(x)))  # Handle rows with zero variance
    }
    (x - mean(x)) / sd(x)  # Z-score normalization
  }))
  
  # Prepare annotation for heatmap columns (cells)
  annotation_col <- cluster_annotations  # Use the ordered cluster annotations for the cells
  
  
  # Generate the heatmap
  pheatmap(
    expression_matrix,scale="row",
    cluster_rows = FALSE,       # Do not cluster genes (rows)
    cluster_cols = to_cl,       # Do not cluster cells (columns)
    show_colnames = FALSE,      # Hide column names (cells)
    treeheight_row = 0,         # No dendrogram for rows
    cellwidth = 0.25, 
    cellheight = 10, 
    annotation_col = annotation_col,  # Add cluster annotations to columns
    color = colorRampPalette(color_palette)(50),  # Custom color scale,
    treeheight_col = F
  )
}


#Create a Seurat object with this scRNA seq data:
library(Seurat)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(ggtrace)
rds_path<-paste0(res_path,"seurat_MG.rds")
if (file.exists(rds_path)) {
  seurat_MG<-readRDS(rds_path)
}else{
  seurat_T_cellmicroglia<-readRDS(paste0(data_loc, "/Sort_seq_unprocessed.rds"))
  ##subset the wells that passed the QC:
  seurat_T_cellmicroglia<-subset(seurat_T_cellmicroglia, ERCC_well=="pass")
  # Define the mapping from the Library to the cells_in_plate labels
  plate_mapping <- c(
    "PMC-AW-s006" = "PMP d21 from PO",
    "PMC-AW-s007" = "PMP d21 from DMGO",
    "PMC-AW-s008" = "CAR from DMGO",
    "PMC-AW-s010" = "CAR from DMGO + PMP",
    "PMC-AW-s012" = "Unexposed CAR",
    "PMC-AW-s014" = "Unexposed PMP",
    "PMC-AW-s016" = "NCAM+ CAR from DMGO",
    "PMC-AW-s018" = "NCAM+ CAR from DMGO + PMP",
    "PMC-AW-s021" = "PMP d21 from PO",
    "PMC-AW-s022" = "PMP d21 from DMGO",
    "PMC-AW-s023" = "CAR from DMGO",
    "PMC-AW-s025" = "CAR from DMGO + PMP",
    "PMC-AW-s026" = "PMP from DMGO + CAR",
    "PMC-AW-s027" = "NCAM+ CAR from DMGO",
    "PMC-AW-s029" = "NCAM+ CAR from DMGO + PMP",
    "PMC-AW-s030" = "PMP from DMGO + NCAM+ CAR",
    "PMC-AW-s031" = "Unexposed CAR",
    "PMC-AW-s033" = "Unexposed PMP"
  )
  
  # Create the new column based on the "Library" metadata
  seurat_T_cellmicroglia@meta.data$cells_in_plate <- plate_mapping[seurat_T_cellmicroglia@meta.data$Library]
  
  DimPlot(seurat_T_cellmicroglia , group.by = "Library", split.by = "cells_in_plate")
  

  # Create the new column based on the "Library" metadata
  seurat_T_cellmicroglia@meta.data$cells_in_plate <- plate_mapping[seurat_T_cellmicroglia@meta.data$Library]
  
  seurat_T_cellmicroglia<-FindClusters(seurat_T_cellmicroglia, resolution = 0.15)
  DimPlot(seurat_T_cellmicroglia , group.by ="cells_in_plate")
  
  DimPlot(seurat_T_cellmicroglia )
  
  
  ##split the CAR and the Microglia data, we split based on plate and on cluster, since there are some cells that are contaminations from the organoid
  seurat_MG<-subset(seurat_T_cellmicroglia, cells_in_plate%in%c("Unexposed PMP",
                                                                "PMP d21 from PO",
                                                                "PMP d21 from DMGO",
                                                                "PMP from DMGO + CAR",
                                                                "PMP from DMGO + NCAM+ CAR")&seurat_clusters %in%c("6","5", "1"))
  
  seurat_CAR<-subset(seurat_T_cellmicroglia, cells_in_plate%in%c("CAR from DMGO",
                                                                 "CAR from DMGO + PMP",
                                                                 "Unexposed CAR",
                                                                 "NCAM+ CAR from DMGO",
                                                                 "NCAM+ CAR from DMGO + PMP")&seurat_clusters %in%c("2","0", "4","3")
  )
  
  
  
  DimPlot(seurat_MG )
  DimPlot(seurat_CAR)
  
  saveRDS(seurat_CAR, paste0(res_path,"seurat_CAR.rds"))
  saveRDS(seurat_MG , paste0(res_path,"seurat_MG.rds"))
  
}


###Analysis exposed non exposed
seurat_MG <- SetIdent(object = seurat_MG, value = "cells_in_plate")
DimPlot(seurat_MG )
Idents(seurat_MG) <- factor(Idents(seurat_MG), levels = c("Unexposed PMP", "PMP d21 from PO", "PMP d21 from DMGO","PMP from DMGO + CAR","PMP from DMGO + NCAM+ CAR"))


##Plot markers of interest
##plot only PMP without CAR
seurat_MG_noCAR<-subset(seurat_MG,cells_in_plate%in%c("Unexposed PMP",
                                                      "PMP d21 from PO",
                                                      "PMP d21 from DMGO"))
# Heatmap 1 genes
heatmap1_genes <- c("SPI1", "CTCF", "IRF1", "IRF2", "IRF3", "IRF8", "IRF9", "RUNX1", "RUNX2", "JUN", "JUNB", "JUND", 
                    "FOS", "FOSB", "FOSL2", "ATF4", "BATF", "BATF2", "BATF3", "CEBPA", "CEBPB", "CEBPG", 
                    "MEF2A", "MEF2B", "MEF2C", "MEF2D", "SMAD3", "MAF", "MAF1", "MAFB", "MAFF", "MAFG", "MAFK")


##select only the genes that are upregulated in PMP in PO and or PMP in DMGO vs unexposed PMP


markers_PMP_exposure<-FindMarkers(seurat_MG,ident.1=c("PMP d21 from PO"),      
                                  ident.2="Unexposed PMP" )

markers_PMP_exposure$cluster<-"PMP d21 from PO"
markers_PMP_exposure$gene<-rownames(markers_PMP_exposure)
# Write to Excel file
#write_xlsx(markers_PMP_exposure,paste0(res_path,"DEG_PMP exposure vs non exposure", Sys.Date(), ".xlsx"))
markers_PMP_exposure_up<-subset(markers_PMP_exposure, p_val_adj<0.05&avg_log2FC>0)

heatmap1_genes<-heatmap1_genes[heatmap1_genes%in%markers_PMP_exposure_up$gene]

pdf(paste0(res_path, "/PMP_BRO_DMGO_plot1", Sys.Date(),".pdf"), width=15, height = 15)
generate_heatmap(seurat_MG_noCAR, idents= c("Unexposed PMP", "PMP d21 from PO"), features=microglia_genes_Nat )

dev.off()


##Plot markers of mouse microglia from science paper: this Science paper https://www.science.org/doi/full/10.1126/science.aad8670 on mouse microglia 
library(readxl)
micrglia_Science <- read_excel(paste0(data_loc,"/aad8670-matcovitch-natan-sm.tables-s1-to-s10_Science.xlsx"), 
                               sheet = "Supp. Table 1", skip = 1)

micrglia_Science <- micrglia_Science %>%
  mutate(Gene = toupper(Gene))

MG_embryo<-subset(micrglia_Science , Cluster%in%c(1,2,3))$Gene
MG_birth<-subset(micrglia_Science , Cluster%in%c(4,5))$Gene
MG_adult<-subset(micrglia_Science , Cluster%in%c(6,7))$Gene

seurat_MG_noCAR <- AddModuleScore(
  object = seurat_MG_noCAR,
  features = list(MG_embryo),  # Genes for the current state
  name =  "MG_embryo")

seurat_MG_noCAR <- AddModuleScore(
  object = seurat_MG_noCAR,
  features = list(MG_birth),  # Genes for the current state
  name =  "MG_birth")

seurat_MG_noCAR <- AddModuleScore(
  object = seurat_MG_noCAR,
  features = list(MG_adult),  # Genes for the current state
  name =  "MG_adult")


pdf(paste0(res_path, "/PMP_BRO_DMGO_diff_combinations", Sys.Date(),".pdf"), width=15, height = 15)

data<-DotPlot(seurat_MG_noCAR, features = c("MG_adult1","MG_birth1","MG_embryo1"),scale=T, cols=c("deepskyblue", "coral2")) +coord_flip()+
  ggtitle(paste("Signatures Science paper 2016")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
df<-data[["data"]]

df$id <- factor(df$id, 
                levels = c( "PMP d21 from DMGO", "PMP d21 from PO" , "Unexposed PMP"))  # Replace with your actual order
ggplot(df, aes(x = features.plot, y =id, fill = avg.exp.scaled)) + 
  geom_tile(color="grey2") + 
  scale_fill_gradient2(low = "deepskyblue", high = "coral2", mid = "white", midpoint = 0.0) + 
  theme_minimal() + 
  coord_flip() + 
  ggtitle("Signatures Science paper 2016 not scaled") + 
  theme(asepct.ratio=1,axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.text.y = element_text(size = 6)) +  # Adjust y-axis labels if needed
  scale_y_discrete(limits = rev(levels(df$id)))  # Order clusters by the factor levels

dev.off()




## Adnrade et al markers: https://www.nature.com/articles/s41467-024-52096-w#MOESM1
Andrade_Nat_Com_2024 <- as.data.frame(read_excel(paste0(data_loc,"/Andrade_Nat_Com_2024_41467_2024_52096_MOESM9_ESM.xlsx"), 
                                                 sheet = "S. Fig2C"))

Andrade_Nat_Com_2024 <- Andrade_Nat_Com_2024%>%  ###select the top overexpressed markers
  filter(avg_log2FC > 1 | avg_log2FC < -1)

unique_clusters<-unique(Andrade_Nat_Com_2024$cluster)

# Ensure your upregulated and downregulated genes are in a named list
genes_per_cluster_up <- subset(Andrade_Nat_Com_2024 ,avg_log2FC > 0)
genes_per_cluster_down <- subset(Andrade_Nat_Com_2024 ,avg_log2FC < 0)

# Adding combined module scores for each cluster
for (clusters in unique_clusters) {
  # Get the upregulated and downregulated genes for this cluster
  upregulated_genes <- subset(genes_per_cluster_up, cluster==clusters)$gene
  downregulated_genes <- subset(genes_per_cluster_down, cluster==clusters)$gene
  
  # Compute scores separately for upregulated and downregulated genes
  # Check if there are any upregulated genes
  if (length(upregulated_genes) > 0) {
    # Compute module score for upregulated genes
    seurat_MG_noCAR <- AddModuleScore(
      object = seurat_MG_noCAR,
      features = list(upregulated_genes),
      name = paste0(clusters, "_Upregulated")
    )
  } else {
    # Assign zeros for the upregulated module score
    seurat_MG_noCAR@meta.data[[paste0(cluster, "_Upregulated1")]] <- 0
  }
  
  # Check if there are any downregulated genes
  if (length(downregulated_genes) > 0) {
    # Compute module score for downregulated genes
    seurat_MG_noCAR <- AddModuleScore(
      object = seurat_MG_noCAR,
      features = list(downregulated_genes),
      name = paste0(clusters, "_Downregulated")
    )
  } else {
    # Assign zeros for the downregulated module score
    seurat_MG_noCAR@meta.data[[paste0(clusters, "_Downregulated1")]] <- 0
  }
  
  # Combine the scores by subtracting downregulated from upregulated
  combined_score_name <- paste0(clusters, "_CombinedScore")
  seurat_MG_noCAR@meta.data[[combined_score_name]] <-
    seurat_MG_noCAR@meta.data[[paste0(clusters, "_Upregulated1")]] -
    seurat_MG_noCAR@meta.data[[paste0(clusters, "_Downregulated1")]]
  
  # Remove intermediate scores (optional)
  seurat_MG_noCAR@meta.data[[paste0(clusters, "_Upregulated1")]] <- NULL
  seurat_MG_noCAR@meta.data[[paste0(clusters, "_Downregulated1")]] <- NULL
}

# Inspect the final metadata to ensure combined scores are present
##plot only PMP without CAR

andrade_sig<-paste0(unique_clusters, "_CombinedScore")

pdf(paste0(res_path, "/Mg markers from Andrade Paper", Sys.Date(),".pdf"), width=15, height = 15)


VlnPlot(seurat_MG_noCAR, idents = c("Unexposed PMP","PMP d21 from PO"),
        features=c("Mg_TAM_CombinedScore","Mo_TAM_CombinedScore"), flip = T,stack = TRUE, fill.by = "feature",same.y.lims = T, alpha = 0.5)+
  scale_y_continuous(limits = c(-0.3, 0.3)) +  # Apply consistent y-axis limits
  theme(aspect.ratio = 0.5)+
  stat_summary(
    fun = median,  # Function to calculate the median
    geom = "point",  # Use points to display the medians
    color = "yellow",  # Color of the points
    size = 2  # Size of the points
  )+  geom_jitter(
    aes(fill = "grey"),  # Adjust color by group if needed
    size = 0.5,          # Size of points
    width = 0.1,         # Slight jitter to avoid overlap
    height = 0           # Avoid vertical jitter
  )+
  geom_hline(
    yintercept = 0,  # Add horizontal line at y = 0
    linetype = "dashed",  # Dashed line for distinction
    color = "gray40",  # Gray color for the line
    size = 0.5  # Thickness of the line
  )


# Example of usage:
t_test_results_unexp_PO <- perform_ttest(seurat_MG_noCAR, ident_values=c("Unexposed PMP", "PMP d21 from PO"), features=c("Mg_TAM_CombinedScore", "Mo_TAM_CombinedScore"))

write.csv(t_test_results_unexp_PO, file = paste0(res_path, "/t_test_results_unexp_PO", Sys.Date(),".csv"))

VlnPlot(seurat_MG_noCAR, idents = c("PMP d21 from PO","PMP d21 from DMGO"),
        features=c("Mg_TAM_CombinedScore","IFN_Mg_TAM_CombinedScore" ,        
                   "Phago_Lipid_Mg_TAM_CombinedScore", "Hypoxic_Mg_TAM_CombinedScore" ), flip = T,stack = TRUE, fill.by = "feature",group.by=,same.y.lims = T, alpha = 0.5)+
  scale_y_continuous(limits = c(-0.4, 0.6)) +  # Apply consistent y-axis limits
  theme(aspect.ratio = 0.5)+
  stat_summary(
    fun = median,  # Function to calculate the median
    geom = "point",  # Use points to display the medians
    color = "yellow",  # Color of the points
    size = 2  # Size of the points
  )+geom_jitter(
    aes(fill = "grey"),  # Adjust color by group if needed
    size = 0.5,          # Size of points
    width = 0.1,         # Slight jitter to avoid overlap
    height = 0           # Avoid vertical jitter
  )+
  geom_hline(
    yintercept = 0,  # Add horizontal line at y = 0
    linetype = "dashed",  # Dashed line for distinction
    color = "gray40",  # Gray color for the line
    size = 0.5  # Thickness of the line
  )


t_test_results_PO_DMGO <- perform_ttest(seurat_MG_noCAR, ident_values=c("PMP d21 from PO","PMP d21 from DMGO"), features=c("Mg_TAM_CombinedScore","IFN_Mg_TAM_CombinedScore" ,        
                                                                                                                           "Phago_Lipid_Mg_TAM_CombinedScore", "Hypoxic_Mg_TAM_CombinedScore" ))
write.csv(t_test_results_PO_DMGO, file = paste0(res_path, "/t_test_results_PO_DMGO", Sys.Date(),".csv"))




immune_genes_to_plot <- rev(c(
  "CCL2", "CCL20", "CCL3", "CCL3L1", "CCL3L3", "CCL4", 
  "CCL4L1", "CCL4L2", "CXCL2", "CXCL3", "IL8", "PLAUR", 
  "CD81", "CX3CR1", "FCGBP", "HAVCR2", "KLF2", "LGALS9", 
  "S100A9", "SIRPA", "TGFB1"
))

DotPlot(seurat_MG_noCAR, idents= c("PMP d21 from DMGO","PMP d21 from PO"), features =immune_genes_to_plot , cols=c("deepskyblue", "coral2"))+coord_flip()+ theme(aspect.ratio=1,axis.text.x = element_text(angle = 90, hjust = 1)) +  # Optionally adjust the text angle if needed
  ggtitle("Microglia /Immune suppressive Genes from Andrade paper all populations")  # Add title for clarity




dev.off()





##Matche the populations of SORTseq CART with and without MG to the n3 CART dataset where we identified different t cell types
seurat_T_cellBROn3.2<-readRDS(paste0(data_loc, "T_Cells_BRO_n3_processed.rds"))

# Step 2: Find transfer anchors between the two datasets
anchors <- FindTransferAnchors(
  query = seurat_CAR,         # Reference dataset (source of labels)
  reference  = seurat_T_cellBROn3.2,  # Query dataset (where labels are transferred to)
  dims = 1:30                    # Use the appropriate dimensions for the analysis
)

# Step 3: Transfer the `Idents` or other metadata from reference to query
seurat_CAR<- TransferData(
  anchorset = anchors,
  query = seurat_CAR,
  reference = seurat_T_cellBROn3.2,
  refdata = Idents(seurat_T_cellBROn3.2),  # Transfer cluster labels
  dims = 1:30
)

# Step 4: Set the transferred identities as the new `Idents` in the query object
Idents(seurat_CAR) <- seurat_CAR$predicted.id
seurat_CAR$seurat_clusters<-seurat_CAR@active.ident

cell_type_distribution <- seurat_CAR@meta.data %>%
  group_by(cells_in_plate) %>%
  group_by(seurat_clusters, cells_in_plate) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  # Calculate total cell count per treatment
  group_by(cells_in_plate) %>%
  mutate(total_cells = sum(cell_count)) %>%
  # Calculate percentage for each cluster within treatment
  mutate(percentage = (cell_count / total_cells) * 100) %>%
  ungroup()  # Optional, to remove grouping after calculations

# Plot the barplot using ggplot
pdf(paste0(res_path, "Matching of Sortseq CAR to n3 dataset", Sys.Date(),".pdf"), height = 15, width = 10)

ggplot(cell_type_distribution, aes(fill = seurat_clusters, y = percentage, x = cells_in_plate)) +
  geom_bar(stat = "identity", position = "stack") +  # 'fill' makes it proportional
  labs(x = "Cluster", y = "Proportion", fill = "Cell Type") +
  theme_minimal() +theme(aspect.ratio=1,axis.text.x = element_text(angle = 90 ,hjust =0, vjust=0.5))+
  ggtitle("Proportion of treatment in Each Cluster")#+facet_grid(sequencing_batch~., scales = "free_y")

dev.off()

saveRDS(seurat_CAR,paste0(data_loc,"Sortseq_CAR.rds") )
saveRDS(seurat_MG_noCAR,paste0(data_loc,"Sortseq_MG.rds") )






