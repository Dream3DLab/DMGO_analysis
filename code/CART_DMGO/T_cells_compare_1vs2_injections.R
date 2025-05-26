# Set seed for reproducibility
set.seed(123)

# Define paths
data_loc <- "D:/SURF_2/Shared/Dream3DLab (Groupfolder)/1.Projects/BRO_BrainOrganoids/3.Analysis/Github repo/Data/"
res_path <- "D:/SURF_2/Shared/Dream3DLab (Groupfolder)/1.Projects/BRO_BrainOrganoids/3.Analysis/Github repo/Repeat_output/"

# Create results directory if it doesn't exist
if (!dir.exists(res_path)) {
  dir.create(res_path, recursive = TRUE)
}

# Load libraries
library(openxlsx)
library(ggplot2)
library(Seurat)
library(dplyr)

# Load Seurat objects
seurat_n3 <- readRDS(paste0(data_loc, "T_Cells_BRO_n3_processed.rds"))
seurat_sortseq <- readRDS(paste0(data_loc, "Sortseq_CART_predicted_n3.rds"))

# Assign clusters to metadata
seurat_n3$seurat_clusters <- seurat_n3@active.ident

# Add early/late treatment label
meta_n3 <- seurat_n3@meta.data
meta_n3$treatment <- ifelse(meta_n3$Sample %in% c("X8_v30_IPSC_CAR"),
                            paste(meta_n3$treatment, "Late"),
                            paste(meta_n3$treatment, "Early"))

# Calculate cell type proportions for BRO n3 data
distribution_n3 <- meta_n3 %>%
  group_by(seurat_clusters, sequencing_batch, treatment) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(treatment, sequencing_batch) %>%
  mutate(total_cells = sum(cell_count),
         percentage = (cell_count / total_cells) * 100) %>%
  ungroup() %>%
  mutate(unique_exp = interaction(treatment, sequencing_batch))

# Calculate cell type proportions for sort-seq data
distribution_sortseq <- seurat_sortseq@meta.data %>%
  group_by(seurat_clusters, cells_in_plate, Library) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(cells_in_plate, Library) %>%
  mutate(total_cells = sum(cell_count),
         percentage = (cell_count / total_cells) * 100) %>%
  ungroup() %>%
  mutate(unique_exp = interaction(cells_in_plate, Library))

# Filter sortseq data to include only relevant CAR condition
distribution_sortseq <- subset(distribution_sortseq, cells_in_plate == "CAR from DMGO")

# Combine both datasets
combined_distribution <- rbind(
  distribution_n3[, c("seurat_clusters", "cell_count", "total_cells", "percentage", "unique_exp")],
  distribution_sortseq[, c("seurat_clusters", "cell_count", "total_cells", "percentage", "unique_exp")]
)

# Subset to selected experimental conditions
selected_samples <- c(
  "Bulk Early.1", "Bulk Late.1", "Bulk Early.3",
  "CAR from DMGO.PMC-AW-s008", "CAR from DMGO.PMC-AW-s023"
)

filtered_distribution <- subset(combined_distribution, unique_exp %in% selected_samples)

# Add grouped labels
filtered_distribution$groups <- case_when(
  filtered_distribution$unique_exp == "Bulk Late.1" ~ "Late",
  filtered_distribution$unique_exp %in% c("Bulk Early.1", "Bulk Early.3") ~ "Early_bulk",
  TRUE ~ "Early_sortseq"
)

# Normalize percentages per group
normalized_distribution <- filtered_distribution %>%
  group_by(groups, seurat_clusters) %>%
  summarise(percentage = mean(percentage), .groups = "drop") %>%
  group_by(groups) %>%
  mutate(total = sum(percentage),
         norm_percentage = (percentage / total) * 100)

# Plotting
pdf(file = paste0(res_path, "/Split_seq_exp_early_late_", Sys.Date(), ".pdf"))

ggplot(normalized_distribution, aes(x = groups, y = norm_percentage, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Proportion of Treatment in Each Cluster",
       x = "Experimental Group",
       y = "Normalized Proportion (%)",
       fill = "Cell Type") +
  theme_minimal() +
  theme(aspect.ratio = 1,
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))

dev.off()