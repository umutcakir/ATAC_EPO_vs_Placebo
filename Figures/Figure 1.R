library(Seurat)
library(Signac)
library(ggplot2)
library(scplotter)
library(tidyr)
set.seed(0)

ident_order <- c("A1", "A2", "A4", "A6", "B7", "B8", "B9", "B11")

# General function to save ggplot objects in high resolution
save_highres_plot <- function(plot, filename, width = 12, height = 10, dpi = 1200, format = "png") {
  # Determine file extension
  ext <- tolower(format)
  
  # Ensure filename ends with the correct extension
  if (!grepl(paste0("\\.", ext, "$"), filename)) {
    filename <- paste0(filename, ".", ext)
  }
  
  # Choose appropriate device based on format
  ggsave(
    filename = filename,
    plot = plot,
    device = ext,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    bg = "white"  # white background avoids transparency issues
  )
}

palette_umap <- c(
  "Oligodendrocytes" = "maroon",      
  "Pyramidal neurons" = "darkcyan",      
  "Interneurons" = "gold",     
  "Intermediate cells" = "navy",    
  "Dentate gyrus neurons" = "turquoise",   
  "Astrocytes" = "blue",   
  "Microglia" = "darkorange",  
  "Endothelial cells" = "deepskyblue",     
  "Pericytes" = "purple",    
  "Ependymal cells" = "grey",     
  "Neuroimmune cells" = "violet"
)



# 1B-1C Start
ATAC <- readRDS("../ATAC_analysis_with_combined_peaks/ATAC_with_sharedUMAP_and_imputedRNA.rds")
# Create the UMAP cluster plot with labels
p <- CellDimPlot(
  ATAC,
  group_by = "seurat_clusters",
  reduction = "umap",
  theme = "theme_blank",
  label = TRUE,
  label_insitu = TRUE,
  label_repel = TRUE
)

# Tabulate and sort predicted.id counts
cell_order <- names(sort(table(ATAC@meta.data$predicted.id), decreasing = TRUE))

# Save the plot using your general function
save_highres_plot(p, "Figure1/Figure1B_seurat_clusters_label_TRUE", format = "png")

# Tabulate and sort predicted.id counts
cell_order <- names(sort(table(ATAC@meta.data$predicted.id), decreasing = TRUE))

# Reassign predicted.id as a factor with the new order
ATAC@meta.data$predicted.id <- factor(ATAC@meta.data$predicted.id, levels = cell_order)

p <- CellDimPlot(
  ATAC,
  group_by = "predicted.id",
  reduction = "umap",
  theme = "theme_blank",
  palcolor = palette_umap,
  label = FALSE,
  label_insitu = FALSE,
  label_repel = FALSE
)

# Save the plot using your general function
save_highres_plot(p, "Figure1/Figure1E_11_clusters", format = "png")

p <- CellDimPlot(
  ATAC,
  group_by = "predicted.id",
  reduction = "umap",
  theme = "theme_blank",
  palcolor = palette_umap,
  label = TRUE,
  label_insitu = TRUE,
  label_repel = TRUE
)

# Save the plot using your general function
save_highres_plot(p, "Figure1/Figure1E_11_clusters_with_labels", format = "png")
# 1B-C Complete 

#S1C Start
ATAC <- readRDS("../ATAC_analysis_with_combined_peaks/ATAC_with_sharedUMAP_and_imputedRNA.rds")
library(ggplot2)
library(patchwork)
ATAC@meta.data$orig.ident <- factor(
  ATAC@meta.data$orig.ident,
  levels = ident_order
)

# Panel A: TSS Enrichment Distribution by Sample
p1 <- ggplot(ATAC@meta.data, aes(x = TSS.enrichment, fill = orig.ident)) +
  geom_density(alpha = 0.6) +
  #geom_vline(xintercept = 2, linetype = "dashed", color = "red") +
  theme_classic() +
  labs(title = "TSS Enrichment by Sample", x = "TSS Enrichment", y = "Density") +
  theme(legend.position = "top")

# Panel B: Nucleosome Signal Distribution by Sample
p2 <- ggplot(ATAC@meta.data, aes(x = nucleosome_signal, fill = orig.ident)) +
  geom_density(alpha = 0.6) +
  #geom_vline(xintercept = 4, linetype = "dashed", color = "red") +
  theme_classic() +
  labs(title = "Nucleosome Signal by Sample", x = "Nucleosome Signal", y = "Density") +
  theme(legend.position = "top")

# Panel C: ENCODE Blacklist Fraction by Sample
p3 <- ggplot(ATAC@meta.data, aes(x = blacklist_ratio, fill = orig.ident)) +
  geom_density(alpha = 0.6) +
  #geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  theme_classic() +
  labs(title = "ENCODE Blacklist Fraction by Sample", x = "Blacklist Ratio", y = "Density") +
  theme(legend.position = "top")

# Combine the three panels
figure_S1C <- p1 / p2 / p3 + plot_annotation(tag_levels = "A")

# Show the figure
figure_S1C
save_highres_plot(figure_S1C, "Figure1/FigureS1C_TSS_Enrichment_Nucleosome_Signal_Blacklist_Fraction")
#S1C Complete 


# 1D Start
# Set peaks as the default assay if not already
DefaultAssay(ATAC) <- "peaks"
Idents(ATAC) <- "predicted.id"  # Use predicted identities for grouping

# Marker loci for each cell type
marker_genes <- c(
  "Slc1a2",      # Astrocytes
  "Prox1",       # Dentate gyrus neurons
  "Slc17a7",     # Neuronal marker (vGlut2)
  "Plp1",        # Oligodendrocytes
  "Mbp",         # Oligodendrocytes
  "Mog",         # Oligodendrocytes
  "Pdgfra"       # OPCs / Oligodendrocyte-lineage
)

# Plot and arrange in 2 columns
p <- CoveragePlot(
  object = ATAC,
  region = marker_genes,
  annotation = TRUE,
  group.by = "predicted.id",
  extend.upstream = 2000,
  extend.downstream = 2000,
  ncol = 2
)

# Save the figure
save_highres_plot(p, "Figure1/Figure1D_MarkerLoci_ChromatinAccessibility", width = 12, height = 10, dpi = 600, format = "png")
# 1D Complete

# 1E Start
anchors <- readRDS("../ATAC_analysis_with_combined_peaks/anchors_from_FindIntegrationAnchors.RDS")
# Integrate data
integrated <- IntegrateData(anchorset = anchors)
# Dimensionality reduction on integrated data
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, dims = 1:30)
# Plot co-embedding
DimPlot(integrated, group.by = "dataset_source", raster = FALSE)
coembed <- integrated
coembed$global_identity <- ifelse(!is.na(coembed$global_identity), coembed$global_identity, coembed$predicted.id)

# Generate UMAP plot
p <- DimPlot(integrated, group.by = "dataset_source", raster = FALSE, label = FALSE) +
  scale_color_manual(
    values = c("ATAC_source" = "maroon", "RNA_source" = "steelblue4"),
    labels = c("ATAC_source" = "ATAC", "RNA_source" = "RNA")
  ) +
  labs(color = NULL)  # Optional: remove legend title

# Save with your general function
save_highres_plot(p, "Figure1/Figure1C_Coembedding_UMAP_ATAC_RNA_together", format = "png")

coembed@meta.data$global_identity <- factor(coembed@meta.data$global_identity, levels = cell_order)
# Generate UMAP plot
p <- DimPlot(coembed, group.by = "global_identity", raster = FALSE, label = TRUE, cols = palette_umap) +
  labs(color = NULL)  # Optional: remove legend title

# Save with your general function
save_highres_plot(p, "Figure1/Figure1D_Coembedding_UMAP_ATAC_RNA_together_with_identity_label_TRUE", format = "png")

p <- DimPlot(coembed, group.by = "global_identity", raster = FALSE, label = FALSE, cols = palette_umap) +
  labs(color = NULL)  # Optional: remove legend title

# Save with your general function
save_highres_plot(p, "Figure1/Figure1D_Coembedding_UMAP_ATAC_RNA_together_with_identity_label_FALSE", format = "png")



# Calculate overall score
overall_score <- ATAC@meta.data$prediction.score.max
overall_df <- data.frame(
  predicted.id = "Overall",
  prediction.score.max = overall_score
)

# Extract actual data and append overall
meta_df <- ATAC@meta.data %>%
  dplyr::select(predicted.id, prediction.score.max)

# Combine both
plot_df <- bind_rows(meta_df, overall_df)

# Add "Overall" to end of levels
cell_order_with_overall <- c(cell_order, "Overall")
plot_df$predicted.id <- factor(plot_df$predicted.id, levels = cell_order_with_overall)

# Add black color for Overall
palette_umap_with_overall <- c(palette_umap, "Overall" = "black")

# Plot
p <- ggplot(plot_df, aes(x = predicted.id, y = prediction.score.max, fill = predicted.id)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8, color = NA) +
  geom_boxplot(width = 0.1, outlier.size = 0.4, color = "black", alpha = 0.5) +
  geom_hline(yintercept = median(overall_score), linetype = "dashed", color = "black", size = 0.6) +
  # annotate("text", x = 0.5, y = median(overall_score) + 0.02,
  #          label = "median: 0.98", hjust = 0, size = 4, color = "black") +
  scale_fill_manual(values = palette_umap_with_overall) +
  labs(
    x = NULL,
    y = "Prediction score"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 18, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    plot.title = element_text(size = 16, face = "bold", color = "black", hjust = 0.5),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    panel.grid.major = element_line(size = 0.2, color = "gray85"),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 40)
  )

save_highres_plot(p, "Figure1/Figure1F_Prediction_Score_distribution", format = "png", height = 10)



# 1E Complete
ATAC_backup <- ATAC

DefaultAssay(ATAC) <- "peaks"
Idents(ATAC) <- "predicted.id"

cluster_markers <- FindAllMarkers(ATAC, only.pos = TRUE, min.pct = 0.55, logfc.threshold = 0.55, test.use = "LR", latent.vars = 'nCount_peaks')

cluster_markers <- FindAllMarkers(ATAC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

library(dplyr)
cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top10

heatmap <- DoHeatmap(
  object = ATAC,
  features = top10$gene,
  group.by = "predicted.id",
  slot = "data",
  group.colors = palette_umap,
  label = FALSE
) +
  scale_fill_gradientn(colors = c("navy", "yellow", "red2")) +
  ggtitle("Top Variable Peaks Across All Cells") +
  theme(
    axis.text.y = element_blank(),        # hides cell labels
    axis.text.x = element_blank(),        # hides peak labels
    axis.ticks.y = element_blank(),       # remove ticks
    axis.ticks.x = element_blank(),       # remove ticks
    axis.title.x = element_blank(),       # remove axis title
    axis.title.y = element_blank(),       # remove axis title
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

save_highres_plot(heatmap, "Figure1/Figure1G_heatmap_for_peaks", width = 14)


DoHeatmap(ATAC, features = top10$gene, slot = "data") + NoLegend()

 

# Set correct default assay
DefaultAssay(ATAC) <- "peaks"

#top_peaks <- FindTopFeatures(ATAC, min.cutoff = "10")

# Select top variable peaks
ATAC <- FindTopFeatures(ATAC)

top_peaks <- head(VariableFeatures(ATAC), 50000) 

# Set identity class
Idents(ATAC) <- "predicted.id"
  
# Heatmap  
heatmap <- DoHeatmap(
  object = ATAC,
  features = top_peaks,
  group.by = "predicted.id",
  slot = "data",
  group.colors = palette_umap,
  label = FALSE
) +
  scale_fill_gradientn(colors = c("navy", "yellow", "red2")) +
  ggtitle("Top Variable Peaks Across All Cells") +
  theme(
    axis.text.y = element_blank(),        # hides cell labels
    axis.text.x = element_blank(),        # hides peak labels
    axis.ticks.y = element_blank(),       # remove ticks
    axis.ticks.x = element_blank(),       # remove ticks
    axis.title.x = element_blank(),       # remove axis title
    axis.title.y = element_blank(),       # remove axis title
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


save_highres_plot(heatmap, "Figure1/Figure1G_heatmap_for_peaks")


library(ComplexHeatmap)

# Get top peaks (e.g., by variance)
top_peaks <- head(VariableFeatures(ATAC), 10000)

# Extract normalized accessibility matrix
mat <- GetAssayData(ATAC, slot = "data")[top_peaks, ]

# Optional: scale (z-score) across cells
scaled_mat <- t(scale(t(as.matrix(mat))))  # rows = peaks, cols = cells

scaled_mat <- as.matrix(mat)

# Subset to reasonable number of cells (optional)
set.seed(42)
cell_subset <- sample(colnames(scaled_mat), 300)
scaled_mat <- scaled_mat[, cell_subset]

cell_types <- ATAC$predicted.id[cell_subset]

col_annot <- HeatmapAnnotation(
  CellType = cell_types,
  col = list(CellType = palette_umap),
  show_annotation_name = FALSE
)

heatmap <- Heatmap(
  matrix = scaled_mat,
  name = "Accessibility",
  top_annotation = col_annot,
  cluster_rows = TRUE,
  cluster_columns = FALSE,  # optional: disable to keep ladder effect
  show_row_names = FALSE,
  show_column_names = FALSE,
  column_split = cell_types,  # optional: separate by predicted.id
  col = colorRamp2(c(-2, 0, 2), c("navy", "yellow", "red2")),
  heatmap_legend_param = list(title = "Accessibility")
)



# Genes to visualize
genes_to_plot <- c("Slc1a2", "Slc17a7", "Plp1", "Mog")

# Generate plots
coverage_plots <- lapply(genes_to_plot, function(gene) {
  CoveragePlot(
    object = ATAC,
    region = gene,
    assay = "peaks",
    group.by = "predicted.id",
    extend.upstream = 2000,
    extend.downstream = 2000,
    annotation = TRUE,
    peaks = TRUE,
    links = FALSE
  ) + ggtitle(gene)
})

# Combine and display
library(patchwork)
wrap_plots(coverage_plots, ncol = 2)


# Generate the publication-quality plot for Slc1a2
Slc1a2 <- CoveragePlot(
  object = ATAC,
  region = "chr2-102650000-102670000",
  assay = "peaks",
  group.by = "predicted.id",
  extend.upstream = 2000,
  extend.downstream = 2000,
  annotation = TRUE,
  peaks = TRUE,
  links = FALSE,
  cols = palette_umap  # assign your color palette here
) +
  ggtitle("Slc1a2") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  ) & scale_fill_manual(values = palette_umap)

save_highres_plot(Slc1a2, "Figure1/Figure1H_Slc1a2_11_celltypes_coverage_plot")

Slc17a7 <- CoveragePlot(
  object = ATAC,
  region = "chr7-45165000-45175000",
  assay = "peaks",
  group.by = "predicted.id",
  extend.upstream = 2000,
  extend.downstream = 500,
  annotation = TRUE,
  peaks = TRUE,
  links = FALSE,
  cols = palette_umap  # assign your color palette here
) +
  ggtitle("Slc17a7") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  ) & scale_fill_manual(values = palette_umap)

Slc17a7

save_highres_plot(Slc17a7, "Figure1/Figure1H_Slc17a7_11_celltypes_coverage_plot")


Plp1 <- CoveragePlot(
  object = ATAC,
  region = "chrX-136820000-136825000",
  assay = "peaks",
  group.by = "predicted.id",
  extend.upstream = 500,
  extend.downstream = 500,
  annotation = TRUE,
  peaks = TRUE,
  links = FALSE,
  cols = palette_umap  # assign your color palette here
) +
  ggtitle("Plp1") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  ) & scale_fill_manual(values = palette_umap)

Plp1

save_highres_plot(Plp1, "Figure1/Figure1H_Plp1_11_celltypes_coverage_plot")
 
Mog <- CoveragePlot(
  object = ATAC,
  region = "chr17-37020000-37025000",
  assay = "peaks",
  group.by = "predicted.id",
  extend.upstream = 1000,
  extend.downstream = 1000,
  annotation = TRUE,
  peaks = TRUE,
  links = FALSE,
  cols = palette_umap  # assign your color palette here
) +
  ggtitle("Mog") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black", face = "bold"),
    strip.text = element_text(size = 12, face = "bold", color = "black"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 10, 10, 10)
  ) & scale_fill_manual(values = palette_umap)

Mog

save_highres_plot(Mog, "Figure1/Figure1H_Mog_11_celltypes_coverage_plot")

