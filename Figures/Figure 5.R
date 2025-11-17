library(Seurat)
library(Signac)
library(ggplot2)
library(scplotter)
library(tidyr)
library(dplyr)
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
 
anchors <- readRDS("E:/Projects/ATAC-EPO/ATAC_analysis_with_individual_peaks/Pyramidal_Lineage/anchors_from_FindIntegrationAnchors.RDS")
# Integrate data
integrated <- IntegrateData(anchorset = anchors)
# Dimensionality reduction on integrated data
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated)
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, dims = 1:30)
# Plot co-embedding
DimPlot(integrated, group.by = "dataset_source", raster = FALSE)
#coembed <- integrated
#coembed$global_identity <- ifelse(!is.na(coembed$global_identity), coembed$global_identity, coembed$predicted.id)

p <- DimPlot(integrated, group.by = "dataset_source", raster = FALSE, label = FALSE) +
  scale_color_manual(
    values = c("ATAC_source" = "maroon", "RNA_source" = "steelblue4"),
    labels = c("ATAC_source" = "ATAC", "RNA_source" = "RNA")
  ) +
  labs(color = NULL)  # Optional: remove legend title

save_highres_plot(p, filename = "Figure5/ATAC_RNA_pyramidal_UMAP")


ATAC <- readRDS("E:/Projects/ATAC-EPO/ATAC_analysis_with_individual_peaks/Pyramidal_Lineage/ATAC_pyramidal_curated.rds")
set.seed(0)
ATAC_pyramidal <- ATAC
before_harmony <- DimPlot(ATAC_pyramidal, label = TRUE) + ggtitle("scATAC-seq (Unintegrated)")
library(harmony)
ATAC_pyramidal <- RunHarmony(object = ATAC_pyramidal, group.by.vars = 'orig.ident', reduction.use = 'lsi', assay.use = 'peaks', project.dim = FALSE)
ATAC_pyramidal <- RunUMAP(ATAC_pyramidal, dims = 2:30, reduction = 'harmony')
after_harmony <- DimPlot(ATAC_pyramidal, label = TRUE) + ggtitle("scATAC-seq (Integrated using Harmony)")
before_harmony | after_harmony

# Rename the levels
ATAC_pyramidal@meta.data$final_celltype <- gsub(
  pattern = "New_formed_Ventral_CA3_serotonin",
  replacement = "Newly_formed_Ventral_CA3_serotonin",
  x = ATAC_pyramidal@meta.data$final_celltype
)

ATAC_pyramidal@meta.data$final_celltype <- gsub(
  pattern = "New_formed_CA3_superficial_Mbp",
  replacement = "Newly_formed_CA3_superficial_Mbp",
  x = ATAC_pyramidal@meta.data$final_celltype
)

# Remove underscores
ATAC_pyramidal$final_celltype <- gsub("_", " ", ATAC_pyramidal$final_celltype)

ATAC_pyramidal@meta.data$final_celltype <- gsub(
  pattern = "chandelier cell",
  replacement = "Chandelier cell",
  x = ATAC_pyramidal@meta.data$final_celltype
)

ATAC_pyramidal@meta.data$final_celltype <- gsub(
  pattern = "Deep layer thin-CA3",
  replacement = "Deep layer thin CA3",
  x = ATAC_pyramidal@meta.data$final_celltype
)

ATAC_pyramidal@meta.data$final_celltype <- gsub(
  pattern = "New formed deep Foxp2",
  replacement = "Newly formed deep Foxp2",
  x = ATAC_pyramidal@meta.data$final_celltype
)

ATAC_pyramidal@meta.data$final_celltype <- gsub(
  pattern = "Newly formed migrating-deep Sox5",
  replacement = "Newly formed migrating deep Sox5",
  x = ATAC_pyramidal@meta.data$final_celltype
)

Idents(ATAC_pyramidal) <- "final_celltype"
DimPlot(ATAC_pyramidal, label = TRUE) + ggtitle("scATAC-seq (Integrated using Harmony)")



new_order <- c(
  # Group 1: Newly formed
  "Newly formed migrating superficial",
  "Newly formed deep Foxp2",
  "Newly formed migrating serotonin firing Cux2 ",
  "Newly formed migrating serotonin firing Dok5",
  "Newly formed migrating deep Sox5",
  "Newly formed migrating serotonin firing Unc5d",
  "Newly formed Ventral CA3 serotonin",
  "Newly formed CA3 superficial Mbp",
  "Ventral CA1 migrating",
  
  # Group 2: Mature
  "CA1 dorsal",
  "CA1 retrohippo",
  "CA3 dorsal",
  "Deep layer thin CA3",
  "CA2",
  
  # Group 3: Other (zero in helper)
  "Chandelier cell",
  "Foxp2 Bnc2 prelineage",
  "Cajal like"
)

ATAC_pyramidal@meta.data$final_celltype <- factor(ATAC_pyramidal$final_celltype, levels = new_order)

# Load
library(paletteer)

# Get 17 colors from the "alphabet2" palette
celltype_colors <- DiscretePalette(n = 17, palette = "alphabet2")

celltype_colors[celltype_colors == "#E2E2E2"] <- "#194A8D"  # A rich blue-gray 

# Name them according to your cell type levels
names(celltype_colors) <- levels(ATAC_pyramidal$final_celltype)

# Use in DimPlot
p <- DimPlot(
  ATAC_pyramidal,
  group.by = "final_celltype",
  cols = celltype_colors,
  label = TRUE, raster = FALSE
) + labs(color = NULL)  # Optional: remove legend title





save_highres_plot(p, filename = "Figure5/pyramidal_UMAP")
saveRDS(celltype_colors, file = "Figure5/celltype_colors.RDS")

ATAC_pyramidal@meta.data <- ATAC_pyramidal@meta.data %>%
  mutate(
    condition = case_when(
      startsWith(orig.ident, "A") ~ "EPO",
      startsWith(orig.ident, "B") ~ "Placebo",
      TRUE ~ NA_character_
    )
  )




ATAC <- readRDS("Data/Pyramidal/EPO_vs_Placebo/Newly_formed_vs_Mature_ATAC_individual_peaks_for_Pyramidal_neurons_mature_vs_newly_formed.RDS")
head(ATAC@meta.data)
library(dplyr)
ATAC@meta.data <- ATAC@meta.data %>%
  mutate(
    condition = case_when(
      startsWith(type, "A") ~ "EPO",
      startsWith(type, "B") ~ "Placebo",
      TRUE ~ NA_character_
    )
  )

# Create combined cluster+condition identity
ATAC@meta.data$mature_vs_newlyformed_status_condition <- paste0(ATAC$mature_vs_newlyformed_status, "_", ATAC@meta.data$condition)
Idents(ATAC) <- "mature_vs_newlyformed_status_condition"

# Unique clusters
clusters <- unique(ATAC$mature_vs_newlyformed_status)


# Fix factor labels
ATAC$mature_vs_newlyformed_status <- factor(
  ATAC$mature_vs_newlyformed_status,
  levels = c("newly_formed", "mature"),
  labels = c("newly formed", "mature")
)

# Plot with custom colors
p <- DimPlot(
  ATAC,
  group.by = "mature_vs_newlyformed_status",
  label = FALSE, raster = FALSE,
  split.by = "condition",
  cols = c("newly formed" = "red", "mature" = "blue")
) + labs(color = NULL)


save_highres_plot(p, filename = "Figure5/pyramidal_UMAP_without_labels_newly_formed_vs_mature")

 
# Use in DimPlot
p <- DimPlot(
  ATAC,
  group.by = "mature_vs_newlyformed_status",
  label = TRUE, raster = FALSE, split.by = "condition"
) + labs(color = NULL)  # Optional: remove legend title




# Use in DimPlot
p <- DimPlot(
  ATAC_pyramidal,
  group.by = "final_celltype",
  cols = celltype_colors,
  label = FALSE, raster = FALSE
) + labs(color = NULL)  # Optional: remove legend title

save_highres_plot(p, filename = "Figure5/pyramidal_UMAP_without_labels")

# Use in DimPlot
p <- DimPlot(
  ATAC_pyramidal,
  group.by = "final_celltype",
  cols = celltype_colors,
  label = FALSE,
  raster = FALSE
)

p <- LabelClusters(plot = p, id = "final_celltype", repel = TRUE) +
  NoLegend()



save_highres_plot(p, filename = "Figure5/pyramidal_UMAP_no_legend")



library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(scales)

# Load the DARs CSV
DARs <- read_csv("Data/Pyramidal/Without_Treatment_Separation/Mature_vs_NewlyFormed/DA_peaks_combined_pyramidal_celltypes_Without_treatment_separation_filtered_annotated_for_genes_and_TEs.csv")

# Filter for adjusted p-value < 0.01 and cluster = "Cluster_mature"
filtered_DARs <- DARs %>%
  filter(p_val_adj < 0.01, cluster == "Cluster_newly_formed")
# Filter for adjusted p-value < 0.01 and cluster = "Cluster_mature"
filtered_DARs <- filtered_DARs %>%
  filter(pct.1 >= 0.10 | pct.2 >= 0.10)


# Correct TE overlaps by excluding low-confidence classes or families
exclude_classes <- c("Low_complexity", "Satellite", "rRNA", "tRNA", "scRNA", "snRNA", 
                     "RNA", "RC", "srpRNA", "DNA?", "LTR?", "LINE?", "RC?", "SINE?")

filtered_DARs <- filtered_DARs %>%
  mutate(overlapping_TE = ifelse(
    overlapping_TE & (
      TE_class %in% exclude_classes |
        grepl("\\?$", TE_family)
    ),
    FALSE,
    overlapping_TE
  ))

# Annotate direction of regulation
filtered_DARs <- filtered_DARs %>%
  mutate(
    direction = ifelse(avg_log2FC > 0, "Up", "Down")
  )

# Create summary for All, Promoter, TE overlapping
summary_all <- filtered_DARs %>%
  summarise(
    Total = n(),
    Promoter = sum(overlapping_promoter),
    TE = sum(overlapping_TE)
  ) %>%
  pivot_longer(cols = everything(), names_to = "Category", values_to = "Count") %>%
  mutate(Direction = "All")

# Breakdown by direction (Up/Down)
summary_by_direction <- filtered_DARs %>%
  group_by(direction) %>%
  summarise(
    Total = n(),
    Promoter = sum(overlapping_promoter),
    TE = sum(overlapping_TE)
  ) %>%
  pivot_longer(cols = -direction, names_to = "Category", values_to = "Count") %>%
  dplyr::rename(Direction = direction)

# Combine into one dataframe for plotting
plot_data <- bind_rows(summary_all, summary_by_direction)

plot_data <- plot_data %>%
  mutate(Category = recode(Category,
                           "Total" = "Total DARs",
                           "Promoter" = "Gene-TSS DARs",
                           "TE" = "TE-derived DARs"
  ))

# ✅ Factor ordering
plot_data$Category <- factor(plot_data$Category, levels = c("Total DARs", "Gene-TSS DARs", "TE-derived DARs"))
plot_data$Direction <- factor(plot_data$Direction, levels = c("All", "Up", "Down"))

# ✅ Color palette (match your preferred order)
color_palette <- c("All" = "#999999", "Down" = "steelblue1", "Up" = "hotpink2")

# ✅ Plot
p <- ggplot(plot_data, aes(x = Category, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = Count),
            position = position_dodge(width = 0.7),
            vjust = -0.4,
            size = 4,
            color = "black") +
  scale_fill_manual(values = color_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), labels = comma) +
  labs(
    title = "Chromatin Accessibility Differences (Newly Formed vs Mature)",
    subtitle = "Distribution of Total, Promoter-Overlapping, and TE-Derived DARs",
    x = NULL,
    y = "DARs (Newly Formed vs Mature Neurons)",
    fill = "Regulation"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10)),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 13),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
 
p

save_highres_plot(p, filename = "Figure5/Chromatin Accessibility Differences (Newly formed vs Mature Neurons)", width = 6)









library(dplyr)
library(ggplot2)
library(readr)

# Load DARs
DARs <- read_csv("Data/Pyramidal/Without_Treatment_Separation/Mature_vs_NewlyFormed/DA_peaks_combined_pyramidal_celltypes_without_treatment_separation.csv")

# Filter for cluster "Cluster_newly_formed", expression ≥10%, padj < 0.01
volcano_data <- DARs %>%
  filter(cluster == "Cluster_newly_formed") %>%
  filter(pct.1 >= 0.10 | pct.2 >= 0.10) %>%
  mutate(
    log10padj = -log10(p_val_adj),
    significance = ifelse(p_val_adj < 0.01, "Significant", "Not Significant")
  )
 
# Cap infinite values
volcano_data$log10padj[is.infinite(volcano_data$log10padj)] <-
  max(volcano_data$log10padj[is.finite(volcano_data$log10padj)], na.rm = TRUE)

# Volcano plot
p <- ggplot(volcano_data, aes(x = avg_log2FC, y = log10padj, color = significance)) +
  geom_point(size = 1.3, alpha = 0.7) +
  scale_color_manual(values = c("Significant" = "#D73027", "Not Significant" = "#CCCCCC")) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black", linewidth = 0.4) +
  #geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", color = "blue", linewidth = 0.4) +
  labs(
    title = "Volcano Plot of DARs in Newly Formed Pyramidal Cells",
    subtitle = "Peaks expressed in ≥10% of cells in at least one condition",
    x = expression(log[2]~"Fold Change"),
    y = expression(-log[10]~"Adjusted P-value"),
    color = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black", size = 0.25),
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 14),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )

save_highres_plot(p, filename = "Figure5/Chromatin Accessibility Differences (Newly formed vs Mature Neurons) volcano plot")









library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

# Load the DARs CSV
DARs <- read_csv("Data/Pyramidal/Without_Treatment_Separation/Mature_vs_NewlyFormed/DA_peaks_combined_pyramidal_celltypes_Without_treatment_separation_filtered_annotated_for_genes_and_TEs.csv")
DARs$overlapping_promoter <- NA
DARs$gene_id <- NA
DARs$gene_symbol_old <- DARs$gene_symbol
DARs$gene_symbol <- NA
head(DARs)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# === Parse summit into chr/start/end ===
DARs <- DARs %>%
  separate(summit, into = c("chr", "range"), sep = ":") %>%
  separate(range, into = c("start", "end"), sep = "-") %>%
  mutate(
    start = as.integer(start),
    end = as.integer(end),
    summit_label = paste0(chr, ":", start, "-", end)
  )

# === Convert to GRanges ===
dar_summits <- GRanges(
  seqnames = DARs$chr,
  ranges = IRanges(
    start = floor((DARs$start + DARs$end) / 2),
    end = floor((DARs$start + DARs$end) / 2)
  )
)

# === Load transcript TSSs from txdb ===
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
transcripts_gr <- transcripts(txdb)

# Map TXID → GENEID
tx2gene <- AnnotationDbi::select(
  txdb,
  keys = as.character(mcols(transcripts_gr)$tx_id),
  keytype = "TXID",
  columns = "GENEID"
)
transcripts_gr$gene_id <- tx2gene$GENEID[match(mcols(transcripts_gr)$tx_id, tx2gene$TXID)]

# Map GENEID → SYMBOL
gene_id_to_symbol <- mapIds(
  org.Mm.eg.db,
  keys = transcripts_gr$gene_id,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)
transcripts_gr$symbol <- gene_id_to_symbol

# === Define TSS as 1bp GRanges and expand ±1kb for promoter ===
tss_gr <- GRanges(
  seqnames = seqnames(transcripts_gr),
  ranges = IRanges(
    start = ifelse(strand(transcripts_gr) == "+", start(transcripts_gr), end(transcripts_gr)),
    width = 1
  ),
  strand = strand(transcripts_gr),
  gene_id = as.character(transcripts_gr$gene_id),
  symbol = as.character(transcripts_gr$symbol)
)

promoter_gr <- resize(tss_gr, width = 2001, fix = "center")

# === Find promoter overlaps only ===
overlaps <- findOverlaps(dar_summits, promoter_gr)
DARs$overlapping_promoter <- FALSE
DARs$overlapping_promoter[queryHits(overlaps)] <- TRUE

# === Assign gene_id and symbol only for overlapping promoters ===
DARs$gene_id[queryHits(overlaps)] <- mcols(promoter_gr)$gene_id[subjectHits(overlaps)]
DARs$gene_symbol[queryHits(overlaps)] <- mcols(promoter_gr)$symbol[subjectHits(overlaps)]

# === Done — Optional preview ===
DARs %>% dplyr::select(summit_label, gene_id, gene_symbol, overlapping_promoter) %>% head()

sum(DARs$gene_symbol_old == DARs$gene_symbol, na.rm = TRUE)
mean(DARs$gene_symbol_old == DARs$gene_symbol, na.rm = TRUE) * 100

library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
# Filter for promoter-DARs with increased accessibility in newly formed neurons
promoter_dars_up_new <- DARs %>%
  filter(
    p_val_adj < 0.01,
    overlapping_promoter == TRUE,
    cluster == "Cluster_newly_formed",
    avg_log2FC > 0.0,
    !is.na(gene_symbol)
  )

# Unique gene symbols
genes_up_new <- unique(promoter_dars_up_new$gene_symbol)

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(genes_up_new, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

library(ggplot2)
library(dplyr)
genes <- entrez_ids$ENTREZID
# Combine results from all ontologies (BP, CC, MF) if you haven't already
ego_BP <- enrichGO(genes, OrgDb = org.Mm.eg.db, ont = "BP", keyType      = "ENTREZID",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
ego_CC <- enrichGO(genes, OrgDb = org.Mm.eg.db, ont = "CC", keyType      = "ENTREZID",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)
ego_MF <- enrichGO(genes, OrgDb = org.Mm.eg.db, ont = "MF", keyType      = "ENTREZID",
                   pAdjustMethod = "BH",
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

# Convert to data frame and annotate source
df_BP <- as.data.frame(ego_BP)[1:20, ] %>% mutate(ontology = "biological process")
df_CC <- as.data.frame(ego_CC)[1:20, ] %>% mutate(ontology = "cellular component")
df_MF <- as.data.frame(ego_MF)[1:20, ] %>% mutate(ontology = "molecular function")

# Combine all
go_df <- bind_rows(df_BP, df_CC, df_MF)

# Format factor levels for correct order
go_df$Description <- factor(go_df$Description, levels = rev(go_df$Description))

# Set ontology colors (classic scheme)
ontology_colors <- c(
  "biological process" = "#000080",     # navy
  "cellular component" = "#006400",     # dark green
  "molecular function" = "#DAA520"      # goldenrod
)

# Plot
p <- ggplot(go_df, aes(x = Description, y = Count, fill = ontology)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(
    aes(label = round(-log10(p.adjust), 1)),
    hjust = -0.2, size = 3.5, color = "black"
  ) +
  scale_fill_manual(values = c(
    "biological process" = "#1f78b4",
    "cellular component" = "#33a02c",
    "molecular function" = "#ff7f00"
  )) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title = "GO Enrichment of Promoter DARs (Upregulated in Newly Formed Neurons)",
    subtitle = expression("Numbers represent " * -log[10]~"adjusted p-values"),
    x = NULL,
    y = "Number of Genes",
    fill = "Ontology"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.title = element_text(color = "black", face = "bold"),
    plot.title = element_text(color = "black", face = "bold", size = 16, hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(color = "black", size = 12, hjust = 0.5),
    panel.grid = element_blank(),
    legend.position = "right"
  )

save_highres_plot(p, filename = "Figure5/GO analysis for Upregulated in Newly Formed Neurons")
 
# GO enrichment (Biological Process)
ego <- enrichGO(
  gene         = entrez_ids$ENTREZID,
  OrgDb        = org.Mm.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Barplot of top GO terms
barplot(ego, showCategory = 20, title = "GO Enrichment: Newly Formed Promoter DARs", font.size = 12)



library(tidyverse)

# Load DARs table
DARs <- read.csv("Data/Pyramidal/EPO_vs_Placebo/Newly_formed_vs_Mature_DA_peaks_combined_pyramidal_celltypes_EPO_vs_Placebo_separation_filtered_annotated_for_genes_and_TEs.csv")
DARs$cluster <- paste0("Cluster_", DARs$cell_stage)
# Filter for significant DARs
DARs_sig <- DARs %>%
  filter(p_val_adj < 0.01)

# Annotate direction
DARs_sig <- DARs_sig %>%
  mutate(direction = case_when(
    avg_log2FC > 0 ~ "Up upon EPO",
    avg_log2FC < 0 ~ "Down upon EPO",
    TRUE ~ NA_character_
  ))

# Summarize and rename clusters
dar_percent <- DARs_sig %>%
  filter(!is.na(direction)) %>%
  group_by(cluster = case_when(
    cluster == "Cluster_newly_formed" ~ "Newly Formed",
    cluster == "Cluster_mature" ~ "Mature",
    TRUE ~ cluster
  ), direction) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(
    total = sum(n),
    percent = n / total * 100,
    label = paste0(round(percent, 1), "% (", n, ")")
  ) %>%
  ungroup()

# Plot
p <- ggplot(dar_percent, aes(x = cluster, y = percent, fill = direction)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = 0.7) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 4.5, fontface = "bold", color = "black") +
  scale_fill_manual(
    values = c("Up upon EPO" = "hotpink2",  # Green
               "Down upon EPO" = "steelblue1"),  # Orange
    name = "Accessibility Change",
    guide = guide_legend(reverse = TRUE)
  ) +
  labs(
    title = "Proportion of Significant DARs per Cluster",
    x = NULL,
    y = "Percentage of Significant DARs"
  ) +
  coord_flip() +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black", face = "bold"),
    axis.title = element_text(color = "black", face = "bold"),
    plot.title = element_text(face = "bold", size = 16, color = "black"),
    legend.title = element_text(face = "bold", color = "black"),
    legend.text = element_text(color = "black"),
    panel.grid = element_blank(),
    legend.position = "top"
  )


save_highres_plot(p, filename = "Figure5/Proportion of Significant DARs per Cluster", height = 3)


  



DARs <- read.csv("Data/Pyramidal/EPO_vs_Placebo/Newly_formed_vs_Mature_DA_peaks_combined_pyramidal_celltypes_EPO_vs_Placebo_separation_filtered_annotated_for_genes_and_TEs.csv")
ATAC <- readRDS("Data/Pyramidal/EPO_vs_Placebo/Newly_formed_vs_Mature_ATAC_individual_peaks_for_Pyramidal_neurons_mature_vs_newly_formed.RDS")
head(ATAC@meta.data)
library(dplyr)
ATAC@meta.data <- ATAC@meta.data %>%
  mutate(
    condition = case_when(
      startsWith(type, "A") ~ "EPO",
      startsWith(type, "B") ~ "Placebo",
      TRUE ~ NA_character_
    )
  )

# Create combined cluster+condition identity
ATAC@meta.data$mature_vs_newlyformed_status_condition <- paste0(ATAC$mature_vs_newlyformed_status, "_", ATAC@meta.data$condition)
Idents(ATAC) <- "mature_vs_newlyformed_status_condition"

# Unique clusters
clusters <- unique(ATAC$mature_vs_newlyformed_status)


da_peaks_per_cluster_list <- list()
set.seed(0)
for (clust in clusters) {
  message(sprintf("Processing cluster: %s", clust))
  
  epo_group <- paste0(clust, "_EPO")
  placebo_group <- paste0(clust, "_Placebo")
  
  if (all(c(epo_group, placebo_group) %in% Idents(ATAC))) {
    n_epo <- sum(Idents(ATAC) == epo_group)
    n_placebo <- sum(Idents(ATAC) == placebo_group)
    
    if (n_epo >= 3 && n_placebo >= 3) {
      message("Running DA for cluster: ", clust)
      
      da_peaks <- FindMarkers(
        object = ATAC,
        ident.1 = epo_group,
        ident.2 = placebo_group,
        only.pos = FALSE,
        min.pct = 0.01,
        test.use = "LR",
        latent.vars = "nCount_peaks"
      )
      
      if (nrow(da_peaks) > 0) {
        da_peaks$comparison <- paste0(epo_group, "_vs_", placebo_group)
        da_peaks$peak_id <- rownames(da_peaks)
        da_peaks_per_cluster_list[[clust]] <- da_peaks
      } else {
        message("No DA peaks found for cluster: ", clust)
      }
    } else {
      message("Skipping cluster ", clust, ": not enough cells (EPO = ", n_epo, ", Placebo = ", n_placebo, ")")
    }
  } else {
    message("Skipping cluster ", clust, ": missing group identities")
  }
}

saveRDS(da_peaks_per_cluster_list, file = "Figure5/da_peaks_per_cluster_list_file_0point01_newlyformed")

saveRDS(da_peaks_per_cluster_list, file = "Figure5/da_peaks_per_cluster_list_file_0point01_mature")


# Read in
df_newlyformed <- readRDS("Figure5/da_peaks_per_cluster_list_file_0point01_newlyformed")
df_mature      <- readRDS("Figure5/da_peaks_per_cluster_list_file_0point01_mature")

# Extract the data frames from the lists and add a group column
newlyformed_df <- df_newlyformed$newly_formed
newlyformed_df$cell_stage <- "newly_formed"

mature_df <- df_mature$mature
mature_df$cell_stage <- "mature"

# Merge into a single data frame
merged_df <- rbind(newlyformed_df, mature_df)

# Check result
head(merged_df)

newly_formed_and_mature_EPO_vs_Placebo_sig <- merged_df[merged_df$p_val_adj < 0.01,]
table(newly_formed_and_mature_EPO_vs_Placebo_sig$comparison)

head(newly_formed_and_mature_EPO_vs_Placebo_sig)

da_peaks_combined <- newly_formed_and_mature_EPO_vs_Placebo_sig

# Reorder columns: peak right after "row index" (rownames are numeric now)
da_peaks_combined$peak <- da_peaks_combined$peak_id
da_peaks_combined <- da_peaks_combined[, c("peak", setdiff(names(da_peaks_combined), "peak"))]

da_peaks_combined_seurat_level <- da_peaks_combined #even though it is called seurat_level, it is at 11 level. 

peak_annotation_df <- read_csv("E:/Projects/ATAC-EPO/ATAC_analysis_with_individual_peaks/ATAC_peak_annotation_for_genes_and_TEs_for_individual_peak_file_+1_-1kbfrompromoters.csv")  # change if your annotation file is named differently

# Step 2: Convert ':' to '-' in full_peak to match ...1 column
peak_annotation_df <- peak_annotation_df %>%
  mutate(full_peak_hyphen = gsub(":", "-", full_peak))

# Step 3: Join on peak coordinates
combined_df <- da_peaks_combined_seurat_level %>%
  left_join(peak_annotation_df, by = c("peak" = "full_peak_hyphen"))

# Step 4: Filter by adjusted p-value
filtered_df <- combined_df %>%
  filter(p_val_adj < 0.01)

# Step 5: Save to CSV
filtered_df <- cbind(rownames = rownames(filtered_df), filtered_df)[, c("rownames", setdiff(names(filtered_df), "rownames"))]
da_peaks_combined_seurat_level <- cbind(rownames = rownames(da_peaks_combined_seurat_level), da_peaks_combined_seurat_level)[, c("rownames", setdiff(names(da_peaks_combined_seurat_level), "rownames"))]
write_csv(filtered_df, "Figure5/DA_peaks_combined_pyramidal_celltypes_EPO_vs_Placebo_separation_filtered_annotated_for_genes_and_TEs.csv")
saveRDS(filtered_df, "Figure5/DA_peaks_combined_pyramidal_celltypes_EPO_vs_Placebo_separation_filtered_annotated_for_genes_and_TEs.RDS")
saveRDS(merged_df, "Figure5/DA_peaks_combined_pyramidal_celltypes_EPO_vs_Placebo_separation.RDS")
#write_csv(da_peaks_combined_seurat_level, file = "DA_peaks_combined_pyramidal_cellclusters_EPO_vs_Placebo.csv")
#saveRDS(da_peaks_per_cluster_list, file = "DA_peaks_per_cluster_only.pos.FALSE_pyramidal_cellclusters_EPO_vs_Placebo.rds")





library(tidyverse)
library(pheatmap)


all_da_peaks <- da_peaks_per_cluster_list %>%
  imap_dfr(~ .x %>% mutate(cluster = .y, peak_id = rownames(.x)))

# Create wide-format matrix with avg_log2FC per group
heatmap_matrix_all <- all_da_peaks %>%
  filter(cluster %in% c("mature", "newly_formed")) %>%
  group_by(peak_id, cluster) %>%
  summarize(avg_log2FC = mean(avg_log2FC, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = cluster, values_from = avg_log2FC, values_fill = 0)

# Convert to matrix
heatmap_mat_all <- heatmap_matrix_all %>%
  column_to_rownames("peak_id") %>%
  as.matrix()

# Final check: should be two columns
print(dim(heatmap_mat_all))  # rows = peaks, cols = 2 ("mature", "newly_formed")

# Cap values at ±2
heatmap_mat_all_clipped <- heatmap_mat_all
heatmap_mat_all_clipped[heatmap_mat_all_clipped > 2] <- 2
heatmap_mat_all_clipped[heatmap_mat_all_clipped < -2] <- -2

# Rename columns (remove underscores) and reorder
col_order <- c("newly formed", "mature")

colnames(heatmap_mat_all_clipped) <- gsub("_", " ", colnames(heatmap_mat_all_clipped))
heatmap_mat_all_clipped <- heatmap_mat_all_clipped[, col_order]


library(grid)

pheatmap(
  mat = heatmap_mat_all_clipped,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-2, 2, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize = 12,
  main = "Differential Accessibility - log2FC",
  show_rownames = FALSE,
  show_colnames = TRUE,
  border_color = NA,
  cellwidth = NA,
  cellheight = NA,
  legend = TRUE,
  fontsize_col = 12,
  fontsize_row = 0
)

# Add right-side label
grid.text(
  "Significant peak in at least one cluster",
  x = unit(0.99, "npc"),
  y = 0.5,
  rot = 270,
  gp = gpar(fontsize = 12)
)

# Set the legend title by assigning a name to the matrix
attr(heatmap_mat_all_clipped, "name") <- "log2FC"


heatmap_mat_all_clipped <- heatmap_mat_all_clipped[, c("newly formed", "mature")]
# Save as PNG
png("Figure5/DA_peaks_heatmap.png", width = 1200*5, height = 1600*5, res = 1200)

# Plot + label
ht <- Heatmap(
  heatmap_mat_all_clipped,
  name = "log2FC",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  column_title = "Differential Accessibility - log2FC",
  column_names_gp = gpar(fontsize = 12),
  row_dend_side = "left",  # explicitly place the row dendrogram
  heatmap_legend_param = list(
    title_position = "lefttop-rot",
    legend_direction = "vertical",
    legend_side = "left"
  )
)

# Draw the heatmap with left legend
draw(ht, heatmap_legend_side = "left")

# Finish
dev.off()




# 1. Save peak order from log2FC heatmap
peak_order <- rownames(heatmap_mat_all_clipped)

# 2. Build wide-format matrix of padj values
padj_matrix_all <- all_da_peaks %>%
  filter(cluster %in% c("mature", "newly_formed")) %>%
  group_by(peak_id, cluster) %>%
  summarize(p_val_adj = mean(p_val_adj, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = cluster, values_from = p_val_adj, values_fill = 1)

# 3. Convert to matrix and reorder
padj_mat <- padj_matrix_all %>%
  column_to_rownames("peak_id") %>%
  as.matrix()

padj_mat_ordered <- padj_mat[peak_order, ]

# 4. Transform to -log10(padj), cap at 20
padj_mat_log10 <- -log10(padj_mat_ordered)
padj_mat_log10[is.infinite(padj_mat_log10)] <- 20
padj_mat_log10[padj_mat_log10 > 20] <- 20

# 5. Clean column names and reorder
colnames(padj_mat_log10) <- gsub("_", " ", colnames(padj_mat_log10))
padj_mat_log10 <- padj_mat_log10[, c("newly formed", "mature")]

# 6. Save plot as PNG
png("Figure5/DA_peaks_padj_heatmap.png", width = 1200*5, height = 1600*5, res = 1200)

# 7. Draw heatmap
# Define the color function
col_fun <- colorRamp2(c(0, 10, 20), c("white", "orange", "red"))

# Create heatmap object
ht_padj <- Heatmap(
  padj_mat_log10,
  name = "-log10(padj)",
  col = col_fun,
  cluster_rows = FALSE,              # preserve row order
  cluster_columns = FALSE,           
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 12),
  heatmap_legend_param = list(
    title = "-log10(padj)",
    legend_direction = "vertical",
    title_position = "lefttop-rot"
  ),
  column_title = "-log10(Adjusted P-value)"
)

# Draw with legend on the left
draw(ht_padj, heatmap_legend_side = "right")

# 9. Close PNG device
dev.off()


library(ggplot2)
library(dplyr)





DARs <- read.csv("Data/Pyramidal/EPO_vs_Placebo/Newly_formed_vs_Mature_DA_peaks_combined_pyramidal_celltypes_EPO_vs_Placebo_separation_filtered_annotated_for_genes_and_TEs.csv")
DARs$cluster <- DARs$cell_stage
ATAC <- readRDS("Data/Pyramidal/EPO_vs_Placebo/Newly_formed_vs_Mature_ATAC_individual_peaks_for_Pyramidal_neurons_mature_vs_newly_formed.RDS")
head(ATAC@meta.data)
library(dplyr)
ATAC@meta.data <- ATAC@meta.data %>%
  mutate(
    condition = case_when(
      startsWith(type, "A") ~ "EPO",
      startsWith(type, "B") ~ "Placebo",
      TRUE ~ NA_character_
    )
  )

# Create combined cluster+condition identity
ATAC@meta.data$mature_vs_newlyformed_status_condition <- paste0(ATAC$mature_vs_newlyformed_status, "_", ATAC@meta.data$condition)
Idents(ATAC) <- "mature_vs_newlyformed_status_condition"

# Unique clusters
clusters <- unique(ATAC$mature_vs_newlyformed_status)


da_peaks_per_cluster_list <- list()
set.seed(0)
for (clust in clusters) {
  message(sprintf("Processing cluster: %s", clust))
  
  epo_group <- paste0(clust, "_EPO")
  placebo_group <- paste0(clust, "_Placebo")
  
  if (all(c(epo_group, placebo_group) %in% Idents(ATAC))) {
    n_epo <- sum(Idents(ATAC) == epo_group)
    n_placebo <- sum(Idents(ATAC) == placebo_group)
    
    if (n_epo >= 3 && n_placebo >= 3) {
      message("Running DA for cluster: ", clust)
      
      da_peaks <- FindMarkers(
        object = ATAC,
        ident.1 = epo_group,
        ident.2 = placebo_group,
        only.pos = FALSE,
        min.pct = 0.00,
        logfc.threshold = 0.00,
        test.use = "LR",
        latent.vars = "nCount_peaks",
        features = unique(DARs$peak)
      )
      
      if (nrow(da_peaks) > 0) {
        da_peaks$comparison <- paste0(epo_group, "_vs_", placebo_group)
        da_peaks$peak_id <- rownames(da_peaks)
        da_peaks_per_cluster_list[[clust]] <- da_peaks
      } else {
        message("No DA peaks found for cluster: ", clust)
      }
    } else {
      message("Skipping cluster ", clust, ": not enough cells (EPO = ", n_epo, ", Placebo = ", n_placebo, ")")
    }
  } else {
    message("Skipping cluster ", clust, ": missing group identities")
  }
}

saveRDS(da_peaks_per_cluster_list, file = "Figure5/Volcano_plot_DARs_EPO_vs_Placebo.RDS")

library(tidyverse)
library(pheatmap)
 

all_da_peaks <- da_peaks_per_cluster_list %>%
  imap_dfr(~ .x %>% mutate(cluster = .y, peak_id = rownames(.x)))




readRDS("Figure5/da_peaks_per_cluster_list_file_0point01_newlyformed")$newly_formed -> newlyformed_df
readRDS("Figure5/da_peaks_per_cluster_list_file_0point01_mature")$mature -> mature_df
colnames(newlyformed_df)
colnames(mature_df)

# Add an origin column (so you know from which group they came)
newlyformed_df$cluster <- "newly_formed"
mature_df$cluster <- "mature"

# Combine
combined_df <- bind_rows(newlyformed_df, mature_df)
all_da_peaks <- combined_df
all_da_peaks <- all_da_peaks[all_da_peaks$pct.1 > 0.10 | all_da_peaks$pct.2 > 0.10, ]

# Prepare and clean cluster names
volcano_data <- all_da_peaks %>%
  filter(cluster %in% c("mature", "newly_formed")) %>%
  mutate(
    cluster = recode(cluster, "newly_formed" = "newly formed"),
    cluster = factor(cluster, levels = c("newly formed", "mature")),
    log10padj = -log10(p_val_adj),
    significance = ifelse(p_val_adj < 0.01, "Significant", "Not Significant")
  )

# Cap infinite padj values
volcano_data$log10padj[is.infinite(volcano_data$log10padj)] <- 
  max(volcano_data$log10padj[is.finite(volcano_data$log10padj)], na.rm = TRUE)

# Plot
p <- ggplot(volcano_data, aes(x = avg_log2FC, y = log10padj, color = significance)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_manual(values = c("Significant" = "#D73027", "Not Significant" = "#CCCCCC")) +
  facet_wrap(~ cluster, scales = "free_y") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey35", linewidth = 0.4) +   # <- log2FC > 0.5
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "grey35", linewidth = 0.4) + # <- optional: log2FC < -0.5
  labs(
    title = "Differential Accessibility by Cell Type",
    subtitle = "Peaks expressed in ≥10% of cells in at least one condition",
    x = expression(log[2]~"Fold Change"),
    y = expression(-log[10]~"Adjusted P-value"),
    color = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black", size = 0.25),
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(face = "bold", size = 16),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24)
  )

print(p)


save_highres_plot(p, filename = "Figure5/Volcanoplot2", width = 16)

 




DARs <- read.csv("Data/Pyramidal/EPO_vs_Placebo/DA_peaks_combined_pyramidal_cellclusters_EPO_vs_Placebo.csv")
head(DARs)

library(ggplot2)
library(dplyr)

# Prepare the data
volcano_data <- DARs %>%
  mutate(
    cluster = gsub("Cluster_", "", cluster),
    cluster = gsub("_", " ", cluster),
    cluster = factor(cluster, levels = c("newly formed", "mature")),
    log10padj = -log10(p_val_adj),
    significance = ifelse(p_val_adj < 0.01, "Significant", "Not Significant")
  )

# Cap -log10(padj) values to avoid Inf
volcano_data$log10padj[is.infinite(volcano_data$log10padj)] <- 
  max(volcano_data$log10padj[is.finite(volcano_data$log10padj)], na.rm = TRUE)

# Plot
p <- ggplot(volcano_data, aes(x = avg_log2FC, y = log10padj, color = significance)) +
  geom_point(size = 1.2, alpha = 0.7) +
  scale_color_manual(values = c("Significant" = "#D73027", "Not Significant" = "#CCCCCC")) +
  facet_wrap(~ cluster, scales = "free_y") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black", linewidth = 0.4) +
  labs(
    title = "Differential Accessibility by Cell Type",
    x = "log2 Fold Change",
    y = "-log10 Adjusted P-value",
    color = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey90", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black", size = 0.25),
    legend.position = "top",
    legend.justification = "left",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    strip.text = element_text(face = "bold", size = 16),
    axis.title.x = element_text(size = 24),
    axis.title.y = element_text(size = 24),
  )

# Show plot
print(p)

save_highres_plot(p, filename = "Figure5/Volcanoplot", width = 16)



library(dplyr)
library(tidyr)
library(pheatmap)

# Cleanly combine DA peaks per cluster
combined_da_peaks <- bind_rows(
  lapply(names(da_peaks_per_cluster_list), function(clust) {
    df <- da_peaks_per_cluster_list[[clust]]
    df$peak <- rownames(df)
    df$cluster <- clust
    return(df)
  })
)








## Load required libraries
library(tidyverse)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

DARs <- read.csv("Data/Pyramidal/EPO_vs_Placebo/Newly_formed_vs_Mature_DA_peaks_combined_pyramidal_celltypes_EPO_vs_Placebo_separation_filtered_annotated_for_genes_and_TEs.csv")
head(DARs)


DARs_nf <- DARs %>%
  filter(grepl("newly_formed", comparison, ignore.case = TRUE)) %>%
  filter(!is.na(summit)) %>%
  separate(summit, into = c("chr", "start", "end"), sep = "[:-]", convert = TRUE) %>%
  mutate(chr = ifelse(grepl("^chr", chr), chr, paste0("chr", chr)))

# Step 2: Convert summit to GRanges
summit_gr <- makeGRangesFromDataFrame(DARs_nf,
                                      seqnames.field = "chr",
                                      start.field = "start",
                                      end.field = "end",
                                      keep.extra.columns = TRUE)

# Step 3: Get TSS coordinates from mm10 TxDb
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
tss_gr <- promoters(transcripts(txdb), upstream = 0, downstream = 1)

# Step 4: Find distance to nearest TSS
nearest_hits <- distanceToNearest(summit_gr, tss_gr)
DARs_nf$distance_to_TSS <- NA
DARs_nf$distance_to_TSS[queryHits(nearest_hits)] <- mcols(nearest_hits)$distance

# Step 5: Categorize by distance to TSS (±100 kb window)
DARs_nf$DAR_category <- case_when(
  DARs_nf$distance_to_TSS <= 1000 ~ "Promoter-proximal (1 kb)",
  DARs_nf$distance_to_TSS <= 10000 ~ "Distal-only (1–10 kb)",
  TRUE ~ "Intergenic (>10 kb)"
)

# Set factor level order
DARs_nf$DAR_category <- factor(DARs_nf$DAR_category,
                               levels = rev(c("Promoter-proximal (1 kb)",
                                          "Distal-only (1–10 kb)",
                                          "Intergenic (>10 kb)")))

# Step 6: Count per category
summary_df <- DARs_nf %>%
  count(DAR_category)

# Step 7: Horizontal stacked bar plot
p <- ggplot(summary_df, aes(x = "DARs", y = n, fill = DAR_category)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  geom_text(aes(label = n),
            position = position_stack(vjust = 0.5),
            fontface = "bold", size = 5, color = "black") +
  scale_fill_manual(values = c(
    "Promoter-proximal (1 kb)" = "#1f78b4",
    "Distal-only (1–10 kb)" = "#b2df8a",
    "Intergenic (>10 kb)" = "#999999"
  ),  guide = guide_legend(reverse = TRUE)) +
  labs(
    title = "Genomic Distribution of DARs in Newly Formed Neurons",
    subtitle = "Summit distance to nearest TSS",
    x = NULL,
    y = "Number of DARs",
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 13),
    legend.position = "top"
  ) +
  coord_flip()

save_highres_plot(p, filename = "Figure5/Summit distance to nearest TSS (mm10) using ±10 kb classification", height = 3)








library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

# Load DARs file
DARs <- read_csv("Data/Pyramidal/EPO_vs_Placebo/Newly_formed_vs_Mature_DA_peaks_combined_pyramidal_celltypes_EPO_vs_Placebo_separation_filtered_annotated_for_genes_and_TEs.csv")
DARs$cluster <- paste0("Cluster_", DARs$cell_stage)
DARs$overlapping_promoter_old <- DARs$overlapping_promoter
DARs$overlapping_promoter <- NA
# Filter for significant DARs
DARs_filtered <- DARs %>%
  filter(p_val_adj < 0.01)

#DARs_filtered <- DARs_filtered %>%
#  filter(pct.1 >= 0.10 | pct.2 >= 0.10)

# Clean TE overlaps by removing low-confidence classes/families
exclude_classes <- c("Low_complexity", "Satellite", "rRNA", "tRNA", "scRNA", "snRNA", 
                     "RNA", "RC", "srpRNA", "DNA?", "LTR?", "LINE?", "RC?", "SINE?")

DARs_filtered <- DARs_filtered %>%
  mutate(overlapping_TE = ifelse(
    overlapping_TE & (TE_class %in% exclude_classes | grepl("\\?$", TE_family)),
    FALSE,
    overlapping_TE
  ))

head(DARs_filtered)


library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# 1. Get ±1 kb promoter regions
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoter_regions <- promoters(transcripts(txdb), upstream = 1000, downstream = 1000)

# 2. Extract peak coordinates and create GRanges for DAR summits
peak_split <- stringr::str_split_fixed(DARs_filtered$peak, "-", 3)

dar_summits <- GRanges(
  seqnames = peak_split[, 1],
  ranges = IRanges(
    start = floor((as.numeric(peak_split[, 2]) + as.numeric(peak_split[, 3])) / 2),
    end = floor((as.numeric(peak_split[, 2]) + as.numeric(peak_split[, 3])) / 2)
  )
)

# 3. Check overlaps between DAR summits and promoter regions
overlap_idx <- queryHits(findOverlaps(dar_summits, promoter_regions))
DARs_filtered$overlapping_promoter <- FALSE
DARs_filtered$overlapping_promoter[unique(overlap_idx)] <- TRUE


library(dplyr)
# 4. Preview result
#head(DARs_filtered %>% select(peak, overlapping_promoter))

# Create clean cluster labels
DARs_filtered <- DARs_filtered %>%
  mutate(Neuron_Type = case_when(
    cluster == "Cluster_newly_formed" ~ "Newly Formed",
    cluster == "Cluster_mature" ~ "Mature",
    TRUE ~ cluster
  ))

# Summarize counts
summary_counts <- DARs_filtered %>%
  group_by(Neuron_Type) %>%
  summarise(
    `Total DARs` = n(),
    `Gene-TSS DARs` = sum(overlapping_promoter),
    `TE-derived DARs` = sum(overlapping_TE),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = -Neuron_Type, names_to = "Category", values_to = "Count")

# Factor ordering
summary_counts$Category <- factor(summary_counts$Category,
                                  levels = c("Total DARs", "Gene-TSS DARs", "TE-derived DARs"))
summary_counts$Neuron_Type <- factor(summary_counts$Neuron_Type,
                                     levels = c("Newly Formed", "Mature"))

# Color palette
color_palette <- c("Newly Formed" = "#e31a1c", "Mature" = "#1f78b4")

# Plot
p <- ggplot(summary_counts, aes(x = Category, y = Count, fill = Neuron_Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, color = "black") +
  geom_text(aes(label = Count),
            position = position_dodge(width = 0.8),
            vjust = -0.5,
            size = 4,
            fontface = "bold",
            color = "black") +
  scale_fill_manual(values = color_palette) +
  labs(
    title = "EPO-Induced Chromatin Remodeling",
    subtitle = "Comparison of DAR Types by Neuron Maturation State",
    x = NULL,
    y = "Number of DARs",
    fill = "Neuron Type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 13),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 13, face = "bold"),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p

save_highres_plot(p, filename = "Figure5/Comparison of DAR Types by Neuron Maturation State", width = 6, height = 16)










# Load libraries
library(tidyverse)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ggplot2)

# 1. Load and parse DARs
DARs <- read_csv("Data/Pyramidal/EPO_vs_Placebo/Newly_formed_vs_Mature_DA_peaks_combined_pyramidal_celltypes_EPO_vs_Placebo_separation_filtered_annotated_for_genes_and_TEs.csv")
DARs$cluster <- paste0("Cluster_", DARs$cell_stage)
DARs <- DARs %>%
  separate(summit, into = c("chr", "range"), sep = ":") %>%
  separate(range, into = c("start", "end"), sep = "-") %>%
  mutate(
    start = as.integer(start),
    end = as.integer(end),
    summit_label = paste0(chr, ":", start, "-", end)
  )

dar_gr <- GRanges(
  seqnames = DARs$chr,
  ranges = IRanges(start = DARs$start, end = DARs$end),
  summit_label = DARs$summit_label
)

# 2. Build clean gene list
manual_genes <- c("Egr1", "Egr2", "Egr3", "Ascl1", "Bcl11a", "Nrg3",
                  "Dnmt3a", "Setd5", "Septin4", "Actn1", "Kalrn", "Ntrk3",
                  "Per2", "Neto1", "Cep95", "Ryr1", "Zbtb47")

dar_genes <- unique(DARs$gene_symbol)

all_genes <- unique(c(manual_genes, dar_genes))
all_genes <- all_genes[!is.na(all_genes)]
all_genes <- all_genes[!grepl("Rik$", all_genes, ignore.case = TRUE)]

# 3. Map to TSS using TxDb + OrgDb
target_entrez <- mapIds(org.Mm.eg.db, keys = all_genes,
                        column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
gene_gr <- genes(txdb)[names(genes(txdb)) %in% target_entrez]
gene_gr$symbol <- mapIds(org.Mm.eg.db,
                         keys = names(gene_gr),
                         column = "SYMBOL",
                         keytype = "ENTREZID",
                         multiVals = "first")

# 4. Create ±100 kb window
tss_gr <- resize(gene_gr, width = 1, fix = "start")
window_gr <- resize(tss_gr, width = 200001, fix = "center")

# 5. Find overlaps
hits <- findOverlaps(window_gr, dar_gr)

gene_symbols <- gene_gr$symbol[queryHits(hits)]
summits <- mcols(dar_gr)$summit_label[subjectHits(hits)]
tss_coords <- start(tss_gr)[queryHits(hits)]
summit_coords <- start(dar_gr)[subjectHits(hits)]
distances <- summit_coords - tss_coords
dar_rows <- subjectHits(hits)

# 6. Create annotation table
peak_table <- tibble(
  gene = gene_symbols,
  summit = summits,
  distance_to_TSS = distances
) %>%
  left_join(DARs %>% dplyr::select(summit_label, avg_log2FC, p_val_adj),
            by = c("summit" = "summit_label")) %>%
  mutate(
    significance = ifelse(p_val_adj < 0.05, "FDR < 0.05", "ns"),
    is_distal = abs(distance_to_TSS) > 20000,
    plot_distance = case_when(
      distance_to_TSS > 20000 ~ 22000,
      distance_to_TSS < -20000 ~ -22000,
      TRUE ~ distance_to_TSS
    ),
    distance_label = case_when(
      distance_to_TSS > 20000 ~ ">+20 kb",
      distance_to_TSS < -20000 ~ "<-20 kb",
      TRUE ~ NA_character_
    ),
    gene = factor(gene, levels = sort(unique(gene)))
  )

# 7. Plot with jitter for distal points
ggplot() +
  # background lines from TSS to peaks
  geom_segment(data = peak_table, aes(x = 0, xend = plot_distance, y = gene, yend = gene), color = "gray70", linewidth = 0.3) +
  
  # proximal peaks (≤ ±20 kb)
  geom_point(
    data = peak_table %>% filter(!is_distal),
    aes(x = plot_distance, y = gene, color = avg_log2FC, shape = significance),
    size = 3.2
  ) +
  
  # distal peaks (> ±20 kb) — jittered to separate overlapping log2FC
  geom_point(
    data = peak_table %>% filter(is_distal),
    aes(x = plot_distance, y = gene, color = avg_log2FC, shape = significance),
    position = position_jitter(width = 1000, height = 0),
    size = 3.2
  ) +
  
  # labels for distal bins
  geom_text_repel(
    data = peak_table %>% filter(is_distal),
    aes(x = plot_distance, y = gene, label = distance_label),
    size = 3,
    color = "black",
    fontface = "italic",
    direction = "y",
    nudge_y = 0.25,
    segment.size = 0.2,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = 100
  ) +
  
  scale_color_gradient2(
    low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0,
    name = expression(log[2]~"FC")
  ) +
  scale_shape_manual(values = c("FDR < 0.05" = 19, "ns" = 1), name = "Significance") +
  scale_x_continuous(
    breaks = c(-20000, -10000, 0, 10000, 20000),
    limits = c(-25000, 25000),
    labels = scales::comma_format(scale = 1)
  ) +
  labs(
    title = "DAR Summits Near TSS (±100 kb Search, ±20 kb Display)",
    subtitle = "Jittered distal points show distinct log2FC beyond ±20 kb",
    x = "Distance from TSS (bp)",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),
    axis.text.y = element_text(face = "bold", size = 11),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5)
  )

# Optional: save output
# ggsave("Final_Lollipop_DARs.pdf", width = 11, height = 7)





# Load libraries
library(tidyverse)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ggtext)
 
# 1. Load and parse DARs
DARs <- read_csv("Data/Pyramidal/EPO_vs_Placebo/Newly_formed_vs_Mature_DA_peaks_combined_pyramidal_celltypes_EPO_vs_Placebo_separation_filtered_annotated_for_genes_and_TEs.csv")
DARs$gene_symbol_old <- DARs$gene_symbol
DARs$overlapping_promoter_old <- DARs$overlapping_promoter
DARs$gene_symbol <- NA
DARs$overlapping_promoter <- NA

# 2. Parse summit coordinates into GRanges
DARs <- DARs %>%
  separate(summit, into = c("chr", "range"), sep = ":") %>%
  separate(range, into = c("start", "end"), sep = "-") %>%
  mutate(
    start = as.integer(start),
    end = as.integer(end),
    summit_label = paste0(chr, ":", start, "-", end)
  )

dar_summits <- GRanges(
  seqnames = DARs$chr,
  ranges = IRanges(
    start = floor((DARs$start + DARs$end) / 2),
    end = floor((DARs$start + DARs$end) / 2)
  )
)

# 3. Load transcript annotations and map to gene symbols
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
transcripts_gr <- transcripts(txdb)

# Map TXID → GENEID
tx2gene <- AnnotationDbi::select(
  txdb,
  keys = as.character(mcols(transcripts_gr)$tx_id),
  keytype = "TXID",
  columns = "GENEID"
)
transcripts_gr$gene_id <- tx2gene$GENEID[match(mcols(transcripts_gr)$tx_id, tx2gene$TXID)]

# Map GENEID → SYMBOL
gene_id_to_symbol <- mapIds(org.Mm.eg.db,
                            keys = transcripts_gr$gene_id,
                            column = "SYMBOL",
                            keytype = "ENTREZID",
                            multiVals = "first")
transcripts_gr$symbol <- gene_id_to_symbol

# 4. Define TSS points and ±10 kb promoter regions
tss_gr <- GRanges(
  seqnames = seqnames(transcripts_gr),
  ranges = IRanges(
    start = ifelse(strand(transcripts_gr) == "+", start(transcripts_gr), end(transcripts_gr)),
    width = 1
  ),
  strand = strand(transcripts_gr),
  gene_id = as.character(transcripts_gr$gene_id),
  symbol = as.character(transcripts_gr$symbol)
)

promoter_regions <- resize(tss_gr, width = 20001, fix = "center")  # ±10 kb

# 5. Nearest TSS → assign gene symbol
nearest_idx <- nearest(dar_summits, tss_gr)
DARs$gene_symbol <- mcols(tss_gr)$symbol[nearest_idx]

# 6. Overlap with promoter regions
overlap_hits <- queryHits(findOverlaps(dar_summits, promoter_regions))
DARs$overlapping_promoter <- FALSE
DARs$overlapping_promoter[unique(overlap_hits)] <- TRUE

# 7. Confirm output
DARs %>% dplyr::select(peak, summit_label, gene_symbol, overlapping_promoter) %>% head()


dar_gr <- GRanges(
  seqnames = DARs$chr,
  ranges = IRanges(start = DARs$start, end = DARs$end),
  summit_label = DARs$summit_label
)

# 2. Manual genes + all DAR genes
manual_genes <- c("Egr1", "Egr2", "Egr3", "Ascl1", "Bcl11a", "Nrg3",
                  "Dnmt3a", "Setd5", "Septin4", "Actn1", "Kalrn", "Ntrk3",
                  "Per2", "Neto1", "Cep95", "Ryr1", "Zbtb47")

dar_genes <- unique(DARs$gene_symbol)
all_genes <- unique(c(manual_genes, dar_genes))
all_genes <- all_genes[!is.na(all_genes)]
all_genes <- all_genes[!grepl("Rik$", all_genes, ignore.case = TRUE)]
all_genes <- all_genes[!grepl("^Gm\\d+", all_genes, ignore.case = TRUE)]


# 3. Map gene symbols to transcript-level TSS
transcripts_gr <- transcripts(txdb)

# Map TXID to GENEID
tx2gene <- AnnotationDbi::select(
  txdb,
  keys = as.character(mcols(transcripts_gr)$tx_id),
  keytype = "TXID",
  columns = "GENEID"
)
transcripts_gr$gene_id <- tx2gene$GENEID[match(mcols(transcripts_gr)$tx_id, tx2gene$TXID)]

# Map GENEID to SYMBOL
gene_id_to_symbol <- mapIds(org.Mm.eg.db,
                            keys = transcripts_gr$gene_id,
                            column = "SYMBOL",
                            keytype = "ENTREZID",
                            multiVals = "first")
transcripts_gr$symbol <- gene_id_to_symbol

# Filter only for target genes
transcripts_gr <- transcripts_gr[transcripts_gr$symbol %in% all_genes]

# 4. Define TSS and ±10 kb window
tss_gr <- GRanges(
  seqnames = seqnames(transcripts_gr),
  ranges = IRanges(
    start = ifelse(strand(transcripts_gr) == "+", start(transcripts_gr), end(transcripts_gr)),
    width = 1
  ),
  strand = as.character(strand(transcripts_gr)),
  symbol = as.character(transcripts_gr$symbol)
)

window_gr <- resize(tss_gr, width = 20001, fix = "center")  # ±10 kb

# 5. Overlap peaks with gene windows
hits <- findOverlaps(window_gr, dar_gr)

gene_symbols <- mcols(tss_gr)$symbol[queryHits(hits)]
summits <- mcols(dar_gr)$summit_label[subjectHits(hits)]
tss_coords <- start(tss_gr)[queryHits(hits)]
summit_coords <- start(dar_gr)[subjectHits(hits)]
distances <- summit_coords - tss_coords

# 6. Create base peak table with filtering
peak_table <- tibble(
  gene = gene_symbols,
  summit = summits,
  distance_to_TSS = distances
) %>%
  left_join(DARs %>% dplyr::select(summit_label, avg_log2FC, p_val_adj),
            by = c("summit" = "summit_label")) %>%
  filter(!is.na(avg_log2FC)) %>%
  group_by(gene) %>%
  mutate(min_dist = min(abs(distance_to_TSS))) %>%
  filter(abs(distance_to_TSS) == min_dist) %>%  # Keep only summit(s) with minimum distance
  summarise(
    distance_to_TSS = mean(distance_to_TSS),   # If ties, average their distance
    avg_log2FC = mean(avg_log2FC),
    p_val_adj = mean(p_val_adj),
    .groups = "drop"
  ) %>%
  mutate(
    significance = ifelse(p_val_adj < 0.01 & abs(avg_log2FC) >= 0.50,
                          "FDR < 0.01 & log2FC ≥ 0.50", "NS"),
    capped_log2FC = pmax(pmin(avg_log2FC, 2), -2),
    is_manual = gene %in% manual_genes,
    gene_label = ifelse(is_manual, paste0("<b>", gene, "</b>"), gene)
  ) %>%
  filter(significance != "NS")

# 7. Split into manual and other gene facets
manual_table <- peak_table %>%
  filter(is_manual) %>%
  arrange(gene_label) %>%
  mutate(
    gene_label = factor(gene_label, levels = rev(unique(gene_label))),
    panel = "Neurogenic, Epigenetic, Synaptic Genes"
  )

other_table <- peak_table %>%
  filter(!is_manual) %>%
  arrange(gene_label) %>%
  mutate(
    gene_label = factor(gene_label, levels = rev(unique(gene_label))),
    block = ceiling(row_number() / ceiling(nrow(.) / 8)),  # 8 blocks
    panel = paste("Other Gene Panel", block)
  ) %>%
  dplyr::select(-block)

# Combine both tables
peak_table_final <- bind_rows(manual_table, other_table)

# Plot with 3 columns
p <- ggplot(peak_table_final, aes(x = distance_to_TSS, y = gene_label)) +
  geom_segment(aes(x = 0, xend = distance_to_TSS, yend = gene_label),
               color = "gray70", linewidth = 0.3) +
  geom_point(aes(color = capped_log2FC), shape = 19, size = 3.2) +
  scale_color_gradient2(
    low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0,
    limits = c(-2, 2),
    name = expression(log[2]~"FC")
  ) +
  scale_x_continuous(
    breaks = c(-10000, -5000, 0, 5000, 10000),
    limits = c(-10000, 10000),
    labels = scales::comma_format(scale = 1)
  ) +
  facet_wrap(~ panel, scales = "free_y", ncol = 3) +  # ⬅️ 3 columns = 3x3
  labs(
    title = "Significant DAR Summits Within ±10 kb of TSS",
    subtitle = "FDR < 0.01 and |log2FC| ≥ 0.50 (EPO vs Placebo)",
    x = "Distance from TSS (bp)",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 9, face = "bold"),  # ⬅️ reduced facet label size
    panel.spacing.y = unit(2, "lines"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),
    axis.text.y = element_markdown(size = 6),
    axis.text.x = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),  # ⬅️ global title
    plot.subtitle = element_text(size = 11, hjust = 0.5)
  )


p

save_highres_plot(p, "Figure5/Significant DAR Summits Within ±10 kb of TSS", width = 16, height = 12)

 
