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

# Tabulate and sort predicted.id counts
cell_order <- names(sort(table(ATAC@meta.data$predicted.id), decreasing = TRUE))

# Reassign predicted.id as a factor with the new order
ATAC@meta.data$predicted.id <- factor(ATAC@meta.data$predicted.id, levels = cell_order)



# Load motif enrichment data
motifs <- read.csv("Data/WithoutTreatment_Separation/motif_enrichment_all_clusters_without_treatment_separation.csv")

# Load libraries
library(dplyr)
library(tidyr)
library(pheatmap)

# Load motif enrichment data
motifs <- read.csv("Data/WithoutTreatment_Separation/motif_enrichment_combined_by_cluster.csv")

# Step 1: Identify motifs that pass the threshold in at least one cluster
motif_to_keep <- motifs %>%
  filter(fold.enrichment > 1.5, p.adjust < 1e-10) %>%
  pull(motif.name) %>%
  unique()

# Step 2: Keep ALL clusters for those motifs
motifs_filt <- motifs %>%
  filter(motif.name %in% motif_to_keep) %>%
  mutate(neg_log10_padj = -log10(p.adjust))

# Step 3: Prepare matrix: motif.name × cluster
motif_matrix <- motifs_filt %>%
  select(motif.name, cluster, neg_log10_padj) %>%
  pivot_wider(names_from = cluster, values_from = neg_log10_padj, values_fill = 0)

# Step 4: Convert to matrix and fix rownames
mat <- as.data.frame(motif_matrix)
rownames(mat) <- mat$motif.name
mat <- mat[, setdiff(colnames(mat), "motif.name")]

# Remove clusters with very few data points
mat <- mat[, !colnames(mat) %in% c(
  "Cluster_Intermediate cells", 
  "Cluster_Pericytes", 
  "Cluster_Ependymal cells", 
  "Cluster_Neuroimmune cells"
)]

# Step 5: Cap high values for visual scaling
mat <- as.matrix(mat)
mat <- pmin(mat, 100)  # or 300 if desired

# Step 6: Set legend scale starting from 10
breaks <- seq(10, 100, length.out = 100)

mat[is.na(mat)] <- 0

df_long <- mat %>%
  as.data.frame() %>%
  rownames_to_column("motif") %>%
  pivot_longer(
    cols = -motif,
    names_to = "cell_type",  # directly name it here!
    values_to = "value"
  )

df_long <- df_long %>%
  arrange(motif) %>%
  mutate(motif_group = paste0("Group_", ceiling(as.numeric(factor(motif)) / 65)))

library(stringr)
df_long$cell_type <- str_remove(df_long$cell_type, "^Cluster_")
df_long$cell_type <- factor(df_long$cell_type, levels = rev(cell_order))

p <- ggplot(df_long, aes(x = motif, y = cell_type, fill = value)) +
  geom_tile(color = "white", size = 0.1) +
  scale_fill_gradient(
    low = "white",  # light orange
    high = "chocolate4", # dark brown
    name = "-log10(p.adj)"
  ) +
  facet_wrap(~motif_group, scales = "free_x", ncol = 1) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12, face = "bold"),  # <-- bold y-axis
    strip.text = element_blank(),  # hides facet labels
    panel.grid = element_blank()
  ) +
  labs(x = NULL, y = NULL)

save_highres_plot(p, "Figure2/Figure_2A_MotifEnrichment_Heatmap.png", format = "png", height = 12, width = 18)


# Step 7: Plot heatmap
p <- pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("white", "red3"))(99),
  breaks = breaks,
  main = expression("-log"[10]~"(adjusted p-value), filtered by fold.enrichment > 1.5 & p.adj < 1e-10"),
  fontsize_row = 8,
  fontsize_col = 10,
  legend = TRUE,
  legend_labels = NULL,
  legend_breaks = NULL
)

# Step 8: Save high-resolution heatmap
save_highres_plot(p, "Figure2/Figure_2AS3_MotifEnrichment_Heatmap.png", format = "png", width = 20, height = 20)
# S3/2A Complete

 
# 2B Start 
ATAC <- readRDS("../ATAC_analysis_with_individual_peaks/ATAC_with_individual_peaks_annotated_filtered_motifs_added.RDS")
DefaultAssay(ATAC) <- "peaks"
da_peaks_by_cluster <- readRDS("Data/WithoutTreatment_Separation/DA_peaks_by_cluster_without_treatment_separation_for_11levels.rds")
motif_results_by_cluster <- list()

for (clust_name in names(da_peaks_by_cluster)) {
  da_table <- da_peaks_by_cluster[[clust_name]]
   
  sig_peaks <- rownames(da_table[da_table$p_val_adj < 1e-50, ])
  sig_peaks <- sig_peaks[grepl("^chr", sig_peaks)]  # keep only reference
  
  if (length(sig_peaks) > 3) {
    message("Running motif enrichment for ", clust_name)
    
    motif_enrichment <- FindMotifs(
      object = ATAC,
      features = sig_peaks
    )
    
    motif_results_by_cluster[[clust_name]] <- motif_enrichment
  } else {
    message("Skipping ", clust_name, " (too few peaks)")
  }
}

library(dplyr)

combined_motifs <- bind_rows(
  lapply(names(motif_results_by_cluster), function(clust) {
    df <- motif_results_by_cluster[[clust]]
    df$cluster <- clust  # add cluster column
    return(df)
  })
)

write.csv(combined_motifs, "Data/WithoutTreatment_Separation/motif_enrichment_combined_by_cluster.csv", row.names = FALSE)



# Step 1: Identify motifs that are significant in at least one cluster
motif_to_keep <- combined_motifs %>%
  filter(fold.enrichment > 0, pvalue < 1) %>%  # pvalue < 1 to exclude missing/invalid entries
  pull(motif.name) %>%
  unique()

# Step 2: Keep all clusters for those motifs
combined_motifs_filtered <- combined_motifs %>%
  filter(motif.name %in% motif_to_keep) %>%
  mutate(neg_log10_pval = -log10(pvalue))

# Step 3: Pivot to wide matrix
motif_matrix <- combined_motifs_filtered %>%
  select(motif.name, cluster, neg_log10_pval) %>%
  pivot_wider(names_from = cluster, values_from = neg_log10_pval, values_fill = 0)

# Step 4: Prepare matrix with correct rownames
mat <- as.data.frame(motif_matrix)
rownames(mat) <- mat$motif.name
mat <- mat[, setdiff(colnames(mat), "motif.name")]

# Step 5: Cap extreme values
mat <- pmin(mat, 300)

# Step 6: Set color scale breaks
breaks <- seq(10, 300, length.out = 100)

# Step 7: Replace any NA with 0
mat[is.na(mat)] <- 0

colnames(mat) <- gsub("Cluster_", "", colnames(mat))

library(pheatmap)
# Step 8: Generate heatmap without row labels
p <- pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("white", "darkred"))(99),
  breaks = breaks,
  main = expression("-log"[10]~"(p-value), enriched motifs per cluster"),
  fontsize_row = 8,
  fontsize_col = 10,
  border_color = NA,
  labels_row = NA
)

# Step 9: Save high-resolution image
save_highres_plot(p, "Figure2/Figure_Go_to_Supplementary_2B_MotifEnrichment_Heatmap.png", format = "png")


# Step 1: Find motifs that pass threshold in at least one cluster
motif_to_keep <- combined_motifs %>%
  filter(fold.enrichment > 1.5, p.adjust < 1e-200) %>%
  pull(motif.name) %>%
  unique()

# Step 2: For those motifs, keep all p.adjust values across all clusters
combined_motifs_filtered <- combined_motifs %>%
  filter(motif.name %in% motif_to_keep) %>%
  mutate(neg_log10_padj = -log10(p.adjust))

# Step 3: Pivot to wide format: motif.name × cluster
motif_matrix <- combined_motifs_filtered %>%
  select(motif.name, cluster, neg_log10_padj) %>%
  pivot_wider(names_from = cluster, values_from = neg_log10_padj, values_fill = 0)

# Step 4: Prepare matrix and fix rownames
mat <- as.data.frame(motif_matrix)
rownames(mat) <- mat$motif.name
mat <- mat[, setdiff(colnames(mat), "motif.name")]

# Step 5: Cap large values (e.g., >300)
mat <- pmin(mat, 300)

# Step 6: Define color scale
breaks <- seq(10, 300, length.out = 100)

# Step 7: Replace any NA with 0 (for safety)
mat[is.na(mat)] <- 0

colnames(mat) <- gsub("Cluster_", "", colnames(mat))

# Step 8: Plot heatmap
p <- pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("white", "darkred"))(99),
  breaks = breaks,
  main = expression("-log"[10]~"(adjusted p-value), enriched motifs per cluster, filtered by fold.enrichment > 1.5, p.adjust < 1e-200"),
  fontsize_row = 8,
  fontsize_col = 10,
  border_color = NA
)

# Step 9: Save high-resolution image
save_highres_plot(p, "Figure1/Figure_2B_MotifEnrichment_Heatmap.png", format = "png", width = 16, height = 16)
 
# 2B Complete 

# 2C Start

# Load libraries
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPseeker)
library(dplyr)
library(ggplot2)
library(tibble)
library(forcats)

# Load TxDb
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Load ATAC peaks
ATAC <- readRDS("../ATAC_analysis_with_individual_peaks/ATAC_with_individual_peaks_annotated_filtered_motifs_added.RDS")
atac_gr <- granges(ATAC)
standard_chroms <- standardChromosomes(atac_gr)
atac_gr <- atac_gr[seqnames(atac_gr) %in% standard_chroms]
peak <- resize(atac_gr, width = 1, fix = "center")

# Load and preprocess TE annotations
te_df <- read.delim("Data/mm10_RepeatMasker/mm10_RepeatMasker.tsv", header = TRUE)
exclude_classes <- c("Low_complexity", "Satellite", "rRNA", "tRNA", "scRNA", "snRNA", 
                     "RNA", "RC", "srpRNA", "DNA?", "LTR?", "LINE?", "RC?", "SINE?")
te_df <- subset(te_df, !(repClass %in% exclude_classes))
# Remove families with "?" in the name
te_df <- te_df[!grepl("\\?", te_df$repFamily), ]

te_gr <- GRanges(
  seqnames = te_df$genoName,
  ranges = IRanges(start = te_df$genoStart + 1, end = te_df$genoEnd),
  strand = "*"
)
te_gr <- te_gr[seqnames(te_gr) %in% standardChromosomes(te_gr)]
te_gr <- te_gr[seqnames(te_gr) %in% standardChromosomes(te_gr)]
te_gr_center <- resize(te_gr, width = 1, fix = "center")


# Annotate with ChIPseeker
peak_anno <- annotatePeak(peak, TxDb = txdb, verbose = TRUE)
te_anno <- annotatePeak(te_gr_center, TxDb = txdb, verbose = TRUE)
overlap_idx <- which(countOverlaps(peak, te_gr) > 0)
peak_in_te_anno <- annotatePeak(peak[overlap_idx], TxDb = txdb, verbose = TRUE)

# Annotation summaries from ChIPseeker
df_all_peaks <- as.data.frame(peak_anno@annoStat) %>% 
  mutate(group = "all_peaks")

df_all_te <- as.data.frame(te_anno@annoStat) %>%
  mutate(group = "all_TE")

df_peaks_in_te <- as.data.frame(peak_in_te_anno@annoStat) %>%
  mutate(group = "peaks_in_TE")

# Combine
anno_df <- bind_rows(df_all_te, df_all_peaks, df_peaks_in_te)

# Rename for clarity
colnames(anno_df) <- c("annotation", "percent", "group")

# Get order of annotations based on decreasing percent in the 'all_TE' group
anno_order <- anno_df %>%
  filter(group == "all_TE") %>%
  arrange(desc(percent)) %>%
  pull(annotation)

# Apply the new annotation order as factor levels
anno_df$annotation <- factor(anno_df$annotation, levels = rev(anno_order))

# Color scheme (update as needed for additional categories like Promoter (2-3kb))
annotation_colors <- c(
  "Promoter (<=1kb)" = "#de2d26",
  "Promoter (1-2kb)" = "#fdae6b",
  "Promoter (2-3kb)" = "#fee6ce",
  "5' UTR" = "#a1d99b",
  "3' UTR" = "#238b45",
  "1st Exon" = "#9e79b3",
  "Other Exon" = "#d6a0cd",
  "1st Intron" = "#fdd9a0",
  "Other Intron" = "#a67c52",
  "Downstream (<=300)" = "#9ecae1",
  "Distal Intergenic" = "#3182bd"
)

# Rename group labels for display
anno_df$group <- recode(anno_df$group,
                        all_TE = "All TE",
                        all_peaks = "All Peaks",
                        peaks_in_TE = "Peaks in TE"
)

# Set the desired display order
anno_df$group <- factor(anno_df$group, levels = c("All TE", "All Peaks", "Peaks in TE"))

# Plot
library(ggplot2)
library(forcats)

p <- ggplot(anno_df, aes(x = fct_rev(group), y = percent, fill = annotation)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = annotation_colors, drop = FALSE) +
  coord_flip() +
  labs(y = "Percentage (%)", x = NULL, fill = NULL) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.5, "cm"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(fill = guide_legend(reverse = TRUE, nrow = 2, byrow = TRUE))

save_highres_plot(p, filename = "Figure2/peak_summits_annotation", height = 5)



readRDS("../ATAC_analysis_with_individual_peaks/ATAC_with_individual_peaks_annotated_filtered_motifs_added.RDS") -> alldata

granges(alldata[["peaks"]]) -> all_peaks 

summits <- GRanges(
  seqnames = seqnames(all_peaks),
  ranges = IRanges(
    start = start(all_peaks) + floor(width(all_peaks) / 2),
    width = 1
  )
)

# Keep only chromosomes starting with "chr"
summits <- summits[grepl("^chr", seqnames(summits))]

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes <- genes(txdb)

te_df <- read.delim("Data/mm10_RepeatMasker/mm10_RepeatMasker.tsv", header = TRUE)
exclude_classes <- c("Low_complexity", "Satellite", "rRNA", "tRNA", "scRNA", "snRNA", 
                     "RNA", "RC", "srpRNA", "DNA?", "LTR?", "LINE?", "RC?", "SINE?")
 
# Exclude repClass AND repFamily if they match any in exclude_classes
te_df <- te_df %>%
  filter(
    !repClass %in% exclude_classes,
    !repFamily %in% exclude_classes,
    !grepl("\\?", repFamily)
  )

te_bed <- GRanges(
  seqnames = te_df$genoName,
  ranges = IRanges(start = te_df$genoStart + 1, end = te_df$genoEnd),  # BED 0-based
  strand = te_df$strand,
  TE_name = te_df$repName,
  TE_class = te_df$repClass,
  TE_family = te_df$repFamily
)


hits_genes <- findOverlaps(summits, genes)
hits_tes <- findOverlaps(summits, te_bed)

# Assign overlap status
summits$overlap <- "None"
summits$overlap[queryHits(hits_genes)] <- "Gene"
summits$overlap[queryHits(hits_tes)] <- ifelse(
  summits$overlap[queryHits(hits_tes)] == "Gene", "Both", "TE")

library(ChIPseeker)
library(org.Mm.eg.db)
peak_anno <- annotatePeak(summits, TxDb = txdb)

# Add region annotation
summits$gene_region <- peak_anno@anno$annotation

te_hits <- findOverlaps(summits, te_bed)
summits$TE_class <- NA
summits$TE_family <- NA
summits$TE_name <- NA

summits$TE_class[queryHits(te_hits)] <- mcols(te_bed)$TE_class[subjectHits(te_hits)]
summits$TE_family[queryHits(te_hits)] <- mcols(te_bed)$TE_family[subjectHits(te_hits)]
summits$TE_name[queryHits(te_hits)] <- mcols(te_bed)$TE_name[subjectHits(te_hits)]


# Reorder
overlap_table <- table(summits$overlap)
desired_order <- c("Gene", "TE", "Both", "None")
overlap_table <- overlap_table[desired_order]

# Convert table to dataframe for ggplot
df_overlap <- as.data.frame(overlap_table)
colnames(df_overlap) <- c("Category", "Count")


# Load libraries
library(eulerr)
library(grid)

# Extract counts
gene_only <- df_overlap$Count[df_overlap$Category == "Gene"]
te_only   <- df_overlap$Count[df_overlap$Category == "TE"]
both      <- df_overlap$Count[df_overlap$Category == "Both"]
none      <- df_overlap$Count[df_overlap$Category == "None"]

# Create euler object
fit <- euler(c(
  Gene = gene_only,
  TE = te_only,
  "Gene&TE" = both
))

# Total number
total <- gene_only + te_only + both + none

# Open PNG device
png("Figure2/Venn_diagram.png", width = 6, height = 6, units = "in", res = 1200)

# === Your full plotting code below ===
# Plot Venn
plot(
  fit,
  fills = list(fill = c("skyblue", "orange"), alpha = 0.7),
  labels = list(font = 2, cex = 1.6),
  quantities = list(cex = 1.8, font = 2, formatter = function(x) format(x, big.mark = ",")),
  edges = TRUE
)

# Draw total frame
grid.rect(
  x = 0.5, y = 0.5, width = 1.05, height = 1.05,
  gp = gpar(col = "black", fill = NA, lwd = 2)
)

# Add formatted "Total" label
grid.text(
  label = paste0("Total: ", format(total, big.mark = ",")),
  x = unit(0.02, "npc"), y = unit(0.97, "npc"),
  just = "left", gp = gpar(fontsize = 16, fontface = "bold")
)

# === End of plotting ===

# Close PNG device
dev.off()



library(dplyr)
library(ggplot2)
library(forcats)

# Subset summits overlapping TE or Both
te_peaks <- summits[summits$overlap %in% c("TE", "Both")]
te_df <- as.data.frame(te_peaks)

# Count TE_family × TE_class
te_counts <- te_df %>%
  count(TE_class, TE_family, name = "Count") %>%
  filter(!is.na(TE_class)) %>%
  arrange(desc(Count))

# Order TE_class by total abundance
class_order <- te_counts %>%
  group_by(TE_class) %>%
  summarise(total = sum(Count)) %>%
  arrange(desc(total)) %>%
  pull(TE_class)

te_counts$TE_class <- factor(te_counts$TE_class, levels = class_order)

# Define color palette
te_colors <- c(
  "DNA" = "#E41A1C",
  "LINE" = "#FF7F00",
  "LTR" = "#4DAF4A",
  "SINE" = "#377EB8",
  "Simple_repeat" = "#984EA3",
  "Other" = "#A65628",
  "Unknown" = "#999999"
)[class_order]

# Create legend labels with counts
class_counts <- te_counts %>%
  group_by(TE_class) %>%
  summarise(total = sum(Count)) %>%
  arrange(match(TE_class, class_order))

legend_labels <- paste0(class_counts$TE_class, " (", format(class_counts$total, big.mark = ",", trim = TRUE), ")")
names(legend_labels) <- class_counts$TE_class  # ✅ important to name the vector

# Plot
p <- ggplot(te_counts, aes(x = fct_reorder(TE_family, Count), y = Count, fill = TE_class)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = te_colors,
    labels = legend_labels,
    drop = FALSE,
    name = "TE Class"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.025))) +
  labs(
    title = "Distribution of TE Families in ATAC Peak Summits",
    x = "TE Family",
    y = "Number of Overlapping Summits"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_text(face = "bold", size = 18, margin = margin(r = 2)),
    axis.text.x = element_text(size = 18, face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 18),
    legend.position = "right",
    panel.grid.major.y = element_blank()
  )

save_highres_plot(p, filename = "Figure2/Distribution_of_TE_Families_in_ATAC_Peak_Summits", width = 20, height = 9)


library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPseeker)
library(org.Mm.eg.db)
library(AnnotationDbi) 

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Step 1: Annotate summits using ChIPseeker
peak_anno <- annotatePeak(summits, TxDb = txdb, annoDb = "org.Mm.eg.db")

# Step 2: Extract annotation dataframe
anno_df <- as.data.frame(peak_anno)

# Map Entrez gene ID to gene type and gene symbol
anno_df$geneSymbol <- mapIds(
  org.Mm.eg.db,
  keys = anno_df$geneId,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

anno_df$geneType <- mapIds(
  org.Mm.eg.db,
  keys = anno_df$geneId,
  column = "GENETYPE",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Match by rownames — ensure 1:1 mapping
stopifnot(length(summits) == nrow(anno_df))

mcols(summits)$geneId     <- anno_df$geneId
mcols(summits)$geneSymbol <- anno_df$geneSymbol
mcols(summits)$geneType   <- anno_df$geneType

# Load libraries
library(ggplot2)
library(ggbreak)
library(dplyr)
library(forcats)
library(RColorBrewer)

# Prepare data
gene_df <- summits[summits$overlap %in% c("Gene", "Both")] %>%
  as.data.frame() %>%
  mutate(geneType = ifelse(is.na(geneType), "unknown", geneType))

# Count and sort gene types
gene_counts <- gene_df %>%
  count(geneType, name = "Count") %>%
  arrange(desc(Count))

# Reorder geneType as factor for both fill and legend
gene_counts$geneType <- factor(gene_counts$geneType, levels = gene_counts$geneType)

# Legend labels with formatted counts
legend_labels <- paste0(levels(gene_counts$geneType), " (", format(gene_counts$Count, big.mark = ",", trim = TRUE), ")")
names(legend_labels) <- levels(gene_counts$geneType)

# Color palette
gene_colors <- RColorBrewer::brewer.pal(min(nrow(gene_counts), 8), "Dark2")
names(gene_colors) <- levels(gene_counts$geneType)

# Plot
p <- ggplot(gene_counts, aes(x = fct_reorder(geneType, Count), y = Count, fill = geneType)) +
  geom_col(color = "black", width = 0.8) +
  coord_flip(clip = "off") +
  scale_y_break(c(15000, 110000), scales = 0.5) +
  scale_fill_manual(
    values = gene_colors,
    name = "Gene Type",
    labels = legend_labels,
    breaks = levels(gene_counts$geneType)  # ensure correct order
  ) +
  labs(
    title = "Distribution of Gene Types in ATAC Peak Summits",
    x = "Gene Type",
    y = "Number of Overlapping Summits"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_text(face = "bold", size = 16),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 14),
    legend.position = "right",
    panel.grid.major.y = element_blank()
  )

# Show plot
p



save_highres_plot(p, "Figure2/gene_type_distribution.png", width = 25, height = 6)

 

library(regioneR)
library(GenomicRanges)
library(BiocParallel)
#register(MulticoreParam(workers = 12))  # adjust to your CPU count
set.seed(0)

# Load ATAC peaks
ATAC <- readRDS("../ATAC_analysis_with_individual_peaks/ATAC_with_individual_peaks_annotated_filtered_motifs_added.RDS")
atac_gr <- granges(ATAC)
standard_chroms <- standardChromosomes(atac_gr)
atac_gr <- atac_gr[seqnames(atac_gr) %in% standard_chroms]
peak <- resize(atac_gr, width = 1, fix = "center")

# Load RepeatMasker with header
te_df <- read.delim("Data/mm10_RepeatMasker/mm10_RepeatMasker.tsv", header = TRUE, stringsAsFactors = FALSE)

# Exclude unwanted TE classes
exclude_classes <- c("Low_complexity", "Satellite", "rRNA", "tRNA", "scRNA", "snRNA",
                     "RNA", "RC", "srpRNA", "DNA?", "LTR?", "LINE?", "RC?", "SINE?")
te_df <- subset(te_df, !(repClass %in% exclude_classes))

# Keep only standard chromosomes (chr1–19, X, Y)
te_df <- te_df[te_df$genoName %in% paste0("chr", c(1:19, "X", "Y")), ]

# Convert to GRanges
TEs <- GRanges(
  seqnames = te_df$genoName,
  ranges = IRanges(start = te_df$genoStart + 1, end = te_df$genoEnd),
  strand = te_df$strand,
  repName = te_df$repName,
  repClass = te_df$repClass,
  repFamily = te_df$repFamily
)

peaks_gr <- peak
TEs_chr <- TEs

# Unique TE names
repNames <- unique(mcols(TEs)$repName)
tf_df <- data.frame(
TE_name = repNames,
Pvalue = NA_real_,
Observed = NA_integer_,
Expected = NA_real_,
OE_ratio = NA_real_,
stringsAsFactors = FALSE
)

n_total <- length(repNames)
set.seed(0)
for (i in seq_along(repNames)) {
rn <- repNames[i]

cat(sprintf("Analyzing %s (%d of %d)\n", rn, i, n_total))  # progress message
fam_gr <- TEs_chr[mcols(TEs_chr)$repName == rn]

pt <- permTest(
  A = fam_gr,
  B = peaks_gr,
  ntimes = 10,
  randomize.function = resampleRegions,
  universe = TEs_chr,
  evaluate.function = numOverlaps,
  verbose = FALSE,
  alternative = "greater",
  force.parallel = TRUE  
)

obs <- pt$numOverlaps$observed
exp <- mean(pt$numOverlaps$permuted)

tf_df$Pvalue[i]   <- pt$numOverlaps$pval
tf_df$Observed[i] <- obs
tf_df$Expected[i] <- round(exp, 2)
tf_df$OE_ratio[i] <- ifelse(!is.na(exp) && exp > 0, round(obs / exp, 3), NA)
}

 

# Save final merged result
write.table(tf_df, "Figure2/Combined_TE_enrichment_by_repName.txt", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(tf_df, "Figure2/Combined_TE_enrichment_by_repName.rds")

cat("✅ Enrichment testing using repName (with O/E ratio) is complete.\n")
 

library(dplyr)
library(ggplot2)
library(forcats)

# Step 1: Map repName to repClass
rep_info <- te_df %>%
  distinct(repName, repClass)

# Step 2: Merge with enrichment results
plot_df <- tf_df %>%
  left_join(rep_info, by = c("TE_name" = "repName")) %>%
  mutate(
    repClass = ifelse(is.na(repClass), "Unknown", repClass),
    neg_log10_pval = ifelse(!is.na(Pvalue) & Pvalue > 0, -log10(Pvalue), 0),
    OE_ratio = ifelse(is.na(OE_ratio), 0, OE_ratio)
  )

# Step 3: Define desired TE class order
class_order <- c("LTR", "SINE", "LINE", "Simple_repeat", "DNA", "Other", "Unknown")
plot_df$repClass <- factor(plot_df$repClass, levels = class_order)

# Step 4: Define color palette in the correct order
te_colors <- c(
  "LTR"           = "#4DAF4A",
  "SINE"          = "#377EB8",
  "LINE"          = "#FF7F00",
  "Simple_repeat" = "#984EA3",
  "DNA"           = "#E41A1C",
  "Other"         = "#A65628",
  "Unknown"       = "#999999"
)

# Step 4: Select labeled TEs
top_labels <- plot_df %>%
  filter(Pvalue < 1e-10, OE_ratio > 15) %>%
  slice_max(order_by = neg_log10_pval, n = 10)

# Step 5: Plot
p <- ggplot(plot_df, aes(x = OE_ratio, y = neg_log10_pval, color = repClass)) +
  geom_point(alpha = 0.4, size = 1.75) +
  geom_text_repel(
    data = top_labels,
    aes(label = TE_name),
    size = 2.8,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.size = 0.2,
    max.overlaps = 100,
    min.segment.length = 0
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = 2.85, linetype = "dashed", color = "grey50") +
  scale_x_log10(
    limits = c(0.5, NA),
    breaks = c(0.5, 1, 2, 5, 10, 20, 50, 100),
    name = "Observed / Expected (log scale)"
  ) +
  scale_color_manual(values = te_colors, name = "TE Class") +
  labs(
    y = expression(-log[10](p-value)),
    title = "Transposon Enrichment in ATAC Peaks"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 14)
  )

p

save_highres_plot(p, filename = "Figure2/Transposon_enrichment_1_v2", width = 12)


# Define significance threshold
signif_threshold <- 3.0

# Filter significant TEs
signif_counts_df <- plot_df %>%
  filter(neg_log10_pval >= signif_threshold)

# Identify TEs to label (OE ratio > 5)
label_df <- signif_counts_df %>%
  filter(OE_ratio > 5)

# TE class order and colors
class_order <- c("LTR", "SINE", "LINE", "Simple_repeat", "DNA", "Other", "Unknown")
signif_counts_df$repClass <- factor(signif_counts_df$repClass, levels = class_order)
label_df$repClass <- factor(label_df$repClass, levels = class_order)


# Plot
p2 <- ggplot() +
  # Faded background points
  geom_point(
    data = signif_counts_df,
    aes(x = Expected, y = Observed),
    color = "grey75", size = 1.2, alpha = 0.3
  ) +
  # Highlighted points
  geom_point(
    data = label_df,
    aes(x = Expected, y = Observed, color = repClass),
    size = 2.2, alpha = 0.9
  ) +
  # Labels
  geom_text_repel(
    data = label_df,
    aes(x = Expected, y = Observed, label = TE_name, color = repClass),
    size = 2.8,
    box.padding = 0.3,
    point.padding = 0.25,
    segment.size = 0.3,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  # Diagonal line (Observed = Expected)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  # Scales
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = te_colors, name = "TE Class") +
  labs(
    title = "Observed vs Expected Counts for Significant Transposons",
    x = "Expected Count (log10)",
    y = "Observed Count (log10)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 14)
  )

p2


save_highres_plot(p2, filename = "Figure2/Transposon_enrichment_2", width = 12)
 
p2 <- ggplot() +
  # Faded background points
  geom_point(
    data = signif_counts_df,
    aes(x = Expected, y = Observed),
    color = "grey75", size = 1.8, alpha = 0.3
  ) +
  # Highlighted points
  geom_point(
    data = label_df,
    aes(x = Expected, y = Observed, color = repClass),
    size = 3.0, alpha = 0.9
  ) +
  # Diagonal line (Observed = Expected)
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  # Scales
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = te_colors, name = "TE Class") +
  labs(
    title = "Observed vs Expected Counts for Significant Transposons",
    x = "Expected Count (log10)",
    y = "Observed Count (log10)"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 14)
  )

p2

save_highres_plot(p2, filename = "Figure2/Transposon_enrichment_2_version2", width = 12)



# Load required libraries
library(tidyr)
library(dplyr)

# Load the combined enrichment results
combined_df <- readRDS("Combined_TE_enrichment_by_repName.rds")

# Pivot into matrix form: rows = TE_name, columns = TF, values = Pvalue
pval_matrix <- combined_df %>%
  select(TE_name, TF, Pvalue) %>%
  pivot_wider(names_from = TF, values_from = Pvalue)

# Convert to data.frame and keep TE_name as first column (no rownames)
pval_matrix_df <- as.data.frame(pval_matrix)

# Save the table with TE_name as first column (no row names)
write.table(pval_matrix_df, "TE_pvalue_matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(pval_matrix_df, "TE_pvalue_matrix.rds")

cat("✅ P-value matrix saved as properly formatted text and RDS.\n")





# Load required libraries
library(rtracklayer)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPseeker)
library(ggplot2)
library(dplyr)
library(forcats)

# Set TxDb object
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Load ATAC object
ATAC <- readRDS("../ATAC_analysis_with_individual_peaks/ATAC_with_individual_peaks_annotated_filtered.RDS")

# Load and combine cCRE BED files
bed_dir <- "Data/cCREs_mousebrain_catlasv1/"
bed_files <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE)
gr_list <- lapply(bed_files, function(file) {
  message("Importing: ", file)
  import(file, format = "BED")
})
combined_gr <- do.call(c, gr_list)
combined_gr_unique <- unique(combined_gr)

# Prepare ATAC peaks
atac_gr <- granges(ATAC)
standard_chroms <- standardChromosomes(atac_gr)
atac_gr <- atac_gr[seqnames(atac_gr) %in% standard_chroms]

# Find overlaps
overlaps <- findOverlaps(atac_gr, combined_gr_unique, ignore.strand = TRUE)
overlapping_idx <- unique(queryHits(overlaps))

# Split peaks
atac_overlapping <- atac_gr[overlapping_idx]
atac_nonoverlapping <- atac_gr[-overlapping_idx]

# Annotate with ChIPseeker
anno_overlapping <- annotatePeak(atac_overlapping, TxDb = txdb, verbose = FALSE)
anno_nonoverlapping <- annotatePeak(atac_nonoverlapping, TxDb = txdb, verbose = FALSE)

# Extract annotation summaries
df_overlapping <- as.data.frame(anno_overlapping@annoStat) %>%
  mutate(group = "Overlapping")

df_nonoverlapping <- as.data.frame(anno_nonoverlapping@annoStat) %>%
  mutate(group = "Not Overlapping")

# Combine and clean
df_anno <- bind_rows(df_overlapping, df_nonoverlapping)
colnames(df_anno) <- c("annotation", "percent", "group")

# Set annotation order by descending abundance in "Overlapping"
anno_order <- df_anno %>%
  filter(group == "Overlapping") %>%
  arrange(desc(percent)) %>%
  pull(annotation)

df_anno$annotation <- factor(df_anno$annotation, levels = anno_order)
df_anno$group <- factor(df_anno$group, levels = c("Overlapping", "Not Overlapping"))

# Plot
p1 <- ggplot(df_anno, aes(x = annotation, y = percent, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("Overlapping" = "black", "Not Overlapping" = "gray80")) +
  labs(
    title = "Genomic Annotation of ATAC Peaks by cCRE Overlap",
    x = "Genomic Annotation",
    y = "Percentage of Peaks",
    fill = "Group"  # Legend title
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 16),
    axis.text.y = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 16),  # Increase legend title
    legend.text = element_text(size = 14),
    legend.position = "top"
  )


# Count total overlapping and non-overlapping peaks
total_peaks <- length(atac_gr)
num_overlapping <- length(atac_overlapping)
num_nonoverlapping <- length(atac_nonoverlapping)

df_overlap_counts <- data.frame(
  Category = c("Overlapping", "Not Overlapping"),
  Count = c(num_overlapping, num_nonoverlapping)
)

# Plot: Horizontal bar chart
p1.1 <- ggplot(df_overlap_counts, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  geom_text(
    aes(label = format(Count, big.mark = ",")),
    hjust = -0.1, fontface = "bold", size = 5
  ) +
  scale_fill_manual(values = c("Overlapping" = "black", "Not Overlapping" = "gray80")) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Number of ATAC Peaks",
    title = "ATAC Peaks Overlapping with Mouse Brain cCREs"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(face = "bold", size = 16),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold", size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "none"
  ) +
  expand_limits(y = max(df_overlap_counts$Count) * 1.15)


save_highres_plot(p1, filename = "Figure2/ATAC_cCREs_overlap_1", width = 16)
save_highres_plot(p1.1, filename = "Figure2/ATAC_cCREs_overlap_2", height = 4)

