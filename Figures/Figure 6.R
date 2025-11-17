# Load libraries
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(readr)
library(dplyr)
library(stringr)

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


# Load the data
df <- read_csv("Data/Pyramidal/EPO_vs_Placebo/DA_peaks_combined_pyramidal_mature_vs_newly_formed_celltypes_EPO_vs_Placebo_separation_filtered_annotated_for_genes_and_TEs.csv")
head(df)
library(dplyr)
df <- df %>%
  mutate(
    cluster = case_when(
      comparison == "newly_formed_EPO_vs_newly_formed_Placebo" ~ "Cluster_newly_formed",
      comparison == "mature_EPO_vs_mature_Placebo" ~ "Cluster_mature",
      TRUE ~ NA_character_
    )
  )
dim(df)
library(dplyr)

std_chroms <- paste0("chr", c(1:19, "X", "Y", "M"))

df <- df %>%
  mutate(chrom = sub("[-:].*$", "", peak)) %>%  # get chr from 'peak'
  filter(chrom %in% std_chroms) %>%
  select(-chrom)  # remove helper column if not needed

dim(df)
# Define exclusion list
exclude_classes <- c("Low_complexity", "Satellite", "rRNA", "tRNA", "scRNA", "snRNA",  
                     "RNA", "RC", "srpRNA", "DNA?", "LTR?", "LINE?", "RC?", "SINE?")

# Modify 'overlapping_TE' in-place without removing rows
df <- df %>%
  mutate(overlapping_TE = ifelse(
    TE_class %in% exclude_classes |
      TE_family %in% exclude_classes |
      grepl("\\?", TE_family),
    FALSE,
    overlapping_TE
  ))

# Extract summit position (format: chr:start-end → we use 'start')
df <- df %>%
  mutate(
    chrom = str_extract(summit, "^chr[^:]+"),
    summit_pos = as.integer(str_extract(summit, "(?<=:)[0-9]+"))
  )

# Create GRanges object for all summits
summit_gr <- GRanges(
  seqnames = df$chrom,
  ranges = IRanges(start = df$summit_pos, end = df$summit_pos),
  strand = "*"
)

# Get ±10 kb promoter regions from mm10 gene TSSs
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoter_regions <- promoters(transcripts(txdb), upstream = 10000, downstream = 10000)

# Find overlaps
hits <- findOverlaps(summit_gr, promoter_regions)

# Annotate full data with summit proximity
df$within_10kb_promoter_summit <- FALSE
df$within_10kb_promoter_summit[queryHits(hits)] <- TRUE


library(AnnotationDbi)
library(org.Mm.eg.db)  # For mapping gene IDs to symbols

# Extract matching rows from both df and promoter_regions
summit_df_hits <- df[queryHits(hits), ]
promoter_hits <- promoter_regions[subjectHits(hits)]

# Get gene IDs from promoter regions
gene_ids <- mcols(promoter_hits)$tx_id

# Map transcript ID → gene ID
tx2gene <- AnnotationDbi::select(txdb, keys = as.character(gene_ids), columns = "GENEID", keytype = "TXID")

# Add gene symbols using org.Mm.eg.db
gene_info <- AnnotationDbi::select(org.Mm.eg.db,
                    keys = unique(tx2gene$GENEID),
                    columns = c("SYMBOL"),
                    keytype = "ENTREZID")

# Merge mappings back
tx2gene <- left_join(tx2gene, gene_info, by = c("GENEID" = "ENTREZID"))

# Associate each summit with its mapped gene(s)
summit_gene_table <- tibble(
  row_index = queryHits(hits),
  summit_id = df$peak[queryHits(hits)],
  tx_id = as.character(mcols(promoter_hits)$tx_id)  # <— FIX here
) %>%
  left_join(tx2gene, by = c("tx_id" = "TXID")) %>%
  group_by(row_index, summit_id) %>%
  summarise(promoter_genes = paste(unique(SYMBOL), collapse = ", "), .groups = "drop")


# Add this back into the main dataframe
df$promoter_gene_symbols <- NA
df$promoter_gene_symbols[summit_gene_table$row_index] <- summit_gene_table$promoter_genes



View(df)
library(tidyr)
fig6B_df <- df %>%
  filter(overlapping_TE == TRUE) %>%
  group_by(cluster) %>%
  summarise(
    Proximal = sum(within_10kb_promoter_summit),
    Distal = n() - Proximal,
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(Proximal, Distal),
               names_to = "Proximity", values_to = "Count") %>%
  filter(cluster %in% c("Cluster_newly_formed", "Cluster_mature")) %>%
  mutate(
    cluster = recode(cluster,
                     "Cluster_newly_formed" = "Newly Formed Neurons",
                     "Cluster_mature" = "Mature Neurons"),
    Proximity = recode(Proximity,
                       "Proximal" = "Within ±10kb",
                       "Distal" = "Distal (>10kb)")
  )

p <- ggplot(fig6B_df, aes(x = cluster, y = Count, fill = Proximity)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("Within ±10kb" = "#1f78b4", "Distal (>10kb)" = "#999999")) +
  theme_minimal(base_size = 16) +
  labs(
    title = "Proximity of TE-DARs to Promoters",
    x = NULL, y = "Number of TE-Associated DARs"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title.y = element_text(face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.title = element_blank(),
    legend.position = "top",
    panel.grid = element_blank()
  )

p <- p + geom_text(aes(label = Count),
              position = position_stack(vjust = 0.5),
              size = 5,
              color = "white",
              fontface = "bold")

p

save_highres_plot(p, "Figure6/Proximity of TE-DARs to Promoters", height = 3)


library(dplyr)
library(ggplot2)

# Define TE class colors
te_colors <- c(
  "DNA"           = "#E41A1C",
  "LINE"          = "#FF7F00",
  "LTR"           = "#4DAF4A",
  "SINE"          = "#377EB8",
  "Simple_repeat" = "#984EA3",
  "Other"         = "#A65628",
  "Unknown"       = "#999999"
)

# Filter TE-only and reassign matching family/class to Unknown
te_only <- df %>%
  filter(overlapping_TE == TRUE) %>%
  mutate(
    same_fc = TE_family == TE_class & !TE_class %in% c("Simple_repeat", "Other"),
    TE_class  = if_else(same_fc, "Unknown", TE_class),
    TE_family = if_else(same_fc, "Unknown", TE_family)
  ) %>%
  select(-same_fc)

# Count TE_family × TE_class
te_fam_class_counts <- te_only %>%
  count(TE_family, TE_class, sort = TRUE)

# Combine counts for TE_class = "Unknown" if multiple families overlap
class_counts <- te_fam_class_counts %>%
  group_by(TE_class) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  mutate(
    legend_label = paste0(TE_class, " (", total, ")"),
    legend_label = factor(legend_label, levels = legend_label)
  )

# Merge class info
te_fam_class_counts <- te_fam_class_counts %>%
  left_join(class_counts, by = "TE_class")

# Order TE_family levels for plotting
te_fam_class_counts <- te_fam_class_counts %>%
  mutate(TE_family = factor(TE_family, levels = unique(TE_family[order(n)])))

# Set fill color mapping (named by legend label in correct order)
color_map <- setNames(te_colors[class_counts$TE_class], class_counts$legend_label)

# Final plot
p <- ggplot(te_fam_class_counts, aes(x = TE_family, y = n, fill = legend_label)) +
  geom_col(color = "black") +
  geom_text(aes(label = n), hjust = -0.25, size = 3.2, fontface = "bold", color = "black") +
  coord_flip() +
  scale_fill_manual(
    name = "TE Class (Count)",
    values = color_map
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribution of TE Families in Overlapping TE-DARs",
    x = "TE Family",
    y = "Count"
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5, color = "black"),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.position = "right",
    panel.grid = element_blank()
  )  +
  expand_limits(y = max(te_fam_class_counts$n) * 1.1)

p

save_highres_plot(p, filename = "Figure6/Distribution of TE Families in Overlapping TE-DARs", height = 6)


# Load libraries
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(dplyr)
library(stringr)
library(ggplot2)
library(readr)


# Filter to TE-overlapping DARs only
te_dars <- df %>% dplyr::filter(overlapping_TE == TRUE)

te_dars <- te_dars %>%
  filter(pct.1 >= 0.10 | pct.2 >= 0.10)

# Extract summit chromosome and position
te_dars <- te_dars %>%
  mutate(
    chrom = str_extract(summit, "^chr[^:]+"),
    summit_pos = as.integer(str_extract(summit, "(?<=:)[0-9]+"))
  )

# Create GRanges object for summit positions
summit_gr <- GRanges(
  seqnames = te_dars$chrom,
  ranges = IRanges(start = te_dars$summit_pos, end = te_dars$summit_pos),
  strand = "*"
)

# Use UCSC-based annotation to get promoters
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoters_10kb <- promoters(transcripts(txdb), upstream = 10000, downstream = 10000)

# Prepare overlap metadata
tx_ids <- mcols(promoters_10kb)$tx_id
strands <- as.character(strand(promoters_10kb))
all_transcripts <- transcripts(txdb)  # full transcript coordinates for TSS extraction

# Find overlaps between summit and promoters
overlap_hits <- findOverlaps(summit_gr, promoters_10kb)

# Extract info from overlaps
overlapping_tx_ids <- tx_ids[subjectHits(overlap_hits)]
overlapping_strands <- strands[subjectHits(overlap_hits)]

# Get actual TSS position from transcripts
tss_pos <- ifelse(
  overlapping_strands == "+",
  start(all_transcripts)[subjectHits(overlap_hits)],
  end(all_transcripts)[subjectHits(overlap_hits)]
)

# Calculate signed distance (positive = upstream, negative = downstream)
signed_distance_to_TSS <- te_dars$summit_pos[queryHits(overlap_hits)] - tss_pos
signed_distance_to_TSS <- ifelse(overlapping_strands == "-", -signed_distance_to_TSS, signed_distance_to_TSS)

# Create main tibble with overlaps
gene_match_df <- tibble(
  row_index = queryHits(overlap_hits),
  tx_id = overlapping_tx_ids,
  strand = overlapping_strands,
  peak = te_dars$peak[queryHits(overlap_hits)],
  full_peak = te_dars$full_peak[queryHits(overlap_hits)],
  summit = te_dars$summit[queryHits(overlap_hits)],
  avg_log2FC = te_dars$avg_log2FC[queryHits(overlap_hits)],
  summit_pos = te_dars$summit_pos[queryHits(overlap_hits)],
  chrom = te_dars$chrom[queryHits(overlap_hits)],
  signed_distance_to_TSS = signed_distance_to_TSS
)

# Map transcript ID → gene ID
tx2gene <- AnnotationDbi::select(txdb,
                  keys = as.character(unique(gene_match_df$tx_id)),
                  keytype = "TXID",
                  columns = "GENEID")

# Map gene ID → gene symbol
gene_info <- AnnotationDbi::select(org.Mm.eg.db,
                    keys = unique(tx2gene$GENEID),
                    columns = "SYMBOL",
                    keytype = "ENTREZID")

# Merge gene mappings
gene_match_df <- gene_match_df %>%
  mutate(tx_id = as.character(tx_id)) %>%
  left_join(tx2gene, by = c("tx_id" = "TXID")) %>%
  left_join(gene_info, by = c("GENEID" = "ENTREZID")) %>%
  dplyr::rename(gene_symbol = SYMBOL)

# Keep only the closest TE-DAR to TSS per gene
shortest_df <- gene_match_df %>%
  dplyr::filter(!is.na(gene_symbol)) %>%
  group_by(gene_symbol) %>%
  slice_min(order_by = abs(signed_distance_to_TSS), n = 1, with_ties = FALSE) %>%
  ungroup()

# Add combined label for plot
shortest_df <- shortest_df %>%
  mutate(label = paste0(gene_symbol, " (", full_peak, ")"))

library(ggtext)

shortest_df <- shortest_df %>%
  mutate(
    label = paste0("(", full_peak, ") <b>", gene_symbol, "</b>")
  )

p <- ggplot(shortest_df, aes(x = signed_distance_to_TSS, y = reorder(label, signed_distance_to_TSS))) +
  geom_col(aes(fill = pmin(pmax(avg_log2FC, -2), 2)), width = 0.6, color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_fill_gradient2(
    low = "#4575b4", mid = "white", high = "#d73027",
    midpoint = 0, limits = c(-2, 2),
    name = expression(log[2]~Fold~Change~"(EPO vs Placebo)"),
    guide = guide_colorbar(barwidth = 1.5, barheight = 8)
  ) +
  scale_x_continuous(
    labels = scales::comma_format(scale = 1),
    limits = c(-10000, 10000),
    breaks = seq(-10000, 10000, 5000)
  ) +
  labs(
    title = "Filtered TE-DARs (≥10% Cells Expressing) Within ±10 kb of TSS",
    subtitle = "Bar direction reflects TE-DAR position relative to TSS (upstream/downstream)",
    x = "Signed Distance to TSS (bp)",
    y = "Peak (Gene)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = ggtext::element_markdown(size = 10),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10)),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

save_highres_plot(p, filename = "Figure6/TE-DARs Within ±10 kb of TSS ")
 


# Load libraries
library(readr)
library(dplyr)
library(ggplot2)

# Load data
te_df <- read_tsv("Figure6/Combined_TE_enrichment_by_repName_for_Public_ChIP_Seq_data.txt")

# Define exclusion list
exclude_classes <- c("Low_complexity", "Satellite", "rRNA", "tRNA", "scRNA", "snRNA",
                     "RNA", "RC", "srpRNA", "DNA?", "LTR?", "LINE?", "RC?", "SINE?")


# Calculate -log10(pvalue)
te_df <- te_df %>%
  mutate(log10_p = -log10(ifelse(Pvalue == 0, 1e-300, Pvalue)))

# Count significantly enriched TE families per TF (p < 0.001)
tf_summary <- te_df %>%
  dplyr::filter(Pvalue < 0.001) %>%
  group_by(TF) %>%
  summarise(Enriched_TEs = n_distinct(TE_name)) %>%
  arrange(desc(Enriched_TEs))

# Plot
p <- ggplot(tf_summary, aes(x = reorder(TF, Enriched_TEs), y = Enriched_TEs)) +
  geom_col(fill = "#1f78b4", color = "black", width = 0.7) +
  geom_text(
    aes(label = Enriched_TEs),
    hjust = -0.15, size = 4, fontface = "bold", color = "black"
  ) +
  coord_flip(clip = "off") +
  theme_minimal(base_size = 16) +
  labs(
    title = "Number of Significantly Enriched TE Families per TF",
    x = NULL,
    y = "Number of Enriched TE Families"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 17),
    axis.title.y = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.text = element_text(color = "black", size = 13),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 30, 10, 10)  # allow room for text labels
  ) +
  expand_limits(y = max(tf_summary$Enriched_TEs) * 1.15)

save_highres_plot(p, filename = "Figure6/Number of Significantly Enriched TE Families per TF for ChIP-peaks", height = 4)


# Load data
te_df <- read_tsv("Figure6/Combined_TE_enrichment_by_repName_for_Public_ChIP_Seq_data.txt")

# Prepare binary matrix
binary_matrix <- te_df %>%
  dplyr::filter(!is.na(TE_name)) %>%
  mutate(significant = ifelse(Pvalue < 0.001, 1, 0)) %>%
  dplyr::select(TE_name, TF, significant) %>%
  distinct() %>%
  pivot_wider(names_from = TF, values_from = significant, values_fill = 0)

# Convert to matrix
mat <- as.data.frame(binary_matrix)
rownames(mat) <- mat$TE_name
mat <- mat[, -1]

# Sort rows by total TF presence (descending)
mat <- mat[order(rowSums(mat), decreasing = TRUE), ]

# Plot heatmap (no row clustering)
p <- pheatmap(
  mat,
  color = c("white", "black"),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  legend_breaks = c(0, 1),
  legend_labels = c("Not Enriched", "Enriched"),
  fontsize = 12,
  fontsize_row = 8,
  fontsize_col = 10,
  main = "TF Enrichment Across TE Families (p < 0.001)",
  show_rownames = FALSE,
  show_colnames = TRUE,
  angle_col = 90
)

save_highres_plot(p, filename = "Figure6/heatmap_for_TFs", width = 6.3)


# Load data
te_df <- read_tsv("Figure6/Combined_TE_enrichment_by_repName_for_Public_ChIP_Seq_data.txt")

# Create binary enrichment matrix
binary_matrix <- te_df %>%
  dplyr::filter(!is.na(TE_name)) %>%
  mutate(significant = ifelse(Pvalue < 0.001, 1, 0)) %>%
  dplyr::select(TE_name, TF, significant) %>%
  distinct() %>%
  pivot_wider(names_from = TF, values_from = significant, values_fill = 0)

# Convert to matrix and set row names
mat <- as.data.frame(binary_matrix)
rownames(mat) <- mat$TE_name
mat <- mat[, -1]

# Filter to TE families enriched in ≤ 8 TFs
mat_zoomed <- mat[rowSums(mat) >= 8, ]

# Sort rows from high to low
mat_zoomed <- mat_zoomed[order(rowSums(mat_zoomed), decreasing = TRUE), ]

mat_zoomed <- mat_zoomed[!grepl("\\)n$", rownames(mat_zoomed)), ]

# Plot heatmap
p <- pheatmap(
  mat_zoomed,
  color = c("white", "black"),
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  border_color = NA,
  legend_breaks = c(0, 1),
  legend_labels = c("Not Enriched", "Enriched"),
  fontsize = 16,
  fontsize_row = 16,
  fontsize_col = 16,
  main = "TE Families Significantly (p < 0.001) Enriched in ≥ 8 TFs",
  show_rownames = TRUE,
  show_colnames = TRUE,
  angle_col = 90
)

save_highres_plot(p, filename = "Figure6/TE Families Enriched in at least 8 TFs")




# Load data
te_df <- read_tsv("Figure6/Combined_TE_enrichment_by_repName_for_Public_ChIP_Seq_data.txt")

# Count how many TFs each TE is significantly enriched in (p < 0.001)
te_tf_counts <- te_df %>%
  dplyr::filter(Pvalue < 0.001, !is.na(TE_name)) %>%
  distinct(TE_name, TF) %>%
  count(TE_name, name = "TF_count")

# Histogram: how many TEs have N significant TF overlaps
ggplot(te_tf_counts, aes(x = TF_count)) +
  geom_histogram(binwidth = 1, fill = "#4daf4a", color = "black") +
  scale_x_continuous(breaks = 1:max(te_tf_counts$TF_count)) +
  labs(
    title = "Distribution of TE Families by TF Enrichment",
    x = "Number of TFs a TE Family is Significantly Enriched In",
    y = "Number of TE Families"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )


library(tidyr)
library(dplyr)
library(readr)
library(pheatmap)

# Load the data
te_df <- read_tsv("Figure6/Combined_TE_enrichment_by_repName_for_Public_ChIP_Seq_data.txt")

# Create matrix of -log10(p-values)
pval_matrix <- te_df %>%
  dplyr::filter(!is.na(TE_name)) %>%
  mutate(
    pval_capped = ifelse(Pvalue == 0, 1e-10, Pvalue),  # just in case
    log10_p = -log10(pval_capped)
  ) %>%
  dplyr::select(TE_name, TF, log10_p) %>%
  distinct() %>%
  pivot_wider(names_from = TF, values_from = log10_p, values_fill = 0)

# Convert to matrix
mat <- as.data.frame(pval_matrix)
rownames(mat) <- mat$TE_name
mat <- mat[, -1]

# Filter TE families enriched in ≥2 TFs (p < 0.001 → log10 > 3)
#mat <- mat[rowSums(mat > 3) >= 2, ]

# Sort rows by total enrichment signal
mat <- mat[order(rowSums(mat), decreasing = TRUE), ]

# Use soft log scale (since max is probably ~3–4)
p <- pheatmap(
  mat,
  color = colorRampPalette(c("white", "#fddede", "red2"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA,
  legend = TRUE,
  fontsize = 11,
  fontsize_row = 8,
  fontsize_col = 10,
  main = expression("-log"[10]*"(p-value) of TE Family Enrichment per TF"),
  show_rownames = FALSE,
  show_colnames = TRUE,
  angle_col = 90
)

save_highres_plot(p, filename =  "Figure6/heatmap_for_TFs_with_pvalues", width = 6.3)








library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)

# Load data
te_df <- read_tsv("Figure6/Combined_TE_enrichment_by_repName_for_Public_ChIP_Seq_data.txt")

# Count how many TFs each TE is significantly enriched in (p < 0.001)
te_tf_counts <- te_df %>%
  dplyr::filter(Pvalue < 0.001, !is.na(TE_name)) %>%
  distinct(TE_name, TF) %>%
  count(TE_name, name = "TF_count")

# Histogram: how many TE families fall into each TF count bucket
ggplot(te_tf_counts, aes(x = TF_count)) +
  geom_histogram(binwidth = 1, fill = "#d73027", color = "black") +
  scale_x_continuous(breaks = 1:max(te_tf_counts$TF_count)) +
  labs(
    title = "Distribution of TE Families by Number of TFs with Significant Enrichment",
    x = "Number of TFs (p < 0.001)",
    y = "Number of TE Families"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )







library(ggplot2)
library(dplyr)
library(readr)

# Load data
te_df <- read_tsv("Figure6/Combined_TE_enrichment_by_repName_for_Public_ChIP_Seq_data.txt")

# Count TFs per TE family
te_tf_counts <- te_df %>%
  dplyr::filter(Pvalue < 0.001, !is.na(TE_name)) %>%
  distinct(TE_name, TF) %>%
  count(TE_name, name = "TF_count")

# Count number of TE families in each TF_count
distribution_df <- te_tf_counts %>%
  count(TF_count, name = "Num_TEs")

# Lollipop plot
ggplot(distribution_df, aes(x = TF_count, y = Num_TEs)) +
  geom_segment(aes(xend = TF_count, yend = 0), color = "grey40", size = 0.8) +
  geom_point(color = "#d73027", size = 4) +
  geom_text(aes(label = Num_TEs), vjust = -0.8, fontface = "bold", size = 4) +
  scale_x_continuous(breaks = 1:max(distribution_df$TF_count)) +
  labs(
    title = "Distribution of TE Families by Number of TFs with Significant Enrichment",
    subtitle = "p < 0.001 for enrichment across neurogenic TFs",
    x = "Number of TFs a TE Family is Enriched In",
    y = "Number of TE Families"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  ) +
  expand_limits(y = max(distribution_df$Num_TEs) * 1.15)








library(readr)
library(dplyr)
library(tidyr)
library(ComplexUpset)
library(ggplot2)

# Load your TE enrichment data
te_df <- read_tsv("Figure6/Combined_TE_enrichment_by_repName_for_Public_ChIP_Seq_data.txt")

# Create binary matrix of TE × TF (significant if p < 0.001)
binary_df <- te_df %>%
  dplyr::filter(Pvalue < 0.001, !is.na(TE_name)) %>%
  distinct(TE_name, TF) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = TF, values_from = value, values_fill = 0)


# Convert to data.frame with binary columns (required format)
binary_df <- as.data.frame(binary_df)

# Plot the UpSet chart for TF columns
ComplexUpset::upset(
  binary_df,
  intersect = colnames(binary_df)[-1],  # All TFs
  name = "TE Family Overlap",
  base_annotations = list(
    'Intersection size' = intersection_size(
      text = list(vjust = -0.5, size = 3)
    )
  ),
  width_ratio = 0.2,
  sort_sets = 'descending',
  themes = upset_modify_themes(
    list(
      'intersect_size' = theme(axis.text.x = element_text(angle = 45, hjust = 1)),
      'set_size' = theme(axis.text.y = element_text(size = 10))
    )
  )
)






library(ggplot2)
library(dplyr)
library(readr)

# Load TE enrichment data
te_df <- read_tsv("Figure6/Combined_TE_enrichment_by_repName_for_Public_ChIP_Seq_data.txt")

# Count number of TFs each TE family is enriched in (p < 0.001)
te_tf_counts <- te_df %>%
  dplyr::filter(Pvalue < 0.001, !is.na(TE_name)) %>%
  distinct(TE_name, TF) %>%
  count(TE_name, name = "TF_count")

# Histogram data: How many TE families fall into each TF count
hist_data <- te_tf_counts %>%
  count(TF_count, name = "Num_TE_families")

# Plot histogram (as bar chart)
p <- ggplot(hist_data, aes(x = TF_count, y = Num_TE_families)) +
  geom_col(fill = "#fc8d62", color = "black", width = 0.7) +
  geom_text(aes(label = Num_TE_families), vjust = -0.5, fontface = "bold", size = 3.8) +
  scale_x_continuous(breaks = seq(min(hist_data$TF_count), max(hist_data$TF_count), 1)) +
  labs(
    title = "Distribution of TE Families by Number of TFs with Significant Enrichment",
    subtitle = "Threshold: p < 0.001",
    x = "Number of TFs in Which a TE Family Is Significantly Enriched",
    y = "Number of TE Families"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
    
  ) +
  expand_limits(y = max(hist_data$Num_TE_families) * 1.15)


save_highres_plot(p, filename = "Figure6/Number of TFs in Which a TE Family Is Significantly Enriched", height = 5)



library(dplyr)
library(ggplot2)

# Define TE class colors
te_colors <- c(
  "DNA"           = "#E41A1C",
  "LINE"          = "#FF7F00",
  "LTR"           = "#4DAF4A",
  "SINE"          = "#377EB8",
  "Simple_repeat" = "#984EA3",
  "Other"         = "#A65628",
  "Unknown"       = "#999999"
)

# Filter to overlapping TEs
te_only <- df %>% dplyr::filter(overlapping_TE == TRUE)

# Count TE_name × TE_class
te_name_class_counts <- te_only %>%
  count(TE_name, TE_class, sort = TRUE)

# Compute TE_class totals for legend
class_counts <- te_name_class_counts %>%
  group_by(TE_class) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  mutate(
    legend_label = paste0(TE_class, " (", total, ")"),
    legend_label = factor(legend_label, levels = legend_label)
  )

# Merge class info
te_name_class_counts <- te_name_class_counts %>%
  left_join(class_counts, by = "TE_class")

# Order TE names
te_name_class_counts <- te_name_class_counts %>%
  mutate(TE_name = factor(TE_name, levels = unique(TE_name[order(n)])))

# Set fill color mapping
color_map <- setNames(te_colors[class_counts$TE_class], class_counts$legend_label)

# Final plot
p <- ggplot(te_name_class_counts, aes(x = TE_name, y = n, fill = legend_label)) +
  geom_col(color = "black") +
  geom_text(aes(label = n), hjust = -0.25, size = 3.2, fontface = "bold", color = "black") +
  coord_flip() +
  scale_fill_manual(
    name = "TE Class (Count)",
    values = color_map
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribution of TE Names in Overlapping TE-DARs",
    x = "TE Name",
    y = "Count"
  ) +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5, color = "black"),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  expand_limits(y = max(te_name_class_counts$n) * 1.1)

# Optional save command (if you have a custom function for export)
# save_highres_plot(p, filename = "Figure6/Distribution of TE Names in Overlapping TE-DARs", height = 4)

print(p)



# Load libraries
library(readr)
library(dplyr)
library(ggplot2)

# Load the data
df <- read_csv("Data/Pyramidal/EPO_vs_Placebo/DA_peaks_combined_pyramidal_mature_vs_newly_formed_celltypes_EPO_vs_Placebo_separation_filtered_annotated_for_genes_and_TEs.csv")
df <- df %>%
  mutate(
    same_fc = TE_family == TE_class & !TE_class %in% c("Simple_repeat", "Other"),
    TE_class  = if_else(same_fc, "Unknown", TE_class),
    TE_family = if_else(same_fc, "Unknown", TE_family)
  ) %>%
  select(-same_fc)

head(df)
library(dplyr)
df <- df %>%
  mutate(
    cluster = case_when(
      comparison == "newly_formed_EPO_vs_newly_formed_Placebo" ~ "Cluster_newly_formed",
      comparison == "mature_EPO_vs_mature_Placebo" ~ "Cluster_mature",
      TRUE ~ NA_character_
    )
  )
dim(df)
library(dplyr)

std_chroms <- paste0("chr", c(1:19, "X", "Y", "M"))

df <- df %>%
  mutate(chrom = sub("[-:].*$", "", peak)) %>%  # get chr from 'peak'
  filter(chrom %in% std_chroms) %>%
  select(-chrom)  # remove helper column if not needed

dim(df)



# Define exclusion list
exclude_classes <- c("Low_complexity", "Satellite", "rRNA", "tRNA", "scRNA", "snRNA",  
                     "RNA", "RC", "srpRNA", "DNA?", "LTR?", "LINE?", "RC?", "SINE?")

# Modify 'overlapping_TE' in-place without removing rows
df <- df %>%
  mutate(overlapping_TE = ifelse(
    TE_class %in% exclude_classes |
      TE_family %in% exclude_classes |
      grepl("\\?", TE_family),
    FALSE,
    overlapping_TE
  ))

# Extract summit position (format: chr:start-end → we use 'start')
df <- df %>%
  mutate(
    chrom = str_extract(summit, "^chr[^:]+"),
    summit_pos = as.integer(str_extract(summit, "(?<=:)[0-9]+"))
  )

# Create GRanges object for all summits
summit_gr <- GRanges(
  seqnames = df$chrom,
  ranges = IRanges(start = df$summit_pos, end = df$summit_pos),
  strand = "*"
)

# Get ±10 kb promoter regions from mm10 gene TSSs
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoter_regions <- promoters(transcripts(txdb), upstream = 10000, downstream = 10000)

# Find overlaps
hits <- findOverlaps(summit_gr, promoter_regions)

# Annotate full data with summit proximity
df$within_10kb_promoter_summit <- FALSE
df$within_10kb_promoter_summit[queryHits(hits)] <- TRUE


library(AnnotationDbi)
library(org.Mm.eg.db)  # For mapping gene IDs to symbols

# Extract matching rows from both df and promoter_regions
summit_df_hits <- df[queryHits(hits), ]
promoter_hits <- promoter_regions[subjectHits(hits)]

# Get gene IDs from promoter regions
gene_ids <- mcols(promoter_hits)$tx_id

# Map transcript ID → gene ID
tx2gene <- AnnotationDbi::select(txdb, keys = as.character(gene_ids), columns = "GENEID", keytype = "TXID")

# Add gene symbols using org.Mm.eg.db
gene_info <- AnnotationDbi::select(org.Mm.eg.db,
                    keys = unique(tx2gene$GENEID),
                    columns = c("SYMBOL"),
                    keytype = "ENTREZID")

# Merge mappings back
tx2gene <- left_join(tx2gene, gene_info, by = c("GENEID" = "ENTREZID"))

# Associate each summit with its mapped gene(s)
summit_gene_table <- tibble(
  row_index = queryHits(hits),
  summit_id = df$peak[queryHits(hits)],
  tx_id = as.character(mcols(promoter_hits)$tx_id)  # <— FIX here
) %>%
  left_join(tx2gene, by = c("tx_id" = "TXID")) %>%
  group_by(row_index, summit_id) %>%
  summarise(promoter_genes = paste(unique(SYMBOL), collapse = ", "), .groups = "drop")


# Add this back into the main dataframe
df$promoter_gene_symbols <- NA
df$promoter_gene_symbols[summit_gene_table$row_index] <- summit_gene_table$promoter_genes


# Filter to TE-overlapping DARs and categorize proximity

te_family_proximity <- df %>%
  dplyr::filter(overlapping_TE == TRUE, !is.na(TE_family)) %>%
  mutate(Proximity = ifelse(within_10kb_promoter_summit, "Within ±10kb", "Distal (>10kb)")) %>%
  count(TE_family, Proximity, sort = TRUE)

# Optional: limit to top N TE families for cleaner plot (e.g., top 20)
top_families <- te_family_proximity %>%
  group_by(TE_family) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  slice_max(order_by = total, n = 200) %>%
  pull(TE_family)

# Filter to top families (optional)
te_family_proximity <- te_family_proximity %>%
  dplyr::filter(TE_family %in% top_families)

# Order TE_family factor by total count
te_family_proximity <- te_family_proximity %>%
  group_by(TE_family) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  mutate(TE_family = factor(TE_family, levels = (unique(TE_family[order(total)]))))

# Plot
p <- ggplot(te_family_proximity, aes(x = TE_family, y = n, fill = Proximity)) +
  geom_col(color = "black", width = 0.7) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold", color = "black") +
  coord_flip() +
  scale_fill_manual(
    values = c("Within ±10kb" = "#1f78b4", "Distal (>10kb)" = "#999999"),
    breaks = c("Within ±10kb", "Distal (>10kb)")   # 👈 legend order
  ) +
  labs(
    title = "Promoter Proximity of TE-DARs by TE Family",
    x = "TE Family",
    y = "Count of TE-DARs"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid = element_blank()
    
  )

# Display the plot
print(p)

save_highres_plot(p, "Figure6/Promoter Proximity of TE-DARs by TE Family Graph", height = 6)












# Load required packages
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(dplyr)
library(readr)
library(stringr)

# Read your input file (adjust path if needed)
df <- read_csv("Data/Pyramidal/EPO_vs_Placebo/DA_peaks_combined_pyramidal_mature_vs_newly_formed_celltypes_EPO_vs_Placebo_separation_filtered_annotated_for_genes_and_TEs.csv")

# Modify 'overlapping_TE' in-place without removing rows
df <- df %>%
  mutate(overlapping_TE = ifelse(
    TE_class %in% exclude_classes |
      TE_family %in% exclude_classes |
      grepl("\\?", TE_family),
    FALSE,
    overlapping_TE
  ))

# Extract summit position
df <- df %>%
  mutate(
    chrom = str_extract(summit, "^chr[^:]+"),
    summit_pos = as.integer(str_extract(summit, "(?<=:)[0-9]+"))
  )

# Filter to TE-overlapping DARs
te_dars <- df %>% filter(overlapping_TE == TRUE)

# Create GRanges object for summits
summit_gr <- GRanges(
  seqnames = te_dars$chrom,
  ranges = IRanges(start = te_dars$summit_pos, end = te_dars$summit_pos),
  strand = "*"
)

# Define promoter regions ±10kb
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoters_10kb <- promoters(transcripts(txdb), upstream = 10000, downstream = 10000)

# Find overlaps between summits and promoter regions
overlap_hits <- findOverlaps(summit_gr, promoters_10kb)
tx_ids <- mcols(promoters_10kb)$tx_id[subjectHits(overlap_hits)]

# Map transcript IDs to gene IDs and gene symbols
tx2gene <- AnnotationDbi::select(txdb, keys = as.character(tx_ids), columns = "GENEID", keytype = "TXID")
gene_info <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = unique(tx2gene$GENEID),
                                   columns = c("SYMBOL"),
                                   keytype = "ENTREZID")

# Combine mappings
tx2gene <- tx2gene %>% left_join(gene_info, by = c("GENEID" = "ENTREZID"))
gene_symbols <- unique(tx2gene$SYMBOL)
gene_symbols <- gene_symbols[!is.na(gene_symbols)]

# Perform GO enrichment with SYMBOLs
ego_symbol <- enrichGO(
  gene          = gene_symbols,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# View top terms
head(ego_symbol)

# Plot (optional)
barplot(ego_symbol, showCategory = 50, title = "GO Enrichment for Genes Near TE-bound DARs (Gene Symbols)")
 
























library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(stringr)
library(tidyr)
set.seed(0)



# Load the data
df <- read_csv("Data/Pyramidal/EPO_vs_Placebo/DA_peaks_combined_pyramidal_mature_vs_newly_formed_celltypes_EPO_vs_Placebo_separation_filtered_annotated_for_genes_and_TEs.csv")
head(df)
library(dplyr)
df <- df %>%
  mutate(
    cluster = case_when(
      comparison == "newly_formed_EPO_vs_newly_formed_Placebo" ~ "Cluster_newly_formed",
      comparison == "mature_EPO_vs_mature_Placebo" ~ "Cluster_mature",
      TRUE ~ NA_character_
    )
  )
dim(df)
library(dplyr)

std_chroms <- paste0("chr", c(1:19, "X", "Y", "M"))

df <- df %>%
  mutate(chrom = sub("[-:].*$", "", peak)) %>%  # get chr from 'peak'
  filter(chrom %in% std_chroms) %>%
  select(-chrom)  # remove helper column if not needed

dim(df)
# Define exclusion list
exclude_classes <- c("Low_complexity", "Satellite", "rRNA", "tRNA", "scRNA", "snRNA",  
                     "RNA", "RC", "srpRNA", "DNA?", "LTR?", "LINE?", "RC?", "SINE?")

# Modify 'overlapping_TE' in-place without removing rows
df <- df %>%
  mutate(overlapping_TE = ifelse(
    TE_class %in% exclude_classes |
      TE_family %in% exclude_classes |
      grepl("\\?", TE_family),
    FALSE,
    overlapping_TE
  ))

da_peaks_df <- df
da_peaks_df <- da_peaks_df[da_peaks_df$p_val_adj < 0.01,]

# Extract chr/start/end from 'peak' column
peak_split <- str_split_fixed(da_peaks_df$peak, "-", 3)
da_peaks_df <- da_peaks_df %>%
  mutate(
    seqnames = peak_split[,1],
    start = as.numeric(peak_split[,2]),
    end = as.numeric(peak_split[,3])
  )

# Function to count overlaps
count_overlaps <- function(chip_gr, da_gr) {
  overlaps <- countOverlaps(da_gr, chip_gr)
  sum(overlaps > 0)
}

# Load all TF peak files
tf_files <- list.files("Figure6/final_bed_files/final_bed_files/", pattern = "\\.bed$", full.names = TRUE)
tf_names <- sub("\\.bed$", "", basename(tf_files))
tf_gr_list <- setNames(lapply(tf_files, import), tf_names)

# Combine all TFs into one GRanges object
all_tfs_combined <- do.call(c, unname(tf_gr_list))

# Initialize result lists
per_tf_results <- list()
combined_results <- list()

# Loop through each cluster
for (clust in unique(da_peaks_df$cluster)) {
  
  # Subset DA peaks
  subset_df <- da_peaks_df %>% filter(cluster == clust)
  
  # Convert to GRanges
  da_gr <- GRanges(
    seqnames = subset_df$seqnames,
    ranges = IRanges(start = subset_df$start, end = subset_df$end)
  )
  
  total_da <- length(da_gr)
  
  # Overlap with each TF
  tf_overlap_counts <- sapply(tf_gr_list, count_overlaps, da_gr = da_gr)
  
  per_tf_results[[clust]] <- data.frame(
    cluster = clust,
    TF = names(tf_overlap_counts),
    Overlapping_DA_Peaks = tf_overlap_counts,
    Total_DA_Peaks = total_da
  )
  
  # Overlap with any TF
  overlap_any_tf <- countOverlaps(da_gr, all_tfs_combined)
  combined_results[[clust]] <- data.frame(
    cluster = clust,
    Total_DA_Peaks = total_da,
    DA_Peaks_Overlapping_Any_TF = sum(overlap_any_tf > 0),
    Percent_Overlap_Any_TF = round(100 * sum(overlap_any_tf > 0) / total_da, 2)
  )
}

# Merge all results
final_tf_overlap <- bind_rows(per_tf_results) %>%
  mutate(Percent_Overlap_Per_TF = round(100 * Overlapping_DA_Peaks / Total_DA_Peaks, 2))

final_combined_overlap <- bind_rows(combined_results)

# Save to CSV
write.csv(final_tf_overlap, "Figure6/DA_peak_overlap_per_TF_per_cluster.csv", row.names = FALSE)
write.csv(final_combined_overlap, "Figure6/DA_peak_overlap_with_any_TF_per_cluster.csv", row.names = FALSE)

# Print summaries
cat("✅ Per TF Overlaps:\n")
print(final_tf_overlap)

cat("\n✅ Combined TF Overlaps:\n")
print(final_combined_overlap)










library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(stringr)
library(tidyr)
set.seed(0)



# Load the data
df <- read_csv("Data/Pyramidal/EPO_vs_Placebo/DA_peaks_combined_pyramidal_mature_vs_newly_formed_celltypes_EPO_vs_Placebo_separation_filtered_annotated_for_genes_and_TEs.csv")
head(df)
library(dplyr)
df <- df %>%
  mutate(
    cluster = case_when(
      comparison == "newly_formed_EPO_vs_newly_formed_Placebo" ~ "Cluster_newly_formed",
      comparison == "mature_EPO_vs_mature_Placebo" ~ "Cluster_mature",
      TRUE ~ NA_character_
    )
  )
dim(df)
library(dplyr)

std_chroms <- paste0("chr", c(1:19, "X", "Y", "M"))

df <- df %>%
  mutate(chrom = sub("[-:].*$", "", peak)) %>%  # get chr from 'peak'
  filter(chrom %in% std_chroms) %>%
  select(-chrom)  # remove helper column if not needed

dim(df)
# Define exclusion list
exclude_classes <- c("Low_complexity", "Satellite", "rRNA", "tRNA", "scRNA", "snRNA",  
                     "RNA", "RC", "srpRNA", "DNA?", "LTR?", "LINE?", "RC?", "SINE?")

# Modify 'overlapping_TE' in-place without removing rows
df <- df %>%
  mutate(overlapping_TE = ifelse(
    TE_class %in% exclude_classes |
      TE_family %in% exclude_classes |
      grepl("\\?", TE_family),
    FALSE,
    overlapping_TE
  ))

da_peaks_df <- df

da_peaks_df <- da_peaks_df %>%
  mutate(overlapping_TE = ifelse(
    TE_class %in% exclude_classes |
      TE_family %in% exclude_classes |
      grepl("\\?", TE_family),
    FALSE,
    overlapping_TE
  ))

da_peaks_df <- da_peaks_df[da_peaks_df$overlapping_TE == TRUE,]

da_peaks_df <- da_peaks_df[da_peaks_df$p_val_adj < 0.01,]

# Extract chr/start/end from 'peak' column
peak_split <- str_split_fixed(da_peaks_df$peak, "-", 3)
da_peaks_df <- da_peaks_df %>%
  mutate(
    seqnames = peak_split[,1],
    start = as.numeric(peak_split[,2]),
    end = as.numeric(peak_split[,3])
  )

# Function to count overlaps
count_overlaps <- function(chip_gr, da_gr) {
  overlaps <- countOverlaps(da_gr, chip_gr)
  sum(overlaps > 0)
}

# Load all TF peak files
tf_files <- list.files("Figure6/final_bed_files/final_bed_files/", pattern = "\\.bed$", full.names = TRUE)
tf_names <- sub("\\.bed$", "", basename(tf_files))
tf_gr_list <- setNames(lapply(tf_files, import), tf_names)

# Combine all TFs into one GRanges object
all_tfs_combined <- do.call(c, unname(tf_gr_list))

# Initialize result lists
per_tf_results <- list()
combined_results <- list()

# Loop through each cluster
for (clust in unique(da_peaks_df$cluster)) {
  
  # Subset DA peaks
  subset_df <- da_peaks_df %>% filter(cluster == clust)
  
  # Convert to GRanges
  da_gr <- GRanges(
    seqnames = subset_df$seqnames,
    ranges = IRanges(start = subset_df$start, end = subset_df$end)
  )
  
  total_da <- length(da_gr)
  
  # Overlap with each TF
  tf_overlap_counts <- sapply(tf_gr_list, count_overlaps, da_gr = da_gr)
  
  per_tf_results[[clust]] <- data.frame(
    cluster = clust,
    TF = names(tf_overlap_counts),
    Overlapping_DA_Peaks = tf_overlap_counts,
    Total_DA_Peaks = total_da
  )
  
  # Overlap with any TF
  overlap_any_tf <- countOverlaps(da_gr, all_tfs_combined)
  combined_results[[clust]] <- data.frame(
    cluster = clust,
    Total_DA_Peaks = total_da,
    DA_Peaks_Overlapping_Any_TF = sum(overlap_any_tf > 0),
    Percent_Overlap_Any_TF = round(100 * sum(overlap_any_tf > 0) / total_da, 2)
  )
}

# Merge all results
final_tf_overlap <- bind_rows(per_tf_results) %>%
  mutate(Percent_Overlap_Per_TF = round(100 * Overlapping_DA_Peaks / Total_DA_Peaks, 2))

final_combined_overlap <- bind_rows(combined_results)

# Save to CSV
write.csv(final_tf_overlap, "Figure6/TE-DA_peak_overlap_per_TF_per_cluster.csv", row.names = FALSE)
write.csv(final_combined_overlap, "Figure6/TE-DA_peak_overlap_with_any_TF_per_cluster.csv", row.names = FALSE)

# Print summaries
cat("✅ Per TF Overlaps:\n")
print(final_tf_overlap)

cat("\n✅ Combined TF Overlaps:\n")
print(final_combined_overlap)



















library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(stringr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(openxlsx)

# === Load gene annotations ===
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes_tx <- genes(txdb)

promoter_sites <- GRanges(
  seqnames = seqnames(genes_tx),
  ranges = IRanges(
    start = ifelse(strand(genes_tx) == "+", start(genes_tx), end(genes_tx)),
    end = ifelse(strand(genes_tx) == "+", start(genes_tx), end(genes_tx))
  ),
  strand = strand(genes_tx),
  gene_id = names(genes_tx)
)

# === Load DA peaks ===
da_peaks_df 

peak_split <- str_split_fixed(da_peaks_df$peak, "-", 3)
da_peaks_df <- da_peaks_df %>%
  mutate(
    seqnames = peak_split[, 1],
    start = as.numeric(peak_split[, 2]),
    end = as.numeric(peak_split[, 3])
  )

# === Load ChIP-seq TF peaks ===
tf_files <- list.files("Figure6/final_bed_files/final_bed_files/", pattern = "\\.bed$", full.names = TRUE)
tf_names <- sub("\\.bed$", "", basename(tf_files))

# Safely import all TF GRanges
tf_gr_list <- lapply(tf_files, function(f) {
  gr <- try(import(f), silent = TRUE)
  if (inherits(gr, "try-error")) return(NULL)
  return(gr)
})
valid_idx <- !sapply(tf_gr_list, is.null)
tf_gr_list <- setNames(tf_gr_list[valid_idx], tf_names[valid_idx])

# === Initialize output ===
go_results_list <- list()
pdf("Figure6/GO_enrichment_DApeaks_perCluster_perTF.pdf", width = 14, height = 10)
wb <- createWorkbook()

# === Loop over clusters and TFs ===
for (clust in unique(da_peaks_df$cluster)) {
  message("🔬 Cluster: ", clust)
  
  cluster_df <- da_peaks_df %>% filter(cluster == clust)
  da_gr <- GRanges(seqnames = cluster_df$seqnames, ranges = IRanges(cluster_df$start, cluster_df$end))
  
  for (tf in names(tf_gr_list)) {
    tf_gr <- tf_gr_list[[tf]]
    
    message("   ➤ TF: ", tf)
    
    # Overlap
    overlap_hits <- findOverlaps(da_gr, tf_gr)
    if (length(overlap_hits) == 0) {
      message("      ⚠️ No overlaps")
      next
    }
    
    da_overlap <- da_gr[queryHits(overlap_hits)]
    
    # Find nearest promoters
    nearest_idx <- nearest(da_overlap, promoter_sites)
    nearest_promoters <- promoter_sites[nearest_idx]
    gene_ids <- unique(mcols(nearest_promoters)$gene_id)
    gene_ids <- gene_ids[!is.na(gene_ids)]
    
    if (length(gene_ids) < 10) {
      message("      ⚠️ Not enough genes (", length(gene_ids), ") for GO")
      next
    }
    
    # GO enrichment
    ego <- enrichGO(
      gene = gene_ids,
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.01,
      readable = TRUE
    )
    
    if (!is.null(ego) && !is.null(ego@result) && nrow(ego@result) > 0) {
      # Try-catch in case pairwise_termsim fails for edge cases
      try({
        ego <- pairwise_termsim(ego)
        
        # Save to result list
        short_clust <- gsub("Cluster_", "", clust)
        short_clust <- gsub(" ", "_", short_clust)
        tag <- paste(short_clust, tf, sep = "__")
        tag <- substr(tag, 1, 31)
        go_results_list[[tag]] <- ego
        
        # Plot
        print(barplot(ego, showCategory = 30, title = paste(clust, "-", tf)))
        print(dotplot(ego, showCategory = 30, title = paste(clust, "-", tf)))
        print(emapplot(ego, showCategory = 30, layout = "kk") + ggtitle(paste(clust, "-", tf)))
        
        # Save to Excel sheet
        addWorksheet(wb, tag)
        writeData(wb, sheet = tag, ego@result)
      }, silent = TRUE)
    } else {
      message("      ⚠️ No significant GO terms")
    }
  }
}

dev.off()

# Save workbook
saveWorkbook(wb, "Figure6/GO_enrichment_results_perCluster_perTF.xlsx", overwrite = TRUE)
message("✅ GO enrichment done for all cluster × TF combinations")






library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(forcats)

# Load data
df <- read_csv("Figure6/DA_peak_overlap_per_TF_per_cluster.csv")

# === Clean and format data ===
df <- df %>%
  mutate(
    # Remove underscores and prefix
    cluster = str_replace_all(cluster, "_", " "),
    cluster = str_remove(cluster, regex("^Cluster\\s*", ignore_case = TRUE)),
    cluster = str_to_title(cluster),  # Capitalize first letter of each word (e.g., "Newly Formed")
    
    TF = str_replace_all(TF, "_", " "),
    TF = toupper(TF),  # Uppercase all TF names
    
    cluster_TF = paste(cluster, TF, sep = " | "),  # Construct axis label
    
    # Assign group
    cluster_lower = tolower(cluster),
    type = case_when(
      str_detect(cluster_lower, "new|immature|nf|ascl1|nrg3") ~ "Newly Formed",
      str_detect(cluster_lower, "mature|ca1|ca3|dentate|dg|ctip2|deep") ~ "Mature",
      TRUE ~ "Other"
    )
  ) %>%
  filter(type != "Other") %>%
  mutate(
    type = factor(type, levels = c("Newly Formed", "Mature")),
    cluster_TF = fct_reorder(cluster_TF, Percent_Overlap_Per_TF),
    label = paste0(Overlapping_DA_Peaks, "/", Total_DA_Peaks)
  )

# === Create plot ===
p <- ggplot(df, aes(x = Percent_Overlap_Per_TF, y = cluster_TF, fill = type)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_text(aes(label = label), hjust = -0.1, size = 3.5) +
  facet_wrap(~type, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("Newly Formed" = "orange", "Mature" = "red2")) +
  labs(
    title = "% Overlap of DARs with TFs",
    x = "Percent Overlap (%)",
    y = "Cluster | TF"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    strip.text = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    legend.position = "none",
    panel.grid.major.y = element_blank()
  ) +
  xlim(0, max(df$Percent_Overlap_Per_TF) + 10)

save_highres_plot(p, filename = "Figure6/Percent Overlap of DA Peaks with Transcription Factors", height = 6)


library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(forcats)

# === Load TF plot to match the order ===
df_tf <- read_csv("Figure6/DA_peak_overlap_per_TF_per_cluster.csv") %>%
  mutate(
    cluster = str_replace_all(cluster, "_", " "),
    cluster = str_remove(cluster, regex("^Cluster\\s*", ignore_case = TRUE)),
    cluster = str_to_title(cluster),
    TF = str_replace_all(TF, "_", " "),
    TF = toupper(TF),
    cluster_TF = paste(cluster, TF, sep = " | "),
    cluster_lower = tolower(cluster),
    type = case_when(
      str_detect(cluster_lower, "new|immature|nf|ascl1|nrg3") ~ "Newly Formed",
      str_detect(cluster_lower, "mature|ca1|ca3|dentate|dg|ctip2|deep") ~ "Mature",
      TRUE ~ "Other"
    )
  ) %>%
  filter(type != "Other") %>%
  mutate(cluster_TF = fct_reorder(cluster_TF, Percent_Overlap_Per_TF))

ordered_levels <- levels(df_tf$cluster_TF)

# === Load TE data ===
df_te <- read_csv("Figure6/TE-DA_peak_overlap_per_TF_per_cluster.csv") %>%
  mutate(
    cluster = str_replace_all(cluster, "_", " "),
    cluster = str_remove(cluster, regex("^Cluster\\s*", ignore_case = TRUE)),
    cluster = str_to_title(cluster),
    TF = str_replace_all(TF, "_", " "),
    TF = toupper(TF),
    cluster_TF = paste(cluster, TF, sep = " | "),
    cluster_lower = tolower(cluster),
    type = case_when(
      str_detect(cluster_lower, "new|immature|nf|ascl1|nrg3") ~ "Newly Formed",
      str_detect(cluster_lower, "mature|ca1|ca3|dentate|dg|ctip2|deep") ~ "Mature",
      TRUE ~ "Other"
    )
  ) %>%
  filter(type != "Other") %>%
  mutate(
    type = factor(type, levels = c("Newly Formed", "Mature")),
    cluster_TF = factor(cluster_TF, levels = ordered_levels),
    label = paste0(Overlapping_DA_Peaks, "/", Total_DA_Peaks)
  )

# === Plot ===
p <- ggplot(df_te, aes(x = Percent_Overlap_Per_TF, y = cluster_TF, fill = type)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_text(aes(label = label), hjust = 1.1, size = 3.5) +
  facet_wrap(~type, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("Newly Formed" = "orange", "Mature" = "red2")) +
  scale_x_continuous(
    trans = "reverse",
    breaks = c(0, 20, 40, 60),
    expand = expansion(mult = c(0.05, 0))
  ) +
  scale_y_discrete(position = "right") +
  labs(
    title = "% Overlap of TE-derived DARs with TFs",
    x = "Percent Overlap (%)",
    y = "Cluster | TF"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    strip.text = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    legend.position = "none",
    panel.grid.major.y = element_blank()
  )

save_highres_plot(p, filename = "Figure6/Percent Overlap of DA Peaks with TEs with TFs", height = 6, width = 10)

 


library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# Load and tag datasets
df_all <- read_csv("Figure6/DA_peak_overlap_with_any_TF_per_cluster.csv") %>%
  mutate(DAR_type = "All DARs")

df_te <- read_csv("Figure6/TE-DA_peak_overlap_with_any_TF_per_cluster.csv") %>%
  mutate(DAR_type = "TE-derived DARs")

# Combine and clean
df_combined <- bind_rows(df_all, df_te) %>%
  mutate(
    cluster = str_replace_all(cluster, "_", " "),
    Cluster = case_when(
      str_detect(tolower(cluster), "new") ~ "Newly Formed",
      str_detect(tolower(cluster), "mature") ~ "Mature",
      TRUE ~ NA_character_
    ),
    Cluster = factor(Cluster, levels = c("Newly Formed", "Mature")),
    DAR_type = factor(DAR_type, levels = c("All DARs", "TE-derived DARs")),
    label = paste0(DA_Peaks_Overlapping_Any_TF, "/", Total_DA_Peaks)
  )

# Plot with horizontal bars
p <- ggplot(df_combined, aes(x = Percent_Overlap_Any_TF, y = DAR_type, fill = Cluster)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.9),
    height = 0.4
  ) +
  geom_text(
    aes(label = label),
    position = position_dodge(width = 0.9),
    hjust = -0.1,
    size = 4
  ) +
  scale_fill_manual(values = c("Newly Formed" = "orange", "Mature" = "red2")) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Percent of DAR Types Overlapping Any TF Peak",
    x = "Percent Overlap (%)",
    y = "DAR Type",
    fill = "Cluster"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "top"
  )


save_highres_plot(p, filename = "Figure6/Percent of DAR Types Overlapping Any TF Peak", height = 5, width = 18)


p <- ggplot(df_combined, aes(x = Percent_Overlap_Any_TF, y = DAR_type, fill = Cluster)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.9),
    height = 0.4
  ) +
  geom_text(
    aes(label = label),
    position = position_dodge(width = 0.9),
    hjust = 1.1,  # aligns text inside bars since bars go right to left
    size = 4
  ) +
  scale_fill_manual(values = c("Newly Formed" = "orange", "Mature" = "red2")) +
  scale_x_reverse(
    limits = c(100, 0),
    breaks = c(100, 75, 50, 25, 0),
    expand = expansion(mult = c(0.05, 0))
  ) +
  scale_y_discrete(position = "right") +
  labs(
    title = "Percent of DAR Types Overlapping Any TF Peak",
    x = "Percent Overlap (%)",
    y = "DAR Type",
    fill = "Cluster"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "top"
  )

save_highres_plot(p, filename = "Figure6/Percent of DAR Types Overlapping Any TF Peak_y_axis_right", height = 4, width = 18)














library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(stringr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(openxlsx)

# === Load gene annotations ===
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes_tx <- transcripts(txdb)

promoter_sites <- GRanges(
  seqnames = seqnames(genes_tx),
  ranges = IRanges(
    start = ifelse(strand(genes_tx) == "+", start(genes_tx), end(genes_tx)),
    end = ifelse(strand(genes_tx) == "+", start(genes_tx), end(genes_tx))
  ),
  strand = strand(genes_tx),
  gene_id = names(genes_tx)
)

# === Load DA peaks ===
da_peaks_df 

peak_split <- str_split_fixed(da_peaks_df$peak, "-", 3)
da_peaks_df <- da_peaks_df %>%
  mutate(
    seqnames = peak_split[, 1],
    start = as.numeric(peak_split[, 2]),
    end = as.numeric(peak_split[, 3])
  )

# === Load ChIP-seq TF peaks ===
tf_files <- list.files("Figure6/final_bed_files/final_bed_files/", pattern = "\\.bed$", full.names = TRUE)
tf_names <- sub("\\.bed$", "", basename(tf_files))

# Safely import all TF GRanges
tf_gr_list <- lapply(tf_files, function(f) {
  gr <- try(import(f), silent = TRUE)
  if (inherits(gr, "try-error")) return(NULL)
  return(gr)
})
valid_idx <- !sapply(tf_gr_list, is.null)
tf_gr_list <- setNames(tf_gr_list[valid_idx], tf_names[valid_idx])

# === Initialize output ===
go_results_list <- list()
pdf("Figure6/GO_enrichment_DApeaks_perCluster_perTF.pdf", width = 14, height = 10)
wb <- createWorkbook()

# === Loop over clusters and TFs ===
for (clust in unique(da_peaks_df$cluster)) {
  message("🔬 Cluster: ", clust)
  
  cluster_df <- da_peaks_df %>% filter(cluster == clust)
  da_gr <- GRanges(seqnames = cluster_df$seqnames, ranges = IRanges(cluster_df$start, cluster_df$end))
  
  for (tf in names(tf_gr_list)) {
    tf_gr <- tf_gr_list[[tf]]
    
    message("   ➤ TF: ", tf)
    
    # Overlap
    overlap_hits <- findOverlaps(da_gr, tf_gr)
    if (length(overlap_hits) == 0) {
      message("      ⚠️ No overlaps")
      next
    }
    
    da_overlap <- da_gr[queryHits(overlap_hits)]
    
    # Find nearest promoters
    nearest_idx <- nearest(da_overlap, promoter_sites)
    nearest_promoters <- promoter_sites[nearest_idx]
    gene_ids <- unique(mcols(nearest_promoters)$gene_id)
    gene_ids <- gene_ids[!is.na(gene_ids)]
    
    if (length(gene_ids) < 10) {
      message("      ⚠️ Not enough genes (", length(gene_ids), ") for GO")
      next
    }
    
    # GO enrichment
    ego <- enrichGO(
      gene = gene_ids,
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.01,
      readable = TRUE
    )
    
    if (!is.null(ego) && !is.null(ego@result) && nrow(ego@result) > 0) {
      # Try-catch in case pairwise_termsim fails for edge cases
      try({
        ego <- pairwise_termsim(ego)
        
        # Save to result list
        short_clust <- gsub("Cluster_", "", clust)
        short_clust <- gsub(" ", "_", short_clust)
        tag <- paste(short_clust, tf, sep = "__")
        tag <- substr(tag, 1, 31)
        go_results_list[[tag]] <- ego
        
        # Plot
        print(barplot(ego, showCategory = 30, title = paste(clust, "-", tf)))
        print(dotplot(ego, showCategory = 30, title = paste(clust, "-", tf)))
        print(emapplot(ego, showCategory = 30, layout = "kk") + ggtitle(paste(clust, "-", tf)))
        
        # Save to Excel sheet
        addWorksheet(wb, tag)
        writeData(wb, sheet = tag, ego@result)
      }, silent = TRUE)
    } else {
      message("      ⚠️ No significant GO terms")
    }
  }
}

dev.off()

# Save workbook
saveWorkbook(wb, "Figure6/GO_enrichment_results_perCluster_perTF.xlsx", overwrite = TRUE)
message("✅ GO enrichment done for all cluster × TF combinations")













library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(stringr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(openxlsx)
library(AnnotationDbi)

# === Load gene annotations using transcripts ===
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes_tx <- transcripts(txdb)

# Map transcript IDs to gene IDs
valid_txids <- as.character(mcols(genes_tx)$tx_id)

tx2gene <- AnnotationDbi::select(
  txdb,
  keys = valid_txids,
  keytype = "TXID",
  columns = "GENEID"
)


# Attach matched gene_id to genes_tx
genes_tx$gene_id <- tx2gene$GENEID[match(mcols(genes_tx)$tx_id, tx2gene$TXID)]

# Create promoter GRanges (1 bp TSS)
promoter_sites <- GRanges(
  seqnames = seqnames(genes_tx),
  ranges = IRanges(
    start = ifelse(strand(genes_tx) == "+", start(genes_tx), end(genes_tx)),
    end = ifelse(strand(genes_tx) == "+", start(genes_tx), end(genes_tx))
  ),
  strand = strand(genes_tx),
  gene_id = genes_tx$gene_id
)

# === Load DA peaks ===
# Replace this with your actual data loading line
da_peaks_df <- read_csv("Data/Pyramidal/EPO_vs_Placebo/DA_peaks_combined_pyramidal_mature_vs_newly_formed_celltypes_EPO_vs_Placebo_separation_filtered_annotated_for_genes_and_TEs.csv")

# Modify 'overlapping_TE' in-place without removing rows
da_peaks_df <- da_peaks_df %>%
  mutate(overlapping_TE = ifelse(
    TE_class %in% exclude_classes |
      TE_family %in% exclude_classes |
      grepl("\\?", TE_family),
    FALSE,
    overlapping_TE
  ))

da_peaks_df <- da_peaks_df[da_peaks_df$overlapping_TE == TRUE, ]

peak_split <- str_split_fixed(da_peaks_df$peak, "-", 3)
da_peaks_df <- da_peaks_df %>%
  mutate(
    seqnames = peak_split[, 1],
    start = as.numeric(peak_split[, 2]),
    end = as.numeric(peak_split[, 3])
  )

# === Load ChIP-seq TF peaks ===
tf_files <- list.files("Figure6/final_bed_files/final_bed_files/", pattern = "\\.bed$", full.names = TRUE)
tf_names <- sub("\\.bed$", "", basename(tf_files))

# Safely import all TF GRanges
tf_gr_list <- lapply(tf_files, function(f) {
  gr <- try(import(f), silent = TRUE)
  if (inherits(gr, "try-error")) return(NULL)
  return(gr)
})
valid_idx <- !sapply(tf_gr_list, is.null)
tf_gr_list <- setNames(tf_gr_list[valid_idx], tf_names[valid_idx])

# === Initialize output ===
go_results_list <- list()
pdf("Figure6/GO_enrichment_DApeaks_perCluster_perTF.pdf", width = 14, height = 10)
wb <- createWorkbook()

# === Loop over clusters and TFs ===
for (clust in unique(da_peaks_df$cluster)) {
  message("🔬 Cluster: ", clust)
  
  cluster_df <- da_peaks_df %>% filter(cluster == clust)
  da_gr <- GRanges(seqnames = cluster_df$seqnames, ranges = IRanges(cluster_df$start, cluster_df$end))
  
  for (tf in names(tf_gr_list)) {
    tf_gr <- tf_gr_list[[tf]]
    
    message("   ➤ TF: ", tf)
    
    # Overlap
    overlap_hits <- findOverlaps(da_gr, tf_gr)
    if (length(overlap_hits) == 0) {
      message("      ⚠️ No overlaps")
      next
    }
    
    da_overlap <- da_gr[queryHits(overlap_hits)]
    
    # Find nearest promoters
    nearest_idx <- nearest(da_overlap, promoter_sites)
    nearest_promoters <- promoter_sites[nearest_idx]
    gene_ids <- unique(mcols(nearest_promoters)$gene_id)
    gene_ids <- gene_ids[!is.na(gene_ids)]
    
    if (length(gene_ids) < 10) {
      message("      ⚠️ Not enough genes (", length(gene_ids), ") for GO")
      next
    }
    
    # GO enrichment
    ego <- enrichGO(
      gene = gene_ids,
      OrgDb = org.Mm.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.01,
      readable = TRUE
    )
    
    if (!is.null(ego) && !is.null(ego@result) && nrow(ego@result) > 0) {
      try({
        ego <- pairwise_termsim(ego)
        
        # Save results
        short_clust <- gsub("Cluster_", "", clust)
        short_clust <- gsub(" ", "_", short_clust)
        tag <- paste(short_clust, tf, sep = "__")
        tag <- substr(tag, 1, 31)
        go_results_list[[tag]] <- ego
        
        # Plot
        print(barplot(ego, showCategory = 30, title = paste(clust, "-", tf)))
        print(dotplot(ego, showCategory = 30, title = paste(clust, "-", tf)))
        print(emapplot(ego, showCategory = 30, layout = "kk") + ggtitle(paste(clust, "-", tf)))
        
        # Save to Excel
        addWorksheet(wb, tag)
        writeData(wb, sheet = tag, ego@result)
      }, silent = TRUE)
    } else {
      message("      ⚠️ No significant GO terms")
    }
  }
}

dev.off()
saveWorkbook(wb, "Figure6/GO_enrichment_results_perCluster_perTF.xlsx", overwrite = TRUE)
message("✅ GO enrichment done for all cluster × TF combinations")















library(dplyr)
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(stringr)
library(GenomeInfoDb)
library(tidyr)

# === Load and filter TE-derived DARs ===
da_peaks_df <- read_csv(
  "Data/Pyramidal/EPO_vs_Placebo/DA_peaks_combined_pyramidal_mature_vs_newly_formed_celltypes_EPO_vs_Placebo_separation_filtered_annotated_for_genes_and_TEs.csv",
  show_col_types = FALSE
)
da_peaks_df <- da_peaks_df[da_peaks_df$p_val_adj < 0.01, ]
da_peaks_df <- da_peaks_df[da_peaks_df$avg_log2FC > 0, ]
da_peaks_df <- da_peaks_df %>% dplyr::select(1:12)  # keep first 12 columns

# --- Exclusions for classes (and we’ll drop families containing '?') ---
exclude_classes <- c(
  "Low_complexity","Satellite","rRNA","tRNA","scRNA","snRNA",
  "RNA","RC","srpRNA","DNA?","LTR?","LINE?","RC?","SINE?"
)

# ========= Read RepeatMasker TSV and FILTER on import (NO +1 CONVERSION) =========
rmsk_tsv <- "Data/mm10_RepeatMasker/mm10_RepeatMasker.tsv"  # or .tsv.gz
rmsk_df <- readr::read_tsv(rmsk_tsv, show_col_types = FALSE)

# Keep only rows where class NOT excluded and family does NOT contain "?"
rmsk_df <- rmsk_df %>%
  mutate(
    repClass  = as.character(repClass),
    repFamily = as.character(repFamily),
    repFamily_nafix = ifelse(is.na(repFamily), "", repFamily)
  ) %>%
  filter(!(repClass %in% exclude_classes | grepl("\\?", repFamily_nafix))) %>%
  select(-repFamily_nafix)

# Build GRanges DIRECTLY from TSV (no coordinate shift)
rmsk_gr <- GRanges(
  seqnames = rmsk_df$genoName,
  ranges   = IRanges(start = rmsk_df$genoStart, end = rmsk_df$genoEnd),
  strand   = rmsk_df$strand
)

mcols(rmsk_gr)$repName   <- rmsk_df$repName
mcols(rmsk_gr)$repClass  <- rmsk_df$repClass
mcols(rmsk_gr)$repFamily <- rmsk_df$repFamily


# === Convert ALL peaks to GRanges (do NOT prune this; keep same length as df) ===
stopifnot("peak" %in% colnames(da_peaks_df))
peak_split <- str_split_fixed(da_peaks_df$peak, "-", 3)
peak_gr <- GRanges(
  seqnames = peak_split[, 1],
  ranges   = IRanges(
    start = as.numeric(peak_split[, 2]),
    end   = as.numeric(peak_split[, 3])
  )
)

# Harmonize seqlevel style but DON'T drop any peaks
ref_style <- seqlevelsStyle(rmsk_gr)[1]
seqlevelsStyle(peak_gr) <- ref_style

# Work on a subset copy that lives on common chromosomes; keep an index map
common_seqs <- intersect(seqlevels(peak_gr), seqlevels(rmsk_gr))
idx_common  <- which(as.character(seqnames(peak_gr)) %in% common_seqs)
peak_gr_sub <- peak_gr[idx_common]
rmsk_gr_sub <- keepSeqlevels(rmsk_gr, common_seqs, pruning.mode = "coarse")

# ===== Summit (midpoint) then extend ±250 bp (width = 501) =====
peak_summit_sub <- resize(peak_gr_sub, width = 1,    fix = "center")
peak_win500_sub  <- resize(peak_summit_sub, width = 501, fix = "center")

# (safety: clamp starts to >=1; keeps width if any would become <1)
st <- start(peak_win500_sub); st[st < 1] <- 1L; start(peak_win500_sub) <- st

# Overlap the ±500 windows with filtered RepeatMasker
hits <- findOverlaps(peak_win500_sub, rmsk_gr_sub, ignore.strand = TRUE)

# ===== Build overlap flag for FULL dataframe length (avoid length mismatch) =====
overlap_flag <- rep(FALSE, length(peak_gr))
if (length(hits) > 0) {
  overlap_flag[idx_common[unique(queryHits(hits))]] <- TRUE
}
da_peaks_df$overlapping_TE <- overlap_flag

# ===== (Optional) TE annotations per peak using the ±500 window =====
if (length(hits) > 0) {
  ol_bp <- width(pintersect(ranges(peak_win500_sub)[queryHits(hits)],
                            ranges(rmsk_gr_sub)[subjectHits(hits)]))
  
  ov_tbl <- tibble::tibble(
    row       = idx_common[queryHits(hits)],  # map back to original row index
    repName   = mcols(rmsk_gr_sub)$repName[subjectHits(hits)],
    repClass  = mcols(rmsk_gr_sub)$repClass[subjectHits(hits)],
    repFamily = mcols(rmsk_gr_sub)$repFamily[subjectHits(hits)],
    overlap_bp = ol_bp
  )
  
  collapsed <- ov_tbl %>%
    dplyr::group_by(row) %>%
    dplyr::summarise(
      TE_name   = paste(unique(stats::na.omit(repName)),   collapse = ";"),
      TE_class  = paste(unique(stats::na.omit(repClass)),  collapse = ";"),
      TE_family = paste(unique(stats::na.omit(repFamily)), collapse = ";"),
      .groups = "drop"
    )
  
  best <- ov_tbl %>%
    dplyr::group_by(row) %>%
    dplyr::slice_max(overlap_bp, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(row,
                     TE_name_top   = repName,
                     TE_class_top  = repClass,
                     TE_family_top = repFamily,
                     TE_overlap_bp = overlap_bp
    )
  
  ann <- dplyr::full_join(collapsed, best, by = "row")
  
  da_peaks_df <- da_peaks_df %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    dplyr::left_join(ann, by = "row") %>%
    dplyr::select(-row)
} else {
  da_peaks_df <- da_peaks_df %>%
    dplyr::mutate(
      TE_name = NA_character_, TE_class = NA_character_, TE_family = NA_character_,
      TE_name_top = NA_character_, TE_class_top = NA_character_, TE_family_top = NA_character_,
      TE_overlap_bp = NA_integer_
    )
}

# Define TE-derived DARs as those overlapping via the ±500bp window
te_derived_dars_df <- da_peaks_df %>% dplyr::filter(overlapping_TE)
te_derived_dars_gr <- peak_gr[da_peaks_df$overlapping_TE]

# Summits for TE-derived DARs (still 1bp center if you need them later)
te_summits <- resize(te_derived_dars_gr, width = 1, fix = "center")

cat("All peaks:", length(peak_gr),
    "| TE-derived DARs (summit ±100bp):", length(te_derived_dars_gr), "\n")


# === Load TF peaks ===
tf_files <- list.files("Figure6/final_bed_files/final_bed_files/", pattern = "\\.bed$", full.names = TRUE)
tf_names <- sub("\\.bed$", "", basename(tf_files))
tf_gr_list <- setNames(lapply(tf_files, import), tf_names)

# === Combine all TF peaks ===
all_tfs_gr <- do.call(c, unname(tf_gr_list))

# === Filter TE-DARs overlapping any TF ===
overlap_hits <- findOverlaps(te_derived_dars_gr, all_tfs_gr)
te_dars_with_tf <- te_derived_dars_gr[queryHits(overlap_hits)]

# === Use summits of TE-DARs ===
te_summits <- resize(te_dars_with_tf, width = 1, fix = "center")

# === Get transcript annotations ===
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes_tx <- transcripts(txdb)

# === Map transcripts to genes ===
tx2gene <- AnnotationDbi::select(
  txdb,
  keys = as.character(mcols(genes_tx)$tx_id),
  keytype = "TXID",
  columns = "GENEID"
)
genes_tx$gene_id <- tx2gene$GENEID[match(mcols(genes_tx)$tx_id, tx2gene$TXID)]

# === Define TSS sites ===
tss_sites <- GRanges(
  seqnames = seqnames(genes_tx),
  ranges = IRanges(
    start = ifelse(strand(genes_tx) == "+", start(genes_tx), end(genes_tx)),
    end = ifelse(strand(genes_tx) == "+", start(genes_tx), end(genes_tx))
  ),
  strand = strand(genes_tx),
  gene_id = genes_tx$gene_id
)

# === Extend TSS ±10 kb ===
tss_extended <- resize(tss_sites, width = 20001, fix = "center")

# === Find TSSs overlapping TE-DAR summits ===
tss_hits <- findOverlaps(tss_extended, te_summits)
nearby_gene_ids <- unique(mcols(tss_extended)$gene_id[queryHits(tss_hits)])
nearby_gene_ids <- nearby_gene_ids[!is.na(nearby_gene_ids)]

# === Convert to gene symbols ===
gene_symbols <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = nearby_gene_ids,
  columns = "SYMBOL",
  keytype = "ENTREZID"
)

unique_gene_symbols <- unique(gene_symbols$SYMBOL)
cat("✅ Number of nearby genes:", length(unique_gene_symbols), "\n")
head(unique_gene_symbols)
writeClipboard(paste(unique_gene_symbols, collapse = "\n"))



# Ensure symbols are unique and non-missing
gene_symbols_clean <- unique(na.omit(unique_gene_symbols))

# --- Run GO enrichment (unchanged) ---
ego_BP <- enrichGO(gene = gene_symbols_clean, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                   ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
ego_CC <- enrichGO(gene = gene_symbols_clean, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                   ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
ego_MF <- enrichGO(gene = gene_symbols_clean, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
                   ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)


library(dplyr)
library(tidytext)   # for reorder_within / scale_y_reordered
library(ggplot2)
library(scales)

# helper: parse "a/b" to numeric ratio
.parse_ratio <- function(x) {
  x <- as.character(x)
  sapply(strsplit(x, "/"), function(y) as.numeric(y[1]) / as.numeric(y[2]))
}

as_go_df <- function(ego, label, top_n = 30) {
  df <- as.data.frame(ego)
  if (is.null(df) || nrow(df) == 0) return(tibble())
  df <- df %>%
    mutate(
      ontology      = label,
      GeneRatio_num = .parse_ratio(GeneRatio),
      negLog10Padj  = -log10(p.adjust)
    )
  top_n <- min(top_n, nrow(df))  # <-- compute outside data-masking
  df %>% arrange(desc(GeneRatio_num)) %>% slice_head(n = top_n)
}

go_df <- bind_rows(
  as_go_df(ego_BP, "Biological Process", top_n = 20),
  as_go_df(ego_CC, "Cellular Component", top_n = 20),
  as_go_df(ego_MF, "Molecular Function", top_n = 20)
)

# guard in case nothing enriched
if (nrow(go_df) == 0) stop("No enriched GO terms at the given threshold.")

# order terms by GeneRatio within each ontology for plotting
go_df <- go_df %>%
  mutate(Description = reorder_within(Description, GeneRatio_num, ontology))

# cleaner dot plot (size = Count, color = -log10 q)
p <- ggplot(go_df, aes(x = GeneRatio_num, y = Description)) +
  geom_point(aes(size = Count, color = negLog10Padj)) +
  facet_wrap(~ ontology, scales = "free_y", nrow = 1) +
  scale_y_reordered() +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_viridis_c(name = expression(-log[10]~"q"), option = "plasma", direction = -1) +
  guides(size = guide_legend(title = "Gene count")) +
  labs(
    title = "GO enrichment for genes near TE-DARs",
    x = "Gene ratio", y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 17),
    strip.text  = element_text(size = 18, face = "bold"),  # ⟵ bigger facet titles
    # optional: soften the strip background
    # strip.background = element_rect(fill = "grey95", color = NA),
    legend.position = "right"
  )

p

save_highres_plot(p, filename = "Figure6/GO_enrichment_DARs_overlapping_with_TFs", width = 24, height = 8)



library(dplyr)
library(AnnotationDbi)
library(org.Mm.eg.db)

# Table of peaks (index aligns with te_summits / te_dars_with_tf)
peaks_tbl <- tibble(
  peak_idx = seq_along(te_dars_with_tf),
  peak     = paste0(as.character(seqnames(te_dars_with_tf)),
                    "-", start(te_dars_with_tf), "-", end(te_dars_with_tf))
)

# Use existing overlaps between TSS±10kb windows (tss_extended) and peak summits (te_summits)
# We now compute distance to the ACTUAL TSS (tss_sites) and choose the closest per peak.
if (length(tss_hits) > 0) {
  # indices for tss_sites and te_summits
  tss_idx  <- queryHits(tss_hits)     # rows in tss_extended (same order as tss_sites)
  peak_idx <- subjectHits(tss_hits)   # rows in te_summits (and te_dars_with_tf)
  
  # distance from summit (1bp) to TSS (1bp)
  dist_bp <- distance(tss_sites[tss_idx], te_summits[peak_idx])
  
  hits_df <- tibble(
    peak_idx = peak_idx,
    tss_idx  = tss_idx,
    gene_id  = mcols(tss_sites)$gene_id[tss_idx],
    dist_bp  = dist_bp
  ) %>%
    filter(!is.na(gene_id))
  
  # keep the single closest TSS per peak
  closest_df <- hits_df %>%
    group_by(peak_idx) %>%
    slice_min(dist_bp, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # map Entrez -> SYMBOL
  sym_map <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = unique(closest_df$gene_id),
                                   keytype = "ENTREZID",
                                   columns = "SYMBOL")
  
  peak_gene_df <- closest_df %>%
    left_join(sym_map, by = c("gene_id" = "ENTREZID")) %>%
    left_join(peaks_tbl, by = "peak_idx") %>%
    transmute(
      peak,
      closest_gene_symbol = SYMBOL,
      dist_to_TSS_bp      = dist_bp
    ) %>%
    arrange(dist_to_TSS_bp)
  
} else {
  peak_gene_df <- peaks_tbl %>%
    transmute(peak, closest_gene_symbol = NA_character_, dist_to_TSS_bp = NA_integer_)
}

# inspect / save
head(peak_gene_df)
# readr::write_csv(peak_gene_df, "TE_DAR_peaks_closest_TSS_gene.csv")

library(dplyr)
library(readr)

# Load original df (as you wrote)
da_peaks_df <- read_csv(
  "Data/Pyramidal/EPO_vs_Placebo/DA_peaks_combined_pyramidal_mature_vs_newly_formed_celltypes_EPO_vs_Placebo_separation_filtered_annotated_for_genes_and_TEs.csv",
  show_col_types = FALSE
)

# Ensure character type for joining key
da_peaks_df <- da_peaks_df %>% mutate(peak = as.character(peak))

# Make the peak->closest gene mapping UNIQUE by 'peak'
# (requires 'peak_gene_df' from the previous step with columns:
#  peak, closest_gene_symbol, dist_to_TSS_bp)
peak_gene_df_unique <- peak_gene_df %>%
  mutate(peak = as.character(peak)) %>%
  filter(!is.na(closest_gene_symbol)) %>%
  distinct(peak, .keep_all = TRUE)

# Join onto the original df
da_peaks_df <- da_peaks_df %>%
  left_join(peak_gene_df_unique, by = "peak")

# Quick summary
cat("Peaks annotated with a closest TSS gene:",
    sum(!is.na(da_peaks_df$closest_gene_symbol)), "of", nrow(da_peaks_df), "\n")

View(da_peaks_df)

da_peaks_df_with_gene <- da_peaks_df %>%
  filter(!is.na(closest_gene_symbol) & str_trim(closest_gene_symbol) != "")

# quick check
cat("Rows kept:", nrow(da_peaks_df_with_gene), "of", nrow(da_peaks_df), "\n")

da_peaks_df_with_peak_up <- da_peaks_df_with_gene %>%
  dplyr::filter(!is.na(avg_log2FC) & avg_log2FC > 0)

head(da_peaks_df_with_peak_up)


library(dplyr)
library(stringr)
library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)

# --- 1) Convert peaks to GRanges ---
peak_split <- str_split_fixed(da_peaks_df_with_peak_up$peak, "-", 3)
peak_gr <- GRanges(
  seqnames = peak_split[, 1],
  ranges   = IRanges(start = as.numeric(peak_split[, 2]),
                     end   = as.numeric(peak_split[, 3]))
)

# --- 2) Import all BED files ---
bed_files <- list.files("Figure6/final_bed_files/final_bed_files/",
                        pattern = "\\.bed$", full.names = TRUE)
stopifnot(length(bed_files) > 0)

bed_names <- sub("\\.bed$", "", basename(bed_files))
bed_list <- setNames(lapply(bed_files, import), bed_names)

# Harmonize seqlevels
ref_style <- seqlevelsStyle(bed_list[[1]])[1]
seqlevelsStyle(peak_gr) <- ref_style
bed_list <- lapply(bed_list, function(gr) { seqlevelsStyle(gr) <- ref_style; gr })

# --- 3) Check overlaps per BED file ---
overlap_flag <- rep(FALSE, length(peak_gr))
overlap_names <- vector("list", length(peak_gr))  # list to collect TF names

for (i in seq_along(bed_list)) {
  gr <- bed_list[[i]]
  hits <- findOverlaps(peak_gr, gr, ignore.strand = TRUE)
  
  if (length(hits) > 0) {
    overlap_flag[unique(queryHits(hits))] <- TRUE
    # Append file name to each overlapping peak
    for (q in unique(queryHits(hits))) {
      overlap_names[[q]] <- c(overlap_names[[q]], names(bed_list)[i])
    }
  }
}

# Collapse TF names into a semicolon-separated string
overlap_names_str <- sapply(overlap_names, function(x) {
  if (is.null(x)) NA_character_ else paste(unique(x), collapse = ";")
})

# --- 4) Add columns to dataframe ---
da_peaks_df_with_peak_up <- da_peaks_df_with_peak_up %>%
  mutate(
    overlaps_any_bed = overlap_flag,
    bed_files_overlap = overlap_names_str
  )

# quick check
table(da_peaks_df_with_peak_up$overlaps_any_bed)
head(da_peaks_df_with_peak_up %>% select(peak, overlaps_any_bed, bed_files_overlap))


library(dplyr)
library(ggplot2)
library(forcats)
 
da_peaks_df_with_peak_up_filtered <- da_peaks_df_with_peak_up %>%
  filter(pct.1 >= 0.10 | pct.2 >= 0.10)

da_peaks_df_with_peak_up_filtered <- da_peaks_df_with_peak_up_filtered[da_peaks_df_with_peak_up_filtered$dist_to_TSS_bp<10000,]


# 1) Average log2FC for duplicate peaks, cap at 4
heat_df <- da_peaks_df_with_peak_up_tf_filtered %>%
  group_by(peak) %>%
  summarise(avg_log2FC = mean(avg_log2FC, na.rm = TRUE), .groups = "drop") %>%
  mutate(avg_log2FC_capped = pmin(avg_log2FC, 4)) %>%
  arrange(desc(avg_log2FC_capped))

# 2) Single-column heatmap (white → red, 0–4)
ggplot(heat_df, aes(x = "avg_log2FC",
                    y = fct_reorder(peak, avg_log2FC_capped, .desc = TRUE),
                    fill = avg_log2FC_capped)) +
  geom_tile() +
  scale_fill_gradient(
    name = "avg log2FC",
    low = "white", high = "red",
    limits = c(0, 4), oob = scales::squish
  ) +
  labs(
    title = "Upregulated peaks (values > 4 capped at 4)",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# Optional: save dataframe
# readr::write_csv(heat_df, "peak_log2FC_single_column_capped.csv")



library(readr)
library(dplyr)
library(stringr)
library(GenomicRanges)
library(AnnotationDbi)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(rtracklayer)
library(pheatmap)

exclude_classes <- c(
  "Low_complexity","Satellite","rRNA","tRNA","scRNA","snRNA",
  "RNA","RC","srpRNA","DNA?","LTR?","LINE?","RC?","SINE?"
)

# === Load and filter TE-derived DARs ===
da_peaks_df <- read_csv("Data/Pyramidal/EPO_vs_Placebo/DA_peaks_combined_pyramidal_mature_vs_newly_formed_celltypes_EPO_vs_Placebo_separation_filtered_annotated_for_genes_and_TEs.csv")
da_peaks_df <- da_peaks_df[da_peaks_df$p_val_adj < 0.01,]
da_peaks_df <- da_peaks_df[da_peaks_df$avg_log2FC > 0,]
da_peaks_df <- da_peaks_df[da_peaks_df$pct.1 > 0.10 | da_peaks_df$pct.2 > 0.10,]

# Recalculate overlapping_TE flag
da_peaks_df <- da_peaks_df %>%
  mutate(overlapping_TE = ifelse(
    TE_class %in% exclude_classes |
      TE_family %in% exclude_classes |
      grepl("\\?", TE_family),
    FALSE,
    overlapping_TE
  ))

# Keep only TE-derived DARs
da_peaks_df <- da_peaks_df[da_peaks_df$overlapping_TE == TRUE, ]

# === Convert to GRanges and keep peak names ===
peak_split <- str_split_fixed(da_peaks_df$peak, "-", 3)
te_derived_dars_gr <- GRanges(
  seqnames = peak_split[, 1],
  ranges = IRanges(
    start = as.numeric(peak_split[, 2]),
    end = as.numeric(peak_split[, 3])
  ),
  peak_name = da_peaks_df$peak
)

# === Load TF peaks ===
tf_files <- list.files("Figure6/final_bed_files/final_bed_files/", pattern = "\\.bed$", full.names = TRUE)
tf_gr_list <- setNames(lapply(tf_files, import), sub("\\.bed$", "", basename(tf_files)))

# Combine all TF peaks
all_tfs_gr <- do.call(c, unname(tf_gr_list))

# === Filter TE-DARs overlapping any TF ===
overlap_hits <- findOverlaps(te_derived_dars_gr, all_tfs_gr)
te_dars_with_tf <- te_derived_dars_gr[queryHits(overlap_hits)]

# Use summits of TE-DARs
te_summits <- resize(te_dars_with_tf, width = 1, fix = "center")

# === Get transcript annotations ===
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
genes_tx <- transcripts(txdb)

# Map transcripts to genes
tx2gene <- AnnotationDbi::select(
  txdb,
  keys = as.character(mcols(genes_tx)$tx_id),
  keytype = "TXID",
  columns = "GENEID"
)
genes_tx$gene_id <- tx2gene$GENEID[match(mcols(genes_tx)$tx_id, tx2gene$TXID)]

# Define TSS sites
tss_sites <- GRanges(
  seqnames = seqnames(genes_tx),
  ranges = IRanges(
    start = ifelse(strand(genes_tx) == "+", start(genes_tx), end(genes_tx)),
    end = ifelse(strand(genes_tx) == "+", start(genes_tx), end(genes_tx))
  ),
  strand = strand(genes_tx),
  gene_id = genes_tx$gene_id
)

# Extend TSS ±10 kb
tss_extended <- resize(tss_sites, width = 20001, fix = "center")

# === Find TSSs overlapping TE-DAR summits ===
tss_hits <- findOverlaps(tss_extended, te_summits)

# Peak–gene mapping
peak_gene_df <- data.frame(
  peak = mcols(te_summits)$peak_name[subjectHits(tss_hits)],
  gene_id = mcols(tss_extended)$gene_id[queryHits(tss_hits)]
)
peak_gene_df <- peak_gene_df[!is.na(peak_gene_df$gene_id), ]

# Convert to gene symbols
gene_symbols <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = unique(peak_gene_df$gene_id),
  columns = "SYMBOL",
  keytype = "ENTREZID"
)

# Merge symbols with peaks
peak_gene_df <- merge(
  peak_gene_df,
  gene_symbols,
  by.x = "gene_id",
  by.y = "ENTREZID",
  all.x = TRUE
)
peak_gene_df <- unique(peak_gene_df)

cat("✅ Peaks with at least one nearby gene:", length(unique(peak_gene_df$peak)), "\n")
cat("✅ Unique nearby genes:", length(unique(peak_gene_df$SYMBOL)), "\n")

# === Merge with DA peaks to get log2FC ===
peak_gene_df2 <- merge(
  peak_gene_df,
  da_peaks_df[, c("peak", "avg_log2FC")],
  by = "peak",
  all.x = TRUE
)

# Replace NA with 0 and cap between 0–2
peak_gene_df2$avg_log2FC[is.na(peak_gene_df2$avg_log2FC)] <- 0
peak_gene_df2$avg_log2FC <- pmin(pmax(peak_gene_df2$avg_log2FC, 0), 2)

# Create row labels
peak_gene_df2$row_label <- paste0(peak_gene_df2$peak, " - ", peak_gene_df2$SYMBOL, "")

# If multiple genes per peak → average log2FC
heat_df <- peak_gene_df2 %>%
  group_by(row_label) %>%
  summarise(log2FC = mean(avg_log2FC, na.rm = TRUE)) %>%
  ungroup()

# Convert to matrix
heat_mat <- as.matrix(heat_df$log2FC)
rownames(heat_mat) <- heat_df$row_label

library(RColorBrewer)

# Remove unwanted row
heat_mat <- heat_mat[!rownames(heat_mat) %in% "chr15-79798080-79799685 - Npcd", , drop = FALSE]
heat_mat <- heat_mat[!rownames(heat_mat) %in% "chr4-41028935-41030559 - Aqp7", , drop = FALSE]


# Define a cleaner color palette
atac_colors <- colorRampPalette(c("white", brewer.pal(9, "Reds")[3:9]))(100)

# Publication-quality ATAC heatmap
atac_pheat <- pheatmap(
  heat_mat,
  cluster_rows   = TRUE,
  cluster_cols   = FALSE,
  color          = atac_colors,
  legend         = TRUE,
  legend_labels  = c("Low", "High"),
  border_color   = NA,
  main           = "TE-derived Upregulated DARs \n(≥10% Cells Expressing) Overlapping with TFs",
  fontsize       = 10,        # general text size
  fontsize_row   = 12,         # smaller row labels for many peaks
  fontsize_col   = 10,
  angle_col      = 45,        # angled column labels for readability
  treeheight_row = 12,        # smaller tree height for compactness
  show_rownames  = TRUE,      # show row labels
  show_colnames  = TRUE
)

library(Seurat)
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(pheatmap)

# === Load RNA object ===
RNA <- readRDS("../ATAC_analysis_with_combined_peaks/RNA_final_(without_A3_A5_B10_B12).RDS")

# === Set identities ===
RNA@meta.data$global_identity_treatment <- paste0(
  RNA@meta.data$global_identity, "_", RNA@meta.data$treatment
)
Idents(RNA) <- "global_identity_treatment"

# === Genes of interest: same as in ATAC heatmap ===
valid_genes <- intersect(unique(peak_gene_df2$SYMBOL), rownames(RNA))

# === Run FindMarkers for each cell type ===
all_identities <- sort(unique(RNA@meta.data$global_identity))

marker_list <- map(all_identities, function(ct) {
  
  id_epo     <- paste0(ct, "_Erythropoietin")
  id_placebo <- paste0(ct, "_Placebo")
  
  fm <- FindMarkers(
    object   = RNA,
    ident.1  = id_epo,
    ident.2  = id_placebo,
    features = valid_genes,
    test.use = "bimod",
    logfc.threshold = 0,
    min.pct = 0.01
  )
  
  fm$gene       <- rownames(fm)
  fm$cell_type  <- ct
  fm
})

markers_all <- bind_rows(marker_list)

# === Build wide matrices for log2FC and p-values ===
fc_mat <- markers_all %>%
  select(cell_type, gene, avg_log2FC) %>%
  pivot_wider(names_from = cell_type, values_from = avg_log2FC) %>%
  column_to_rownames("gene") %>%
  as.matrix()

p_mat <- markers_all %>%
  select(cell_type, gene, p_val_adj) %>%
  pivot_wider(names_from = cell_type, values_from = p_val_adj) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Replace missing values
fc_mat[is.na(fc_mat)] <- 0
p_mat[is.na(p_mat)]   <- 1

# === Cap log2FC values between -1 and 1 ===
fc_mat[fc_mat >  1] <-  1
fc_mat[fc_mat < -1] <- -1

# === Star matrix for significance ===
star_mat <- ifelse(p_mat < 0.01, "*", "")

# === Match row order with ATAC heatmap ===
# === Extract ATAC gene order from pheatmap object ===
# atac_pheat is the object returned when you plotted ATAC heatmap
atac_gene_order <- sub(".* - ", "", rownames(heat_mat)[atac_pheat$tree_row$order])

# === Keep only genes present in RNA matrix ===
valid_gene_order <- atac_gene_order[atac_gene_order %in% rownames(fc_mat)]

# === Reorder RNA matrices to match ATAC heatmap visual order ===
fc_mat_ordered <- fc_mat[match(valid_gene_order, rownames(fc_mat)), , drop = FALSE]
star_mat_ordered <- star_mat[match(valid_gene_order, rownames(star_mat)), , drop = FALSE]

# === Plot RNA-seq heatmap ===
RNA_pheatmap <- pheatmap(
  fc_mat_ordered,
  cluster_rows   = FALSE, # keep ATAC row order
  cluster_cols   = TRUE,
  color          = colorRampPalette(c("blue", "white", "red"))(100),
  main           = "RNA-seq: EPO vs Placebo (avg_log2FC, ★ = adj P < 0.01)",
  display_numbers = star_mat_ordered,
  number_color   = "black",
  fontsize_row   = 8,
  fontsize_number = 8
)


save_highres_plot(atac_pheat, filename = "Figure6/atac_pheat_for_TE-derived upregulated DARs overlapping with TFs", width = 5)
save_highres_plot(RNA_pheatmap, filename = "Figure6/RNA_pheat_for_TE-derived upregulated DARs overlapping with TFs", width = 7)








 

RNA <- readRDS("../ATAC_analysis_with_combined_peaks/RNA_final_(without_A3_A5_B10_B12).RDS")
library(Seurat)
library(dplyr)
library(purrr)      # map()
library(tidyr)      # pivot_wider()
library(tibble)     # column_to_rownames()
library(pheatmap)

# ---------------------------------------------------------------------------
# 1) Set the cell identities (already done earlier)
RNA@meta.data$global_identity_treatment <- paste0(
  RNA@meta.data$global_identity, "_", RNA@meta.data$treatment
)
Idents(RNA) <- "global_identity_treatment"

# ---------------------------------------------------------------------------
# 2) Loop through each global_identity and run FindMarkers -------------------
all_identities <- sort(unique(RNA@meta.data$global_identity))
available_genes <- rownames(RNA)
valid_genes <- intersect(unique(unique_gene_symbols), available_genes)

marker_list <- map(all_identities, function(ct) {
  
  id_epo     <- paste0(ct, "_Erythropoietin")
  id_placebo <- paste0(ct, "_Placebo")
  
  fm <- FindMarkers(
    object   = RNA,
    ident.1  = id_epo,
    ident.2  = id_placebo,
    features = valid_genes,   # TE-near genes
    test.use = "bimod",
    logfc.threshold = 0,      # keep every gene; we’ll decide significance later
    min.pct = 0.01
  )
  
  fm$gene       <- rownames(fm)
  fm$cell_type  <- ct        # tag with the originating identity
  fm
})

markers_all <- bind_rows(marker_list)

# ---------------------------------------------------------------------------
# 2) Wide matrices: log2FCs  and  P-values -----------------------------------
fc_mat <- markers_all %>%
  select(cell_type, gene, avg_log2FC) %>%
  pivot_wider(names_from = cell_type, values_from = avg_log2FC) %>%
  column_to_rownames("gene") %>%
  as.matrix()

p_mat <- markers_all %>%
  select(cell_type, gene, p_val_adj) %>%
  pivot_wider(names_from = cell_type, values_from = p_val_adj) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Replace missing values
fc_mat[is.na(fc_mat)] <- 0
p_mat[is.na(p_mat)]   <- 1   # treat as non-significant

# ---------------------------------------------------------------------------
# 3) Cap log2FC to ±2 --------------------------------------------------------
fc_mat[fc_mat >  1] <-  1
fc_mat[fc_mat < -1] <- -1

# ---------------------------------------------------------------------------
# 4) Build a “star” matrix for display_numbers -------------------------------
star_mat <- ifelse(p_mat < 0.01, "*", "")   # put star where adj-p < 0.01

# ---------------------------------------------------------------------------
# 5) Heat-map ----------------------------------------------------------------
p <- pheatmap(
  fc_mat,
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  color          = colorRampPalette(c("blue", "white", "red"))(100),
  main           = "EPO vs Placebo (avg_log2FC, ★ = adj P < 0.01)",
  display_numbers = star_mat,
  number_color   = "black",
  fontsize_row   = 12,
  fontsize_number = 12
)

p



RNA <- readRDS("../ATAC_analysis_with_combined_peaks/RNA_final_(without_A3_A5_B10_B12).RDS")
library(Seurat)
library(dplyr)
library(purrr)      # map()
library(tidyr)      # pivot_wider()
library(tibble)     # column_to_rownames()
library(pheatmap)

# ---------------------------------------------------------------------------
# 1) Set the cell identities (already done earlier)
RNA@meta.data$global_identity_treatment <- paste0(
  RNA@meta.data$global_identity, "_", RNA@meta.data$treatment
)
Idents(RNA) <- "global_identity_treatment"

# ---------------------------------------------------------------------------
# 2) Loop through each global_identity and run FindMarkers -------------------
all_identities <- sort(unique(RNA@meta.data$global_identity))
available_genes <- rownames(RNA)
valid_genes <- intersect(unique(da_peaks_df_with_peak_up_filtered$closest_gene_symbol), available_genes)
 
marker_list <- map(all_identities, function(ct) {
  
  id_epo     <- paste0(ct, "_Erythropoietin")
  id_placebo <- paste0(ct, "_Placebo")
  
  fm <- FindMarkers(
    object   = RNA,
    ident.1  = id_epo,
    ident.2  = id_placebo,
    features = valid_genes,   # TE-near genes
    test.use = "bimod",
    logfc.threshold = 0,      # keep every gene; we’ll decide significance later
    min.pct = 0.01
  )
  
  fm$gene       <- rownames(fm)
  fm$cell_type  <- ct        # tag with the originating identity
  fm
})

markers_all <- bind_rows(marker_list)

# ---------------------------------------------------------------------------
# 2) Wide matrices: log2FCs  and  P-values -----------------------------------
fc_mat <- markers_all %>%
  select(cell_type, gene, avg_log2FC) %>%
  pivot_wider(names_from = cell_type, values_from = avg_log2FC) %>%
  column_to_rownames("gene") %>%
  as.matrix()

p_mat <- markers_all %>%
  select(cell_type, gene, p_val_adj) %>%
  pivot_wider(names_from = cell_type, values_from = p_val_adj) %>%
  column_to_rownames("gene") %>%
  as.matrix()

# Replace missing values
fc_mat[is.na(fc_mat)] <- 0
p_mat[is.na(p_mat)]   <- 1   # treat as non-significant

# ---------------------------------------------------------------------------
# 3) Cap log2FC to ±2 --------------------------------------------------------
fc_mat[fc_mat >  1] <-  1
fc_mat[fc_mat < -1] <- -1

# ---------------------------------------------------------------------------
# 4) Build a “star” matrix for display_numbers -------------------------------
star_mat <- ifelse(p_mat < 0.01, "*", "")   # put star where adj-p < 0.01

# ---------------------------------------------------------------------------
# 5) Heat-map ----------------------------------------------------------------
p <- pheatmap(
  fc_mat,
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  color          = colorRampPalette(c("blue", "white", "red"))(100),
  main           = "EPO vs Placebo (avg_log2FC, ★ = adj P < 0.01)",
  display_numbers = star_mat,
  number_color   = "black",
  fontsize_row   = 12,
  fontsize_number = 12
)

save_highres_plot(p, filename = "Figure6/EPOvsPlacebo_TEproxGenes_log2FC_padj0.01_globalIdentities_heatmap", width = 8)













RNA <- readRDS("../ATAC_analysis_with_combined_peaks/RNA_final_(without_A3_A5_B10_B12).RDS")


head(RNA@meta.data)

RNA@meta.data$global_identity_treatment <- paste0(RNA@meta.data$global_identity, "_", RNA@meta.data$treatment)

Idents(RNA) <- "global_identity_treatment"
available_genes <- rownames(RNA)
valid_genes <- intersect(unique(da_peaks_df_with_peak_up$closest_gene_symbol), available_genes)
markers <- FindMarkers(
  object = RNA,
  test.use = "bimod",
  features = valid_genes,
  ident.1 = "Pyramidal neurons_Erythropoietin",
  ident.2 = "Pyramidal neurons_Placebo",
  logfc.threshold = 0,
  min.pct = 0
)

unique(te_dars_with_tf)

te_dars_unique <- unique(te_dars_with_tf)

# Turn the GRanges coordinates back into the same "chr:start-end" string
# that da_peaks_df$peak uses.
peak_ids <- paste0(seqnames(te_dars_unique), "-", start(te_dars_unique), "-", end(te_dars_unique))

# --- 2. Pull log2FCs and average duplicates ---------------------------------
heat_df <- da_peaks_df %>%
  filter(peak %in% peak_ids) %>%            # keep only the TF-overlapping peaks
  group_by(peak) %>%                       # collapse any duplicates
  summarise(avg_log2FC = mean(avg_log2FC, na.rm = TRUE), .groups = "drop")

# --- 2b. Cap avg_log2FC between –2 and +2 -----------------------------------
heat_df <- heat_df %>% 
  mutate(avg_log2FC = pmax(pmin(avg_log2FC,  2),  0))   # upper-cap at +2, lower-cap at –2


# --- 3. Build the matrix (now guaranteed unique row names) ------------------
heat_mat <- as.matrix(heat_df %>% column_to_rownames("peak"))
colnames(heat_mat) <- "avg_log2FC"

# --- 4. Optional ordering ----------------------------------------------------
heat_mat <- heat_mat[order(heat_mat[, 1], decreasing = TRUE), , drop = FALSE]

# --- 5. Plot -----------------------------------------------------------------
p <- pheatmap(
  heat_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color        = colorRampPalette(c("pink", "red2"))(100),
  main         = "TE-derived upregulated DARs overlapping any TF\n(The nearby genes are checked by snRNA-Seq)"
)

save_highres_plot(p, filename = "Figure6/TE-derived upregulated DARs overlapping any TF", width = 8)
markers
 
FindMarkers(RNA, test.use = "bimod", features = as.character(unique_gene_symbols), ident.1 = "Erythropoietin", ident.2 = "Placebo")
FindMarkers(RNA, test.use = "bimod", ident.1 = "Erythropoietin", ident.2 = "Placebo")

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(dplyr)
library(ggtext)

# Ensure symbols are unique and non-missing
gene_symbols_clean <- unique(na.omit(unique_gene_symbols))

# Run GO enrichment for all 3 ontologies
ego_BP <- enrichGO(
  gene          = gene_symbols_clean,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_CC <- enrichGO(
  gene          = gene_symbols_clean,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "CC",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

ego_MF <- enrichGO(
  gene          = gene_symbols_clean,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",
  ont           = "MF",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
# Extract top terms (20 per ontology) and label ontology type
df_BP <- as.data.frame(ego_BP)[1:20, ] %>% mutate(ontology = "Biological Process")
df_CC <- as.data.frame(ego_CC)[1:20, ] %>% mutate(ontology = "Cellular Component")
df_MF <- as.data.frame(ego_MF)[1:19, ] %>% mutate(ontology = "Molecular Function")
 
# Combine and clean
go_df <- bind_rows(df_BP, df_CC, df_MF)

# Format factor levels for y-axis ordering within each ontology
go_df <- go_df %>%
  group_by(ontology) %>%
  mutate(Description = factor(Description, levels = rev(unique(Description)))) %>%
  ungroup()

# Plot
p <- ggplot(go_df, aes(x = Description, y = Count, fill = ontology)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(
    aes(label = round(-log10(p.adjust), 1)),
    hjust = -0.2, size = 3.5, color = "black"
  ) +
  facet_wrap(~ ontology, scales = "free_y", nrow = 1) +  # 3 side-by-side panels
  scale_fill_manual(values = c(
    "Biological Process" = "#1f78b4",
    "Cellular Component" = "#33a02c",
    "Molecular Function" = "#ff7f00"
  )) +
  coord_flip() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(
    title = "GO Enrichment for Genes Near DARs Overlapping with TFs",
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
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none",  # hide legend since each facet already shows ontology
    plot.title = element_text(color = "black", face = "bold", size = 16, hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(color = "black", size = 12, hjust = 0.5),
    panel.grid = element_blank()
  )

# Save
save_highres_plot(p, filename = "Figure6/GO_enrichment_DARs_overlapping_with_TFs", width = 24, height = 8)








library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

# Ensure symbols are valid
gene_symbols_clean <- unique(na.omit(unique_gene_symbols))

# Run GO enrichment directly with SYMBOLs
ego <- enrichGO(
  gene          = gene_symbols_clean,
  OrgDb         = org.Mm.eg.db,
  keyType       = "SYMBOL",         # ✅ use SYMBOL directly
  ont           = "BP",             # Biological Process
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Plot results
if (!is.null(ego) && nrow(ego@result) > 0) {
  print(dotplot(ego, showCategory = 40, title = "GO Enrichment for Genes Near TE-derived DARs"))
} else {
  message("⚠️ No significant GO terms found.")
}
  







library(readr)   # fast CSV reader
library(dplyr)   # data-manipulation verbs

# 1. Read the file
peaks <- read_csv(
  "Data/Pyramidal/EPO_vs_Placebo/Newly_formed_vs_Mature_DA_peaks_combined_pyramidal_celltypes_EPO_vs_Placebo_separation_filtered_annotated_for_genes_and_TEs.csv",
  show_col_types = FALSE        # suppress column-type messages
)

# 2. Filter for the desired comparison
nf_only <- peaks %>%
  filter(comparison == "newly_formed_EPO_vs_newly_formed_Placebo")

library(dplyr)

nf_TE_consistent <- nf_only %>% 
  group_by(TE_family) %>% 
  # keep a TE_name group only if all log2FCs are either >0 or <0
  filter(all(avg_log2FC > 0) | all(avg_log2FC < 0)) %>% 
  ungroup()

# quick sanity-check
nf_TE_consistent %>% 
  count(TE_name, name = "rows_per_TE") %>% 
  print(n = 100)






















readRDS("../ATAC_analysis_with_combined_peaks/ATAC_with_sharedUMAP_and_imputedRNA.rds") -> ATAC_test

ATAC <- readRDS("../ATAC_analysis_with_combined_peaks/ATAC_with_sharedUMAP_and_imputedRNA.rds")
# Tabulate and sort predicted.id counts
cell_order <- names(sort(table(ATAC@meta.data$predicted.id), decreasing = TRUE))
ATAC@meta.data$predicted.id <- factor(ATAC@meta.data$predicted.id, levels = cell_order)

DefaultAssay(ATAC) <- "peaks"
library(Signac)
library(GenomicRanges)
library(ggplot2)

# 1. Define the top peak region
region <- "chr9-53300903-53303235"


# 4. Plot with TE annotation
p <- CoveragePlot(
  object = ATAC,
  region = region,
  features = NULL,
  extend.upstream = 1000,
  extend.downstream = 1000,
  group.by = "predicted.id",  # or your preferred grouping
  peaks = TRUE, 
  assay = "peaks",
  #idents = "Pyramidal neurons"
)

# 5. Show plot
print(p)

library(AnnotationHub)
ah = AnnotationHub()
annotations <- ah[["AH99012"]]
exclude_classes <- c("Low_complexity", "Satellite", "rRNA", "tRNA", "scRNA", "snRNA", 
                     "RNA", "RC", "srpRNA", "DNA?", "LTR?", "LINE?", "RC?", "SINE?") #Simple_repeat is specific to here
annotations <- annotations[
  !(annotations$repClass %in% exclude_classes | grepl("\\?$", annotations$repClass))
]
strand(annotations) <- "*"

head(annotations)
genome(annotations) <- "mm10"
# Set required metadata for TE annotations
annotations$type <- "exon"  # Needed for drawing in CoveragePlot
annotations$gene_name <- paste0(annotations$repName, "_TE_", seq_along(annotations))
annotations$gene_id <- annotations$gene_name  # Can be the same for simplicity
annotations$tx_id <- NA  # Optional if missing
annotations$gene_biotype <- annotations$repClass  # Optional but useful

# Set genome if missing
GenomeInfoDb::genome(annotations) <- "mm10"

common_seqlevels <- intersect(
  GenomeInfoDb::seqlevels(annotations),
  GenomeInfoDb::seqlevels(Annotation(ATAC))
)

annotations <- GenomeInfoDb::keepSeqlevels(
  annotations,
  value = common_seqlevels,
  pruning.mode = "coarse"
)

combined_annotations <- c(Annotation(ATAC), annotations)
Annotation(ATAC) <- combined_annotations

# 4. Plot with TE annotation
CoveragePlot(
  object = ATAC,
  region = region,
  features = NULL,
  extend.upstream = 1000,
  extend.downstream = 1000,
  group.by = "predicted.id",  # or your preferred grouping
  peaks = TRUE, 
  assay = "peaks",
  #idents = "Pyramidal neurons"
)


regions_full <- c(
  "chr15-79798080-79799685",  # Nptxr
  "chr4-154956456-154959454", # Hes5
  "chr4-154956456-154959454", # Pank4
  "chr4-41028935-41030559",   # Nfx1
  "chr11-57520630-57523177",  # Fam114a2
  "chr11-57520630-57523177",  # Mfap3
  "chr1-171259815-171264537", # B4galt3
  "chr12-111666114-111669956",# Ckb
  "chr5-110597356-110599454", # Galnt9
  "chr7-133767706-133769095", # Dhx32
  "chr7-133767706-133769095", # Fank1
  "chr6-115890559-115892228", # Ift122
  "chr17-27655274-27658086",  # Pacsin1
  "chr3-122478475-122481415"  # Bcar3
)

CoverageBrowser(ATAC, region, assay = "peaks", sep = c("-", "-"))

# Loop through regions and make coverage plots
for (region in regions_full) {
  p <- CoveragePlot(
    object = ATAC,
    region = region,
    features = NULL,
    extend.upstream = 1000,
    extend.downstream = 1000,
    group.by = "predicted.id",  # adjust grouping if needed
    peaks = TRUE,
    assay = "peaks"
    # idents = "Pyramidal neurons"  # optional
  ) & scale_fill_manual(values = palette_umap)
  # Save each plot
  save_highres_plot(p,
    filename = paste0("Figure6/igvplots/+-1kb region/coverage_", gsub("[:-]", "_", region)),
    width = 10,
    height = 5
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



Ckb <- CoveragePlot(
  object = ATAC,
  region = "chr12-111666114-111669956",
  features = NULL,
  extend.upstream = 2000,
  extend.downstream = 2000,
  group.by = "predicted.id",  # adjust grouping if needed
  peaks = TRUE,
  assay = "peaks",
  cols = palette_umap  
  # idents = "Pyramidal neurons"  # optional
) & scale_fill_manual(values = palette_umap)

save_highres_plot(Ckb, filename = "Figure6/igvplots/Ckb")

Galnt9 <- CoveragePlot(
  object = ATAC,
  region = "chr5-110597356-110599454",
  features = NULL,
  extend.upstream = 2000,
  extend.downstream = 2000,
  group.by = "predicted.id",  # adjust grouping if needed
  peaks = TRUE,
  assay = "peaks"
  # idents = "Pyramidal neurons"  # optional
) & scale_fill_manual(values = palette_umap)

save_highres_plot(Galnt9, filename = "Figure6/igvplots/Galnt9")

Nfx1 <- CoveragePlot(
  object = ATAC,
  region = "chr4-41028935-41030559",
  features = NULL,
  extend.upstream = 5000,
  extend.downstream = 1000,
  group.by = "predicted.id",  # adjust grouping if needed
  peaks = TRUE,
  assay = "peaks"
  # idents = "Pyramidal neurons"  # optional
) & scale_fill_manual(values = palette_umap)

save_highres_plot(Nfx1, filename = "Figure6/igvplots/Nfx1")


Bcar <- CoveragePlot(
  object = ATAC,
  region = "chr3-122478475-122481415",
  features = NULL,
  extend.upstream = 2000,
  extend.downstream = 2000,
  group.by = "predicted.id",  # adjust grouping if needed
  peaks = TRUE,
  assay = "peaks"
  # idents = "Pyramidal neurons"  # optional
) & scale_fill_manual(values = palette_umap)

save_highres_plot(Bcar, filename = "Figure6/igvplots/Bcar")


Pacsin1 <- CoveragePlot(
  object = ATAC,
  region = "chr17-27655274-27658086",
  features = NULL,
  extend.upstream = 1000,
  extend.downstream = 2000,
  group.by = "predicted.id",  # adjust grouping if needed
  peaks = TRUE,
  assay = "peaks"
  # idents = "Pyramidal neurons"  # optional
) & scale_fill_manual(values = palette_umap)
  
save_highres_plot(Pacsin1, filename = "Figure6/igvplots/Pacsin1")


 
Ift22 <- CoveragePlot(
  object = ATAC,
  region = "chr6-115890559-115892228",
  features = NULL,
  extend.upstream = 500,
  extend.downstream = 500,
  group.by = "predicted.id",  # adjust grouping if needed
  peaks = TRUE,
  assay = "peaks"
  # idents = "Pyramidal neurons"  # optional
)  & scale_fill_manual(values = palette_umap)

save_highres_plot(Ift22, filename = "Figure6/igvplots/Ift22")
  

B4galt3 <- CoveragePlot(
  object = ATAC,
  region = "chr1-171259815-171264537",
  features = NULL,
  extend.upstream = 500,
  extend.downstream = 10000,
  group.by = "predicted.id",  # adjust grouping if needed
  peaks = TRUE,
  assay = "peaks"
  # idents = "Pyramidal neurons"  # optional
) & scale_fill_manual(values = palette_umap)

save_highres_plot(B4galt3, filename = "Figure6/igvplots/B4galt3")



Mfap3 <- CoveragePlot(
  object = ATAC,
  region = "chr11-57520630-57523177",
  features = NULL,
  extend.upstream = 1000,
  extend.downstream = 1000,
  group.by = "predicted.id",  # adjust grouping if needed
  peaks = TRUE,
  assay = "peaks"
  # idents = "Pyramidal neurons"  # optional
) & scale_fill_manual(values = palette_umap)

save_highres_plot(Mfap3, filename = "Figure6/igvplots/Mfap3")












ATAC_newly_formed_mature <- readRDS("Data/Pyramidal/EPO_vs_Placebo/Newly_formed_vs_Mature_ATAC_individual_peaks_for_Pyramidal_neurons_mature_vs_newly_formed.RDS")



ATAC_11_lineages <- readRDS("../ATAC_analysis_with_combined_peaks/ATAC_with_sharedUMAP_and_imputedRNA.rds")
DefaultAssay(ATAC_11_lineages) <- "peaks"
ATAC_11_lineages@meta.data <- ATAC_11_lineages@meta.data %>%
  mutate(
    condition = case_when(
      startsWith(orig.ident, "A") ~ "EPO",
      startsWith(orig.ident, "B") ~ "Placebo",
      TRUE ~ NA_character_
    )
  )
ATAC_11_lineages
head(ATAC_11_lineages)


# ── 0) packages ─────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)

# ── 1) pick columns & optional QC filter ────────────────────
meta <- ATAC_11_lineages@meta.data

sample_col   <- if ("orig.ident"   %in% colnames(meta)) "orig.ident"   else "dataset"
celltype_col <- if ("predicted.id" %in% colnames(meta)) "predicted.id" else "seurat_clusters"

meta_use <- meta %>%
  { if ("is__cell_barcode" %in% names(.)) filter(., is__cell_barcode == 1) else . } %>%
  { if ("excluded_reason" %in% names(.))   filter(., is.na(excluded_reason) | excluded_reason == 0) else . } %>%
  select(sample = all_of(sample_col), celltype = all_of(celltype_col)) %>%
  mutate(sample  = as.factor(sample),
         celltype = as.factor(celltype))

# ── 2) counts and normalized composition per animal ────────
comp_long <- meta_use %>%
  count(sample, celltype, name = "n") %>%
  group_by(sample) %>%
  mutate(total_cells = sum(n),
         prop        = n / total_cells,      # 0–1
         pct         = 100 * prop,           # 0–100
         per_1k      = 1000 * prop) %>%      # per 1000 cells (sometimes handy)
  ungroup()

# nice cell-type order (by overall abundance)
ct_order <- comp_long %>%
  group_by(celltype) %>% summarise(mean_prop = mean(prop)) %>%
  arrange(desc(mean_prop)) %>% pull(celltype)
comp_long <- comp_long %>%
  mutate(celltype = factor(celltype, levels = ct_order),
         sample   = fct_inorder(sample))

# ── 3) wide tables (optional saves) ─────────────────────────
# proportions (rows = sample, cols = cell types), sums to 1 per row
comp_prop_wide <- comp_long %>%
  select(sample, celltype, prop) %>%
  tidyr::pivot_wider(names_from = celltype, values_from = prop, values_fill = 0)

# raw counts
comp_count_wide <- comp_long %>%
  select(sample, celltype, n) %>%
  pivot_wider(names_from = celltype, values_from = n, values_fill = 0)

write.csv(comp_prop_wide, "extra_analysis/composition/composition_proportions_by_animal.csv", row.names = FALSE)
write.csv(comp_count_wide, "extra_analysis/composition/composition_counts_by_animal.csv", row.names = FALSE)


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

# ── 4) plots ────────────────────────────────────────────────
# (a) stacked bar: fraction per cell type within each animal
# 2) Stacked bar with your palette, legend & stack follow ct_order
p_bar <- ggplot(comp_long, aes(x = sample, y = prop, fill = celltype)) +
  geom_col(width = 0.85) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(
    values = palette_umap,   # maps by name
    breaks = ct_order,       # legend order = your ct_order
    drop   = TRUE,
    na.value = "grey70"
  ) +
  labs(x = "Animal", y = "Composition (fraction of cells)", fill = "Cell type",
       title = "Per-animal cell-type composition (normalized within animal)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

# ggsave("composition_stacked_bar.pdf", p_bar, width = 9, height = 5)

save_highres_plot(p_bar, filename = "extra_analysis/composition/composition_stacked_bar")

# (b) heatmap: animals × cell types, values are proportions
p_heat <- ggplot(comp_long, aes(x = celltype, y = sample, fill = prop)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue", labels = percent_format(accuracy = 1)) +
  labs(x = "Cell type", y = "Animal", fill = "Fraction",
       title = "Per-animal cell-type composition") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))
# ggsave("composition_heatmap.pdf", p_heat, width = 10, height = 5)

save_highres_plot(p_heat, filename = "extra_analysis/composition/composition_heatmap")














 

set.seed(0)
ATAC_11_lineages <- readRDS("../ATAC_analysis_with_combined_peaks/ATAC_with_sharedUMAP_and_imputedRNA.rds")
DefaultAssay(ATAC_11_lineages) <- "peaks"
ATAC_11_lineages@meta.data <- ATAC_11_lineages@meta.data %>%
  mutate(
    condition = case_when(
      startsWith(orig.ident, "A") ~ "EPO",
      startsWith(orig.ident, "B") ~ "Placebo",
      TRUE ~ NA_character_
    )
  )







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



meta <- ATAC_11_lineages@meta.data

# identify columns
sample_col   <- if ("orig.ident"   %in% colnames(meta)) "orig.ident"   else "dataset"
celltype_col <- if ("predicted.id" %in% colnames(meta)) "predicted.id" else "seurat_clusters"
cond_col <- if ("Condition" %in% colnames(meta)) {
  "Condition"
} else if ("condition" %in% colnames(meta)) {
  "condition"
} else {
  stop("No Condition/condition column found in metadata.")
}

# keep everything; just select needed columns
meta_use <- meta %>%
  dplyr::select(sample    = dplyr::all_of(sample_col),
                celltype  = dplyr::all_of(celltype_col),
                condition = dplyr::all_of(cond_col)) %>%
  dplyr::mutate(dplyr::across(c(sample, celltype, condition), as.factor))

# ── per-animal composition ─────────────────────────────────
comp_long <- meta_use %>%
  dplyr::count(sample, celltype, name = "n") %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(
    total_cells = sum(n),
    prop        = n / total_cells,
    pct         = 100 * prop,
    per_1k      = 1000 * prop
  ) %>%
  dplyr::ungroup()

# abundance-based order (your rule)
ct_order <- comp_long %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(mean_prop = mean(prop), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(mean_prop)) %>%
  dplyr::pull(celltype)

comp_long <- comp_long %>%
  dplyr::mutate(
    celltype = factor(celltype, levels = ct_order),
    sample   = forcats::fct_inorder(sample)
  )

# add condition per sample (no exclusions)
sample_condition <- meta_use %>% dplyr::select(sample, condition) %>% dplyr::distinct()
comp_long <- dplyr::left_join(comp_long, sample_condition, by = "sample")

# wide tables (animals × cell types)
comp_prop_wide <- comp_long %>%
  dplyr::select(sample, celltype, prop) %>%
  tidyr::pivot_wider(names_from = celltype, values_from = prop, values_fill = 0)

comp_count_wide <- comp_long %>%
  dplyr::select(sample, celltype, n) %>%
  tidyr::pivot_wider(names_from = celltype, values_from = n, values_fill = 0)

dir.create("extra_analysis/composition", recursive = TRUE, showWarnings = FALSE)
write.csv(comp_prop_wide, "extra_analysis/composition/composition_proportions_by_animal.csv", row.names = FALSE)
write.csv(comp_count_wide, "extra_analysis/composition/composition_counts_by_animal.csv", row.names = FALSE)

cond_colors <- c("EPO" = "red", "Placebo" = "blue")

# ── plots: per-animal ──────────────────────────────────────
p_bar <- ggplot2::ggplot(comp_long, ggplot2::aes(x = sample, y = prop, fill = celltype)) +
  ggplot2::geom_col(width = 0.85) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_fill_manual(values = palette_umap, breaks = ct_order, drop = TRUE) +
  ggplot2::labs(x = "Animal", y = "Composition (fraction of cells)",
                fill = "Cell type",
                title = "Per-animal cell-type composition (normalized within animal)") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(legend.position = "right")


save_highres_plot(p_bar, filename = "extra_analysis/composition/composition_stacked_bar")
# ggsave("extra_analysis/composition/composition_stacked_bar.pdf", p_bar, width = 9, height = 5)

p_heat <- ggplot2::ggplot(comp_long, ggplot2::aes(x = celltype, y = sample, fill = prop)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient(low = "white", high = "steelblue",
                               labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(x = "Cell type", y = "Animal", fill = "Fraction",
                title = "Per-animal cell-type composition") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))

save_highres_plot(p_heat, filename = "extra_analysis/composition/composition_heatmap")
# ggsave("extra_analysis/composition/composition_heatmap.pdf", p_heat, width = 10, height = 5)

# ── per-condition composition (EPO vs Placebo) ─────────────
comp_long_cond <- meta_use %>%
  dplyr::count(condition, celltype, name = "n") %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(
    total_cells = sum(n),
    prop        = n / total_cells,
    pct         = 100 * prop,
    per_1k      = 1000 * prop
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    celltype  = factor(celltype, levels = ct_order),
    condition = forcats::fct_inorder(condition)   # preserves metadata order
  )

# wide tables (conditions × cell types)
comp_prop_wide_cond <- comp_long_cond %>%
  dplyr::select(condition, celltype, prop) %>%
  tidyr::pivot_wider(names_from = celltype, values_from = prop, values_fill = 0)

comp_count_wide_cond <- comp_long_cond %>%
  dplyr::select(condition, celltype, n) %>%
  tidyr::pivot_wider(names_from = celltype, values_from = n, values_fill = 0)

dir.create("extra_analysis/composition_condition", recursive = TRUE, showWarnings = FALSE)
write.csv(comp_prop_wide_cond,  "extra_analysis/composition_condition/composition_proportions_by_condition.csv", row.names = FALSE)
write.csv(comp_count_wide_cond, "extra_analysis/composition_condition/composition_counts_by_condition.csv",      row.names = FALSE)

# bars: stack by cell type, colored outline by condition (EPO=red, Placebo=blue)
p_bar_cond <- ggplot2::ggplot(comp_long_cond,
                              ggplot2::aes(x = condition, y = prop, fill = celltype)) +
  ggplot2::geom_col(width = 0.85) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_fill_manual(values = palette_umap, breaks = ct_order, drop = TRUE) +
  ggplot2::labs(x = "Condition", y = "Composition (fraction of cells)",
                fill = "Cell type",
                title = "Per-condition cell-type composition (normalized within condition)") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(legend.position = "right")

save_highres_plot(p_bar_cond, filename = "extra_analysis/composition_condition/composition_stacked_bar_condition")
# ggsave("extra_analysis/composition_condition/composition_stacked_bar_condition.pdf", p_bar_cond, width = 9, height = 5)

# heatmap: conditions × cell types
p_heat_cond <- ggplot2::ggplot(comp_long_cond, ggplot2::aes(x = celltype, y = condition, fill = prop)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient(low = "white", high = "steelblue",
                               labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(x = "Cell type", y = "Condition", fill = "Fraction",
                title = "Per-condition cell-type composition") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))

save_highres_plot(p_heat_cond, filename = "extra_analysis/composition_condition/composition_heatmap_condition")
# ggsave("extra_analysis/composition_condition/composition_heatmap_condition.pdf", p_heat_cond, width = 10, height = 5)
















set.seed(0)
ATAC_newly_formed_mature <- readRDS("Data/Pyramidal/EPO_vs_Placebo/Newly_formed_vs_Mature_ATAC_individual_peaks_for_Pyramidal_neurons_mature_vs_newly_formed.RDS")
DefaultAssay(ATAC_newly_formed_mature) <- "peaks"
ATAC_newly_formed_mature@meta.data <- ATAC_newly_formed_mature@meta.data %>%
  mutate(
    condition = case_when(
      startsWith(orig.ident, "A") ~ "EPO",
      startsWith(orig.ident, "B") ~ "Placebo",
      TRUE ~ NA_character_
    )
  )


meta <- ATAC_newly_formed_mature@meta.data

# identify columns
sample_col   <- if ("orig.ident"   %in% colnames(meta)) "orig.ident"   else "dataset"
celltype_col <- if ("mature_vs_newlyformed_status" %in% colnames(meta)) "mature_vs_newlyformed_status" else "mature_vs_newlyformed_status"
cond_col <- if ("Condition" %in% colnames(meta)) {
  "Condition"
} else if ("condition" %in% colnames(meta)) {
  "condition"
} else {
  stop("No Condition/condition column found in metadata.")
}

# keep everything; just select needed columns
meta_use <- meta %>%
  dplyr::select(sample    = dplyr::all_of(sample_col),
                celltype  = dplyr::all_of(celltype_col),
                condition = dplyr::all_of(cond_col)) %>%
  dplyr::mutate(dplyr::across(c(sample, celltype, condition), as.factor))

# ── per-animal composition ─────────────────────────────────
comp_long <- meta_use %>%
  dplyr::count(sample, celltype, name = "n") %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(
    total_cells = sum(n),
    prop        = n / total_cells,
    pct         = 100 * prop,
    per_1k      = 1000 * prop
  ) %>%
  dplyr::ungroup()

# abundance-based order (your rule)
ct_order <- comp_long %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(mean_prop = mean(prop), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(mean_prop)) %>%
  dplyr::pull(celltype)

comp_long <- comp_long %>%
  dplyr::mutate(
    celltype = factor(celltype, levels = ct_order),
    sample   = forcats::fct_inorder(sample)
  )

# add condition per sample (no exclusions)
sample_condition <- meta_use %>% dplyr::select(sample, condition) %>% dplyr::distinct()
comp_long <- dplyr::left_join(comp_long, sample_condition, by = "sample")

# wide tables (animals × cell types)
comp_prop_wide <- comp_long %>%
  dplyr::select(sample, celltype, prop) %>%
  tidyr::pivot_wider(names_from = celltype, values_from = prop, values_fill = 0)

comp_count_wide <- comp_long %>%
  dplyr::select(sample, celltype, n) %>%
  tidyr::pivot_wider(names_from = celltype, values_from = n, values_fill = 0)

dir.create("extra_analysis/composition", recursive = TRUE, showWarnings = FALSE)
write.csv(comp_prop_wide, "extra_analysis/composition/composition_proportions_by_animal.csv", row.names = FALSE)
write.csv(comp_count_wide, "extra_analysis/composition/composition_counts_by_animal.csv", row.names = FALSE)

cond_colors <- c("EPO" = "red", "Placebo" = "blue")
palette_umap <- c("newly_formed" = "#e31a1c", "mature" = "#1f78b4")


# ── plots: per-animal ──────────────────────────────────────
p_bar <- ggplot2::ggplot(comp_long, ggplot2::aes(x = sample, y = prop, fill = celltype)) +
  ggplot2::geom_col(width = 0.85) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_fill_manual(values = palette_umap, breaks = ct_order, drop = TRUE) +
  ggplot2::labs(x = "Animal", y = "Composition (fraction of cells)",
                fill = "Cell type",
                title = "Per-animal cell-type composition (normalized within animal)") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(legend.position = "right")


save_highres_plot(p_bar, filename = "extra_analysis/composition/composition_stacked_bar")
# ggsave("extra_analysis/composition/composition_stacked_bar.pdf", p_bar, width = 9, height = 5)

p_heat <- ggplot2::ggplot(comp_long, ggplot2::aes(x = celltype, y = sample, fill = prop)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient(low = "white", high = "steelblue",
                               labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(x = "Cell type", y = "Animal", fill = "Fraction",
                title = "Per-animal cell-type composition") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))

save_highres_plot(p_heat, filename = "extra_analysis/composition/composition_heatmap")
# ggsave("extra_analysis/composition/composition_heatmap.pdf", p_heat, width = 10, height = 5)

# ── per-condition composition (EPO vs Placebo) ─────────────
comp_long_cond <- meta_use %>%
  dplyr::count(condition, celltype, name = "n") %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(
    total_cells = sum(n),
    prop        = n / total_cells,
    pct         = 100 * prop,
    per_1k      = 1000 * prop
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    celltype  = factor(celltype, levels = ct_order),
    condition = forcats::fct_inorder(condition)   # preserves metadata order
  )

# wide tables (conditions × cell types)
comp_prop_wide_cond <- comp_long_cond %>%
  dplyr::select(condition, celltype, prop) %>%
  tidyr::pivot_wider(names_from = celltype, values_from = prop, values_fill = 0)

comp_count_wide_cond <- comp_long_cond %>%
  dplyr::select(condition, celltype, n) %>%
  tidyr::pivot_wider(names_from = celltype, values_from = n, values_fill = 0)

dir.create("extra_analysis/composition_condition", recursive = TRUE, showWarnings = FALSE)
write.csv(comp_prop_wide_cond,  "extra_analysis/composition_condition/composition_proportions_by_condition.csv", row.names = FALSE)
write.csv(comp_count_wide_cond, "extra_analysis/composition_condition/composition_counts_by_condition.csv",      row.names = FALSE)

# bars: stack by cell type, colored outline by condition (EPO=red, Placebo=blue)
p_bar_cond <- ggplot2::ggplot(comp_long_cond,
                              ggplot2::aes(x = condition, y = prop, fill = celltype)) +
  ggplot2::geom_col(width = 0.85) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_fill_manual(values = palette_umap, breaks = ct_order, drop = TRUE) +
  ggplot2::labs(x = "Condition", y = "Composition (fraction of cells)",
                fill = "Cell type",
                title = "Per-condition cell-type composition (normalized within condition)") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(legend.position = "right")

save_highres_plot(p_bar_cond, filename = "extra_analysis/composition_condition/composition_stacked_bar_condition")
# ggsave("extra_analysis/composition_condition/composition_stacked_bar_condition.pdf", p_bar_cond, width = 9, height = 5)

# heatmap: conditions × cell types
p_heat_cond <- ggplot2::ggplot(comp_long_cond, ggplot2::aes(x = celltype, y = condition, fill = prop)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient(low = "white", high = "steelblue",
                               labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(x = "Cell type", y = "Condition", fill = "Fraction",
                title = "Per-condition cell-type composition") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))

save_highres_plot(p_heat_cond, filename = "extra_analysis/composition_condition/composition_heatmap_condition")
# ggsave("extra_analysis/composition_condition/composition_heatmap_condition.pdf", p_heat_cond, width = 10, height = 5)















set.seed(0)
ATAC_newly_formed_mature <- readRDS("Data/Pyramidal/EPO_vs_Placebo/Newly_formed_vs_Mature_ATAC_individual_peaks_for_Pyramidal_neurons_mature_vs_newly_formed.RDS")
DefaultAssay(ATAC_newly_formed_mature) <- "peaks"
ATAC_newly_formed_mature@meta.data <- ATAC_newly_formed_mature@meta.data %>%
  mutate(
    condition = case_when(
      startsWith(orig.ident, "A") ~ "EPO",
      startsWith(orig.ident, "B") ~ "Placebo",
      TRUE ~ NA_character_
    )
  )


meta <- ATAC_newly_formed_mature@meta.data

# identify columns
sample_col   <- if ("orig.ident"   %in% colnames(meta)) "orig.ident"   else "dataset"
celltype_col <- if ("transfered_predicted.id2" %in% colnames(meta)) "transfered_predicted.id2" else "transfered_predicted.id2"
cond_col <- if ("Condition" %in% colnames(meta)) {
  "Condition"
} else if ("condition" %in% colnames(meta)) {
  "condition"
} else {
  stop("No Condition/condition column found in metadata.")
}

# keep everything; just select needed columns
meta_use <- meta %>%
  dplyr::select(sample    = dplyr::all_of(sample_col),
                celltype  = dplyr::all_of(celltype_col),
                condition = dplyr::all_of(cond_col)) %>%
  dplyr::mutate(dplyr::across(c(sample, celltype, condition), as.factor))

# ── per-animal composition ─────────────────────────────────
comp_long <- meta_use %>%
  dplyr::count(sample, celltype, name = "n") %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(
    total_cells = sum(n),
    prop        = n / total_cells,
    pct         = 100 * prop,
    per_1k      = 1000 * prop
  ) %>%
  dplyr::ungroup()

# abundance-based order (your rule)
ct_order <- comp_long %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(mean_prop = mean(prop), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(mean_prop)) %>%
  dplyr::pull(celltype)

comp_long <- comp_long %>%
  dplyr::mutate(
    celltype = factor(celltype, levels = ct_order),
    sample   = forcats::fct_inorder(sample)
  )

# add condition per sample (no exclusions)
sample_condition <- meta_use %>% dplyr::select(sample, condition) %>% dplyr::distinct()
comp_long <- dplyr::left_join(comp_long, sample_condition, by = "sample")

# wide tables (animals × cell types)
comp_prop_wide <- comp_long %>%
  dplyr::select(sample, celltype, prop) %>%
  tidyr::pivot_wider(names_from = celltype, values_from = prop, values_fill = 0)

comp_count_wide <- comp_long %>%
  dplyr::select(sample, celltype, n) %>%
  tidyr::pivot_wider(names_from = celltype, values_from = n, values_fill = 0)

dir.create("extra_analysis/composition", recursive = TRUE, showWarnings = FALSE)
write.csv(comp_prop_wide, "extra_analysis/composition/composition_proportions_by_animal.csv", row.names = FALSE)
write.csv(comp_count_wide, "extra_analysis/composition/composition_counts_by_animal.csv", row.names = FALSE)

cond_colors <- c("EPO" = "red", "Placebo" = "blue")
#palette_umap <- c("newly_formed" = "#e31a1c", "mature" = "#1f78b4")

set.seed(123)  # for reproducibility
celltype_colors <- DiscretePalette(n = length(unique(comp_long$celltype)), palette = "alphabet2")
celltype_colors[celltype_colors == "#E2E2E2"] <- "#194A8D"  # A rich blue-gray 
palette_umap <- celltype_colors

# ── plots: per-animal ──────────────────────────────────────
p_bar <- ggplot2::ggplot(comp_long, ggplot2::aes(x = sample, y = prop, fill = celltype)) +
  ggplot2::geom_col(width = 0.85) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_fill_manual(values = palette_umap, breaks = ct_order, drop = TRUE) +
  ggplot2::labs(x = "Animal", y = "Composition (fraction of cells)",
                fill = "Cell type",
                title = "Per-animal cell-type composition (normalized within animal)") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(legend.position = "right")


save_highres_plot(p_bar, filename = "extra_analysis/composition/composition_stacked_bar")
# ggsave("extra_analysis/composition/composition_stacked_bar.pdf", p_bar, width = 9, height = 5)

p_heat <- ggplot2::ggplot(comp_long, ggplot2::aes(x = celltype, y = sample, fill = prop)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient(low = "white", high = "steelblue",
                               labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(x = "Cell type", y = "Animal", fill = "Fraction",
                title = "Per-animal cell-type composition") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))

save_highres_plot(p_heat, filename = "extra_analysis/composition/composition_heatmap")
# ggsave("extra_analysis/composition/composition_heatmap.pdf", p_heat, width = 10, height = 5)

# ── per-condition composition (EPO vs Placebo) ─────────────
comp_long_cond <- meta_use %>%
  dplyr::count(condition, celltype, name = "n") %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(
    total_cells = sum(n),
    prop        = n / total_cells,
    pct         = 100 * prop,
    per_1k      = 1000 * prop
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    celltype  = factor(celltype, levels = ct_order),
    condition = forcats::fct_inorder(condition)   # preserves metadata order
  )

# wide tables (conditions × cell types)
comp_prop_wide_cond <- comp_long_cond %>%
  dplyr::select(condition, celltype, prop) %>%
  tidyr::pivot_wider(names_from = celltype, values_from = prop, values_fill = 0)

comp_count_wide_cond <- comp_long_cond %>%
  dplyr::select(condition, celltype, n) %>%
  tidyr::pivot_wider(names_from = celltype, values_from = n, values_fill = 0)

dir.create("extra_analysis/composition_condition", recursive = TRUE, showWarnings = FALSE)
write.csv(comp_prop_wide_cond,  "extra_analysis/composition_condition/composition_proportions_by_condition.csv", row.names = FALSE)
write.csv(comp_count_wide_cond, "extra_analysis/composition_condition/composition_counts_by_condition.csv",      row.names = FALSE)

# bars: stack by cell type, colored outline by condition (EPO=red, Placebo=blue)
p_bar_cond <- ggplot2::ggplot(comp_long_cond,
                              ggplot2::aes(x = condition, y = prop, fill = celltype)) +
  ggplot2::geom_col(width = 0.85) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_fill_manual(values = palette_umap, breaks = ct_order, drop = TRUE) +
  ggplot2::labs(x = "Condition", y = "Composition (fraction of cells)",
                fill = "Cell type",
                title = "Per-condition cell-type composition (normalized within condition)") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(legend.position = "right")

save_highres_plot(p_bar_cond, filename = "extra_analysis/composition_condition/composition_stacked_bar_condition")
# ggsave("extra_analysis/composition_condition/composition_stacked_bar_condition.pdf", p_bar_cond, width = 9, height = 5)

# heatmap: conditions × cell types
p_heat_cond <- ggplot2::ggplot(comp_long_cond, ggplot2::aes(x = celltype, y = condition, fill = prop)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient(low = "white", high = "steelblue",
                               labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(x = "Cell type", y = "Condition", fill = "Fraction",
                title = "Per-condition cell-type composition") +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 40, hjust = 1))

save_highres_plot(p_heat_cond, filename = "extra_analysis/composition_condition/composition_heatmap_condition")
# ggsave("extra_analysis/composition_condition/composition_heatmap_condition.pdf", p_heat_cond, width = 10, height = 5)


# For ChIP-Seq: 

library(tidyverse)

# ---- read the raw lines ----
lines <- read_lines("positive_negative_repeats.txt")

# ---- parse with regex ----
df <- tibble(raw = lines) %>%
  extract(
    raw,
    into = c("family", "TF", "positive", "negative"),
    regex = "^(\\S+) \\+ (\\S+): (\\d+) positive, (\\d+) negative repeats$",
    convert = TRUE
  )


write_tsv(df, "positive_negative_repeats.tsv")
cat("\nSaved table: positive_negative_repeats.tsv\n")
