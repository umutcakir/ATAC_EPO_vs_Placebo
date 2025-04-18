alldata <- readRDS("ATAC_with_individual_peaks_annotated_filtered_motifs_added.RDS")
set.seed(0)
library(Signac)
library(Seurat)
library(JASPAR2020)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)

Idents(alldata) <- "seurat_clusters"

print("Motifs are added") 

clusters <- unique(Idents(alldata))

library(future)
library(future.apply)
library(progressr)

# Set up parallel backend (adjust workers as needed)
plan("multisession", workers = 10)

options(future.globals.maxSize = 20 * 1024^3)  # 20 GB

with_progress({
  p <- progressor(along = clusters)
  
  da_peaks_list <- future_lapply(clusters, function(clust) {
    p(sprintf("Processing cluster %s", clust))
    
    da <- FindMarkers(
      object = alldata,
      ident.1 = clust,
      only.pos = TRUE,
      min.pct = 0.1,
      test.use = "LR",
      latent.vars = "nCount_peaks"
    )
    
    if (nrow(da) > 0) {
      da$cluster <- paste0("Cluster_", clust)
      da$peak_coordinate <- rownames(da)
      return(da)
    } else {
      return(NULL)
    }
  }, future.seed = TRUE)
})

# Clean NULLs from list
da_peaks_list <- Filter(Negate(is.null), da_peaks_list)
da_peaks_by_cluster <- da_peaks_list

# Combine into a single data frame
combined_da_peaks <- do.call(rbind, da_peaks_list)

# Save to CSV
write.csv(combined_da_peaks, "da_peaks_all_clusters_without_treatment_separation.csv", row.names = TRUE)

saveRDS(da_peaks_by_cluster, file = "DA_peaks_by_cluster_without_treatment_separation.rds")



motif_results_by_cluster <- list()

for (clust_name in names(da_peaks_by_cluster)) {
  da_table <- da_peaks_by_cluster[[clust_name]]
  
  sig_peaks <- rownames(da_table[da_table$p_val_adj < 0.05, ])
  sig_peaks <- sig_peaks[grepl("^chr", sig_peaks)]  # keep only reference
  
  if (length(sig_peaks) > 10) {
    message("Running motif enrichment for ", clust_name)
    
    motif_enrichment <- FindMotifs(
      object = alldata,
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

# Save to single CSV
write.csv(combined_motifs, "motif_enrichment_all_clusters_without_treatment_separation.csv", row.names = FALSE)

saveRDS(motif_results_by_cluster, file = "Motif_enrichment_by_cluster_without_treatment_separation.rds")



