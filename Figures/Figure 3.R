
dat <- read.delim("Diff_peak_scATACseq_EPO_Pl.bed", header=F)
row.names(dat) <- paste(dat[,1], dat[,2], dat[,3], sep="_")
dat <- dat[,c(4,6)]
dat[,1] <- -log10(dat[,1])

colnames(dat) <- c("log10P.adjusted", "avg_logFC")
attach(dat)

  png("Volcano_EPO_PL_DEATACseq_singlecell.png",   width = 10.7, height = 14.8, units = "cm", res = 600, pointsize = 12)
plot(avg_log2FC*-1, logp, pch=19, cex=1, 
                     ylab="-log10 p-value", xlab="Log2 Fold Change", xlim=c(-2,2),ylim=c(0,300),
 col=ifelse(avg_log2FC >= 0 & logp > 2,
                     "darkslateblue", ifelse(avg_log2FC < 0 & logp > 2,
 "darkred", "bisque4")))
abline(h = 3, col = "grey20", lty = 1, lwd = 1)
abline(v = c(-0.25,0.25), col = "grey20", lty = 2, lwd = 2)
 dev.off()

 dat <- read.csv("MotifsBvsA.epo.csv", header=T, sep=",")

 
  png("Volcano_EPO_PL_Motifs_DEATACseq_singlecell.png",   width = 10.7, height = 14.8, units = "cm", res = 600, pointsize = 12)
plot(fold.enrichment, logP, pch=19, cex=1, 
                     ylab="-log10 p-value", xlab="Log2 Fold Change", xlim=c(1,4),ylim=c(0,180),
 col=ifelse(fold.enrichment >= 0 & logP > 2,
                     "darkred", ifelse(fold.enrichment < 0 & logP > 2,
 "blue", "bisque4")))
abline(h = 3, col = "grey20", lty = 1, lwd = 1)
abline(v = c(0,0.25), col = "grey20", lty = 2, lwd = 2)
 dev.off()

savehistory("ATAC_seq_analysis_EPO_Riki.Rhistory")

load("ATAC_seq_analysis_EPO_Riki.RData")
ls()
head(dat)
head(data)
dim(data)
list.files()
da <- read.delim("EPO_Diff_peaks_over_TSS_genes.bed", header=F, stringsAsFactors=F)

da <- da[!duplicated(da[,7]),]
row.names(da) <- da[,7]

da <- da[,c(4:5)]

ge <- read.delim("DARs/DEG_EPO_Pl_Pyr_snRNAseq_hippo.tsv", row.names=1, header=T, stringsAsFactors=F)

dat <- merge(da, ge, by="row.names", all=F)

ka <- read.delim("Sorted_Diff_peak_scATACseq_EPO_Pl_neighbor_genes.bed", header=F, stringsAsFactors=F)

row.names(ka) <- paste(ka$V1, ka$V2, ka$V3, sep = "-")
ka$names <- paste(ka$V1, ka$V2, ka$V3, sep = "-")

ka <- ka[which(!duplicated(ka$names)),]

row.names(ka) <- ka$names
gdata <- merge(data, ka, by="row.names", all=F)

dat <- gdata[,c(1,3,6,16)]

attach(dat)
 png("Scatterplot_EPO_PL_Motifs_DEATACseq_DERNAseq.png",   width = 12.7, height = 12.8, units = "cm", res = 600, pointsize = 12)
plot(avg_log2FC.x, avg_log2FC.y, pch=19, cex=1, 
                     ylab="ATAC_Log2 Fold Change", xlab="RNA Log2 Fold Change", ylim=c(-1.7,1.7),xlim=c(-0.70,0.70),
 col=ifelse(avg_log2FC.x > 0 & avg_log2FC.y < 0,
                     "grey", ifelse(avg_log2FC.x > 0 & avg_log2FC.x > 0,
                              "red2", ifelse(avg_log2FC.x < 0 & avg_log2FC.y < 0,
 "blue", "grey")))
)
 dev.off()
 
la <- read.delim("EPO_Diff_peaks_over_TE_overlapping_with_TSS.bed", header=F)

ka$geens <- paste(ka$V6, ka$V7, ka$V8, sep = "-")
colnames(ka)

kas <- ka[,c(4,5,9,11,12,13)]

kas <- kas[!duplicated(kas[,3]),]
row.names(kas) <- kas[,3]

kas <- kas[,-c(1,2)]

kas <- kas[,-1]
head(kas)
colnames(kas) <- c("Distance","TEs","genes")

te$TEs <- paste(te$V1, te$V2, te$V3, sep = "-")

tes <- te[!duplicated(te$TEs),]

man <- merge(kas, tes, by="TEs", all=F)

row.names(man) <- man[,1]

man <- man[,c(1,2,3,8,9,10)]

manu <- merge(man, dat, by="row.names", all=F)

df <- manu[,c(6,11,3,5,10,9)]

colnames(df) <- c("TEs","Genes","Distance","ATAC_FDR","RNA_FDR","RNA_log_FC")

write.table(df, "TE_neighbor_genes_DAR_overall_merged_data.tsv", row.names=T, col.names=NA, sep="\t", dec=".", quote=F) 

df <- read.table("TE_neighbor_genes_DAR_overall_merged_data.tsv", header=T, stringsAsFactors=F, row.names=1)

library(dplyr)

library(ggplot2)

 png("Barplot_TEs_DAR_next_to_genes.png",   width = 12.7, height = 15.8, units = "cm", res = 600, pointsize = 12)
df %>%
  filter(ATAC_FDR < 0.05) %>%
  mutate(Family = sub("_.*", "", TEs)) %>%
  count(Family) %>%
  mutate(Frequency = n / sum(n)) %>%
  top_n(20, Frequency) %>%
  ggplot(aes(x = reorder(Family, Frequency), y = Frequency)) +
  geom_col(fill = "firebrick3") +
  coord_flip() +
  labs(x = "TE family", y = "Proportion among EPO-responsive DARs") +
  theme_classic(base_size = 14)
dev.off()

df$log_RNA  <- -log10(df$RNA_FDR)
df$log_ATAC <- -log10(df$ATAC_FDR)

 png("BarPlot_TE_DEATACseq_DERNAseq.png",   width = 12.7, height = 12.8, units = "cm", res = 600, pointsize = 12)
df %>%
  filter(ATAC_FDR < 0.05) %>%
  mutate(Family = sub("_.*", "", TEs)) %>%
  count(Family) %>%
  mutate(Frequency = n / sum(n)) %>%
  top_n(20, Frequency) %>%
  ggplot(aes(x = reorder(Family, Frequency), y = Frequency)) +
  geom_col(fill = "firebrick3") +
  coord_flip() +
  labs(x = "TE family", y = "Proportion among EPO-responsive DARs") +
  theme_classic(base_size = 14)
dev.off()


ggplot(df, aes(x = Distance/1000, y = RNA_log_FC, color = log_ATAC)) +
  geom_point(alpha = 0.6) +
  scale_color_gradient(low = "lightgrey", high = "darkred") +
  labs(x = "Distance from TE to gene (kb)", y = "Gene log2 FC", color = "Accessibility") +
  theme_classic(base_size = 14)
dev.off()
list.files()
df2 <- read.delim("Pairwise_ATAC_RNA_DEG_DAR_gene_names.tsv", header=T, row.names=1)
head(df2)
q()
savehistory("plotting_Supple6.Rhistory")

