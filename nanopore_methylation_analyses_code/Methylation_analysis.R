library(NanoMethViz)
library(Rsamtools)
library(dplyr)

setwd("/shared/ucl/depts/mottlab/Ziming_analysis/Nanopore_methylation")
directory <- "/shared/ucl/depts/mottlab/Ziming_analysis/Nanopore_methylation/"

####### visualization with NanoMethViz #######################
# import methylation data
directory <- "/shared/ucl/depts/mottlab/Ziming_analysis/Nanopore_methylation/"
methy_tabix <- file.path(directory, "methy_data_new.bgz")
samples <- c("Can", "Col")
methy_calls <- c("Can_colref_per_read_modified_base_calls.txt.gz", "Col_per_read_modified_base_calls.txt.gz")
create_tabix_file(methy_calls, methy_tabix, samples)
read.table(gzfile("methy_data_new.bgz"), col.names = methy_col_names(), nrows = 6)

group <- c("Can", "Col")
samples <- c("Can", "Col")
sample_anno <- data.frame(samples, group, stringsAsFactors = FALSE)
sample_anno
sample_anno <- sample_anno %>%
  as.data.frame() %>%
  dplyr::rename(
    sample = samples
  )

# import annotation
extraCols_name <- c(strandy="character", gene_idy="character", transcript_idy="character", symbol="character")
anno_exon <- rtracklayer::import("Col_all_anno.bed", format = "BED", extraCols = extraCols_name)
head(anno_exon)

anno <- anno_exon %>%
  as.data.frame() %>%
  dplyr::select(gene_idy, seqnames, strandy, start, end, transcript_idy, symbol) %>%
  dplyr::rename(
    chr = seqnames,
    gene_id = gene_idy,
    transcript_id = transcript_idy,
    strand = strandy,
    symbol = symbol
  ) %>%
  dplyr::select(gene_id, chr, strand, start, end, transcript_id, symbol)

head(anno)

nmeth_results <- NanoMethResult("methy_data_new.bgz", sample_anno, anno)
nmeth_results

# to plot methylation heat map
pdf(file="rRNA_region.pdf")
plot_grange_heatmap(nmeth_results, GenomicRanges::GRanges("Chr2_RagTag_polished:0-22472"))+ggtitle("rRNA_Chr2")
plot_grange_heatmap(nmeth_results, GenomicRanges::GRanges("Chr4_RagTag_polished:0-22472"))+ggtitle("rRNA_Chr4")
dev.off()
ggsave("rRNA_Chr4.png", width = 5, height = 5)


################ perform DML Aand DMR analysis #####################
library(DSS)
require(bsseq)
library(tidyverse)

#### Nanopore methylation analysis
bss <- methy_to_bsseq(nmeth_results)
bss
dmlTesm_nanopore.sm <- DMLtest(bss, group1="Can", group2="Col", smoothing=TRUE)
dmls_nanopore <- callDML(dmlTesm_nanopore.sm, p.threshold=0.01)
class(dmls_nanopore) <- c("data.frame","DMLtest")
dmrs <- callDMR(dmls_nanopore, p.threshold=0.05)
write.csv(dmrs,"DMRs_Can_Col_p0.05.csv")
write.csv(dmls_nanopore,"DMLs_Can_Col_p0.01.csv")
save.image("Nanopore_NanoMethViz_DSS_data.RData")

lmr <- bsseq_to_log_methy_ratio(bss)

# convert nanometh result to edgeR format
edger_mat <- methy_to_edger(nmeth_results)
head (edger_mat)

###################### visualization of DMR for Nanopore methylation data
# loading saved results from previous bsseq analysis
bsseq_dmr <- read.table(
  system.file(package = "NanoMethViz", "dmr_subset.tsv.gz"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)
plot_gene(nmeth_results, anno_regions = bsseq_dmr, heatmap = TRUE)

pdf(file="rRNA_Chr2_several.pdf")
plot_gene(nmeth_results, "18s_Chr2_5", heatmap = TRUE) + ggtitle("18s_Chr2_5")
plot_gene(nmeth_results, "Col_AT2G01020_1_rRNA_5S", heatmap = TRUE) + ggtitle("Col_AT2G01020_1_rRNA_5.8S")
plot_gene(nmeth_results, "25s_Chr2_9", heatmap = TRUE) + ggtitle("25s_Chr2_9")
plot_gene(nmeth_results, "Col_AT2G03855_1_miRNA_ath-MIR5642b", heatmap = TRUE) + ggtitle("Col_AT2G03855_1_miRNA_ath-MIR5642b")
plot_gene(nmeth_results, "18s_Chr2_6 ", heatmap = TRUE) + ggtitle("18s_Chr2_6 ")
dev.off()

# to compare with the NOR2 and NOR4, we also look at the rRNA array in Chr1 (should be the 5s)
pdf(file="rRNA_Chr1_array.pdf")
plot_grange_heatmap(nmeth_results, GenomicRanges::GRanges("Chr1_RagTag_polished:2410550-4010729"))+ggtitle("rRNA_Chr1")
dev.off()

save.image(file="methylation.RData") 
load(file="methylation.RData")


