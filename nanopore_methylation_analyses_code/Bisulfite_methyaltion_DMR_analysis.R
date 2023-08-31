library(DSS)
require(bsseq)
library(tidyverse)

#### bisulfite analysis (http://www.bioconductor.org/packages/devel/bioc/vignettes/DSS/inst/doc/DSS.html)
Can_3_raw <- read.table(file.path("/shared/ucl/depts/mottlab/Ziming_analysis/Bisulfite_methylation", "Can_3.multiple.deduplicated.bedGraph.gz.bismark.zero.cov"), col.names=c("chr", "start", "end", "methylation_percent", "count_methylated", "count_unmethylated"))
Can_4_raw <- read.table(file.path("/shared/ucl/depts/mottlab/Ziming_analysis/Bisulfite_methylation", "Can_4.multiple.deduplicated.bedGraph.gz.bismark.zero.cov"), col.names=c("chr", "start", "end", "methylation_percent", "count_methylated", "count_unmethylated"))
Col_1_raw <- read.table(file.path("/shared/ucl/depts/mottlab/Ziming_analysis/Bisulfite_methylation", "Col_1.multiple.deduplicated.bedGraph.gz.bismark.zero.cov"), col.names=c("chr", "start", "end", "methylation_percent", "count_methylated", "count_unmethylated"))
Col_2_raw <- read.table(file.path("/shared/ucl/depts/mottlab/Ziming_analysis/Bisulfite_methylation", "Col_2.multiple.deduplicated.bedGraph.gz.bismark.zero.cov"), col.names=c("chr", "start", "end", "methylation_percent", "count_methylated", "count_unmethylated"))

Can_3 <- Can_3_raw %>% 
  select(chr, end, count_methylated, count_unmethylated) %>%
  mutate(total_reads = count_methylated + count_unmethylated) %>%
  select(chr, end, total_reads, count_methylated) %>%
  filter(total_reads > 10) %>%
  rename("pos" = "end", "N" = "total_reads", "X" = "count_methylated")

Can_4 <- Can_4_raw %>% 
  select(chr, end, count_methylated, count_unmethylated) %>%
  mutate(total_reads = count_methylated + count_unmethylated) %>%
  select(chr, end, total_reads, count_methylated) %>%
  filter(total_reads > 10) %>%
  rename("pos" = "end", "N" = "total_reads", "X" = "count_methylated")

Col_1 <- Col_1_raw %>% 
  select(chr, end, count_methylated, count_unmethylated) %>%
  mutate(total_reads = count_methylated + count_unmethylated) %>%
  select(chr, end, total_reads, count_methylated) %>%
  filter(total_reads > 10) %>%
  rename("pos" = "end", "N" = "total_reads", "X" = "count_methylated")

Col_2 <- Col_2_raw %>% 
  select(chr, end, count_methylated, count_unmethylated) %>%
  mutate(total_reads = count_methylated + count_unmethylated) %>%
  select(chr, end, total_reads, count_methylated) %>%
  filter(total_reads > 10) %>%
  rename("pos" = "end", "N" = "total_reads", "X" = "count_methylated")


BS_bisulfite <- makeBSseqData(list(Can_3, Can_4, Col_1, Col_2), c("Can3","Can4", "Col1", "Col2") )


dmlTest.sm <- DMLtest(BS_bisulfite, group1=c("Can3", "Can4"), group2=c("Col1", "Col2"), smoothing=TRUE)
dmls_p0.001 <- callDML(dmlTest.sm, p.threshold=0.001)

# By default, the test is based on the null hypothesis that the difference in methylation levels is 0. Alternatively, users can specify a threshold for difference. For example, to detect loci with difference greater than 0.1, do:
# dmls2 <-  callDML(dmlTest, delta=0.1, p.threshold=0.001) did not run
class(dmls_p0.001) <- c("data.frame", "DMLtest") # need to assign the DMLtest to this class object, otherwise errors
dmrs_p0.01_dml0.001 <- callDMR(dmls_p0.001, p.threshold=0.01)
# DMRs are sorted by areaStat, which is defined in bsseq as the sum of the test statistics of all CpG sites within the DMR.
write.table(dmls_p0.001, file = "Bisulfite_DML_results_p0.001", row.names = FALSE, sep = "\t")
write.table(dmrs_p0.01_dml0.001, file = "Bisulfite_DMR_results_dmlp0.001_dmrp0.01", sep = "\t",row.names = FALSE)
save.image("Bisulfite_DSS_data.RData")
load("Bisulfite_DSS_data.RData") # if need to re-use it

####################### use methylKit for some other analysis #################
library(methylKit) # in linux change to another environment
file.list = list( "/shared/ucl/depts/mottlab/Ziming_analysis/Bisulfite_methylation/Can_3.multiple.deduplicated.bedGraph.gz.bismark.zero.cov",
                  "/shared/ucl/depts/mottlab/Ziming_analysis/Bisulfite_methylation/Can_4.multiple.deduplicated.bedGraph.gz.bismark.zero.cov",
                  "/shared/ucl/depts/mottlab/Ziming_analysis/Bisulfite_methylation/Col_1.multiple.deduplicated.bedGraph.gz.bismark.zero.cov",
                  "/shared/ucl/depts/mottlab/Ziming_analysis/Bisulfite_methylation/Col_2.multiple.deduplicated.bedGraph.gz.bismark.zero.cov")
CpG.bismark.cov=methRead(file.list, sample.id = list("Can3","Can4","Col1","Col2"), assembly = "Col", treatment = c(1,1,0,0), context = "CpG", pipeline = "bismarkCoverage")

# perform some basic statistics
sink("MethylKit_basic_statistics_bisulfite.log")
cat("basic statistics\n\n")
cat("# percent methylation statistics for sample Can3 \n")

pdf("MethylKit_basic_statistics.pdf")
getMethylationStats(CpG.bismark.cov[[1]],plot=TRUE,both.strands=FALSE) #for Can3
getCoverageStats(CpG.bismark.cov[[1]],plot=TRUE,both.strands=FALSE) #for Can3

# filter the data by coverage
cat("\n # filter data based on coverage \n")
filtered.CpG.bismark.cov=filterByCoverage(CpG.bismark.cov,lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)

#merge data together
cat("\n # merge all samples together \n")
meth=unite(filtered.CpG.bismark.cov, destrand=FALSE) # the unite data will only contain loci in all samples
# one could relax the stringency by the following comand
# meth.min=unite(myobj,min.per.group=1L) (not run)
head(meth)

cat("\n # perform correlation between samples \n")
getCorrelation(meth,plot=TRUE) # correlation between samples
cat("\n # cluster samples by their methylation profile \n")
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE) # cluster sample basing on methylation profiles
PCASamples(meth) # PCA plot
sink()
dev.off()

# differentially methylated region/locus test
sink("MethylKit_Differentially_methylated_locus&Regions.log")
cat("# perform differentially methylation region test \n")
myDiff=calculateDiffMeth(meth)

# get hyper methylated bases (above normal/average)
myDiff25p.hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")
# get hypo methylated bases (below normal/average)
myDiff25p.hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")
# get all differentially methylated bases
myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)
# visualize the  distribution of hypo/hyper-methylated bases/regions per chromosome 
diffMethPerChr(myDiff,plot=FALSE,qvalue.cutoff=0.01, meth.cutoff=25)
sink()

save.image("MethylKit_data.RData")














