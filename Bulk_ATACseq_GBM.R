#set your working directory and load packages
library(ChrAccR)
library(uwot)
library(BSgenome.Hsapiens.UCSC.hg19)
library(pheatmap)
library(ggplot2)
setwd("~/")
#Configure analysis: A properly curated sample annotation table is essential for the entire analysis.
AnnotFn <- file.path("GBM_ATAC_info.txt")
bamDir <- file.path("bam")
Annot <- read.table(AnnotFn, sep="\t", header=TRUE, stringsAsFactors=FALSE, fill=FALSE, strip.white=TRUE)
# add a column that ChrAccR can use to find the correct bam file for each sample
Annot[,"bamFilenameFull"] <- file.path(bamDir, Annot[,"bamFilename"])

# set a couple of analysis parameters and reporting preferences. These options can be set using the setConfigElement function 
setConfigElement("filteringCovgReqSamples", 0.25)
setConfigElement("filteringCovgCount", 4L)
setConfigElement("filteringSexChroms", TRUE)
setConfigElement("cleanMem", TRUE)
setConfigElement("doPeakCalling", TRUE)
setConfigElement("annotationColumns", c("tissueType", "donor"))
setConfigElement("annotationMinGroupSize", 1)
setConfigElement("annotationPeakGroupColumn", "tissueType")
setConfigElement("annotationPeakGroupAgreePerc", 0.25)
cleanMem(iter.gc = 6L)

theme_set(muRtools::theme_nogrid())

#Running the analysis; note: could provide specific regions of interest in regionSets
run_atac("GBM_ChrAccR", "bamFilenameFull", Annot, genome="hg19", sampleIdCol="sampleID", regionSets=NULL)
#load raw data
dsa_raw <- readRDS("~/GBM_ChrAccR/data/dsATAC_filtered/ds.rds")
#working with raw
Fetal <- getSampleAnnot(dsa_raw)[,"sampleType"] == "Fetal"
dsa_gbm <- removeSamples(dsa_raw, Fetal)
GBM11M <- getSampleAnnot(dsa_gbm)[,"sampleID"] == "GBM11M"
dsa_gbm <- removeSamples(dsa_gbm, GBM11M)

new_dsa <- addSampleAnnotCol(dsa_gbm, "lobe", lobe)
new_dsa <- addSampleAnnotCol(new_dsa, "recurrent", recurrent)
new_dsa <- addSampleAnnotCol(new_dsa, "EGFR", EGFR)
new_dsa <- addSampleAnnotCol(new_dsa, "PDGFRA", PDGFRA)
new_dsa <- addSampleAnnotCol(new_dsa, "CDKN2A_B", CDKN2A_B)
new_dsa <- addSampleAnnotCol(new_dsa, "PTEN", PTEN)
new_dsa <- addSampleAnnotCol(new_dsa, "IDH1", IDH1)
new_dsa <- addSampleAnnotCol(new_dsa, "MGMT", MGMT)
new_dsa <- addSampleAnnotCol(new_dsa, "CHR1P_LOSS", CHR1P_LOSS)
new_dsa <- addSampleAnnotCol(new_dsa, "CHR19Q_LOSS", CHR19Q_LOSS)
new_dsa <- addSampleAnnotCol(new_dsa, "CHR7P_GAIN", CHR7P_GAIN)
new_dsa <- addSampleAnnotCol(new_dsa, "CHR10Q_LOSS", CHR10Q_LOSS)

dsa_rpkm <- transformCounts(new_dsa, method="RPKM", regionTypes=c(".peaks.cons"))
dsa_l2rpkm <- transformCounts(dsa_rpkm, method="log2", regionTypes=c(".peaks.cons"))
dsa_ql2rpkm <- transformCounts(dsa_l2rpkm, method="quantile", regionTypes=c(".peaks.cons"))

counts_norm <- as.data.frame(dsa_ql2rpkm@counts$.peaks.cons)
counts_norm$variance = apply(counts_norm, 1, var)
counts_norm_ordered<- counts_norm[order(counts_norm$variance, decreasing = TRUE),]
counts_norm_ordered <- counts_norm_ordered[1:500,]
counts_norm_ordered <- counts_norm_ordered[,-31]
counts_norm_ordered_matrix <- as.matrix(counts_norm_ordered)
heatmap.2(counts_norm_ordered_matrix,
          col=viridis,
          trace="none", 
          main="500 most variable peaks",
          scale="row",
          key = TRUE)


reportDir <- file.path(".", "ChrAccR_reports_exploratory_updated3")
setConfigElement("annotationColumns", c("tissueType", "recurrent", "EGFR", "PDGFRA", "CDKN2A_B", "PTEN", "IDH1", "MGMT", "CHR1P_LOSS", "CHR19Q_LOSS", "CHR7P_GAIN", "CHR10Q_LOSS"))
createReport_exploratory(dsa_ql2rpkm, reportDir)

#try MDS plots too
dsa_qnorm <- transformCounts(dsa_gbm, method="quantile")
dsa_logqnorm <- transformCounts(dsa_qnorm, method="log10")
cm_norm <- getCounts(dsa_logqnorm, ".peaks.cons")

counts <- edgeR::DGEList(cm_norm)

#Multidimensional scaling plot because I do not have replicates for some time points, so the chraccr package won't cluseter based on age/day
#make sure row names are in same order as column names were for cm
sampleinfo <- read.delim2("GBM_ATAC_info2.txt", header = TRUE, sep = "\t")

col.tissue <- c("tumor" = "#2930e6", "margin"="#fab82a")[sampleinfo$tissueType]
col.recurrent <- c("margin" = "#2930e6", "N"="#080807", "Y"="#f7aa05")[sampleinfo$recurrent]
col.EGFR <- c("margin" = "#2930e6", "N"="#080807", "Y"="#f7aa05")[sampleinfo$EGFR]
col.PDGFRA <- c("margin" = "#2930e6", "N"="#080807", "Y"="#f7aa05")[sampleinfo$PDGFRA]
col.CDKN2A_B <- c("margin" = "#2930e6", "N"="#080807", "Y"="#f7aa05")[sampleinfo$CDKN2A_B]
col.PTEN <- c("margin" = "#2930e6", "N"="#080807", "Y"="#f7aa05")[sampleinfo$PTEN]
col.IDH1 <- c("margin" = "#2930e6", "N"="#080807", "Y"="#f7aa05")[sampleinfo$IDH1]
col.MGMT_me <- c("margin" = "#2930e6", "N"="#080807", "Y"="#f7aa05")[sampleinfo$MGMT_me]
col.1p_loss <- c("margin" = "#2930e6", "N"="#080807", "Y"="#f7aa05")[sampleinfo$X1p_loss]
col.19q_loss <- c("margin" = "#2930e6", "N"="#080807", "Y"="#f7aa05")[sampleinfo$X19q_loss]
col.7p_gain <- c("margin" = "#2930e6", "N"="#080807", "Y"="#f7aa05")[sampleinfo$X7p_gain]
col.10q_loss <- c("margin" = "#2930e6", "N"="#080807", "Y"="#f7aa05")[sampleinfo$X10q_loss]
col.lobe <- c("frontal" = "#2930e6", "temporal"="#fab82a", "occipital", "temporoparietal", "parietal", "NP")[sampleinfo$tissueType]

data.frame(sampleinfo$PDGFRA,col.PDGFRA)
pdf("GBM PCA_PDGFRA.pdf")
plotMDS(counts,col=col.PDGFRA, pch = 19)
legend("bottomright",
       fill=c("#2930e6","#080807", "#f7aa05"),legend=levels(sampleinfo$PDGFRA))
dev.off()


################################################Find DA peaks
##################Universal. Approach: make comparisons based on PCA clustering results.
##There are clear differences between margin, IDH-WT tumor, and IDH-mutant tumor, so make 3 comparisons (margin vs IDH-WT, margin vs IDH-mutant, IDH-WT vs IDH-mutant)

#make consensus peaksets
grouping <- as.vector(sampleinfo$IDH1)
grouping_updated <- grouping[-c(11:16)]

peakGrl_perSample <- readRDS("~/GBM_ChrAccR/data/peakGrl_perSample.rds")
peakGrl_perSample_1<- peakGrl_perSample[-3]
peakGrl_perSample_updated<- peakGrl_perSample_1[-32]
peakGrl_perSample_updated<- peakGrl_perSample_updated[-31]
peakGrl_perSample_updated_3<- peakGrl_perSample_updated[-c(11:16)]#take out recurrent, but keep IDHMT (updated: 4/20/22)

consPeaks <- getConsensusPeakSet(
  peakGrl_perSample_updated,
  mode = "no_by_score",
  grouping = grouping,
  groupAgreePerc = 0.75, #keep peaks, so long as 75% of the samples in a group have that peak
  groupConsSelect = TRUE, 
  scoreCol = "score",
  keepOvInfo = TRUE
)

#updated: 4/20/22
consPeaks_updated <- getConsensusPeakSet(
  peakGrl_perSample_updated_3,
  mode = "no_by_score",
  grouping = grouping_updated,
  groupAgreePerc = 0.75, #keep peaks, so long as 75% of the samples in a group have that peak
  groupConsSelect = TRUE, 
  scoreCol = "score",
  keepOvInfo = TRUE
)
#filter all peaks to only include consensus peaks
cons_peaks <- filterByGRanges(new_dsa,consPeaks, method = "white")
cons_peaks_3 <- filterByGRanges(new_dsa,consPeaks_updated, method = "white") #updated: 4/20/22
#remove recurrent samples
GBM21M <- getSampleAnnot(cons_peaks)[,"sampleID"] == "GBM21M"
cons_peaks <- removeSamples(cons_peaks, GBM21M)
GBM21T <- getSampleAnnot(cons_peaks)[,"sampleID"] == "GBM21T"
cons_peaks <- removeSamples(cons_peaks, GBM21T)
GBM22M <- getSampleAnnot(cons_peaks)[,"sampleID"] == "GBM22M"
cons_peaks <- removeSamples(cons_peaks, GBM22M)
GBM22T <- getSampleAnnot(cons_peaks)[,"sampleID"] == "GBM22T"
cons_peaks <- removeSamples(cons_peaks, GBM22T)
GBM24M <- getSampleAnnot(cons_peaks)[,"sampleID"] == "GBM24M"
cons_peaks <- removeSamples(cons_peaks, GBM24M)
GBM24T <- getSampleAnnot(cons_peaks)[,"sampleID"] == "GBM24T"
cons_peaks <- removeSamples(cons_peaks, GBM24T)
saveRDS(cons_peaks, file = "MT&WT_cons_peaks.rds")
saveRDS(cons_peaks_3, file = "updated_cons_peaks.rds")


reportDir <- file.path(".", "ChrAccR_analysis_DA_lessStrict")
setConfigElement("differentialColumns", "IDH1")
createReport_differential(cons_peaks, reportDir)
createReport_differential(cons_peaks_3, reportDir) #updated: 4/20/22

####margin Vs IDHWT
DA_mVSwt<- read.table("ChrAccR_analysis_DA_lessStrict/differential_data/diffTab_1_peaksCons.tsv", header= TRUE)
DA_mVSwt <- as.data.frame(DA_mVSwt)
DA_mVSwt <- DA_mVSwt %>% drop_na(padj) #get rid of peaks w/ NAs in padj row
sig_DA_mVSwt <- DA_mVSwt[DA_mVSwt$padj < 0.05,] #take signif diff peaks
margin_IDHWT <- sig_DA_mVSwt[sig_DA_mVSwt$log2FoldChange > 1,] #take those up in margin
IDHWT_margin <- sig_DA_mVSwt[sig_DA_mVSwt$log2FoldChange < -1,] #take those up in tumor

with(DA_mVSwt, plot(log2BaseMean, log2FoldChange, pch=20, cex=0.5 ,main="Differentially Accessible Peaks (margin vs IDH1WT)", xlim=c(2,12), ylim=c(-5,5)))
with(subset(DA_mVSwt, padj <.05 & log2FoldChange < -1), points(log2BaseMean, log2FoldChange, pch=20, cex=0.5, col="#f51836"))
with(subset(DA_mVSwt, padj <.05 & log2FoldChange > 1), points(log2BaseMean, log2FoldChange, pch=20, cex=0.5, col="#f51836"))

#####margin Vs IDHMT
DA_mVSmt<- read.table("ChrAccR_analysis_DA_lessStrict/differential_data/diffTab_2_peaksCons.tsv", header= TRUE)
DA_mVSmt <- as.data.frame(DA_mVSmt)
DA_mVSmt <- DA_mVSmt %>% drop_na(padj) #get rid of peaks w/ NAs in padj row
sig_DA_mVSmt <- DA_mVSmt[DA_mVSmt$padj < 0.05,] #take signif diff peaks
margin_IDHMT <- sig_DA_mVSmt[sig_DA_mVSmt$log2FoldChange > 1,] #take those up in margin
IDHMT_margin <- sig_DA_mVSmt[sig_DA_mVSmt$log2FoldChange < -1,] #take those up in tumor

with(DA_mVSmt, plot(log2BaseMean, log2FoldChange, pch=20, cex=0.5 ,main="Differentially Accessible Peaks (margin vs IDH1MT)", xlim=c(2,12), ylim=c(-5,5)))
with(subset(DA_mVSmt, padj <.05 & log2FoldChange < -1), points(log2BaseMean, log2FoldChange, pch=20, cex=0.5, col="#f51836"))
with(subset(DA_mVSmt, padj <.05 & log2FoldChange > 1), points(log2BaseMean, log2FoldChange, pch=20, cex=0.5, col="#f51836"))

#####IDHWT Vs IDHMT
DA_wtVSmt<- read.table("ChrAccR_analysis_DA_lessStrict/differential_data/diffTab_3_peaksCons.tsv", header= TRUE)
DA_wtVSmt <- as.data.frame(DA_wtVSmt)
DA_wtVSmt <- DA_wtVSmt %>% drop_na(padj) #get rid of peaks w/ NAs in padj row
sig_DA_wtVSmt <- DA_wtVSmt[DA_wtVSmt$padj < 0.05,] #take signif diff peaks
IDHWT_IDHMT <- sig_DA_wtVSmt[sig_DA_wtVSmt$log2FoldChange > 1,] #take those up in margin
IDHMT_IDHWT <- sig_DA_wtVSmt[sig_DA_wtVSmt$log2FoldChange < -1,] #take those up in tumor

with(DA_wtVSmt, plot(log2BaseMean, log2FoldChange, pch=20, cex=0.5 ,main="Differentially Accessible Peaks (IDH1WT vs IDH1MT)", xlim=c(2,14), ylim=c(-8,8)))
with(subset(DA_wtVSmt, padj <.05 & log2FoldChange < -1), points(log2BaseMean, log2FoldChange, pch=20, cex=0.5, col="#f51836"))
with(subset(DA_wtVSmt, padj <.05 & log2FoldChange > 1), points(log2BaseMean, log2FoldChange, pch=20, cex=0.5, col="#f51836"))

#UPDATED:4/20/22
#load the third peakCons tsv file for Wt vs Mt differences
DA_wtVSmt_2<- read.table("ChrAccR_analysis_DA_4/20updates/differential_data/diffTab_3_peaksCons.tsv", header= TRUE)
DA_mVSmwt_2<- read.table("ChrAccR_analysis_DA_4/20updates/differential_data/diffTab_1_peaksCons.tsv", header= TRUE)
DA_wtVSmt_2 <- as.data.frame(DA_wtVSmt_2)
DA_mVSwt_2 <- as.data.frame(DA_mVSwt_2)
DA_wtVSmt_2 <- DA_wtVSmt_2 %>% drop_na(padj) #get rid of peaks w/ NAs in padj row
DA_mVSwt_2 <- DA_mVSwt_2 %>% drop_na(padj) 
sig_DA_wtVSmt_2 <- DA_wtVSmt_2[DA_wtVSmt_2$padj < 0.05,] #take signif diff peaks
sig_DA_mVSwt_2 <- DA_mVSwt_2[DA_mVSwt_2$padj < 0.05,]
IDHWT_IDHMT_2 <- sig_DA_wtVSmt_2[sig_DA_wtVSmt_2$log2FoldChange > 2,] #take those up in margin
M_IDHWT_2 <- sig_DA_mVSwt_2[sig_DA_mVSwt_2$log2FoldChange > 2,]
IDHMT_IDHWT_2 <- sig_DA_wtVSmt_2[sig_DA_wtVSmt_2$log2FoldChange < -2,] #take those up in tumor
IDHWT_M_2 <- sig_DA_mVSwt_2[sig_DA_mVSwt_2$log2FoldChange < -2,]
IDHWT_IDHMT_2<- cbind(IDHWT_IDHMT_2, direction='IDHWT_gain')
IDHMT_IDHWT_2<- cbind(IDHMT_IDHWT_2, direction='IDHWT_loss')
write.csv(IDHWT_IDHMT_2, "IDHWT_gain_peaks.csv")
write.csv(IDHMT_IDHWT_2, "IDHWT_loss_peaks.csv")

##Annotate DA peaks with nearest genes
library(readxl)
library(biomaRt)
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
library(ChIPseeker)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(GenomicRanges)
#clusterProfiler is super annoying to install. You need an older version of rvcheck, so first:
devtools::install_version('rvcheck',version='0.1.8')
library(rvcheck)
#then install clusterProfiler, but make sure not to update rvcheck
BiocManager::install("clusterProfiler")
library(clusterProfiler)
##make DA regions into GRanges lists
IDHMT_IDHWT_gr <- makeGRangesFromDataFrame(IDHMT_IDHWT,
                                           keep.extra.columns=TRUE,
                                           ignore.strand=FALSE,
                                           seqinfo=NULL,
                                           seqnames.field=c("chrom"),
                                           start.field="chromStart",
                                           end.field=c("chromEnd"),
                                           strand.field="strand",
                                           starts.in.df.are.0based=FALSE)
IDHWT_IDHMT_gr <- makeGRangesFromDataFrame(IDHWT_IDHMT,
                                           keep.extra.columns=TRUE,
                                           ignore.strand=FALSE,
                                           seqinfo=NULL,
                                           seqnames.field=c("chrom"),
                                           start.field="chromStart",
                                           end.field=c("chromEnd"),
                                           strand.field="strand",
                                           starts.in.df.are.0based=FALSE)
IDHWT_margin_gr <- makeGRangesFromDataFrame(IDHWT_margin,
                                            keep.extra.columns=TRUE,
                                            ignore.strand=FALSE,
                                            seqinfo=NULL,
                                            seqnames.field=c("chrom"),
                                            start.field="chromStart",
                                            end.field=c("chromEnd"),
                                            strand.field="strand",
                                            starts.in.df.are.0based=FALSE)
margin_IDHWT_gr <- makeGRangesFromDataFrame(margin_IDHWT,
                                            keep.extra.columns=TRUE,
                                            ignore.strand=FALSE,
                                            seqinfo=NULL,
                                            seqnames.field=c("chrom"),
                                            start.field="chromStart",
                                            end.field=c("chromEnd"),
                                            strand.field="strand",
                                            starts.in.df.are.0based=FALSE)
IDHMT_margin_gr <- makeGRangesFromDataFrame(IDHMT_margin,
                                            keep.extra.columns=TRUE,
                                            ignore.strand=FALSE,
                                            seqinfo=NULL,
                                            seqnames.field=c("chrom"),
                                            start.field="chromStart",
                                            end.field=c("chromEnd"),
                                            strand.field="strand",
                                            starts.in.df.are.0based=FALSE)
margin_IDHMT_gr <- makeGRangesFromDataFrame(margin_IDHMT,
                                            keep.extra.columns=TRUE,
                                            ignore.strand=FALSE,
                                            seqinfo=NULL,
                                            seqnames.field=c("chrom"),
                                            start.field="chromStart",
                                            end.field=c("chromEnd"),
                                            strand.field="strand",
                                            starts.in.df.are.0based=FALSE)

#find TSS sites near peaks
peakAnno_IDHMT_IDHWT = annotatePeak(IDHMT_IDHWT_gr, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoDF_IDHMT_IDHWT<- as.data.frame(peakAnno_IDHMT_IDHWT@anno)
#only include if +/-1000 bp
peakAnno_promoter_IDHMT_IDHWT <- peakAnnoDF_IDHMT_IDHWT[peakAnnoDF_IDHMT_IDHWT$distanceToTSS < 2000,]
peakAnno_promoter_IDHMT_IDHWT <- peakAnno_promoter_IDHMT_IDHWT[peakAnno_promoter_IDHMT_IDHWT$distanceToTSS > -1000,]
write.csv(as.data.frame(peakAnno_promoter_IDHMT_IDHWT), "GO_promoter_IDHMT_IDHWT.csv")

peakAnno_IDHWT_IDHMT = annotatePeak(IDHWT_IDHMT_gr, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoDF_IDHWT_IDHMT<- as.data.frame(peakAnno_IDHWT_IDHMT@anno)
peakAnno_promoter_IDHWT_IDHMT <- peakAnnoDF_IDHWT_IDHMT[peakAnnoDF_IDHWT_IDHMT$distanceToTSS < 2000,]
peakAnno_promoter_IDHWT_IDHMT <- peakAnno_promoter_IDHWT_IDHMT[peakAnno_promoter_IDHWT_IDHMT$distanceToTSS > -1000,]
write.csv(as.data.frame(peakAnno_promoter_IDHWT_IDHMT), "GO_promoter_IDHWT_IDHMT.csv")

peakAnno_IDHWT_margin = annotatePeak(IDHWT_margin_gr, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoDF_IDHWT_margin<- as.data.frame(peakAnno_IDHWT_margin@anno)
peakAnno_promoter_IDHWT_margin <- peakAnnoDF_IDHWT_margin[peakAnnoDF_IDHWT_margin$distanceToTSS < 2000,]
peakAnno_promoter_IDHWT_margin <- peakAnno_promoter_IDHWT_margin[peakAnno_promoter_IDHWT_margin$distanceToTSS > -1000,]
write.csv(as.data.frame(peakAnno_promoter_IDHWT_margin), "GO_promoter_IDHWT_margin.csv")

peakAnno_margin_IDHWT = annotatePeak(margin_IDHWT_gr, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoDF_margin_IDHWT<- as.data.frame(peakAnno_margin_IDHWT@anno)
peakAnno_promoter_margin_IDHWT <- peakAnnoDF_margin_IDHWT[peakAnnoDF_margin_IDHWT$distanceToTSS < 2000,]
peakAnno_promoter_margin_IDHWT <- peakAnno_promoter_margin_IDHWT[peakAnno_promoter_margin_IDHWT$distanceToTSS > -1000,]
write.csv(as.data.frame(peakAnno_promoter_margin_IDHWT), "GO_promoter_margin_IDHWT.csv")

peakAnno_IDHMT_margin = annotatePeak(IDHMT_margin_gr, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoDF_IDHMT_margin<- as.data.frame(peakAnno_IDHMT_margin@anno)
peakAnno_promoter_IDHMT_margin <- peakAnnoDF_IDHMT_margin[peakAnnoDF_IDHMT_margin$distanceToTSS < 2000,]
peakAnno_promoter_IDHMT_margin <- peakAnno_promoter_IDHMT_margin[peakAnno_promoter_IDHMT_margin$distanceToTSS > -1000,]
write.csv(as.data.frame(peakAnno_promoter_IDHMT_margin), "GO_promoter_IDHMT_margin.csv")

peakAnno_margin_IDHMT = annotatePeak(margin_IDHMT_gr, tssRegion=c(-1000, 1000), TxDb=txdb, annoDb="org.Hs.eg.db")
peakAnnoDF_margin_IDHMT<- as.data.frame(peakAnno_margin_IDHMT@anno)
peakAnno_promoter_margin_IDHMT <- peakAnnoDF_margin_IDHMT[peakAnnoDF_margin_IDHMT$distanceToTSS < 2000,]
peakAnno_promoter_margin_IDHMT <- peakAnno_promoter_margin_IDHMT[peakAnno_promoter_margin_IDHMT$distanceToTSS > -1000,]
write.csv(as.data.frame(peakAnno_promoter_margin_IDHMT), "GO_promoter_margin_IDHMT.csv")

##GO analysis; *note that I have not subsetted to only use those peaks that are in the promoter regions
ego_margin_IDHMT <- enrichGO(as.data.frame(peakAnno_promoter_margin_IDHMT)$geneId, 
                             OrgDb = org.Hs.eg.db,             
                             ont = "BP", 
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05, 
                             readable = TRUE)
dotplot(ego_margin_IDHMT, showCategory=10)
write.csv(as.data.frame(ego_margin_IDHMT), "GO_margin_IDHMT.csv")

ego_IDHMT_margin <- enrichGO(as.data.frame(peakAnno_promoter_IDHMT_margin)$geneId, 
                             OrgDb = org.Hs.eg.db,             
                             ont = "BP", 
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05, 
                             readable = TRUE)
dotplot(ego_IDHMT_margin, showCategory=10)
write.csv(as.data.frame(ego_IDHMT_margin), "GO_IDHMT_margin.csv")

ego_margin_IDHWT <- enrichGO(as.data.frame(peakAnno_promoter_margin_IDHWT)$geneId, 
                             OrgDb = org.Hs.eg.db,             
                             ont = "BP", 
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05, 
                             readable = TRUE)
dotplot(ego_margin_IDHWT, showCategory=10)
write.csv(as.data.frame(ego_margin_IDHWT), "GO_margin_IDHWT.csv")


ego_IDHWT_margin <- enrichGO(as.data.frame(peakAnno_promoter_IDHWT_margin)$geneId, 
                             OrgDb = org.Hs.eg.db,             
                             ont = "BP", 
                             pAdjustMethod = "BH", 
                             qvalueCutoff = 0.05, 
                             readable = TRUE)
dotplot(ego_IDHWT_margin, showCategory=10)
write.csv(as.data.frame(ego_IDHWT_margin), "GO_IDHWT_margin.csv")

ego_IDHWT_IDHMT <- enrichGO(as.data.frame(peakAnno_IDHWT_IDHMT)$geneId, 
                            OrgDb = org.Hs.eg.db,             
                            ont = "BP", 
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05, 
                            readable = TRUE)
dotplot(ego_IDHWT_IDHMT, showCategory=10)
write.csv(as.data.frame(ego_IDHWT_IDHMT), "GO_IDHWT_IDHMT.csv")

ego_IDHMT_IDHWT <- enrichGO(as.data.frame(peakAnno_IDHMT_IDHWT)$geneId, 
                            OrgDb = org.Hs.eg.db,             
                            ont = "BP", 
                            pAdjustMethod = "BH", 
                            qvalueCutoff = 0.05, 
                            readable = TRUE)
dotplot(ego_IDHMT_IDHWT, showCategory=10)
write.csv(as.data.frame(ego_IDHMT_IDHWT), "GO_IDHMT_IDHWT.csv")

####################################################
#generate DA heatmap with stemness scores
#read in counts and coords
ATAC_GBM_peaks<- read.csv("C:/Users/caitl/Documents/third rotation/GBM/ATAC/CPM_peaks_coords.csv", header = T)
ATAC_GBM_peak_coords<- read.csv("C:/Users/caitl/Documents/third rotation/GBM/ATAC/new_coords.csv", header = T)
rownames(ATAC_GBM_peak_coords) <- ATAC_GBM_peak_coords$X
rownames(ATAC_GBM_peaks) <- ATAC_GBM_peak_coords$X

ATAC_GBM_peak_coords <- ATAC_GBM_peak_coords %>% drop_na(padj) #get rid of peaks w/ NAs in padj row
sig_da <- ATAC_GBM_peak_coords[ATAC_GBM_peak_coords$padj < 0.05,]
margin_IDHWT <- sig_da[sig_da$log2FoldChange > 2,] #take those up in margin 2282 peaks
IDHWT_margin <- sig_da[sig_da$log2FoldChange < -2,]
sig_da <- sig_da %>% filter(log2FoldChange > 2 | log2FoldChange < -2)

ATAC_GBM_peaks <- ATAC_GBM_peaks[rownames(ATAC_GBM_peaks) %in% rownames(sig_da),c(11,4,17,19,7,12,13,6,3,15,10,18,14,5,16,2)]

ATAC_matrix <- as.matrix(ATAC_GBM_peaks)


#make heatmap and grab row order
da_heat<- pheatmap::pheatmap(ATAC_matrix, scale = "row", cluster_cols = FALSE)
row.order <- da_heat$tree_row$order

library(ComplexHeatmap)
library(RColorBrewer)
#just ATAC heatmap

ATAC_matrix<- ATAC_matrix[row.order,]

#####calculate "stemness" score and add column label
setwd("~/third rotation/GBM/RNA")
#load data
feature_counts <- read_xlsx("GBM feature counts.xlsx")
#removing duplicates genes 
df1<- feature_counts[!duplicated(feature_counts$Geneid), ]
#remove rows for which there is no approved gene symbol (df2), make gene name the row name (df2), remove the gene name column (df3).
df2 <- df1[!is.na(df1$Geneid), ]
df2 <- as.data.frame(df2)
row.names(df2) <- df2$Geneid
df3 <- df2[,c(11,7,21,22,9,15,17,8,6,19,10,38,35,26,37,24)]
#read in dataframe with stem genes
stem_genes <- read_xlsx("stem genes.xlsx", sheet = 2)
#Pull feturecounts data only for these genes
df3_l2cpm<- cpm(df3, log = F)
df4 <- df3_l2cpm[row.names(df3_l2cpm) %in% stem_genes$`Stem genes`,]
stem_GBM_matrix <- as.matrix(df4)
zscore_ColSum<- as.data.frame(colSums(stem_GBM_matrix))
scaled_stem_GBM<- as.data.frame(scale(zscore_ColSum))
library(circlize)
col_fun = colorRamp2(c(-4, 0, 4), c( "red", "white", "blue"))
column_ha = HeatmapAnnotation(foo1 = scaled_stem_GBM$`colSums(stem_GBM_matrix)`, col = list(foo = col_fun))

#ATAC heatmap with stemness column label and astro gene labels
#flipped the order of cols so tumor is on left and margin on right
ComplexHeatmap::pheatmap(ATAC_matrix, name = "mat",cluster_rows = FALSE, cluster_cols = FALSE,
                         color =colorRampPalette(brewer.pal(n=9,name="Blues"), bias = 0.6)(100), scale = "row", top_annotation = column_ha)
