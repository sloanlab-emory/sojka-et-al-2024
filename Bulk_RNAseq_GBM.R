#######################compare IDHwt to IDHmt
#remove organoid, fetal/adult, recurrent, and all margin

##try with DESEQ2
library(DESeq2)
library(ggplot2)

countData <- read.csv('C:/Users/caitl/Documents/third rotation/GBM/RNA/tumor_counts.csv', header = TRUE, sep = ",")
metaData <- read.csv('C:/Users/caitl/Documents/third rotation/GBM/RNA/tumor_meta.csv', header = TRUE, sep = ",")
#data with recurrent (not used for the IDH comparison)
countData<- read.csv("C:/Users/caitl/Documents/third rotation/GBM/RNA/tumor_counts2.csv", header = TRUE, sep = ",")
metaData<- read.csv("C:/Users/caitl/Documents/third rotation/GBM/RNA/tumor_meta2.csv", header = TRUE, sep = ",")
countData <- countData[,-c(12,14,16,19)]
metaData <- metaData[-c(11,13,15,18),]
#data with all tumor and margin
countData <- read.csv("C:/Users/caitl/Documents/third rotation/GBM/RNA/GBM_counts2.csv", header = T, sep = ",")
dds <- DESeqDataSetFromMatrix(countData=df1, 
                              colData=GBM_meta, 
                              design=~Group, tidy = TRUE)


df1<- countData[!duplicated(countData$Gene), ]
dds <- DESeqDataSetFromMatrix(countData=df1, 
                              colData=metaData, 
                              design=~IDH, tidy = TRUE)
df2 <- as.data.frame(df1)
row.names(df2) <- df2$Gene
df3 <- df2[,-1]

#Filtering count data:
#Step1: filter out the genes for which there are no reads, or there are inconsistent reads across replicate samples, or there are low reads. Chose a normalization technique: CPM (counts per million), for example (GBM_CPM). We filter using CPM values rather than counts because they account for differences in sequencing depth between samples. 
GBM_CPM <- edgeR::cpm(df3)
head(GBM_CPM)

#Step2: Our filtering rule is to keep transcripts that have CPM > X in at least two samples. X is the cutoff value (in CPM) that will consistently yield a minimum transcript count. We need to first determine the scientifically relevant minimum transcript count. As a general rule, a good threshold can be chosen for a CPM value that corresponds to a count of 10 (for 10-20M reads/sample depth). This will vary based on the depth of sequence; for low depth sequence <1M, try a count of 1. CPM = counts per million, or how many counts would I get for a gene if the sample had a library size of 1M. So for a library size of 1M, 1 count = 1 CPM. For a library size of 10M, 10 counts = 1 CPM. For a library size of 20M, 10 counts = 0.5 CPM.

#Step 3: Next, impose the threshold (0.5 is an example value, replace with determined X value).
GBM_thresh <- GBM_CPM > 0.5

#Step 4: Identify/subset the genes for which CPM > 0.5 in at least two samples. Use DGElist to convert conunts.keep df back to an oject (GBM_count) for ease of analysis moving forward.
GBM_keep <- rowSums(GBM_thresh) >= 2
GBM_counts.keep <- df3[GBM_keep,]
GBM_counts.keep<- cpm(GBM_counts.keep, log = T)
pcaCoord <- muRtools::getDimRedCoords.pca(t(GBM_counts.keep), components = c(1, 2))
metaData2 <- metaData[,c(1,7)]
rownames(metaData2) <- metaData2$Sample
metaData2 <- as.data.frame(metaData2[,-1])
metaData2 <- as.data.frame(metaData2)
muRtools::getDimRedPlot(pcaCoord, annot=metaData2, colorCol="category", shapeCol = NULL, addLabels = T)

dds <- DESeq(dds)
res <- results(dds)
head(results(dds, tidy=TRUE))
summary(res) 
res <- res[order(res$padj),]
head(res)
IDHmtvswt_astro<- as.data.frame(res)

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="IDH")
plotPCA(vsdata, intgroup="category")

plotCounts(dds, gene="PTTG1", intgroup="IDH")
plotCounts(dds, gene="TOP2A", intgroup="IDH")
plotCounts(dds, gene="ETNPPL", intgroup="IDH")
plotCounts(dds, gene="ALDOC", intgroup="IDH")


#look at organoid modules within DESEQ results
IDHmtvswt_astro$SYMBOL <- rownames(IDHmtvswt_astro)
signifgenes_IDHmtVsIDHwt <- IDHmtvswt_astro  %>% filter(padj<.05 & abs(log2FoldChange)>2)
write.csv(signifgenes_IDHmtVsIDHwt, file = "signifgenes_IDH1mtVsIDHwtTumor.csv")
signifgenes <- IDHmtvswt_astro %>% filter(padj < 0.05) #1550 genes
signifgenes_down_astro <- signifgenes  %>% filter(log2FoldChange < -2) #890
signifgenes_up_astro <- signifgenes  %>% filter(log2FoldChange > 2) #211; 
write.csv(signifgenes_down_astro, file = "signifgenes_down_astro.csv")
write.csv(signifgenes_up_astro, file = "signifgenes_up_astro.csv")
#save IDH DESEQ R objects
save(countData, metaData, df1, dds, res, IDHmtvswt, vsdata, signifgenes, signifgenes_down, signifgenes_up, file="IDH_DESeq.RData")


down_ann <- select(org.Hs.eg.db,keys=rownames(signifgenes_down),keytype = "SYMBOL" ,columns=c("ENTREZID","SYMBOL","GENENAME"))
head(down_ann)
all_ann <- select(org.Hs.eg.db,keys=rownames(IDHmtvswt_astro),keytype = "SYMBOL" ,columns=c("ENTREZID","SYMBOL","GENENAME"))
allgenes <- as.character(all_ann$ENTREZID)
down_GO <- as.character(down_ann$ENTREZID)
ego_down <- enrichGO(gene          = down_GO,
                     universe      = allgenes,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     readable      = TRUE)
dotplot(ego_down, showCategory=20, orderBy = "count")
up_ann <- select(org.Hs.eg.db,keys=rownames(signifgenes_up),keytype = "SYMBOL" ,columns=c("ENTREZID","SYMBOL","GENENAME"))
up_GO <- as.character(up_ann$ENTREZID)
ego_up <- enrichGO(gene          = up_GO,
                   universe      = allgenes,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   readable      = TRUE)
dotplot(ego_up, showCategory=20, orderBy = "count")




#genes near diff 5hmc that appear to be DE  in bulk tissue RNA-seq
library(readxl)
IDHup_5hmc<- read_xlsx("C:/Users/caitl/Documents/third rotation/5HMC and 5MC collab/Copy of PositiveCorrelate_Genebody5hmC_GeneExpression.xlsx", sheet = 1, col_names = T)
as.vector(intersect(signifgenes_up$SYMBOL, IDHup_5hmc$gene))#overlap of Yangping's list with my list
#"SMAD9"    "NEGR1"    "RTN1"     "PDE8B"    "FAM171A1" "CSMD3"    "DLGAP1"  
# "SSTR2"    "TOX" 
as.vector(intersect(blu_mod$gene_id, IDHup_5hmc$gene))#overlap with Middle/Late module
#"NEGR1"    "SORCS3"   "SOX6"     "LTBP3"    "PRKCA"    "ZHX3"     "PHACTR3"  "ZBTB20"   "EPHB1"    "PDE8B"   
#"MEF2C"    "TRAF3IP2" "RORB"     "GRIA3"
as.vector(intersect(t_mod$gene_id, IDHup_5hmc$gene))#overlap with Early module
#"SYDE2"   "RHOU"    "RHOBTB1" "SLIT1"   "AFAP1L2" "PLEKHA5" "ANKS1B"  "SRRM4"   "PITPNM3" "SSTR2"   "RBFOX3" 
#"PXDN"    "AFF3"    "CNTN4"   "MSX1"    "CDH18"   "SLIT3"   "KHDRBS2" "BNC2"    "SLC24A2" "FOCAD"   "DNAI1" 
as.vector(intersect(br_mod$gene_id, IDHup_5hmc$gene))#overlap with Middle module
#"RTN1"  "HAGHL" "TIMP3" "NRG1"  "NTNG2"
as.vector(intersect(g_mod$gene_id, IDHup_5hmc$gene))#overlap with Late module
#"FAM171A1" "TMEM38A"  "TNIK"     "ADRA1A"   "AR"

#genes near diff 5hmc
library(readxl)
IDHup_5hmc<- read_xlsx("C:/Users/caitl/Documents/third rotation/5HMC and 5MC collab/Copy of Genelist.xlsx", sheet = 1, col_names = T)
as.vector(intersect(signifgenes_up$SYMBOL, IDHup_5hmc$gene))#overlap of Yangping's list with my list

#genes that overlap with fetal/adult astro genes
Fastro_genes <- na.omit(cell_type$F.astro)
Aastro_genes <- na.omit(cell_type$A.astro)
as.vector(intersect(Fastro_genes, IDHup_5hmc$gene))
#"NRG1"  "PXDN"  "LTBP3" "DNAI1" "DGKK" 
as.vector(intersect(Aastro_genes, IDHup_5hmc$gene))
#"TIMP3"    "TNIK"     "PRKCA"    "RORB"     "RFX4"     "FAM171A1"
#"FOXO1"    "NRG3"     "ZHX3"     "C12orf49" "EPHB1"    "YAP1"    
#"TRAF3IP2" "TMEM38A" 



#genes that are down and in various modules
down_blue <- as.vector(intersect(signifgenes_down_astro$SYMBOL, blu_mod$gene_id)) #54
down_turq <- as.vector(intersect(signifgenes_down_astro$SYMBOL, t_mod$gene_id)) #145
down_brown <- as.vector(intersect(signifgenes_down_astro$SYMBOL, br_mod$gene_id)) #54
down_green <- as.vector(intersect(signifgenes_down_astro$SYMBOL, g_mod$gene_id)) #32
down_yellow <- as.vector(intersect(signifgenes_down_astro$SYMBOL, y_mod$gene_id)) #125
up_blue <- as.vector(intersect(signifgenes_up_astro$SYMBOL, blu_mod$gene_id)) #46
up_turq <- as.vector(intersect(signifgenes_up_astro$SYMBOL, t_mod$gene_id)) #18
up_brown <- as.vector(intersect(signifgenes_up_astro$SYMBOL, br_mod$gene_id)) #7
up_green <- as.vector(intersect(signifgenes_up_astro$SYMBOL, g_mod$gene_id)) #13
up_yellow <- as.vector(intersect(signifgenes_up_astro$SYMBOL, y_mod$gene_id)) #2
save(down_blue, down_brown, down_green, down_turq, down_yellow, up_blue, up_brown, up_green, up_turq, up_yellow, file = "IDH_DEG_moduleOverlap.RData")
write.csv(down_blue, file = "blue_down.csv")
write.csv(down_turq, file = "turq_down.csv")
write.csv(down_brown, file = "brown_down.csv")
write.csv(down_green, file = "green_down.csv")
write.csv(down_yellow, file = "yellow_down.csv")

write.csv(up_blue, file = "blue_up.csv")
write.csv(up_turq, file = "turq_up.csv")
write.csv(up_brown, file = "brown_up.csv")
write.csv(up_green, file = "green_up.csv")
write.csv(up_yellow, file = "yellow_up.csv")

up_Fastro <- as.vector(intersect(signifgenes_up_astro$SYMBOL, Fastro_genes)) #only 3 genes
up_Aastro <- as.vector(intersect(signifgenes_up_astro$SYMBOL, Aastro_genes)) #43 genes
down_Fastro <- as.vector(intersect(signifgenes_down_astro$SYMBOL, Fastro_genes)) #162 genes
down_Aastro <- as.vector(intersect(signifgenes_down_astro$SYMBOL, Aastro_genes)) #28 genes
write.csv(up_Fastro, file = "up_Fastro.csv")
write.csv(up_Aastro, file = "up_Aastro.csv")
write.csv(down_Fastro, file = "down_Fastro.csv")
write.csv(down_Aastro, file = "down_Aastro.csv")
write.csv(Fastro_genes, file = "Fastro_genes.csv")
write.csv(Aastro_genes, file = "Aastro_genes.csv")
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes IDH1mt vs IDHwt Fetal/Adult genes",xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% down_Fastro | rownames(res) %in% up_Fastro), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, rownames(res) %in% down_Aastro | rownames(res) %in% up_Aastro), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))

down_all <- as.vector(intersect(signifgenes_down$SYMBOL, submod$gene_id))
up_all <- as.vector(intersect(signifgenes_up$SYMBOL, submod$gene_id))
with(lrt.res$table, plot(logFC, -log10(PValue), pch=20, main="DE genes Tumor vs Margin Module overlap", xlim=c(-12,12), ylim=c(0,65)))
with(subset(lrt.res$table, FDR<.05 & logFC < -2), points(logFC, -log10(PValue), pch=20, col="#C1C7C9"))
with(subset(lrt.res$table, FDR<.05 & logFC>2), points(logFC, -log10(PValue), pch=20, col="#C1C7C9"))
with(subset(lrt.res$table, SYMBOL %in% down_all | SYMBOL %in% up_all), points(logFC, -log10(PValue), pch=20, col="#549b96"))
#lots of difference between tumor and margin, but roughly 1/3 of DE genes fall into one of my astro maturation modules ()
#math for above: ((down_all + up_all)/(signifgenes_down + signifgenes_up)) 

#let's see if this is specific to this set of genes
website_data<- read_xlsx("C:/Users/caitl/Documents/third rotation/Copy of website_data_4_2019 (version 1).xlsx", sheet = "FC>1", col_names = TRUE)
website_sample<- website_data[sample(nrow(website_data), 3944, replace = FALSE, prob = NULL),]
down_website <- as.vector(intersect(signifgenes_down$SYMBOL, website_sample$genes)) #77, 67, 65
up_website <- as.vector(intersect(signifgenes_up$SYMBOL, website_sample$genes)) #28, 25, 21 
#105 total out of (695+403)=1098, so roughly 9.5%, 8%, 8%

website_data<- read_xlsx("C:/Users/caitl/Documents/third rotation/Copy of website_data_4_2019 (version 1).xlsx", sheet = "Sheet6", col_names = TRUE)
fetal_website <- website_data[website_data$fetal>2,1]
adult_website <- website_data[website_data$adult>2,1]
website_both <- unique(c(fetal_website$...1, adult_website$...1))
website_both2<- sample(website_both, 3944, replace = FALSE, prob = NULL)
down_website <- as.vector(intersect(signifgenes_down$SYMBOL, website_both2)) #107, 98, 104
up_website <- as.vector(intersect(signifgenes_up$SYMBOL, website_both2)) #35, 31, 35
#13%, 11.7%, 12.7%



#Make volcano plots for following modules: blue, turquoise, yellow, brown, green

#blue
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes IDHmt vs wt Middle/Late Module", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% down_blue | rownames(res) %in% up_blue), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))

#turquoise
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes IDHmt vs wt Early Module", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% down_turq | rownames(res) %in% up_turq), points(log2FoldChange, -log10(pvalue), pch=20, col="#30D5C8"))


#Brown
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes IDHmt vs wt Middle Module", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% down_brown | rownames(res) %in% up_brown), points(log2FoldChange, -log10(pvalue), pch=20, col="brown"))


#Green
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes IDHmt vs wt Late Module", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% down_green | rownames(res) %in% up_green), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))


#Yellow
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes IDHmt vs wt Early/Middle Module", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% down_yellow | rownames(res) %in% up_yellow), points(log2FoldChange, -log10(pvalue), pch=20, col="#FFD700"))


#look for overlap with 5hmC genes
IDHmut_up_5hmc_GB<- read_xlsx("C:/Users/caitl/Documents/third rotation/5HMC and 5MC collab/IDHmut_vs_IDHwt_Dif5mC_Dif5hmC_GeneList_updated.xlsx", col_names = TRUE)
IDHmut_up_5hmc_GB <- as.vector(IDHmut_up_5hmc_GB$ABCA1)
IDHmut_down_5hmc_GB<- read_xlsx("C:/Users/caitl/Documents/third rotation/5HMC and 5MC collab/IDHmut_vs_IDHwt_Dif5mC_Dif5hmC_GeneList_updated.xlsx", col_names = TRUE, sheet = 2)
IDHmut_down_5hmc_GB <- as.vector(IDHmut_down_5hmc_GB$CREB5)
IDHmut_up_5hmc_TSS<- read_xlsx("C:/Users/caitl/Documents/third rotation/5HMC and 5MC collab/IDHmut_vs_IDHwt_Dif5mC_Dif5hmC_GeneList_updated.xlsx", col_names = TRUE, sheet = 3)
IDHmut_up_5hmc_TSS <- as.vector(IDHmut_up_5hmc_TSS$ABLIM1)

IDHmut_up_5mc_genebody<- read_xlsx("C:/Users/caitl/Documents/third rotation/5HMC and 5MC collab/IDHmut_vs_IDHwt_Dif5mC_Dif5hmC_GeneList_updated.xlsx", col_names = TRUE, sheet = 5)
IDHmut_up_5mc_genebody <- as.vector(IDHmut_up_5mc_genebody$AATK)
IDHmut_down_5mc_genebody<- read_xlsx("C:/Users/caitl/Documents/third rotation/5HMC and 5MC collab/IDHmut_vs_IDHwt_Dif5mC_Dif5hmC_GeneList_updated.xlsx", col_names = TRUE, sheet = 6)
IDHmut_down_5mc_genebody <- as.vector(IDHmut_down_5mc_genebody$ARHGAP39)
IDHmut_up_5mc_TSSpromoter<- read_xlsx("C:/Users/caitl/Documents/third rotation/5HMC and 5MC collab/IDHmut_vs_IDHwt_Dif5mC_Dif5hmC_GeneList_updated.xlsx", col_names = TRUE, sheet = 7)
IDHmut_up_5mc_TSSpromoter <- as.vector(IDHmut_up_5mc_TSSpromoter$ABHD12B)
IDHmut_down_5mc_TSSpromoter<- read_xlsx("C:/Users/caitl/Documents/third rotation/5HMC and 5MC collab/IDHmut_vs_IDHwt_Dif5mC_Dif5hmC_GeneList_updated.xlsx", col_names = TRUE, sheet = 8)
IDHmut_down_5mc_TSSpromoter <- as.vector(IDHmut_down_5mc_TSSpromoter$EGFR)

#try FC of 1 as well
signifgenes_down <- signifgenes %>% filter(log2FoldChange < -1) 
signifgenes_up <- signifgenes %>% filter(log2FoldChange > 1) 

down_IDHmut_down_5hmc_GB <- as.vector(intersect(signifgenes_down$SYMBOL, IDHmut_down_5hmc_GB))
down_IDHmut_down_5mc_GB <- as.vector(intersect(signifgenes_down$SYMBOL, IDHmut_down_5mc_genebody))
down_IDHmut_down_5mc_TSS <- as.vector(intersect(signifgenes_down$SYMBOL, IDHmut_down_5mc_TSSpromoter))

up_IDHmut_up_5hmc_GB <- as.vector(intersect(signifgenes_up$SYMBOL, IDHmut_up_5hmc_GB))
up_IDHmut_up_5hmc_TSS <- as.vector(intersect(signifgenes_up$SYMBOL, IDHmut_up_5hmc_TSS))
up_IDHmut_up_5mc_GB <- as.vector(intersect(signifgenes_up$SYMBOL, IDHmut_up_5mc_genebody))
up_IDHmut_up_5mc_TSS <- as.vector(intersect(signifgenes_up$SYMBOL, IDHmut_up_5mc_TSSpromoter))

down_IDHmut_up_5mc_GB <- as.vector(intersect(signifgenes_down$SYMBOL, IDHmut_up_5mc_genebody))#35/557
down_IDHmut_up_5mc_TSS <- as.vector(intersect(signifgenes_down$SYMBOL, IDHmut_up_5mc_TSSpromoter))#20/219



#overlap with Diff 5hmC GB (TSS didn't have much overlap)
#tried FC of 2 and 1
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="IDH DE genes overlap with Diff 5hmC", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% down_IDHmut_down_5hmc_GB | rownames(res) %in% up_IDHmut_up_5hmc_GB), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, rownames(res) %in% up_IDHmut_up_5hmc_TSS), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="IDH DE genes overlap with Diff 5mC", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% down_IDHmut_down_5mc_GB | rownames(res) %in% up_IDHmut_up_5mc_GB), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, rownames(res) %in% down_IDHmut_down_5mc_TSS | rownames(res) %in% up_IDHmut_up_5mc_TSS), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))

#table of all overlapping genes
IDH_DE_5mc_overlap<- cbind(up_IDHmut_up_5hmc_GB, up_IDHmut_up_5hmc_TSS, up_IDHmut_up_5mc_GB, up_IDHmut_up_5mc_TSS, down_IDHmut_down_5hmc_GB, down_IDHmut_down_5mc_GB, down_IDHmut_down_5mc_TSS)
write.csv(IDH_DE_5mc_overlap, file = "IDH_DE_5mc_overlap.csv")
save(IDHmut_up_5hmc_GB, IDHmut_up_5hmc_TSS, IDHmut_up_5mc_genebody, IDHmut_up_5mc_TSSpromoter, IDHmut_down_5hmc_GB, IDHmut_down_5mc_genebody, IDHmut_down_5mc_TSSpromoter, file = "Diff_5hmC_genes.RData")

####5hmc genes that overlap with  module genes
IDHmut_up_5hmc_GB_blue <- as.vector(intersect(IDHmut_up_5hmc_GB, blu_mod$gene_id)) #62
IDHmut_up_5hmc_GB_turq <- as.vector(intersect(IDHmut_up_5hmc_GB, t_mod$gene_id)) #66
IDHmut_up_5hmc_GB_brown <- as.vector(intersect(IDHmut_up_5hmc_GB, br_mod$gene_id)) #34
IDHmut_up_5hmc_GB_green <- as.vector(intersect(IDHmut_up_5hmc_GB, g_mod$gene_id)) #21
IDHmut_up_5hmc_GB_yellow <- as.vector(intersect(IDHmut_up_5hmc_GB, y_mod$gene_id))#27

IDHmut_up_5hmc_TSS_blue <- as.vector(intersect(IDHmut_up_5hmc_TSS, blu_mod$gene_id)) #40
IDHmut_up_5hmc_TSS_turq <- as.vector(intersect(IDHmut_up_5hmc_TSS, t_mod$gene_id)) #23
IDHmut_up_5hmc_TSS_brown <- as.vector(intersect(IDHmut_up_5hmc_TSS, br_mod$gene_id)) #22
IDHmut_up_5hmc_TSS_green <- as.vector(intersect(IDHmut_up_5hmc_TSS, g_mod$gene_id)) #11
IDHmut_up_5hmc_TSS_yellow <- as.vector(intersect(IDHmut_up_5hmc_TSS, y_mod$gene_id))#6

IDHmut_up_5hmc_GB_blue <- as.vector(intersect(IDHmut_up_5hmc_GB, blu_mod$gene_id)) #62
IDHmut_up_5hmc_GB_turq <- as.vector(intersect(IDHmut_up_5hmc_GB, t_mod$gene_id)) #66
IDHmut_up_5hmc_GB_brown <- as.vector(intersect(IDHmut_up_5hmc_GB, br_mod$gene_id)) #34
IDHmut_up_5hmc_GB_green <- as.vector(intersect(IDHmut_up_5hmc_GB, g_mod$gene_id)) #21
IDHmut_up_5hmc_GB_yellow <- as.vector(intersect(IDHmut_up_5hmc_GB, y_mod$gene_id))#27

####5mc genes that overlap with  module genes
IDHmut_up_5mc_genebody_blue <- as.vector(intersect(IDHmut_up_5mc_genebody, blu_mod$gene_id)) #38
IDHmut_up_5mc_genebody_turq <- as.vector(intersect(IDHmut_up_5mc_genebody, t_mod$gene_id)) #47
IDHmut_up_5mc_genebody_brown <- as.vector(intersect(IDHmut_up_5mc_genebody, br_mod$gene_id)) #39
IDHmut_up_5mc_genebody_green <- as.vector(intersect(IDHmut_up_5mc_genebody, g_mod$gene_id)) #15
IDHmut_up_5mc_genebody_yellow <- as.vector(intersect(IDHmut_up_5mc_genebody, y_mod$gene_id))#4

down_IDHmut_up_5mc_GB_blue <- as.vector(intersect(down_IDHmut_up_5mc_GB, blu_mod$gene_id)) #3
down_IDHmut_up_5mc_GB_turq <- as.vector(intersect(down_IDHmut_up_5mc_GB, t_mod$gene_id)) #10/35
#[1] "FLNC"     "NES"      "DES"      "LAMB1"    "TCEA3"    "ISYNA1"   "RARRES2"  "LDHA"     "MARVELD1" "PODXL" 
down_IDHmut_up_5mc_GB_brown <- as.vector(intersect(down_IDHmut_up_5mc_GB, br_mod$gene_id)) #5/35
#[1] "DMRTA2" "NXPH4"  "PDGFA"  "SSBP4"  "MXRA7"
down_IDHmut_up_5mc_GB_green <- as.vector(intersect(down_IDHmut_up_5mc_GB, g_mod$gene_id)) #0
down_IDHmut_up_5mc_GB_yellow <- as.vector(intersect(down_IDHmut_up_5mc_GB, y_mod$gene_id))#1/35

#######look for overlap with reactive gene sets
library(readxl)
reactive<- read_xlsx("C:/Users/caitl/Documents/third rotation/astro_reactivity_IDH/Astrocyte_Genesets.xlsx", sheet = 1, col_names = T)
as.vector(intersect(signifgenes_up_astro$SYMBOL, reactive$`Inflammatory Genes`)) #0 genes
as.vector(intersect(signifgenes_down_astro$SYMBOL, reactive$`Inflammatory Genes`)) #5 gene
as.vector(intersect(signifgenes_up_astro$SYMBOL, reactive$`Alternative Activation`))#1 genes
as.vector(intersect(signifgenes_down_astro$SYMBOL, reactive$`Alternative Activation`))#20 genes
as.vector(intersect(signifgenes_up_astro$SYMBOL, reactive$Progenitor))#2 genes
as.vector(intersect(signifgenes_down_astro$SYMBOL, reactive$Progenitor))#4 genes
as.vector(intersect(signifgenes_up_astro$SYMBOL, reactive$Mature))#9 genes
as.vector(intersect(signifgenes_down_astro$SYMBOL, reactive$Mature))#0 genes


with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes IDHmt vs wt Inflammatory genes", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% reactive$`Inflammatory Genes`), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes IDHmt vs wt Alternative activation genes", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% reactive$`Alternative Activation`), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes IDHmt vs wt progenitor genes", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% reactive$Progenitor), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes IDHmt vs wt mature genes", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% reactive$Mature), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))


#############try IDHwt vs IDHmt in bulk data from Yangping
##try with DESEQ2
library(DESeq2)
library(ggplot2)
setwd("~/third rotation/5HMC and 5MC collab")

#data with recurrent (not used for the IDH comparison)
countData<- read.csv("RNA_counts2.csv", header = TRUE, sep = ",")
metaData<- read.csv("RNA_meta.csv", header = TRUE, sep = ",")

df1<- countData[!duplicated(countData$Gene), ]
dds <- DESeqDataSetFromMatrix(countData=df1, 
                              colData=metaData, 
                              design=~IDH, tidy = TRUE)
#if you get an error saing that the count data isn't all intergers, format cells in excel to numeric with 0 sig figs

dds <- DESeq(dds)
res <- results(dds)
head(results(dds, tidy=TRUE))
summary(res) 
res <- res[order(res$padj),]
head(res)
IDHmtvswt_bulk<- as.data.frame(res)
IDHmtvswt_bulk$gene <- rownames(IDHmtvswt_bulk)

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="IDH")
plotPCA(vsdata, intgroup="category")

plotCounts(dds, gene="PTTG1", intgroup="IDH")
plotCounts(dds, gene="TOP2A", intgroup="IDH")
plotCounts(dds, gene="ETNPPL", intgroup="IDH")
plotCounts(dds, gene="ALDOC", intgroup="IDH")

signifgenes <- IDHmtvswt_bulk %>% filter(padj < 0.05) #2678 genes
signifgenes_down_bulk <- signifgenes %>% filter(log2FoldChange < -2) #1897 genes up in wt tumor
signifgenes_up_bulk <- signifgenes %>% filter(log2FoldChange > 2) #741 genes up in IDHmt tumor

signifgenes_up$genes <- rownames(signifgenes_up)
signifgenes_down$genes <- rownames(signifgenes_down)
signifgenes$genes <- rownames(signifgenes)

#Genes that are up in the IDH tumor bulk data (padj < 0.05; FC > 2) and have a gain of 5hmc in the IDH mut tumor
genes_5hmc_gain <- read_xlsx("Copy of Genelist.xlsx", sheet = 1, col_names = T)
as.vector(intersect(signifgenes_up$genes, genes_5hmc_gain$gene))
# [1] "NEGR1"      "ARPP21"     "RBFOX3"     "MTAP"       "FOXP2"      "CDH18"      "TRAF3IP2"   "PHACTR3"   
# [9] "AFF3"       "MKLN1"      "TRHDE"      "DLGAP1-AS1" "GRIN1"      "PITPNM3"    "PCDH15"     "SMAD9"     
# [17] "CNTLN"      "KLLN"       "DLEU1"      "CSMD3"      "C12orf49"

##Genes that are up in the wt tumor bulk data (padj < 0.05; FC < -2) and have a gain of 5hmc in the IDH mut tumor
as.vector(intersect(signifgenes_down$genes, genes_5hmc_gain$gene))
#[1] "LDHA"         "SERPINH1"     "FSTL1"        "VIM"          "SFRP4"        "TNFRSF1A"     "VCL"         
#[8] "SMIM3"        "TNFRSF12A"    "GLIS3"        "PTX3"         "TMSB4X"       "MYH9"         "RHOJ"        
#[15] "VASP"         "CD9"          "TGIF1"        "ZNF521"       "CDH11"        "FOSL2"        "RHBDF1"      
#[22] "WNT5A"        "ACSS3"        "SAMD9L"       "FNDC3B"       "EFEMP1"       "IQGAP1"       "STBD1"       
#[29] "AK4"          "LGALS1"       "IER5L"        "CLU"          "EPHA2"        "TMEM87B"      "GPC4"        
#[36] "TMEM221"      "PYGL"         "RAB34"        "SH2D4A"       "MRC2"         "FN1"          "CLCF1"       
#[43] "WWTR1"        "NDRG1"        "CSRP2"        "SEMA6D"       "CEBPD"        "CYFIP1"       "VAT1"        
#[50] "BBC3"         "GJA1"         "TGFB2"        "PREX1"        "DLC1"         "HSPB1"        "LYPD1"       
#[57] "ARNTL"        "LAMA4"        "MT2A"         "TMSB10"       "ITGA7"        "TMOD1"        "GAP43"       
#[64] "PLEKHF2"      "ANXA5"        "LMO2"         "KLF9"         "HEG1"         "EHD4"         "CXCR4"       
#[71] "NES"          "MAN2A1"       "ESYT1"        "ZC3HAV1"      "ABCA1"        "EYA2"         "ADAMTS1"     
#[78] "MGLL"         "SYTL2"        "FBXL7"        "LOC100126784" "CHST14"       "PLXNB2"       "SEC24D"      
#[85] "FAM184A"      "JAG1"         "PXDC1"        "SOX9"         "LHX2"         "MAP3K5"       "B3GNT2"      
#[92] "PDE4B"        "ROBO2"        "FHL3"         "KLF16"        "NRP2"         "SSBP4"        "RASSF3"      
#[99] "ITGA2"        "SLC35G2"      "ITGA6"        "CEP112"       "RRAS2"        "SLC37A1"      "EPS15"       
#[106] "TTYH3"        "PIEZO1"       "CDC42EP4"     "COBL"         "RIN1"         "ARHGEF6"      "YBX3"        
#[113] "EPAS1"        "CMTM3"        "PTK2B"        "SORCS1"       "DAP"          "NID1"  

genes_5hmc_loss <- read_xlsx("Copy of Genelist.xlsx", sheet = 2, col_names = T)
as.vector(intersect(signifgenes_down$genes, genes_5hmc_loss$gene))
#[1] "SEC61G" "RGS6"   "HOXB3"  "SEMA3A" "MEOX2"  "PSPH" 


####intersection of bulk and astro-specific differentially expressed genes
up_overlap<- as.vector(intersect(signifgenes_up_bulk$gene, signifgenes_up_astro$SYMBOL)) #277 overlap
#[1] "GNAL"        "NEGR1"       "ZBTB16"      "GABBR1"      "PCDH11X"     "PCBP3"       "CACNG2"      "FERMT1"      "TPTEP1"     
#[10] "FAM222A-AS1" "KSR2"        "WSCD2"       "ATOH8"       "CDR1"        "KCNN3"       "LINC00925"   "ACSL6"       "KIAA1755"   
#[19] "LOC729732"   "SLITRK1"     "ACBD7"       "SMOC1"       "FAM57B"      "CBLN1"       "CSMD3"       "FAM13C" 
down_overlap<- as.vector(intersect(signifgenes_down_bulk$gene, signifgenes_down_astro$SYMBOL))
write.csv(up_overlap, file = "up_overlap.csv")
write.csv(down_overlap, file = "down_overlap.csv")





########################compare IDH1mt to margin
library(DESeq2)
library(ggplot2)

countData <- read.csv('C:/Users/caitl/Documents/third rotation/GBM/RNA/IDHonly_counts.csv', header = TRUE, sep = ",")
metaData <- read.csv('C:/Users/caitl/Documents/third rotation/GBM/RNA/IDHonly_meta.csv', header = TRUE, sep = ",")

df1<- countData[!duplicated(countData$Gene), ]

dds <- DESeqDataSetFromMatrix(countData=df1, 
                              colData=metaData, 
                              design=~1, tidy = TRUE)
design(dds) <- ~ Patient + Tissue


dds <- DESeq(dds)
res <- results(dds)
head(results(dds, tidy=TRUE))
summary(res) 
res <- res[order(res$padj),]
head(res)
IDHmtvsM<- as.data.frame(res)
signifgenes_IDHmtVsMargin <- IDHmtvsM  %>% filter(padj<.05 & abs(log2FoldChange)>2)
write.csv(signifgenes_IDHmtVsMargin, file = "signifgenes_IDHmtVsMargin.csv")
write.csv(IDHmtvsM, "IDHmtvsM.csv")
saveRDS(IDHmtvsM, file = "DESEQ_IDH1mtVSmargin.Rds")

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="IDH1mt tumor vs Margin", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
#number of genes up in margin and IDH1mt: 
signifgenes <- IDHmtvsM %>% filter(padj < 0.05) #1550 genes
signifgenes_down <- signifgenes %>% filter(log2FoldChange < -2) #749 genes up in margin
signifgenes_up <- signifgenes %>% filter(log2FoldChange > 2) #146 genes up in IDH1mt tumor

Aastro_genes <- na.omit(cell_type$A.astro)
Fastro_genes <- na.omit(cell_type$F.astro)

as.vector(intersect(Fastro_genes, rownames(signifgenes_down)))
#[1] "LHX2"     "EFNB2"    "MOXD1"    "CENPV"    "LDLR"     "FABP5"    "RAI14"    "DPYSL5"   "TYMS"     "LRIG3"   
#[11] "PLAGL1"   "ITGA2"    "CDH4"     "EPHB2"    "LLGL1"    "DFNB31"   "PCDHGA12" "DEPDC1B"  "PCDHGA11" "PDLIM4"  
#[21] "LRRC61"   "PLEKHF1"
margin_adult_genes<- as.vector(intersect(Aastro_genes, rownames(signifgenes_down))) #of the genes up in margin, a lot are "mature" astro genes
saveRDS(margin_adult_genes, file = "margin_adult_genes.Rds")
#[1] "S100B"     "S100A13"   "CSRP1"     "GJB6"      "S100A1"    "PTGDS"     "SLC39A12"  "WIF1"      "FGFR2"    
#[10] "ACAA2"     "PAQR8"     "CLDN10"    "HOPX"      "SLCO1C1"   "GABRA2"    "ELOVL2"    "LGI1"      "SASH1"    
#[19] "DIO2"      "FGFR3"     "GRM3"      "PPP1R1B"   "LRRC16A"   "MT1E"      "LRRC3B"    "MTMR10"    "ALDH1A1"  
#[28] "GPC5"      "GSTM1"     "HHATL"     "SLC48A1"   "GSTM5"     "MAOB"      "CABLES1"   "SELENBP1"  "EMX2OS"   
#[37] "HSD17B6"   "MYO6"      "RGN"       "NQO1"      "RGS20"     "LGALS3"    "CYP2J2"    "EMX2"      "ME1"      
#[46] "RYR3"      "HEPH"      "CHI3L1"    "EGLN3"     "IQCA1"     "ADORA2B"   "MEGF10"    "MYLK"      "PRSS35"   
#[55] "C3orf70"   "SLC13A5"   "FXYD1"     "DLC1"      "IDI2-AS1"  "OAF"       "MT1M"      "ITPKB"     "ARRB1"    
#[64] "FBXO2"     "NTSR2"     "RERG"      "MT1F"      "ABCD2"     "MT1G"      "CYP7B1"    "C16orf89"  "SLC7A10"  
#[73] "SGCD"      "PYGM"      "ACSS3"     "MYBPC1"    "PIR"       "GPR143"    "RLBP1"     "FBLN1"     "SYT9"     
#[82] "PPP2R2C"   "CDS1"      "GNA14"     "CPNE6"     "DNAH7"     "FAS"       "EPHX2"     "ACACB"     "CNTNAP3"  
#[91] "AIFM3"     "KIAA0930"  "GRIN2C"    "CHST1"     "LINC00092" "PCDHGA3"   "ARSD"      "SORCS2"    "TMEM229A" 
#[100] "EFHC2"     "ARSF"      "TMPRSS3"   "AQP1"      "BDH1"      "SLC13A3"   "LGI4"      "BHMT2"     "TMEM132C" 
#[109] "ACOT11"    "RAPGEF3"   "ADORA1"    "PYROXD2"   "DBX2"      "HSPB2"     "HDHD3"     "ANKRD35"   "RHCG"     
#[118] "ASB4"      "SOD3"      "TRPV3"     "PI16"      "VASN"

###genes to plot separately: WIF1, HOPX, GABRA2?, S100B, GJB6, ALDH1A1

vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="Tissue", ntop=20000)
plotPCA(vsdata, intgroup="Patient")


plotCounts(dds, gene="PTTG1", intgroup="IDH")
plotCounts(dds, gene="TOP2A", intgroup="IDH")
plotCounts(dds, gene="ETNPPL", intgroup="IDH")
plotCounts(dds, gene="ALDOC", intgroup="IDH")

#try looking at reactive gene sets
reactive<- read_xlsx("C:/Users/caitl/Documents/third rotation/astro_reactivity_IDH/Astrocyte_Genesets.xlsx", sheet = 1, col_names = T)
as.vector(intersect(signifgenes_up$SYMBOL, reactive$`Inflammatory Genes`)) #0 genes
as.vector(intersect(signifgenes_down$SYMBOL, reactive$`Inflammatory Genes`)) #5 gene
as.vector(intersect(signifgenes_up$SYMBOL, reactive$`Alternative Activation`))#1 genes
as.vector(intersect(signifgenes_down$SYMBOL, reactive$`Alternative Activation`))#20 genes
as.vector(intersect(signifgenes_up$SYMBOL, reactive$Progenitor))#2 genes
as.vector(intersect(signifgenes_down$SYMBOL, reactive$Progenitor))#4 genes
as.vector(intersect(signifgenes_up$SYMBOL, reactive$Mature))#9 genes
as.vector(intersect(signifgenes_down$SYMBOL, reactive$Mature))#0 genes


with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes IDHmt vs margin Inflammatory genes", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% reactive$`Inflammatory Genes`), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes IDHmt vs margin Alternative activation genes", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% reactive$`Alternative Activation`), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes IDHmt vs margin progenitor genes", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% reactive$Progenitor), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes IDHmt vs margin mature genes", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% reactive$Mature), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))


#######################

########################compare IDHwt to margin
#remove organoid, fetal/adult, recurrent, and IDHmt

##try with DESEQ2
library(DESeq2)
library(ggplot2)
setwd("~/third rotation/GBM/RNA")
countData <- read.csv('GBM_counts.csv', header = TRUE, sep = ",")
metaData <- read.csv('GBM_meta.csv', header = TRUE, sep = ",")
df1<- countData[!duplicated(countData$Geneid), ]
dds <- DESeqDataSetFromMatrix(countData=df1, 
                              colData=metaData, 
                              design=~1, tidy = TRUE)
#<https://support.bioconductor.org/p/84241/>
design(dds) <- ~ Patient + Tissue

dds <- DESeq(dds)
res <- results(dds)
head(results(dds, tidy=TRUE))
summary(res) 
res <- res[order(res$padj),]
head(res)
TumorvsMargin<- as.data.frame(res)
signifgenes <- TumorvsMargin %>% filter(padj < 0.05 & abs(log2FoldChange)>2)#1333
signifgenes_down <- signifgenes %>% filter(log2FoldChange < -2) #777
signifgenes_up <- signifgenes %>% filter(log2FoldChange > 2)
write.csv(signifgenes, file = "DESeq_IDHwtVsMargin.csv")

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-15,15)))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


library(clusterProfiler)
library(Rgraphviz)
library(DOSE)

library(org.Hs.eg.db)
down_ann <- select(org.Hs.eg.db,keys=rownames(signifgenes_down),keytype = "SYMBOL" ,columns=c("ENTREZID","SYMBOL","GENENAME"))
head(down_ann)
all_ann <- select(org.Hs.eg.db,keys=rownames(TumorvsMargin),keytype = "SYMBOL" ,columns=c("ENTREZID","SYMBOL","GENENAME"))
allgenes <- as.character(all_ann$ENTREZID)
down_GO <- as.character(down_ann$ENTREZID)
ego_down <- enrichGO(gene          = down_GO,
                     universe      = allgenes,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     readable      = TRUE)
dotplot(ego_down, showCategory=20, orderBy = "count")
up_ann <- select(org.Hs.eg.db,keys=rownames(signifgenes_up),keytype = "SYMBOL" ,columns=c("ENTREZID","SYMBOL","GENENAME"))
up_GO <- as.character(up_ann$ENTREZID)
ego_up <- enrichGO(gene          = up_GO,
                   universe      = allgenes,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   readable      = TRUE)
dotplot(ego_up, showCategory=20, orderBy = "count")

down_kegg <- enrichKEGG(gene         = down_GO,
                        organism     = 'hsa',
                        pvalueCutoff = 0.05)
up_kegg<- enrichKEGG(gene         = up_GO,
                     organism     = 'hsa',
                     pvalueCutoff = 0.05)
dotplot(down_kegg, showCategory=20)

library(ReactomePA)
react_down<- enrichPathway(gene=down_GO, pvalueCutoff = 0.05, readable=TRUE)
react_up<- enrichPathway(gene=up_GO, pvalueCutoff = 0.05, readable=TRUE)


#DEGs heatmap
rownames(df1) <- df1$Geneid
df2 <- df1[,-1]
df2<- cpm(df2, log = T)
df3 <- df2[rownames(signifgenes),c(15,6,17,2,21,19,13,11,4,8,9,7,18,5,1,14,12,10,20,16,3)]
library(ComplexHeatmap)
ComplexHeatmap::pheatmap(as.matrix(df3), name = "mat",cluster_rows = TRUE, cluster_cols = FALSE,
                         color =colorRampPalette(brewer.pal(n=9,name="Blues"), bias = 0.8)(100), scale = "row")


vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="Tissue")
plotPCA(vsdata, intgroup="Patient")

plotCounts(dds, gene="PTTG1", intgroup="Tissue")
plotCounts(dds, gene="ETNPPL", intgroup="Tissue")
plotCounts(dds, gene="ALDOC", intgroup="Tissue")
plotCounts(dds, gene="NFIA", intgroup="Tissue")
plotCounts(dds, gene="NR3C2", intgroup="Tissue")
plotCounts(dds, gene="OLIG2", intgroup="Tissue")
plotCounts(dds, gene="OLIG1", intgroup="Tissue")
plotCounts(dds, gene="SOX10", intgroup="Tissue")
plotCounts(dds, gene="RXRG", intgroup="Tissue")
plotCounts(dds, gene="SOX3", intgroup="Tissue")



#look at organoid modules within DESEQ results
TumorvsMargin$SYMBOL <- rownames(TumorvsMargin)
signifgenes <- TumorvsMargin %>% filter(padj < 0.05) #1550 genes
signifgenes_down <- signifgenes %>% filter(log2FoldChange < -2) #777
signifgenes_up <- signifgenes %>% filter(log2FoldChange > 2) #556

#save IDH DESEQ R objects
save(countData, metaData, df1, dds, res, TumorvsMargin, vsdata, signifgenes, signifgenes_down, signifgenes_up, file="GBM_DESeq.RData")
#genes that are down and in various modules
down_blue <- as.vector(intersect(signifgenes_down$SYMBOL, blu_mod$gene_id)) #72
down_turq <- as.vector(intersect(signifgenes_down$SYMBOL, t_mod$gene_id)) #99
down_brown <- as.vector(intersect(signifgenes_down$SYMBOL, br_mod$gene_id)) #35
down_green <- as.vector(intersect(signifgenes_down$SYMBOL, g_mod$gene_id)) #29
down_yellow <- as.vector(intersect(signifgenes_down$SYMBOL, y_mod$gene_id)) #31
up_blue <- as.vector(intersect(signifgenes_up$SYMBOL, blu_mod$gene_id)) #27
up_turq <- as.vector(intersect(signifgenes_up$SYMBOL, t_mod$gene_id)) #107
up_brown <- as.vector(intersect(signifgenes_up$SYMBOL, br_mod$gene_id)) #12
up_green <- as.vector(intersect(signifgenes_up$SYMBOL, g_mod$gene_id)) #33
up_yellow <- as.vector(intersect(signifgenes_up$SYMBOL, y_mod$gene_id)) #25
save(down_blue, down_brown, down_green, down_turq, down_yellow, up_blue, up_brown, up_green, up_turq, up_yellow, file = "GBM_DEG_moduleOverlap.RData")
##calculating enrichment of these modules
#can be done two ways? Try both
#option 1: of the early module genes that are DE, % in tumor and % in margin
#option 2: of the genes up in margin, what % is in early module (and vice versa for tumor)
#option 1. first, add the corresponding down/up
#down_blue + up_blue = 99; turq = 206; brown = 47; green = 62; yellow = 56
tissue <- rep(c("margin", "tumor"), times=5)
module_total <- rep(c(99,206,47,62,56), each=2)
module_values <- c(-72,27,-99,107,-35,12,-29,33,-31,25)
modules <- rep(c("Middle_Late", "Early", "Middle", "Late", "Early_Middle"), each=2)
fake_groups <- rep(c("a","b","c","d","e"), each=2)
tissue_mod_df <- data.frame(tissue, module_total, module_values, modules, fake_groups)
tissue_mod_df <- tissue_mod_df %>% mutate(freq=(module_values/module_total)*100)

#option 2. what are the total genes up in margin and tumor?
#up in margin = 777; up in tumor = 556
DE_tissue <- rep(c(777,556), times=5)
tissue_mod_df2 <- data.frame(tissue, DE_tissue, module_values, modules)
tissue_mod_df2 <- tissue_mod_df2 %>% mutate(freq=(module_values/DE_tissue)*100)

#option 2 is ugly (values are small), so let's go with option 1 for now.
#make a diverging bar chart for each module of the % module genes enriched
library(ggplot2)
ggplot(tissue_mod_df, aes(x = fake_groups, y = freq)) +
  geom_bar(stat = "identity") + scale_y_continuous(limits = c(-100, 100)) +
  coord_flip() +
  facet_wrap(factor(modules, levels = unique(modules))~.)




down_all <- as.vector(intersect(signifgenes_down$SYMBOL, submod$gene_id))
up_all <- as.vector(intersect(signifgenes_up$SYMBOL, submod$gene_id))
with(lrt.res$table, plot(logFC, -log10(PValue), pch=20, main="DE genes Tumor vs Margin Module overlap", xlim=c(-12,12), ylim=c(0,65)))
with(subset(lrt.res$table, FDR<.05 & logFC < -2), points(logFC, -log10(PValue), pch=20, col="#C1C7C9"))
with(subset(lrt.res$table, FDR<.05 & logFC>2), points(logFC, -log10(PValue), pch=20, col="#C1C7C9"))
with(subset(lrt.res$table, SYMBOL %in% down_all | SYMBOL %in% up_all), points(logFC, -log10(PValue), pch=20, col="#549b96"))
#lots of difference between tumor and margin, but roughly 1/3 of DE genes fall into one of my astro maturation modules ()
#math for above: ((down_all + up_all)/(signifgenes_down + signifgenes_up)) 

#let's see if this is specific to this set of genes
website_data<- read_xlsx("C:/Users/caitl/Documents/third rotation/Copy of website_data_4_2019 (version 1).xlsx", sheet = "FC>1", col_names = TRUE)
website_sample<- website_data[sample(nrow(website_data), 3944, replace = FALSE, prob = NULL),]
down_website <- as.vector(intersect(signifgenes_down$SYMBOL, website_sample$genes)) #77, 67, 65
up_website <- as.vector(intersect(signifgenes_up$SYMBOL, website_sample$genes)) #28, 25, 21 
#105 total out of (695+403)=1098, so roughly 9.5%, 8%, 8%

website_data<- read_xlsx("C:/Users/caitl/Documents/third rotation/Copy of website_data_4_2019 (version 1).xlsx", sheet = "Sheet6", col_names = TRUE)
fetal_website <- website_data[website_data$fetal>2,1]
adult_website <- website_data[website_data$adult>2,1]
website_both <- unique(c(fetal_website$...1, adult_website$...1))
website_both2<- sample(website_both, 3944, replace = FALSE, prob = NULL)
down_website <- as.vector(intersect(signifgenes_down$SYMBOL, website_both2)) #107, 98, 104
up_website <- as.vector(intersect(signifgenes_up$SYMBOL, website_both2)) #35, 31, 35
#13%, 11.7%, 12.7%



#Make volcano plots for following modules: blue, turquoise, yellow, brown, green

#blue
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes T vs M Middle/Late Module"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% down_blue | rownames(res) %in% up_blue), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))

#turquoise
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes T vs M Early Module"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% down_turq | rownames(res) %in% up_turq), points(log2FoldChange, -log10(pvalue), pch=20, col="#30D5C8"))


#Brown
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes T vs M Middle Module"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% down_brown | rownames(res) %in% up_brown), points(log2FoldChange, -log10(pvalue), pch=20, col="brown"))


#Green
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes T vs M Late Module"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% down_green | rownames(res) %in% up_green), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))


#Yellow
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes T vs M Early/Middle Module"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% down_yellow | rownames(res) %in% up_yellow), points(log2FoldChange, -log10(pvalue), pch=20, col="#FFD700"))


#overlay fetal and mature astro genes on volcano
Fastro_genes <- na.omit(cell_type$F.astro)
Aastro_genes <- na.omit(cell_type$A.astro)

down_fetal <- as.vector(intersect(signifgenes_down$SYMBOL, Fastro_genes))
down_adult <- as.vector(intersect(signifgenes_down$SYMBOL, Aastro_genes))
up_fetal <- as.vector(intersect(signifgenes_up$SYMBOL, Fastro_genes))
up_adult <- as.vector(intersect(signifgenes_up$SYMBOL, Aastro_genes))

sp_down_adult <- as.vector(c("GJB6", "WIF1", "SLC1A2", "PTGDS", "SLC14A1", "S100A1", "FAM171B", "GABRA2", "FGFR2", "GRM3"))
sp_up_fetal <- as.vector(c("KIF23", "HIST1H2AL", "HIST1H2BH", "TNC", "CD24", "GPX3", "PTX3", "VCAM1", "IGF2BP2", "DTL"))

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes T vs M Fetal/Adult genes"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% down_fetal | rownames(res) %in% up_fetal), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, rownames(res) %in% sp_up_fetal), text(log2FoldChange, -log10(padj), labels=subset(rownames(res), rownames(res) %in% sp_up_fetal), cex=0.4, pos=3))
with(subset(res, rownames(res) %in% sp_up_fetal), points(log2FoldChange, -log10(pvalue), pch=20, col="purple"))
with(subset(res, rownames(res) %in% down_adult | rownames(res) %in% up_adult), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, rownames(res) %in% sp_down_adult), text(log2FoldChange, -log10(padj), labels=subset(rownames(res), rownames(res) %in% sp_down_adult), cex=0.4, pos=3))
with(subset(res, rownames(res) %in% sp_down_adult), points(log2FoldChange, -log10(pvalue), pch=20, col="purple"))

save(Fastro_genes,Aastro_genes,down_fetal,down_adult,up_fetal,up_adult, file = "GBM_DESEQ_feta_adult_overlap.RData")


#try looking at reactive gene sets
reactive<- read_xlsx("C:/Users/caitl/Documents/third rotation/astro_reactivity_IDH/Astrocyte_Genesets.xlsx", sheet = 1, col_names = T)
as.vector(intersect(signifgenes_up$SYMBOL, reactive$`Inflammatory Genes`)) #0 genes
as.vector(intersect(signifgenes_down$SYMBOL, reactive$`Inflammatory Genes`)) #5 gene
as.vector(intersect(signifgenes_up$SYMBOL, reactive$`Alternative Activation`))#1 genes
as.vector(intersect(signifgenes_down$SYMBOL, reactive$`Alternative Activation`))#20 genes
as.vector(intersect(signifgenes_up$SYMBOL, reactive$Progenitor))#2 genes
as.vector(intersect(signifgenes_down$SYMBOL, reactive$Progenitor))#4 genes
as.vector(intersect(signifgenes_up$SYMBOL, reactive$Mature))#9 genes
as.vector(intersect(signifgenes_down$SYMBOL, reactive$Mature))#0 genes


with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes IDHwt vs margin Inflammatory genes", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% reactive$`Inflammatory Genes`), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes IDHwt vs margin Alternative activation genes", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% reactive$`Alternative Activation`), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes IDHwt vs margin progenitor genes", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% reactive$Progenitor), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="DE genes IDHwt vs margin mature genes", xlim=c(-15,15), ylim=c(0,30)))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="#C1C7C9"))
with(subset(res, rownames(res) %in% reactive$Mature), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))

###############correlate GBM and organoid samples using only module genes
setwd("~/third rotation/Organoid maturation timeline/HiSeq/RNA-seq")

early_top100<- read.csv("C:/Users/caitl/Documents/third rotation/Organoid maturation timeline/HiSeq/RNA-seq/WGCNA/early_top100.csv", header = T)
earlyMiddle_top100<- read.csv("C:/Users/caitl/Documents/third rotation/Organoid maturation timeline/HiSeq/RNA-seq/WGCNA/earlyMiddle_top100.csv", header = T)
middle_top100<- read.csv("C:/Users/caitl/Documents/third rotation/Organoid maturation timeline/HiSeq/RNA-seq/WGCNA/Middle_top100.csv", header = T)
middleLate_top100<- read.csv("C:/Users/caitl/Documents/third rotation/Organoid maturation timeline/HiSeq/RNA-seq/WGCNA/MiddleLate_top100.csv", header = T)
late_top100<- read.csv("C:/Users/caitl/Documents/third rotation/Organoid maturation timeline/HiSeq/RNA-seq/WGCNA/Late_top100.csv", header = T)
top_module <- rbind(early_top100,earlyMiddle_top100,middle_top100,middleLate_top100,late_top100)
top_module<- top_module[,-1]

early_top50<- read.csv("C:/Users/caitl/Documents/third rotation/Organoid maturation timeline/HiSeq/RNA-seq/WGCNA/early_top50.csv", header = T)
earlyMiddle_top50<- read.csv("C:/Users/caitl/Documents/third rotation/Organoid maturation timeline/HiSeq/RNA-seq/WGCNA/earlyMiddle_top50.csv", header = T)
middle_top50<- read.csv("C:/Users/caitl/Documents/third rotation/Organoid maturation timeline/HiSeq/RNA-seq/WGCNA/Middle_top50.csv", header = T)
middleLate_top50<- read.csv("C:/Users/caitl/Documents/third rotation/Organoid maturation timeline/HiSeq/RNA-seq/WGCNA/MiddleLate_top50.csv", header = T)
late_top50<- read.csv("C:/Users/caitl/Documents/third rotation/Organoid maturation timeline/HiSeq/RNA-seq/WGCNA/Late_top50.csv", header = T)
top_module <- rbind(early_top50,earlyMiddle_top50,middle_top50,middleLate_top50,late_top50)
top_module<- top_module[,-1]

all_module <- read.csv("C:/Users/caitl/Documents/third rotation/Organoid maturation timeline/HiSeq/RNA-seq/WGCNA/org_RNAseq_gene_modules.csv", header = T)

library(readxl)
cts <- read_xlsx("wgcna input.xlsx", col_names = TRUE) 
cts<- cts[!duplicated(cts$Geneid), ]

#remove rows for which there is no approved gene symbol (df2), make gene name the row name (df2), remove the gene name column (df3).
cts <- cts[!is.na(cts$`C4 d550`), ]
cts2 <- as.data.frame(cts)
rownames(cts2) <- cts2$Geneid

#load data
library(dplyr)
feature_counts <- read_xlsx("~/third rotation/GBM/RNA/GBM feature counts.xlsx")
#removing duplicates genes 
df1<- feature_counts[!duplicated(feature_counts$Geneid), ]
#remove rows for which there is no approved gene symbol (df2), make gene name the row name (df2), remove the gene name column (df3).
df2 <- df1[!is.na(df1$Geneid), ]
df2 <- as.data.frame(df2)
row.names(df2) <- df2$Geneid
#tpm before filtering columns
df2 <- df2[,-c(1,2)]
#normalize all using TPM method
df2_tpm <- tpm3(df2, gene_length)


df3 <- df2_tpm[,-c(1:3,38)]
#create a tumor object and a margin object (only IDHwt)
df_t <- df3[,c(1:6,10,12,14,16,17)]
df_m <- df3[,c(19:23,27,28,30,32:33)]

df3 <- df2[,-c(2,12:14,29:31)]
df4 <- df3[,-c(12,14,16,19,27,29,32)]
gbm_filtered<- df4 %>% filter(df4$Geneid %in% top_module)
gbm_filtered <- gbm_filtered[,-c(2:4,26)]


###corelate with early, middle, late
org_cts<- read_xlsx("C:/Users/caitl/Documents/third rotation/Organoid maturation timeline/HiSeq/RNA-seq/master_featurecounts.xlsx", col_names = T)
org_cts<- org_cts[!duplicated(org_cts$Geneid), ]
org_cts <- org_cts[!is.na(org_cts$Geneid), ]
gene_length <- as.vector(org_cts$Length)
org_cts <- as.data.frame(org_cts)
rownames(org_cts) <- org_cts$Geneid
org_cts <- org_cts[,-c(1,2)]
#normalize all using TPM method
tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}
org_tpm <- tpm3(org_cts, gene_length)
org_tpm <- as.data.frame(org_tpm)

##filter tpm matrices to only include module genes
org_tpm$Geneid <- rownames(org_tpm)
org_tpm_filtered<- org_tpm %>% filter(org_tpm$Geneid %in% df_m_filtered$Geneid)

df_m <- as.data.frame(df_m)
df_m$Geneid <- rownames(df_m)
df_m_filtered<- df_m %>% filter(df_m$Geneid %in% all_module$gene_id)

df_t <- as.data.frame(df_t)
df_t$Geneid <- rownames(df_t)
df_t_filtered<- df_t %>% filter(df_t$Geneid %in% all_module$gene_id)
#remove gene id column
org_tpm_filtered <- org_tpm_filtered[,-c(21,22)]
df_m_filtered <- df_m_filtered[,-11]
df_t_filtered <- df_t_filtered[,-12]


df_t_filtered<- df_t %>% filter(df_t$Geneid %in% top_module)
df_t_filtered <- df_t_filtered[,-1]
df_m_filtered<- df_m %>% filter(df_m$Geneid %in% top_module)
df_m_filtered <- df_m_filtered[,-1]

merged_data <- merge(gbm_filtered, cts2, by.x = "Geneid", by.y = "Geneid")
rownames(merged_data) <- merged_data$Geneid
merged_data <- merged_data[,-1]
merged_matrix <- as.matrix(merged_data)

#merge tumor and margin matrices with separate early, middle, and late matrices
early_df <- org_tpm_filtered[,c(1:7)]
middle_df <- org_tpm_filtered[,c(8:15)]
late_df <- org_tpm_filtered[,c(16:20)]

#to perform correlation, take rowmeans for each df
early_means <- rowMeans(early_df)
middle_means <- rowMeans(middle_df)
late_means <- rowMeans(late_df)
tumor_means <- rowMeans(df_t_filtered)
margin_means <- rowMeans(df_m_filtered)

# Step 2: Perform Spearman correlation analysis
correlation_result <- cor.test(early_means, tumor_means, method = "spearman")

#spearman correlation that compares all samples
library(Hmisc)
cor_result_merged <- round(cor(merged_matrix, method="spearman"),2)

##get rid of unnecessary columns and make that the cor_matrix
library(gplots)
heatmap.2(as.matrix(cor_result_merged[-c(22:42),-c(1:21)]), dendrogram= "none",Colv = F,Rowv = F, trace = NULL,scale = "column",col=colorRampPalette(rev(brewer.pal(n=11,name="RdBu")))(100))
heatmap.2(as.matrix(cor_result_merged[-c(22:42),-c(1:21)]), dendrogram= "none",Colv = F,Rowv = F, trace = NULL,scale = "column",col=colorRampPalette(rev(brewer.pal(n=11,name="RdBu")), bias=1.2)(100))


###########################anaylsis between Steven's data and my tumor data
countData<- read.csv("C:/Users/caitl/Documents/third rotation/GBM/RNA/all_GBM_samples.csv", header = TRUE, sep = ",")
countData2<- read.csv("C:/Users/caitl/Documents/third rotation/feature_counts_human_astros.csv", header = TRUE, sep = ",")
countData3 <- countData2[,-2]
metaData<- read.csv("C:/Users/caitl/Documents/third rotation/GBM/RNA/all_tissue_meta.csv", header = TRUE, sep = ",")

#vst normalization
df1<- countData[!duplicated(countData$Gene), ]
df2<- countData3[!duplicated(countData3$Geneid), ]
tissue_df<- merge(df1, df2, by.x="Gene", by.y="Geneid")
dds <- DESeqDataSetFromMatrix(countData=tissue_df, 
                              colData=metaData, 
                              design=~Tissue, tidy = TRUE)
dds <- DESeq(dds)
vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="Tissue")


#log2CPM
df3 <- as.data.frame(tissue_df)
row.names(df3) <- df3$Gene
df3 <- df3[,-1]
#Filtering count data:
#Step1: filter out the genes for which there are no reads, or there are inconsistent reads across replicate samples, or there are low reads. Chose a normalization technique: CPM (counts per million), for example (GBM_CPM). We filter using CPM values rather than counts because they account for differences in sequencing depth between samples. 
GBM_CPM <- edgeR::cpm(df3)
head(GBM_CPM)
#Step2: Our filtering rule is to keep transcripts that have CPM > X in at least two samples. X is the cutoff value (in CPM) that will consistently yield a minimum transcript count. We need to first determine the scientifically relevant minimum transcript count. As a general rule, a good threshold can be chosen for a CPM value that corresponds to a count of 10 (for 10-20M reads/sample depth). This will vary based on the depth of sequence; for low depth sequence <1M, try a count of 1. CPM = counts per million, or how many counts would I get for a gene if the sample had a library size of 1M. So for a library size of 1M, 1 count = 1 CPM. For a library size of 10M, 10 counts = 1 CPM. For a library size of 20M, 10 counts = 0.5 CPM.
#Step 3: Next, impose the threshold (0.5 is an example value, replace with determined X value).
GBM_thresh <- GBM_CPM > 0.5
#Step 4: Identify/subset the genes for which CPM > 0.5 in at least two samples. Use DGElist to convert conunts.keep df back to an oject (GBM_count) for ease of analysis moving forward.
GBM_keep <- rowSums(GBM_thresh) >= 2
GBM_counts.keep <- df3[GBM_keep,]
GBM_counts.keep2<- cpm(GBM_counts.keep, log = T)
pcaCoord <- muRtools::getDimRedCoords.pca(t(GBM_counts.keep2), components = c(2))
muRtools::getDimRedPlot(pcaCoord, annot=metaData, colorCol="Tissue", shapeCol = "Something", addLabels = F, ptSize = 3)

#do this one more time after filtering for only mature and immature genes to see if PCA looks even better
Aastro_genes<- as.vector(Aastro_genes)
Fastro_genes<- as.vector(Fastro_genes)
astro_genes <- c(Fastro_genes, Aastro_genes)
GBM_counts.keep3 <- as.data.frame(GBM_counts.keep2)
GBM_counts.keep3<- GBM_counts.keep3 %>% filter(rownames(GBM_counts.keep3) %in% astro_genes)
pcaCoord <- muRtools::getDimRedCoords.pca(t(GBM_counts.keep3), components = c(2))
muRtools::getDimRedPlot(pcaCoord, annot=metaData, colorCol="Tissue", shapeCol = "Something", addLabels = F, ptSize = 3)

#make bar plots to show
mature_genes_OI <- as.vector(c("WIF1", "GJB6", "CLDN10", "FGFR2", "PTGDS", "ELOVL2"))
GBM_counts.keep3 <- as.data.frame(GBM_counts.keep2)
GBM_counts.keep3 <- GBM_counts.keep3[,c(20,23,24,27,28,33:52)]
GBM_counts.keep4<- GBM_counts.keep3 %>% filter(rownames(GBM_counts.keep3) %in% Aastro_genes)
GBM_counts.keep5<- GBM_counts.keep3 %>% filter(rownames(GBM_counts.keep3) %in% margin_adult_genes)
GBM_counts.keep6<- GBM_counts.keep3 %>% filter(rownames(GBM_counts.keep3) %in% mature_genes_OI)
signifgenes <- DESEQ_IDH1mtVSmargin %>% filter(padj < 0.05) 
signifgenes_down <- signifgenes %>% filter(log2FoldChange < -2) 
signifgenes_up <- signifgenes %>% filter(log2FoldChange > 2)
GBM_counts.keep7<- GBM_counts.keep3 %>% filter(rownames(GBM_counts.keep3) %in% rownames(signifgenes_up))
GBM_counts.keep6$Gene <- rownames(GBM_counts.keep6)
GBM_counts.keep4$Gene <- rownames(GBM_counts.keep4)
GBM_counts.keep5$Gene <- rownames(GBM_counts.keep5)
GBM_counts.keep7$Gene <- rownames(GBM_counts.keep7)
library(tidyr)
long_GBM_counts<-GBM_counts.keep6 %>%
  pivot_longer(!Gene, names_to = "sample", values_to = "count")
new_value <- as.vector(c("IDH1_tumor", "margin", "IDH1_tumor", "margin", "IDH1_tumor", "margin", "IDH1_tumor", rep("adult", times=12), rep("fetal", times=6)))
new_value <- rep(new_value, times=6)
long_GBM_counts$tissue <- as.factor(new_value)
long_GBM_counts2 <- long_GBM_counts %>%
  group_by(Gene, tissue) %>%
  mutate(tissue_mean = mean(count), tissue_SD = sd(count))
long_GBM_counts2$tissue <- factor(long_GBM_counts2$tissue, levels = c("adult", "margin", "IDH1_tumor", "fetal"))
long_GBM_counts3 <- long_GBM_counts2[,c(1,4,5,6)]
long_GBM_counts3 <- unique(long_GBM_counts3)

library(ggplot2)
ggplot(long_GBM_counts3, aes(x = tissue, y = tissue_mean)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 25)) + geom_errorbar(aes(ymin=tissue_mean-tissue_SD, ymax=tissue_mean+tissue_SD), width=.2) +
  facet_wrap(factor(Gene, levels = unique(Gene))~.)

#dotplot version
ggplot(long_GBM_counts2, aes(x = tissue, y = count)) +
  geom_point(stat = "identity") +
  theme(axis.text.x = element_text(angle = 25)) +
  facet_wrap(factor(Gene, levels = unique(Gene))~.)

#overall mean of margin and mature genes
#first, adult astro genes
long_GBM_counts<-GBM_counts.keep4 %>%
  pivot_longer(!Gene, names_to = "sample", values_to = "count")
new_value <- as.vector(c("IDH1_tumor", "margin", "IDH1_tumor", "margin", "IDH1_tumor", "margin", "IDH1_tumor", rep("adult", times=12), rep("fetal", times=6)))
new_value <- rep(new_value, times=917)
long_GBM_counts$tissue <- as.factor(new_value)
long_GBM_counts2 <- long_GBM_counts %>%
  group_by(tissue) %>%
  mutate(tissue_mean = mean(count), tissue_SD = sd(count))
long_GBM_counts2$tissue <- factor(long_GBM_counts2$tissue, levels = c("adult", "margin", "IDH1_tumor", "fetal"))
long_GBM_counts3 <- long_GBM_counts2[,c(4,5,6)]
long_GBM_counts3 <- unique(long_GBM_counts3)
ggplot(long_GBM_counts3, aes(x = tissue, y = tissue_mean)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 25)) + geom_errorbar(aes(ymin=tissue_mean-tissue_SD, ymax=tissue_mean+tissue_SD), width=.2)

#now, with margin adult genes
long_GBM_counts<-GBM_counts.keep5 %>%
  pivot_longer(!Gene, names_to = "sample", values_to = "count")
new_value <- as.vector(c("IDH1_tumor", "margin", "IDH1_tumor", "margin", "IDH1_tumor", "margin", "IDH1_tumor", rep("adult", times=12), rep("fetal", times=6)))
new_value <- rep(new_value, times=122)
long_GBM_counts$tissue <- as.factor(new_value)
long_GBM_counts2 <- long_GBM_counts %>%
  group_by(tissue) %>%
  mutate(tissue_mean = mean(count), tissue_SE = sd(count)/sqrt(length(count)))
long_GBM_counts2$tissue <- factor(long_GBM_counts2$tissue, levels = c("adult", "margin", "IDH1_tumor", "fetal"))
long_GBM_counts3 <- long_GBM_counts2[,c(4,5,6)]
long_GBM_counts3 <- unique(long_GBM_counts3)
ggplot(long_GBM_counts3, aes(x = tissue, y = tissue_mean)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 25)) + geom_errorbar(aes(ymin=tissue_mean-tissue_SE, ymax=tissue_mean+tissue_SE), width=.2)

#now, with all upregulated margin genes
long_GBM_counts<-GBM_counts.keep7 %>%
  pivot_longer(!Gene, names_to = "sample", values_to = "count")
new_value <- as.vector(c("IDH1_tumor", "margin", "IDH1_tumor", "margin", "IDH1_tumor", "margin", "IDH1_tumor", rep("adult", times=12), rep("fetal", times=6)))
new_value <- rep(new_value, times=146)#747 for margin; #146 for tumor
long_GBM_counts$tissue <- as.factor(new_value)
long_GBM_counts2 <- long_GBM_counts %>%
  group_by(tissue) %>%
  mutate(tissue_mean = mean(count), tissue_SE = sd(count)/sqrt(length(count)))
long_GBM_counts2$tissue <- factor(long_GBM_counts2$tissue, levels = c("adult", "margin", "IDH1_tumor", "fetal"))
long_GBM_counts3 <- long_GBM_counts2[,c(4,5,6)]
long_GBM_counts3 <- unique(long_GBM_counts3)
long_GBM_counts4 <- long_GBM_counts3[-c(3,4),]
long_GBM_counts5 <- long_GBM_counts3[-c(1,2),]
ggplot(long_GBM_counts4, aes(x = tissue, y = tissue_mean)) +
  geom_bar(stat = "identity") + ylim(0,5)+
  theme(axis.text.x = element_text(angle = 25)) + geom_errorbar(aes(ymin=tissue_mean-tissue_SE, ymax=tissue_mean+tissue_SE), width=.2)
#now do the stats, compare expn of margin genes in adult vs fetal and do for tumor genes
long_GBM_counts6<- long_GBM_counts2[long_GBM_counts2$tissue == "adult" | long_GBM_counts2$tissue == "fetal",]
long_GBM_counts7 <- long_GBM_counts6 %>%
  group_by(sample) %>%
  mutate(tissue_mean = mean(count), tissue_SE = sd(count)/sqrt(length(count)))
long_GBM_counts8 <- long_GBM_counts7[,c(2,4,5,6)]
long_GBM_counts8 <- unique(long_GBM_counts8)
adult<- as.vector(long_GBM_counts8[1:12,3])
fetal<- as.vector(long_GBM_counts8[13:18,3])
t.test(adult,fetal,var.equal = T)
#Two Sample t-test for the margin genes
#data:  adult and fetal
#t = 14.92, df = 16, p-value = 8.281e-11
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  1.797659 2.393116
#sample estimates:
#  mean of x mean of y 
#3.375797  1.280409 

#Two Sample t-test for the tumor genes
#data:  adult and fetal
#t = -2.0279, df = 16, p-value = 0.05956
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.49229989  0.01091767
#sample estimates:
#  mean of x mean of y 
#1.131650  1.372341

###QC: show enrichment of astrocyte genes in tumor and margin panned samples
library(readxl)
library(dplyr)
library(RColorBrewer)
cell_markers <- read_xlsx("C:/Users/caitl/Documents/cell_markers.xlsx", sheet = 2)
cell_markers <- as.vector(cell_markers$markers)
df2<- df1[df1$Geneid %in% cell_markers,]
rownames(df2) <- df2$Geneid
df2 <- df2[,-1]
df2 <- df2[cell_markers,]
df3 <- scale(df2)
library(gplots)
heatmap.2(as.matrix(df2), dendrogram= "none",Colv = F,Rowv = F, trace = NULL,scale = "row",col=colorRampPalette(rev(brewer.pal(n=11,name="RdBu")))(100))


GBM_data <- read_excel("C:/Users/caitl/Documents/third rotation/GBM/RNA/GBM feature counts.xlsx", col_names = TRUE)
cell_markers2 <- read_xlsx("C:/Users/caitl/Documents/cell_markers.xlsx", sheet = 3)

GBM_Fastro<- filtered_GBM_data[cell_markers2$`F astro`,]
GBM_Aastro<- filtered_GBM_data[cell_markers2$`A astro`,]
GBM_Olig<- filtered_GBM_data[cell_markers2$Oligo,]
GBM_Neuron<- filtered_GBM_data[cell_markers2$Neuron,]
GBM_Micro<- filtered_GBM_data[cell_markers2$Micro,]
GBM_Endo<- filtered_GBM_data[cell_markers2$Endo,]

mean_Fastro<- colMeans(GBM_Fastro,na.rm=TRUE)
mean_Aastro<- colMeans(GBM_Aastro,na.rm=TRUE)
mean_Olig<- colMeans(GBM_Olig,na.rm=TRUE)
mean_Neuron<- colMeans(GBM_Neuron,na.rm=TRUE)
mean_Micro<- colMeans(GBM_Micro,na.rm=TRUE)
mean_Endo<- colMeans(GBM_Endo,na.rm=TRUE)


##try to plot mean expression of astro, olig, neuron etc genes in tumor and margin
countData<- read_excel("C:/Users/caitl/Documents/third rotation/GBM/RNA/GBM feature counts.xlsx", col_names = TRUE)
#vst normalization
df2<- countData[!duplicated(countData$Geneid), ]
df3 <- as.data.frame(df2)
row.names(df3) <- df3$Geneid
df3 <- df3[,-c(1:5,40)]
#Filtering count data:
#Step1: filter out the genes for which there are no reads, or there are inconsistent reads across replicate samples, or there are low reads. Chose a normalization technique: CPM (counts per million), for example (GBM_CPM). We filter using CPM values rather than counts because they account for differences in sequencing depth between samples. 
GBM_CPM <- edgeR::cpm(df3)
head(GBM_CPM)
#Step2: Our filtering rule is to keep transcripts that have CPM > X in at least two samples. X is the cutoff value (in CPM) that will consistently yield a minimum transcript count. We need to first determine the scientifically relevant minimum transcript count. As a general rule, a good threshold can be chosen for a CPM value that corresponds to a count of 10 (for 10-20M reads/sample depth). This will vary based on the depth of sequence; for low depth sequence <1M, try a count of 1. CPM = counts per million, or how many counts would I get for a gene if the sample had a library size of 1M. So for a library size of 1M, 1 count = 1 CPM. For a library size of 10M, 10 counts = 1 CPM. For a library size of 20M, 10 counts = 0.5 CPM.
#Step 3: Next, impose the threshold (0.5 is an example value, replace with determined X value).
GBM_thresh <- GBM_CPM > 0.5
#Step 4: Identify/subset the genes for which CPM > 0.5 in at least two samples. Use DGElist to convert conunts.keep df back to an oject (GBM_count) for ease of analysis moving forward.
GBM_keep <- rowSums(GBM_thresh) >= 2
GBM_counts.keep <- df3[GBM_keep,]
GBM_counts.keep2<- cpm(GBM_counts.keep, log = T)

cell_markers2 <- read_xlsx("C:/Users/caitl/Documents/cell_markers.xlsx", sheet = 4)
GBM_counts.keep3 <- as.data.frame(GBM_counts.keep2)
GBM_counts.keep3<- GBM_counts.keep3 %>% filter(rownames(GBM_counts.keep3) %in% cell_markers2$Gene)#filter normalized count mtx to only include master list genes
GBM_counts.keep3$Gene <- rownames(GBM_counts.keep3)

#add tissue annotation
long_GBM_counts<-GBM_counts.keep3 %>%
  pivot_longer(!Gene, names_to = "sample", values_to = "count")
new_value <- as.vector(c(rep("tumor", times=18), rep("margin", times=16)))
new_value <- rep(new_value, times=285)
long_GBM_counts$tissue <- as.factor(new_value)
long_GBM_counts<- merge(long_GBM_counts, cell_markers2, by="Gene")
long_GBM_counts2 <- long_GBM_counts %>%
  group_by(Type,tissue) %>%
  mutate(tissue_mean = mean(count), tissue_SD = sd(count))
long_GBM_counts2$tissue <- factor(long_GBM_counts2$tissue, levels = c("tumor", "margin"))
long_GBM_counts3 <- long_GBM_counts2[,c(4,5,6,7)]
long_GBM_counts3 <- unique(long_GBM_counts3)
long_GBM_counts4 <- as.data.frame(long_GBM_counts3)
long_GBM_counts4$concat <- paste(long_GBM_counts4$tissue, long_GBM_counts4$Type)
ggplot(long_GBM_counts4, aes(x = concat, y = tissue_mean)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 25)) + geom_errorbar(aes(ymin=tissue_mean-tissue_SD, ymax=tissue_mean+tissue_SD), width=.2)