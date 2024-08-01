library(tidyverse)     # tidyverse will pull in ggplot2, readr, other useful libraries
library(magrittr)      # provides the %>% operator
library(WGCNA)
library(dplyr)
library(ggplot2)
library(readxl)

#Load Count Data
setwd("~/third rotation/Organoid maturation timeline/HiSeq/RNA-seq")
cts <- read_xlsx("wgcna input.xlsx", col_names = TRUE) 
cts<- cts[!duplicated(cts$Geneid), ]

#remove rows for which there is no approved gene symbol (df2), make gene name the row name (df2), remove the gene name column (df3).
cts <- cts[!is.na(cts$`C4 d550`), ]

#Make a pivot table separated by gene name 
col_sel = names(cts)[-1]
mdata <- cts %>%
  tidyr::pivot_longer(., col = all_of(col_sel)) %>%   # The dot is the the input data, magrittr tutorial
  mutate(group = gsub(".*(d)","\\1", name) %>%  gsub("_HEPA","", .)   # Get the shorter treatment names
  )

##########################################################
##DESeq
#Normalize Counts with DESeq
##########################################################
library(DESeq2)
#Prepare DESeq input, which is expecting a matrix of integers.
de_input = as.matrix(cts[,-1])
row.names(de_input) = cts$Geneid
de_input[1:5,1:10]

meta_df <- data.frame( Sample = names(cts[-1])) %>%
  mutate(
    Type = gsub(".*(d)","\\1", Sample) %>% gsub("_HEPA","", .)
  )

dds <- DESeqDataSetFromMatrix(round(de_input),
                              meta_df,
                              design = ~Type)
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds)

wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)
summary(rv_wpn)

q75_wpn <- quantile( rowVars(wpn_vsd), .75)  # <= original
expr_normalized <- wpn_vsd[ rv_wpn > q75_wpn, ]
q85_wpn <- quantile( rowVars(wpn_vsd), .85)  # 3973 genes *ended up using this cutoff
expr_normalized <- wpn_vsd[ rv_wpn > q85_wpn, ]
q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= changed to 95 quantile to reduce dataset
expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]

expr_normalized[1:5,1:10]

expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(
    Gene_id = row.names(expr_normalized)
  ) %>%
  pivot_longer(-Gene_id)

expr_normalized_df %>% ggplot(., aes(x = name, y = value)) +
  geom_violin() +
  geom_point() +
  theme_bw() +
  theme(
    axis.text.x = element_text( angle = 90)
  ) +
  ylim(0, NA) +
  labs(
    title = "Normalized and 95 quantile Expression",
    x = "treatment",
    y = "normalized expression"
  )

#####################################
#Now let's transpose the data and prepare the dataset for WGCNA.
################################################################
input_mat = t(expr_normalized)

input_mat[1:5,1:10]

#We can see now that the rows = treatments and columns = gene probes. We're ready to start WGCNA. 
#A correlation network will be a complete network (all genes are connected to all other genes). 
#Ergo we will need to pick a threshhold value (if correlation is below threshold, remove the edge).

#To do that, WGCNA will try a range of soft thresholds and create a diagnostic plot. This step will 
#take several minutes so feel free to run and get coffee.

allowWGCNAThreads()          # allow multi-threading (optional)

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")
#Pick a soft threshold power near the curve of the plot. For me, maybe 10, 12, or 14.
#Now we can create the network using the blockwiseModules command. The blockwiseModule 
#may take a while to run, since it is constructing the TOM (topological overlap matrix) 
#and several other steps. 
picked_power = 12
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 200, #updated from 60 on 5/20/22
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.325, #default was 0.25, 0.3 for more modules, 0.35 for fewer; updated from 0.35 on 5/20/22
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
cor <- temp_cor     # Return cor function to original namespace

##looking at the modules
# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )


netwk$colors[netwk$blockGenes[[1]]]
table(netwk$colors)





#Relate Module (cluster) Assignments to Treatment Groups
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)
write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")

module_df$colors <- as.character(module_df$colors)

module_df$colors[module_df$colors == "turquoise"] <- "Early"
module_df$colors[module_df$colors == "yellow"] <- "Early_Middle"
module_df$colors[module_df$colors == "brown"] <- "Middle"
module_df$colors[module_df$colors == "green"] <- "Late"
module_df$colors[module_df$colors == "blue"] <- "Middle_Late"

module_df2<-module_df[!(module_df$colors=="grey"),]

write.csv(module_df2, "org_RNAseq_gene_modules.csv", col.names = TRUE)

#we need to figure out which modules are associated with each 
#trait/treatment group. WGCNA will calcuate an Eigangene (hypothetical central gene) 
#for each module, so it easier to determine if modules are associated with different treatments.

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)
MEs0$treatment <- factor(MEs0$treatment, levels = MEs0$treatment)


# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "darkblue",
    high = "darkred",
    mid = "white",
    space = "Lab",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")


mME %>% ggplot(., aes(x=treatment, y=name)) +
  geom_tile(aes(fill = value),colour = "white") +
  theme_bw() +
  scale_fill_distiller(palette = "RdBu", direction = -1) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")

mME %>% ggplot(., aes(x=treatment, y=name)) +
  geom_tile(aes(fill = value),colour = "white") +
  theme_bw() +
  scale_fill_distiller(palette = "BuPu", direction = 1) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")


####pull out genes for further analysis
# pick out a few modules of interest here
modules_of_interest = c("blue", "turquoise", "brown", "green", "yellow")#most up to date

modules_of_interest = c("blue", "turquoise", "brown", "red","green", "yellow", "black")
modules_of_interest = c("blue", "turquoise", "brown", "yellow") #updated: 5/9/2022

# Pull out list of genes in that module
submod_early = module_df %>%
  subset(colors %in% modules_of_interest)
submod_middle = module_df %>%
  subset(colors %in% modules_of_interest)
submod_late = module_df %>%
  subset(colors %in% modules_of_interest)
save(submod_early, submod_middle, submod_late, file = "binned_submods.RData")

submod = module_df %>%
  subset(colors %in% modules_of_interest)


row.names(module_df) = module_df$gene_id
subexpr = expr_normalized[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df$name <- factor(submod_df$name, levels = unique(submod_df$name))

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")


######################################################
#Look at various modules in GBM data
#####################################################
blue_mod = c("blue")
turquoise_mod = c("turquoise")
brown_mod = c("brown")
yellow_mod = c("yellow") #added 5/9/2022
green_mod = c("green")

# Pull out list of genes in that module
blu_mod<- module_df[module_df$colors == "blue",] #Middle/Late = 926 genes
t_mod<- module_df[module_df$colors == "turquoise",] #Early = 1134 genes
br_mod<- module_df[module_df$colors == "brown",] #Middle = 861 genes
y_mod<- module_df[module_df$colors == "yellow",] #Early/Middle = 656 genescai
g_mod<- module_df[module_df$colors == "green",] #Late = 367 genes

###########preparing GBM data
library(Glimma)
library(tidyverse)
library(factoextra)
library(viridis)
library(cluster)
library(clustertend)
library(edgeR)
library(limma)
library(gplots)
library(RColorBrewer)
setwd("~/third rotation/GBM/RNA")
feature_counts <- read_csv("feature_counts_raw.csv")
#removing duplicates genes
df1<- feature_counts[!duplicated(feature_counts$Geneid), ]
#remove rows for which there is no approved gene symbol (df2), make gene name the row name (df2), remove the gene name column (df3).
df2 <- df1[!is.na(df1$Geneid), ]
df2 <- as.data.frame(df2)
row.names(df2) <- df2$Geneid
df3 <- df2[,-1]
#Filtering count data:
#Step1: filter out the genes for which there are no reads, or there are inconsistent reads across replicate samples, or there are low reads. Chose a normalization technique: CPM (counts per million), for example (GBM_CPM). We filter using CPM values rather than counts because they account for differences in sequencing depth between samples.
GBM_CPM <- edgeR::cpm(df3)
head(GBM_CPM)
GBM_thresh <- GBM_CPM > 0.5
#Step 4: Identify/subset the genes for which CPM > 0.5 in at least two samples. Use DGElist to convert conunts.keep df back to an oject (GBM_count) for ease of analysis moving forward.
GBM_keep <- rowSums(GBM_thresh) >= 2
GBM_counts.keep <- df3[GBM_keep,]
GBM_counts.keep<- cpm(GBM_counts.keep)
GBM_counts <- edgeR::DGEList(GBM_counts.keep)
sampleinfo <- read.delim2("meta_all.txt", header = TRUE, sep = "\t")
GBM_counts.keep2<- GBM_counts.keep[,-c(1:28, 40:45, 48, 51, 52, 55, 56, 61, 62)]
GBM_counts <- edgeR::DGEList(GBM_counts.keep2)
sampleinfo_tissue<- sampleinfo[-c(1:28, 40:45, 48, 51, 52, 55, 56, 61, 62),]
sampleinfo_tissue$patient <- c("A","A","B","B","C","C","D","D","E","F","F","G","G","H","H","I","I","J","J","K","K")
###########try new DE method (12/13/21)
Patient <- factor(paste(sampleinfo_tissue$patient,sep="."))
Tissue <- factor(paste(sampleinfo_tissue$age,sep="."))
design <- model.matrix(~ 0+Patient+Tissue)
rownames(design) <- colnames(GBM_counts)
y <- estimateDisp(GBM_counts, design, robust=TRUE)
y$common.dispersion
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
library(org.Hs.eg.db)
ann <- select(org.Hs.eg.db,keys=rownames(lrt),keytype = "SYMBOL" ,columns=c("ENTREZID","SYMBOL","GENENAME"))
head(ann)
ann2 <- ann %>% distinct(SYMBOL, .keep_all = TRUE)
table(ann2$ENTREZID==rownames(lrt))
lrt$genes <- ann2
lrt.res<- topTags(lrt, n="Inf", sort.by = "PValue")
signifgenes <- lrt.res$table %>% filter(FDR < 0.05)
signifgenes_down <- signifgenes %>% filter(logFC < -2)
signifgenes_up <- signifgenes %>% filter(logFC > 2)

#genes that are down and in various modules
down_blue <- as.vector(intersect(signifgenes_down$SYMBOL, blu_mod$gene_id)) #63; 5/20/22: 67
down_turq <- as.vector(intersect(signifgenes_down$SYMBOL, t_mod$gene_id)) #90; 5/20/22: 89
down_brown <- as.vector(intersect(signifgenes_down$SYMBOL, br_mod$gene_id)) #23; 5/20/22: 35
down_red <- as.vector(intersect(signifgenes_down$SYMBOL, r_mod$gene_id)) #12
down_green <- as.vector(intersect(signifgenes_down$SYMBOL, g_mod$gene_id)) #22; 5/20/22: 22
down_black <- as.vector(intersect(signifgenes_down$SYMBOL, bla_mod$gene_id)) #0
down_yellow <- as.vector(intersect(signifgenes_down$SYMBOL, y_mod$gene_id)) #27; 5/20/22: 24
up_blue <- as.vector(intersect(signifgenes_up$SYMBOL, blu_mod$gene_id)) #16; 5/20/22: 19
up_turq <- as.vector(intersect(signifgenes_up$SYMBOL, t_mod$gene_id)) #96; 5/20/22: 96
up_brown <- as.vector(intersect(signifgenes_up$SYMBOL, br_mod$gene_id)) #10; 5/20/22: 8
up_red <- as.vector(intersect(signifgenes_up$SYMBOL, r_mod$gene_id)) #2
up_green <- as.vector(intersect(signifgenes_up$SYMBOL, g_mod$gene_id)) #6; 5/20/22: 27
up_black <- as.vector(intersect(signifgenes_up$SYMBOL, bla_mod$gene_id)) #1
up_yellow <- as.vector(intersect(signifgenes_up$SYMBOL, y_mod$gene_id)) #30; 5/20/22: 10


#Make volcano plots for following modules: blue, turquoise, yellow, brown, green

#blue
with(lrt.res$table, plot(logFC, -log10(PValue), pch=20, main="DE genes Tumor vs Margin Middle/Late Module", xlim=c(-12,12), ylim=c(0,65)))
with(subset(lrt.res$table, FDR<.05 & logFC < -2), points(logFC, -log10(PValue), pch=20, col="#C1C7C9"))
with(subset(lrt.res$table, FDR<.05 & logFC>2), points(logFC, -log10(PValue), pch=20, col="#C1C7C9"))
with(subset(lrt.res$table, SYMBOL %in% down_blue | SYMBOL %in% up_blue), points(logFC, -log10(PValue), pch=20, col="blue"))

#turquoise
with(lrt.res$table, plot(logFC, -log10(PValue), pch=20, main="DE genes Tumor vs Margin Early Module", xlim=c(-12,12), ylim=c(0,65)))
with(subset(lrt.res$table, FDR<.05 & logFC < -2), points(logFC, -log10(PValue), pch=20, col="#C1C7C9"))
with(subset(lrt.res$table, FDR<.05 & logFC>2), points(logFC, -log10(PValue), pch=20, col="#C1C7C9"))
with(subset(lrt.res$table, SYMBOL %in% down_turq | SYMBOL %in% up_turq), points(logFC, -log10(PValue), pch=20, col="#30D5C8"))

#Brown
with(lrt.res$table, plot(logFC, -log10(PValue), pch=20, main="DE genes Tumor vs Margin Middle Module", xlim=c(-12,12), ylim=c(0,65)))
with(subset(lrt.res$table, FDR<.05 & logFC < -2), points(logFC, -log10(PValue), pch=20, col="#C1C7C9"))
with(subset(lrt.res$table, FDR<.05 & logFC>2), points(logFC, -log10(PValue), pch=20, col="#C1C7C9"))
with(subset(lrt.res$table, SYMBOL %in% down_brown | SYMBOL %in% up_brown), points(logFC, -log10(PValue), pch=20, col="brown"))

#red
with(lrt.res$table, plot(logFC, -log10(PValue), pch=20, main="DE genes Tumor vs Margin Middle Module", xlim=c(-12,12), ylim=c(0,65)))
with(subset(lrt.res$table, FDR<.05 & logFC < -2), points(logFC, -log10(PValue), pch=20, col="#C1C7C9"))
with(subset(lrt.res$table, FDR<.05 & logFC>2), points(logFC, -log10(PValue), pch=20, col="#C1C7C9"))
with(subset(lrt.res$table, SYMBOL %in% down_red | SYMBOL %in% up_red), points(logFC, -log10(PValue), pch=20, col="red"))

#Green
with(lrt.res$table, plot(logFC, -log10(PValue), pch=20, main="DE genes Tumor vs Margin Late Module", xlim=c(-12,12), ylim=c(0,65)))
with(subset(lrt.res$table, FDR<.05 & logFC < -2), points(logFC, -log10(PValue), pch=20, col="#C1C7C9"))
with(subset(lrt.res$table, FDR<.05 & logFC>2), points(logFC, -log10(PValue), pch=20, col="#C1C7C9"))
with(subset(lrt.res$table, SYMBOL %in% down_green | SYMBOL %in% up_green), points(logFC, -log10(PValue), pch=20, col="green"))

#Black
with(lrt.res$table, plot(logFC, -log10(PValue), pch=20, main="DE genes Tumor vs Margin Middle 2 Module", xlim=c(-12,12), ylim=c(0,65)))
with(subset(lrt.res$table, FDR<.05 & logFC < -2), points(logFC, -log10(PValue), pch=20, col="#C1C7C9"))
with(subset(lrt.res$table, FDR<.05 & logFC>2), points(logFC, -log10(PValue), pch=20, col="#C1C7C9"))
with(subset(lrt.res$table, SYMBOL %in% down_black | SYMBOL %in% up_black), points(logFC, -log10(PValue), pch=20, col="pink"))

#Yellow
with(lrt.res$table, plot(logFC, -log10(PValue), pch=20, main="DE genes Tumor vs Margin Early/Middle Module", xlim=c(-12,12), ylim=c(0,65)))
with(subset(lrt.res$table, FDR<.05 & logFC < -2), points(logFC, -log10(PValue), pch=20, col="#C1C7C9"))
with(subset(lrt.res$table, FDR<.05 & logFC>2), points(logFC, -log10(PValue), pch=20, col="#C1C7C9"))
with(subset(lrt.res$table, SYMBOL %in% down_yellow | SYMBOL %in% up_yellow), points(logFC, -log10(PValue), pch=20, col="#FFD700"))



####################correlate GBM and organoids

#load data
feature_counts <- read_xlsx("~/third rotation/GBM/RNA/GBM feature counts.xlsx")
#removing duplicates genes 
df1<- feature_counts[!duplicated(feature_counts$Geneid), ]
#remove rows for which there is no approved gene symbol (df2), make gene name the row name (df2), remove the gene name column (df3).
df2 <- df1[!is.na(df1$Geneid), ]
df2 <- as.data.frame(df2)
row.names(df2) <- df2$Geneid
df3 <- df2[,-c(1:2)]

gbm_green<- df3 %>% filter(rownames(df3) %in% g_mod$gene_id)
gbm_red<- df3 %>% filter(rownames(df3) %in% y_mod$gene_id)
gbm_brown<- df3 %>% filter(rownames(df3) %in% br_mod$gene_id)
gbm_turq<- df3 %>% filter(rownames(df3) %in% t_mod$gene_id)
gbm_blue<- df3 %>% filter(rownames(df3) %in% blu_mod$gene_id)
gbm_filtered<- df3 %>% filter(rownames(df3) %in% submod$gene_id)
submod2<- submod %>% filter(submod$gene_id %in% rownames(gbm_filtered))
# split by a vector specifying rowgroups
library(ComplexHeatmap)
pseudoNormCounts <- cpm(gbm_filtered, prior.count = 1)
scaled_mat = t(scale(t(as.matrix(pseudoNormCounts))))#scale by row
scaled_mat_2 <- scaled_mat[!rowSums(!is.finite(scaled_mat)),]
submod2<- submod2 %>% filter(submod2$gene_id %in% rownames(scaled_mat_2))
#split heatmap by module
Heatmap(scaled_mat_2, split = submod2$colors, cluster_rows = TRUE,cluster_columns = FALSE, 
        row_names_gp = gpar(fontsize = 7))
#don't split, cluster rows, and color label each row by module
Heatmap(scaled_mat_2, cluster_rows = TRUE,cluster_columns = FALSE, 
        row_names_gp = gpar(fontsize = 7)) +
  Heatmap(submod2$colors, name = "module", width = unit(5, "mm"))


#remove IDH1 mutants and recurrants and re-order the columns
updated_gbm_count_df <- gbm_count_df2[,-c(12,14,16,18,21,24,30)] #removing weird tumors, kept the paired margin samples for now
updated_gbm_count_df <- updated_gbm_count_df[,c(2,3,5,6,8,9,10,14,15,18,20,22,1,4,7,11,12,13,16,17,19,21,23)]

#re-order columns by decreasing STEM score (from RNA-seq data) with all samples
updated_gbm_count_df <- gbm_filtered[,c(5,8,4,14,2,12,10,6,3,7,17,1,16,18,9,15,13,11,22,33,23,30,25,28,20,21,32,19,34,27,24,26,29,31)]
#the above, but without IDH mutant and recurrent
updated_gbm_count_df <- gbm_filtered[,c(5,4,14,2,12,10,6,3,17,1,16,22,33,23,30,25,28,20,21,32,19,34,27,29,31)]

##########correlate GBM and organoid samples using only module genes
cts2 <- as.data.frame(cts)
rownames(cts2) <- cts2$Geneid
cts2<- cts2[,-1]#?
cts2_filtered<- cts2 %>% filter(rownames(cts2) %in% rownames(gbm_filtered))#?
df3 <- df2[,-c(2:5,40)]
df3 <- df2[,-c(2,12:14,29:31)]
df4 <- df3[,-c(12,14,16,19,27,29,32)]
gbm_filtered<- df3 %>% filter(df3$Geneid %in% submod$gene_id)

merged_data <- merge(gbm_filtered, cts2, by.x = "Geneid", by.y = "Geneid")
rownames(merged_data) <- merged_data$Geneid
merged_data <- merged_data[,-1]
merged_matrix <- as.matrix(merged_data)

#scale data since it was collected differently
scaled_merged <- scale(merged_matrix) # scale and center columns

#spearman correlation that compares all samples
library(Hmisc)
cor_result_merged <- round(cor(scaled_merged, method="spearman"),2)

##get rid of unnecessary columns and make that the cor_matrix
library(gplots)
heatmap.2(as.matrix(cor_result_merged[-c(33:55),-c(1:32)]), dendrogram= "none",Colv = F,Rowv = F, trace = NULL,scale = "column",col=bluered(50))
heatmap.2(as.matrix(cor_result_merged[-c(33:55),-c(1:32)]), dendrogram= "none",Colv = F,Rowv = F, trace = NULL,scale = "column",col=colorRampPalette(rev(brewer.pal(n=11,name="RdBu")))(100))


########scoring GBM samples for each module
gbm_green<- df3 %>% filter(rownames(df3) %in% g_mod$gene_id)
gbm_red<- df3 %>% filter(rownames(df3) %in% r_mod$gene_id)
gbm_brown<- df3 %>% filter(rownames(df3) %in% br_mod$gene_id)
gbm_turq<- df3 %>% filter(rownames(df3) %in% t_mod$gene_id)
gbm_blue<- df3 %>% filter(rownames(df3) %in% bl_mod$gene_id)

colSum_green <- colSums(gbm_green)
colSum_red <- colSums(gbm_red)
colSum_brown <- colSums(gbm_brown)
colSum_turq <- colSums(gbm_turq)
colSum_blue <- colSums(gbm_blue)

colSum_combined <- rbind(colSum_red, colSum_turq, colSum_green, colSum_brown, colSum_blue)
rownames(colSum_combined)<- c("red", "turq", "green", "brown", "blue")
scaled_colSum_combined <- scale(colSum_combined)
heatmap.2(as.matrix(scaled_colSum_combined), scale = "row", Rowv = F)


#########################################################################

#finding the top X genes in each module. adapted from the "chooseTopHubInEachModule" function <https://rdrr.io/cran/WGCNA/src/R/dendrogramAdjustmentFunctions.R>
#for every module, (1) calculate the network adjacency for each gene (see: <https://rdrr.io/cran/WGCNA/man/adjacency.html>),(2) calculate the rowSums of the adjacency, 
# (3) find the max 10 genes 
#steps:

# (1) get gene expression matrix where samples are rows and genes are columns
load("C:/Users/caitl/Documents/third rotation/Organoid maturation timeline/HiSeq/RNA-seq/WGCNA/WGCNA_2.RData")
# input_mat is this format, but is normalized; hopefully it is correct to use the normalized data instead of raw
input_df <- as.data.frame(t(input_mat))
# (2) group by genes in the same module
turq_input<- input_df %>% dplyr::filter(rownames(input_df) %in% t_mod$gene_id)
turq_input <- as.matrix(t(turq_input))

yellow_input<- input_df %>% dplyr::filter(rownames(input_df) %in% y_mod$gene_id)
yellow_input <- as.matrix(t(yellow_input))

brown_input<- input_df %>% dplyr::filter(rownames(input_df) %in% br_mod$gene_id)
brown_input <- as.matrix(t(brown_input))

blue_input<- input_df %>% dplyr::filter(rownames(input_df) %in% blu_mod$gene_id)
blue_input <- as.matrix(t(blue_input))

green_input<- input_df %>% dplyr::filter(rownames(input_df) %in% g_mod$gene_id)
green_input <- as.matrix(t(green_input))

# (3) calculate adj and figure out necessary power value (generates something like a correlation matrix of adj values w/ genes as columns and rows)
turq_adj = adjacency(turq_input,power = 12, corOptions = list(use = "p"))
yellow_adj = adjacency(yellow_input,power = 12, corOptions = list(use = "p"))
brown_adj = adjacency(brown_input,power = 12, corOptions = list(use = "p"))
blue_adj = adjacency(blue_input,power = 12, corOptions = list(use = "p"))
green_adj = adjacency(green_input,power = 12, corOptions = list(use = "p"))

# (4) find indices of genes with the highest (or other #) adj values using: tail(order(rowSums(adj)), 10)
turq_adj_top10<- tail(order(rowSums(turq_adj)), 10)
turq_adj_top25<- tail(order(rowSums(turq_adj)), 25)
turq_adj_top50<- tail(order(rowSums(turq_adj)), 50)
turq_adj_top100<- tail(order(rowSums(turq_adj)), 100)

yellow_adj_top10<- tail(order(rowSums(yellow_adj)), 10)
yellow_adj_top25<- tail(order(rowSums(yellow_adj)), 25)
yellow_adj_top50<- tail(order(rowSums(yellow_adj)), 50)
yellow_adj_top100<- tail(order(rowSums(yellow_adj)), 100)

brown_adj_top10<- tail(order(rowSums(brown_adj)), 10)
brown_adj_top25<- tail(order(rowSums(brown_adj)), 26)
brown_adj_top50<- tail(order(rowSums(brown_adj)), 51)
brown_adj_top100<- tail(order(rowSums(brown_adj)), 101)

blue_adj_top10<- tail(order(rowSums(blue_adj)), 10)
blue_adj_top25<- tail(order(rowSums(blue_adj)), 25)
blue_adj_top50<- tail(order(rowSums(blue_adj)), 51)
blue_adj_top100<- tail(order(rowSums(blue_adj)), 101)

green_adj_top10<- tail(order(rowSums(green_adj)), 10)
green_adj_top25<- tail(order(rowSums(green_adj)), 25)
green_adj_top50<- tail(order(rowSums(green_adj)), 50)
green_adj_top100<- tail(order(rowSums(green_adj)), 100)

# (5) get gene names associated with those indices
turq_adj_top10<- colnames(turq_adj)[turq_adj_top10]
# [1] "RSPO1"    "ARL4A"    "CBLN1"    "GALNT5"   "CRABP2"   "C14orf39" "TPBG"     "UNC5C"    "WNT5A"    "WNT10B"
turq_adj_top25<- colnames(turq_adj)[turq_adj_top25]
turq_adj_top50<- colnames(turq_adj)[turq_adj_top50]
turq_adj_top100<- colnames(turq_adj)[turq_adj_top100]
write.csv(turq_adj_top10, file = "early_top10.csv")
write.csv(turq_adj_top25, file = "early_top25.csv")
write.csv(turq_adj_top50, file = "early_top50.csv")
write.csv(turq_adj_top100, file = "early_top100.csv")

yellow_adj_top10<- colnames(yellow_adj)[yellow_adj_top10]
yellow_adj_top25<- colnames(yellow_adj)[yellow_adj_top25]
yellow_adj_top50<- colnames(yellow_adj)[yellow_adj_top50]
yellow_adj_top100<- colnames(yellow_adj)[yellow_adj_top100]
write.csv(yellow_adj_top10, file = "earlyMiddle_top10.csv")
write.csv(yellow_adj_top25, file = "earlyMiddle_top25.csv")
write.csv(yellow_adj_top50, file = "earlyMiddle_top50.csv")
write.csv(yellow_adj_top100, file = "earlyMiddle_top100.csv")

brown_adj_top10<- colnames(brown_adj)[brown_adj_top10]
brown_adj_top25<- colnames(brown_adj)[brown_adj_top25]
brown_adj_top50<- colnames(brown_adj)[brown_adj_top50]
brown_adj_top100<- colnames(brown_adj)[brown_adj_top100]
write.csv(brown_adj_top10, file = "Middle_top10.csv")
write.csv(brown_adj_top25, file = "Middle_top25.csv")
write.csv(brown_adj_top50, file = "Middle_top50.csv")
write.csv(brown_adj_top100, file = "Middle_top100.csv")

blue_adj_top10<- colnames(blue_adj)[blue_adj_top10]
blue_adj_top25<- colnames(blue_adj)[blue_adj_top25]
blue_adj_top50<- colnames(blue_adj)[blue_adj_top50]
blue_adj_top100<- colnames(blue_adj)[blue_adj_top100]
write.csv(blue_adj_top10, file = "MiddleLate_top10.csv")
write.csv(blue_adj_top25, file = "MiddleLate_top25.csv")
write.csv(blue_adj_top50, file = "MiddleLate_top50.csv")
write.csv(blue_adj_top100, file = "MiddleLate_top100.csv")

green_adj_top10<- colnames(green_adj)[green_adj_top10]
green_adj_top25<- colnames(green_adj)[green_adj_top25]
green_adj_top50<- colnames(green_adj)[green_adj_top50]
green_adj_top100<- colnames(green_adj)[green_adj_top100]
write.csv(green_adj_top10, file = "Late_top10.csv")
write.csv(green_adj_top25, file = "Late_top25.csv")
write.csv(green_adj_top50, file = "Late_top50.csv")
write.csv(green_adj_top100, file = "Late_top100.csv")

stripchart(v$E["CD45",]~DAY,vertical=TRUE,las=2,cex.axis=0.8,pch=16,cex=1.3,col=nice.col,method="jitter",ylab="Normalised log2 expression",main="CD45")