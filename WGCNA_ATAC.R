library(tidyverse)     # tidyverse will pull in ggplot2, readr, other useful libraries
library(magrittr)      # provides the %>% operator
library(WGCNA)
library(dplyr)
library(ggplot2)
library(readxl)

#Make count matrix into the proper format
count_df <- as.data.frame(dsa@counts$.peaks.cons)
cts<- cbind(all_peaks2$name, count_df)
colnames(cts)[1]<- paste("Geneid")
colnames(cts)[7]<- paste("C4_D150")
#did the above steps in my instance and exported the cts object
cts <- readRDS("C:/Users/caitl/Documents/third rotation/Organoid maturation timeline/HiSeq/ATAC-seq/ChrAccR/org_count_df.rds")
#Make a pivot table separated by gene name 
cts<-readRDS("/Users/stevensloan/Downloads/cts.rds") #ignore this step
col_sel = names(cts)[-1]
mdata <- cts %>%
  tidyr::pivot_longer(., col = all_of(col_sel)) %>%   # The dot is the the input data, magrittr tutorial
  mutate(group = gsub(".*(D)","\\1", name) %>%  gsub("_hepa","", .)   # Get the shorter treatment names
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
    Type = gsub(".*(D)","\\1", Sample) %>% gsub("_hepa","", .)
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
q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= changed to 95 quantile to reduce dataset
expr_normalized <- wpn_vsd[ rv_wpn > q95_wpn, ]
q95_wpn <- quantile( rowVars(wpn_vsd), .90)  # <= changed to 90 quantile to reduce dataset
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
#Now let’s transpose the data and prepare the dataset for WGCNA.
################################################################
input_mat = t(expr_normalized)

input_mat[1:5,1:10]

#We can see now that the rows = treatments and columns = gene probes. We’re ready to start WGCNA. 
#A correlation network will be a complete network (all genes are connected to all other genes). 
#Ergo we will need to pick a threshhold value (if correlation is below threshold, remove the edge).

#To do that, WGCNA will try a range of soft thresholds and create a diagnostic plot. This step will 
#take several minutes so feel free to run and get coffee.

allowWGCNAThreads()          # allow multi-threading (optional)

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 40, by = 2))

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



picked_power = 30
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 4,
                          pamRespectsDendro = F,
                          #detectCutHeight = 0.75,
                          minModuleSize = 100,
                          maxBlockSize = 4000,#change to 1000
                          detectCutHeight = 0.995,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.1, #default was 0.25
                          
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

###############################################################
### STEVEN ADDED THIS PART TO RE-SPLIT THE ORIGINAL NETWORK ###
###############################################################

netwk2<-recutBlockwiseTrees(input_mat,netwk$goodSamples,netwk$goodGenes,
                    netwk$blocks,netwk$TOMFiles,netwk$dendrograms,
                    corType = "pearson",
                    networkType = "unsigned",
                    deepSplit = 3,
                    detectCutHeight = 0.95, minModuleSize = 700, #minModuleSize = 100 for top 5%
                    maxCoreScatter = NULL, minGap = NULL,
                    maxAbsCoreScatter = NULL, minAbsGap = NULL,
                    pamStage = TRUE, pamRespectsDendro = TRUE,
                    #minCoreKME = 0.5, minCoreKMESize = minModuleSize/3,
                    #minKMEtoStay = 0.3,
                    reassignThreshold = 1e-6,
                    mergeCutHeight = 0.15, impute = TRUE, #mergecutheight= 0.1 for top 5%
                    trapErrors = FALSE, numericLabels = TRUE,
                    verbose = 3)

mergedColors = labels2colors(netwk2$colors)

# Check to see now how many modules there are
netwk2$colors[netwk2$blockGenes[[1]]]
table(netwk2$colors)

#modules: 

###############################################################
###############################################################
###############################################################

#Relate Module (cluster) Assignments to Treatment Groups
module_df <- data.frame(
  gene_id = names(netwk2$colors),
  colors = labels2colors(netwk2$colors)
)
write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")

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
    low = "blue",
    high = "red",
    mid = "white",
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
modules_of_interest = c("blue", "turquoise", "brown", "grey", "yellow")
#grey module is late module with 1,136 peaks
#turquoise module is mid/late module with 15,159 peaks
#brown module is middle module with 821 peaks
#blue module is early module with 13,386 peaks

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id
subexpr = expr_normalized[submod$gene_id,]

colnames(subexpr) <- factor(colnames(subexpr), levels = colnames(subexpr))

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df$name <- factor(submod_df$name, levels = submod_df$name[1:20])

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

blue_mod<- module_df[module_df$colors == "blue",] #9,798 "late" for top 90
turq_mod<- module_df[module_df$colors == "turquoise",] #12,830 "early" for top 90
brown_mod<- module_df[module_df$colors == "brown",] #5,374 "middle/late" for top 90
grey_mod<- module_df[module_df$colors == "grey",] #1,679 "early/middle" for top 90
yellow_mod<- module_df[module_df$colors == "yellow",] #821 "middle" for top 90

save(blue_mod, turq_mod, brown_mod, grey_mod, yellow_mod, file = "ATAC_WGCNA_mods.RData")
