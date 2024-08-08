
library(Seurat)
setwd("/home/sloanlab")
```



```{r}
#load data
data_GSM6720852 <- Read10X(data.dir = '/home/ubuntu/data/CP/GSM6720852')
data_GSM6720854 <- Read10X(data.dir = '/home/ubuntu/data/CP/GSM6720854')
data_GSM6720856 <- Read10X(data.dir = '/home/ubuntu/data/CP/GSM6720856')
data_GSM6720858 <- Read10X(data.dir = '/home/ubuntu/data/CP/GSM6720858')
data_GSM6720860 <- Read10X(data.dir = '/home/ubuntu/data/CP/GSM6720860')
data_GSM6720862 <- Read10X(data.dir = '/home/ubuntu/data/CP/GSM6720862')
data_GSM6720864 <- Read10X(data.dir = '/home/ubuntu/data/CP/GSM6720864')
data_GSM6720866 <- Read10X(data.dir = '/home/ubuntu/data/CP/GSM6720866')
data_GSM6720868 <- Read10X(data.dir = '/home/ubuntu/data/CP/GSM6720868')
data_GSM6720870 <- Read10X(data.dir = '/home/ubuntu/data/CP/GSM6720870')
data_GSM6720872 <- Read10X(data.dir = '/home/ubuntu/data/CP/GSM6720872')
data_GSM6720874 <- Read10X(data.dir = '/home/ubuntu/data/CP/GSM6720874')
data_GSM6720876 <- Read10X(data.dir = '/home/ubuntu/data/CP/GSM6720876')
data_GSM6720878 <- Read10X(data.dir = '/home/ubuntu/data/CP/GSM6720878')
data_GSM6720880 <- Read10X(data.dir = '/home/ubuntu/data/CP/GSM6720880')


sdata_GSM6720852 <-  CreateSeuratObject(data_GSM6720852, project = "GSM6720852")
sdata_GSM6720854 <-  CreateSeuratObject(data_GSM6720854, project = "GSM6720854")
sdata_GSM6720856 <-  CreateSeuratObject(data_GSM6720856, project = "GSM6720856")
sdata_GSM6720858 <-  CreateSeuratObject(data_GSM6720858, project = "GSM6720858")
sdata_GSM6720860 <-  CreateSeuratObject(data_GSM6720860, project = "GSM6720860")
sdata_GSM6720862 <-  CreateSeuratObject(data_GSM6720862, project = "GSM6720862")
sdata_GSM6720864 <-  CreateSeuratObject(data_GSM6720864, project = "GSM6720864")
sdata_GSM6720866 <-  CreateSeuratObject(data_GSM6720866, project = "GSM6720866")
sdata_GSM6720868 <-  CreateSeuratObject(data_GSM6720868, project = "GSM6720868")
sdata_GSM6720870 <-  CreateSeuratObject(data_GSM6720870, project = "GSM6720870")
sdata_GSM6720872 <-  CreateSeuratObject(data_GSM6720872, project = "GSM6720872")
sdata_GSM6720874 <-  CreateSeuratObject(data_GSM6720874, project = "GSM6720874")
sdata_GSM6720876 <-  CreateSeuratObject(data_GSM6720876, project = "GSM6720876")
sdata_GSM6720878 <-  CreateSeuratObject(data_GSM6720878, project = "GSM6720878")
sdata_GSM6720880 <-  CreateSeuratObject(data_GSM6720880, project = "GSM6720880")

# add metadata
sdata_GSM6720852$type = "GSM6720852"
sdata_GSM6720854$type = "GSM6720854"
sdata_GSM6720856$type = "GSM6720856"
sdata_GSM6720858$type = "GSM6720858"
sdata_GSM6720860$type = "GSM6720860"
sdata_GSM6720862$type = "GSM6720862"
sdata_GSM6720864$type = "GSM6720864"
sdata_GSM6720866$type = "GSM6720866"
sdata_GSM6720868$type = "GSM6720868"
sdata_GSM6720870$type = "GSM6720870"
sdata_GSM6720872$type = "GSM6720872"
sdata_GSM6720874$type = "GSM6720874"
sdata_GSM6720876$type = "GSM6720876"
sdata_GSM6720878$type = "GSM6720878"
sdata_GSM6720880$type = "GSM6720880"

rm(data_GSM6720852, data_GSM6720854, data_GSM6720856, data_GSM6720858, data_GSM6720860, data_GSM6720862, data_GSM6720864, data_GSM6720866, data_GSM6720868, data_GSM6720870, data_GSM6720872, data_GSM6720874, data_GSM6720876, data_GSM6720878, data_GSM6720880)

# run garbage collect to free up memory
gc()

# Merge datasets into one single seurat object
alldata <- merge(sdata_GSM6720852, c(sdata_GSM6720854, sdata_GSM6720856, sdata_GSM6720858, sdata_GSM6720860, sdata_GSM6720862, sdata_GSM6720864, sdata_GSM6720866, sdata_GSM6720868, sdata_GSM6720870, sdata_GSM6720872, sdata_GSM6720874, sdata_GSM6720876, sdata_GSM6720878, sdata_GSM6720880), add.cell.ids = c("GSM6720852", "GSM6720854", "GSM6720856", "GSM6720858",
                                                                                                                                                                                                                                                                                                                   "GSM6720860", "GSM6720862", "GSM6720864", "GSM6720866", "GSM6720868", "GSM6720870", "GSM6720872", "GSM6720874", "GSM6720876", "GSM6720878", "GSM6720880"))

# remove all objects that will not be used.
rm(sdata_GSM6720852, sdata_GSM6720854, sdata_GSM6720856, sdata_GSM6720858, sdata_GSM6720860, sdata_GSM6720862, sdata_GSM6720864, sdata_GSM6720866, sdata_GSM6720868, sdata_GSM6720870, sdata_GSM6720872, sdata_GSM6720874, sdata_GSM6720876, sdata_GSM6720878, sdata_GSM6720880)

# run garbage collect to free up memory
gc()

#save alldata
saveRDS(alldata, file = "alldata.rds")
```



```{r}
#filter using the parameters from Ramos et al GEO meta information. (1) Remove doublets (they did with DoubletFinder) and (2) low quality nuceli (< 400 unique genes, < 1,000 UMI counts, or > 15% mitochondrial genes)

#first, calculate percent MT genes
grep("^MT-",rownames(alldata@assays$RNA@counts),value = TRUE)
PercentageFeatureSet(alldata,pattern="^MT-") -> alldata$percent.MT
#then, filter low quality nuclei
subset(
  alldata,
  nFeature_RNA>400 & 
    nCount_RNA > 1000 & 
    percent.MT < 15
) -> alldata_filt

#filter for cells that are in the published meta data
#1 load CP metadata
#2 extract rownames of alldata, use grep to only include nucleotide sequence, attach back as a column?
alldata_cells <- as.vector(row.names(alldata_filt@meta.data))
alldata_cells<- gsub(".*_", "", alldata_cells)
alldata_cells<- gsub("\\-.*","", alldata_cells)
x <- as.data.frame(alldata_filt@meta.data$orig.ident)
x$orig.ident <- x$`alldata_filt@meta.data$orig.ident`
x<- x %>% mutate(sample_id =
                   case_when(orig.ident == "GSM6720852" ~ "9C", 
                             orig.ident == "GSM6720854" ~ "24C",
                             orig.ident == "GSM6720856" ~ "26C",
                             orig.ident == "GSM6720858" ~ "11C",
                             orig.ident == "GSM6720860" ~ "56C",
                             orig.ident == "GSM6720862" ~ "34C",
                             orig.ident == "GSM6720864" ~ "30C",
                             orig.ident == "GSM6720866" ~ "63C",
                             orig.ident == "GSM6720868" ~ "23C",
                             orig.ident == "GSM6720870" ~ "62C",
                             orig.ident == "GSM6720872" ~ "64C",
                             orig.ident == "GSM6720874" ~ "3C",
                             orig.ident == "GSM6720876" ~ "60C",
                             orig.ident == "GSM6720878" ~ "10C",
                             orig.ident == "GSM6720880" ~ "6C",)
)
#having issues with updating the vctrs package and loading dplyr, so did this mutate bit on my computer and uploaded
x_updated<- read.csv("x_update.csv", header = T)
#combine sample and cell name into one column
cell <- as.vector(paste(alldata_cells, x_updated$sample_id, sep="_"))

#make into a vector and add to metadata as well
alldata2<- AddMetaData(alldata_filt, alldata_cells, col.name = "cell_name")#not necessary
alldata2<- AddMetaData(alldata2, x_updated$sample_id, col.name = "sample_id")#not necessary
alldata2<- AddMetaData(alldata2, cell, col.name = "cell")

#3 do the same for the CP metadata (make column with only nuc sequence)
published_meta_CP <- read.csv("/home/ubuntu/GSE217511_CorticalPlate_Seuratmetadata.csv", header = T)
published_meta_CP$cell_name<- gsub("\\-.*","", published_meta_CP$X)
#combine sample and cell name into one column
published_meta_CP$cell <- paste(published_meta_CP$cell_name, published_meta_CP$sample, sep="_")

#4 use subset to filter for cells that were included in the meta data
alldata_filt2<- subset(x = alldata2, subset = cell %in% published_meta_CP$cell) 

#5 merge metadata together, so that you have associated cell types and other info from published metadata
#make smaller version of published meta with only important columns
smaller_pub_meta<- published_meta_CP[,c(5,6,8:12,14)]
#merge with alldata_filt2 meta data using AddMetaData function
alldata_filt2@meta.data$row_names <- rownames(alldata_filt2@meta.data) 
merged_meta<- merge(alldata_filt2@meta.data, smaller_pub_meta, by= "cell")
rownames(merged_meta) <- merged_meta$row_names
alldata_filt3 <- AddMetaData(object = alldata_filt2, metadata = merged_meta)

#save this seurat object for moving forward. It now is filtered for the same cells that were included in the published analysis and includes the meta data info (including cell classification) that was used in the publication
saveRDS(alldata_filt3, file = "alldata_filt3.rds")

###############################################
##try integrating CP dataset looking at maturation modules across stage ("Age")
library(Seurat)
library(SeuratData)
library(patchwork)

#1) go back to seurat and add Age to meta information
#merge with alldata_filt3 meta data using AddMetaData function
extra_info<- data.frame(orig.ident = c("GSM6720852", "GSM6720854", "GSM6720856", "GSM6720858", "GSM6720860", "GSM6720862", "GSM6720864", "GSM6720866", "GSM6720868", "GSM6720870", "GSM6720872", "GSM6720874", "GSM6720876", "GSM6720878", "GSM6720880"), 
                        GW = c("17", "18", "19", "19.7", "20.3", "21.3", "22", "23.3", "23.7", "24.7", "25.8", "26", "32.8", "33.2", "38.3"),
                        Age = c("17-20", "17-20", "17-20", "17-20", "20-24", "20-24", "20-24", "20-24", "20-24", "24-28", "24-28", "24-28", "32-39", "32-39", "32-39"))

alldata_filt3@meta.data$row_names <- rownames(alldata_filt3@meta.data) 
merged_meta<- merge(alldata_filt3@meta.data, extra_info, by= "orig.ident")
rownames(merged_meta) <- merged_meta$row_names
alldata_filt4 <- AddMetaData(object = alldata_filt3, metadata = merged_meta)

#follow "Performing integration on datasets normalized with SCTransform" <https://satijalab.org/seurat/articles/integration_introduction.html>

#2) split merged Seurat object by age bracket (*note: this the un-normalized data)
alldata_filt4.list <- SplitObject(alldata_filt4, split.by = "Age")

#3)then, apply SCT normalization to each Seurat object individually
alldata_filt4.list.list <- lapply(X = alldata_filt4.list, FUN = SCTransform)

#follow the above vingette (**be sure normalization method is set to SCT)
features <- SelectIntegrationFeatures(object.list = alldata_filt4.list.list, nfeatures = 3000)
alldata_filt4.list.list <- PrepSCTIntegration(object.list = alldata_filt4.list.list, anchor.features = features)
alldata_filt4.anchors <- FindIntegrationAnchors(object.list = alldata_filt4.list.list, normalization.method = "SCT",
                                                anchor.features = features)
saveRDS(alldata_filt4.anchors, file = "alldata_filt4.anchors.rds")
alldata_filt4.combined.sct <- IntegrateData(anchorset = alldata_filt4.anchors, normalization.method = "SCT")
saveRDS(alldata_filt4.combined.sct, file = "alldata_filt4.combined.sct.rds")

alldata_filt4.combined.sct2 <- RunPCA(alldata_filt4.combined.sct, verbose = FALSE)
alldata_filt4.combined.sct2 <- RunUMAP(alldata_filt4.combined.sct2, reduction = "pca", dims = 1:30)
p1 <- DimPlot(alldata_filt4.combined.sct2, reduction = "umap", group.by = "Age", shuffle = T, seed = 1)
p2 <- DimPlot(alldata_filt4.combined.sct2, reduction = "umap", group.by = "celltypes", label = TRUE,
              repel = TRUE)
pdf(file = "integrated_clusters.pdf", width = 11, height = 8.5)
p1 + p2
dev.off()

#plot EOMES expression for reviewer
DefaultAssay(alldata_filt4.combined.sct2) <- "RNA"
alldata_filt4.combined.sct2_2<- NormalizeData(alldata_filt4.combined.sct2)
VlnPlot(alldata_filt4.combined.sct2_2, features = "EOMES", group.by = "celltypes", assay = "RNA")
FeaturePlot(alldata_filt4.combined.sct2_2, features = "EOMES", slot = "data", order = T)
#Ridgeplot not great for subtle differences
library(ggridges, lib.loc = "/usr/local/lib/R/site-library")
RidgePlot(alldata_filt4.combined.sct2_2, features = "EOMES", group.by = "celltypes", assay = "RNA", slot = "data") + geom_density_ridges2(bandwidth = 0.5)

#make another version of umap with celltypes where all clusters are gray except for AC and gIPC
my_cols <- c("L2/3 CPN"="#c4c2c2",  "AC"="#F8766D","nIPC"="#c4c2c2","UD"="#c4c2c2","TAC"="#c4c2c2","IN"="#c4c2c2","SPN"="#c4c2c2","L6 CPN"="#c4c2c2","L4/5a CPN"="#c4c2c2", "MG"="#c4c2c2","gIPC"="#ABA300","CRN"="#c4c2c2","BVC"="#c4c2c2","OPC"="#c4c2c2")

p2<- DimPlot(alldata_filt4.combined.sct2, reduction = "umap", group.by = "celltypes", label = TRUE,
             repel = TRUE, cols = my_cols)
pdf(file = "integrated_clusters_12_5.pdf", width = 11, height = 8.5)
p1 + p2
dev.off()

#subset to only include glia clusters
alldata_filt4.combined.sct2_subset1<- subset(x = alldata_filt4.combined.sct2, subset = celltypes %in% c("AC", "TAC", "gIPC", "OPC"))
#try susetting for only AC and gIPC clusters
alldata_filt4.combined.sct2_subset2<- subset(x = alldata_filt4.combined.sct2, subset = celltypes %in% c("AC", "gIPC"))

#split by "Age" to look at relevant glia clusters across age
pdf(file = "integrated_glia_Age_11_7.pdf", width = 11, height = 8.5)
DimPlot(alldata_filt4.combined.sct2_subset, reduction = "umap", split.by = "Age", group.by = "celltypes", shuffle = T, seed = 1)
dev.off()

#try not splitting by Age (updated:11/30/23)
#edit 12/5/32: make colors the same as in large UMAP
my_cols2 <- c("AC"="#F8766D","gIPC"="#ABA300","OPC"="#a58aff","TAC"="#fc61d7")
my_cols3 <- c("17-20"="#c6dbef","20-24"="#6caed6","24-28"="#2370b5","32-39"="#07316b")

p1 <- DimPlot(alldata_filt4.combined.sct2_subset1, reduction = "umap", group.by = "Age", shuffle = T, seed = 1, cols = my_cols3)
p2 <- DimPlot(alldata_filt4.combined.sct2_subset1, reduction = "umap", group.by = "celltypes", label = TRUE,
              repel = TRUE, cols = my_cols2)
pdf(file = "integrated_clusters_subset_12_5.pdf", width = 11, height = 8.5)
p1 + p2
dev.off()


#try moving outliars to center of plot
alldata_filt4.combined.sct2_subset[["umap_new"]] <- CollapseEmbeddingOutliers(alldata_filt4.combined.sct2_subset,
                                                                              reduction = "umap", reduction.key = 'umap_', outlier.sd = 0.5)

pdf(file = "integrated_glia_Age_CollapseEmbeddingOutliers.pdf", width = 11, height = 8.5)
DimPlot(alldata_filt4.combined.sct2_subset, reduction = "umap_new", split.by = "Age", group.by = "celltypes", shuffle = T, seed = 1)
dev.off()
#--> this didn't really take care of the stray points

####to plot genes we need to switch from integrated to SCT
#note: while it is suggested to use log-normalized RNA for plotting individual genes (see https://github.com/satijalab/seurat/issues/4082), I will be calculating module scores. I tried using both SCT and log-normalized RNA assays as input for calculating module scores and the patterns don't change. The min and max module scores barely changed. Could probably do either. If plotting individual genes, use log-normalized RNA (RNA@data) because we applied SCT to samples individually, then integrated.
DefaultAssay(alldata_filt4.combined.sct2_subset1) <- "SCT"
DefaultAssay(alldata_filt4.combined.sct2_subset1) <- "RNA"
alldata_filt4.combined.sct2_subset1<- NormalizeData(alldata_filt4.combined.sct2_subset1)
#feature plots of maturation module scores umaps split by age
early<- as.list(org_RNAseq_gene_modules[org_RNAseq_gene_modules$colors== "Early",1])
early_middle<- as.list(org_RNAseq_gene_modules[org_RNAseq_gene_modules$colors== "Early_Middle",1])
middle<- as.list(org_RNAseq_gene_modules[org_RNAseq_gene_modules$colors== "Middle",1])
middle_late<- as.list(org_RNAseq_gene_modules[org_RNAseq_gene_modules$colors== "Middle_Late",1])
late<- as.list(org_RNAseq_gene_modules[org_RNAseq_gene_modules$colors== "Late",1])

###top module genes
library(readxl)
top10 <- read_xlsx("topMods.xlsx", sheet = 1, col_names = F)
top25 <- read_xlsx("topMods.xlsx", sheet = 2, col_names = F)
top50 <- read_xlsx("topMods.xlsx", sheet = 3, col_names = F)

org_top10<- org_RNAseq_gene_modules[org_RNAseq_gene_modules$gene_id %in% top10$...1,]
org_top25<- org_RNAseq_gene_modules[org_RNAseq_gene_modules$gene_id %in% top25$...1,]
org_top50<- org_RNAseq_gene_modules[org_RNAseq_gene_modules$gene_id %in% top50$...1,]

early<- as.list(org_top50[org_top50$colors== "Early",1])
early_middle<- as.list(org_top50[org_top50$colors== "Early_Middle",1])
middle<- as.list(org_top50[org_top50$colors== "Middle",1])
middle_late<- as.list(org_top50[org_top50$colors== "Middle_Late",1])
late<- as.list(org_top50[org_top50$colors== "Late",1])

####
#11/7/23: added nbin cutoff of 12 after subsetting to only include AC and gIPC because without nIPC and OPC, we lose cells with reads and no longer have enough cells with reads for module genes to use the default cutoff of 24
object_early <- AddModuleScore(object = alldata_filt4.combined.sct2_subset, features = early, name = "early_score", nbin = 12)
object_early_middle <- AddModuleScore(object = alldata_filt4.combined.sct2_subset, features = early_middle, name = "early_middle_score", nbin = 12)
object_middle <- AddModuleScore(object = alldata_filt4.combined.sct2_subset, features = middle, name = "middle_score", nbin = 12)
object_middle_late <- AddModuleScore(object = alldata_filt4.combined.sct2_subset, features = middle_late, name = "middle_late_score", nbin = 12)
object_late <- AddModuleScore(object = alldata_filt4.combined.sct2_subset, features = late, name = "late_score", nbin = 12)

#11/30/23: try module scores for AC, gIPC, OPC, and TAC
object_early <- AddModuleScore(object = alldata_filt4.combined.sct2_subset1, features = early, name = "early_score")
object_early_middle <- AddModuleScore(object = alldata_filt4.combined.sct2_subset1, features = early_middle, name = "early_middle_score")
object_middle <- AddModuleScore(object = alldata_filt4.combined.sct2_subset1, features = middle, name = "middle_score")
object_middle_late <- AddModuleScore(object = alldata_filt4.combined.sct2_subset1, features = middle_late, name = "middle_late_score")
object_late <- AddModuleScore(object = alldata_filt4.combined.sct2_subset1, features = late, name = "late_score")

#11/30/23: try module scores for just AC and gIPC. had to change bin size for fewer cells for SCT. For RNA, leave as default
object_all <- AddModuleScore(object = alldata_filt4.combined.sct2_subset2, features = early, name = "early_score", nbin = 12)
object_all <- AddModuleScore(object = object_all, features = early_middle, name = "early_middle_score", nbin = 12)
object_all <- AddModuleScore(object = object_all, features = middle, name = "middle_score", nbin = 12)
object_all <- AddModuleScore(object = object_all, features = middle_late, name = "middle_late_score", nbin = 12)
object_all <- AddModuleScore(object = object_all, features = late, name = "late_score", nbin = 12)

#wanted to keep color scale the same across Age plots with "keep.scale" option in feature plots, but this option does not work if you try to customize the color palette
pdf(file = "integrated_glia_Age_early_module.pdf", width = 15, height = 8.5)
FeaturePlot(object = object_early, features = 'early_score1', split.by = "Age", order = T, keep.scale = "all") & scale_color_viridis_c() & theme(legend.position = "right")
dev.off()
#try another approach for these graphs
library(tidyverse)
library(patchwork)
library(viridis)
library(Seurat)
library(scCustomize)
library(RColorBrewer)

pal <- viridis(n = 10, option = "D")

#the package says this about the na_cutoff "scCustomize contains a parameter called na_cutoff which tells the function which values to plot as background. By default this is set to value that means background is treated as 0 or below. Depending on what feature, assay, or value you are interested in this parameter should be modified appropriately.For instance if plotting module score which contains negative values you will probably want to remove the cutoff value entirely to avoid misconstruing results."

pdf(file = "integrated_glia_Age_early_module.pdf", width = 15, height = 8.5)
e<- FeaturePlot_scCustom(seurat_object = object_early, features = "early_score1", split.by = "Age", order = T, na_cutoff = NA)
dev.off()

pdf(file = "integrated_glia_Age_early_middle_module.pdf", width = 15, height = 8.5)
FeaturePlot_scCustom(seurat_object = object_early_middle, features = "early_middle_score1", split.by = "Age", order = T, na_cutoff = NA)
dev.off()

pdf(file = "integrated_glia_Age_middle_module.pdf", width = 15, height = 8.5)
FeaturePlot_scCustom(seurat_object = object_middle, features = "middle_score1", split.by = "Age", order = T, na_cutoff = NA)
dev.off()

pdf(file = "integrated_glia_Age_middle_late_module.pdf", width = 15, height = 8.5)
FeaturePlot_scCustom(seurat_object = object_middle_late, features = "middle_late_score1", split.by = "Age", order = T, na_cutoff = NA)
dev.off()

pdf(file = "integrated_glia_Age_late_module.pdf", width = 15, height = 8.5)
FeaturePlot_scCustom(seurat_object = object_late, features = "late_score1",  split.by = "Age", order = T, na_cutoff = NA)
dev.off()

#new featureplots that are not split by Age and all have the same scale
#Get the range of scores across all gene sets
score_range <- c(
  min(c(object_early@meta.data$early_score1, object_early_middle@meta.data$early_middle_score1,
        object_middle@meta.data$middle_score1, object_middle_late@meta.data$middle_late_score1,
        object_late@meta.data$late_score1)),
  max(c(object_early@meta.data$early_score1, object_early_middle@meta.data$early_middle_score1,
        object_middle@meta.data$middle_score1, object_middle_late@meta.data$middle_late_score1,
        object_late@meta.data$late_score1))
)
#for RNA: [1] -0.1047000  0.4553764
#for SCT: [1] -0.07956151  0.33834284

ylgnbu_palette <- brewer.pal(9, "YlGnBu")

FeaturePlot(object_early, features = "early_score1", cols = ylgnbu_palette, 
            min.cutoff = score_range[1], max.cutoff = score_range[2], order = T)

FeaturePlot(object_early_middle, features = "early_middle_score1", cols = ylgnbu_palette, min.cutoff = score_range[1], max.cutoff = score_range[2], order = T)

FeaturePlot(object_middle, features = "middle_score1", cols = ylgnbu_palette, 
            min.cutoff = score_range[1], max.cutoff = score_range[2], order = T)

FeaturePlot(object_middle_late, features = "middle_late_score1", cols = ylgnbu_palette, 
            min.cutoff = score_range[1], max.cutoff = score_range[2], order = T)

FeaturePlot(object_late, features = "late_score1", cols = ylgnbu_palette, 
            min.cutoff = score_range[1], max.cutoff = score_range[2], order = T)


##try specific genes
FeaturePlot(object = alldata_filt4.combined.sct2_subset1, features = c("EGFR", "OLIG2", "OLIG1", "SOX9"), order = T) & scale_color_viridis_c() & theme(legend.position = "right")

#diff colors
FeaturePlot(object = alldata_filt4.combined.sct2_subset1, features = c("EGFR", "OLIG2", "OLIG1", "SOX9"), order = T, cols = ylgnbu_palette) & theme(legend.position = "right")


#try violin plots because of the outlier issue with featureplots
pdf(file = "integrated_glia_Age_early_module_Vln.pdf", width = 15, height = 8.5)
VlnPlot_scCustom(seurat_object = object_early, features = "early_score1", split.by = "celltypes", group.by = "Age", ggplot_default_colors = T)
dev.off()

pdf(file = "integrated_glia_Age_early_middle_module_Vln.pdf", width = 15, height = 8.5)
VlnPlot_scCustom(seurat_object = object_early_middle, features = "early_middle_score1", split.by = "celltypes", group.by = "Age", ggplot_default_colors = T)
dev.off()

pdf(file = "integrated_glia_Age_middle_module_Vln.pdf", width = 15, height = 8.5)
VlnPlot_scCustom(seurat_object = object_middle, features = "middle_score1", split.by = "celltypes", group.by = "Age", ggplot_default_colors = T)
dev.off()

pdf(file = "integrated_glia_Age_middle_late_module_Vln.pdf", width = 15, height = 8.5)
VlnPlot_scCustom(seurat_object = object_middle_late, features = "middle_late_score1", split.by = "celltypes", group.by = "Age", ggplot_default_colors = T)
dev.off()

pdf(file = "integrated_glia_Age_late_module_Vln.pdf", width = 15, height = 8.5)
VlnPlot_scCustom(seurat_object = object_late, features = "late_score1", split.by = "celltypes", group.by = "Age", ggplot_default_colors = T)
dev.off()

##updated 11/30/23: make violin plots showing module score across age bin
#only include AC and gIPC data (alldata_filt4.combined.sct2_subset2), but combine them in plot, and make scales 
#across the five plots the same.

#Get the range of scores across all gene sets

object_all$Age <- factor(object_all$Age, levels = c("17-20", "20-24", "24-28", "32-39"))

VlnPlot(object_all, features = c("early_score1", "early_middle_score1", "middle_score1", "middle_late_score1", "late_score1"), group.by = "Age", stack = T, same.y.lims = T, ncol = 1, flip = T)




##plot maturation TFs
maturation_TFs <- c("OTX2", "ZIC4", "ZIC1", "ZIC3", "NHLH1", "NEUROD1", "PKNOX2", "RXRG", "HEY1", "PAX3", "LMX1A", "MSX2", "SOX15", "POU3F3", "TFAP2C", "POU3F4", "FOXG1", "LHX2", "POU3F2", "MYCN", "SOX8", "OLIG1", "OLIG2", "EOMES", "ASCL1", "SOX21", "RFX4", "PRRX1", "NR3C2")

pdf(file = "maturation_TFs_feature.pdf", width = 15, height = 8.5)
p5<- FeaturePlot_scCustom(seurat_object = alldata_filt4.combined.sct2_subset, features = maturation_TFs, split.by = "Age", order = T, na_cutoff = NA)
dev.off()


###plot maturation modules in whole UMAP
#attempted this with all genes and top 50. The patterns look nearly identical, so went with all genes
object_early <- AddModuleScore(object = alldata_filt4.combined.sct2, features = early, name = "early_score")
object_early_middle <- AddModuleScore(object = alldata_filt4.combined.sct2, features = early_middle, name = "early_middle_score")
object_middle <- AddModuleScore(object = alldata_filt4.combined.sct2, features = middle, name = "middle_score")
object_middle_late <- AddModuleScore(object = alldata_filt4.combined.sct2, features = middle_late, name = "middle_late_score")
object_late <- AddModuleScore(object = alldata_filt4.combined.sct2, features = late, name = "late_score")


pdf(file = "integrated_allClusters_early_module.pdf", width = 15, height = 8.5)
FeaturePlot_scCustom(seurat_object = object_early, features = "early_score1", order = T, na_cutoff = NA)
dev.off()

pdf(file = "integrated_allClusters_early_middle_module.pdf", width = 15, height = 8.5)
FeaturePlot_scCustom(seurat_object = object_early_middle, features = "early_middle_score1", order = T, na_cutoff = NA)
dev.off()

pdf(file = "integrated_allClusters_middle_module.pdf", width = 15, height = 8.5)
FeaturePlot_scCustom(seurat_object = object_middle, features = "middle_score1", order = T, na_cutoff = NA)
dev.off()

pdf(file = "integrated_allClusters_middle_late_module.pdf", width = 15, height = 8.5)
FeaturePlot_scCustom(seurat_object = object_middle_late, features = "middle_late_score1", order = T, na_cutoff = NA)
dev.off()

pdf(file = "integrated_allClusters_late_module.pdf", width = 15, height = 8.5)
FeaturePlot_scCustom(seurat_object = object_late, features = "late_score1",  order = T, na_cutoff = NA)
dev.off()




####find genes enriched in AS and gIPCs

#ignore this because it is integrated
alldata_filt4.combined.sct3<- PrepSCTFindMarkers(alldata_filt4.combined.sct2)
saveRDS(alldata_filt4.combined.sct3, file = "alldata_filt4.combined.sct3.rds")

alldata_filt3 <- readRDS("~/alldata_filt3.rds")
library(Seurat)
NormalizeData(alldata_filt3, normalization.method = "LogNormalize", scale.factor = 10000) -> norm_data
#find markers of astrocytes
Idents(norm_data) <- "celltypes"
AC.markers <- FindMarkers(norm_data, ident.1 = "AC")
gIPC.markers <- FindMarkers(norm_data, ident.1 = "gIPC")


##################################
#pseudobulk 
#not split by age
# Assuming you have already loaded Seurat and have a Seurat object named 'seuratObj'

# Install and load necessary libraries
# install.packages("Matrix")
# install.packages("dplyr")
# install.packages("Biobase")
# install.packages("Matrix.utils")
library(Matrix)
library(dplyr)
library(Biobase)
library(Matrix.utils)

pseudo_age <- AggregateExpression(alldata_filt4, assays = "RNA", return.seurat = T, group.by = c("GW", "celltypes"))

pseudo_age_cts <- AggregateExpression(alldata_filt4, assays = "RNA", return.seurat = T, group.by = c("GW", "celltypes"), slot = "counts")

pseudo_age_df <- as.data.frame(pseudo_age@assays$RNA@data)
pseudo_age_cts_df <- as.data.frame(pseudo_age_cts@assays$RNA@data)

pseudo <- AggregateExpression(alldata_filt4, assays = "RNA", return.seurat = T, group.by = c("celltypes"))
pseudo_cts <- AggregateExpression(alldata_filt4, assays = "RNA", return.seurat = T, group.by = c("celltypes"), slot = "counts")

pseudo_df <- as.data.frame(pseudo@assays$RNA@data)
pseudo_cts_df <- as.data.frame(pseudo_cts@assays$RNA@data)

pseudo_bin <- AggregateExpression(alldata_filt4, assays = "RNA", return.seurat = T, group.by = c("Age", "celltypes"))
pseudo_bin_cts <- AggregateExpression(alldata_filt4, assays = "RNA", return.seurat = T, group.by = c("Age", "celltypes"), slot = "counts")

pseudo_bin_df <- as.data.frame(pseudo_bin@assays$RNA@data)
pseudo_bin_cts_df <- as.data.frame(pseudo_bin_cts@assays$RNA@data)

write.csv(pseudo_age_cts_df, file = "pseudo_age_cts_df.csv")
write.csv(pseudo_bin_cts_df, file = "pseudo_bin_cts_df.csv")
write.csv(pseudo_cts_df, file = "pseudo_cts_df.csv")

#first, subset by cell type
alldata_filt4.AC<- subset(x = alldata_filt4, subset = celltypes %in% c("AC"))

alldata_filt4.gIPC<- subset(x = alldata_filt4, subset = celltypes %in% c("gIPC"))

# Function to generate pseudobulk count matrix
generatePseudobulkMatrix <- function(seurat_obj) {
  # Get the count matrix from Seurat object
  counts_mat <- t(seurat_obj@assays$RNA@counts)
  
  # Aggregate counts to get pseudobulk counts
  pseudobulk_counts <- Matrix::colSums(counts_mat)
  
  # Create a data frame with pseudobulk counts
  pseudobulk_df <- as.data.frame(pseudobulk_counts)
  colnames(pseudobulk_df) <- "PseudobulkCounts"
  
  return(pseudobulk_df)
}

# Generate pseudobulk count matrix
pseudobulk_matrix_all_AC <- generatePseudobulkMatrix(alldata_filt4.AC)#raw AC
pseudobulk_matrix_all_gIPC <- generatePseudobulkMatrix(alldata_filt4.gIPC)

# Print first few rows of the pseudobulk count matrix
head(pseudobulk_matrix)

#generate different pseudobulk matrices for age bins

# Function to generate pseudobulk count matrix for specific groups
generatePseudobulkMatrixByGroup <- function(seurat_obj, group_column, group_value) {
  # Get the count matrix from Seurat object
  counts_mat <- t(seurat_obj@assays$RNA@counts)
  
  # Get cell IDs for the specific group value
  cells_in_group <- which(seurat_obj@meta.data[[group_column]] == group_value)
  
  # Subset the count matrix for the specific group
  counts_subset <- counts_mat[, cells_in_group]
  
  # Aggregate counts to get pseudobulk counts
  pseudobulk_counts <- Matrix::colSums(counts_subset)
  
  # Create a data frame with pseudobulk counts
  pseudobulk_df <- as.data.frame(pseudobulk_counts)
  colnames(pseudobulk_df) <- paste("PseudobulkCounts_", group_value, sep = "")
  
  return(pseudobulk_df)
}

## Generate pseudobulk count matrix for time bins
pseudobulk_matrix_17_20_AC <- generatePseudobulkMatrixByGroup(alldata_filt4.AC, "Age", "17-20")

pseudobulk_matrix_20_24_AC <- generatePseudobulkMatrixByGroup(alldata_filt4.AC, "Age", "20-24")

pseudobulk_matrix_24_28_AC <- generatePseudobulkMatrixByGroup(alldata_filt4.AC, "Age", "24-28")

pseudobulk_matrix_32_39_AC <- generatePseudobulkMatrixByGroup(alldata_filt4.AC, "Age", "32-39")

#now for gIPC
pseudobulk_matrix_17_20_gIPC <- generatePseudobulkMatrixByGroup(alldata_filt4.gIPC, "Age", "17-20")

pseudobulk_matrix_20_24_gIPC <- generatePseudobulkMatrixByGroup(alldata_filt4.gIPC, "Age", "20-24")

pseudobulk_matrix_24_28_gIPC <- generatePseudobulkMatrixByGroup(alldata_filt4.gIPC, "Age", "24-28")

pseudobulk_matrix_32_39_gIPC <- generatePseudobulkMatrixByGroup(alldata_filt4.gIPC, "Age", "32-39")

write.csv(pseudobulk_matrix_17_20_AC, file = "AC_pseudo_17_20.csv")

write.csv(pseudobulk_matrix_20_24_AC, file = "AC_pseudo_20_24.csv")

write.csv(pseudobulk_matrix_24_28_AC, file = "AC_pseudo_24_28.csv")

write.csv(pseudobulk_matrix_32_39_AC, file = "AC_pseudo_32_39.csv")

write.csv(pseudobulk_matrix_17_20_gIPC, file = "gIPC_pseudo_17_20.csv")

write.csv(pseudobulk_matrix_20_24_gIPC, file = "gIPC_pseudo_20_24.csv")

write.csv(pseudobulk_matrix_24_28_gIPC, file = "gIPC_pseudo_24_28.csv")

write.csv(pseudobulk_matrix_32_39_gIPC, file = "gIPC_pseudo_32_39.csv")

write.csv(pseudobulk_matrix_all_AC, file = "AC_pseudo.csv")

write.csv(pseudobulk_matrix_all_gIPC, file = "gIPC_pseudo.csv")

######################################################################
alldata_filt4.combined.sct <- readRDS("~/alldata_filt4.combined.sct.rds")


alldata_filt4.combined.sct2 <- RunPCA(alldata_filt4.combined.sct, verbose = FALSE)
alldata_filt4.combined.sct2 <- RunUMAP(alldata_filt4.combined.sct2, reduction = "pca", dims = 1:30)
p1 <- DimPlot(alldata_filt4.combined.sct2, reduction = "umap", group.by = "Age", shuffle = T, seed = 1)
p2 <- DimPlot(alldata_filt4.combined.sct2, reduction = "umap", group.by = "celltypes", label = TRUE,
              repel = TRUE)
pdf(file = "integrated_clusters.pdf", width = 11, height = 8.5)
p1 + p2
dev.off()

#subset to only include glia clusters
alldata_filt4.combined.sct2_subset1<- subset(x = alldata_filt4.combined.sct2, subset = celltypes %in% c("AC", "TAC", "gIPC", "OPC"))


### 7/23/24: Update for reviewer
#Subcluster SCT-normalized alldata_filt4.combined.sct2_subset1 and re-plot module scores

scfp <- FindNeighbors(alldata_filt4.combined.sct2_subset1, graph.name = "test", dims = 1:10)
scfp <- FindClusters(scfp, graph.name = "test", resolution = 0.5, algorithm = 1, verbose = TRUE)
scfp<- RunPCA(scfp, verbose = FALSE)
scfp <- RunUMAP(scfp, reduction = "pca", dims = 1:30)
DimPlot(scfp, reduction = "umap",  group.by = "celltypes", label = TRUE,
        repel = TRUE, cols = my_cols2)
DimPlot(scfp, reduction = "umap",  group.by = "test_res.0.5", label = TRUE,
        repel = TRUE)

# Check current cluster identities
levels(Idents(scfp))

# Rename clusters
new_cluster_ids <- c("AS_1", "gIPC_1", "OPC_1", "AS_2", "TAC_1", "AS_3", "OPC_2", "AS_4", "gIPC_2", "TAC_2", "AC_5", "OPC_3")
names(new_cluster_ids) <- levels(scfp)
scfp2 <- RenameIdents(scfp, new_cluster_ids)
DimPlot(scfp2, reduction = "umap", label = TRUE,
        +         repel = TRUE, shuffle = T, seed = 1)
DimPlot(scfp2, reduction = "umap",  group.by = "new_cluster_ids", label = TRUE,
        repel = TRUE)
DimPlot(scfp2, reduction = "umap",  group.by = "celltypes", label = TRUE,
        repel = TRUE, cols = my_cols2)
DimPlot(scfp2, reduction = "umap",  group.by = "Age", label = TRUE,
        repel = TRUE, shuffle = T, seed = 1)
saveRDS(scfp2, file = "scfp2.rds")


###
#plot module scores
DefaultAssay(scfp2) <- "RNA"
scfp2<- NormalizeData(scfp2)

early<- as.list(org_RNAseq_gene_modules[org_RNAseq_gene_modules$colors== "Early",1])
early_middle<- as.list(org_RNAseq_gene_modules[org_RNAseq_gene_modules$colors== "Early_Middle",1])
middle<- as.list(org_RNAseq_gene_modules[org_RNAseq_gene_modules$colors== "Middle",1])
middle_late<- as.list(org_RNAseq_gene_modules[org_RNAseq_gene_modules$colors== "Middle_Late",1])
late<- as.list(org_RNAseq_gene_modules[org_RNAseq_gene_modules$colors== "Late",1])

object_early <- AddModuleScore(object = scfp2, features = early, name = "early_score")
object_early_middle <- AddModuleScore(object = scfp2, features = early_middle, name = "early_middle_score")
object_middle <- AddModuleScore(object = scfp2, features = middle, name = "middle_score")
object_middle_late <- AddModuleScore(object = scfp2, features = middle_late, name = "middle_late_score")
object_late <- AddModuleScore(object = scfp2, features = late, name = "late_score")


#new featureplots that are not split by Age and all have the same scale
#Get the range of scores across all gene sets
score_range <- c(
  min(c(object_early@meta.data$early_score1, object_early_middle@meta.data$early_middle_score1,
        object_middle@meta.data$middle_score1, object_middle_late@meta.data$middle_late_score1,
        object_late@meta.data$late_score1)),
  max(c(object_early@meta.data$early_score1, object_early_middle@meta.data$early_middle_score1,
        object_middle@meta.data$middle_score1, object_middle_late@meta.data$middle_late_score1,
        object_late@meta.data$late_score1))
)
#for RNA: [1] -0.1047000  0.4553764
#for SCT: [1] -0.07956151  0.33834284

ylgnbu_palette <- brewer.pal(9, "YlGnBu")

FeaturePlot(object_early, features = "early_score1", cols = ylgnbu_palette, 
            min.cutoff = score_range[1], max.cutoff = score_range[2], order = T)

FeaturePlot(object_early_middle, features = "early_middle_score1", cols = ylgnbu_palette, min.cutoff = score_range[1], max.cutoff = score_range[2], order = T)

FeaturePlot(object_middle, features = "middle_score1", cols = ylgnbu_palette, 
            min.cutoff = score_range[1], max.cutoff = score_range[2], order = T)

FeaturePlot(object_middle_late, features = "middle_late_score1", cols = ylgnbu_palette, 
            min.cutoff = score_range[1], max.cutoff = score_range[2], order = T)

FeaturePlot(object_late, features = "late_score1", cols = ylgnbu_palette, 
            min.cutoff = score_range[1], max.cutoff = score_range[2], order = T)

#######################################################################
#sankey plot to compare WGCNA modules to NMF groups

library(plotly)
library(dplyr)

#below is for 5 NMF groups
#create a basic sankey
p <- plot_ly(
  type = "sankey",
  orientation = "h",
  
  #each element is a node here...V1 is node 0,
  #V2 is node 1,
  #other is node 10
  node = list(
    label = c("V1", "V2", "V3", 
              "V4", "V5",
              "early","early_middle","middle","middle_late", "late"),
    line = list(
      color = "black",
      width = 0.5
    )
  ),
  
  link = list(
    #All channels become the sources...so, nodes 0-4
    source = c(0, 0,0,0,0,
               1,1,1,1,1,
               2,2,2,2,2,
               3,3,3,3,3,
               4,4,4,4,4),
    #Landing pages become the target...so, nodes 5-10
    target = c(5,6,7,8,9,
               5,6,7,8,9,
               5,6,7,8,9,
               5,6,7,8,9,
               5,6,7,8,9),
    #Assigning values between nodes...
    value =  c(2,0,1,78,41,
               5,40,16,10,0,
               49,11,14,0,0,
               2,3,2,54,30,
               1,9,60,39,0
    )
  )
) %>% 
  layout(
    title = "maturation module enrichment in NMF groups",
    font = list(
      size = 10
    )
  )
p

#saving as a vectorized file type (only could figure out svg)
library(reticulate)
reticulate::install_miniconda()
reticulate::conda_install('r-reticulate', 'python-kaleido')
reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
reticulate::use_miniconda('r-reticulate')
library(plotly)
p %>%
  config(
    toImageButtonOptions = list(
      format = "svg",
      filename = "myplot",
      width = 600,
      height = 700
    )
  )
#when the plot window pops up, hover over the top right corner of plot and click
#the camera icon to download plot. it should now download as an svg instead of png


#########################################################################
#correlate organoid and pseudobulk Ramos data

#read in pseudobulk data
age_pseudo<- read.csv(file = "C:/Users/caitl/Documents/third rotation/pseudotime/pseudo_age_cts_df.csv")
age_pseudo<- read.csv(file = "C:/Users/caitl/Documents/third rotation/pseudotime/pseudo_cts_df.csv")

#read in organoid count matrix
library(readxl)
org_cts<- read_xlsx("C:/Users/caitl/Documents/third rotation/Organoid maturation timeline/HiSeq/RNA-seq/master_featurecounts.xlsx", col_names = T)

#filter age_pseudo to only keep genes that are highly expressed in gIPC
gIPC_genes<- read.csv("C:/Users/caitl/Documents/third rotation/pseudotime/gIPC_markers.csv")
gIPC_genes_subset <- gIPC_genes[gIPC_genes$p_val_adj < 0.01,]
#sort for the top 50, 100, and 500 genes (sorted by FC)
gIPC_genes_subset <- gIPC_genes_subset[order(gIPC_genes_subset$avg_log2FC, decreasing = T),]
gIPC_genes_subset2 <- gIPC_genes_subset[1:50,]
gIPC_genes_subset2 <- gIPC_genes_subset[1:100,] #ended up using this
gIPC_genes_subset2 <- gIPC_genes_subset[1:500,]

library(dplyr)
age_pseudo <- age_pseudo[age_pseudo$X %in% gIPC_genes_subset2$X,]
rownames(age_pseudo) <- age_pseudo$X

# Filter columns containing 'gIPC' in their names
gIPC_age_pseudo <- age_pseudo %>% 
  select(contains("gIPC"))
gIPC_age_pseudo <- as.data.frame(gIPC_age_pseudo)

org_cts <- as.data.frame(org_cts)
org_cts <- org_cts[org_cts$Geneid %in% rownames(gIPC_age_pseudo),]
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

#merge the two and make into a matrix
gIPC_age_pseudo <- as.data.frame(gIPC_age_pseudo)
gIPC_age_pseudo$Gene_id <- rownames(gIPC_age_pseudo)

org_tpm <- as.data.frame(org_tpm)
org_tpm$Gene_id <- rownames(org_tpm)

merged_data <- merge(org_tpm, gIPC_age_pseudo, by = "Gene_id")

rownames(merged_data) <- merged_data$Gene_id
merged_data <- merged_data[,-1]
merged_data_mtx <- as.matrix(merged_data)


library(Hmisc)
cor_result <- round(cor(merged_data_mtx, method="spearman"),2)
a <- cor_result[-c(1:21),-22]
scaled <- as.data.frame(scale(a))
scaled2 <- as.vector(scaled)
#min: -1.833302
#max: 1.629602
write.csv(scaled2, file = "hCS_gIPC_cor.csv")

remotes::install_github("jmw86069/colorjam")
library(colorjam)
scaled$sample<- rownames(scaled)
x <- jamba::nameVector(scaled, y = scaled$sample)
jamba::showColors(vals2colorLevels(x, col = "BuPu", divergent = T))

#check out <https://jmw86069.github.io/colorjam/index.html> to make own color palette

########################

#try correlation with AS split by age
age_pseudo<- read.csv(file = "C:/Users/caitl/Documents/third rotation/pseudotime/pseudo_bin_cts_df.csv")

#read in organoid count matrix
library(readxl)
org_cts<- read_xlsx("C:/Users/caitl/Documents/third rotation/Organoid maturation timeline/HiSeq/RNA-seq/master_featurecounts.xlsx", col_names = T)

#filter age_pseudo to only keep genes that are highly expressed in gIPC
AC_genes<- read.csv("C:/Users/caitl/Documents/third rotation/pseudotime/AC_markers.csv")
AC_genes_subset <- AC_genes[AC_genes$p_val_adj < 0.01,]
#sort for the top 50, 100, and 500 genes (sorted by FC)
AC_genes_subset <- AC_genes_subset[order(AC_genes_subset$avg_log2FC, decreasing = T),]
AC_genes_subset2 <- AC_genes_subset[1:50,]
AC_genes_subset2 <- AC_genes_subset[1:100,]#try this first
AC_genes_subset2 <- AC_genes_subset[1:500,]

library(dplyr)
age_pseudo <- age_pseudo[age_pseudo$X %in% AC_genes_subset2$X,]
rownames(age_pseudo) <- age_pseudo$X

# Filter columns containing 'AC' in their names
AC_age_pseudo <- age_pseudo %>% 
  select(contains("AC"))
AC_age_pseudo <- AC_age_pseudo %>% 
  select(-contains("TAC"))


org_cts <- as.data.frame(org_cts)
org_cts <- org_cts[org_cts$Geneid %in% rownames(AC_age_pseudo),]
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

#merge the two and make into a matrix
AC_age_pseudo <- as.data.frame(AC_age_pseudo)
AC_age_pseudo$Gene_id <- rownames(AC_age_pseudo)

org_tpm <- as.data.frame(org_tpm)
org_tpm$Gene_id <- rownames(org_tpm)

merged_data <- merge(org_tpm, AC_age_pseudo, by = "Gene_id")

rownames(merged_data) <- merged_data$Gene_id
merged_data <- merged_data[,-1]
merged_data_mtx <- as.matrix(merged_data)


library(Hmisc)
cor_result <- round(cor(merged_data_mtx, method="spearman"),2)
a <- cor_result[-c(1:21),-c(22:25)]

#use same color palette as for gIPCs
library(RColorBrewer)
heatmap.2(t(as.matrix(a)), dendrogram= "none",Colv = F, Rowv = F,colRow = F, trace = NULL,scale = "column",col=brewer.pal(n = 9, name = "BuPu"))
write.csv(t(as.matrix(a)), file = "hCS_AC_cor.csv")
