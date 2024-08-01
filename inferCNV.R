#CNV signal correlation with tissue mean signal. For tumor, cells that correlate well
#with tissue mean are more likely to be malignant. Also the reverse (cells that don't correlate with tissue mean)
#for cells in margin tissue

library(infercnv)

#import extracted objects from ArchR for inferCNV input
cellMeta_all2 <- readRDS("~/cellMeta_all2.rds")
expMat_all_Assay2 <- readRDS("~/expMat_all_Assay2.rds")
genePos.all <- readRDS("~/genePos.all.rds")

expMat_all_Assay2_2<- as.matrix(expMat_all_Assay2)

colnames(cellMeta_all2) <- NULL

genePos.all2<- genePos.all %>% select(-c(width, symbol))
colnames(genePos.all2) <- NULL

seRNA.GBM.2 <- readRDS("~/seRNA.GBM.2.rds")
library(SummarizedExperiment)


#try running with just keep clusters
cells_interest <- c("Neo_1", "Trans_1", "Trans_2", "Neo_2", "Trans_3", "Trans_4", "ASC_1", "ASC_2", "ASC_3", "ASC_4", "Micro_1", "Macro_1")
cellMeta_all3<- cellMeta_all2 %>% filter(clu05_short_label %in% cells_interest)
seRNA.GBM.2_filt <- seRNA.GBM.2[,rownames(cellMeta_all3)]
seRNA.GBM.2_filt <- seRNA.GBM.2_filt@assays@data$counts

inCNV_obj_keep = CreateInfercnvObject(raw_counts_matrix=seRNA.GBM.2_filt, 
                                      annotations_file=cellMeta_all3,
                                      delim="\t",
                                      gene_order_file=genePos.all2,
                                      ref_group_names=c("Micro_1", "Macro_1")) #reference cell pop

inCNV_obj_keep = infercnv::run(inCNV_obj_keep,
                               cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                               out_dir= "inCNVdir_keep",  # dir is auto-created for storing outputs
                               cluster_by_groups=T,   # cluster
                               denoise=T, noise_logistic=TRUE, # turns gradient filtering on 
                               HMM=T, debug = TRUE, diagnostics = TRUE, HMM_type="i3")



infercnv_obj_ASC = readRDS('inCNVdir_keep/preliminary.infercnv_obj')
saveRDS(infercnv_obj_ASC, file = "infercnv_obj_ASC.rds")



#CNV data
CNV_df <- as.data.frame(infercnv_obj_ASC@expr.data)

#1)subset matrix into sub-matrices for cells from each tissue sample
tumor_1<- CNV_df[grepl("tumor_1",colnames(CNV_df))]
tumor_2<- CNV_df[grepl("tumor_2",colnames(CNV_df))]
tumor_3<- CNV_df[grepl("tumor_3",colnames(CNV_df))]

margin_1<- CNV_df[grepl("margin_1",colnames(CNV_df))]
margin_2<- CNV_df[grepl("margin_2",colnames(CNV_df))]


#2)Calculate tissue average (row mean)
mean_tumor_1<- rowMeans(tumor_1)
mean_tumor_2<- rowMeans(tumor_2)
mean_tumor_3<- rowMeans(tumor_3)
mean_margin_1<- rowMeans(margin_1)
mean_margin_2<- rowMeans(margin_2)

#3)lappy, for every column (cell) perform correlation with the above tissue average and 
#extract correlation coefficient
correlation_coefficients_tumor1 <- numeric(length(tumor_1))
correlation_coefficients_tumor1 <- apply(tumor_1, 2, function(col) cor(col, mean_tumor_1))
cor_tumor1 <- as.data.frame(correlation_coefficients_tumor1)

correlation_coefficients_tumor2 <- numeric(length(tumor_2))
correlation_coefficients_tumor2 <- apply(tumor_2, 2, function(col) cor(col, mean_tumor_2))
cor_tumor2 <- as.data.frame(correlation_coefficients_tumor2)

correlation_coefficients_tumor3 <- numeric(length(tumor_3))
correlation_coefficients_tumor3 <- apply(tumor_3, 2, function(col) cor(col, mean_tumor_3))
cor_tumor3 <- as.data.frame(correlation_coefficients_tumor3)

correlation_coefficients_margin1 <- numeric(length(margin_1))
correlation_coefficients_margin1 <- apply(margin_1, 2, function(col) cor(col, mean_margin_1))
cor_margin1 <- as.data.frame(correlation_coefficients_margin1)

correlation_coefficients_margin2 <- numeric(length(margin_2))
correlation_coefficients_margin2 <- apply(margin_2, 2, function(col) cor(col, mean_margin_2))
cor_margin2 <- as.data.frame(correlation_coefficients_margin2)

#make a tumor df and a separate margin df
#first combine correlation dfs and normalized CNV signal
comb_tumor1 <- merge(cor_tumor1, scaled_data, by = 'row.names')
colnames(comb_tumor1) <- c("cell", "cor_coef", "mean_of_squares")
comb_tumor2 <- merge(cor_tumor2, scaled_data, by = 'row.names')
colnames(comb_tumor2) <- c("cell", "cor_coef", "mean_of_squares")
comb_tumor3 <- merge(cor_tumor3, scaled_data, by = 'row.names')
colnames(comb_tumor3) <- c("cell", "cor_coef", "mean_of_squares")

comb_tumor <- rbind(comb_tumor1, comb_tumor2, comb_tumor3)

comb_margin1 <- merge(cor_margin1, scaled_data, by = 'row.names')
colnames(comb_margin1) <- c("cell", "cor_coef", "mean_of_squares")
comb_margin2 <- merge(cor_margin2, scaled_data, by = 'row.names')
colnames(comb_margin2) <- c("cell", "cor_coef", "mean_of_squares")

comb_margin <- rbind(comb_margin1, comb_margin2)

#plot distributions
comb_tumor_2 <- merge(comb_tumor, cellMeta_all2%>% select(clu05_short_label), by.x = "cell", by.y = 'row.names')
comb_margin_2 <- merge(comb_margin, cellMeta_all2%>% select(clu05_short_label), by.x = "cell", by.y = 'row.names')
library(ggplot2)
ggplot(comb_tumor_2, aes(x=clu05_short_label, y=cor_coef, fill=clu05_short_label)) + #geom_hline(yintercept=0.370, linetype="dashed", color = "red")+
  geom_violin() + theme(axis.text.x = element_text(angle = 90)) + 
  geom_jitter()
ggplot(comb_margin_2, aes(x=clu05_short_label, y=cor_coef, fill=clu05_short_label)) + #geom_hline(yintercept=0.370, linetype="dashed", color = "red")+
  geom_violin() + theme(axis.text.x = element_text(angle = 90)) + 
  geom_jitter()


#set 0.4 as a correlation cutoff. 
#if a value is larger for tumor dataset, then yes (neoplastic), if less than, then no
#for margin dataset it will be the opposite rule
comb_tumor$malignancy <- ifelse(comb_tumor$cor_coef >= 0.4, "yes", "no")
comb_margin$malignancy <- ifelse(comb_margin$cor_coef >= 0.4, "no", "yes")
#combine into one dataframe
comb_all<- rbind(comb_margin, comb_tumor)
comb_final <- comb_all[,c(1,2,4)]
#export so it can be added as metadata to sn-multiome

###################
#add to UMAP

#load UMAP coordinates
CNV_coord<- read.table(file = "C:/Users/caitl/Documents/third rotation/inferCNV/UMAP3_position .txt")
CNV_label<- read.table(file = "C:/Users/caitl/Documents/third rotation/inferCNV/cell_data_CNV.txt")

#attach column with "yes", "no"
CNV_coord2 <- merge(CNV_coord, CNV_label, by = "row.names")
#make plot
library(ggplot2)

# Create a scatter plot with points colored by the "Label" variable
ggplot(CNV_coord2, aes(x = CNV_coord2$allcells.combined.Ham.UMAP_Dimension_1, y = CNV_coord2$allcells.combined.Ham.UMAP_Dimension_2, color = CNV_coord2$CNV2)) +
  geom_point() +
  scale_color_manual(values = c("yes" = "#306fab", "no" = "#9f9f9f", "missing" = "#d7d7d7")) +
  ggtitle("CNV score")+
  theme_minimal()