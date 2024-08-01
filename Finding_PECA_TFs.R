setwd("~/third rotation/Organoid maturation timeline/HiSeq/TimeReg/binned_PECA_comparisons")
#load TF-TG score output from PECA comparisons
early_late <- read.table("early_late/early_specific_module.txt", header = TRUE)
late_early <- read.table("early_late/late_specific_module.txt", header = TRUE)
early_middle <- read.table("early_middle/early_specific_module.txt", header = TRUE)
middle_early <- read.table("early_middle/middle_specific_module.txt", header = TRUE)
middle_late <- read.table("middle_late/middle_specific_module.txt", header = TRUE)
late_middle <- read.table("middle_late/late_specific_module.txt", header = TRUE)

#subset to only include those with a substantial FC 
early_late_TFs<- as.data.frame(table(early_late$TF))
early_late_TFs <- early_late_TFs[order(early_late_TFs$Freq, decreasing = TRUE),]
top_early_late_TFs <- early_late_TFs[1:5, ]
early_late_top <- early_late %>% filter(Fold > 1.5)#683/5386

late_early_TFs<- as.data.frame(table(late_early$TF))
late_early_TFs <- late_early_TFs[order(late_early_TFs$Freq, decreasing = TRUE),]
top_late_early_TFs <- late_early_TFs[1:10, ]
late_early_top <- late_early %>% filter(Fold > 1.5) #169/3096

early_middle_TFs<- as.data.frame(table(early_middle$TF))
early_middle_TFs <- early_middle_TFs[order(early_middle_TFs$Freq, decreasing = TRUE),]
top_early_middle_TFs <- early_middle_TFs[1:5, ]
early_middle_top <- early_middle %>% filter(Fold > 1.5) #325/5417

middle_early_TFs<- as.data.frame(table(middle_early$TF))
middle_early_TFs <- middle_early_TFs[order(middle_early_TFs$Freq, decreasing = TRUE),]
top_middle_early_TFs <- middle_early_TFs[1:5, ]
middle_early_top <- middle_early %>% filter(Fold > 1.5)#435/7360

middle_late_TFs<- as.data.frame(table(middle_late$TF))
middle_late_TFs <- middle_late_TFs[order(middle_late_TFs$Freq, decreasing = TRUE),]
top_middle_late_TFs <- middle_late_TFs[1:5, ]
middle_late_top <- middle_late %>% filter(Fold > 1.5) #186/5984

late_middle_TFs<- as.data.frame(table(late_middle$TF))
late_middle_TFs <- late_middle_TFs[order(late_middle_TFs$Freq, decreasing = TRUE),]
top_late_middle_TFs <- late_middle_TFs[1:10, ]
late_middle_top <- late_middle %>% filter(Fold > 1.5)#130/2512

###note: tried adjusting FC to 1.5 and 1; 1.5 didn't make a huge difference in # of TGs, but lowering to 1 did
#for example RORA went from 2 to 3 TGs when going from FC of 2 to FC of 1.5, but when going to FC of 1, there are
#now ~100 TGs

#make two dfs: (1) for all top TF-TG pairs and (2) a comprehensive list of all TFs in the top TF-TG df
all_top <- rbind(early_late_top, late_early_top, early_middle_top, middle_early_top, middle_late_top, late_middle_top)

#step 1, with df1 ("all_top"), add a column indicating which module the TG is in from the organoid
#WGCNA RNA-seq analysis, then (step 2) group by TF and look at how many TGs for each are in each module. 
load("C:/Users/caitl/Documents/third rotation/Organoid maturation timeline/HiSeq/RNA-seq/WGCNA_FINAL.RData")

submod<-read.csv("C:/Users/caitl/Documents/third rotation/Organoid maturation timeline/HiSeq/RNA-seq/org_RNAseq_gene_modules.csv")

#1
fancy_top <- all_top %>%
  left_join(submod, by= c("TG"="gene_id"))
write.csv(fancy_top, "top_TF_TG.csv", col.names = TRUE)

fancy_top <- fancy_top[!is.na(fancy_top$colors),]

#2
TF_mod_representation<- fancy_top %>%
  group_by(TF) %>%
  dplyr::count(colors) %>%
  mutate(module_freq = n/sum(n), total_TG = sum(n)) 

TF_mod_representation <- TF_mod_representation %>% filter(total_TG > 5) #must have at least 5 TGs;  FC>1.5
TF_mod_representation <- TF_mod_representation[!is.na(TF_mod_representation$colors), ]
TF_mod_representation <- TF_mod_representation[order(TF_mod_representation$total_TG, decreasing = TRUE),]
TFs <- as.data.frame(unique(TF_mod_representation$TF))#Final TFs
colnames(TFs) <- "TF"
write.csv(TFs, "TFs_from_PECA_FC1.5_FINAL.csv", col.names = TRUE)

TF_mod_representation$TF <- factor(TF_mod_representation$TF)
TF_mod_representation$colors <- factor(TF_mod_representation$colors)
#add a column called "Label" which combines the TF name with the total number of TGs in parentheses
#this will serve as the title for each pie chart
TF_mod_representation$Label<- paste0(TF_mod_representation$TF, " ","(",TF_mod_representation$total_TG, ")")

#make a pie chart for each TF showing the representation of each module in the respective TGs
#modified code from: <https://www.geeksforgeeks.org/create-multiple-pie-charts-using-ggplot2-in-r/#:~:text=For%20building%20a%20Pie%20Chart%20in%20R%2C%20we,to%20use%20an%20additional%20method%20named%20facet_grid%20%28%29>
library(ggplot2)
ggplot(data=TF_mod_representation, aes(x=" ", y=module_freq, group=colors, colour=colors, fill=colors)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) + 
  facet_wrap(factor(Label, levels = unique(Label))~.) +theme_void()

#try making stacked bar chart that could go next to heatmap
#load the heatmap row order and TF labels for bar chart order
load("C:/Users/caitl/Documents/third rotation/Organoid maturation timeline/HiSeq/TimeReg/pecaTFs_heat_info.RData")
pecaTFs_heat_rlabels <- as.data.frame(pecaTFs_heat_rlabels)
pecaTFs_heat_rorder<- as.vector(pecaTFs_heat_rorder)
pecaTFs_heat_rlabels_clean<- gsub("\\_.*", "", pecaTFs_heat_rlabels$pecaTFs_heat_rlabels)
pecaTFs_heat_rlabels_clean<- as.data.frame(pecaTFs_heat_rlabels_clean)
pecaTFs_heat_rlabels_2 <- cbind(pecaTFs_heat_rlabels,pecaTFs_heat_rlabels_clean)
pecaTFs_heat_rlabels_2 <- pecaTFs_heat_rlabels_2[pecaTFs_heat_rorder,] #now we have a cleaned up list of our TFs in the same order they were on the heatmap
pecaTFs_heat_rlabels_2 <- as.data.frame(pecaTFs_heat_rlabels_2)
colnames(pecaTFs_heat_rlabels_2) <- c("TF_encode", "TF")
pecaTFs_heat_rlabels_total<- pecaTFs_heat_rlabels_2 %>%
  group_by(TF_encode) %>%
  left_join(TF_mod_representation, by= c("TF"="TF"))

pecaTFs_heat_rlabels_total$TF_encode <- factor(pecaTFs_heat_rlabels_total$TF_encode, levels = rev(pecaTFs_heat_rlabels_2$TF_encode))

ggplot(pecaTFs_heat_rlabels_total, aes(x=TF_encode, y=module_freq, group=colors, colour=colors, fill=colors)) +
  geom_bar(stat="identity", width = .7, colour="black", lwd=0.1) +
  coord_flip() +
  labs(y="", x="")

