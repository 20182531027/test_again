# load packages
rm(list=ls())
library(stringr)
library(gridExtra)
library(Seurat)
library(patchwork)
library(cowplot)
library(ggplot2)
library(harmony)
library(dplyr)

setwd("D:\\course\\Training_data\\Jev_scRNA_seq")
scRNA.T <- readRDS("scaled_annotated_jev_scRNA_seq.RDS")

scRNA.T@meta.data$virus_raw_counts <- scRNA.T@assays[["RNA"]]@layers[["counts"]][row.names(scRNA.T) %in% "jev-p3",match(row.names(scRNA.T@meta.data),colnames(scRNA.T))]
scRNA.T@meta.data$virus_normal <- scRNA.T@assays[["RNA"]]@layers[["data"]][row.names(scRNA.T) %in% "jev-p3",match(row.names(scRNA.T@meta.data),colnames(scRNA.T))]

# draw  virus rate density plot in
p1 <- ggplot(scRNA.T@meta.data, aes(x = virus_normal, colour = celltypes)) +
  geom_density() +
  geom_vline(xintercept = 5, color = "red", linetype = "dashed") + # 在x=5处添加红色虚线
  theme_classic()
ggsave("scaled_virus_reads_distribution_among_celltypes.pdf",p1,width = 7,height = 5,dpi=300)

ggplot(scRNA.T@meta.data, aes(x = virus_normal, colour = Group)) +
  geom_density() +
  geom_vline(xintercept = 5, color = "red", linetype = "dashed") + # 在x=5处添加红色虚线
  theme_classic()

p2 <- ggplot(scRNA.T@meta.data, aes(x = virus_raw_counts, colour = celltypes)) +
  geom_density() +
  xlim(0, 1000) + # 设置x轴的显示范围为0到1000
  geom_vline(xintercept = 500, color = "red", linetype = "dashed") + # 在x=5处添加红色虚线
  theme_classic()
ggsave("raw_virus_counts_distribution_among_celltypes.pdf",p2,width = 7,height = 5,dpi=300)


ggplot(scRNA.T@meta.data, aes(x = virus_raw_counts, colour = Group)) +
  geom_density() +
  xlim(0, 1000) + # 设置x轴的显示范围为0到1000
  geom_vline(xintercept = 500, color = "red", linetype = "dashed") + # 在x=5处添加红色虚线
  theme_classic()
# Because in the paper, JEV positive cell with virus counts bigger than 500, so we set 500 as the threshold

scRNA.T@meta.data$virus_info <- "negative"
scRNA.T@meta.data[scRNA.T@meta.data$virus_raw_counts >= 500,]$virus_info <- "positive"

table(scRNA.T$Group,scRNA.T$virus_info) # in this threshhold, healthy group have no virus infection

cluster_source_counts <- table(scRNA.T@meta.data$virus_info,scRNA.T@meta.data$celltypes)

cluster_source_counts_df <- as.data.frame(cluster_source_counts)
colnames(cluster_source_counts_df) <- c("virus_info", "subtypes", "Count")
cluster_source_proportions <- prop.table(cluster_source_counts, margin = 2)
Cellratio <- as.data.frame(cluster_source_proportions)
cluster_source_proportions <- prop.table(cluster_source_counts, margin = 2)
Cellratio <- as.data.frame(cluster_source_proportions)
colnames(Cellratio) <- c("virus_rate", "subtypes", "Freq")

#Cellratio$celltypes<-factor(Cellratio$celltypes,levels = c("Oligodentricytes", "Meg3_Neurons", "Immune_cells", "Tac2_Neuron", "Oxt_Neurons", "Astrocytes", "Ctxn3_cells"))
#theme(axis.text.x = element_text(angle = 45, hjust = 1))
Cellratio1 <- Cellratio[Cellratio$virus_rate=="positive",]
Cellratio1 <- Cellratio1 %>% arrange(desc(Freq))
Cellratio1$subtypes <- as.character(Cellratio1$subtypes)
Cellratio$subtypes <- as.character(Cellratio$subtypes)
Cellratio$subtypes <- factor(Cellratio$subtypes,levels = Cellratio1$subtypes)

p3 <- ggplot(Cellratio) + 
  geom_bar(aes(x =subtypes, y= Freq, fill = virus_rate),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x = 'sample', y = 'Ratio', axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18)) +
  scale_fill_manual(values = my_palette[2:1],name="virus_info")+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=15))+
  theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20))+
  theme(axis.text.y = element_text(size=15))+
  theme(legend.text = element_text(size = 18), legend.title = element_text(size = 18))
ggsave("Fig.6b.virus_rate_among_celltypes.pdf",p3,width = 7,height = 5,dpi = 300)





