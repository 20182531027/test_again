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
library(clusterProfiler)
# BiocManager::install("BiocNeighbors")
# devtools::install_github("sqjin/CellChat")
library(CellChat)

# Ensure that when beginning, the directory only has "GSE237915_RAW.tar" file, and depressed this file, then continue the next pipeline.

######################1.prepare the correct format ###########################################
# set work path that store the deprssed "GSE237915_RAW.tar"
setwd("C:/Users/25988/Desktop/GSE237915_RAW") # the directory is where "GSE237915_RAW.tar" exists.
## List the files starting with GSM and ending with gz in the current folder
files<-list.files('./', pattern = '^GSM.*\\.gz$')

## Obtain grouping information for 8 samples, i.e. the relationship between sample GSM number and grouping

df <- data.frame(sample_id = character(), condition = character(), stringsAsFactors = FALSE)
for(i in files){
  split_name <- strsplit(i, "__")[[1]][1]
  # Extract "sample_id" and "condition" information
  sample_id <- str_split(split_name,"_",n=2)[[1]][1]
  condition <- str_split(split_name, "_",n=2)[[1]][2]
  # Add information to the data.frame
  df <- rbind(df, data.frame(sample_id = sample_id, condition = condition, stringsAsFactors = FALSE))
}
df<-unique(df)
df$replicate<-c("1","2","3","1","2","1","1","2")

## Organize each sample and place the "barcodes. tsv. gz", "features. tsv. gz", and "matrix. mtx. gz" of each sample into the sample folder
fs = list.files('./',pattern = '^GSM.*\\.gz$')
samples <- substr(fs,1,10)
lapply( unique(samples), function(x){
  y = fs[grepl(x,fs)]
  folder = paste0('./',strsplit(y[1],split = '_')[[1]][1])
  #creat folder
  dir.create(folder,recursive = T)
  #Rename the subfolders and move them to the corresponding folders
  file.rename(paste0('./',y[1]),file.path(folder,"barcodes.tsv.gz"))
  file.rename(paste0('./',y[2]),file.path(folder,"features.tsv.gz"))
  file.rename(paste0('./',y[3]),file.path(folder,"matrix.mtx.gz"))
})



#############################2.read and Quality control sample by loop #################################################
# due to the more samples take in many storage, I just choose 2 healthy samples, 2 severe symptom samples
# read samples by a loop
df<-df[c(1,3,7,8),]
write.csv(df,"0.sample_info.csv")
sample_list <- list()
plot_list <- list()
for(i in 1:nrow(df)){
  sample_name <- paste(df[i,2],df[i,3],sep="_")
  sample.data <- Read10X(data.dir = df[i,1],gene.column = 2) # you need make some changes
  # Initialize the Seurat object with the raw (non-normalized data).
  sampleRDS <- CreateSeuratObject(counts = sample.data, project = sample_name, min.cells = 3, min.features = 200)
  sampleRDS[["percent.mt"]] <- PercentageFeatureSet(sampleRDS, pattern = "^mt-")
  sampleRDS$Group <- df[i,2]
  p1 <- VlnPlot(sampleRDS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  plot_list[[sample_name]] <- p1
  # filter low quality cells
  sampleRDS <- subset(sampleRDS, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 5)
  sample_list[[sample_name]] <- sampleRDS
}
py<-plot_grid(plotlist = plot_list, nrow = 2, ncol = 2,labels = "AUTO")
ggsave("1.all_sample_QC_violin.pdf",py,width=10,height=8,dpi=300)

## merge sample
hxj <- length(sample_list)
myrds <- merge(sample_list[[1]], y = sample_list[c(2:hxj)])
saveRDS(myrds,"raw_merged_all_sample.RDS")




##############################3.Batch Effect remove and dimention reduction #######################################################
myrds <- NormalizeData(myrds, normalization.method = "LogNormalize", scale.factor = 1e4) 
myrds <- FindVariableFeatures(myrds, selection.method = 'vst', nfeatures = 2000)
"jev-p3" %in% VariableFeatures(myrds)
## remove the gene of the JEV to avoid the influence for cell cluster
VariableFeatures(myrds) <- VariableFeatures(myrds)[!grepl("jev-p3",VariableFeatures(myrds))]
scRNA.T <- ScaleData(myrds, rownames(myrds))

"jev-p3" %in% VariableFeatures(object = scRNA.T)
scRNA.T <- RunPCA(scRNA.T, features = VariableFeatures(object = scRNA.T)) 
ElbowPlot(scRNA.T)
DimPlot(scRNA.T, reduction = "pca")
scRNA.T <- RunUMAP(scRNA.T, dims = 1:15, reduction = "pca")
### Batch correction not carried out
scRNA.T <- RunUMAP(scRNA.T,reduction = "pca", dims = 1:15, reduction.name = "umap_naive")
p2 <- DimPlot(scRNA.T, reduction = "umap_naive",group.by = "orig.ident")


scRNA.T <- RunHarmony(scRNA.T, "orig.ident", assay.use = "RNA", max.iter.harmony = 20,reduction.save = "harmony")
scRNA.T <- RunUMAP(scRNA.T, dims = 1:15, reduction = "harmony")
p3 <- DimPlot(scRNA.T, reduction = "harmony",group.by = "orig.ident")+ggtitle("harmoy_batch_remove")
p <- plot_grid(p2+NoLegend(), p3, labels = "AUTO")
ggsave("2.batch_effect_before_and_after_harmony.pdf",p,width=10,height = 5,dpi=300)




###############################4. Find cluster and cell annotation ########################################################
ElbowPlot(scRNA.T)
scRNA.T <- FindNeighbors(scRNA.T, reduction = "harmony", dims = 1:15)
scRNA.T <- FindClusters(scRNA.T, resolution = 0.05, algorithm = 1)
head(scRNA.T@meta.data)
p4 <- DimPlot(scRNA.T, reduction = "umap", group.by = c("orig.ident","seurat_clusters"), label = TRUE, raster = FALSE)+NoLegend()
ggsave("3.0cluster_sample_umap.pdf",p4,width = 12,height =5,dpi=300)
Idents(scRNA.T)<-"seurat_clusters"
DefaultAssay(scRNA.T)
scRNA.T<-JoinLayers(scRNA.T)
table(scRNA.T$seurat_clusters)
ob1 <- subset(x = scRNA.T, downsample = 100)
test_0.5 <- FindAllMarkers(scRNA.T,logfc.threshold = 0.25, only.pos = TRUE)
test_0.5 <- test_0.5 %>% arrange(cluster, desc(avg_log2FC))
write.csv(test_0.5,"cluster_markers_gene.csv")
test_0.5  %>% group_by(cluster) %>%top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(ob1,features = top10$gene,label = T)

p5 <- DotPlot(scRNA.T,features=unique(top10$gene))+RotatedAxis()
ggsave("cluster_marker_gene_exp_dotplot.pdf",p5,width = 20,height = 8,dpi=300)

p6 <- DotPlot(scRNA.T,features=c("Ly6c2",
                                          "Cx3cr1",
                                          "Snap25",
                                          "Nkg7","Cd3e",
                                          "Plp1",
                                          "Slc1a2","Olig1","Olig2",
                                          "Cldn5","Ly6c1",
                                          "Penk",
                                          "Vtn",
                                          "Tubb2b","Nnat",
                                          "Hbb-bt"))+RotatedAxis()

ggsave("celltypes_marker_gene_exp_dotplot.pdf",p6,width = 20,height = 8,dpi=300)

levels(scRNA.T)
new.cluster.ids <- c("Mon_Mac","Mic_Mac","Neurons","NK_T","OLG","Ast_OPCs","End","Neurons","Mur","Ependymal_cell","erythrocyte")
names(new.cluster.ids) <- levels(scRNA.T)
scRNA.T <- RenameIdents(scRNA.T, new.cluster.ids)

head(scRNA.T@meta.data)
DimPlot(scRNA.T, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

cluster_labels <- c(0,1,2,3,4,5,6,7,8,9,10)
cell_types <- c("Mon_Mac","Mic_Mac","Neurons","NK_T","OLG","Ast_OPCs","End","Neurons","Mur","Ependymal_cell","erythrocyte")
qiqi<-setNames(cell_types, cluster_labels)
scRNA.T@meta.data$celltypes <- qiqi[scRNA.T@meta.data$seurat_clusters]
p7 <- DimPlot(scRNA.T, reduction = "umap", group.by = c("orig.ident","celltypes"), label = TRUE, raster = FALSE,label.size = 5)+NoLegend()
ggsave("3.celltypes_annotated_umap.pdf",p7,width = 12,height =5,dpi=300)
p8 <- DimPlot(scRNA.T, reduction = "umap", group.by = c("seurat_clusters","celltypes"), label = TRUE, raster = FALSE,label.size = 5)+NoLegend()
ggsave("3.1cluster_celltypes_annotated_umap.pdf",p8,width = 12,height =5,dpi=300)
saveRDS(scRNA.T,"scaled_annotated_jev_scRNA_seq.RDS")

p8.1 <- FeaturePlot(scRNA.T, c("Ly6c2",
                    "Cx3cr1",
                    "Snap25",
                    "Nkg7","Cd3e",
                    "Plp1",
                    "Slc1a2","Olig1","Olig2",
                    "Cldn5","Ly6c1",
                    "Penk",
                    "Vtn",
                    "Tubb2b","Nnat",
                    "Hbb-bt","jev-p3"),ncol =6,reduction = "umap")
ggsave("3.2.important_gene_exp_umap.pdf",p8.1,width = 20,height =9,dpi=300)




################################5. Cell ratio variation among Groups#############################################################

DimPlot(scRNA.T, reduction = "umap", group.by = c("celltypes"), label = TRUE, raster = FALSE,label.size = 5,split.by = "Group")+NoLegend()

cluster_source_counts <- table(scRNA.T@meta.data$celltypes,scRNA.T@meta.data$Group)
cluster_source_counts_df <- as.data.frame(cluster_source_counts)
colnames(cluster_source_counts_df) <- c("celltypes", "Sample_Source", "Count")
cluster_source_proportions <- prop.table(cluster_source_counts, margin = 2)
Cellratio <- as.data.frame(cluster_source_proportions)
cluster_source_proportions <- prop.table(cluster_source_counts, margin = 2)
Cellratio <- as.data.frame(cluster_source_proportions)
colnames(Cellratio) <- c("celltypes", "Sample_Source", "Freq")
Cellratio$sample_source<-factor(Cellratio$Sample_Source,levels=c("Healthy","severe_symptom"))
p9 <- ggplot(Cellratio) + 
  geom_bar(aes(x =sample_source, y= Freq, fill = celltypes),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='sample',y = 'Ratio')+
  # scale_fill_manual(values = cols)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=12))+
  theme(axis.text.y = element_text(size=12))

ggsave("4.celltypes_ratio_variation_barplot.pdf",p9,width = 4,height =5,dpi=300)



################################6. Gene expression Variation among Groups && Pathway enrichment analysis##################################################################
Idents(scRNA.T)<-"Group"
table(scRNA.T$Group)
ob1 <- subset(x = scRNA.T, downsample = 1000)
test_df <- FindAllMarkers(ob1,logfc.threshold = 0.25, only.pos = TRUE)
# test_cluster_markers <- FindMarkers(ob1, ident.1 = "Healthy", logfc.threshold = 0.25, test.use = "roc")
test_df <- test_df %>% arrange(cluster, desc(avg_log2FC))
test_df  %>% group_by(cluster) %>%top_n(n = 20, wt = avg_log2FC) -> top20
p10 <- DoHeatmap(ob1,features = top20$gene,label = F)
ggsave("5.group_DEG_heatmap.pdf",p10,width = 8,height =5,dpi=300)
write.csv(test_df,"group_DEGs.csv")

# test_df<-read.csv("../group_DEGs.csv",row.names = 1)

GO_database <- 'org.Mm.eg.db' #GOAnalyze the specified species, please refer to the species abbreviation index table for details http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
#BiocManager::install("org.Mm.eg.db")
#KEGG_database <- 'mmu' #KEGG analysis specifies species, please refer to the species abbreviation index table for details http://www.genome.jp/kegg/catalog/org_list.html

#Set up enriched genes (severe_symptom)
gene<-test_df[["gene"]][test_df$cluster=="severe_symptom"]
gene_C <- bitr(gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database) #To change from another name to symbol, the decimal point needs to be removed


#GO analysis
GO_up<-enrichGO( gene_C$ENTREZID,
                 OrgDb = GO_database,
                 keyType = "ENTREZID",#Set the type of gene ID to be read
                 ont = "BP",#(ont = "ALL", therefore it includes three parts: Biological Process, Cellular Component, and Molecular Function）
                 pvalueCutoff = 0.05,#Set p-value threshold
                 qvalueCutoff = 0.05,#Set q-value threshold
                 readable = T)

#Set up enriched genes (Healthy)
gene<-test_df[["gene"]][test_df$cluster=="Healthy"]
gene_C <- bitr(gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)#To change from another name to symbol, the decimal point needs to be removed

#GO analysis
GO_down<-enrichGO( gene_C$ENTREZID,
                   OrgDb = GO_database,
                   keyType = "ENTREZID",#Set the type of gene ID to be read
                   ont = "BP",#(ont = "ALL", therefore it includes three parts: Biological Process, Cellular Component, and Molecular Function）
                   pvalueCutoff = 0.05,#Set p-value threshold
                   qvalueCutoff = 0.05,#Set q-value threshold
                   readable = T)


#Extract data.frame
df_up=GO_up@result
df_down=GO_down@result
top="p.adj"
#Extract data
if(top == "count"|top =="Count"){
  df_up=df_up[order(df_up$Count,decreasing = T),][1:20,]
  df_down=df_down[order(df_down$Count,decreasing = T),][1:20,]
}else if(top == "p.adj"){
  df_up=df_up[order(df_up$p.adjust),][1:20,]
  df_down=df_down[order(df_down$p.adjust),][1:20,]
}else if(top == "random"){#random Self designated pathway
  df_up=df_up[order(df_up$p.adjust),][c(1,10,21,28,29,33,35,36,48,11,43,60,154,189,439,52,56,62,99,654),]
  df_down=df_down[order(df_down$p.adjust),][1:20,]
}else{
  print("inpute correct 'orderby' parameter")
  
}

#Add label
df_up$type="up"
df_down$type="down"
#Merge two data.frame
df_GO=rbind(df_up,df_down)

## merely draw upregulated gene pathway
df_GO$`-Log(q-value)` <- -log(df_GO$qvalue)


test = df_GO %>% filter(type=="up") 
test_sorted <- test %>% arrange((`-Log(q-value)`))
test_sorted$Description<-factor(test_sorted$Description,levels = test_sorted$Description)
bp <- ggplot() +
  geom_bar(data = test_sorted,
           aes(x = `-Log(q-value)`, y = Description),
           width= 0.5, #Width adjustment
           stat= 'identity',fill="#ffc55b") +theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # panel.border = element_blank(),
    axis.title.x = element_text(size = 15),  # Set the font size of the X-axis title
    axis.title.y = element_text(size = 15),  # Set the font size of the Y-axis title
    axis.text.x = element_text(size = 12),   # Set the font size of the X-axis scale label
    axis.text.y = element_text(size = 12)    # Set the font size of the Y-axis scale label
  )+ggtitle("20 biological pathway upregulated in severe symptom group")+scale_x_continuous(expand = c(0,0))
ggsave("Fig.3c_upregulated_pathway_in_severe_symptom_group.pdf",width = 14,height = 5,dpi = 300)

test = df_GO %>% filter(type=="down") 
test_sorted <- test %>% arrange((`-Log(q-value)`))
test_sorted$Description<-factor(test_sorted$Description,levels = test_sorted$Description)
bp <- ggplot() +
  geom_bar(data = test_sorted,
           aes(x = `-Log(q-value)`, y = Description),
           width= 0.5, #Width adjustment
           stat= 'identity',fill="#ffc55b") +theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        # panel.border = element_blank(),
        axis.title.x = element_text(size = 15),  # Set the font size of the X-axis title
        axis.title.y = element_text(size = 15),  # Set the font size of the Y-axis title
        axis.text.x = element_text(size = 12),   # Set the font size of the X-axis scale label
        axis.text.y = element_text(size = 12)    # Set the font size of the Y-axis scale label
  )+ggtitle("20 biological pathway upregulated in severe symptom group")+scale_x_continuous(expand = c(0,0))
ggsave("Fig.3c_downregulated_pathway_in_severe_symptom_group.pdf",width = 14,height = 5,dpi = 300)



#Processing data.frame
# Create a function to convert fractional form to decimal form
convert_fraction_to_decimal <- function(fraction) {
  parts <- strsplit(fraction, "/")[[1]]  # Split a string into numerator and denominator
  numerator <- as.numeric(parts[1])      # Extract numerator and convert them into numbers
  denominator <- as.numeric(parts[2])    # Extract denominator and convert them into numbers
  decimal <- numerator / denominator      # Calculate the decimal value
  return(decimal)
}
# use GeneRatio function on the column to convert fractional form to decimal point form
df_GO$GeneRatio <- sapply(df_GO$GeneRatio, convert_fraction_to_decimal)
#Adjust the down ratio to a negative number for easier drawing
df_GO$GeneRatio <- ifelse(df_GO$type == "down", -df_GO$GeneRatio, df_GO$GeneRatio)
# Simplify the number represented by Scientific notation to the number before two decimal places, which is convenient for drawing
df_GO$p.adjust <- sprintf("%.2e", df_GO$p.adjust)


#Draw and display the up and down adjustment results simultaneously
library(tidyverse)

ggplot(df_GO,aes(reorder(Description, GeneRatio),GeneRatio,fill=type))+
  geom_col()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        # legend.title = element_blank(),
        axis.text = element_text(color="black",size=10),
        axis.line.x = element_line(color='black'),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size = 20, hjust = 0.5)
        # legend.position = 'none'
  )+
  coord_flip()+
  #Draw the middle dividing line
  geom_segment(aes(y=0, yend=0,x=0,xend=40.5))+
  #Add GO note location
  geom_text(data = df_GO[which(df_GO$GeneRatio > 0),], aes(x = Description, y = -0.005,
                                                           label = Description),  hjust = 1, size = 4) +
  geom_text(data = df_GO[which(df_GO$GeneRatio < 0),], aes(x = Description, y = 0.005, 
                                                           label = Description),  hjust = 0, size = 4) +
  #Add Padj position
  geom_text(data = df_GO[which(df_GO$GeneRatio>0),],aes(label=p.adjust),
            hjust=-0.1, size=4, color='red')+
  geom_text(data = df_GO[which(df_GO$GeneRatio<0),],aes(label=p.adjust),
            hjust=1.1, size=4, color="red")+
  scale_fill_manual(values = c("#1084A4","#8D4873"))+
  #Control the range of x and y axes
  scale_x_discrete(expand = expansion(mult = c(0,0)))+
  ##Adjust the scope according to the actual situation
  ylim(-0.15, 0.15)+
  labs(x='', y='GeneRatio')+
  ggtitle("severe_symptom vs Healthy")

ggsave("5.1DEGs_enrich_GO_pathway.pdf",width =15,height = 9,dpi = 300,limitsize = FALSE)





###################################7. cell cell communication analysis ##################################
rm(sampleRDS,sample.data,sample_list)
scRNA.T <- readRDS("scaled_annotated_jev_scRNA_seq.RDS")
dir.create("cellchat")
setwd("cellchat")
sample_list<-list()
unique(scRNA.T$Group)
for (i in unique(scRNA.T$Group)){
  samplerds <- subset(scRNA.T,Group==i)
  # saveRDS(samplerds,paste(i,"jev_scRNAseq.RDS",sep="_"))
  data.input <- GetAssayData(samplerds, layer = 'data')
  meta <- samplerds@meta.data[,c("orig.ident","celltypes")]
  colnames(meta) <-  c("group","celltypes")
  table(meta$celltypes)
  # meta$celltypes <- gsub(" cells|-cells", "", meta$celltypes)
  table(meta$celltypes)
  # Adipocytes           B      CD4+ T      CD8+ T Endothelial Fibroblasts 
  #        244         667        1007         547          28         762 
  #  Monocytes Neutrophils          NK 
  #        311         378          87 
  # Suggest sorting celltypes in advance~ 
  identical(rownames(meta),colnames(data.input))
  # celltype_order <- c("CD4+ T","CD8+ T","B","Adipocytes","Endothelial",
  #                     "Fibroblasts","Monocytes","Neutrophils","NK")
  meta$celltypes <- as.factor(meta$celltypes)
  # Sort according to the order of meta$celltypes
  ordered_indices <- order(meta$celltypes)
  # Sort "meta" and "data.input"
  meta <- meta[ordered_indices, ]
  data.input <- data.input[, ordered_indices]
  identical(rownames(meta),colnames(data.input))
  
  # creat cellchat object
  cellchat <- createCellChat(object = data.input, 
                             meta = meta, 
                             group.by = "celltypes")
  CellChatDB<-CellChatDB.mouse
  CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact",
                                                    "Non-protein Signaling"), key = "annotation") # use Secreted Signaling
  # set the used database in the object
  cellchat@DB <- CellChatDB.use
  #Preprocessing the expression data for cell-cell communication analysis
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  #future::plan("multiprocess", workers = 4) # do parallel
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.mouse)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  saveRDS(cellchat,paste(i,"jev_scRNAseq_cellchat.RDS",sep="_"))
  sample_list[[i]] <- cellchat
}



sample_list <- list()
## Healthy group
j = "Healthy"
cellchat <-readRDS("Healthy_jev_scRNAseq_cellchat.RDS")
nPatterns = 2
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
p1 <- netAnalysis_river(cellchat, pattern = c("outgoing"))
p2 <- netAnalysis_river(cellchat, pattern = "incoming")
p3 <- p1|p2
ggsave(paste("Fig.4c_",j,"_group_cellchat_river_plot.pdf",sep=""),p3,width =15,height = 12,dpi = 300,limitsize = FALSE)
sample_list[[j]] <- cellchat

## severe symptom
j="severe_symptom"
cellchat <-readRDS("severe_symptom_jev_scRNAseq_cellchat.RDS")
nPatterns = 2
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
p1 <- netAnalysis_river(cellchat, pattern = c("outgoing"))
p2 <- netAnalysis_river(cellchat, pattern = "incoming")
p3 <- p1|p2
ggsave(paste("Fig.4d_",j,"_group_cellchat_river_plot.pdf",sep = ""),p3,width =15,height = 12,dpi = 300,limitsize = FALSE)

sample_list[[j]] <- cellchat

for (i in names(sample_list)){
  cellchat <- sample_list[[i]]
  pdf(paste(i,"circle_plot.pdf",sep = "_"),width = 12,height = 5)
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  dev.off()
  
}

## merged cellchat result
levels(sample_list[[1]]@meta$celltypes)
levels(sample_list[[2]]@meta$celltypes)
group.new = levels(sample_list[[2]]@idents)
sample_list[[1]] <- liftCellChat(sample_list[[1]], group.new)
## Integrate cellchat result to compare between groups
cellchat <- mergeCellChat(sample_list, add.names = names(sample_list),cell.prefix = T)
saveRDS(cellchat,"merged_sample_cellchat.RDS")


png("spweight.png",width = 30,height = 15,units = "cm",res = 200)
weight.max <- getMaxWeight(sample_list, attribute = c("idents","weight"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(sample_list)) {
  netVisual_circle(sample_list[[i]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 6, title.name = paste0("Weight of interactions - ", names(sample_list)[i]))
}
dev.off()

png("spnumber.png",width = 30,height = 15,units = "cm",res = 200)
weight.max <- getMaxWeight(sample_list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(sample_list)) {
  netVisual_circle(sample_list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 6, title.name = paste0("Number of interactions - ", names(sample_list)[i]))
}
dev.off()


gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))+
  # theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 16))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")+
  # theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 16))
gg1 + gg2
ggsave("Fig.4a_b_comparasion_bar_plot_for_number_weight.pdf",gg1 + gg2,width = 10,height = 6,dpi = 300)


pdf("Fig.4a_b_comparasion_shell_plot_for_cell_cell_interaction.pdf", width=10,height = 6)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

pe_0 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
pe <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
ggsave("Fig.4e_comparasion_cell_cell_communication_pathway.pdf", pe_0+pe, width=10,height =15,dpi=300)


## look the interaction for Neurons
png("Fig.4i_j_Neurons+Cells+recevier.png", width = 30,height = 20,units = "cm",res = 200)
par(mfrow = c(1,2))
for (i in 1:length(sample_list)){
  cellchat1=sample_list[[i]]
  groupSize=as.numeric(table(cellchat1@idents))
  name=names(sample_list)[[i]]
  netVisual_circle(cellchat1@net$weight, vertex.weight = groupSize,
                   weight.scale = T, label.edge= F,
                   #title.name = "Interaction number",
                   targets.use = c('Neurons'),
                   title.name = paste0("Neurons","_as_receiver in ",name))
}
dev.off()


png("Fig.4f_g_Neurons+Cells+sender.png", width = 30,height = 20,units = "cm",res = 200)
par(mfrow = c(1,2))
for (i in 1:length(sample_list)){
  cellchat1=sample_list[[i]]
  groupSize=as.numeric(table(cellchat1@idents))
  name=names(sample_list)[[i]]
  netVisual_circle(cellchat1@net$weight, vertex.weight = groupSize,
                   weight.scale = T, label.edge= F,
                   #title.name = "Interaction number",
                   sources.use = c('Neurons'),
                   title.name = paste0("Neurons","_as_sender in ",name))
}
dev.off()

# ligand-receptor analysis
p1<-netVisual_bubble(cellchat, sources.use = c("Neurons"),  
                     #signaling = c("COLLAGEN","CLDN","GALECTIN","LAMININ","GAP","FN1","PSAP","NCAM","PECAM1","PDGF","VEGF","ADGRG","NOTCH","FGF","PTN"),
                     comparison = c(1:2), angle.x = 45)
ggsave("Fig.4h_Neurons_ligands_receptor_bubble_plot.pdf",p1, width = 5,height = 15,dpi = 300)

p2<-netVisual_bubble(cellchat, targets.use = c("Neurons"),  
                     #signaling = c("COLLAGEN","CLDN","GALECTIN","LAMININ","GAP","FN1","PSAP","NCAM","PECAM1","PDGF","VEGF","ADGRG","NOTCH","FGF","PTN"),
                     comparison = c(1:2), angle.x = 45)
ggsave("Fig.4k_Neurons_ligands_receptor_bubble_plot_receptor.pdf",p2, width = 5,height = 15,dpi = 300)






