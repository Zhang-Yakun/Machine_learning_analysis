suppressMessages(require(Seurat))
suppressMessages(require(ggplot2))
suppressMessages(require(cowplot))
suppressMessages(require(scater))
suppressMessages(require(scran))
suppressMessages(require(BiocParallel))
suppressMessages(require(BiocNeighbors))
library(dplyr)
library(stringr)
library(data.table)
library(tidyverse)

# Seurat ------------------------------------------------

setwd("E:/Rwork/lung-GEO-8/immunothreapy/")

# 读取数据 --------------------------------------------------------------------

new_counts <- read.table(file="GSE146100_NormData.txt",header = T)
rownames(new_counts)<-new_counts$Gene
new_counts<-new_counts[,-1]

mydata <- CreateSeuratObject(counts = new_counts, min.cells = 3)
# 批量创建seurat对象 

mydata@assays

# 运行标准流程并进行可视化 ------------------------------------------------------------
setwd("e:/Rwork/GSE144945/Integrate_Anchor")

########  VariableFeature

mydata<-FindVariableFeatures(mydata,selection.method="vst",nfeatures=2000)
top10<-head(VariableFeatures(mydata),10)
plot1<-VariableFeaturePlot(mydata)
#散点图，X轴：平均表达值；Y轴：标准化后的方差，显示前2000个高变基因


###### Scale & PCA

mydata <- ScaleData(mydata, verbose = FALSE)
mydata <- RunPCA(mydata, npcs = 30, verbose = FALSE)

VizDimplot<-VizDimLoadings(mydata,dims = 1:2,reduction="pca")

dir.create("PCA")
#ggsave("PCA/PC_1&2_genes.pdf",plot = VizDimplot)

PCAplot<-DimPlot(mydata,reduction = "pca",group.by = "orig.ident")
#可视化PCA聚类
ggsave("PCA/PCA_batch.pdf",plot = PCAplot)

PCA_heatmap<-DimHeatmap(mydata,dims = 1:15,cells = 500,balanced = TRUE)
#PC1:15的基因在500个细胞中表达值的展示
#ggsave("PCA/PC15_heatmap.pdf",plot = PCA_heatmap)


###### clustering

mydata<-FindNeighbors(mydata,dims=1:20)
mydata<-FindClusters(mydata,resolution = 0.5)

head(Idents(mydata),5)
#查看每个细胞对应的聚类数ID

head(mydata@meta.data)
table(mydata@meta.data$seurat_clusters)
#查看每个cluster中有多少细胞，cluster—0中的细胞数最多

mydata<-RunTSNE(mydata,dims=1:20)
mydata<-RunUMAP(mydata,dims = 1:20)

tsne_plot<-DimPlot(mydata,reduction = "tsne",group.by = "orig.ident")
umap_plot<-DimPlot(mydata,reduction = "umap",group.by = "orig.ident")

dir.create("Clustering")
ggsave("Clustering/tsne.pdf",plot = tsne_plot)
ggsave("Clustering/umap.pdf",plot = umap_plot)

tsne_plot2<-DimPlot(mydata,reduction = "tsne",label = T,label.size=5)
tsne_plot3<-DimPlot(mydata,reduction = "tsne",split.by = "orig.ident",group.by = "orig.ident",ncol = 5)

ggsave("Clustering/tsne2.pdf",plot = tsne_plot2)
ggsave("Clustering/tsne3.pdf",plot = tsne_plot3)

umap_plot2<-DimPlot(mydata,reduction = "umap",label = T,label.size=5)
umap_plot3<-DimPlot(mydata,reduction = "umap",split.by = "orig.ident",group.by = "orig.ident",ncol = 5)

ggsave("Clustering/umap2.pdf",plot = umap_plot2)
ggsave("Clustering/umap3.pdf",plot = umap_plot3)


#######UMAP
mydata <- RunUMAP(mydata, reduction = "pca", dims = 1:20)
p1 <- DimPlot(mydata, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(mydata, reduction = "umap", label = TRUE, repel = TRUE) +
  NoLegend()
plot_grid(p1, p2)

########STNE
mydata <- RunTSNE(mydata, reduction = "pca", dims = 1:20)
p3 <- DimPlot(mydata, reduction = "tsne", group.by = "orig.ident")
p4 <- DimPlot(mydata, reduction = "tsne", label = TRUE, repel = TRUE) +
  NoLegend()
plot_grid(p3, p4)

saveRDS(mydata,file = "mydata_cluster.rds")

# DEG analysis ---------------------------------------------------------

mydata<-readRDS("mydata_cluster.rds")

#识别差异基因
scRNA.pos.markers<-FindAllMarkers(mydata,only.pos=FALSE,min.pct = 0.25,logfc.threshold = 0.25)
scRNA.up.marker<-FindAllMarkers(mydata,only.pos=TRUE,min.pct = 0.25,logfc.threshold = 0.25)
#only.pos=TRUE表示只保留上调的差异基因

write.csv(scRNA.pos.markers,file = "Allmarker20.csv")
write.csv(scRNA.up.marker,file = "Upmarker20.csv")
#该markers文件用于后边的细胞类型注释

head(scRNA.pos.markers)

Top_pos_genes<-scRNA.up.marker%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
#找出各个cluster中的TOP差异基因，按logFC值排序，n=2展示各个cluster的TOP2

# 展示候选marker基因的表达量 --------------------------------------------------------
#从上一部分“输出marker基因”输出的cell_type_markers.csv文件里选出每类p_val_adj排名第一的基因，例如第II类(class为1)的ENSG00000205364。

library("RColorBrewer")
display.brewer.all(type = "qual")

top<-Top_pos_genes$gene%>%unique()

FeaturePlot(mydata, reduction = "tsne",
            features = top[1:9], 
            ncol  = 3, #画在2列里
            cols = c("#66C2A5","grey","#E41A1C"), #用三种颜色，分别展示低表达、中表达和高表达
            #do.return =T, #返回ggplot2对象，便于组图
            #no.legend = F, #显示图例
            pt.size = 1
) 
ggsave("pos_marker_C0-C8.pdf")

## 细胞类型手动注释Cell types annotation ------------------------------------------------------------------

#合并数据，并进行预处理
#GG <- merge(G1, G2)
#rm(G1)
#rm(G2)
#对合并数据重新进行细胞自动注释
cellAnnotation <- function(obj,markerList,assay='SCT',slot='data'){
  markergene <- unique(do.call(c,markerList))
  Idents(obj) <- 'seurat_clusters'
  cluster.averages <- AverageExpression(obj,assays=assay, slot = slot,features=markergene, return.seurat = TRUE)
  if(assay=='SCT'){
    scale.data <- cluster.averages@assays$SCT@data
  }else{
    scale.data <- cluster.averages@assays$RNA@data
  }
  print(scale.data)
  cell_score <- sapply(names(markergeneList),function(x){
    tmp <- scale.data[rownames(scale.data)%in%markergeneList[[x]],]
    if(is.matrix(tmp)){
      if(nrow(tmp)>=2){
        res <- apply(tmp,2,max)
        return(res)
      }else{
        return(rep(-2,ncol(tmp)))
      }
    }else{
      return(tmp)
    }
  })
  print(cell_score)
  celltypeMap <- apply(cell_score,1,function(x){
    colnames(cell_score)[which(x==max(x))]
  },simplify = T)
  obj@meta.data$cellType_auto <- plyr::mapvalues(x = obj@active.ident, from = names(celltypeMap), to = celltypeMap)
  return(obj)
}
lymphocyte <- c('CD3D','CD3E','CD79A','MS4A1','MZB1')
myeloid <- c('CD68','CD14','TPSAB1' , 'TPSB2','CD1E','CD1C','LAMP3', 'IDO1')
EOC <- c('EPCAM','KRT19','CD24')
fibo_gene <- c('DCN','FAP','COL1A2')
endo_gene <- c('PECAM1','VWF')
markergeneList <- list(lymphocyte=lymphocyte,myeloid=myeloid,Epi=EOC,fibo=fibo_gene,endo=endo_gene)
mydata_anno <- cellAnnotation(obj=mydata,assay = 'RNA',markerList=markergeneList)

d1<-DimPlot(mydata_anno, reduction = "tsne", group.by = "cellType_auto")
#可视化注释后的数据
d2<-DimPlot(mydata, reduction = "tsne", label = TRUE, repel = TRUE) +NoLegend()

d1+d2

##展示细胞簇的基因表达情况
genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A','TNFRSF18','SLC4A10','IL7R',
                   #T cell
                   'CD19', 'CD79A', 'MS4A1' ,
                   #B cell
                   'IGHG1', 'MZB1', 'SDC1',
                   #Plasma cell
                   'CD68', 'CD163', 'CD14', 'JCHAIN',
                   #macrophage
                   'TPSAB1' , 'TPSB2',
                   #mast cell
                   'RCVRN','FPR1' , 'ITGAM' ,
                   #dentritic cell
                   'C1QA',  'C1QB',
                   #alveolar macrophage
                   'S100A9', 'S100A8','CD86', 'MMP19',
                   #epithelial cell
                   'LAMP3', 'IDO1','IDO2',
                   #alveolar type II cell
                   'CD1E','CD1C',
                   #dentritic cell
                   'KRT86','GNLY',
                   #NK/cytotoxic T cell
                   'FGF7','MME','ACTA2','GFPT2','CNN1','CNN2',
                   #muscularis cell
                   'DCN', 'LUM',  'GSN' ,
                   #stromal cell
                   'FAP','FN1','THY1','COL1A1','COL3A1', 
                   #stromal cell
                   'PECAM1', 'VWF',
                   #vascular cell
                   'EPCAM' , 'KRT19', 'PROM1', 'CD24','MKI67',
                   #cycling cell
                   'ALDH1A1', 'ALDH1A3','ITGA4','ITGA6'
                   #cancer stem cell
                   )
DotPlot(mydata_anno,group.by = 'seurat_clusters', features = unique(genes_to_check),cluster.idents = T) + coord_flip()

#根据基因表达图，手动注释细胞类型
B_P_cell=c(6,17)
Mast_cell=c(13)
T_cell=c(0,1,2,3,19)
Epithelial=c(14,10,18,9)
Endothelial =c(16)
Myeloid=c(8,7,4,5)
Fibro_cell=c(15)
NK_cell=c(11,12)
table(c(B_P_cell,Mast_cell,T_cell,Epithelial,Endothelial,Myeloid,Fibro_cell,NK_cell))
print(setdiff(0:20,c(B_P_cell,Mast_cell,T_cell,Epithelial,Endothelial,Myeloid,Fibro_cell,NK_cell)))
current.cluster.ids <- c(B_P_cell,Mast_cell,T_cell,Epithelial,Endothelial,Myeloid,Fibro_cell,NK_cell)
new.cluster.ids <- c(rep("B cell",length(B_P_cell)),
                     rep("Mast cell",length(Mast_cell)),
                     rep("T cell",length(T_cell)),
                     rep("Epithelial cell",length(Epithelial)),
                     rep("Endothelial cell",length(Endothelial)),
                     rep("Myeloid cell",length(Myeloid)),
                     rep("Fibrocyte",length(Fibro_cell)),
                     rep("NK cell",length(NK_cell)))

mydata@meta.data$pre_cellType <- plyr::mapvalues(x = mydata$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)
mydata<- subset(mydata, pre_cellType %in% c("T cell","B cell","Myeloid cell","Epithelial cell","Fibrocyte",
                                            "NK cell","Endothelial cell","Mast cell"))

# 细胞浸润差异分析 ----------------------------------------------------------------

mydata<-readRDS("D:/workplace/Immunotherapy/mydata_annotation.rds")

#展示不同组织细胞类型的UMAP图，可以发现两个组织之间的细胞浸润是存在差异的，于是我们针对这个现象进行统计分析
#首先进行CellMarker的展示，证明细胞类型之间具有明显的区分度

Idents(mydata) <- 'pre_cellType'
cellType_markerGene <- FindAllMarkers(mydata,logfc.threshold = 0.5,only.pos = T,test.use = 'roc',min.pct = 0.2)
cellType_markerGene <- subset(cellType_markerGene,!grepl(pattern = 'RP[LS]',gene))
cellType_markerGene <- subset(cellType_markerGene,!grepl(pattern = 'MT-',gene))
top5 <- cellType_markerGene %>% group_by(cluster) %>% top_n(n = 10, wt = myAUC)
print(top5)
B_P_gene <- c('MS4A1','CD79A','IGHG1','JCHAIN')
T_gene <- c('CD3D','CD3E')
Epithelial_gene <- c('EPCAM','CD24')
Endothelial_gene <- c('PECAM1','VWF')
Myeloid_gene <- c('CD68','CD14')
Mast_gene <- c('TPSAB1','TPSB2')
Fibro_gene <- c('DCN','LUM')
NK_gene<-c('KLRD1','GNLY')
cellMarker <- c(B_P_gene,T_gene,Epithelial_gene,Endothelial_gene,Myeloid_gene,
                Mast_gene,Fibro_gene,NK_gene)
mydata$pre_cellType <- factor(mydata$pre_cellType,levels=c("B cell","T cell","Epithelial cell","Endothelial cell","Myeloid cell","Mast cell","Fibrocyte","NK cell"))
DotPlot(mydata,group.by = 'pre_cellType',cols = c('blue','red'), features = unique(cellMarker),cluster.idents = F) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

d1<-DimPlot(mydata, group.by = 'pre_cellType',reduction = "tsne", label = TRUE, repel = TRUE,
            cols =c("#00a8e1","#99cc00","#e30039","#fcd300","#800080","#00994e","#ff6600","#db00c2"))
d2<-DimPlot(mydata, group.by = 'seurat_clusters',reduction = "tsne", label = TRUE, repel = TRUE) 
d1+d2

saveRDS(mydata,file = "mydata_annotation.rds")

d3<-DimPlot(mydata, group.by = 'orig.ident',reduction = "tsne", label = TRUE, repel = TRUE) 
d4<-DimPlot(mydata, group.by = 'treatment',reduction = "tsne", label = TRUE, repel = TRUE) 
d3+d4


###比较分析
library("ggsignif")
library("ggpubr")

cell.Paired <- aggregate(mydata$treatment, list(mydata$pre_cellType, mydata$treatment, mydata$orig.ident, mydata$seurat_clusters), length)

#R
cell.pro.Paired.R <- subset(cell.Paired, Group.2 == "Responded")
cell.pro.Paired.R$proportion <- cell.pro.Paired.R$x/sum(cell.pro.Paired.R$x)
  
#NR
cell.pro.Paired.N <- subset(cell.Paired, Group.2 == "Non-responded")
cell.pro.Paired.N$proportion <- cell.pro.Paired.N$x/sum(cell.pro.Paired.N$x)

cell.pro.Paired <- rbind(cell.pro.Paired.R, cell.pro.Paired.N)

head(cell.Paired)
head(cell.pro.Paired)
compare_means(proportion ~ Group.2, data = cell.pro.Paired, group.by = "Group.1")
options(repr.plot.height = 6, repr.plot.width = 7)
ggboxplot(cell.pro.Paired, x = "Group.1", y = "proportion",
          color = "Group.2", palette = "jco", 
          add = "jitter")+# palette可以按照期刊选择相应的配色，如"npg"等
  stat_compare_means(aes(group = Group.2), label = "p.signif")+
  theme(panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_blank(),
        text=element_text(size = 12, family = "serif"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))
### W1 VS W2
cell.pro.Paired1<-subset(cell.pro.Paired,cell.pro.Paired$Group.3%in%c("W1","W2"))
compare_means(proportion ~ Group.2, data = cell.pro.Paired1, group.by = "Group.1")
options(repr.plot.height = 6, repr.plot.width = 7)
ggboxplot(cell.pro.Paired1, x = "Group.1", y = "proportion",
          color = "Group.2", palette = "jco", 
          add = "jitter")+# palette可以按照期刊选择相应的配色，如"npg"等
  stat_compare_means(aes(group = Group.2), label = "p.signif")+
  theme(panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_blank(),
        text=element_text(size = 12, family = "serif"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

### W3 VS W2
cell.pro.Paired3<-subset(cell.pro.Paired,cell.pro.Paired$Group.3%in%c("W3","W2"))
compare_means(proportion ~ Group.2, data = cell.pro.Paired3, group.by = "Group.1")
options(repr.plot.height = 6, repr.plot.width = 7)
ggboxplot(cell.pro.Paired3, x = "Group.1", y = "proportion",
          color = "Group.2", palette = "jco", 
          add = "jitter")+# palette可以按照期刊选择相应的配色，如"npg"等
  stat_compare_means(aes(group = Group.2), label = "p.signif")+
  theme(panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_blank(),
        text=element_text(size = 12, family = "serif"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

####
cell.pro.Paired1<-subset(cell.pro.Paired,cell.pro.Paired$Group.3%in%c("W1","W2"))
compare_means(proportion ~ Group.2, data = cell.pro.Paired1, group.by = "Group.4")
options(repr.plot.height = 6, repr.plot.width = 7)
ggboxplot(cell.pro.Paired1, x = "Group.4", y = "proportion",
          color = "Group.2", palette = "jco", 
          add = "jitter")+# palette可以按照期刊选择相应的配色，如"npg"等
  stat_compare_means(aes(group = Group.2), label = "p.signif")+
  theme(panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_blank(),
        text=element_text(size = 12, family = "serif"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

# 细胞数量饼状图 -----------------------------------------------------------------

SingleCell <- mydata
SingleCell@meta.data$pre_cellType <- factor(SingleCell@meta.data$pre_cellType,
                                                 levels = c("T cell","B cell","Myeloid cell","Epithelial cell","Fibrocyte",
                                                            "Endothelial cell","Mast cell","NK cell"))
cellNum.pie <- table(SingleCell@meta.data$pre_cellType) %>% as.data.frame()
colnames(cellNum.pie) <- c("Group","Value")
cellNum.pie$Group <- factor(cellNum.pie$Group, levels = c("T cell","B cell","Myeloid cell","Epithelial cell","Fibrocyte",
                                                          "Endothelial cell","Mast cell","NK cell"))
cellNum.pie
library(RColorBrewer)
brewer.pal.info
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

cellType_col<-c("#00a8e1","#99cc00","#e30039","#fcd300","#800080","#00994e","#ff6600","#db00c2")
library(ggforce)
options(repr.plot.height = 6, repr.plot.width = 7.5)
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+#去除没用的ggplot背景，坐标轴
  xlab("")+ylab('')+#添加颜色
  scale_fill_manual(values = cellType_col)+
  geom_arc_bar(data=cellNum.pie,
               stat = "pie",
               aes(x0=0,y0=0,r0=0,r=2,
                   amount=Value,fill=Group)
  )+#饼图
  annotate("text",x=1.4,y=1.6,label="6457",angle=-48)+
  annotate("text",x=-1.2,y=-1.8,label="602",angle=-40)+
  annotate("text",x=-1.9,y=-1,label="2394",angle=-65)+
  annotate("text",x=-1.5,y=1.6,label="1034",angle=20)+
  annotate("text",x=-1.1,y=1.8,label="158",angle=15)+
  annotate("text",x=-1,y=1.9,label="145",angle=0)+
  annotate("text",x=-0.7,y=1.9,label="230",angle=0)+
  annotate("text",x=-0.4,y=2.1,label="592",angle=0)

# HallMarker功能差异分析 ----------------------------------------------------

#获取人类的hall marker基因集
library(msigdbr)
h.human <- msigdbr(species = "Homo sapiens", category = "H")
h.names <- unique(h.human$gs_name)
h.sets <- vector("list", length=length(h.names))
names(h.sets) <- h.names
for(i in names(h.sets)){
  h.sets[[i]] <- subset(h.human, gs_name == i, "gene_symbol") %>% unlist(.)
}
#特异细胞基因集在通路hall marker活性分析
map_names = function(seur=NULL, names=NULL){
  
  # Map "names" to genes or features in a Seurat object
  # ---------------------------------------------------
  # seur = seurat object
  # names = list of names (genes or features)
  # returns list(genes=c(GENES), feats=(FEATS))
  
  # Initialize variables
  names = as.list(names)
  genes = c()
  feats = c()
  
  # Get data and metadata
  data = seur@assays[['RNA']]@data
  meta = seur@meta.data
  
  # Map names
  if(!is.null(names)){
    genes = sapply(names, function(a){intersect(a, rownames(data))}, simplify=F)
    feats = sapply(names, function(a){intersect(a, colnames(meta))}, simplify=F)
  }
  
  # Filter genes and feats
  genes = genes[lengths(genes) > 0]
  feats = feats[lengths(feats) > 0]
  
  # Fix gene names
  if(length(genes) > 0){
    if(is.null(names(genes))){names(genes) = sapply(genes, paste, collapse='.')}
  }
  
  # Fix feat names
  if(length(feats) > 0){
    if(is.null(names(feats))){names(feats) = sapply(feats, paste, collapse='.')}
  }
  
  return(list(genes=genes, feats=feats))
}


score_cells = function(seur=NULL, names=NULL, combine_genes='mean', groups=NULL, group_stat='mean', cells.use=NULL){
  
  # Score genes and features across cells and optionally aggregate
  # The steps are:
  # - calculate mean expression across genes (combine_genes = 'sum', 'mean', 'scale', 'scale2')
  # - calculate mean expression within each cell type (groups, group_stat = 'mean', 'alpha', 'mu')
  
  # Fix input arguments and get data for name mapping
  data = seur@assays[['RNA']]@data
  meta = seur@meta.data
  if(!is.null(groups)){groups = setNames(groups, colnames(seur@assays[['RNA']]@data))}
  scores = NULL
  
  # Map genes and feats
  res = map_names(seur=seur, names=names)
  genes = res$genes
  feats = res$feats
  genes.use = unique(do.call(c, genes))
  
  # Subset cells
  if(!is.null(cells.use)){
    data = data[,cells.use,drop=F]
    meta = meta[cells.use,,drop=F]
    if(!is.null(groups)){groups = groups[cells.use]}
  }
  
  group_genes = function(x, method){
    
    # combine expression data across genes within a signature
    # x = [genes x cells] matrix
    # method = 'sum', 'mean', 'scale'
    # returns [genes x cells] or [1 x cells] matrix
    
    if(nrow(x) == 1){return(x[1,,drop=F])}
    if(method == 'sum'){
      t(colSums(x))
    } else if(method == 'mean'){
      t(colMeans(x, na.rm=T))
    } else if(method == 'scale'){
      x = t(scale(t(x)))
      t(colMeans(x, na.rm=T))
    } else if(method == 'scale2'){
      x = t(scale(t(x), center=F))
      t(colMeans(x, na.rm=T))
    } else if(method == 'none'){
      x
    } else {
      stop('Error: invalid combine_genes method')
    }
  }
  
  group_cells = function(x, groups, method){
    
    # combine expression data across cells
    # x = [genes x cells] matrix
    # group_stat = 'alpha', 'mu', or 'mean'
    # returns [genes x groups] matrix
    
    if(is.null(groups)){return(x)}
    if(method %in% c('n', 'sum')){
      if(method == 'n'){x = x > 0}
      x = t(data.frame(aggregate(t(x), list(groups), sum, na.rm=T), row.names=1))
    } else {
      if(method == 'alpha'){x = x > 0}
      if(method == 'mu'){x[x == 0] = NA}
      x = t(data.frame(aggregate(t(x), list(groups), mean, na.rm=T), row.names=1))
    }
    x[is.na(x)] = 0
    x
  }
  
  # Calculate scores
  names.use = unique(c(names(genes), names(feats)))
  
  # Speed improvements (fast indexing for flat structures)
  name_map = sapply(names.use, function(a) c(genes[[a]], feats[[a]]), simplify=F)
  do.flat = all(lengths(name_map) == 1)
  if(do.flat == TRUE){
    genes[['flat']] = do.call(c, genes)
    feats[['flat']] = do.call(c, feats)
    names.iter = 'flat'
    combine_genes = 'none'
  } else {
    names.iter = names.use
  }
  
  backup = scores
  scores = lapply(names.iter, function(name){
    
    # Combine data and metadata
    if(name %in% names(genes)){
      si = data[genes[[name]],,drop=F]
    } else {
      si = c()
    }
    if(name %in% names(feats)){
      if(is.null(si)){
        si = t(meta[,feats[[name]],drop=F])
      } else {
        si = rBind(si, t(meta[,feats[[name]],drop=F]))
      }
    }
    si = as.matrix(as.data.frame(si))
    si = group_genes(si, method=combine_genes)
    si = group_cells(si, groups=groups, method=group_stat)
    si = data.frame(t(si))
  })
  
  # Collapse scores
  if(do.flat == TRUE){
    scores = scores[[1]][,make.names(name_map[names.use]),drop=F]
  } else {
    do.collapse = all(lapply(scores, ncol) == 1)
    if(do.collapse == TRUE){
      scores = as.data.frame(do.call(cbind, scores))
    }
  }
  
  # Fix names
  names.use = make.names(names.use)
  names(scores) = names.use
  
  # Combine data
  if(!is.null(backup)){
    if(is.data.frame(scores)){
      scores = cbind.data.frame(scores, backup)
    } else {
      scores = c(scores, backup)
    }
  }
  
  if(nrow(scores) == 0){scores = NULL}
  return(scores)
}
AUC_scores <- score_cells(subset(BRCA_SingleCell), names = h.sets)
library(corrplot)
library(pheatmap)
row_AUC <- rownames(AUC_scores)
AUC_scores <- cbind(AUC_scores, cells = row_AUC)
cells_meta <- data.frame(cells = rownames(BRCA_SingleCell@meta.data), pre_cellType = BRCA_SingleCell@meta.data$pre_cellType,
                         TumorType = BRCA_SingleCell@meta.data$TumorType, sampleID = BRCA_SingleCell@meta.data$Patient.id,
                         Dataset = BRCA_SingleCell@meta.data$Dataset)
AUC_scores <- merge(cells_meta, AUC_scores, by = "cells")
AUC_scores <- AUC_scores[,-1]
AUC_BCell_L <- subset(AUC_scores, pre_cellType == "B cell" & TumorType == "Lymph")
AUC_BCell_P <- subset(AUC_scores, pre_cellType == "B cell" & TumorType == "Primary")
AUC_TCell_L <- subset(AUC_scores, pre_cellType == "T cell" & TumorType == "Lymph")
AUC_TCell_P <- subset(AUC_scores, pre_cellType == "T cell" & TumorType == "Primary")
AUC_EOC_L <- subset(AUC_scores, pre_cellType == "Epithelial cell" & TumorType == "Lymph")
AUC_EOC_P <- subset(AUC_scores, pre_cellType == "Epithelial cell" & TumorType == "Primary")
AUC_Fibrocyte_L <- subset(AUC_scores, pre_cellType == "Fibrocyte" & TumorType == "Lymph")
AUC_Fibrocyte_P <- subset(AUC_scores, pre_cellType == "Fibrocyte" & TumorType == "Primary")
AUC_Myeloid_L <- subset(AUC_scores, pre_cellType == "Myeloid cell" & TumorType == "Lymph")
AUC_Myeloid_P <- subset(AUC_scores, pre_cellType == "Myeloid cell" & TumorType == "Primary")
AUC_Endothelial_L <- subset(AUC_scores, pre_cellType == "Endothelial cell" & TumorType == "Lymph")
AUC_Endothelial_P <- subset(AUC_scores, pre_cellType == "Endothelial cell" & TumorType == "Primary")
diff_analysis <- function(x, y){
  p <- c()
  fc <- c()
  for(i in 5:dim(x)[2]){
    pv <- wilcox.test(x[,i],y[,i], exact = FALSE)[3]
    fcv <- ifelse(is.infinite(log(mean(x[,i])/mean(y[,i]), 2)), 0, log(mean(x[,i])/mean(y[,i]), 2))
    fc <- c(fc, fcv)
    p <- c(p, pv)
  }
  fdr <- p.adjust(p, method = "BH") %>% log(., 10) %>% abs(.)
  result <- data.frame(FDR = fdr, LOG2FC = fc)
  return(result)
}
AUC_BCell_diff <- diff_analysis(AUC_BCell_L,AUC_BCell_P)
AUC_TCell_diff <- diff_analysis(AUC_TCell_L,AUC_TCell_P)
AUC_EOC_diff <- diff_analysis(AUC_EOC_L,AUC_EOC_P)
AUC_Fibrocyte_diff <- diff_analysis(AUC_Fibrocyte_L,AUC_Fibrocyte_P)
AUC_Myeloid_diff <- diff_analysis(AUC_Myeloid_L,AUC_Myeloid_P)
AUC_Endothelial_diff <- diff_analysis(AUC_Endothelial_L,AUC_Endothelial_P)
AUC_BCell_diff$pathway <- h.names
AUC_TCell_diff$pathway <- h.names
AUC_EOC_diff$pathway <- h.names
AUC_Fibrocyte_diff$pathway <- h.names
AUC_Myeloid_diff$pathway <- h.names
AUC_Endothelial_diff$pathway <- h.names
AUC_Cells_diff <- cbind(AUC_BCell_diff, pre_cellType = "B cell")
AUC_Cells_diff <- cbind(AUC_TCell_diff, pre_cellType = "T cell") %>% rbind(., AUC_Cells_diff)
AUC_Cells_diff <- cbind(AUC_EOC_diff, pre_cellType = "Epithelial cell") %>% rbind(., AUC_Cells_diff)
AUC_Cells_diff <- cbind(AUC_Fibrocyte_diff, pre_cellType = "Fibrocyte") %>% rbind(., AUC_Cells_diff)
AUC_Cells_diff <- cbind(AUC_Myeloid_diff, pre_cellType = "Myeloid cell") %>% rbind(., AUC_Cells_diff)
AUC_Cells_diff <- cbind(AUC_Endothelial_diff, pre_cellType = "Endothelial cell") %>% rbind(., AUC_Cells_diff)
pathway.name <- AUC_Cells_diff$pathway
pathway.name <- strsplit(pathway.name, split = "HALLMARK_")
pathway.name <- unlist(pathway.name)
pathway.name <- pathway.name[seq(2,length(pathway.name), by = 2)]
AUC_Cells_diff$pathway <- pathway.name
AUC_Cells_diff <- AUC_Cells_diff[AUC_Cells_diff$pathway %in% c("ALLOGRAFT_REJECTION","COAGULATION","COMPLEMENT",
                                                               "IL6_JAK_STAT3_SIGNALING","INFLAMMATORY_RESPONSE",
                                                               "INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE", #Immune
                                                               "BILE_ACID_METABOLISM","CHOLESTEROL_HOMEOSTASIS","FATTY_ACID_METABOLISM",
                                                               "GLYCOLYSIS","HEME_METABOLISM","OXIDATIVE_PHOSPHORYLATION",
                                                               "XENOBIOTIC_METABOLISM", #Metabolism
                                                               "E2F_TARGETS","G2M_CHECKPOINT","MITOTIC_SPINDLE","MYC_TARGETS_V1",
                                                               "MYC_TARGETS_V2", #Proliferation
                                                               "ANGIOGENESIS","EPITHELIAL_MESENCHYMAL_TRANSITION", #Metastasis
                                                               "ANDROGEN_RESPONSE","APOPTOSIS","ESTROGEN_RESPONSE_EARLY",
                                                               "ESTROGEN_RESPONSE_LATE","HEDGEHOG_SIGNALING","HYPOXIA",
                                                               "IL2_STAT5_SIGNALING","KRAS_SIGNALING_DN","KRAS_SIGNALING_UP",
                                                               "MTORC1_SIGNALING","NOTCH_SIGNALING","PI3K_AKT_MTOR_SIGNALING",
                                                               "PROTEIN_SECRETION","REACTIVE_OXYGEN_SPECIES_PATHWAY","TGF_BETA_SIGNALING",
                                                               "TNFA_SIGNALING_VIA_NFKB","UNFOLDED_PROTEIN_RESPONSE",
                                                               "WNT_BETA_CATENIN_SIGNALING"),] #Signaling
AUC_Cells_diff$pathway <- factor(AUC_Cells_diff$pathway,
                                 levels = c("ALLOGRAFT_REJECTION","COAGULATION","COMPLEMENT",
                                            "IL6_JAK_STAT3_SIGNALING","INFLAMMATORY_RESPONSE",
                                            "INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE", #Immune
                                            "BILE_ACID_METABOLISM","CHOLESTEROL_HOMEOSTASIS","FATTY_ACID_METABOLISM",
                                            "GLYCOLYSIS","HEME_METABOLISM","OXIDATIVE_PHOSPHORYLATION",
                                            "XENOBIOTIC_METABOLISM", #Metabolism
                                            "E2F_TARGETS","G2M_CHECKPOINT","MITOTIC_SPINDLE","MYC_TARGETS_V1",
                                            "MYC_TARGETS_V2", #Proliferation
                                            "ANGIOGENESIS","EPITHELIAL_MESENCHYMAL_TRANSITION", #Metastasis
                                            "ANDROGEN_RESPONSE","APOPTOSIS","ESTROGEN_RESPONSE_EARLY",
                                            "ESTROGEN_RESPONSE_LATE","HEDGEHOG_SIGNALING","HYPOXIA",
                                            "IL2_STAT5_SIGNALING","KRAS_SIGNALING_DN","KRAS_SIGNALING_UP",
                                            "MTORC1_SIGNALING","NOTCH_SIGNALING","PI3K_AKT_MTOR_SIGNALING",
                                            "PROTEIN_SECRETION","REACTIVE_OXYGEN_SPECIES_PATHWAY","TGF_BETA_SIGNALING",
                                            "TNFA_SIGNALING_VIA_NFKB","UNFOLDED_PROTEIN_RESPONSE",
                                            "WNT_BETA_CATENIN_SIGNALING"),#Signaling
                                 labels = c("Allograft rejection","Coagulation","Complement",
                                            "IL6 JAK STAT3 signaling","Inflammatory response",
                                            "Interferon alpha response","Interferon gamma response",
                                            "Bile acid metabolism","Cholesterol homeostasis","Fatty acid matebolism",
                                            "Glycolysis","Heme Metabolism","Oxidative phosphorylation",
                                            "Xenobiotic metabolism",
                                            "E2f tagets","G2M checkpoint","Mitotic spindle","MYC targets v1",
                                            "MYC targets v2",
                                            "Angiogenesis","Epithelial mesenchymal transition",
                                            "Androgen response","Apoptosis","Estrogen response early",
                                            "Estrogen response late","Hedgehog signaling","Hypoxia",
                                            "IL2 STAT5 signaling","KRAS signaling dn","KRAS signaling up",
                                            "mTORC1 signaling","Notch signaling","PI3K AKT MTOR signaling",
                                            "Protein secretion","Reactive oxygen species pathway","TGF beta signaling",
                                            "TNFA signaling via NFKB","Unfolded protein response",
                                            "WNT beta catenin signaling"))
options(repr.plot.height = 5.5, repr.plot.width = 16)
ggplot(AUC_Cells_diff, aes(pathway, forcats::fct_rev(pre_cellType), fill = LOG2FC, size = FDR)) +
  geom_point(shape = 21, stroke = 0.1) +
  geom_hline(yintercept = seq(1.5, 5.5, 1), size = 0.5, col = "black") +
  geom_vline(xintercept = seq(1.5, 38.5, 1), size = 0.5, col = "black") +
  scale_x_discrete(position = "top") +
  scale_radius(range = c(2, 6)) +
  scale_fill_gradient2(low = "blue", mid = "white",high = "red",
                       breaks = c(-round(max(abs(AUC_Cells_diff$LOG2FC)),1), 0, round(max(abs(AUC_Cells_diff$LOG2FC)),1)),
                       limits = c(-round(max(abs(AUC_Cells_diff$LOG2FC)),1), round(max(abs(AUC_Cells_diff$LOG2FC)),1))) +
  #scale_fill_gradient(low = "blue" , high = "red") +
  #scale_fill_distiller(palette = "Spectral") +
  #theme_minimal() +
  theme(legend.position = "bottom",
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid.major = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 10),
        axis.text.y = element_text(size = 12)) +
  guides(size = guide_legend(override.aes = list(fill = NA, color = "black", stroke = .25), 
                             label.position = "bottom",
                             title.position = "right", 
                             order = 1),
         fill = guide_colorbar(ticks.colour = NA, title.position = "top", order = 2)) +
  labs(size = "FDR", fill = "logFC", x = NULL, y = NULL)
