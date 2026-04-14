# myeloid cell-4,5,7,8 ----------------------------------------------------
mydata<-readRDS("mydata_annotation.rds")

myeloid<-mydata[,mydata@meta.data$pre_cellType=="Myeloid cell"]
DimPlot(myeloid, reduction = "tsne", group.by = "seurat_clusters",
        cols =c("#89c8e8","#caa2f4","#1f78b4","#6a3d9a"))

####髓系细胞细分类：DC MONOCYTE TAM MDSC
cDC=c("FCER1A","CST3")
pDC=c('IL3RA','GZMB','SERPINF1','ITM2C')
Mono=c('CD14','LYZ','S100A8','S100A9','LST1')
Mac=c('CSF1R','CD68')

ref_markers<-data.frame(CellType=c("cDC","pDC","Mono","Mac"),
                        Markers=c("FCER1A,CST3",
                                  "IL3RA,GZMB,SERPINF1,ITM2C",
                                  "CD14,LYZ,S100A8,S100A9,LST1",
                                  "CSF1R,CD68")
                        )

#####

#cluster-specific up-markers
markers<-fread("E:/Rwork/lung-GEO-8/immunothreapy/Upmarker20.csv")

pbmc.markers<-subset(markers,markers$cluster==c("4","5","7","8"))

######超几何检验#######

p<-data.frame()

for(i in 1:length(unique(pbmc.markers$cluster))){
  
  cluster_gene<-pbmc.markers$gene[pbmc.markers$cluster==unique(pbmc.markers$cluster)[i]]
  #cluster_i中特异上调的基因向量
  
  for(j in 1:dim(ref_markers)[1]){
    
    mark<-strsplit(ref_markers$Markers[j],",")[[1]]
    cell_type<-rep(ref_markers$CellType[j],length(mark))
    ref_df<-data.frame(CellType=cell_type,MarkerGene=mark)
    #cell_type_j中的marker基因
    
    N<-union(strsplit(ref_markers$Markers,",")%>%unlist(),pbmc.markers$gene)%>%length()
    #背景基因数
    q<-intersect(cluster_gene,ref_df$MarkerGene)%>%length()
    #celltype1和cluster1的交集基因数
    m<-length(ref_df$MarkerGene)
    #celltype中的基因数
    n<-N-m
    k<-length(cluster_gene)
    #cluster中的基因数
    
    p_value<-1-phyper(q,m,n,k)
    
    df_p<-data.frame(cluster=unique(pbmc.markers$cluster)[i],
                     celltype=ref_markers$CellType[j],
                     phyper_pvalue=p_value)
    p<-rbind(p,df_p)
    
  }
}

annotation_p<-subset(p,p$phyper_pvalue<0.5)

test<-p %>% group_by(cluster) %>% dplyr::mutate(min_rank(phyper_pvalue))
annotation_top<-subset(test,test$`min_rank(phyper_pvalue)`==1)
annotation_top$cluster<-paste(rep("C",length(annotation_top$cluster)),annotation_top$cluster,sep = "")

annotation_top$celltype<-paste(annotation_top$celltype,annotation_top$cluster,sep = "_")

#####重新定义细胞类型######

myeloid@meta.data$new.cluster.ids<-ifelse(myeloid@meta.data$seurat_clusters==4,"pDC_C4",
         ifelse(myeloid@meta.data$seurat_clusters==5,"Mac_C5",
                ifelse(myeloid@meta.data$seurat_clusters==7,"Mono_C7","Mac_C8")))

Idents(myeloid)
#原来细胞聚类的名称

#####可视化注释后的细胞

DimPlot(myeloid, reduction = "tsne", group.by = "new.cluster.ids",
        cols =c("#caa2f4","#6a3d9a","#1f78b4","#89c8e8"))

saveRDS(myeloid, file = "myeloid_annotation.rds")


# 四个髓系细胞富集分析 --------------------------------------------------------------
library(stringr)
library(dplyr)
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(viridisLite)
library(viridis)
rm(list = ls())

setwd("D:/workplace/Immunotherapy/")
myeloid<-readRDS("myeloid_annotation.rds")

gene<-fread("Upmarker20.csv")%>%as.data.frame()

marker<-subset(gene, cluster %in% 8)
up_marker<-marker$gene

geneID <- bitr(up_marker, fromType = "SYMBOL",
               toType = "ENTREZID", 
               OrgDb = org.Hs.eg.db) 

all.GO <- enrichGO(gene = geneID$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   readable = T)

allgg <- as.data.frame(all.GO)

# 按照pvalue取每个ONTOLOGY的前5
top5 <- allgg %>%  group_by(ONTOLOGY) %>%  arrange(pvalue) %>%  slice_head(n = 5)

df_top5 <- rbind(subset(top5, ONTOLOGY=="BP"), subset(top5, ONTOLOGY=="CC"), subset(top5, ONTOLOGY=="MF"))

df_top5$ONTOLOGY <- factor(df_top5$ONTOLOGY, levels=c('BP', 'CC', 'MF'))
df_top5$Description <- factor(df_top5$Description, levels = rev(df_top5$Description))
options(repr.plot.width=7, repr.plot.height=6.5)

mycol3 <- c('#6BA5CE', '#F5AA5F','#C7533B')
cmap <- c("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo")
p <- ggplot(data = df_top5, aes(x = Count, y = Description, fill=ONTOLOGY)) + 
  geom_bar(width = 0.5,stat = 'identity') + 
  theme_classic() + 
  scale_x_continuous(expand = c(0,0.5)) + 
  scale_fill_manual(values = alpha(mycol3, 0.66))

p <- p + theme(axis.text.y = element_blank()) +
  geom_text(data = df_top5, 
            aes(x = 0.1, y = Description, label = Description),
            size = 4.8, 
            hjust = 0) 
p <- p + geom_text(data = df_top5,
                   aes(x = 0.1, y = Description, label = geneID , color=-log10(pvalue)),
                   size = 4,
                   fontface = 'italic',
                   hjust = 0,
                   vjust = 2.7) +
  scale_colour_viridis(option=cmap[7], direction=-1) 

p <- p + labs(title = 'Enriched top 5 BP, CC and MF') + 
  theme( plot.title = element_text(size = 14, face = 'bold'),
         axis.title = element_text(size = 13),
         axis.text = element_text(size = 11), 
         axis.ticks.y = element_blank())
p


# 4类髓系细胞比例饼图 --------------------------------------------------------------------

library(devtools)
library(dplyr)
library(data.table)
library(Seurat)
library(magrittr)
library(monocle)
options(stringsAsFactors = FALSE) #禁止chr转成factor

setwd("D:/workplace/Immunotherapy/")
myeloid<-readRDS("myeloid_annotation.rds")

myeloid$new.cluster.ids<-as.character(myeloid$new.cluster.ids)
t<-table(myeloid$new.cluster.ids)

#####pie1

pie(as.numeric(t),
    names(t),
    radius = 1.0,
    clockwise = T,
    main="Cell proportion",
    col =c("#caa2f4","#6a3d9a","#1f78b4","#89c8e8"))

numbers=as.numeric(t)
types= names(t)

#####pie2

library(ggplot2)
library(ggforce)

A <- data.frame(numbers, types)

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
  scale_fill_manual(values = c('#E5D2DD', '#476D87', '#F1BB72', '#F3B1A0'))+
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=numbers,fill=types)
  )


# monocle -----------------------------------------------------------------

# 构建monocle对象
#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(myeloid@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = myeloid@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())


# 数据过滤 --------------------------------------------------------------------

HSMM<-monocle_cds
## 归一化 
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
#Filtering low-quality cells
HSMM <- detectGenes(HSMM, min_expr = 3 )
print(head(fData(HSMM)))

expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))
print(head(pData(HSMM)))


# 轨迹构建 -------------------------------------------------------

diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = "~new.cluster.ids")
#大于10个细胞中表达的基：expressed_genes
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01)) 
#显著差异的基因并排序

HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)

####降维 Trajectory step 2: reduce data dimensionality  
HSMM <- reduceDimension(HSMM, max_components = 2,
                        method = 'DDRTree')
##排序 #Trajectory step 3: order cells along the trajectory  
#devtools::load_all("C:/Users/zyk15/Documents/R/win-library/4.2/monocle")

HSMM <- orderCells(HSMM)
plot_cell_trajectory(HSMM, color_by = "new.cluster.ids")

###轨迹
plot_cell_trajectory(HSMM, color_by = "State")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
plot_cell_trajectory(HSMM, color_by = "new.cluster.ids") +
  facet_wrap(~new.cluster.ids, nrow = 1)

saveRDS(HSMM, file = "pseudotime-myeloid.rds")
#HSMM<-readRDS("pseudotime-myeloid.rds")

# #特征基因热图 -----------------------------------------------------------------

expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 10))

diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.01))
##383个拟时序列中差异表达的特征基因

plot.new()
plot_pseudotime_heatmap(HSMM[sig_gene_names[1:100],],
                          num_clusters = 2,
                          cores = 1,
                          show_rownames = T,
                          return_heatmap=F,
                          use_gene_short_name = TRUE)

clusters <- cutree(p$tree_row, k = 2)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)

write.table(clustering,"modelu_cluster2.txt",quote = F,row.names = T)

# 特征基因的功能 -----------------------------------------------------------------

##对热图中的四个基因模块进行功能富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(ggplot2)
library("RColorBrewer")
#display.brewer.all()

col_df<-data.frame(c1=c("#DECAE0","#C9D8EB"),
                   c2=c("#8B5C9E","#52A5C7"))

for(i in 1:2){
  
  gene<-subset(clustering,clustering$Gene_Clusters==i,)%>%rownames()
  ID <- bitr(gene, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
  ego <- enrichGO(
    gene=ID$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "ALL", # 或MF或CC
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = 1
  )
  ego1 <- setReadable(ego, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  ego2 <- ego1[ego1$pvalue < 0.05, asis = T]
  go_res<-ego2@result%>%arrange(p.adjust)
  write.csv(go_res,file=paste("GO_module_",i,".csv",sep = ""))
  
  ggplot(data=go_res[1:10,],aes(x=reorder(as.factor(Description),-log10(pvalue)),y=-log10(pvalue),fill=pvalue))+
    geom_bar(stat = "identity")+
    scale_fill_gradient(low=col_df[i,2], high=col_df[i,1])+
    coord_flip()+
    labs(title=paste("GO-term C",i,sep=""),y="Enrichment Score",x="")
  ggsave(paste("GO_cluster",i,".pdf",sep = ""),width = 15, height = 10)
}


# 细胞衰老功能分析 ------------------------------------------------------------------
sce<-readRDS("myeloid_annotation.rds")
sene_geneset<-clusterProfiler::read.gmt("REACTOME_CELLULAR_SENESCENCE.v2024.1.Hs.gmt")
#msigdb_细胞衰老基因集有197个基因
genes<-sene_geneset$gene
sce<-AddModuleScore(object = sce,
                    features = list(genes),
                    name = "CellSenescence",
                    assay = "RNA", 
                    search = T)

FeaturePlot(sce, reduction = "umap",
            features = "CDKN1A", 
            cols = c("grey","#E2696F","#E41A1C"), #用三种颜色，分别展示低表达、中表达和高表达
            pt.size = 0.3)
FeaturePlot(sce, reduction = "tsne",
            features = "CCL2", 
            cols = c("grey","#E2696F","#E41A1C"), #用三种颜色，分别展示低表达、中表达和高表达
            pt.size = 0.3)
FeaturePlot(sce, reduction = "tsne",
            features = "IL6", 
            cols = c("grey","#E2696F","#E41A1C"), #用三种颜色，分别展示低表达、中表达和高表达
            pt.size = 0.3)
FeaturePlot(sce, reduction = "tsne",
            features = "CellSenescence1", 
            #cols = c("grey","#E2696F","#E41A1C"), #用三种颜色，分别展示低表达、中表达和高表达
            pt.size = 0.3)


# cytotrace2 -----------------------------------genes# cytotrace2 --------------------------------------------------------------


# 本地安装
#remotes::install_local("D:/workplace/Immunotherapy/digitalcytometry-cytotrace2-6e8041f.tar.gz",
                       #subdir = "cytotrace2_r", # 特殊的
                       #upgrade = F,dependencies = T)
library(CytoTRACE2)
library(Seurat)
library(tidyverse)
library(cowplot)

DimPlot(sce)
Idents(sce)<-"new.cluster.ids"
DimPlot(sce, reduction = "umap", group.by = "new.cluster.ids",
        cols =c("#caa2f4","#6a3d9a","#1f78b4","#89c8e8"))
#计算完的结果会自动添加到metadata中
cytotrace2_result_sce<-cytotrace2(sce,
                                    is_seurat=TRUE,
                                    slot_type="counts",
                                    species='human',#默认是小鼠，这里是人的数据所以改为human
                                    seed=123)

#把原来的注释结果提取出来进行可视化
anno<-data.frame(phenotype=sce@meta.data$new.cluster.ids)
rownames(anno)<-colnames(sce)
plots<-plotData(cytotrace2_result=cytotrace2_result_sce,
                  annotation=anno,
                  is_seurat=T)
p1<-plots$CytoTRACE2_UMAP
p2<-plots$CytoTRACE2_Potency_UMAP
p3<-plots$CytoTRACE2_Relative_UMAP
p4<-plots$Phenotype_UMAP
p5<-plots$CytoTRACE2_Boxplot_byPheno
plot_grid(p1,p2,p3,p4,p5,nrow=2)#简单拼图

saveRDS(cytotrace2_result_sce, file = "sce-myeloid.rds")

setwd("D:/workplace/Immunotherapy/")
sce<-readRDS("sce-myeloid.rds")


# C8,c5巨噬细胞 -------------------------------------------------------------

setwd("D:/workplace/Immunotherapy/")
myeloid<-readRDS("myeloid_annotation.rds")
library(data.table)
library(dplyr)
library(stringr)
library(Seurat)
library(ggplot2)

myeloid@meta.data%>%head()
mac58<-myeloid[,myeloid@meta.data$new.cluster.ids==c("Mac_C5","Mac_C8")]
#提出C5,C8
mac58@meta.data$new.cluster.ids%>%unique()

##定义M1(C8),M2(C5)

##展示细胞簇的基因表达情况
genes_to_check = c('CD80', 'CD86', 'CD16','CD32','IL6',
                   #M1
                   'CD206', 'CD163','CD68','TGFB','IL10'
                   #M2
)
DotPlot(mac58,group.by = 'new.cluster.ids', features = unique(genes_to_check),cluster.idents = T) + coord_flip()

#####可视化C5C8巨噬细胞及特征基因集表达

DimPlot(mac58, reduction = "tsne", group.by = "new.cluster.ids",
        cols =c("#caa2f4","#6a3d9a"))

FeaturePlot(mac58, reduction = "tsne",
            features = "CD274", 
            cols = c("grey","#E2696F","#E41A1C"), #用三种颜色，分别展示低表达、中表达和高表达
            pt.size = 0.8)

###计算基因集得分

## M1
M1<-read.table("M1_Polarization.txt")
genes<-M1$V1
mac58<-AddModuleScore(object = mac58,
                    features = list(genes),
                    name = "M1_Polarization",
                    assay = "RNA", 
                    search = T)

FeaturePlot(mac58, reduction = "tsne",
            features = "M1_Polarization1", 
            cols = c("grey","#E2696F","#E41A1C"), #用三种颜色，分别展示低表达、中表达和高表达
            pt.size = 0.8)

## M2
M2<-read.table("M2_Polarization.txt")
genes<-M2$V1
mac58<-AddModuleScore(object = mac58,
                      features = list(genes),
                      name = "M2_Polarization",
                      assay = "RNA", 
                      search = T)

FeaturePlot(mac58, reduction = "tsne",
            features = "M2_Polarization1", 
            cols = c("grey","#E2696F","#E41A1C"), #用三种颜色，分别展示低表达、中表达和高表达
            pt.size = 0.8)

## Anti_inflammatory

anti_inf<-read.table("Anti_inflammatory.txt")
genes<-anti_inf$V1
mac58<-AddModuleScore(object = mac58,
                      features = list(genes),
                      name = "Anti_inflammatory",
                      assay = "RNA", 
                      search = T)

FeaturePlot(mac58, reduction = "tsne",
            features = "Anti_inflammatory1", 
            cols = c("grey","#E2696F","#E41A1C"), #用三种颜色，分别展示低表达、中表达和高表达
            pt.size = 0.8)

## Pro_inflammatory

pro_inf<-read.table("Pro_inflammatory.txt")
genes<-pro_inf$V1
mac58<-AddModuleScore(object = mac58,
                      features = list(genes),
                      name = "Pro_inflammatory",
                      assay = "RNA", 
                      search = T)

FeaturePlot(mac58, reduction = "tsne",
            features = "Pro_inflammatory1", 
            cols = c("grey","#E2696F","#E41A1C"), #用三种颜色，分别展示低表达、中表达和高表达
            pt.size = 0.8)

## Immune_surveillance

immu1<-read.table("Immune_surveillance.txt")
genes<-immu1$V1
mac58<-AddModuleScore(object = mac58,
                      features = list(genes),
                      name = "Immune_surveillance",
                      assay = "RNA", 
                      search = T)

FeaturePlot(mac58, reduction = "tsne",
            features = "Immune_surveillance1", 
            cols = c("grey","#E2696F","#E41A1C"), #用三种颜色，分别展示低表达、中表达和高表达
            pt.size = 0.8)

### Immune_escape

immu2<-read.table("Immune_escape.txt")
genes<-immu2$V1
mac58<-AddModuleScore(object = mac58,
                      features = list(genes),
                      name = "Immune_escape",
                      assay = "RNA", 
                      search = T)

FeaturePlot(mac58, reduction = "tsne",
            features = "Immune_escape1", 
            cols = c("grey","#E2696F","#E41A1C"), #用三种颜色，分别展示低表达、中表达和高表达
            pt.size = 0.8)

#####C5,C8得分比较

head(mac58@meta.data)

saveRDS(mac58,"Mac_c5_c8.rds")
mac58<-readRDS("Mac_c5_c8.rds")
vio_input<-mac58@meta.data
vio_input$new.cluster.ids<-as.factor(vio_input$new.cluster.ids)


ggplot(vio_input, aes(x=new.cluster.ids, y=M1_Polarization1,fill=new.cluster.ids)) + 
  geom_violin(trim=FALSE,color="white") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线图
  scale_fill_manual(values = c("#caa2f4","#6a3d9a"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=20), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=16),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("M1_Polarization")+xlab("") #设置x轴和y轴的标题

### 组间差异p值
wilcox.test(M1_Polarization1~treatment,data = vio_input)[[3]]

t.test(vio_input$M1_Polarization1 ~ vio_input$treatment)[[3]]


# 免疫治疗分组比较 ----------------------------------------------------------------

####可视化

vio_input<-mac58@meta.data
vio_w1_w2<-subset(vio_input,vio_input$orig.ident%in%c("W1","W2"))
vio_w1_w2$treatment<-as.factor(vio_w1_w2$treatment)


ggplot(vio_w1_w2, aes(x=treatment, y=M1_Polarization1,fill=treatment)) + 
  geom_violin(trim=FALSE,color="white") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线图
  scale_fill_manual(values = c("#0c6db3", "#e4b918"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=20), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=16),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("M1_Polarization")+xlab("") #设置x轴和y轴的标题

### 组间差异p值
#wilcox.test(M1_Polarization1~treatment,data = vio_w1_w2)[[3]]

t.test(vio_w1_w2$M1_Polarization1 ~ vio_w1_w2$treatment)[[3]]


# 应答组和非应答组中比较M1/M2的比例 -----------------------------------------------------

table(vio_w1_w2$treatment,vio_w1_w2$new.cluster.ids)

#查看各组细胞数
table(mac58@meta.data$orig.ident,mac58@meta.data$new.cluster.ids)

Cellratio <- prop.table(table(mac58@meta.data$new.cluster.ids,mac58@meta.data$orig.ident)
, margin = 2) #计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)

allcolour=c("#caa2f4","#6a3d9a")

library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = allcolour)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))

#较高的M1/M2比率反映了促炎M1样TAM的主导地位，增强了T细胞对肿瘤的免疫反应并提高了免疫治疗的有效性。
#相反，较低的M1/M2比值表明免疫抑制性M2样TAM相对增加，促进肿瘤免疫逃逸，削弱免疫反应，促进肿瘤耐药性


# 转移，增殖基因与c5,c8的表达相关性分析 ---------------------------------------------------

emt_geneset<-clusterProfiler::read.gmt("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2024.1.Hs.gmt")
#msigdb_emt有200个基因
emt<-emt_geneset$gene

mac58<-readRDS("Mac_c5_c8.rds")

mydata2<-AddModuleScore(object = mac58,
                        features = list(emt),
                        name = "EMT",
                        assay = "RNA", 
                        search = T)

FeaturePlot(mydata2, reduction = "tsne",
            features = "EMT1", 
            pt.size = 0.3)

library(ggplot2)
library(ggpubr)
wilcox.test(mydata2@meta.data$EMT1~mydata2@meta.data$seurat_clusters)

vio_input_emt<-data.frame(EMT=as.numeric(mydata2@meta.data$EMT1),
                          group=as.factor(mydata2@meta.data$seurat_clusters),
                          celltype=as.factor(mydata2@meta.data$pre_cellType))


ggplot(vio_input_emt, aes(x=group, y=EMT, fill=group)) + 
  geom_bar(stat="summary", fun="mean")+
  scale_fill_manual(values = c("#0c6db3", "#e4b918"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  stat_compare_means(method = "t.test")+
  ylab("Expression of EMT")+xlab("") #设置x轴和y轴的标题



####增殖-G2M检查点

G2M_geneset<-clusterProfiler::read.gmt("HALLMARK_G2M_CHECKPOINT.v2024.1.Hs.gmt")
#msigdb_G2M有200个基因
G2M<-G2M_geneset$gene
mac58<-AddModuleScore(object = mac58,
                      features = list(G2M),
                      name = "G2M",
                      assay = "RNA", 
                      search = T)

FeaturePlot(mac58, reduction = "tsne",
            features = "G2M1", 
            pt.size = 0.3)
#19575*583
####

up_marker<-read.csv("Upmarker20.csv")

c5_marker<-subset(up_marker,up_marker$cluster=="5")
c5_marker$gene
#1568

c8_marker<-subset(up_marker,up_marker$cluster=="8")
c8_marker$gene
#588

identical(colnames(mac58@assays$RNA),rownames(mac58@meta.data))

###与EMT得分正相关的C5marker
cor_df<-data.frame()
  
for( i in 1:length(c5_marker$gene)){
  
  g<-c5_marker$gene[i]
  gene_exp<-as.numeric(mac58@assays$RNA[g,])
  cor_res<-cor.test(gene_exp,mac58@meta.data$EMT1)
  cor<-cor_res[[4]]
  p<-cor_res[[3]]
  cor_df2<-data.frame(gene=g,
                      cor=cor,
                      p=p)
  cor_df<-rbind(cor_df,cor_df2)
}

#p<0.05,r>0.4

c5_cor<-subset(cor_df,cor>0.4)
c5_cor<-subset(c5_cor,p<0.05)
#126个显著正相关基因

###与EMT得分负相关的C8marker
cor_df_c8<-data.frame()

for( i in 1:length(c8_marker$gene)){
  
  g<-c8_marker$gene[i]
  gene_exp<-as.numeric(mac58@assays$RNA[g,])
  cor_res<-cor.test(gene_exp,mac58@meta.data$EMT1)
  cor<-cor_res[[4]]
  p<-cor_res[[3]]
  cor_df2<-data.frame(gene=g,
                      cor=cor,
                      p=p)
  cor_df_c8<-rbind(cor_df_c8,cor_df2)
}

#p<0.05,r<-0.4

c8_cor<-subset(cor_df_c8,cor<=-0.4)
#7个显著负相关基因

EMT_gene<-c(c8_cor$gene,c5_cor$gene)
#133个肿瘤转移相关的巨噬细胞marker

write.table(EMT_gene,"EMT_M1_M2_marker.txt",row.names = F,
            col.names = F,quote = F)

saveRDS(mac58,"Mac_c5_c8.rds")

# emt相关性图 -----------------------------------------------------------------
setwd("D:/workplace/Immunotherapy/")
mac58<-readRDS("Mac_c5_c8.rds")

library(ggstatsplot)#加载包
library(ggside)
library(ggplot2)

######C5-EMT相关性

c5_mean<-apply(mac58@assays$RNA[c5_cor$gene,],2,mean)

identical(names(mac58$EMT1),names(c5_mean))

cor_5_emt<-data.frame(emt=as.numeric(mac58$EMT1),
                      c5=as.numeric(c5_mean))

ggscatterstats(data = cor_5_emt ,
               x = emt, 
               y = c5, 
               centrality.para = "mean",
               margins = "both",
               xfill = "#009E73",
               yfill = "#caa2f4",
               marginal.type = "histogram")

######C8-EMT相关性

c8_mean<-apply(mac58@assays$RNA[c8_cor$gene,],2,mean)

identical(names(mac58$EMT1),names(c8_mean))

cor_8_emt<-data.frame(emt=as.numeric(mac58$EMT1),
                      c8=as.numeric(c8_mean))

ggscatterstats(data = cor_8_emt ,
               x = emt, 
               y = c8, 
               centrality.para = "mean",
               margins = "both",
               xfill = "#009E73",
               yfill = "#caa2f4",
               marginal.type = "histogram")


# PD-1组间差异 ----------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(ggpubr)

wilcox.test(as.numeric(mydata2@assays$RNA["PDCD1",])~mydata2@meta.data$treatment)

vio_input_PD1<-data.frame(PDCD1=as.numeric(mydata2@assays$RNA["PDCD1",]),
                          group=as.factor(mydata2@meta.data$treatment),
                          celltype=as.factor(mydata2@meta.data$pre_cellType))


ggplot(vio_input_PD1, aes(x=group, y=PDCD1, fill=group)) + 
  geom_bar(stat="identity")+
  scale_fill_manual(values = c("#0c6db3", "#e4b918"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  stat_compare_means(aes(label=paste0("p = ",after_stat(p.format))))+
  ylab("Expression of PDCD1")+xlab("") #设置x轴和y轴的标题


