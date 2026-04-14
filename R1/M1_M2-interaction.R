
# c5,c8与其他细胞互作比较 ----------------------------------------------------------

#BiocManager::install("BiocNeighbors",force = TRUE)
#devtools::install_github("sqjin/CellChat")
library(CellChat)
library(patchwork)
library(ggalluvial)
library(igraph)
library(dplyr)

setwd("D:/workplace/Immunotherapy/cellchat/")

mydata<-readRDS("mydata_annotation.rds")

#设置硬件参数，8线程
suppressWarnings(suppressMessages(future::plan("multiprocess", workers = 8)))


# 分开两组数据

data.input <-  mydata@assays$RNA # normalized data matrix

meta <- mydata@meta.data
meta$group<-as.factor(meta$treatment)

unique(meta$group)
cell.use1 <- rownames(meta)[meta$group == "Non-responded"] 
cell.use2 <- rownames(meta)[meta$group == "Responded"] 
# 按指定的变量提取细胞

data.input1 <- data.input[,cell.use1]#取出对应细胞,也就是说，data的列名是meta的行名
data.input2 <- data.input[,cell.use2]
data.input2[1:5,1:5]

meta1 = meta[cell.use1, ]#取出对应细胞的meta信息
meta2 = meta[cell.use2, ]
unique(meta2$pre_cellType)#看meta中储存的细胞注释信息，稍后用它作为分组依据

#####创建cellchat对象

cellchat1 <- createCellChat(object = data.input1, meta = meta1, group.by = "pre_cellType")
cellchat2 <- createCellChat(object = data.input2, meta = meta2, group.by = "pre_cellType")

####分别处理这两组样本

#创建celllchat对象，group.by指定通讯间的对象，用meta中的注释作为分组依据
cellchat1 <- addMeta(cellchat1, meta = meta1)
cellchat1 <- setIdent(cellchat1, ident.use = "pre_cellType") # set "labels" as default cell identity
groupSize1 <- as.numeric(table(cellchat1@idents)) #每种细胞的细胞数量

cellchat2 <- addMeta(cellchat2, meta = meta2)
cellchat2 <- setIdent(cellchat2, ident.use = "pre_cellType") # set "labels" as default cell identity
groupSize2 <- as.numeric(table(cellchat2@idents)) #每种细胞的细胞数量

#######设置参考数据库

CellChatDB <- CellChatDB.human
#查看数据库的组成比例
showDatabaseCategory(CellChatDB)
# 查看数据库具体信息
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")#取出相应分类用作分析数据库
cellchat1@DB <- CellChatDB.use#将数据库内容载入cellchat对象中，相当于设置好接下来要参考的数据库
cellchat2@DB <- CellChatDB.use#将数据库内容载入cellchat对象中，相当于设置好接下来要参考的数据库

##########表达数据预处理

cellchat1 <- subsetData(cellchat1)#取出表达数据
cellchat1 <- identifyOverExpressedGenes(cellchat1)#寻找高表达的基因#
cellchat1 <- identifyOverExpressedInteractions(cellchat1)#寻找高表达的通路
cellchat1 <- projectData(cellchat1, PPI.human)#投影到PPI，储存上一步的结果到cellchat@LR$LRsig

cellchat2 <- subsetData(cellchat2)#取出表达数据
cellchat2 <- identifyOverExpressedGenes(cellchat2)#寻找高表达的基因#
cellchat2<- identifyOverExpressedInteractions(cellchat2)#寻找高表达的通路
cellchat2 <- projectData(cellchat2, PPI.human)#投影到PPI，储存上一步的结果到cellchat@LR$LRsig


#####推断受体配体信号网络

## 1.在配受体水平上计算细胞通讯
#计算细胞与细胞之间通信的概率
cellchat1 <- computeCommunProb(cellchat1, raw.use = T)#默认计算方式为#type = "truncatedMean",
cellchat2 <- computeCommunProb(cellchat2, raw.use = T)

##去掉通讯数量很少的细胞，默认cutoff的值为20%，即表达比例在25%以下的基因会被认为是0， trim = 0.1可以调整比例阈值
cellchat1 <- filterCommunication(cellchat1, min.cells = 10)
cellchat2 <- filterCommunication(cellchat2, min.cells = 10)

#将推断的结果提取出来
df.net1 <- subsetCommunication(cellchat1)#将细胞通讯预测结果以数据框的形式取出
df.net2 <- subsetCommunication(cellchat2)#将细胞通讯预测结果以数据框的形式取出

# 在细胞通路水平上计算细胞间通讯
cellchat1 <- computeCommunProbPathway(cellchat1)
cellchat2 <- computeCommunProbPathway(cellchat2)

##计算细胞-细胞聚合通信网络
cellchat1 <- aggregateNet(cellchat1)
cellchat2 <- aggregateNet(cellchat2)

saveRDS(cellchat1, file = "cellchat1_NR.rds")
saveRDS(cellchat2, file = "cellchat2_R.rds")

####合并两组数据
rm(list=ls()) #清空所有变量
options(stringsAsFactors = F) #输入数据不自动转换成因子（防止数据格式错误）

cellchat1<-readRDS("cellchat1_NR.rds")
cellchat1<-updateCellChat(cellchat1)

cellchat2<-readRDS("cellchat2_R.rds")
cellchat2<-updateCellChat(cellchat2)

object.list <- list(NR = cellchat1, R = cellchat2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cellchat

#####比较交互总数和交互强度

gg1 <- compareInteractions(cellchat, show.legend = T, group = c("NR","R"))
gg2 <- compareInteractions(cellchat, show.legend = T, group = c("NR","R"), measure = "weight")
gg1 + gg2


######两个数据集之间的细胞-细胞通信网络中相互作用的差异数量或差异作用强度可以使用圆图来可视化，
##其中红色的（或者蓝色的) 彩色边缘代表与第一个数据集相比第二个数据集中的信号增加（或者减少) 

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

######heatmap

gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2


#####cellchat plot

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


######特定细胞类型之间的通讯-互作数

group.cellType <- c(rep("Myeloid cell", 3), rep("Epithelial cell", 3), rep("Fibrocyte",3))
group.cellType <- factor(group.cellType, levels = c("Myeloid cell", "Epithelial cell", "Fibrocyte"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.cellType)})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

####比较 2D 空间中的传出和传入交互强度可以很容易地识别在不同数据集之间发送或接收信号发生显着变化的细胞群。
dev.new(width = 800, height = 600)
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  object.list[[i]] <- netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP")
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)


######## 识别与一个细胞组相关的信号变化

gg1 <- netAnalysis_signalingChanges_scatter(object.list, idents.use = "Myeloid cell")
gg2 <- netAnalysis_signalingChanges_scatter(object.list, idents.use = "Epithelial cell")


######比较每个信号通路的整体信号流

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

#可以通过简单地比较每个信号通路的信息流来识别保守和特定于上下文的信号通路，
#该信息流由推断网络中所有细胞组对之间的通信概率之和（即网络中的总权重）定义）。
#此条形图可以以堆叠模式绘制，也可以不绘制。根据NR和 R之间推断网络内的整体信息流的差异
#对重要的信号通路进行排序。红色的顶部信号通路在NR中富集，而这些绿色在R中富集。

##########outgoing,

library(ComplexHeatmap)

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 10)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

########incoming

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 10, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 10, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

########overall

ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 10, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 10, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


#########识别上调和下调的信号配体-受体对

netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:4,6:8),  comparison = c(1, 2), angle.x = 45)

######可视化某一信号
dev.new(width = 800, height = 600)
pathways.show <- c("TGFb") 

netVisual_aggregate(object.list[[1]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[1]))
netVisual_aggregate(object.list[[2]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[2]))

###R组信号通路
dev.new(width = 800, height = 600)
pathways.show <- c("PTN") 
netVisual_aggregate(object.list[[2]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list[[2]])))


##
saveRDS(cellchat, file = "cellchat_merge.rds")
saveRDS(object.list, file = "object.list.rds")
