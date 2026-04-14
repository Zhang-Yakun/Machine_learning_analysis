rm(list = ls())

library(dplyr)
library(data.table)

setwd("D:/workplace/Immunotherapy/datasets/cibersort/")

# 细胞比例分组比较 -------------------------------------------------------------------

cell1<-fread("CIBERSORT-Results3.txt")%>%as.data.frame()
cell1<-cell1[,1:23]
colnames(cell1)[1]<-"ID"
cell1$M1M2<-ifelse(cell1$`Macrophages M1`>cell1$`Macrophages M2`,1,0)

group1<-fread("data3_group.csv")%>%as.data.frame()
input1<-inner_join(group1,cell1)
rownames(input1)<-input1$ID
input1<-input1[,-1]

### 4套数据集免疫细胞R,NR比较分析
library("ggsignif")
library("ggpubr")

cell.pro.Paired<-data.frame()

for(i in 1:(dim(input1)[2]-1)){
  cell.Paired <- data.frame(group=input1$Var,
                            celltype=rep(colnames(input1)[i+1],dim(input1)[1]),
                            proportion=input1[,i+1])
  cell.pro.Paired<-rbind(cell.pro.Paired,cell.Paired)
}

head(cell.pro.Paired)
compare_means(proportion ~ group, data = cell.pro.Paired, group.by = "celltype")
options(repr.plot.height = 6, repr.plot.width = 7)
ggboxplot(cell.pro.Paired, x = "celltype", y = "proportion",
          color = "group", palette = "jco", 
          add = "jitter")+# palette可以按照期刊选择相应的配色，如"npg"等
  stat_compare_means(aes(group = group), label = "p.signif")+
  theme(panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_blank(),
        text=element_text(size = 12, family = "serif"),
        axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))


# 计算M1/M2比例 ---------------------------------------------------------------

MM_df<-input1[,c(1,16,17)]
MM_df$ratio<-(MM_df$`Macrophages M1`)/(MM_df$`Macrophages M2`)

MM_df$ratio_group<-ifelse(MM_df$ratio>1,"high_ratio","low_ratio")
table(MM_df$ratio_group)
table(MM_df$Var)


# 细胞比例分组差异分析 --------------------------------------------------------------

###读取基因表达谱

data1<-fread("DATA3.txt")%>%as.data.frame()
#TPM值
rownames(data1)<-data1$GENE
data1<-data1[,-1]

###log2(TPM+1)
log_tpm<-log2(data1+1)


###样本对应

identical(rownames(MM_df),colnames(log_tpm))
group<-MM_df$ratio_group %>% as.factor()

#####limma

library(limma)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(log_tpm , design)
contrast.matrix <- makeContrasts(high_ratio - low_ratio,
                                 levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

tempOutput <- topTable(fit2, coef = 1, n = Inf)
nrDEG1 <- na.omit(tempOutput)  #P值和logFC
nrDEG1$gene <-rownames(nrDEG1)

dif1<-subset(nrDEG1,nrDEG1$P.Value<0.05)
dif2<-subset(dif1,abs(dif1$logFC)>2)

gene_ratio<-unique(dif2$gene)
#645
# R-NR分组差异分析 --------------------------------------------------------------

###样本对应

identical(rownames(MM_df),colnames(log_tpm))
group2<-factor(MM_df$Var,levels = c("Y","N"))

#####limma

library(limma)
design2 <- model.matrix(~ 0 + group2)
colnames(design2) <- levels(group2)
fit <- lmFit(log_tpm , design2)
contrast.matrix <- makeContrasts(Y - N,
                                 levels = design2
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

tempOutput <- topTable(fit2, coef = 1, n = Inf)
nrDEG1 <- na.omit(tempOutput)  #P值和logFC
nrDEG1$gene <-rownames(nrDEG1)

dif1<-subset(nrDEG1,nrDEG1$P.Value<0.05)
dif2<-subset(dif1,abs(dif1$logFC)>2)

gene_response<-unique(dif2$gene)
#179

# 两组差异基因的交集 ---------------------------------------------------------------

diff_sig1<-intersect(gene_response,gene_ratio)
#13
diff_sig3<-intersect(gene_response,gene_ratio)
#47
diff_sig4<-intersect(gene_response,gene_ratio)
#6

write.table(diff_sig1,
            "D:/workplace/Immunotherapy/datasets/M1M2/diff_sig13.txt",
            row.names = F,
            col.names = F,
            quote = F)

# 韦恩图 ---------------------------------------------------------------------

deg <- list()

deg$response_sig<-gene_response
deg$ratio_sig<-gene_ratio

library(VennDiagram)
library(grid)
library(futile.logger)

#自定义颜色
color=c("#6095CE", "#E87D72")
#读入作图文件

p1<-venn.diagram(
  x = deg[1:2],
  category.names = names(deg)[1:2],
  filename = NULL,
  output=TRUE,
  fill = color[1:(length(deg[1:2]))],
  col = 'black', cex = 1,
  cat.col = rep('black', 2),
  fontface = "bold"
)

pdf("D:/workplace/Immunotherapy/datasets/M1M2/dataset2_venn.pdf")
grid.draw(p1)
dev.off()


# 差异基因热图 ------------------------------------------------------------------

model_gene<-fread("D:/workplace/Immunotherapy/datasets/M1M2/diff_sig13.txt",header = F)%>%as.data.frame()

mg<-model_gene$V1

heat_input<-log_tpm[mg,]

group_df<-data.frame(sample=rownames(MM_df),
                     group=MM_df$Var)
group_df<-group_df[order(group_df$group),]

heat_input<-heat_input[,group_df$sample]
identical(colnames(heat_input),group_df$sample)

#####可视化

anno_col = data.frame(
  Group = as.factor(group_df$group)
)

rownames(anno_col) = colnames(heat_input)

n <- t(scale(t(heat_input)))
n[n > 2] <- 2
n[n < -2] <- -2
df <- n
rownames(df)<-rownames(heat_input)

library(pheatmap)
pheatmap(heat_input,cellwidth =12, cellheight = 12, fontsize =10,fontsize_row=10,
         method="pearson", #计算gene或sample之间的相关性的方法，可选"pearson" (default), "kendall", or "spearman"
         scale="row", #为基因做scale
         cluster_rows=T,#为基因做聚类
         cluster_cols=F,#为sample做聚类
         show_colnames=F,show_rownames =T,
         annotation_col = anno_col,
         treeheight_col = "0",#不画树
         border_color = "NA")

# 通路富集图 -------------------------------------------------------------------

# ========== 第一部分：安装和加载必要的包 ==========
# 安装必要的包（如果尚未安装）
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", 
                       "DOSE", "ggplot2", "igraph", "ggraph", "tidygraph",
                       "stringr", "dplyr", "RColorBrewer")

for (pkg in required_packages) {
  if (!require(pkg, quietly = TRUE, character.only = TRUE)) {
    if (pkg %in% c("clusterProfiler", "org.Hs.eg.db", "DOSE", "enrichplot")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  } else {
    library(pkg, character.only = TRUE)
  }
}

# ========== 第二部分：准备基因列表 ==========
# 输入您的13个基因
gene_symbols <- c("CD8A", "GMFG", "AOAH", "ARHGAP9", "LCP1", "CD3E", 
                  "CSF2RB", "CD3D", "GZMB", "HOXC9", "SLAMF7", 
                  "HIST1H4F", "RIMS2")

# 将基因符号转换为Entrez ID（KEGG分析需要）
gene_entrez <- bitr(gene_symbols, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = "org.Hs.eg.db")

# 检查转换结果
cat("成功转换的基因数量:", nrow(gene_entrez), "/", length(gene_symbols), "\n")
print(gene_entrez)

# ========== 第三部分：KEGG通路富集分析 ==========
# 执行KEGG富集分析
KEGG <- read.csv("D:/workplace/Immunotherapy/HSA_KEGG.csv")

kegg_enrich <- 
TERM2GENE <- KEGG[, c("KEGGID", "ENTREZID")]
TERM2NAME <- KEGG[, c("KEGGID", "DESCRIPTION")]

yy <- enricher(gene_entrez$ENTREZID, TERM2GENE = TERM2GENE, TERM2NAME = TERM2NAME, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 1)
head(yy)
yy1 <- setReadable(yy, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
head(yy1)
yy2 <- yy1[yy1$pvalue < 0.05, asis = T]

kegg_enrich<-yy2

  # ========== 第四部分：富集分析可视化 ==========
  
  # 1. 条形图：显示top富集通路
  p1 <- barplot(yy2, 
                showCategory = 10, 
                title = "Top 10 显著富集的KEGG通路",
                font.size = 10,
                col = "p.adjust") + 
    theme(axis.text.y = element_text(size = 10))+
  scale_fill_gradient(low = "#E09F87",high = "#A3BAE2")
c("lightblue", "mistyrose",
  "lightcyan", "lavender")
  # 2. 点图：显示基因计数和p值
  p2 <- dotplot(yy2, 
                showCategory = 10, 
                title = "KEGG通路富集分析",
                color = "p.adjust",
                size = "Count") +
    scale_color_gradient(low = "#FF9999", high = "#99CCFF")

  # 3. 基因-通路网络图（主要可视化）
  # 先提取前10个最显著的通路
  if (nrow(kegg_enrich) > 10) {
    top_pathways <- kegg_enrich[1:10, ]
  } else {
    top_pathways <- kegg_enrich
  }
  
  # 准备网络图数据
  net_data <- as.data.frame(top_pathways)
  
  if (nrow(net_data) > 0) {
    # 创建基因-通路关系数据
    pathway_genes <- strsplit(net_data$geneID, "/")
    
    # 构建边列表
    edges <- data.frame()
    for (i in 1:nrow(net_data)) {
      pathway <- net_data$Description[i]
      genes <- pathway_genes[[i]]
      for (gene in genes) {
        edges <- rbind(edges, data.frame(
          from = pathway,
          to = gene,
          pvalue = -log10(net_data$p.adjust[i])
        ))
      }
    }
    
    # 创建节点列表
    pathway_nodes <- data.frame(
      name = net_data$Description,
      type = "pathway",
      size = net_data$Count,
      pvalue = -log10(net_data$p.adjust)
    )
    
    gene_nodes <- data.frame(
      name = unique(unlist(pathway_genes)),
      type = "gene",
      size = 3,
      pvalue = 1
    )
    
    nodes <- rbind(pathway_nodes, gene_nodes)
    
    # 创建igraph对象
    g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
    
    # 设置节点属性
    V(g)$color <- ifelse(V(g)$type == "pathway", 
                         "#FF9999",  # 通路节点颜色
                         "#99CCFF")  # 基因节点颜色
    
    V(g)$size <- ifelse(V(g)$type == "pathway",
                        V(g)$size * 2,  # 通路节点大小
                        5)  # 基因节点大小
    
    # 设置边属性
    E(g)$width <- edges$pvalue / max(edges$pvalue) * 3
    
    # 可视化网络
    set.seed(123)  # 确保布局可重现
    
    # 方法1：使用igraph的基本绘图
    pdf("KEGG_pathway_network_igraph.pdf", width = 12, height = 10)
    
    plot(g,
         layout = layout_with_fr(g),  # Fruchterman-Reingold布局
         vertex.label = ifelse(V(g)$type == "pathway", 
                               V(g)$name, 
                               NA),  # 只显示通路名称
         vertex.label.cex = ifelse(V(g)$type == "pathway", 0.7, 0),
         vertex.label.color = "black",
         vertex.frame.color = NA,
         edge.color = "gray70",
         main = "KEGG通路-基因网络图")
    
    # 添加图例
    legend("bottomright", 
           legend = c("通路", "基因"), 
           pch = 21,
           col = c("#FF9999", "#99CCFF"),
           pt.bg = c("#FF9999", "#99CCFF"),
           pt.cex = 2,
           cex = 0.8)
    
    dev.off()
    
    # 方法2：使用ggraph绘制更美观的网络图
    if (require("ggraph") && require("tidygraph")) {
      library(ggraph)
      library(tidygraph)
      
      # 转换为tidygraph对象
      tg <- as_tbl_graph(g)
      
      p3 <- ggraph(tg, layout = "fr") + 
        geom_edge_link(aes(width = width), 
                       alpha = 0.5, 
                       color = "gray70") + 
        geom_node_point(aes(color = type, size = size), 
                        alpha = 0.8) + 
        geom_node_text(aes(label = ifelse(type == "pathway", name, "")), 
                       size = 3, 
                       repel = TRUE,
                       max.overlaps = 20) +
        scale_color_manual(values = c("pathway" = "#FF9999", "gene" = "#99CCFF"),
                           labels = c("通路", "基因")) +
        scale_size(range = c(3, 10)) +
        theme_void() +
        theme(
          legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
        ) +
        labs(title = "KEGG通路-基因关联网络",
             color = "节点类型",
             size = "节点大小")
      
      ggsave("KEGG_pathway_network_ggraph.pdf", p3, width = 12, height = 10)
      ggsave("KEGG_pathway_network_ggraph.png", p3, width = 12, height = 10, dpi = 300)
    }
    
    # 4. 富集结果的热图
    # 提取基因在通路中的存在矩阵
    gene_pathway_matrix <- matrix(0, 
                                  nrow = length(unique(unlist(pathway_genes))), 
                                  ncol = nrow(net_data))
    
    rownames(gene_pathway_matrix) <- unique(unlist(pathway_genes))
    colnames(gene_pathway_matrix) <- net_data$Description
    
    for (i in 1:nrow(net_data)) {
      genes_in_pathway <- pathway_genes[[i]]
      gene_pathway_matrix[genes_in_pathway, i] <- 1
    }
    
    # 只保留至少在一个通路中出现的基因
    gene_pathway_matrix <- gene_pathway_matrix[rowSums(gene_pathway_matrix) > 0, ]
    
    if (nrow(gene_pathway_matrix) > 1 && ncol(gene_pathway_matrix) > 1) {
      # 绘制热图
      pheatmap::pheatmap(gene_pathway_matrix,
                         main = "基因-通路关联热图",
                         color = c("white", "darkblue"),
                         cluster_rows = TRUE,
                         cluster_cols = TRUE,
                         show_rownames = TRUE,
                         fontsize_row = 8,
                         filename = "gene_pathway_heatmap.pdf")
    }
    
    cat("\n======= 可视化文件已生成 =======\n")
    cat("1. KEGG_pathway_network_igraph.pdf - 通路网络图(igraph版)\n")
    cat("2. KEGG_pathway_network_ggraph.pdf - 通路网络图(ggraph版)\n")
    cat("3. gene_pathway_heatmap.pdf - 基因-通路关联热图\n")
    
  } else {
    cat("\n没有足够的通路数据生成网络图。\n")
  }
  
  # 显示所有图形
  print(p1)
  print(p2)
  if (exists("p3")) print(p3)
  
} else {
  cat("\n未发现显著富集的KEGG通路。\n")
  cat("可能原因：\n")
  cat("1. 基因数量太少\n")
  cat("2. 这些基因不共同参与已知的KEGG通路\n")
  cat("3. 需要调整富集分析的参数\n")
}

library("enrichplot")
cnetplot(yy2)


# ========== 第五部分：附加功能 - GO富集分析 ==========
# 如果需要，也可以进行GO富集分析
go_enrich <- enrichGO(gene = gene_entrez$ENTREZID,
                      OrgDb = org.Hs.eg.db,
                      keyType = "ENTREZID",
                      ont = "BP",  # 生物过程
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = TRUE)

if (!is.null(go_enrich) && nrow(go_enrich) > 0) {
  cat("\n======= GO富集分析结果 =======\n")
  print(head(go_enrich, 10))
  write.csv(as.data.frame(go_enrich), "GO_enrichment_results.csv", row.names = FALSE)
  
  # 绘制GO富集结果图
  p_go <- barplot(go_enrich, showCategory = 10, title = "Top 10 GO生物过程")
  print(p_go)
  ggsave("GO_enrichment_barplot.pdf", p_go, width = 10, height = 8)
}



# ... ---------------------------------------------------------------------

# ========== 安装和加载必要包 ==========
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c("clusterProfiler", "org.Hs.eg.db", "pheatmap", 
                       "RColorBrewer", "ggplot2", "dplyr")

for (pkg in required_packages) {
  if (!require(pkg, quietly = TRUE, character.only = TRUE)) {
    if (pkg %in% c("clusterProfiler", "org.Hs.eg.db")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# ========== 准备基因列表并进行富集分析 ==========
gene_symbols <- c("CD8A", "GMFG", "AOAH", "ARHGAP9", "LCP1", "CD3E", 
                  "CSF2RB", "CD3D", "GZMB", "HOXC9", "SLAMF7", 
                  "HIST1H4F", "RIMS2")

# 转换为Entrez ID
gene_entrez <- bitr(gene_symbols, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = "org.Hs.eg.db")

# KEGG富集分析
kegg_result<-yy2



# 逻辑回归计算OR值 -------------------------------------------------------------------

library(glmnet)

###提取基因表达谱
autogene<-read.table("12个自噬差异基因.txt",header = T,stringsAsFactors = F,sep = "\t")
allgene<-read.table("GSE16561_expression.txt",header = T,sep = "\t",quote = "")
autogene<-data.frame(autogene[,1])
colnames(autogene)<-c("symbol")
autoexp<-merge(autogene,allgene,by="symbol")
write.table(autoexp,"dif_auto_geneexpr.txt",quote = F,sep="\t",row.names = F)
rownames(autoexp)<-autoexp[,1]
autoexp<-autoexp[,-1]
autoexp<-as.data.frame(t(autoexp))
autoexp$sample<-row.names(autoexp)
clinical<-read.table("clinical0.txt",quote="",sep="\t",header = T)[,1:2]
library(stringr)
group_list=ifelse(str_detect(clinical$title,"Control"),"control","treat")
table(group_list)
group_list = factor(group_list,
                    levels = c("control","treat"))
clinical$status<-group_list
names(clinical)[2]<-"sample"
write.table(clinical[,2:3],"GSE16561_group.txt",sep="\t",quote=F,row.names=F)
auto_cli<-merge(clinical,autoexp,all=F)[,-2]
row.names(auto_cli)<-auto_cli$sample
auto_cli<-auto_cli[,-1]

####单变量Logistic
uni_glm_model<-function(x){#限合结局和变量
  FL<-as.formula(paste("status==","'treat'","~" ,x,sep=""))
  glm1<-glm(FL,data=auto_cli,family ="binomial")#提取所有回归结果放入g1m2中
  glm2<-summary(glm1)
  OR<-round(exp(coef(glm1)),2)
  SE<-glm2$coefficients[,2]#计算CI保留两位小数并合并
  CI5<-round(exp(coef(glm1)-1.96*SE),2)
  CI95<-round(exp(coef(glm1)+1.96*SE),2)
  CI<-paste0(CI5,"-",CI95)#提取P值
  P<-round(glm2$coefficients[,4],4)
  uni_glm_model <- data.frame("Characteristics"=x,
                              "OR" = OR,
                              "CI5"=CI5,
                              "CI95"=CI95,
                              "CI"=CI,
                              "P"=P)[-1,]
  return(uni_glm_model)}
variable.names<- colnames(auto_cli)[c(2:13)];variable.names
uni_logistic<- lapply(variable.names,uni_glm_model)#批量输出结果并合并在一起
uni_logistic
library(plyr)
uni_glm<-ldply(uni_logistic,data.frame);uni_glm
write.table(uni_glm,"uni_logistic.txt",sep="\t",quote=F,row.names=F)
###install.packages("forestplot")
library(forestplot)
pdf("forestplot_all.pdf",width = 12,height = 4)
forestplot(labeltext=as.matrix(uni_glm[,1]), #文本信息
           mean = uni_glm$OR,##HR值
           lower = as.numeric(uni_glm$CI5),##95%置信区间
           upper = as.numeric(uni_glm$CI95),#95%置信区间
           boxsize = 0.1,##大小
           graph.pos=2,#图在表中的列位置
           graphwidth = unit(0.6,"npc"),#图在表中的宽度比例
           fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石,可以更改fpDrawNormalCI；fpDrawCircleCI等
           col=fpColors(box="steelblue", lines="black", zero = "black"),#颜色设置
           lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#置信区间用线宽、高、型
           zero=1,#zero线横坐标
           lwd.zero=2,#zero线宽
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.5), xlab = gpar(cex = 0.7), cex = 0.7),
           grid=T,
           lwd.xaxis=2,#X轴线宽
           xlab="Odds ratio",#X轴标题
           clip=c(0, 50),#边界
           colgap = unit(0.8,"cm"),#列间隙
           new_page = T#是否新页
)
dev.off()
uni_glm<-uni_glm[uni_glm$P<=0.05,]

####Lasso
reg_lasso <- function(inputArr=NULL, x_cols=NULL, y_cols=NULL, out_dir=NULL, reg_family="binomial", return_fitted_values=TRUE, pred_type="link"){
  library(glmnet)
  cv.fit <- cv.glmnet(x=as.matrix(inputArr[, x_cols]), y=inputArr[, y_cols], family=reg_family)
  model_fit <- glmnet(x=as.matrix(inputArr[, x_cols]), y=inputArr[, y_cols], family=reg_family)
  # 输出建模中间结果
  if(!is.null(out_dir)){
    #checkDir(out_dir)
    pdf(paste0(out_dir, "likehood_dev.pdf"))
    plot(cv.fit)
    dev.off()
    pdf(paste0(out_dir, "feat_coef_L1_Norm.pdf"))
    plot(model_fit, label = TRUE)
    dev.off()
  }
  # 提取lambda.min时的系数
  feat_coef <- coef(object=model_fit, s=cv.fit$lambda.min)
  # 提取系数不为0的变量
  ex_idx <- which(feat_coef!=0)
  final_coef <- feat_coef[ex_idx]
  final_feat <- row.names(feat_coef)[ex_idx]
  names(final_coef) <- final_feat
  if(return_fitted_values){
    # 计算拟合值
    risk_score <- predict(model_fit, newx=as.matrix(inputArr[, x_cols]), s=cv.fit$lambda.min, type=pred_type)
    return(list(final_coef=final_coef, model_fit=model_fit, cv.fit=cv.fit, fitted_values=risk_score))
  }else{
    return(list(final_coef=final_coef, model_fit=model_fit, cv.fit=cv.fit))
  }
}
reg_lasso_res <- reg_lasso(inputArr=auto_cli, x_cols=names(auto_cli)[-1], y_cols="status", out_dir="D:/stroke autophagy/", reg_family="binomial", return_fitted_values=FALSE, pred_type="link")
reg_lasso_res$final_coef

reg_lasso_res$model_fit
reg_lasso_res$cv.fit
library(forestplot)
uni_glm<-uni_glm[c(3,4,5,6,7,8,11),]

forestplot(labeltext=as.matrix(uni_glm[,1]), #文本信息
           mean = uni_glm$OR,##HR值
           lower = as.numeric(uni_glm$CI5),##95%置信区间
           upper = as.numeric(uni_glm$CI95),#95%置信区间
           boxsize = 0.8,##大小
           graph.pos=2,#图在表中的列位置
           graphwidth = unit(0.6,"npc"),#图在表中的宽度比例
           fn.ci_norm="fpDrawDiamondCI",#box类型选择钻石,可以更改fpDrawNormalCI；fpDrawCircleCI等
           col=fpColors(box="steelblue", lines="black", zero = "black"),#颜色设置
           lwd.ci=2,ci.vertices.height = 0.1,ci.vertices=TRUE,#置信区间用线宽、高、型
           zero=1,#zero线横坐标
           lwd.zero=2,#zero线宽
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.5), xlab = gpar(cex = 0.7), cex = 0.7),
           grid=T,
           lwd.xaxis=2,#X轴线宽
           xlab="",#X轴标题
           clip=c(-Inf, 3),#边界
           colgap = unit(0.8,"cm"),#列间隙
           new_page = F#是否新页
)


##多因素logistic
FL_multi<-as.formula(paste("status==","'treat'","~" ,"EIF2AK2+GABARAPL2+HSPA8+LAMP2+PRKCQ+RAB33B+TNFSF10",sep=""))
glm1<-glm(FL_multi,data=auto_cli,family ="binomial")#提取所有回归结果放入g1m2中
glm2<-summary(glm1)
OR<-round(exp(coef(glm1)),2)
SE<-glm2$coefficients[,2]#计算CI保留两位小数并合并
CI5<-round(exp(coef(glm1)-1.96*SE),2)
CI95<-round(exp(coef(glm1)+1.96*SE),2)
CI<-paste0(CI5,"-",CI95)#提取P值
P<-round(glm2$coefficients[,4],2)
###
risk_score<-glm2$coefficients[1,1]+glm2$coefficients[2,1]*autoexp[,"EIF2AK2"]+glm2$coefficients[3,1]*autoexp[,"GABARAPL2"]+glm2$coefficients[4,1]*autoexp[,"HSPA8"]+glm2$coefficients[5,1]*autoexp[,"LAMP2"]+glm2$coefficients[6,1]*autoexp[,"PRKCQ"]+glm2$coefficients[7,1]*autoexp[,"TNFSF10"]

risk_score<-data.frame(sample=row.names(autoexp),riskscore=risk_score)
risk_score_merge<-merge(clinical[,2:3],risk_score,all=F)
###risk-score分布图
library(ggpubr)
ggboxplot(risk_score_merge, x = "status", y = "riskscore",
          color = "status",add = "jitter", palette = c( "#D51877", "#6D8F91"))+ 
  stat_compare_means(aes(group = status), method = "t.test", label = "p.format")+ # Add pairwise comparisons p-value
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust =1,size = 8)) 
ggsave("boxplot_risk_score.pdf",width = 12,height = 8)
#模型参数：data.
glm1$coefficients
#线性模型的预测数据：
glm1$linear.predictors
#vip等于1的概率#
glm1$fitted.values
#线性拟合模型的残差
glm1$residuals
