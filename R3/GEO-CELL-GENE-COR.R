setwd("D:/workplace/Immunotherapy/datasets/M1M2/")

###细胞含量
cell_exp<-fread("CIBERSORT-Results3.txt")%>%as.data.frame()
rownames(cell_exp)<-cell_exp$Mixture

cor_input<-cell_exp[,-1]
cor_input<-cor_input[,15:16]
cor_input$M1_M2_ratio<-(cor_input$`Macrophages M1`/cor_input$`Macrophages M2`)

###基因表达
sig57<-read.table("D:/workplace/Immunotherapy/datasets/M1M2/diff_sig13.txt")
sig<-sig57$V1

exp<-fread("DATA3.txt")%>%as.data.frame()
rownames(exp)<-exp$GENE

exp_sig<-exp[sig,-1]

# 提取免疫治疗有效患者 --------------------------------------------------------------

group<-fread("data3_group.csv",header = T)%>%as.data.frame()

g<-subset(group,group$Var=="Y")
g$ID

exp_sig<-exp_sig[,g$ID]
cor_input<-cor_input[g$ID,]
###基因与细胞的相关性

options(stringsAsFactors = FALSE) #禁止chr转成factor

cor_input<-cor_input[colnames(exp_sig),]
#n=1129
identical(rownames(cor_input),colnames(exp_sig))

gene <-rownames(exp_sig)
immuscore <- function(gene){
  y <- as.numeric(exp_sig[gene,])
  colnames <- colnames(cor_input)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- cor.test(as.numeric(cor_input[,x]), y , method="spearman")
    data.frame(gene=gene,immune_cells=x,cor=dd$estimate,p.value=dd$p.value )
  }))
}
##先写一个函数,输入一个基因，即可返回跟免疫基因的相关性、p值

immuscore("CD8A")
#以某基因为例，测试一下函数

res <- do.call(rbind,lapply(gene,immuscore))
head(res)
#基因与免疫细胞的相关性结果

######相关性画图
res$pstar <- ifelse(res$p.value < 0.05,
                    ifelse(res$p.value < 0.01,"**","*"),
                    "")
res$pstar[1:20]
library(ggplot2)
ggplot(res, aes(immune_cells, gene)) + 
  geom_tile(aes(fill = cor), colour = "white",size=1)+
  scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  geom_text(aes(label=pstar),col ="black",size = 5)+
  theme_minimal()+# 不要背景
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(size=10,angle = 45, hjust = 1),# 调整x轴文字
        axis.text.y = element_text(size = 10))+#调整y轴文字
  #调整legen
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))
