
library(BiocManager)
library(parallel)
library(randomForestSRC)
#devtools::install_github("l-magnificence/Mime")
library(Mime1)
library(data.table)
library(dplyr)
library(stringr)

setwd("D:/workplace/Immunotherapy/datasets/")

# 所要研究的gene列表
gene<-fread("D:/workplace/Immunotherapy/datasets/M1M2/diff_sig13.txt",header = F)%>%as.data.frame()
gene<-gene$V1

# 输入免疫治疗肺癌数据集 --------------------------------------------------------------------

# 读取文件，将其存储到数据框中，保留行名和列名


# ###GSE93157 -------------------------------------------------------------

setwd("D:/workplace/Immunotherapy/datasets/GSE93157/")

sam<-fread("GSE93157_info.txt")%>%as.data.frame()

sam_t<-t(sam)%>%as.data.frame()
sam_t<-sam_t[-1,]

sam2<-subset(sam_t,sam_t$V7%in%c("LUNG NON-SQUAMOUS CANCER",
                     "SQUAMOUS LUNG CANCER"))
##n=35
sam_df=data.frame(ID=sam2$V1,
                  Var=sam2$V20)

sam_df$Var<-str_remove(sam_df$Var,"best.resp ")
sam_df$Var<-ifelse(sam_df$Var=="PD","N","Y")
###N:PD=14,Y:NPD=21

####expression

exp<-fread("GSE93157_series_matrix.txt")%>%as.data.frame()
rownames(exp)<-exp$ID_REF
exp<-exp[,-1]

####判断是否log2处理
ex <- exp
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
exprSet <- log2(ex+1)
print("log2 transform finished")}else{print("log2 transform not needed")}

####整合样本和表达

exp_t<-t(exp)%>%as.data.frame()

exp2<-exp_t[sam_df$ID,]
identical(rownames(exp2),sam_df$ID)

exp_gene<-intersect(gene,colnames(exp2))
exp_df<-exp2[,exp_gene]

dataset1<-cbind(sam_df,exp_df)


# GSE111414 ---------------------------------------------------------------

setwd("D:/workplace/Immunotherapy/datasets/GSE111414/")

gse1<-fread("GSE111414_info.txt",stringsAsFactors=FALSE)%>%as.data.frame()
#(PBMCs) from metastatic NSCLC patients with PD-1 blockade
rownames(gse1)<-1:length(gse1$`Patient BS373 timepoint N1`)
gse1<-gse1[,-1]
colnames(gse1)<-gse1[1,]

N1<-seq(1,19,2)
N2<-seq(2,20,2)

df_n1<-gse1[c(7,11,14),N1]
##before the first administration (N1)
df_n2<-gse1[c(7,11,14),N2]
##during treatment (at 4 weeks, N2)

#####TIME:N1
df_n1[3,]%>%as.character()%>%table()
#non-responder(PD)=5,responder(PR)=5
df_n1[3,]<-str_extract(df_n1[3,],"[PD,PR]+")

#####TIME:N2
df_n2[3,]%>%as.character()%>%table()
#non-responder(PD)=5,responder(PR)=5
df_n2[3,]<-str_extract(df_n2[3,],"[PD,PR]+")

df_n20<-cbind(df_n1,df_n2)

sam_df2<-data.frame(ID=colnames(df_n20),
                    Var=as.character(df_n20[1,]))

sam_df2$Var<-str_extract(sam_df2$Var,"[a-z]*")
sam_df2$Var<-ifelse(sam_df2$Var=="responder","Y","N")

# 表达谱 ---------------------------------------------------------------------

exp<-fread("GSE111414_gene_counts.csv")%>%as.data.frame()

# 表达谱预处理count转TPM  --------------------------------------------------------

###多个探针均值作为基因表达值
exp5<-aggregate(exp,by=list(exp$SYMBOL),FUN=mean,na.rm=T) 
exp5<-exp5[,-2]
rownames(exp5)<-exp5$Group.1
exp5<-exp5[,-1]
##count matrix

#expM<-ifelse(expMatrix<0,0,expMatrix) 
expMatrix<-exp5

# 计算基因长度
# 基因的有效长度
eff_length2 <- read.csv("eff_length_symbol.csv", row.names = 1, header = T)
eff_length2$gene_id <- rownames(eff_length2) #基因长度

# read count转TPM
# 从输入数据里提取基因名
feature_ids <- rownames(expMatrix)

# 检查gtf文件和表达量输入文件里基因名的一致性
if (!all(feature_ids %in% rownames(eff_length2))) {
  tbl <- table(feature_ids %in% rownames(eff_length2))
  msg1 <- sprintf("%i gene is shared, %i gene is specified", tbl[[2]], tbl[[1]])
  warning(msg1)
}

if (!identical(feature_ids, rownames(eff_length2))) {
  msg2 <- sprintf("Given GTF file only contain %i gene, but experssion matrix has %i gene", nrow(eff_length2), nrow(expMatrix))
  warning(msg2)
}

expMatrix <- expMatrix[feature_ids %in% rownames(eff_length2), ]
mm <- match(rownames(expMatrix), rownames(eff_length2))
eff_length2 <- eff_length2[mm, ]

if (identical(rownames(eff_length2), rownames(expMatrix))) {
  print("GTF and expression matix now have the same gene and gene in same order")
}

x <- expMatrix / eff_length2$eff_length
expMatrix_tpm <- t(t(x) / colSums(x)) * 1e6
# 查看前三个基因的TPM值
expMatrix_tpm[1:3, 1:5]

log_tpm<-log2(expMatrix_tpm+1)
###log2(TPM+1)

exp_gene2<-intersect(rownames(log_tpm),gene)
exp1<-log_tpm[exp_gene2,]

exp_df2<-t(exp1)%>%as.data.frame()

####对应样本名称

s1<-str_extract(gse1[9,],"S[0-9]*")%>%str_remove("S")
s2<-str_extract(gse1[15,],"N[0-9]*")

s12<-str_c(s1,s2,sep="_")
identical(s12,rownames(exp_df2))
rownames(exp_df2)<-colnames(gse1)
exp_df2<-exp_df2[sam_df2$ID,]
identical(sam_df2$ID,rownames(exp_df2))

dataset2<-cbind(sam_df2,exp_df2)

# ###GSE126044 -------------------------------------------------------------

setwd("D:/workplace/Immunotherapy/datasets/GSE126044/")

sam<-fread("GSE126044_series_matrix.txt")%>%as.data.frame()

sam_t<-t(sam)%>%as.data.frame()
sam_t<-sam_t[-1,]

##n=16
sam_df=data.frame(ID=sam_t$V1,
                  Var=sam_t$V11)

sam_df$Var<-str_remove(sam_df$Var,"patient response: ")
sam_df$Var<-ifelse(sam_df$Var=="non-responder","N","Y")
###N:=11,Y:=5

####expression

exp<-fread("GSE126044_counts.txt")%>%as.data.frame()

###多个探针均值作为基因表达值
exp5<-aggregate(exp,by=list(exp$V1),FUN=mean,na.rm=T) 
exp5<-exp5[,-2]
rownames(exp5)<-exp5$Group.1
exp5<-exp5[,-1]
##count matrix

#expM<-ifelse(expMatrix<0,0,expMatrix) 
expMatrix<-exp5

# 计算基因长度
# 基因的有效长度
eff_length2 <- read.csv("eff_length_symbol.csv", row.names = 1, header = T)
eff_length2$gene_id <- rownames(eff_length2) #基因长度

# read count转TPM
# 从输入数据里提取基因名
feature_ids <- rownames(expMatrix)

# 检查gtf文件和表达量输入文件里基因名的一致性
if (!all(feature_ids %in% rownames(eff_length2))) {
  tbl <- table(feature_ids %in% rownames(eff_length2))
  msg1 <- sprintf("%i gene is shared, %i gene is specified", tbl[[2]], tbl[[1]])
  warning(msg1)
}

if (!identical(feature_ids, rownames(eff_length2))) {
  msg2 <- sprintf("Given GTF file only contain %i gene, but experssion matrix has %i gene", nrow(eff_length2), nrow(expMatrix))
  warning(msg2)
}

expMatrix <- expMatrix[feature_ids %in% rownames(eff_length2), ]
mm <- match(rownames(expMatrix), rownames(eff_length2))
eff_length2 <- eff_length2[mm, ]

if (identical(rownames(eff_length2), rownames(expMatrix))) {
  print("GTF and expression matix now have the same gene and gene in same order")
}

x <- expMatrix / eff_length2$eff_length
expMatrix_tpm <- t(t(x) / colSums(x)) * 1e6
# 查看前三个基因的TPM值
expMatrix_tpm[1:3, 1:5]

log_tpm<-log2(expMatrix_tpm+1)
###log2(TPM+1)

exp_gene3<-intersect(rownames(log_tpm),gene)
exp3<-log_tpm[exp_gene3,]

rownames(sam_t)<-str_remove(rownames(sam_t),"RNA-seq_")

exp3<-exp3[,rownames(sam_t)]
identical(colnames(exp3),rownames(sam_t))

colnames(exp3)<-sam_t$V1
exp_df3<-t(exp3)%>%as.data.frame()

dataset3<-cbind(sam_df,exp_df3)

# ###GSE135222 -------------------------------------------------------------

setwd("D:/workplace/Immunotherapy/datasets/GSE135222/")

sam<-fread("GSE135222_info.txt")%>%as.data.frame()

sam_t<-t(sam)%>%as.data.frame()
sam_t<-sam_t[-1,]

##n=27
sam_df=data.frame(ID=sam_t$V1,
                  Var=sam_t$V13)

sam_df$Var<-str_remove(sam_df$Var,"pfs.time: ")%>%as.numeric()
sam_df$Var<-ifelse(sam_df$Var<=180,"N","Y")
###N:=20,Y:=7

####expression
library(clusterProfiler)
library(org.Hs.eg.db)

exp<-fread("GSE135222_GEO_RNA-seq_omicslab_exp.tsv")%>%as.data.frame()

exp$gene_id<-str_extract(exp$gene_id,"[0-9A-Z]*")

genes<-bitr(exp$gene_id,"ENSEMBL","SYMBOL","org.Hs.eg.db")

colnames(exp)[1]<-"ENSEMBL"

exp4<-inner_join(genes,exp)
exp4<-exp4[,-1]

###多个探针均值作为基因表达值
exp5<-aggregate(exp4,by=list(exp4$SYMBOL),FUN=mean,na.rm=T) 
exp5<-exp5[,-2]
rownames(exp5)<-exp5$Group.1
exp5<-exp5[,-1]
##count matrix

#expM<-ifelse(expMatrix<0,0,expMatrix) 
expMatrix<-exp5

# 计算基因长度
# 基因的有效长度
eff_length2 <- read.csv("eff_length_symbol.csv", row.names = 1, header = T)
eff_length2$gene_id <- rownames(eff_length2) #基因长度

# read count转TPM
# 从输入数据里提取基因名
feature_ids <- rownames(expMatrix)

# 检查gtf文件和表达量输入文件里基因名的一致性
if (!all(feature_ids %in% rownames(eff_length2))) {
  tbl <- table(feature_ids %in% rownames(eff_length2))
  msg1 <- sprintf("%i gene is shared, %i gene is specified", tbl[[2]], tbl[[1]])
  warning(msg1)
}

if (!identical(feature_ids, rownames(eff_length2))) {
  msg2 <- sprintf("Given GTF file only contain %i gene, but experssion matrix has %i gene", nrow(eff_length2), nrow(expMatrix))
  warning(msg2)
}

expMatrix <- expMatrix[feature_ids %in% rownames(eff_length2), ]
mm <- match(rownames(expMatrix), rownames(eff_length2))
eff_length2 <- eff_length2[mm, ]

if (identical(rownames(eff_length2), rownames(expMatrix))) {
  print("GTF and expression matix now have the same gene and gene in same order")
}

x <- expMatrix / eff_length2$eff_length
expMatrix_tpm <- t(t(x) / colSums(x)) * 1e6
# 查看前三个基因的TPM值
expMatrix_tpm[1:3, 1:5]

log_tpm<-log2(expMatrix_tpm+1)
###log2(TPM+1)

exp_gene4<-intersect(rownames(log_tpm),gene)
exp44<-log_tpm[exp_gene4,]

rownames(sam_t)<-str_remove(rownames(sam_t)," ")

exp44<-exp44[,rownames(sam_t)]
identical(colnames(exp44),rownames(sam_t))

colnames(exp44)<-sam_t$V1
exp_df4<-t(exp44)%>%as.data.frame()

identical(sam_df$ID,rownames(exp_df4))

dataset4<-cbind(sam_df,exp_df4)

# # 将3个数据框合并到一个列表中 --------------------------------------------------------


list_train_vali_Data <- list(training = dataset3, validation1 = dataset2)##Dataset1是训练数据集，其他数据集作为验证数据集
save(list_train_vali_Data,file = 'D:/workplace/Immunotherapy/datasets/M1M2/datasets2.Rdata')

setwd("D:/workplace/Immunotherapy/datasets/M1M2/")
load("datasets2.Rdata")

# 模型选择 --------------------------------------------------------------------
#devtools::install_github("bhklab/survcomp")
#devtools::install_github("mixOmicsTeam/mixOmics")
#devtools::install_github("binderh/CoxBoost")
#devtools::install_github("souravc83/fastAdaboost")

library(fastAdaboost)
library(cancerclass)
library(Mime1)
library(mixOmics)
library(survcomp)
library(CoxBoost)

res.ici2 <- ML.Dev.Pred.Category.Sig(train_data = list_train_vali_Data$training,
                                    list_train_vali_Data = list_train_vali_Data,
                                    candidate_genes = gene,
                                    methods = c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass'),
                                    seed = 5201314,
                                    cores_for_parallel = 16
)

####绘制不同方法在不同数据集之间的AUC
auc_vis_category_all(res.ici2,dataset = c("training","validation1"),
                     order= c("training","validation1"))


####绘制不同数据集间具体方法的ROC曲线
plot_list<-list()
methods <- c('nb','svmRadialWeights','rf','kknn','adaboost','LogitBoost','cancerclass')
for (i in methods) {
  plot_list[[i]]<-roc_vis_category(res.ici2,model_name = i,dataset = c("training","validation1"),
                                   order= c("training","validation1"),
                                   anno_position=c(0.4,0.25))
}
aplot::plot_list(gglist=plot_list,ncol=3)


####
res.ici2$sig.gene

write.table(res.ici2$sig.gene,
            "D:/workplace/Immunotherapy/datasets/M1M2/model_sig6.txt",
            row.names = F,
            col.names = F,
            quote = F)

