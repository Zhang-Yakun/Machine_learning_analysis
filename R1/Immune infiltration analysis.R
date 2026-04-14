
library(data.table)
library(dplyr)
library(stringr)

# 输入4个数据集 -----------------------------------------------------------------

###GSE93157

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

exp<-fread("GSE93157_raw_data_values.txt")%>%as.data.frame()
rownames(exp)<-exp$ID_REF
exp<-exp[,-1]
exp2<-exp[,-(66:69)]

####整合样本和表达

colnames(exp2)<-sam_t$V1

exp_count<-exp2[,sam_df$ID]

###COUNTS 转 TPM

#expM<-ifelse(expMatrix<0,0,expMatrix) 
expMatrix<-exp_count

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
expMatrix_tpm[1:3,1:5]
## 17557 * 1129

tpm_df<-as.data.frame(expMatrix_tpm)
tpm_df$GENE<-rownames(tpm_df)

tpm_df2<-tpm_df%>%select("GENE",everything())

write.table(tpm_df2,"DATA1.txt",
            row.names = F, quote = FALSE,sep = "\t")
write.csv(sam_df,"data1_group.csv")

####GSE111414

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

tpm_df<-as.data.frame(expMatrix_tpm)

s1<-str_extract(gse1[9,],"S[0-9]*")%>%str_remove("S")
s2<-str_extract(gse1[15,],"N[0-9]*")

s12<-str_c(s1,s2,sep="_")
identical(s12,colnames(tpm_df))
colnames(tpm_df)<-colnames(gse1)
tpm_df<-tpm_df[,sam_df2$ID]
identical(sam_df2$ID,colnames(tpm_df))

tpm_df$GENE<-rownames(tpm_df)

tpm_df2<-tpm_df%>%select("GENE",everything())

write.table(tpm_df2,"DATA2.txt",
            row.names = F, quote = FALSE,sep = "\t")
write.csv(sam_df2,"data2_group.csv")

####GSE126044

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

rownames(sam_t)<-str_remove(rownames(sam_t),"RNA-seq_")

expMatrix_tpm<-expMatrix_tpm[,rownames(sam_t)]
identical(colnames(expMatrix_tpm),rownames(sam_t))

colnames(expMatrix_tpm)<-sam_t$V1
tpm_df<-expMatrix_tpm%>%as.data.frame()

tpm_df$GENE<-rownames(tpm_df)

tpm_df2<-tpm_df%>%select("GENE",everything())

write.table(tpm_df2,"DATA3.txt",
            row.names = F, quote = FALSE,sep = "\t")
write.csv(sam_df,"data3_group.csv")

####GSE135222

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

rownames(sam_t)<-str_remove(rownames(sam_t)," ")

expMatrix_tpm<-expMatrix_tpm[,rownames(sam_t)]
identical(colnames(expMatrix_tpm),rownames(sam_t))

colnames(expMatrix_tpm)<-sam_t$V1
tpm_df<-t(expMatrix_tpm)%>%as.data.frame()

identical(sam_df$ID,rownames(tpm_df))

tpm_df<-t(tpm_df)%>%as.data.frame()

tpm_df$GENE<-rownames(tpm_df)

tpm_df2<-tpm_df%>%dplyr::select("GENE",everything())

write.table(tpm_df2,"DATA4.txt",
            row.names = F, quote = FALSE,sep = "\t")
write.csv(sam_df,"data4_group.csv")


# 四个GEO数据分别免疫浸润 ------------------------------------------------------

setwd("D:/workplace/Immunotherapy/datasets/cibersort/")

library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)

#####GSE93157

source('Cibersort.R')

result1 <- CIBERSORT('LM22.txt','DATA1.txt', perm = 1000, QN = T)
#perm置换次数=1000，QN分位数归一化=TRUE

result2 <- CIBERSORT('LM22.txt','DATA2.txt', perm = 1000, QN = T)
#perm置换次数=1000，QN分位数归一化=TRUE

result3 <- CIBERSORT('LM22.txt','DATA3.txt', perm = 1000, QN = T)
#perm置换次数=1000，QN分位数归一化=TRUE

result4 <- CIBERSORT('LM22.txt','DATA4.txt', perm = 1000, QN = T)
#perm置换次数=1000，QN分位数归一化=TRUE

# M1,M2细胞比例 

cell_exp1<-fread("CIBERSORT-Results4.txt")%>%as.data.frame()

M1_M2<-data.frame(Macrophages_M1=cell_exp1$`Macrophages M1`,
                  Macrophages_M2=cell_exp1$`Macrophages M2`)
rownames(M1_M2)<-cell_exp1$Mixture

# 比较M1/M2比例在免疫治疗应答组间差异 ----------------------------------------------------

group<-fread("data4_group.csv")%>%as.data.frame()

identical(group$ID,rownames(M1_M2))
M1_M2$group<-group$Var%>%as.factor()
M1_M2$ratio<-M1_M2$Macrophages_M1/M1_M2$Macrophages_M2

t.test(ratio~group,data=M1_M2)
#p=0.3


