setwd("D:/workplace/Immunotherapy/mime_survival/")

library(dplyr)
library(data.table)

sig9<-read.table("signature_9.txt")
gene<-sig9$V1  

marker<-read.table("D:/workplace/Immunotherapy/TCGA/EMT_M1_M2_marker.txt")
marker<-marker$V1

# GSE26939 ----------------------------------------------------------------

setwd("D:/workplace/geneset/geo/GSE26939/")

# 生存数据预处理

sur<-fread("AGE.txt")%>%as.data.frame()
sur<-sur[,-1]

sample=sur[1,]

sample<-c()
OS<-c()
time<-c()

s<-c("survival_status(0='alive',1='dead'): 1","survival_status(0='alive',1='dead'): 0")

for(i in 1:dim(sur)[2]){
  
  x=sur[,i]
  sample[i]<-x[1]
  
  for(j in 1:length(x)){
    if(x[j]%in%s){
      OS[i]<-x[j]
      time[i]<-x[j-1]
    }
  }
}

sur.df<-rbind(sample,OS,time)
colnames(sur.df)<-sur.df[1,]

sur.df2<-t(sur.df)
colnames(sur.df2)<-c("sample","OS","time")
sur.df2<-as.data.frame(sur.df2)

library(stringr)

sur.df2$OS<-sur.df2$OS %>%str_remove("survival_status\\(0\\=\\'alive\\'\\,1\\=\\'dead\\'\\)\\: ")%>%as.numeric()
sur.df2$time<-sur.df2$time%>%str_remove("survival\\_months\\: ")%>%as.numeric()

sur.df3<-sur.df2[!is.na(sur.df2$OS),]
sur.df3<-sur.df3[!is.na(sur.df3$time),]
#去掉生存数据缺失值


##########处理表达谱

exp<-fread("GSE26939_series_matrix.txt")%>%as.data.frame()
#探针表达谱

#####GPL

gpl<-fread('GPL9053.txt',header=T)%>%as.data.frame()
gpl2<-gpl[,c(1,3)]

gpl3<-subset(gpl2,gpl2$ORF%in%marker)
unique(gpl3$ORF)

exp2<-subset(exp,exp$ID_REF%in%gpl3$ID)
colnames(exp2)[1]<-"ID"
exp4<-inner_join(gpl3,exp2)
exp4<-exp4[,-1]

###多个探针均值作为基因表达值
exp5<-aggregate(exp4,by=list(exp4$ORF),FUN=mean,na.rm=T) 
exp5<-exp5[,-2]
exp5<-exp5[,-2]
rownames(exp5)<-exp5$Group.1
exp5<-exp5[,-1]

exp6<-t(exp5)%>%as.data.frame()
exp6<-exp6[sur.df3$sample,]
exp6<-exp6[-1,]
sur.df3<-sur.df3[rownames(exp6),]

identical(sur.df3$sample,rownames(exp6))

####合并表达和生存

sur_exp_26939<-cbind(sur.df3,exp6)
colnames(sur_exp_26939)[1:3]<-c("ID","OS","OS.time")

sur_exp_26939<-sur_exp_26939%>%dplyr::select("ID","OS.time","OS",everything())
####n=114

e<-sur_exp_26939[,-(1:3)]
e2<-e+10
sur_exp_26939<-cbind(sur_exp_26939[,(1:3)],e2)
write.csv(sur_exp_26939,"D:/workplace/Immunotherapy/mime_survival/sur_exp_26939.csv",row.names = F)


# GSE72094 ----------------------------------------------------------------

setwd("D:/workplace/geneset/geo/GSE72094/")

###处理生存

sur<-fread("AGE.txt")%>%as.data.frame()
sur<-sur[,-1]

sample=sur[1,]
OS=sur[18,]
time=sur[19,]

sur.df<-rbind(sample,time,OS)
colnames(sur.df)<-sur.df[1,]
sur.df2<-t(sur.df)
colnames(sur.df2)<-c("ID","OS.time","OS")
sur.df2<-as.data.frame(sur.df2)

library(stringr)

sur.df2$OS<-sur.df2$OS %>%str_remove("vital_status: ")
sur.df2$OS.time<-sur.df2$OS.time %>%str_remove("survival_time_in_days: ")

sur.df2$OS<-ifelse(sur.df2$OS=="Alive",0,ifelse(sur.df2$OS=="Dead",1,NA))

sur.df3<-sur.df2[!is.na(sur.df2$OS),]
sur.df3<-sur.df3[!is.na(sur.df3$OS.time),]
#去掉生存数据缺失值

sur.df3$OS.time<-as.numeric(sur.df3$OS.time)
sur.df3<-sur.df3[!is.na(sur.df3$OS.time),]
#n=398

######处理表达数据

exp<-fread("GSE72094_series_matrix.txt")%>%as.data.frame()
#探针表达谱

#####GPL15048
library(Biobase)
gpl<-fread('GPL15048.txt',header=T)%>%as.data.frame()
gpl2<-gpl[,c(1,4)]

gpl3<-subset(gpl2,gpl2$GeneSymbol%in%marker)

exp2<-subset(exp,exp$ID_REF%in%gpl3$ID)
colnames(exp2)[1]<-"ID"
exp4<-inner_join(gpl3,exp2)

###多个探针均值作为基因表达值
exp5<-aggregate(exp4,by=list(exp4$GeneSymbol),FUN=mean,na.rm=T) 
exp5<-exp5[,-2]
exp5<-exp5[,-2]
rownames(exp5)<-exp5$Group.1
exp5<-exp5[,-1]

exp6<-t(exp5)%>%as.data.frame()
exp6<-exp6[sur.df3$ID,]
identical(sur.df3$ID,rownames(exp6))

sur_exp_72094<-cbind(sur.df3,exp6)
####n=398

write.csv(sur_exp_72094,"D:/workplace/Immunotherapy/mime_survival/sur_exp_72094.csv",row.names = F)


# GSE68571 ----------------------------------------------------------------

setwd("D:/workplace/geneset/geo/GSE68571/")

######生存数据

sur<-fread("AGE.txt")%>%as.data.frame()
sur<-sur[,-1]

sample=colnames(sur)
OS=sur[19,]
time=sur[18,]

sur.df<-rbind(sample,time,OS)
colnames(sur.df)<-sur.df[1,]
sur.df2<-t(sur.df)
colnames(sur.df2)<-c("ID","OS.time","OS")
sur.df2<-as.data.frame(sur.df2)

library(stringr)

sur.df2$OS<-str_sub(sur.df2$OS,25,25)%>%as.numeric()
sur.df2$OS.time<-sur.df2$OS.time %>%str_remove("[followup time (months):]*")%>%as.numeric()

sur.df3<-sur.df2[!is.na(sur.df2$OS),]
sur.df3<-sur.df3[!is.na(sur.df3$OS.time),]
#去掉生存数据缺失值


######表达谱

exp<-fread("GSE68571_series_matrix.txt")%>%as.data.frame()
#探针表达谱

gpl<-fread("GPL80.annot")%>%as.data.frame()
gpl2<-gpl[,c(1,3)]

gpl3<-subset(gpl2,gpl2$`Gene symbol`%in%marker)

exp2<-subset(exp,exp$ID_REF%in%gpl3$ID)
colnames(exp2)[1]<-"ID"
exp4<-inner_join(gpl3,exp2)

###多个探针均值作为基因表达值
exp5<-aggregate(exp4,by=list(exp4$`Gene symbol`),FUN=mean,na.rm=T) 
exp5<-exp5[,-2]
exp5<-exp5[,-2]
rownames(exp5)<-exp5$Group.1
exp5<-exp5[,-1]

exp6<-t(exp5)%>%as.data.frame()
exp6<-exp6[rownames(sur.df3),]
identical(rownames(sur.df3),rownames(exp6))

sur_exp_68571<-cbind(sur.df3,exp6)
####n=86

####判断是否log2处理
ex <- sur_exp_68571[,-(1:3)]
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
exprSet <- log2(ex+1)
print("log2 transform finished")}else{print("log2 transform not needed")}

sur_exp_68571<-cbind(sur_exp_68571[,(1:3)],exprSet)

sur_exp_68571[is.na(sur_exp_68571)]=0

write.csv(sur_exp_68571,"D:/workplace/Immunotherapy/mime_survival/sur_exp_68571.csv",row.names = F)


# GSE68465 ----------------------------------------------------------------

setwd("D:/workplace/geneset/geo/GSE68465/")

######生存数据

sur<-fread("AGE.txt")%>%as.data.frame()
sur<-sur[,-1]

sample=sur[1,]
OS=sur[13,]
time=sur[20,]

sur.df<-rbind(sample,time,OS)
colnames(sur.df)<-sur.df[1,]
sur.df2<-t(sur.df)
colnames(sur.df2)<-c("ID","OS.time","OS")
sur.df2<-as.data.frame(sur.df2)

library(stringr)

sur.df2$OS<-sur.df2$OS %>%str_remove("vital_status: ")
sur.df2$OS.time<-sur.df2$OS.time %>%str_remove("months_to_last_contact_or_death: ")%>%as.numeric()

sur.df2$OS<-ifelse(sur.df2$OS=="Alive",0,ifelse(sur.df2$OS=="Dead",1,NA))

sur.df3<-sur.df2[!is.na(sur.df2$OS),]
sur.df3<-sur.df3[!is.na(sur.df3$OS.time),]
#去掉生存数据缺失值

######表达谱

exp<-fread("GSE68465_series_matrix.txt")%>%as.data.frame()
#探针表达谱

gpl<-fread("GPL96.annot")%>%as.data.frame()
gpl2<-gpl[,c(1,3)]

gpl3<-subset(gpl2,gpl2$`Gene symbol` %in%marker)

exp2<-subset(exp,exp$ID_REF%in%gpl3$ID)
colnames(exp2)[1]<-"ID"
exp4<-inner_join(gpl3,exp2)

###多个探针均值作为基因表达值
exp5<-aggregate(exp4,by=list(exp4$`Gene symbol`),FUN=mean,na.rm=T) 
exp5<-exp5[,-2]
exp5<-exp5[,-2]
rownames(exp5)<-exp5$Group.1
exp5<-exp5[,-1]

exp6<-t(exp5)%>%as.data.frame()
exp6<-exp6[sur.df3$ID,]
identical(sur.df3$ID,rownames(exp6))

sur_exp_68465<-cbind(sur.df3,exp6)
####n=86

####判断是否log2处理
ex <- sur_exp_68465[,-(1:3)]
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
exprSet <- log2(ex+1)
print("log2 transform finished")}else{print("log2 transform not needed")}

sur_exp_68465<-cbind(sur_exp_68465[,(1:3)],exprSet)


write.csv(sur_exp_68465,"D:/workplace/Immunotherapy/mime_survival/sur_exp_68465.csv",row.names = F)


# GSE11969 ----------------------------------------------------------------

setwd("D:/workplace/geneset/geo/GSE11969/")

library(data.table)
library(dplyr)


# 生存数据预处理

sur<-fread("AGE.txt")%>%as.data.frame()
sur<-sur[,-1]

sample=sur[1,]
status=sur[17,]
time=sur[18,]

sur.df<-rbind(sample,time,status)
sur.df<-sur.df[,1:90]
colnames(sur.df)<-sur.df[1,]

sur.df2<-t(sur.df)
colnames(sur.df2)<-c("ID","time","status")
sur.df2<-as.data.frame(sur.df2)

library(stringr)
sur.df2$status<-sur.df2$status %>%str_remove("Status: ")
sur.df2$time<-sur.df2$time%>%str_remove("Survival \\(days\\)\\: ")%>%as.numeric()

sur.df2$status<-ifelse(sur.df2$status=="Alive",0,ifelse(sur.df2$status=="Dead",1,NA))%>%as.numeric()

sur.df3<-sur.df2[!is.na(sur.df2$status),]
sur.df3<-sur.df3[!is.na(sur.df3$time),]
#去掉生存数据缺失值

####处理表达

exp<-fread("GSE11969_series_matrix.txt")%>%as.data.frame()
#探针表达谱

#####GPL

gpl<-fread('GPL7015.txt',header=T)%>%as.data.frame()
gpl2<-gpl[,c(1,6)]

gpl3<-subset(gpl2,gpl2$`Gene symbol`%in%marker)

exp2<-subset(exp,exp$ID_REF%in%gpl3$ID)
colnames(exp2)[1]<-"ID"
exp4<-inner_join(gpl3,exp2)

###多个探针均值作为基因表达值
exp5<-aggregate(exp4,by=list(exp4$`Gene symbol`),FUN=mean,na.rm=T) 
exp5<-exp5[,-2]
exp5<-exp5[,-2]
rownames(exp5)<-exp5$Group.1
exp5<-exp5[,-1]

exp6<-t(exp5)%>%as.data.frame()
exp6<-exp6[sur.df3$ID,]
identical(sur.df3$ID,rownames(exp6))

sur_exp_11969<-cbind(sur.df3,exp6)
colnames(sur_exp_11969)[1:3]<-c("ID","OS.time","OS")
####n=86

e<-sur_exp_11969[,-(1:3)]
e2<-e+10
sur_exp_11969<-cbind(sur_exp_11969[,(1:3)],e2)

write.csv(sur_exp_11969,"D:/workplace/Immunotherapy/mime_survival/sur_exp_11969.csv",row.names = F)


# GSE87340 ----------------------------------------------------------------

setwd("D:/workplace/geneset/geo/GSE87340/")

######生存数据

sur<-fread("AGE.txt")%>%as.data.frame()
sur<-sur[,-1]

sample=sur[1,]
status=sur[17,]
time=sur[18,]

sur.df<-rbind(sample,time,status)
sur.df2<-t(sur.df)
colnames(sur.df2)<-c("sample","time","status")
sur.df2<-as.data.frame(sur.df2)

library(stringr)
sur.df2$status<-sur.df2$status %>%str_remove("survstatus: ")%>%as.numeric()
sur.df2$time<-sur.df2$time %>%str_remove("time to last followup: ")%>%as.numeric()
###

sur.df3<-sur.df2[!is.na(sur.df2$status),]
sur.df3<-sur.df3[!is.na(sur.df3$time),]
#去掉生存数据缺失值

######表达谱

exp<-fread("GSE87340_RPKM_log2.txt")%>%as.data.frame()
#探针表达谱

exp2<-subset(exp,exp$ID%in%marker)

rownames(sur.df3)<-str_remove(rownames(sur.df3),pattern = ".RNA-seq")

rownames(exp2)<-exp2$ID
exp2<-exp2[,-1]
exp6<-t(exp2)%>%as.data.frame()
exp6<-exp6[rownames(sur.df3),]
identical(rownames(sur.df3),rownames(exp6))

sur_exp_87340<-cbind(sur.df3,exp6)
colnames(sur_exp_87340)[1:3]<-c("ID","OS.time","OS")
####n=86

write.csv(sur_exp_87340,"D:/workplace/Immunotherapy/mime_survival/sur_exp_87340.csv",row.names = F)


# GSE31852 ----------------------------------------------------------------
setwd("D:/workplace/geneset/geo/GSE31852/")

######生存数据

sur<-fread("AGE.txt")%>%as.data.frame()
sur<-sur[,-1]

sample<-c()
status<-c()
time<-c()

s<-c("pfsc (1=progressed; 0=not progressed): 1","pfsc (1=progressed; 0=not progressed): 0","progression-free survival status: 1","progression-free survival status: 0")

for(i in 1:dim(sur)[2]){
  
  x=sur[,i]
  sample[i]<-x[1]
  
  for(j in 1:length(x)){
    if(x[j]%in%s){
      status[i]<-x[j]
      time[i]<-x[j+1]
    }
  }
}


sur.df<-rbind(sample,time,status)
colnames(sur.df)<-sur.df[1,]
sur.df2<-t(sur.df)
colnames(sur.df2)<-c("sample","time","status")
sur.df2<-as.data.frame(sur.df2)

sur.df2<-subset(sur.df2,!is.na(sur.df2$status))

status2<-str_split(sur.df2$status,pattern = ":")%>%unlist()
sur.df2$status<-status2[seq(0,length(status2),2)]

time2<-str_split(sur.df2$time,pattern = ":")%>%unlist()
sur.df2$time<-time2[seq(0,length(time2),2)]%>%as.numeric()

sur.df3<-sur.df2[!is.na(sur.df2$status),]
sur.df3<-sur.df3[!is.na(sur.df3$time),]
#去掉生存数据缺失值


######表达谱

#####
exp<-fread("GSE31852_series_matrix.txt")%>%as.data.frame()
#探针表达谱

gpl<-fread("GPL6244.annot")%>%as.data.frame()
gpl2<-gpl[,c(1,3)]

gpl3<-subset(gpl2,gpl2$`Gene symbol`%in%marker)

exp2<-subset(exp,exp$ID_REF%in%gpl3$ID)
colnames(exp2)[1]<-"ID"
exp4<-inner_join(gpl3,exp2)

###多个探针均值作为基因表达值
exp5<-aggregate(exp4,by=list(exp4$`Gene symbol`),FUN=mean,na.rm=T) 
exp5<-exp5[,-2]
exp5<-exp5[,-2]
rownames(exp5)<-exp5$Group.1
exp5<-exp5[,-1]

exp6<-t(exp5)%>%as.data.frame()
exp6<-exp6[rownames(sur.df3),]
identical(rownames(sur.df3),rownames(exp6))

sur_exp_31852<-cbind(sur.df3,exp6)
colnames(sur_exp_31852)[1:3]<-c("ID","OS.time","OS")
####n=86

write.csv(sur_exp_31852,"D:/workplace/Immunotherapy/mime_survival/sur_exp_31852.csv",row.names = F)


# TCGA --------------------------------------------------------------------

setwd("D:/workplace/Immunotherapy/TCGA/")

library(data.table)
library(dplyr)

####表达谱
tcga_exp<-fread("HiSeqV2")%>%as.data.frame()
#log2(norm_count+1)

rownames(tcga_exp)<-tcga_exp$sample
tcga_exp<-tcga_exp[,-1]

exp_df<-t(tcga_exp)%>%as.data.frame()
exp_df<-exp_df[,intersect(colnames(exp_df),marker)]
####生存

sur<-fread("lung-survival.txt")%>%as.data.frame()

sur_df<-data.frame(ID=sur$sample,
                   OS.time=sur$OS.time,
                   OS=sur$OS)

sample<-intersect(sur_df$ID,rownames(exp_df))
#n=1077 表达和生存对应的样本

rownames(sur_df)<-sur_df$ID

survival<-sur_df[sample,]
expre<-exp_df[sample,]

identical(rownames(survival),rownames(expre))

sur_exp_TCGA<-cbind(survival,expre)#n=1077

write.csv(sur_exp_TCGA,"D:/workplace/Immunotherapy/mime_survival/sur_exp_TCGA.csv",row.names = F)

# # 将8个数据框合并到一个列表中 --------------------------------------------------------
setwd("D:/workplace/Immunotherapy/mime_survival/")

list_train_vali_Data <- list(Dataset1 = sur_exp_TCGA, 
                             Dataset2 = sur_exp_26939,
                             Dataset3 = sur_exp_72094,
                             #Dataset4 = sur_exp_68571,
                             Dataset5 = sur_exp_68465,
                             Dataset6 = sur_exp_11969,
                             Dataset7 = sur_exp_87340,
                             Dataset8 = sur_exp_31852
                             )

save(list_train_vali_Data,file = 'datasets8.Rdata')

load("datasets8.Rdata")

# 模型选择 --------------------------------------------------------------------
rm(list = ls())
options(stringsAsFactors = FALSE)
load("D:/workplace/Immunotherapy/mime_survival/datasets8.Rdata")
marker<-read.table("D:/workplace/Immunotherapy/TCGA/EMT_M1_M2_marker.txt")
marker<-marker$V1

library(fastAdaboost)
library(cancerclass)
library(Mime1)
library(mixOmics)
library(survcomp)
library(CoxBoost)


res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data2$Dataset1,
                       list_train_vali_Data = list_train_vali_Data2,
                       unicox.filter.for.candi = T,
                       unicox_p_cutoff = 0.05,
                       candidate_genes = marker,
                       mode = 'all',nodesize =5,seed = 123 )
