
setwd("D:/workplace/Immunotherapy/datasets/M1M2/")

#install.packages("glmnet")
library(glmnet)


# 提取目标基因 ------------------------------------------------------------------

autogene<-read.table("model_sig57.txt",header = F,stringsAsFactors = F,sep = "\t")

# 表达谱 ---------------------------------------------------------------------
load("datasets3.Rdata")
allgene2<-list_train_vali_Data[[1]]

allgene<-allgene2[,-(1:2)]%>%t()%>%as.data.frame()
allgene$symbol<-rownames(allgene)

autogene<-data.frame(autogene$V1)
colnames(autogene)<-c("symbol")

autoexp<-merge(autogene,allgene,by="symbol")
write.table(autoexp,"dif_expr.txt",quote = F,sep="\t",row.names = F)

rownames(autoexp)<-autoexp[,1]
autoexp<-autoexp[,-1]
autoexp<-as.data.frame(t(autoexp))
autoexp$sample<-row.names(autoexp)
clinical<-read.table("clinical0.txt",quote="",sep="\t",header = T)[,1:2]


# 细胞比例矩阵 ------------------------------------------------------------------
library(data.table)
cell<-fread("CIBERSORT-Results3.txt")%>%as.data.frame()

cell$`Macrophages M2`[16]<-0.0000000001

cell_df<-data.frame(sample=cell$Mixture,
                    M1=cell$`Macrophages M1`,
                    M2=cell$`Macrophages M2`,
                    M1M2=(cell$`Macrophages M1`/cell$`Macrophages M2`))
rownames(cell_df)<-cell$Mixture


# 分组信息 --------------------------------------------------------------------
library(stringr)
library(dplyr)

group_list=allgene2[,2]
table(group_list)

identical(autoexp$sample,rownames(allgene2))

group_list = factor(group_list,
                    levels = c("Y","N"))

auto_cli<-data.frame(group=group_list)
rownames(auto_cli)<-rownames(allgene2)
auto_cli<-cbind(auto_cli,autoexp)

auto_cli2<-inner_join(auto_cli,cell_df)
rownames(auto_cli2)<-auto_cli2$sample
auto_cli<-auto_cli2[,-59]
colnames(auto_cli)

# 单变量Logistic -------------------------------------------------------------


uni_glm_model<-function(x){#限合结局和变量
  FL<-as.formula(paste("group ~ " ,x,sep=""))
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
  return(uni_glm_model)
}

variable.names<- colnames(auto_cli)[c(2:dim(auto_cli)[2])];variable.names
uni_logistic<- lapply(variable.names,uni_glm_model)#批量输出结果并合并在一起
uni_logistic
library(plyr)
uni_glm<-ldply(uni_logistic,data.frame);uni_glm
write.table(uni_glm,"uni_logistic.txt",sep="\t",quote=F,row.names=F)

uni_glm<-subset(uni_glm,uni_glm$P<0.05)

# 绘制森林图 -------------------------------------------------------------------


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
reg_lasso_res <- reg_lasso(inputArr=auto_cli, x_cols=names(auto_cli)[-1], y_cols="group", reg_family="binomial", return_fitted_values=FALSE, pred_type="link")
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


#模型预测
#install.packages("pROC")
library(pROC)
clinical[1:39,4]<-0
clinical[40:63,4]<-1#转换变量
pre_logistic<-as.numeric(predict(glm1,type = "response")>0.5)
rocl<-roc(clinical$V4,pre_logistic)
pdf("ROC.pdf",width = 9,height = 9)
plot(rocl, # roc1换为roc2，更改参数可绘制roc2曲
     print.auc=TRUE,print.auc.x=0.5,print.auc.y=0.5, # 图像上输出AUC值,坐标为（x，y）
     auc.polygon=TRUE, auc.polygon.col="skyblue", # 设置ROC曲线下填充色
     max.auc.polygon=TRUE, # 填充整个图像
     grid=c(0.1,0.1), grid.col=c("green", "red"), # 设置间距为0.1，0.2，线条颜色
     #print.thres=TRUE, print.thres.cex=0.8,  # 图像上输出最佳截断值，字体缩放0.8倍
     legacy.axes=TRUE)  # 使横轴从0到1，表示为1-特异度
dev.off()
