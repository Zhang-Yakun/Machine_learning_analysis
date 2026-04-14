
# TCGA-Lung cancer --------------------------------------------------------

setwd("D:/workplace/Immunotherapy/TCGA/")

library(data.table)
library(dplyr)

tcga_exp<-fread("HiSeqV2")%>%as.data.frame()
#log2(norm_count+1)

colnames(tcga_exp)[1]<-"GENE"
rownames(tcga_exp)<-tcga_exp$GENE
tcga_exp<-tcga_exp[,-1]

a<-2^(tcga_exp)-1
#counts matrix

###COUNTS 转 TPM

#expM<-ifelse(expMatrix<0,0,expMatrix) 
expMatrix<-a

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
## 17557 * 1129

tpm_df<-as.data.frame(expMatrix_tpm)
tpm_df$GENE<-rownames(tpm_df)

tpm_df2<-tpm_df%>%select("GENE",everything())

write.table(tpm_df2,"DATA.txt",
            row.names = F, quote = FALSE,sep = "\t")
## tpm (17557 * 1129)

# cibersort ---------------------------------------------------------------
remove(list = ls()) #一键清空
#加载包
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)

setwd("D:/workplace/Immunotherapy/TCGA/Cibersort/")

source('Cibersort.R')

result1 <- CIBERSORT('LM22.txt','DATA.txt', perm = 1000, QN = T)
#perm置换次数=1000，QN分位数归一化=TRUE


# M1,M2细胞比例 ---------------------------------------------------------------

cell_exp<-fread("CIBERSORT-Results_TCGA.txt")%>%as.data.frame()

M1_M2<-data.frame(Macrophages_M1=cell_exp$`Macrophages M1`,
                  Macrophages_M2=cell_exp$`Macrophages M2`)
rownames(M1_M2)<-cell_exp$Mixture

###去掉M1和M2含量为0的样本

mac_df<-subset(M1_M2,M1_M2$Macrophages_M1!=0)
mac_df<-subset(mac_df,mac_df$Macrophages_M2!=0)
##1083非零样本

# TCGA-生存分析 ---------------------------------------------------------------

setwd("D:/workplace/Immunotherapy/TCGA/")

sur<-fread("D:/workplace/Immunotherapy/TCGA/lung-survival.txt")%>%as.data.frame()

sur_df<-data.frame(sample=sur$sample,
                   os=sur$OS,
                   os.time=sur$OS.time)

sample<-intersect(sur_df$sample,rownames(mac_df))
#1032 表达和生存对应的样本

rownames(sur_df)<-sur_df$sample

survival<-sur_df[sample,]
mac_df<-mac_df[sample,]

mac_df$ratio<-(mac_df$Macrophages_M1/mac_df$Macrophages_M2)%>%as.numeric()

mac_df_group<-ifelse(mac_df$ratio>=1,"M1>=M2","M1<M2")
table(mac_df_group)

survival$group<-mac_df_group

####

t.test(mac_df$Macrophages_M1,mac_df$Macrophages_M2)

####生存分析
library(survival)
library(survminer)
library(ggpubr)

fit1 <- survfit(Surv(os.time, os) ~ group, data = survival)

ggsurvplot(fit1,
           pval = TRUE, #添加logrank检验的p值
           conf.int = TRUE, #添加置信区间
           risk.table = TRUE, #添加风险表
           risk.table.col='strata', #根据数据分组为风险表添加颜色
           linetype = 'strata', #不同组别生存曲线线型
           surv.median.line = 'hv', #标注中位生存时间
           ggtheme = theme_bw() #主题风格
           #palette =  c("#00BFC4","#F8766D")
           ) #颜色风格

#####

fit2 <- survfit(Surv(os.time, os) ~ group2, data = survival)

ggsurvplot(fit2,
           pval = TRUE, #添加logrank检验的p值
           conf.int = TRUE, #添加置信区间
           risk.table = TRUE, #添加风险表
           risk.table.col='strata', #根据数据分组为风险表添加颜色
           linetype = 'strata', #不同组别生存曲线线型
           surv.median.line = 'hv', #标注中位生存时间
           ggtheme = theme_bw(), #主题风格
           palette =  c("#00BFC4","#F8766D")
) #颜色风格

#####

library(survival)
library(survminer)
library(ggpubr)

survival2<-subset(survival,survival$group%in%c("M1low & M2high","M1high & M2low"))

table(survival2$group)

fit3 <- survfit(Surv(os.time, os) ~ group, data = survival2)

ggsurvplot(fit3,
           pval = TRUE, #添加logrank检验的p值
           conf.int = TRUE, #添加置信区间
           risk.table = TRUE, #添加风险表
           risk.table.col='strata', #根据数据分组为风险表添加颜色
           linetype = 'strata', #不同组别生存曲线线型
           surv.median.line = 'hv', #标注中位生存时间
           ggtheme = theme_bw() #主题风格
           #palette =  c("#00BFC4","#F8766D")
) #颜色风格


# cox回归筛选风险基因（133genes） ---------------------------------------------------

tpm<-fread("DATA.txt")%>%as.data.frame()
rownames(tpm)<-tpm$GENE
tpm<-tpm[,-1]

log_tpm<-log2(tpm+1)

cox_input<-t(log_tpm)%>%as.data.frame()

cox_input$sample<-rownames(cox_input)
cox_input<-cox_input %>% select("sample",everything())

cox_input<-inner_join(sur_df,cox_input)
rownames(cox_input)<-cox_input$sample

# 多个协变量-cox ---------------------------------------------------------------
gene<-read.table("EMT_M1_M2_marker.txt")
#133 signatures

covariates <- intersect(gene$V1,colnames(cox_input))
#130 signatures

cox_input2<-cbind(cox_input[,2:3],cox_input[,covariates])
#1077sample,130genes


univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(os.time, os)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = cox_input2)})
# 提取数据，并制作数据表格 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
res_df<-as.data.frame(res)

sig<-subset(res_df,res_df$p.value<0.05)

write.csv(sig,"cox_sig9.csv")

#####9个signature的风险比可视化

library(forestplot)
library(grid)
library(checkmate)

# 提取风险比、95%CI、P值等，构建森林图的基本表

# 提取风险比、95%CI、P值等，构建森林图的基本表

# 第一种方法，在之前的cox回归介绍中，我们就是这么提取的！

pvalue <- sig$p.value
HR <- sig$`HR (95% CI for HR)`%>%str_extract("[0-9.]*")%>%as.numeric()
low <- sig$`HR (95% CI for HR)`%>%
  str_extract("\\([0-9.]*")%>%
  str_remove("\\(")%>%
  as.numeric()
  
high <- sig$`HR (95% CI for HR)`%>%
  str_extract("-[0-9.]*")%>%
  str_remove("-")%>%
  as.numeric()

# 整合为数据框，部分调整是为了可视化效果，大家可以自行设定
res_mul_cox <- data.frame(gene=rownames(sig),p.value = pvalue, HR = HR, low = low, high = high, stringsAsFactors = F)

# 给它加一列 HR (95%CI)
res_mul_cox$`HR (95%CI)` = paste0(HR, " (", low, "-", high, ")", sep = "")

################################ forestplot 包 #################################

# 简单的基本森林图
fig <- forestplot(
  res_mul_cox[, c(1,2,6)],   # 选择要在森林图中显示的数据列，第1、5、6列
  mean = res_mul_cox[,3],     # 指定均值数据列（HR），它将显示为森林图的小方块
  lower = res_mul_cox[,4],    # 指定95%置信区间的下限数据列
  upper = res_mul_cox[,5],    # 指定95%置信区间的上限数据列，这些数据将显示为线段穿过方块
  zero = 1,               # 设置零线或参考线为HR=1，这是x轴的垂直线
  boxsize = 0.1,          # 设置小方块的大小
  graph.pos = 2           # 指定森林图应该插入到图形中的位置，这里是第2列
)
fig

####优化森林图

# 更细节地优化一波

# 创建森林图 (forest plot) 并存储在 fig3_1 变量中
fig2 <- forestplot(
  res_mul_cox[, c(1,2,6)],   # 选择要在森林图中显示的数据列，第1、5、6列
  mean = res_mul_cox[,3],     # 指定均值数据列（HR），它将显示为森林图的小方块
  lower = res_mul_cox[,4],    # 指定95%置信区间的下限数据列
  upper = res_mul_cox[,5],    # 指定95%置信区间的上限数据列，这些数据将显示为线段穿过方块
  zero = 1,                        # 均值为1时的参考线，也就是零线
  boxsize = 0.3,                   # 方框的大小
  graph.pos = "right",             # 森林图在右侧显示
  hrzl_lines = list(               # 水平线样式的设置
    "1" = gpar(lty = 1, lwd = 2),  # 均值线
    "2" = gpar(lty = 2),           # 下限和上限之间的虚线
    "10" = gpar(lwd = 2, lty = 1, columns = c(1:4)) # 下限和上限线
  ),
  graphwidth = unit(.25, "npc"),   # 森林图的宽度
  xlab = "\n ", # x轴标签
  fn.ci_norm = "fpDrawDiamondCI",  # 置信区间的形状（这里是钻石形）
  title = "Cox回归森林图",   # 森林图的标题
  col = fpColors(                  # 颜色设置
    box = "blue4",                 # 方框颜色
    lines = "blue4",               # 线条颜色
    zero = "black"                 # 均值为0时的颜色
  )
)
fig2


