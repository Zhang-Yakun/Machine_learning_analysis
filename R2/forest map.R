library(dplyr)
library(data.table)

rm(list = ls())

####20个共有的预后特征
covariates <- read.table("D:/workplace/Immunotherapy/TCGA/M1M2/dataset6_result/cox_dataset6_signature36.txt")
covariates<-covariates$V1

# TCGA-cox ----------------------------------------------------------------

##cox输入数据
tcga<-read.csv("D:/workplace/Immunotherapy/TCGA/M1M2/sur_exp_TCGA.csv")
rownames(tcga)<-tcga$ID
tcga<-tcga[,-1]
tcga2<-tcga%>%dplyr::select(OS,OS.time,everything())
colnames(tcga2)[1:2]<-c("os","os.time")

cox_input2<-tcga2
####cox

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

# 森林图 ---------------------------------------------------------------------

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

###按P值排序，从小到大
res_mul_cox<-res_mul_cox[order(res_mul_cox$p.value),]


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
    "37" = gpar(lwd = 2, lty = 1, columns = c(1:4)) # 下限和上限线
  ),
  graphwidth = unit(.25, "npc"),   # 森林图的宽度
  xlab = "\n ", # x轴标签
  fn.ci_norm = "fpDrawDiamondCI",  # 置信区间的形状（这里是钻石形）
  title = "TCGA Cox regression",   # 森林图的标题
  col = fpColors(                  # 颜色设置
    box = "blue4",                 # 方框颜色
    lines = "blue4",               # 线条颜色
    zero = "black"                 # 均值为0时的颜色
  )
)
fig2


## GEO-cox ----------------------------------------------------------------

##cox输入数据
tcga<-read.csv("D:/workplace/Immunotherapy/TCGA/M1M2/sur_exp_87340.csv")
rownames(tcga)<-tcga$ID
tcga<-tcga[,-1]
tcga2<-tcga%>%dplyr::select(OS,OS.time,everything())
colnames(tcga2)[1:2]<-c("os","os.time")

cox_input2<-tcga2
####cox

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

# 森林图 ---------------------------------------------------------------------

library(forestplot)
library(grid)
library(checkmate)

# 提取风险比、95%CI、P值等，构建森林图的基本表

# 提取风险比、95%CI、P值等，构建森林图的基本表

# 第一种方法，在之前的cox回归介绍中，我们就是这么提取的！

pvalue <- res_df$p.value
HR <- res_df$`HR (95% CI for HR)`%>%str_extract("[0-9.]*")%>%as.numeric()
low <- res_df$`HR (95% CI for HR)`%>%
  str_extract("\\([0-9.]*")%>%
  str_remove("\\(")%>%
  as.numeric()

high <- res_df$`HR (95% CI for HR)`%>%
  str_extract("-[0-9.]*")%>%
  str_remove("-")%>%
  as.numeric()

# 整合为数据框，部分调整是为了可视化效果，大家可以自行设定
res_mul_cox <- data.frame(gene=rownames(res_df),p.value = pvalue, HR = HR, low = low, high = high, stringsAsFactors = F)

# 给它加一列 HR (95%CI)
res_mul_cox$`HR (95%CI)` = paste0(HR, " (", low, "-", high, ")", sep = "")
###按P值排序，从小到大
res_mul_cox<-res_mul_cox[order(res_mul_cox$p.value),]

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
    "37" = gpar(lwd = 2, lty = 1, columns = c(1:4)) # 下限和上限线
  ),
  graphwidth = unit(.25, "npc"),   # 森林图的宽度
  xlab = "\n ", # x轴标签
  fn.ci_norm = "fpDrawDiamondCI",  # 置信区间的形状（这里是钻石形）
  title = "GSE87340 Cox regression",   # 森林图的标题
  col = fpColors(                  # 颜色设置
    box = "blue4",                 # 方框颜色
    lines = "blue4",               # 线条颜色
    zero = "black"                 # 均值为0时的颜色
  )
)
fig2


