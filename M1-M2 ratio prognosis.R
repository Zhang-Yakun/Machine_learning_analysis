library(dplyr)
library(data.table)

setwd("D:/workplace/Immunotherapy/TCGA/Cibersort/")

# M1,M2细胞比例 ---------------------------------------------------------------

cell_exp<-fread("CIBERSORT-Results_TCGA.txt")%>%as.data.frame()

M1_M2<-data.frame(Macrophages_M1=cell_exp$`Macrophages M1`,
                  Macrophages_M2=cell_exp$`Macrophages M2`)
rownames(M1_M2)<-cell_exp$Mixture

###去掉M1和M2含量为0的样本

mac_df<-subset(M1_M2,M1_M2$Macrophages_M1!=0)
mac_df<-subset(mac_df,mac_df$Macrophages_M2!=0)
##1083非零样本

###M1/M2排序
mac_df$ratio<-(mac_df$Macrophages_M1/mac_df$Macrophages_M2)
quantile(mac_df$ratio)

mac_df<-mac_df[order(mac_df$ratio,decreasing = T),]


# M1/M2降序细胞比例 -----------------------------------------------------------

t.test(mac_df$Macrophages_M1,mac_df$Macrophages_M2)

mac_df2<-data.frame(value=c(mac_df$Macrophages_M1,mac_df$Macrophages_M2),
                    cell=rep(c("M1","M2"),c(dim(mac_df)[1],dim(mac_df)[1])))

library(ggplot2)
# 绘制上半部分密度图
p_top <- ggplot(mac_df2, aes(x = value, color = cell, fill =cell)) +
  geom_density() +
  # 让箱子的所有位置都颜色统一，如例文所示
  scale_color_manual(values = c(alpha("#D7C3DF",0.7),alpha("#AE84BE",0.7))) + # 设置透明色
  scale_fill_manual(values = c(alpha("#D7C3DF",0.7),alpha("#AE84BE",0.7))) +
  theme_classic() + # 如果显示采用这一行#4b84b3，#f89f68
  
  xlab("TCGA percentage of M1 & M2") +
  
  ylab(NULL) + 
  theme(legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(size = 12,color = "black"),
        axis.text.y = element_blank(), # 原文不显示纵轴的密度
        #axis.text.y = element_text(size = 12,color = "black"), # 如果要显示采用这一行
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_rug()

#绘制下半部分图
p_bot <-ggplot(mac_df2, aes(value, fill = cell))  + 
  geom_boxplot(aes(col =cell)) + 
  scale_fill_manual(values = c(alpha("#D7C3DF",1),alpha("#AE84BE",1))) + 
  scale_color_manual(values = c(alpha("#D7C3DF",1),alpha("#AE84BE",1))) + 
  xlab(NULL) + 
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        axis.text.x = element_blank(), # 原文不显示箱线图的x轴
        #axis.text.x = element_text(size = 12,color = "black"), # 如要显示箱线图x轴采用这一行
        axis.text.y = element_text(size = 11,color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

#组合图
library(patchwork)
p_top + p_bot +
  plot_layout(nrow = 2)

# 对应样本 --------------------------------------------------------------------

tcga_exp<-fread("D:/workplace/Immunotherapy/TCGA/HiSeqV2")%>%as.data.frame()
#log2(norm_count+1)

rownames(tcga_exp)<-tcga_exp$sample
tcga_exp<-tcga_exp[,-1]

#####生存数据
sur<-fread("D:/workplace/Immunotherapy/TCGA/lung-survival.txt")%>%as.data.frame()

sur_df<-data.frame(sample=sur$sample,
                   os=sur$OS,
                   os.time=sur$OS.time)

sample<-intersect(sur_df$sample,colnames(tcga_exp))
#1077 表达和生存对应的样本

sample2<-intersect(sample,rownames(mac_df))
#1032 表达生存免疫对应样本

rownames(sur_df)<-sur_df$sample

surv_df<-sur_df[sample2,]
expr_df<-tcga_exp[,sample2]
mac_df<-mac_df[sample2,]
#1032总样本数目

# 筛选两组样本 -----------------------------------------------------------------

high_ratio<-subset(mac_df,mac_df$ratio>=quantile(mac_df$ratio,0.75))
#258
low_ratio<-subset(mac_df,mac_df$ratio<=quantile(mac_df$ratio,0.25))
#258
ratio_group<-data.frame(sample=c(rownames(high_ratio),rownames(low_ratio)),
                        group=rep(c("high","low"),c(258,258)))
#516

expr_df<-expr_df[,ratio_group$sample]
surv_df<-surv_df[ratio_group$sample,]
surv_df$group<-ratio_group$group

###输出TCGA数据集
sig4<-read.table("D:/workplace/Immunotherapy/TCGA/M1M2/marker108.txt")
exp_df3<-expr_df[intersect(rownames(expr_df),sig4$V1),]
exp_df4<-t(exp_df3)%>%as.data.frame()

sur_df4<-data.frame(ID=surv_df$sample,
                   OS.time=surv_df$os.time,
                   OS=surv_df$os)

sample4<-intersect(sur_df4$ID,rownames(exp_df4))
#n=1077 表达和生存对应的样本

rownames(sur_df4)<-sur_df4$ID

surv4<-sur_df4[sample4,]
expr4<-exp_df4[sample4,]

identical(rownames(surv4),rownames(expr4))

sur_exp_TCGA<-cbind(surv4,expr4)#n=1077

write.csv(sur_exp_TCGA,"D:/workplace/Immunotherapy/TCGA/M1M2/sur_exp_TCGA.csv",row.names = F)


# 读取表达数据-差异分析 -------------------------------------------------------------

#####limma
group<-rep(c("high","low"),c(258,258))
group<-factor(group,levels = c("high","low"))


library(limma)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(expr_df , design)
contrast.matrix <- makeContrasts(high - low,
                                 levels = design
)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

tempOutput <- topTable(fit2, coef = 1, n = Inf)
nrDEG1 <- na.omit(tempOutput)  #P值和logFC
nrDEG1$gene <-rownames(nrDEG1)

dif1<-subset(nrDEG1,nrDEG1$adj.P.Val<0.05)
dif2<-subset(dif1,abs(dif1$logFC)>1)

gene2<-unique(dif2$gene)
library(stringr)
gene3<-gene2[-which(gene2=="?|729884")]
gene3<-gene3[-which(gene3=="NKX2-1")]
gene<-gene3[-which(gene3=="NKX2-5")]
##1363个差异表达基因

write.csv(dif2,"M1_M2_limma_diff.csv",row.names = F)


# cox回归分析 -----------------------------------------------------------------
library(survival)
library(survminer)
library(ggpubr)

expr_df2<-expr_df[gene,]
cox_input<-t(expr_df2)%>%as.data.frame()

cox_input$sample<-rownames(cox_input)
cox_input<-cox_input %>% select("sample",everything())

cox_input<-inner_join(surv_df,cox_input)
rownames(cox_input)<-cox_input$sample

covariates <- intersect(gene,colnames(cox_input))
#1363 signatures

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
#108

write.csv(sig,"cox_sig_M1M2_res.csv")
write.table(rownames(sig),"marker108.txt",row.names = F,col.names = F,quote = F)
