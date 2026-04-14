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

# TCGA-生存分析 ---------------------------------------------------------------

setwd("D:/workplace/Immunotherapy/TCGA/")

sur<-fread("lung-survival.txt")%>%as.data.frame()

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

mac_df2<-data.frame(value=c(mac_df$Macrophages_M1,mac_df$Macrophages_M2),
                    cell=rep(c("M1","M2"),c(dim(mac_df)[1],dim(mac_df)[1])))

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
