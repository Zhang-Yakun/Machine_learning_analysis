#TCGA-LUAD\M1M2C5C8\中位数\不筛选正常肿瘤差异-直接对肿瘤样本打分分组生存\免疫浸润
#1
setwd("C:\\Users\\123\\Desktop\\RESULT\\TCGA-LUAD\\M1M2C5C8\\中位数\\不筛选正常肿瘤差异-直接对肿瘤样本打分分组生存\\免疫浸润")



LUAD_Tumor<-read.csv("F:\\RESULT\\表达数据\\LUAD肿瘤样本基因表达谱(tpm).csv",row.names=1,check.names=F)



C5_marker<-unlist(read.table("D:\\肺癌极化基因聚类\\C5-marker.txt"))
C8_marker<-unlist(read.table("D:\\肺癌极化基因聚类\\C8-marker.txt"))
M1_Polarization<-unlist(read.table("D:\\肺癌极化基因聚类\\M1_Polarization.txt"))
M2_Polarization<-unlist(read.table("D:\\肺癌极化基因聚类\\M2_Polarization.txt"))
geneA<-intersect(M1_Polarization,C8_marker)
geneB<-intersect(M2_Polarization,C5_marker)
gene_sets <- list(
  M1_C8_Marker = geneA,
  M2_C5_Marker = geneB)


library(GSVA)
library(dplyr)
expr_matrix <- as.matrix(LUAD_Tumor)

gsvaPar <- gsvaParam(
  exprData = expr_matrix,
  geneSets = gene_sets,
  kcdf = "Gaussian"         # 连续型表达数据
)
gsva_scores <- gsva(gsvaPar, verbose = FALSE)
gsva_scores <- as.data.frame(t(gsva_scores))


m1_median <- median(gsva_scores$M1_C8_Marker)
m2_median <- median(gsva_scores$M2_C5_Marker)
gsva_scores$M1_Group <- ifelse(gsva_scores$M1_C8_Marker >= m1_median, "M1_High", "M1_Low")
gsva_scores$M2_Group <- ifelse(gsva_scores$M2_C5_Marker >= m2_median, "M2_High", "M2_Low")
gsva_scores$Phenotype_Median <- paste(gsva_scores$M1_Group, gsva_scores$M2_Group, sep = "_")  # 专属列名
table(gsva_scores$Phenotype_Median)







library(readr)
library(tidyverse)
library(dplyr)
library(data.table)
clin<- fread("D:\\TCGA数据下载\\TCGA-LUAD\\TCGA-LUAD-Tumor\\clinical.tsv",data.table = F)
metadata<-read.csv("D:\\肺癌极化基因聚类\\TCGA肺癌数据分析\\TCGA-LUAD\\LUAD-Tumor\\TCGA_LUAD_sample_info.csv",row.names=1,check.names = F)

yingshe<-data.frame(metadata$sample_id,metadata$TCGA_case_id)
colnames(yingshe)<-c("cases.submitter_id","case_id")

head(yingshe)
clin <- clin[, c(
  "cases.case_id",
  "demographic.vital_status",
  "demographic.days_to_death",
  "diagnoses.days_to_last_follow_up",
  "demographic.age_at_index",
  "demographic.gender",
  "diagnoses.ajcc_pathologic_t",
  "diagnoses.ajcc_pathologic_m",
  "diagnoses.ajcc_pathologic_n",
  "diagnoses.ajcc_pathologic_stage"
)]
head(clin)
dim(yingshe)
dim(clin)
library(dplyr)
library(tidyr)
yingshe_clean <- yingshe %>%
  distinct(cases.submitter_id, .keep_all = TRUE) # 按样本ID去重

# ----------------- 2. 预处理 clin 表 (临床数据) -----------------
# 2.1 强制转换数据类型
# 将随访时间转换为数值型。如果原数据是 "--" 或 "NA"，转换后会变成 NA
clin$diagnoses.days_to_last_follow_up <- as.numeric(clin$diagnoses.days_to_last_follow_up)
clin$demographic.days_to_death <- as.numeric(clin$demographic.days_to_death)

# 2.2 对 clin 表进行预去重
# 同样先按 case_id 去重，防止合并时出现多对多
# 这里先简单去重，保留第一条，具体的生存逻辑留到合并后做
clin_clean <- clin %>%
  distinct(cases.case_id, .keep_all = TRUE)

# ----------------- 3. 合并数据 -----------------
# 现在两个表都是“一对一”了，合并就不会报 many-to-many 警告
clin_merge <- left_join(clin_clean, yingshe_clean, by = c("cases.case_id" = "case_id"))

# ----------------- 4. 智能去重与生存时间计算 -----------------
clin_final <- clin_merge %>%
  # 再次检查去重 (防止万一还有重复)
  # 排序逻辑：
  # 1. 按病人 ID 分组
  # 2. 死亡 (Dead) 的排前面
  # 3. 随访时间 (days_to_last_follow_up) 长的排前面 (注意：此时已是数值型，NA会自动排到最后)
  arrange(
    cases.case_id, 
    desc(demographic.vital_status == "Dead"), 
    desc(diagnoses.days_to_last_follow_up)
  ) %>%
  distinct(cases.case_id, .keep_all = TRUE)

# ----------------- 5. 构建最终生存数据框 -----------------
survival_data <- clin_final %>%
  dplyr::select(
    sample_id = cases.submitter_id,
    case_id   = cases.case_id,
    OS_status = demographic.vital_status,
    follow_up = diagnoses.days_to_last_follow_up,
    death_time = demographic.days_to_death,
    age       = demographic.age_at_index,
    gender    = demographic.gender,
    stage     = diagnoses.ajcc_pathologic_stage
  ) %>%
  mutate(
    # 合并时间：如果是死人，优先用死亡时间；如果是活人，用随访时间
    OS.time = ifelse(!is.na(death_time), death_time, follow_up),
    # 转换状态：Dead -> 1, Alive -> 0
    OS = ifelse(OS_status == "Dead", 1, 0)
  ) %>%
  # 去除没有生存时间的样本
  filter(!is.na(OS.time) & OS.time > 0) %>%
  dplyr::select(-death_time, -follow_up) # 删除中间列

# 查看结果
head(survival_data)
cat("最终样本数:", nrow(survival_data), "\n")
#最终样本数: 443 

#gsva_groups<-read.csv("D:\\肺癌极化基因聚类\\TCGA肺癌数据分析\\只基于M1-M2集合进行打分\\基于肿瘤样本的打分\\TCGA-LUAD\\GSVA得分.csv",row.names = 1,check.names = F)
group_info <- data.frame(
  sample_id = rownames(gsva_scores),
  group = gsva_scores$Phenotype_Median
)

# 将分组信息合并到生存数据中
final_survival <- survival_data %>%
  left_join(group_info, by = "sample_id")

# 筛选出我们关注的两组 (M1_High_M2_Low 和 M1_Low_M2_High)
# 去掉 NA 和其他中间型
final_survival <- final_survival %>%
  filter(group %in% c("M1_High_M2_Low", "M1_Low_M2_High"))

# 检查数据
table(final_survival$group)

# M1_High_M2_Low M1_Low_M2_High 
# 45             49 

library(survival)
library(survminer)

# 1. 构建生存对象
# time = OS.time (天), event = OS (0/1)
surv_obj <- Surv(time = final_survival$OS.time, event = final_survival$OS)

# 2. 拟合曲线
fit <- survfit(surv_obj ~ group, data = final_survival)

# 3. 绘图
p <- ggsurvplot(
  fit, 
  data = final_survival,
  pval = TRUE,              # 显示 P 值
  conf.int = TRUE,          # 显示置信区间
  risk.table = TRUE,        # 显示风险表
  xlab = "Time (Days)",     # X轴标签
  ylab = "Overall Survival Probability", # Y轴标签
  title = "Survival Analysis: M1 vs M2 Phenotype",
  palette = c("#E41A1C", "#377EB8"), # 颜色对应 M1(红) M2(蓝)
  ggtheme = theme_minimal()
)

print(p)



#免疫浸润
expr<-LUAD_Tumor

group_info <- gsva_scores$Phenotype_Median
names(group_info) <- rownames(gsva_scores)
target_groups <- c("M1_High_M2_Low", "M1_Low_M2_High")


target_sample_ids <- rownames(gsva_scores[gsva_scores$Phenotype_Median %in% target_groups, ])

expr_target <- expr[, target_sample_ids]

lm22 <- read.table('./LM22.txt', header = TRUE, row.names = 1, 
                   sep = "\t", check.names = FALSE, quote = "")
lm22_genes <- rownames(lm22)  # 获取基因列表

common_genes <- intersect(rownames(expr_target), lm22_genes)

cat("匹配基因数：", length(common_genes), "\n")
tpm_filtered <- expr_target[common_genes, , drop=FALSE]

group_info <- gsva_scores$Phenotype_Median
names(group_info) <- rownames(gsva_scores)
target_groups <- c("M1_High_M2_Low", "M1_Low_M2_High")
keep_samples <- intersect(colnames(tpm_filtered), names(group_info)[group_info %in% target_groups])
tpm_filtered <- tpm_filtered[, keep_samples, drop=FALSE]
write.table(tpm_filtered, file = "免疫浸润分析表达矩阵(tpm).txt", 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

source("./CIBERSORT.R")
CIBERSORT_results <- CIBERSORT('./LM22.txt', '免疫浸润分析表达矩阵(tpm).txt', 
                               perm = 1000, QN = FALSE)

library(tidyverse)
library(ggpubr)
library(pheatmap)
library(corrplot)
library(ggsci)

CIBERSORT_results <- read.table("CIBERSORT-Results.txt",sep="\t",header=T)
rownames(CIBERSORT_results)<-CIBERSORT_results[,1]
CIBERSORT_results<-CIBERSORT_results[,-1]
head(CIBERSORT_results)
cell_cols <- setdiff(colnames(CIBERSORT_results), c("P-value", "Correlation", "RMSE"))
cell_scores <- CIBERSORT_results[, cell_cols, drop = FALSE]
group_info <- gsva_scores$Phenotype_Median
names(group_info) <- rownames(gsva_scores)

# 确保样本名一致
common_samples <- intersect(rownames(cell_scores), names(group_info))
cell_scores <- cell_scores[common_samples, ]
group_filtered <- group_info[common_samples]

# 创建包含分组的数据框（长格式）
cell_df <- as.data.frame(cell_scores) %>%
  rownames_to_column("sample") %>%
  mutate(Group = group_filtered[sample]) %>%
  pivot_longer(cols = all_of(cell_cols), names_to = "CellType", values_to = "Fraction")





cibersort_scores <- CIBERSORT_results[, 1:22]  # 根据实际列数调整

# 确保样本名一致，并与分组信息对齐
common_samples <- intersect(rownames(cibersort_scores), names(group_info))
cibersort_filtered <- cibersort_scores[common_samples, ]
group_filtered <- group_info[common_samples]

# 合并分组信息
cibersort_df <- cibersort_filtered %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(Group = group_filtered[sample])

cell_types <- colnames(cibersort_filtered)
p_values <- sapply(cell_types, function(ct) {
  # reformulate 自动处理变量名中的空格（相当于 `ct`）
  form <- reformulate("Group", response = ct)
  wilcox.test(form, data = cibersort_df)$p.value
})
p_adj <- p.adjust(p_values, method = "BH")  # 可选：多重检验校正

# 筛选 p < 0.05 的细胞类型（您也可以改用 p_adj < 0.05）
sig_cells <- names(p_values)[p_values < 0.05]
cat("显著差异的细胞类型（p < 0.05）：\n")
print(sig_cells)





# ---------- 3. 转换为长格式，仅保留显著细胞 ----------
plot_data <- cibersort_df %>%
  dplyr::select(sample, Group, all_of(sig_cells)) %>%
  pivot_longer(cols = -c(sample, Group), 
               names_to = "CellType", 
               values_to = "Fraction")

# ---------- 4. 绘制箱线图 ----------
library(ggplot2)
library(ggpubr)

# 确保已加载必要的包
library(ggplot2)
library(ggpubr)

# ---------- 可选：按 p 值排序细胞类型（使图形更美观） ----------
# 计算 p 值排序（假设您已计算出 p_values 向量）
sig_order <- names(sort(p_values[sig_cells]))  # 按 p 值升序排列
plot_data$CellType <- factor(plot_data$CellType, levels = sig_order)

# ---------- 绘制不分面的箱线图 ----------
# ---------- 绘制不分面的箱线图（无散点）----------
p_sig <- ggplot(plot_data, aes(x = CellType, y = Fraction, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.8) +   # 箱线图，无异常点标记
  # geom_jitter 已移除，取消散点
  scale_fill_manual(values = c("M1_High_M2_Low" = "#1E90FF", 
                               "M1_Low_M2_High" = "#FFC0CB")) +
  labs(x = NULL, y = "Relative Proportion", 
       title = "Significant Immune Cell Infiltration Differences") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  stat_compare_means(aes(group = Group), 
                     label = "p.signif",          # 显示星号：ns, *, **, ***
                     method = "wilcox.test",
                     size = 5,
                     label.y.npc = 0.95,          # 标记位置（箱线图上方）
                     hide.ns = TRUE)               # 隐藏不显著的（但此处均为显著，故全显示）

# 查看图形
print(p_sig)

# 保存为 PDF
ggsave("显著免疫细胞两组比较箱线图.pdf", plot = p_sig, width = 14, height = 8)


#差异基因热图
#分组差异分析
library(limma)
library(dplyr)
target_groups <- c("M1_High_M2_Low", "M1_Low_M2_High")
target_samples <- gsva_scores[gsva_scores$Phenotype_Median %in% target_groups, ]
sample_ids <- rownames(target_samples)


group_vector <- factor(
  target_samples$Phenotype_Median,
  levels = c("M1_Low_M2_High", "M1_High_M2_Low")  # 参照组设为 M2 主导，对比 M1 vs M2
)
names(group_vector) <- sample_ids

expr_matrix_log<-log2(as.matrix(expr) + 1)
expr_subset <- expr_matrix_log[, sample_ids, drop = FALSE]

expr_subset <- expr_subset[, names(group_vector)]
design <- model.matrix(~0 + group_vector)
colnames(design) <- c("M1_Low_M2_High", "M1_High_M2_Low")

# 拟合模型
fit <- lmFit(expr_subset, design)

# 对比：M1_High_M2_Low vs M1_Low_M2_High
contrast_matrix <- makeContrasts(
  M1_vs_M2 = M1_High_M2_Low - M1_Low_M2_High,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 提取所有基因结果
deg_results <- topTable(fit2, coef = "M1_vs_M2", number = Inf, adjust.method = "BH")

# 添加显著性标签
deg_results$Significance <- ifelse(
  deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 1,
  ifelse(deg_results$logFC > 0, "Up_in_M1", "Up_in_M2"),
  "Not_Sig"
)

# 筛选显著差异基因
deg_sig <- deg_results %>% filter(Significance != "Not_Sig")
#445个

cat("总基因数:", nrow(deg_results), "\n")
cat("M1 高表达基因数 (Up_in_M1):", sum(deg_results$Significance == "Up_in_M1"), "\n")
cat("M2 高表达基因数 (Up_in_M2):", sum(deg_results$Significance == "Up_in_M2"), "\n")
expr_gene<-expr[rownames(deg_sig),]

expr_log <- log2(as.matrix(expr_gene) + 1)


target_groups <- c("M1_High_M2_Low", "M1_Low_M2_High")

target_samples_info <- gsva_scores[gsva_scores$Phenotype_Median %in% target_groups, ]

target_sample_ids <- intersect(rownames(target_samples_info), colnames(expr_log))

target_sample_ids <- target_sample_ids[order(target_samples_info$Phenotype_Median[match(target_sample_ids, rownames(target_samples_info))])]
expr_plot <- expr_log[, target_sample_ids, drop = FALSE]



annotation_col <- data.frame(
  Group = target_samples_info$Phenotype_Median[match(colnames(expr_plot), rownames(target_samples_info))]
)
rownames(annotation_col) <- colnames(expr_plot)


ann_colors <- list(
  Group = c(
    "M1_High_M2_Low" = "#E41A1C", # 红色代表 M1 主导
    "M1_Low_M2_High" = "#377EB8"  # 蓝色代表 M2 主导
  )
)



pdf("GSVA分组差异基因热图-LUSC.pdf", width = 12, height = 10) # 设置 PDF 尺寸
pheatmap(expr_plot,
         annotation_col = annotation_col,      # 添加顶部注释
         annotation_colors = ann_colors,       # 注释颜色
         scale = "row",                        # 按行标准化 (Z-score)
         clustering_method = "complete",       
         cluster_cols = FALSE,                 # 不聚类列（样本按顺序排）
         cluster_rows = FALSE,                 # 不聚类行（基因按原始顺序排，无树状图）
         show_rownames = T,                # 不显示基因名
         show_colnames = FALSE,                # 不显示样本名
         main = "M1 vs M2 DEGs Expression (Target Samples)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
)
dev.off()



# 1. 从 deg_results 中提取上调/下调基因，并按 logFC 排序取 top30
up_genes <- deg_results %>% 
  filter(Significance == "Up_in_M1") %>% 
  arrange(desc(logFC)) %>%          # 按 logFC 降序（上调最显著在前）
  head(30) %>% 
  rownames()

down_genes <- deg_results %>% 
  filter(Significance == "Up_in_M2") %>% 
  arrange(logFC) %>%                 # 按 logFC 升序（下调最显著在前，负值越小越显著）
  head(30) %>% 
  rownames()

# 合并基因，并定义展示顺序：上调基因在前，下调基因在后（各自保持排序）
selected_genes <- c(up_genes, down_genes)

# 2. 提取这些基因在目标样本中的表达矩阵（已取过 log2(tpm+1)）
expr_plot <- expr_log[selected_genes, target_sample_ids, drop = FALSE]

# 3. 如果需要，按我们设定的基因顺序重新排列行（防止聚类打乱顺序）
# 此处保留 selected_genes 的顺序，并设置 cluster_rows = FALSE
# 同时可对列按分组排序（已在之前完成 target_sample_ids 的排序）

# 4. 绘制热图，增加行高以适应60个基因
pdf("GSVA分组差异基因热图_top30_up_down_LUSC.pdf", width = 12, height = 14)  # 增加高度
pheatmap(expr_plot,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         scale = "row",
         cluster_rows = FALSE,          # 不聚类行，保持我们指定的顺序
         cluster_cols = FALSE,
         show_rownames = TRUE,          # 显示基因名（60个基因可以看清）
         show_colnames = FALSE,
         main = "Top 30 Up/Down-regulated DEGs (M1 vs M2)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         fontsize_row = 10               # 适当缩小行字体，避免重叠
)
dev.off()




##富集分析

# 加载必备包
# 加载必备包
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggpubr)
library(dplyr)

# ---------------------- 1. 定义自定义背景基因集（LUAD表达谱所有基因） ----------------------
bg_genes <- rownames(LUAD_Tumor)  # 你的LUAD肿瘤表达矩阵所有基因
bg_entrez <- mapIds(org.Hs.eg.db, keys = bg_genes, keytype = "SYMBOL", column = "ENTREZID")
bg_entrez <- bg_entrez[!is.na(bg_entrez)]  # 移除无法映射的基因
cat("自定义背景基因数：", length(bg_entrez), "\n")  # 验证背景基因数

# ---------------------- 2. 提取上下调差异基因 ----------------------
up_genes <- rownames(deg_results[deg_results$Significance == "Up_in_M1", ])  # M1高表达
down_genes <- rownames(deg_results[deg_results$Significance == "Up_in_M2", ])  # M2高表达

# 转换为ENTREZ ID
up_entrez <- mapIds(org.Hs.eg.db, keys = up_genes, keytype = "SYMBOL", column = "ENTREZID")
down_entrez <- mapIds(org.Hs.eg.db, keys = down_genes, keytype = "SYMBOL", column = "ENTREZID")
up_entrez <- up_entrez[!is.na(up_entrez)]
down_entrez <- down_entrez[!is.na(down_entrez)]

# ---------------------- 3. GO-BP富集分析（指定自定义背景） ----------------------
go_up <- enrichGO(
  gene = up_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # 生物学过程
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE,
  universe = bg_entrez  # 自定义背景
)
go_down <- enrichGO(
  gene = down_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE,
  universe = bg_entrez
)

# ---------------------- 4. KEGG富集分析（指定自定义背景） ----------------------
kegg_up <- enrichKEGG(
  gene = up_entrez,
  organism = "hsa",  # 人类
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = bg_entrez
)
kegg_down <- enrichKEGG(
  gene = down_entrez,
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = bg_entrez
)

# ---------------------- 5. 数据预处理：计算-log10(p.adjust)并提取前10条 ----------------------
# GO上调（M1）：添加-log10(p.adjust)列
go_up_df <- as.data.frame(go_up) %>% 
  mutate(log10_padj = -log10(p.adjust)) %>%  # 核心：计算-log10(p.adjust)
  arrange(desc(log10_padj)) %>%  # 按显著性降序排序
  head(10)

# GO下调（M2）
go_down_df <- as.data.frame(go_down) %>% 
  mutate(log10_padj = -log10(p.adjust)) %>% 
  arrange(desc(log10_padj)) %>% 
  head(10)

# KEGG上调（M1）
kegg_up_df <- as.data.frame(kegg_up) %>% 
  mutate(log10_padj = -log10(p.adjust)) %>% 
  arrange(desc(log10_padj)) %>% 
  head(10)

# KEGG下调（M2）
kegg_down_df <- as.data.frame(kegg_down) %>% 
  mutate(log10_padj = -log10(p.adjust)) %>% 
  arrange(desc(log10_padj)) %>% 
  head(10)

# ---------------------- 6. 配色方案（匹配你的图） ----------------------
purple_palette <- c("#F9EEF1", "#CFA8BF", "#A873A4", "#7B497F", "#6D70A6")  # M1上调（紫）
blue_palette <- c("#D3E0E8", "#8EA2C5", "#3976AD", "#8787B2", "#5E6087")    # M2上调（蓝）

# ---------------------- 7. GO水平条形图（去边框 + 保留水平内网线） ----------------------
# M1上调GO水平条形图
p_go_up <- ggplot(go_up_df, 
                  aes(y = reorder(Description, log10_padj), x = log10_padj)) +  # X轴=log10_padj
  geom_bar(stat = "identity", aes(fill = Description), color = "white", width = 0.7) +
  scale_fill_manual(values = rep(purple_palette, 2)[1:10]) +  # 循环紫色系
  labs(x = "-log10(Adjusted P-value)", y = "GO Biological Process", 
       title = "GO Enrichment (Up-regulated in M1 Phenotype)") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",  # 隐藏图例
    panel.grid = element_blank() ,	# 隐藏网格线
   panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank()
	)

# M2上调GO水平条形图
p_go_down <- ggplot(go_down_df, 
                    aes(y = reorder(Description, log10_padj), x = log10_padj)) +
  geom_bar(stat = "identity", aes(fill = Description), color = "white", width = 0.7) +
  scale_fill_manual(values = rep(blue_palette, 2)[1:10]) +  # 循环蓝色系
  labs(x = "-log10(Adjusted P-value)", y = "GO Biological Process", 
       title = "GO Enrichment (Up-regulated in M2 Phenotype)") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),
	 panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank()
	
  )

# 拼接GO图
go_combined <- ggarrange(p_go_up, p_go_down, ncol = 1, nrow = 2, heights = c(1, 1))
ggsave("GO_Enrichment_Horizontal_Barplot_log10padj.pdf", go_combined, width = 12, height = 16)

# ---------------------- 8. KEGG水平条形图（X轴=-log10(p.adjust)） ----------------------
# M1上调KEGG水平条形图
p_kegg_up <- ggplot(kegg_up_df, 
                    aes(y = reorder(Description, log10_padj), x = log10_padj)) +
  geom_bar(stat = "identity", aes(fill = Description), color = "white", width = 0.7) +
  scale_fill_manual(values = rep(purple_palette, 2)[1:10]) +
  labs(x = "-log10(Adjusted P-value)", y = "KEGG Pathway", 
       title = "KEGG Enrichment (Up-regulated in M1 Phenotype)") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),
	panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank()
  )

# M2上调KEGG水平条形图
p_kegg_down <- ggplot(kegg_down_df, 
                      aes(y = reorder(Description, log10_padj), x = log10_padj)) +
  geom_bar(stat = "identity", aes(fill = Description), color = "white", width = 0.7) +
  scale_fill_manual(values = rep(blue_palette, 2)[1:10]) +
  labs(x = "-log10(Adjusted P-value)", y = "KEGG Pathway", 
       title = "KEGG Enrichment (Up-regulated in M2 Phenotype)") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),
	panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank()
  )

# 拼接KEGG图
kegg_combined <- ggarrange(p_kegg_up, p_kegg_down, ncol = 1, nrow = 2, heights = c(1, 1))
ggsave("KEGG_Enrichment_Horizontal_Barplot_log10padj.pdf", kegg_combined, width = 12, height = 16)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#2
#TCGA-LUAD-LUSC\\M1M2C5C8\\中位数\\不筛选正常肿瘤差异-直接对肿瘤样本打分分组生存\\免疫浸润


setwd("C:\\Users\\123\\Desktop\\RESULT\\TCGA-LUAD-LUSC\\M1M2C5C8\\中位数\\不筛选正常肿瘤差异-直接对肿瘤样本打分分组生存\\免疫浸润\\")

LUSC_Tumor<-read.csv("F:\\RESULT\\表达数据\\LUSC肿瘤样本基因表达谱(tpm).csv",row.names=1,check.names=F)
LUAD_Tumor<-read.csv("F:\\RESULT\\表达数据\\LUAD肿瘤样本基因表达谱(tpm).csv",row.names=1,check.names=F)
expr<-cbind(LUSC_Tumor,LUAD_Tumor)

C5_marker<-unlist(read.table("D:\\肺癌极化基因聚类\\C5-marker.txt"))
C8_marker<-unlist(read.table("D:\\肺癌极化基因聚类\\C8-marker.txt"))
M1_Polarization<-unlist(read.table("D:\\肺癌极化基因聚类\\M1_Polarization.txt"))
M2_Polarization<-unlist(read.table("D:\\肺癌极化基因聚类\\M2_Polarization.txt"))
geneA<-intersect(M1_Polarization,C8_marker)
geneB<-intersect(M2_Polarization,C5_marker)
gene_sets <- list(
  M1_C8_Marker = geneA,
  M2_C5_Marker = geneB)


library(GSVA)
library(dplyr)
expr_matrix <- as.matrix(expr)

gsvaPar <- gsvaParam(
  exprData = expr_matrix,
  geneSets = gene_sets,
  kcdf = "Gaussian"         # 连续型表达数据
)
gsva_scores <- gsva(gsvaPar, verbose = FALSE)
gsva_scores <- as.data.frame(t(gsva_scores))


m1_median <- median(gsva_scores$M1_C8_Marker)
m2_median <- median(gsva_scores$M2_C5_Marker)
gsva_scores$M1_Group <- ifelse(gsva_scores$M1_C8_Marker >= m1_median, "M1_High", "M1_Low")
gsva_scores$M2_Group <- ifelse(gsva_scores$M2_C5_Marker >= m2_median, "M2_High", "M2_Low")
gsva_scores$Phenotype_Median <- paste(gsva_scores$M1_Group, gsva_scores$M2_Group, sep = "_")  # 专属列名
table(gsva_scores$Phenotype_Median)






library(readr)
library(tidyverse)
library(dplyr)
library(data.table)
clin_LUAD<- fread("D:\\TCGA数据下载\\TCGA-LUAD\\TCGA-LUAD-Tumor\\clinical.tsv",data.table = F)
clin_LUSC<- fread("D:\\TCGA数据下载\\TCGA-LUSC\\TCGA-LUSC-Tumor\\clinical.tsv",data.table = F)
clin<-rbind(clin_LUAD,clin_LUSC)

#提取metadata
LUAD_tumor_sample_info<-read.csv("D:\\肺癌极化基因聚类\\TCGA肺癌数据分析\\TCGA-LUAD\\LUAD-Tumor\\TCGA_LUAD_sample_info.csv",row.names=1,check.names = F)
LUSC_tumor_sample_info<-read.csv("D:\\肺癌极化基因聚类\\TCGA肺癌数据分析\\TCGA-LUSC\\LUSC-Tumor\\TCGA_LUSC_sample_info.csv",row.names=1,check.names = F)
 
metadata<-rbind(LUAD_tumor_sample_info,LUSC_tumor_sample_info)
yingshe<-data.frame(metadata$sample_id,metadata$TCGA_case_id)
colnames(yingshe)<-c("cases.submitter_id","case_id")

head(yingshe)


clin <- clin[, c(
  "cases.case_id",
  "demographic.vital_status",
  "demographic.days_to_death",
  "diagnoses.days_to_last_follow_up",
  "demographic.age_at_index",
  "demographic.gender",
  "diagnoses.ajcc_pathologic_t",
  "diagnoses.ajcc_pathologic_m",
  "diagnoses.ajcc_pathologic_n",
  "diagnoses.ajcc_pathologic_stage"
)]
head(clin)
dim(yingshe)
dim(clin)
library(dplyr)
library(tidyr)

# ----------------- 1. 预处理 yingshe 表 (样本映射表) -----------------
# 防止 yingshe 中一个样本对应多个 case_id (虽然理论上应该是一对一，但清洗一下更保险)
yingshe_clean <- yingshe %>%
  distinct(cases.submitter_id, .keep_all = TRUE) # 按样本ID去重

# ----------------- 2. 预处理 clin 表 (临床数据) -----------------
# 2.1 强制转换数据类型
# 将随访时间转换为数值型。如果原数据是 "--" 或 "NA"，转换后会变成 NA
clin$diagnoses.days_to_last_follow_up <- as.numeric(clin$diagnoses.days_to_last_follow_up)
clin$demographic.days_to_death <- as.numeric(clin$demographic.days_to_death)

# 2.2 对 clin 表进行预去重
# 同样先按 case_id 去重，防止合并时出现多对多
# 这里先简单去重，保留第一条，具体的生存逻辑留到合并后做
clin_clean <- clin %>%
  distinct(cases.case_id, .keep_all = TRUE)

# ----------------- 3. 合并数据 -----------------
# 现在两个表都是“一对一”了，合并就不会报 many-to-many 警告
clin_merge <- left_join(clin_clean, yingshe_clean, by = c("cases.case_id" = "case_id"))

# ----------------- 4. 智能去重与生存时间计算 -----------------
clin_final <- clin_merge %>%
  # 再次检查去重 (防止万一还有重复)
  # 排序逻辑：
  # 1. 按病人 ID 分组
  # 2. 死亡 (Dead) 的排前面
  # 3. 随访时间 (days_to_last_follow_up) 长的排前面 (注意：此时已是数值型，NA会自动排到最后)
  arrange(
    cases.case_id, 
    desc(demographic.vital_status == "Dead"), 
    desc(diagnoses.days_to_last_follow_up)
  ) %>%
  distinct(cases.case_id, .keep_all = TRUE)

# ----------------- 5. 构建最终生存数据框 -----------------
survival_data <- clin_final %>%
  dplyr::select(
    sample_id = cases.submitter_id,
    case_id   = cases.case_id,
    OS_status = demographic.vital_status,
    follow_up = diagnoses.days_to_last_follow_up,
    death_time = demographic.days_to_death,
    age       = demographic.age_at_index,
    gender    = demographic.gender,
    stage     = diagnoses.ajcc_pathologic_stage
  ) %>%
  mutate(
    # 合并时间：如果是死人，优先用死亡时间；如果是活人，用随访时间
    OS.time = ifelse(!is.na(death_time), death_time, follow_up),
    # 转换状态：Dead -> 1, Alive -> 0
    OS = ifelse(OS_status == "Dead", 1, 0)
  ) %>%
  # 去除没有生存时间的样本
  filter(!is.na(OS.time) & OS.time > 0) %>%
  dplyr::select(-death_time, -follow_up) # 删除中间列

# 查看结果
head(survival_data)
cat("最终样本数:", nrow(survival_data), "\n")
#最终样本数: 913 


# 读取之前的 GSVA 分组结果 (假设你之前保存了包含 Phenotype_Median 的文件)
# 这里直接用你代码里的 metadata，它应该包含 sample_id 和分组信息
# 如果 metadata 里没有分组，请读取之前的 "GSVA得分.csv"
# 提取需要的列
group_info <- data.frame(
  sample_id = rownames(gsva_scores),
  group = gsva_scores$Phenotype_Median
)

# 将分组信息合并到生存数据中
final_survival <- survival_data %>%
  left_join(group_info, by = "sample_id")

# 筛选出我们关注的两组 (M1_High_M2_Low 和 M1_Low_M2_High)
# 去掉 NA 和其他中间型
final_survival <- final_survival %>%
  filter(group %in% c("M1_High_M2_Low", "M1_Low_M2_High"))

# 检查数据
table(final_survival$group)

#M1_High_M2_Low M1_Low_M2_High 
#           115            125 

library(survival)
library(survminer)

# 1. 构建生存对象
# time = OS.time (天), event = OS (0/1)
surv_obj <- Surv(time = final_survival$OS.time, event = final_survival$OS)

# 2. 拟合曲线
fit <- survfit(surv_obj ~ group, data = final_survival)

# 3. 绘图
p <- ggsurvplot(
  fit, 
  data = final_survival,
  pval = TRUE,              # 显示 P 值
  conf.int = TRUE,          # 显示置信区间
  risk.table = TRUE,        # 显示风险表
  xlab = "Time (Days)",     # X轴标签
  ylab = "Overall Survival Probability", # Y轴标签
  title = "Survival Analysis: M1 vs M2 Phenotype",
  palette = c("#E41A1C", "#377EB8"), # 颜色对应 M1(红) M2(蓝)
  ggtheme = theme_minimal()
)

print(p)


#免疫浸润


group_info <- gsva_scores$Phenotype_Median
names(group_info) <- rownames(gsva_scores)
target_groups <- c("M1_High_M2_Low", "M1_Low_M2_High")


target_sample_ids <- rownames(gsva_scores[gsva_scores$Phenotype_Median %in% target_groups, ])

expr_target <- expr[, target_sample_ids]

lm22 <- read.table('./LM22.txt', header = TRUE, row.names = 1, 
                   sep = "\t", check.names = FALSE, quote = "")
lm22_genes <- rownames(lm22)  # 获取基因列表

common_genes <- intersect(rownames(expr_target), lm22_genes)

cat("匹配基因数：", length(common_genes), "\n")
tpm_filtered <- expr_target[common_genes, , drop=FALSE]

group_info <- gsva_scores$Phenotype_Median
names(group_info) <- rownames(gsva_scores)
target_groups <- c("M1_High_M2_Low", "M1_Low_M2_High")
keep_samples <- intersect(colnames(tpm_filtered), names(group_info)[group_info %in% target_groups])
tpm_filtered <- tpm_filtered[, keep_samples, drop=FALSE]
write.table(tpm_filtered, file = "免疫浸润分析表达矩阵(tpm).txt", 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

source("./CIBERSORT.R")
CIBERSORT_results <- CIBERSORT('./LM22.txt', '免疫浸润分析表达矩阵(tpm).txt', 
                               perm = 1000, QN = FALSE)

library(tidyverse)
library(ggpubr)
library(pheatmap)
library(corrplot)
library(ggsci)

CIBERSORT_results <- read.table("CIBERSORT-Results.txt",sep="\t",header=T)
rownames(CIBERSORT_results)<-CIBERSORT_results[,1]
CIBERSORT_results<-CIBERSORT_results[,-1]
head(CIBERSORT_results)
cell_cols <- setdiff(colnames(CIBERSORT_results), c("P-value", "Correlation", "RMSE"))
cell_scores <- CIBERSORT_results[, cell_cols, drop = FALSE]
group_info <- gsva_scores$Phenotype_Median
names(group_info) <- rownames(gsva_scores)

# 确保样本名一致
common_samples <- intersect(rownames(cell_scores), names(group_info))
cell_scores <- cell_scores[common_samples, ]
group_filtered <- group_info[common_samples]

# 创建包含分组的数据框（长格式）
cell_df <- as.data.frame(cell_scores) %>%
  rownames_to_column("sample") %>%
  mutate(Group = group_filtered[sample]) %>%
  pivot_longer(cols = all_of(cell_cols), names_to = "CellType", values_to = "Fraction")





cibersort_scores <- CIBERSORT_results[, 1:22]  # 根据实际列数调整

# 确保样本名一致，并与分组信息对齐
common_samples <- intersect(rownames(cibersort_scores), names(group_info))
cibersort_filtered <- cibersort_scores[common_samples, ]
group_filtered <- group_info[common_samples]

# 合并分组信息
cibersort_df <- cibersort_filtered %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(Group = group_filtered[sample])

cell_types <- colnames(cibersort_filtered)
p_values <- sapply(cell_types, function(ct) {
  # reformulate 自动处理变量名中的空格（相当于 `ct`）
  form <- reformulate("Group", response = ct)
  wilcox.test(form, data = cibersort_df)$p.value
})
p_adj <- p.adjust(p_values, method = "BH")  # 可选：多重检验校正

# 筛选 p < 0.05 的细胞类型（您也可以改用 p_adj < 0.05）
sig_cells <- names(p_values)[p_values < 0.05]
cat("显著差异的细胞类型（p < 0.05）：\n")
print(sig_cells)





# ---------- 3. 转换为长格式，仅保留显著细胞 ----------
plot_data <- cibersort_df %>%
  dplyr::select(sample, Group, all_of(sig_cells)) %>%
  pivot_longer(cols = -c(sample, Group), 
               names_to = "CellType", 
               values_to = "Fraction")

# ---------- 4. 绘制箱线图 ----------
library(ggplot2)
library(ggpubr)

# 确保已加载必要的包
library(ggplot2)
library(ggpubr)

# ---------- 可选：按 p 值排序细胞类型（使图形更美观） ----------
# 计算 p 值排序（假设您已计算出 p_values 向量）
sig_order <- names(sort(p_values[sig_cells]))  # 按 p 值升序排列
plot_data$CellType <- factor(plot_data$CellType, levels = sig_order)

# ---------- 绘制不分面的箱线图 ----------
# ---------- 绘制不分面的箱线图（无散点）----------
p_sig <- ggplot(plot_data, aes(x = CellType, y = Fraction, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.8) +   # 箱线图，无异常点标记
  # geom_jitter 已移除，取消散点
  scale_fill_manual(values = c("M1_High_M2_Low" = "#1E90FF", 
                               "M1_Low_M2_High" = "#FFC0CB")) +
  labs(x = NULL, y = "Relative Proportion", 
       title = "Significant Immune Cell Infiltration Differences") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  stat_compare_means(aes(group = Group), 
                     label = "p.signif",          # 显示星号：ns, *, **, ***
                     method = "wilcox.test",
                     size = 5,
                     label.y.npc = 0.95,          # 标记位置（箱线图上方）
                     hide.ns = TRUE)               # 隐藏不显著的（但此处均为显著，故全显示）

# 查看图形
print(p_sig)



library(limma)
library(dplyr)
target_groups <- c("M1_High_M2_Low", "M1_Low_M2_High")
target_samples <- gsva_scores[gsva_scores$Phenotype_Median %in% target_groups, ]
sample_ids <- rownames(target_samples)


group_vector <- factor(
  target_samples$Phenotype_Median,
  levels = c("M1_Low_M2_High", "M1_High_M2_Low")  # 参照组设为 M2 主导，对比 M1 vs M2
)
names(group_vector) <- sample_ids

expr_matrix_log<-log2(as.matrix(expr) + 1)
expr_subset <- expr_matrix_log[, sample_ids, drop = FALSE]

expr_subset <- expr_subset[, names(group_vector)]
design <- model.matrix(~0 + group_vector)
colnames(design) <- c("M1_Low_M2_High", "M1_High_M2_Low")

# 拟合模型
fit <- lmFit(expr_subset, design)

# 对比：M1_High_M2_Low vs M1_Low_M2_High
contrast_matrix <- makeContrasts(
  M1_vs_M2 = M1_High_M2_Low - M1_Low_M2_High,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 提取所有基因结果
deg_results <- topTable(fit2, coef = "M1_vs_M2", number = Inf, adjust.method = "BH")

# 添加显著性标签
deg_results$Significance <- ifelse(
  deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 1,
  ifelse(deg_results$logFC > 0, "Up_in_M1", "Up_in_M2"),
  "Not_Sig"
)

# 筛选显著差异基因
deg_sig <- deg_results %>% filter(Significance != "Not_Sig")
#445个

cat("总基因数:", nrow(deg_results), "\n")
cat("M1 高表达基因数 (Up_in_M1):", sum(deg_results$Significance == "Up_in_M1"), "\n")
cat("M2 高表达基因数 (Up_in_M2):", sum(deg_results$Significance == "Up_in_M2"), "\n")
expr_gene<-expr[rownames(deg_sig),]

expr_log <- log2(as.matrix(expr_gene) + 1)


target_groups <- c("M1_High_M2_Low", "M1_Low_M2_High")

target_samples_info <- gsva_scores[gsva_scores$Phenotype_Median %in% target_groups, ]

target_sample_ids <- intersect(rownames(target_samples_info), colnames(expr_log))

target_sample_ids <- target_sample_ids[order(target_samples_info$Phenotype_Median[match(target_sample_ids, rownames(target_samples_info))])]
expr_plot <- expr_log[, target_sample_ids, drop = FALSE]



annotation_col <- data.frame(
  Group = target_samples_info$Phenotype_Median[match(colnames(expr_plot), rownames(target_samples_info))]
)
rownames(annotation_col) <- colnames(expr_plot)


ann_colors <- list(
  Group = c(
    "M1_High_M2_Low" = "#E41A1C", # 红色代表 M1 主导
    "M1_Low_M2_High" = "#377EB8"  # 蓝色代表 M2 主导
  )
)



pdf("GSVA分组差异基因热图-LUAD-LUSC.pdf", width = 12, height = 10) # 设置 PDF 尺寸
pheatmap(expr_plot,
         annotation_col = annotation_col,      # 添加顶部注释
         annotation_colors = ann_colors,       # 注释颜色
         scale = "row",                        # 按行标准化 (Z-score)
         clustering_method = "complete",       
         cluster_cols = FALSE,                 # 不聚类列（样本按顺序排）
         cluster_rows = FALSE,                 # 不聚类行（基因按原始顺序排，无树状图）
         show_rownames = T,                # 不显示基因名
         show_colnames = FALSE,                # 不显示样本名
         main = "M1 vs M2 DEGs Expression (Target Samples)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
)
dev.off()

up_genes <- deg_results %>% 
  filter(Significance == "Up_in_M1") %>% 
  arrange(desc(logFC)) %>%          # 按 logFC 降序（上调最显著在前）
  head(30) %>% 
  rownames()

down_genes <- deg_results %>% 
  filter(Significance == "Up_in_M2") %>% 
  arrange(logFC) %>%                 # 按 logFC 升序（下调最显著在前，负值越小越显著）
  head(30) %>% 
  rownames()

# 合并基因，并定义展示顺序：上调基因在前，下调基因在后（各自保持排序）
selected_genes <- c(up_genes, down_genes)

# 2. 提取这些基因在目标样本中的表达矩阵（已取过 log2(tpm+1)）
expr_plot <- expr_log[selected_genes, target_sample_ids, drop = FALSE]

# 3. 如果需要，按我们设定的基因顺序重新排列行（防止聚类打乱顺序）
# 此处保留 selected_genes 的顺序，并设置 cluster_rows = FALSE
# 同时可对列按分组排序（已在之前完成 target_sample_ids 的排序）

# 4. 绘制热图，增加行高以适应60个基因
pdf("GSVA分组差异基因热图_top30_up_down_LUAD_LUSC.pdf", width = 12, height = 14)  # 增加高度
pheatmap(expr_plot,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         scale = "row",
         cluster_rows = FALSE,          # 不聚类行，保持我们指定的顺序
         cluster_cols = FALSE,
         show_rownames = TRUE,          # 显示基因名（60个基因可以看清）
         show_colnames = FALSE,
         main = "Top 30 Up/Down-regulated DEGs (M1 vs M2)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         fontsize_row = 6               # 适当缩小行字体，避免重叠
)
dev.off()




##富集分析

# 加载必备包
# 加载必备包
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggpubr)
library(dplyr)

# ---------------------- 1. 定义自定义背景基因集（LUAD表达谱所有基因） ----------------------
bg_genes <- rownames(expr)  # 肿瘤表达矩阵所有基因
bg_entrez <- mapIds(org.Hs.eg.db, keys = bg_genes, keytype = "SYMBOL", column = "ENTREZID")
bg_entrez <- bg_entrez[!is.na(bg_entrez)]  # 移除无法映射的基因
cat("自定义背景基因数：", length(bg_entrez), "\n")  # 验证背景基因数

# ---------------------- 2. 提取上下调差异基因 ----------------------
up_genes <- rownames(deg_results[deg_results$Significance == "Up_in_M1", ])  # M1高表达
down_genes <- rownames(deg_results[deg_results$Significance == "Up_in_M2", ])  # M2高表达

# 转换为ENTREZ ID
up_entrez <- mapIds(org.Hs.eg.db, keys = up_genes, keytype = "SYMBOL", column = "ENTREZID")
down_entrez <- mapIds(org.Hs.eg.db, keys = down_genes, keytype = "SYMBOL", column = "ENTREZID")
up_entrez <- up_entrez[!is.na(up_entrez)]
down_entrez <- down_entrez[!is.na(down_entrez)]

# ---------------------- 3. GO-BP富集分析（指定自定义背景） ----------------------
go_up <- enrichGO(
  gene = up_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # 生物学过程
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE,
  universe = bg_entrez  # 自定义背景
)
go_down <- enrichGO(
  gene = down_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE,
  universe = bg_entrez
)

# ---------------------- 4. KEGG富集分析（指定自定义背景） ----------------------
kegg_up <- enrichKEGG(
  gene = up_entrez,
  organism = "hsa",  # 人类
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = bg_entrez
)
kegg_down <- enrichKEGG(
  gene = down_entrez,
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = bg_entrez
)

# ---------------------- 5. 数据预处理：计算-log10(p.adjust)并提取前10条 ----------------------
# GO上调（M1）：添加-log10(p.adjust)列
go_up_df <- as.data.frame(go_up) %>% 
  mutate(log10_padj = -log10(p.adjust)) %>%  # 核心：计算-log10(p.adjust)
  arrange(desc(log10_padj)) %>%  # 按显著性降序排序
  head(10)

# GO下调（M2）
go_down_df <- as.data.frame(go_down) %>% 
  mutate(log10_padj = -log10(p.adjust)) %>% 
  arrange(desc(log10_padj)) %>% 
  head(10)

# KEGG上调（M1）
kegg_up_df <- as.data.frame(kegg_up) %>% 
  mutate(log10_padj = -log10(p.adjust)) %>% 
  arrange(desc(log10_padj)) %>% 
  head(10)

# KEGG下调（M2）
kegg_down_df <- as.data.frame(kegg_down) %>% 
  mutate(log10_padj = -log10(p.adjust)) %>% 
  arrange(desc(log10_padj)) %>% 
  head(10)

# ---------------------- 6. 配色方案（匹配你的图） ----------------------
purple_palette <- c("#F9EEF1", "#CFA8BF", "#A873A4", "#7B497F", "#6D70A6")  # M1上调（紫）
blue_palette <- c("#D3E0E8", "#8EA2C5", "#3976AD", "#8787B2", "#5E6087")    # M2上调（蓝）

# ---------------------- 7. GO水平条形图（去边框 + 保留水平内网线） ----------------------
# M1上调GO水平条形图
p_go_up <- ggplot(go_up_df, 
                  aes(y = reorder(Description, log10_padj), x = log10_padj)) +  # X轴=log10_padj
  geom_bar(stat = "identity", aes(fill = Description), color = "white", width = 0.7) +
  scale_fill_manual(values = rep(purple_palette, 2)[1:10]) +  # 循环紫色系
  labs(x = "-log10(Adjusted P-value)", y = "GO Biological Process", 
       title = "GO Enrichment (Up-regulated in M1 Phenotype)") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",  # 隐藏图例
    panel.grid = element_blank() ,	# 隐藏网格线
   panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank()
	)

# M2上调GO水平条形图
p_go_down <- ggplot(go_down_df, 
                    aes(y = reorder(Description, log10_padj), x = log10_padj)) +
  geom_bar(stat = "identity", aes(fill = Description), color = "white", width = 0.7) +
  scale_fill_manual(values = rep(blue_palette, 2)[1:10]) +  # 循环蓝色系
  labs(x = "-log10(Adjusted P-value)", y = "GO Biological Process", 
       title = "GO Enrichment (Up-regulated in M2 Phenotype)") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),
	 panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank()
	
  )

# 拼接GO图
go_combined <- ggarrange(p_go_up, p_go_down, ncol = 1, nrow = 2, heights = c(1, 1))
ggsave("GO_Enrichment_Horizontal_Barplot_log10padj.pdf", go_combined, width = 12, height = 16)

# ---------------------- 8. KEGG水平条形图（X轴=-log10(p.adjust)） ----------------------
# M1上调KEGG水平条形图
p_kegg_up <- ggplot(kegg_up_df, 
                    aes(y = reorder(Description, log10_padj), x = log10_padj)) +
  geom_bar(stat = "identity", aes(fill = Description), color = "white", width = 0.7) +
  scale_fill_manual(values = rep(purple_palette, 2)[1:10]) +
  labs(x = "-log10(Adjusted P-value)", y = "KEGG Pathway", 
       title = "KEGG Enrichment (Up-regulated in M1 Phenotype)") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),
	panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank()
  )

# M2上调KEGG水平条形图
p_kegg_down <- ggplot(kegg_down_df, 
                      aes(y = reorder(Description, log10_padj), x = log10_padj)) +
  geom_bar(stat = "identity", aes(fill = Description), color = "white", width = 0.7) +
  scale_fill_manual(values = rep(blue_palette, 2)[1:10]) +
  labs(x = "-log10(Adjusted P-value)", y = "KEGG Pathway", 
       title = "KEGG Enrichment (Up-regulated in M2 Phenotype)") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),
	panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank()
  )

# 拼接KEGG图
kegg_combined <- ggarrange(p_kegg_up, p_kegg_down, ncol = 1, nrow = 2, heights = c(1, 1))
ggsave("KEGG_Enrichment_Horizontal_Barplot_log10padj.pdf", kegg_combined, width = 12, height = 16)




#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#TCGA-LUAD-LUSC\M1M2C5C8\四分位数\筛选正常-肿瘤差异-对肿瘤差异基因表达谱打分
#3

setwd("C:\\Users\\123\\Desktop\\RESULT\\TCGA-LUAD-LUSC\\M1M2C5C8\\四分位数\\筛选正常-肿瘤差异-对肿瘤差异基因表达谱打分\\免疫浸润")

deg_sig1<-read.csv("F:\\RESULT\\差异基因\\LUAD_LUSC_Tumor_vs_Normal_DEGs_Sig.csv",row.names=1,check.names=F)
LUAD_LUSC <- read.csv("D:\\肺癌极化基因聚类\\TCGA肺癌数据分析\\TCGA-LUAD-LUSC\\肿瘤样本基因表达谱(tpm).csv",row.names=1,check.names=F)

expr<- LUAD_LUSC[rownames(deg_sig1),]


C5_marker<-unlist(read.table("D:\\肺癌极化基因聚类\\C5-marker.txt"))
C8_marker<-unlist(read.table("D:\\肺癌极化基因聚类\\C8-marker.txt"))
M1_Polarization<-unlist(read.table("D:\\肺癌极化基因聚类\\M1_Polarization.txt"))
M2_Polarization<-unlist(read.table("D:\\肺癌极化基因聚类\\M2_Polarization.txt"))
geneA<-intersect(M1_Polarization,C8_marker)
geneB<-intersect(M2_Polarization,C5_marker)
gene_sets <- list(
  M1_C8_Marker = geneA,
  M2_C5_Marker = geneB)
  


library(GSVA)
library(dplyr)
expr_matrix <- as.matrix(expr)

gsvaPar <- gsvaParam(
  exprData = expr_matrix,
  geneSets = gene_sets,
  kcdf = "Gaussian"         # 连续型表达数据
)
gsva_scores <- gsva(gsvaPar, verbose = FALSE)
gsva_scores <- as.data.frame(t(gsva_scores))
q1_m1 <- quantile(gsva_scores$M1_C8_Marker, 0.25)
q3_m1 <- quantile(gsva_scores$M1_C8_Marker, 0.75)
q1_m2 <- quantile(gsva_scores$M2_C5_Marker, 0.25)
q3_m2 <- quantile(gsva_scores$M2_C5_Marker, 0.75)

# 根据四分位数分组
gsva_scores$M1_Group <- ifelse(
  gsva_scores$M1_C8_Marker >= q3_m1, "M1_High",
  ifelse(gsva_scores$M1_C8_Marker <= q1_m1, "M1_Low", "M1_Medium")
)
gsva_scores$M2_Group <- ifelse(
  gsva_scores$M2_C5_Marker >= q3_m2, "M2_High",
  ifelse(gsva_scores$M2_C5_Marker <= q1_m2, "M2_Low", "M2_Medium")
)

# 组合表型
gsva_scores$Phenotype_Quartile <- paste(
  gsva_scores$M1_Group, gsva_scores$M2_Group, sep = "_"
)
table(gsva_scores$Phenotype_Quartile)






library(readr)
library(tidyverse)
library(dplyr)
library(data.table)
clin_LUAD<- fread("D:\\TCGA数据下载\\TCGA-LUAD\\TCGA-LUAD-Tumor\\clinical.tsv",data.table = F)
clin_LUSC<- fread("D:\\TCGA数据下载\\TCGA-LUSC\\TCGA-LUSC-Tumor\\clinical.tsv",data.table = F)
clin<-rbind(clin_LUAD,clin_LUSC)

#提取metadata
LUAD_tumor_sample_info<-read.csv("D:\\肺癌极化基因聚类\\TCGA肺癌数据分析\\TCGA-LUAD\\LUAD-Tumor\\TCGA_LUAD_sample_info.csv",row.names=1,check.names = F)
LUSC_tumor_sample_info<-read.csv("D:\\肺癌极化基因聚类\\TCGA肺癌数据分析\\TCGA-LUSC\\LUSC-Tumor\\TCGA_LUSC_sample_info.csv",row.names=1,check.names = F)
 
metadata<-rbind(LUAD_tumor_sample_info,LUSC_tumor_sample_info)
yingshe<-data.frame(metadata$sample_id,metadata$TCGA_case_id)
colnames(yingshe)<-c("cases.submitter_id","case_id")

head(yingshe)


clin <- clin[, c(
  "cases.case_id",
  "demographic.vital_status",
  "demographic.days_to_death",
  "diagnoses.days_to_last_follow_up",
  "demographic.age_at_index",
  "demographic.gender",
  "diagnoses.ajcc_pathologic_t",
  "diagnoses.ajcc_pathologic_m",
  "diagnoses.ajcc_pathologic_n",
  "diagnoses.ajcc_pathologic_stage"
)]
head(clin)
dim(yingshe)
dim(clin)
library(dplyr)
library(tidyr)

# ----------------- 1. 预处理 yingshe 表 (样本映射表) -----------------
# 防止 yingshe 中一个样本对应多个 case_id (虽然理论上应该是一对一，但清洗一下更保险)
yingshe_clean <- yingshe %>%
  distinct(cases.submitter_id, .keep_all = TRUE) # 按样本ID去重

# ----------------- 2. 预处理 clin 表 (临床数据) -----------------
# 2.1 强制转换数据类型
# 将随访时间转换为数值型。如果原数据是 "--" 或 "NA"，转换后会变成 NA
clin$diagnoses.days_to_last_follow_up <- as.numeric(clin$diagnoses.days_to_last_follow_up)
clin$demographic.days_to_death <- as.numeric(clin$demographic.days_to_death)

# 2.2 对 clin 表进行预去重
# 同样先按 case_id 去重，防止合并时出现多对多
# 这里先简单去重，保留第一条，具体的生存逻辑留到合并后做
clin_clean <- clin %>%
  distinct(cases.case_id, .keep_all = TRUE)

# ----------------- 3. 合并数据 -----------------
# 现在两个表都是“一对一”了，合并就不会报 many-to-many 警告
clin_merge <- left_join(clin_clean, yingshe_clean, by = c("cases.case_id" = "case_id"))

# ----------------- 4. 智能去重与生存时间计算 -----------------
clin_final <- clin_merge %>%
  # 再次检查去重 (防止万一还有重复)
  # 排序逻辑：
  # 1. 按病人 ID 分组
  # 2. 死亡 (Dead) 的排前面
  # 3. 随访时间 (days_to_last_follow_up) 长的排前面 (注意：此时已是数值型，NA会自动排到最后)
  arrange(
    cases.case_id, 
    desc(demographic.vital_status == "Dead"), 
    desc(diagnoses.days_to_last_follow_up)
  ) %>%
  distinct(cases.case_id, .keep_all = TRUE)

# ----------------- 5. 构建最终生存数据框 -----------------
survival_data <- clin_final %>%
  dplyr::select(
    sample_id = cases.submitter_id,
    case_id   = cases.case_id,
    OS_status = demographic.vital_status,
    follow_up = diagnoses.days_to_last_follow_up,
    death_time = demographic.days_to_death,
    age       = demographic.age_at_index,
    gender    = demographic.gender,
    stage     = diagnoses.ajcc_pathologic_stage
  ) %>%
  mutate(
    # 合并时间：如果是死人，优先用死亡时间；如果是活人，用随访时间
    OS.time = ifelse(!is.na(death_time), death_time, follow_up),
    # 转换状态：Dead -> 1, Alive -> 0
    OS = ifelse(OS_status == "Dead", 1, 0)
  ) %>%
  # 去除没有生存时间的样本
  filter(!is.na(OS.time) & OS.time > 0) %>%
  dplyr::select(-death_time, -follow_up) # 删除中间列

# 查看结果
head(survival_data)
cat("最终样本数:", nrow(survival_data), "\n")
#最终样本数: 913 


# 读取之前的 GSVA 分组结果 (假设你之前保存了包含 Phenotype_Median 的文件)
# 这里直接用你代码里的 metadata，它应该包含 sample_id 和分组信息
# 如果 metadata 里没有分组，请读取之前的 "GSVA得分.csv"

# 提取需要的列
group_info <- data.frame(
  sample_id = rownames(gsva_scores),
  group = gsva_scores$Phenotype_Quartile
)

# 将分组信息合并到生存数据中
final_survival <- survival_data %>%
  left_join(group_info, by = "sample_id")

# 筛选出我们关注的两组 (M1_High_M2_Low 和 M1_Low_M2_High)
# 去掉 NA 和其他中间型
final_survival <- final_survival %>%
  filter(group %in% c("M1_High_M2_Low", "M1_Low_M2_High"))

# 检查数据
table(final_survival$group)

#M1_High_M2_Low M1_Low_M2_High 
#            15             13 

library(survival)
library(survminer)

# 1. 构建生存对象
# time = OS.time (天), event = OS (0/1)
surv_obj <- Surv(time = final_survival$OS.time, event = final_survival$OS)

# 2. 拟合曲线
fit <- survfit(surv_obj ~ group, data = final_survival)

# 3. 绘图
p <- ggsurvplot(
  fit, 
  data = final_survival,
  pval = TRUE,              # 显示 P 值
  conf.int = TRUE,          # 显示置信区间
  risk.table = TRUE,        # 显示风险表
  xlab = "Time (Days)",     # X轴标签
  ylab = "Overall Survival Probability", # Y轴标签
  title = "Survival Analysis: M1 vs M2 Phenotype",
  palette = c("#E41A1C", "#377EB8"), # 颜色对应 M1(红) M2(蓝)
  ggtheme = theme_minimal()
)



#免疫浸润


group_info <- gsva_scores$Phenotype_Quartile
names(group_info) <- rownames(gsva_scores)
target_groups <- c("M1_High_M2_Low", "M1_Low_M2_High")


target_sample_ids <- rownames(gsva_scores[gsva_scores$Phenotype_Quartile %in% target_groups, ])

expr_target <- expr[, target_sample_ids]

lm22 <- read.table('./LM22.txt', header = TRUE, row.names = 1, 
                   sep = "\t", check.names = FALSE, quote = "")
lm22_genes <- rownames(lm22)  # 获取基因列表

common_genes <- intersect(rownames(expr_target), lm22_genes)

cat("匹配基因数：", length(common_genes), "\n")
tpm_filtered <- expr_target[common_genes, , drop=FALSE]

group_info <- gsva_scores$Phenotype_Quartile
names(group_info) <- rownames(gsva_scores)
target_groups <- c("M1_High_M2_Low", "M1_Low_M2_High")
keep_samples <- intersect(colnames(tpm_filtered), names(group_info)[group_info %in% target_groups])
tpm_filtered <- tpm_filtered[, keep_samples, drop=FALSE]
write.table(tpm_filtered, file = "免疫浸润分析表达矩阵(tpm).txt", 
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

source("./CIBERSORT.R")
CIBERSORT_results <- CIBERSORT('./LM22.txt', '免疫浸润分析表达矩阵(tpm).txt', 
                               perm = 1000, QN = FALSE)

library(tidyverse)
library(ggpubr)
library(pheatmap)
library(corrplot)
library(ggsci)

CIBERSORT_results <- read.table("CIBERSORT-Results.txt",sep="\t",header=T)
rownames(CIBERSORT_results)<-CIBERSORT_results[,1]
CIBERSORT_results<-CIBERSORT_results[,-1]
head(CIBERSORT_results)
cell_cols <- setdiff(colnames(CIBERSORT_results), c("P-value", "Correlation", "RMSE"))
cell_scores <- CIBERSORT_results[, cell_cols, drop = FALSE]
group_info <- gsva_scores$Phenotype_Quartile
names(group_info) <- rownames(gsva_scores)

# 确保样本名一致
common_samples <- intersect(rownames(cell_scores), names(group_info))
cell_scores <- cell_scores[common_samples, ]
group_filtered <- group_info[common_samples]

# 创建包含分组的数据框（长格式）
cell_df <- as.data.frame(cell_scores) %>%
  rownames_to_column("sample") %>%
  mutate(Group = group_filtered[sample]) %>%
  pivot_longer(cols = all_of(cell_cols), names_to = "CellType", values_to = "Fraction")





cibersort_scores <- CIBERSORT_results[, 1:22]  # 根据实际列数调整

# 确保样本名一致，并与分组信息对齐
common_samples <- intersect(rownames(cibersort_scores), names(group_info))
cibersort_filtered <- cibersort_scores[common_samples, ]
group_filtered <- group_info[common_samples]

# 合并分组信息
cibersort_df <- cibersort_filtered %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(Group = group_filtered[sample])

cell_types <- colnames(cibersort_filtered)
p_values <- sapply(cell_types, function(ct) {
  # reformulate 自动处理变量名中的空格（相当于 `ct`）
  form <- reformulate("Group", response = ct)
  wilcox.test(form, data = cibersort_df)$p.value
})
p_adj <- p.adjust(p_values, method = "BH")  # 可选：多重检验校正

# 筛选 p < 0.05 的细胞类型（您也可以改用 p_adj < 0.05）
sig_cells <- names(p_values)[p_values < 0.05]
cat("显著差异的细胞类型（p < 0.05）：\n")
print(sig_cells)





# ---------- 3. 转换为长格式，仅保留显著细胞 ----------
plot_data <- cibersort_df %>%
  dplyr::select(sample, Group, all_of(sig_cells)) %>%
  pivot_longer(cols = -c(sample, Group), 
               names_to = "CellType", 
               values_to = "Fraction")

# ---------- 4. 绘制箱线图 ----------
library(ggplot2)
library(ggpubr)

# 确保已加载必要的包
library(ggplot2)
library(ggpubr)

# ---------- 可选：按 p 值排序细胞类型（使图形更美观） ----------
# 计算 p 值排序（假设您已计算出 p_values 向量）
sig_order <- names(sort(p_values[sig_cells]))  # 按 p 值升序排列
plot_data$CellType <- factor(plot_data$CellType, levels = sig_order)

# ---------- 绘制不分面的箱线图 ----------
# ---------- 绘制不分面的箱线图（无散点）----------
p_sig <- ggplot(plot_data, aes(x = CellType, y = Fraction, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.8) +   # 箱线图，无异常点标记
  # geom_jitter 已移除，取消散点
  scale_fill_manual(values = c("M1_High_M2_Low" = "#1E90FF", 
                               "M1_Low_M2_High" = "#FFC0CB")) +
  labs(x = NULL, y = "Relative Proportion", 
       title = "Significant Immune Cell Infiltration Differences") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold")) +
  stat_compare_means(aes(group = Group), 
                     label = "p.signif",          # 显示星号：ns, *, **, ***
                     method = "wilcox.test",
                     size = 5,
                     label.y.npc = 0.95,          # 标记位置（箱线图上方）
                     hide.ns = TRUE)               # 隐藏不显著的（但此处均为显著，故全显示）

# 查看图形
print(p_sig)












#差异基因热图
#分组差异分析
library(limma)
library(dplyr)
target_groups <- c("M1_High_M2_Low", "M1_Low_M2_High")
target_samples <- gsva_scores[gsva_scores$Phenotype_Quartile %in% target_groups, ]
sample_ids <- rownames(target_samples)


group_vector <- factor(
  target_samples$Phenotype_Quartile,
  levels = c("M1_Low_M2_High", "M1_High_M2_Low")  # 参照组设为 M2 主导，对比 M1 vs M2
)
names(group_vector) <- sample_ids

expr_matrix_log<-log2(as.matrix(expr) + 1)
expr_subset <- expr_matrix_log[, sample_ids, drop = FALSE]

expr_subset <- expr_subset[, names(group_vector)]
design <- model.matrix(~0 + group_vector)
colnames(design) <- c("M1_Low_M2_High", "M1_High_M2_Low")

# 拟合模型
fit <- lmFit(expr_subset, design)

# 对比：M1_High_M2_Low vs M1_Low_M2_High
contrast_matrix <- makeContrasts(
  M1_vs_M2 = M1_High_M2_Low - M1_Low_M2_High,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# 提取所有基因结果
deg_results <- topTable(fit2, coef = "M1_vs_M2", number = Inf, adjust.method = "BH")

# 添加显著性标签
deg_results$Significance <- ifelse(
  deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) > 1,
  ifelse(deg_results$logFC > 0, "Up_in_M1", "Up_in_M2"),
  "Not_Sig"
)

# 筛选显著差异基因
deg_sig <- deg_results %>% filter(Significance != "Not_Sig")
#445个

cat("总基因数:", nrow(deg_results), "\n")
cat("M1 高表达基因数 (Up_in_M1):", sum(deg_results$Significance == "Up_in_M1"), "\n")
cat("M2 高表达基因数 (Up_in_M2):", sum(deg_results$Significance == "Up_in_M2"), "\n")
expr_gene<-expr[rownames(deg_sig),]

expr_log <- log2(as.matrix(expr_gene) + 1)


target_groups <- c("M1_High_M2_Low", "M1_Low_M2_High")

target_samples_info <- gsva_scores[gsva_scores$Phenotype_Quartile %in% target_groups, ]

target_sample_ids <- intersect(rownames(target_samples_info), colnames(expr_log))

target_sample_ids <- target_sample_ids[order(target_samples_info$Phenotype_Quartile[match(target_sample_ids, rownames(target_samples_info))])]
expr_plot <- expr_log[, target_sample_ids, drop = FALSE]



annotation_col <- data.frame(
  Group = target_samples_info$Phenotype_Quartile[match(colnames(expr_plot), rownames(target_samples_info))]
)
rownames(annotation_col) <- colnames(expr_plot)


ann_colors <- list(
  Group = c(
    "M1_High_M2_Low" = "#E41A1C", # 红色代表 M1 主导
    "M1_Low_M2_High" = "#377EB8"  # 蓝色代表 M2 主导
  )
)



pdf("GSVA分组差异基因热图-LUAD-LUSC.pdf", width = 12, height = 10) # 设置 PDF 尺寸
pheatmap(expr_plot,
         annotation_col = annotation_col,      # 添加顶部注释
         annotation_colors = ann_colors,       # 注释颜色
         scale = "row",                        # 按行标准化 (Z-score)
         clustering_method = "complete",       
         cluster_cols = FALSE,                 # 不聚类列（样本按顺序排）
         cluster_rows = FALSE,                 # 不聚类行（基因按原始顺序排，无树状图）
         show_rownames = FALSE,                # 不显示基因名
         show_colnames = FALSE,                # 不显示样本名
         main = "M1 vs M2 DEGs Expression (Target Samples)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100)
)
dev.off()




up_genes <- deg_results %>% 
  filter(Significance == "Up_in_M1") %>% 
  arrange(desc(logFC)) %>%          # 按 logFC 降序（上调最显著在前）
  head(30) %>% 
  rownames()

down_genes <- deg_results %>% 
  filter(Significance == "Up_in_M2") %>% 
  arrange(logFC) %>%                 # 按 logFC 升序（下调最显著在前，负值越小越显著）
  head(30) %>% 
  rownames()

# 合并基因，并定义展示顺序：上调基因在前，下调基因在后（各自保持排序）
selected_genes <- c(up_genes, down_genes)

# 2. 提取这些基因在目标样本中的表达矩阵（已取过 log2(tpm+1)）
expr_plot <- expr_log[selected_genes, target_sample_ids, drop = FALSE]

# 3. 如果需要，按我们设定的基因顺序重新排列行（防止聚类打乱顺序）
# 此处保留 selected_genes 的顺序，并设置 cluster_rows = FALSE
# 同时可对列按分组排序（已在之前完成 target_sample_ids 的排序）

# 4. 绘制热图，增加行高以适应60个基因
pdf("GSVA分组差异基因热图_top30_up_down_LUAD_LUSC.pdf", width = 12, height = 14)  # 增加高度
pheatmap(expr_plot,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         scale = "row",
         cluster_rows = FALSE,          # 不聚类行，保持我们指定的顺序
         cluster_cols = FALSE,
         show_rownames = TRUE,          # 显示基因名（60个基因可以看清）
         show_colnames = FALSE,
         main = "Top 30 Up/Down-regulated DEGs (M1 vs M2)",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         fontsize_row = 10              # 适当缩小行字体，避免重叠
)
dev.off()




##富集分析

# 加载必备包
# 加载必备包
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggpubr)
library(dplyr)

# ---------------------- 1. 定义自定义背景基因集（LUAD表达谱所有基因） ----------------------
bg_genes <- rownames(expr)  # 肿瘤表达矩阵所有基因
bg_entrez <- mapIds(org.Hs.eg.db, keys = bg_genes, keytype = "SYMBOL", column = "ENTREZID")
bg_entrez <- bg_entrez[!is.na(bg_entrez)]  # 移除无法映射的基因
cat("自定义背景基因数：", length(bg_entrez), "\n")  # 验证背景基因数

# ---------------------- 2. 提取上下调差异基因 ----------------------
up_genes <- rownames(deg_results[deg_results$Significance == "Up_in_M1", ])  # M1高表达
down_genes <- rownames(deg_results[deg_results$Significance == "Up_in_M2", ])  # M2高表达

# 转换为ENTREZ ID
up_entrez <- mapIds(org.Hs.eg.db, keys = up_genes, keytype = "SYMBOL", column = "ENTREZID")
down_entrez <- mapIds(org.Hs.eg.db, keys = down_genes, keytype = "SYMBOL", column = "ENTREZID")
up_entrez <- up_entrez[!is.na(up_entrez)]
down_entrez <- down_entrez[!is.na(down_entrez)]

# ---------------------- 3. GO-BP富集分析（指定自定义背景） ----------------------
go_up <- enrichGO(
  gene = up_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # 生物学过程
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE,
  universe = bg_entrez  # 自定义背景
)
go_down <- enrichGO(
  gene = down_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE,
  universe = bg_entrez
)

# ---------------------- 4. KEGG富集分析（指定自定义背景） ----------------------
kegg_up <- enrichKEGG(
  gene = up_entrez,
  organism = "hsa",  # 人类
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = bg_entrez
)
kegg_down <- enrichKEGG(
  gene = down_entrez,
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe = bg_entrez
)

# ---------------------- 5. 数据预处理：计算-log10(p.adjust)并提取前10条 ----------------------
# GO上调（M1）：添加-log10(p.adjust)列
go_up_df <- as.data.frame(go_up) %>% 
  mutate(log10_padj = -log10(p.adjust)) %>%  # 核心：计算-log10(p.adjust)
  arrange(desc(log10_padj)) %>%  # 按显著性降序排序
  head(10)

# GO下调（M2）
go_down_df <- as.data.frame(go_down) %>% 
  mutate(log10_padj = -log10(p.adjust)) %>% 
  arrange(desc(log10_padj)) %>% 
  head(10)

# KEGG上调（M1）
kegg_up_df <- as.data.frame(kegg_up) %>% 
  mutate(log10_padj = -log10(p.adjust)) %>% 
  arrange(desc(log10_padj)) %>% 
  head(10)

# KEGG下调（M2）
kegg_down_df <- as.data.frame(kegg_down) %>% 
  mutate(log10_padj = -log10(p.adjust)) %>% 
  arrange(desc(log10_padj)) %>% 
  head(10)

# ---------------------- 6. 配色方案（匹配你的图） ----------------------
purple_palette <- c("#F9EEF1", "#CFA8BF", "#A873A4", "#7B497F", "#6D70A6")  # M1上调（紫）
blue_palette <- c("#D3E0E8", "#8EA2C5", "#3976AD", "#8787B2", "#5E6087")    # M2上调（蓝）

# ---------------------- 7. GO水平条形图（去边框 + 保留水平内网线） ----------------------
# M1上调GO水平条形图
p_go_up <- ggplot(go_up_df, 
                  aes(y = reorder(Description, log10_padj), x = log10_padj)) +  # X轴=log10_padj
  geom_bar(stat = "identity", aes(fill = Description), color = "white", width = 0.7) +
  scale_fill_manual(values = rep(purple_palette, 2)[1:10]) +  # 循环紫色系
  labs(x = "-log10(Adjusted P-value)", y = "GO Biological Process", 
       title = "GO Enrichment (Up-regulated in M1 Phenotype)") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",  # 隐藏图例
    panel.grid = element_blank() ,	# 隐藏网格线
   panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank()
	)

# M2上调GO水平条形图
p_go_down <- ggplot(go_down_df, 
                    aes(y = reorder(Description, log10_padj), x = log10_padj)) +
  geom_bar(stat = "identity", aes(fill = Description), color = "white", width = 0.7) +
  scale_fill_manual(values = rep(blue_palette, 2)[1:10]) +  # 循环蓝色系
  labs(x = "-log10(Adjusted P-value)", y = "GO Biological Process", 
       title = "GO Enrichment (Up-regulated in M2 Phenotype)") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),
	 panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank()
	
  )

# 拼接GO图
go_combined <- ggarrange(p_go_up, p_go_down, ncol = 1, nrow = 2, heights = c(1, 1))
ggsave("GO_Enrichment_Horizontal_Barplot_log10padj.pdf", go_combined, width = 12, height = 16)

# ---------------------- 8. KEGG水平条形图（X轴=-log10(p.adjust)） ----------------------
# M1上调KEGG水平条形图
p_kegg_up <- ggplot(kegg_up_df, 
                    aes(y = reorder(Description, log10_padj), x = log10_padj)) +
  geom_bar(stat = "identity", aes(fill = Description), color = "white", width = 0.7) +
  scale_fill_manual(values = rep(purple_palette, 2)[1:10]) +
  labs(x = "-log10(Adjusted P-value)", y = "KEGG Pathway", 
       title = "KEGG Enrichment (Up-regulated in M1 Phenotype)") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),
	panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank()
  )

# M2上调KEGG水平条形图
p_kegg_down <- ggplot(kegg_down_df, 
                      aes(y = reorder(Description, log10_padj), x = log10_padj)) +
  geom_bar(stat = "identity", aes(fill = Description), color = "white", width = 0.7) +
  scale_fill_manual(values = rep(blue_palette, 2)[1:10]) +
  labs(x = "-log10(Adjusted P-value)", y = "KEGG Pathway", 
       title = "KEGG Enrichment (Up-regulated in M2 Phenotype)") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),
	panel.border = element_blank(),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank()
  )

# 拼接KEGG图
kegg_combined <- ggarrange(p_kegg_up, p_kegg_down, ncol = 1, nrow = 2, heights = c(1, 1))
ggsave("KEGG_Enrichment_Horizontal_Barplot_log10padj.pdf", kegg_combined, width = 12, height = 16)
