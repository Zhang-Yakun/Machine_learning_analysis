
# lasso -------------------------------------------------------------------


library(glmnet)
expr_matrix<-dataset3[,3:15]%>%as.matrix()
outcome <- factor(dataset3$Var)
cv_fit <- cv.glmnet(x = expr_matrix, 
                    y = outcome,
                    family = "binomial", # 请根据您的结局类型修改
                    alpha = 1,          # 1为LASSO，0为Ridge，中间值为Elastic Net
                    nfolds = 5)         # 5折交叉验证

# 3. 查看最优lambda值 (lambda.min: 误差最小；lambda.1se: 最简单且误差在1个标准差内)
plot(cv_fit) # 绘制误差路径图
print(cv_fit$lambda.min)
print(cv_fit$lambda.1se)

# 4. 提取在最优lambda下系数不为零的基因
# 通常使用 lambda.1se 可以获得更简洁的模型
coef_matrix <- coef(cv_fit, s = "lambda.1se")
selected_genes <- rownames(coef_matrix)[which(coef_matrix != 0)][-1] # 去除截距项
print(paste("筛选出的基因：", paste(selected_genes, collapse=", ")))

# PCA ---------------------------------------------------------------------

expression_scaled <- scale(expr_matrix, center = TRUE, scale = FALSE)

# 2. 执行PCA分析
pca_result <- prcomp(expression_scaled, center = FALSE, scale. = FALSE)

# 查看PCA结果摘要
cat("PCA分析结果摘要:\n")
print(summary(pca_result))

# 3. 提取前6个主成分（相当于6个新特征）
# 每个样本在这些主成分上的得分
pca_scores <- pca_result$x[, 1:6]

# 查看降维后的数据
cat("\n降维后的数据维度（前6个主成分）:", dim(pca_scores), "\n")
cat("主成分名称:", colnames(pca_scores), "\n\n")

# 4. 查看每个主成分的解释方差比例
variance_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
cat("前6个主成分解释的方差比例:\n")
for(i in 1:6) {
  cat(sprintf("PC%d: %.3f (%.1f%%)\n", 
              i, variance_explained[i], variance_explained[i]*100))
}
cat(sprintf("累计解释方差: %.3f (%.1f%%)\n", 
            sum(variance_explained[1:6]), sum(variance_explained[1:6])*100))

# 5. 查看主成分的载荷（loadings），了解每个基因对主成分的贡献
pca_loadings <- pca_result$rotation[, 1:6]
cat("\n前6个主成分的载荷（基因贡献）:\n")
print(round(pca_loadings, 4))

# 6. 可视化：前两个主成分的样本分布
library(ggplot2)

# 创建绘图数据框
plot_data <- data.frame(
  Sample = rownames(pca_scores),
  PC1 = pca_scores[, 1],
  PC2 = pca_scores[, 2],
  Group = factor(rep(c("Group1", "Group2"), each = 8))  # 示例分组，实际使用时请替换
)

# 绘制PCA图
pca_plot <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  geom_text(aes(label = Sample), vjust = -0.8, size = 3) +
  labs(
    title = "PCA Plot: 13 Genes Reduced to 6 Principal Components",
    x = paste0("PC1 (", round(variance_explained[1]*100, 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2]*100, 1), "%)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

print(pca_plot)

# 7. 保存降维后的数据到文件
write.csv(pca_scores, file = "PCA_reduced_data.csv", row.names = TRUE)

# 8. 碎石图：可视化每个主成分的解释方差
scree_data <- data.frame(
  PC = 1:length(variance_explained),
  Variance = variance_explained
)

scree_plot <- ggplot(scree_data[1:10, ], aes(x = factor(PC), y = Variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(aes(group = 1), color = "red", size = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    title = "Scree Plot: Variance Explained by Each Principal Component",
    x = "Principal Component",
    y = "Proportion of Variance Explained"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

print(scree_plot)

# 9. 加载贡献热图：显示每个基因对前6个主成分的贡献
library(pheatmap)

# 绘制载荷热图
pheatmap(
  pca_loadings,
  main = "Gene Loadings on First 6 Principal Components",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize = 10
)

# 输出总结
cat("\n========== PCA降维完成 ==========\n")
cat("原始数据: 16个样本 × 13个基因\n")
cat("降维后数据: 16个样本 × 6个主成分\n")
cat("输出文件: 'PCA_reduced_data.csv' 已保存\n")
