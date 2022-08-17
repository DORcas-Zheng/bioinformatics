#数据降维
# 主成分分析
library("Rtsne")
library("ggplot2")
X <- read.csv("BreastCancer.csv")
Y <- as.matrix(X[,1:3])
#标准化
Y <- scale(Y, center = T, scale = T)
# 
pca.res <- prcomp(Y, center = F, scale = F, rank. = 2)
X$PC1 <- pca.res$x[,1]
X$PC2 <- pca.res$x[,2]
ggplot(X, aes(x=PC1, y=PC2,color=Class)) + geom_point() + theme_bw()

