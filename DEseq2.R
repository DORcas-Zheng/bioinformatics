library(DESeq2)

## 数据预处理
sampleNames <- c("F_1", "F_2", "M_1", "M_2")
# 第一行是命令信息，所以跳过
data <- read.table("C:/Users/99562/Desktop/all_feature.txt", header=TRUE, quote="\t", skip=1)
# 前六列分别是Geneid	Chr	Start	End	Strand	Length
# 我们要的是count数，所以从第七列开始
names(data)[7:10] <- sampleNames
countData <- as.matrix(data[7:10])
rownames(countData) <- data$Geneid
database <- data.frame(name=sampleNames, condition=c("F", "F", "M", "M"))
rownames(database) <- sampleNames

## 设置分组信息并构建dds对象
dds <- DESeqDataSetFromMatrix(countData, colData=database, design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ]

rld <- rlog(dds)
plotPCA(rld, intgroup=c("name","condition"))

## 使用DESeq函数估计离散度，然后差异分析获得res对象
dds <- DESeq(dds)
res <- results(dds)
write.csv(res, "res_des_output.csv")
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata, "all_des_output.csv", row.names=FALSE)

head(res)
summary(res)

plotMA(res, main="DESeq2", ylim=c(-2, 2))

library(ggplot2)

# 这里的resdata也可以用res_des_output.csv这个结果重新导入哦。
# 现在就是用的前面做DESeq的时候的resdata。
resdata$change <- as.factor(
  ifelse(
    resdata$padj<0.01 & abs(resdata$log2FoldChange)>1,
    ifelse(resdata$log2FoldChange>1, "Up", "Down"),
    "NoDiff"
  )
)
valcano <- ggplot(data=resdata, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
  geom_point(alpha=0.8, size=1) + 
  theme_bw(base_size=15) + 
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank()
  ) + 
  ggtitle("DESeq2 Valcano") + 
  scale_color_manual(name="", values=c("red", "green", "black"), limits=c("Up", "Down", "NoDiff")) + 
  geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
  geom_hline(yintercept=-log10(0.01), lty=2, col="gray", lwd=0.5)

valcano

sum(res$padj < 0.1, na.rm=TRUE)

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:500]
nt <- normTransform(dds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[,c("name","condition")])
pdf('heatmap1000.pdf',width = 6, height = 7)
pheatmap(log2.norm.counts, cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)
dev.off()




library(pheatmap)

select <- order(rowMeans(counts(dds, normalized=T)), decreasing=T)[500]
nt <- normTransform(dds)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[, c("name", "condition")])
pheatmap(log2.norm.counts, cluster_rows=F, show_rownames=F, cluster_cols=T, annotation_col=df, fontsize=4)