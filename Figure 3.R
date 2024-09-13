#Figure 3A
library(ConsensusClusterPlus)
expFile="diffGeneExp.txt"

data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=as.matrix(data)

group=sapply(strsplit(colnames(data),"\\_"), "[", 2)
data=data[,group=="treat"]

maxK=9
results=ConsensusClusterPlus(data,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              clusterAlg="pam",
              distance="euclidean",
              seed=123456,
              plot="pdf")

##Set the number of clusters
clusterNum=2
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("m6Acluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$m6Acluster))
cluster$m6Acluster=letter[match(cluster$m6Acluster, uniqClu)]
outTab=cbind(t(data), cluster)
outTab=rbind(ID=colnames(outTab), outTab)
write.table(outTab, file="m6Acluster.txt", sep="\t", quote=F, col.names=F)

#Figure 3B
library(limma)
library(ggplot2)

clusterFile="m6Acluster.txt"

rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
data=rt[,1:(ncol(rt)-1),drop=F]
m6Acluster=as.vector(rt[,ncol(rt)])

data.pca=prcomp(data)
pcaPredict=predict(data.pca)
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], m6Acluster=m6Acluster)
PCA.mean=aggregate(PCA[,1:2], list(m6Acluster=PCA$m6Acluster), mean)

bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
m6aCluCol=bioCol[1:length(levels(factor(m6Acluster)))]

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) {
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
df_ell <- data.frame()
for(g in levels(factor(PCA$m6Acluster))){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$m6Acluster==g,],
                                                   veganCovEllipse(cov.wt(cbind(PC1,PC2),
                                                                          wt=rep(1/length(PC1),length(PC1)))$cov,
                                                                   center=c(mean(PC1),mean(PC2))))), m6Acluster=g))
}

pdf(file="PCA.pdf", height=5, width=6.5)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = m6Acluster)) +
  scale_colour_manual(name="m6Acluster", values =m6aCluCol)+
  theme_bw()+
  theme(plot.margin=unit(rep(1.5,4),'lines'))+
  geom_path(data=df_ell, aes(x=PC1, y=PC2, colour=m6Acluster), size=1, linetype=2)+
  annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$m6Acluster, cex=7)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

#Figure 3C-D
library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)

clusterFile="m6Acluster.txt"

rt=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
rt=rt[order(rt$m6Acluster),]

data=t(rt[,1:(ncol(rt)-1),drop=F])
Type=rt[,ncol(rt),drop=F]

bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
m6aCluCol=bioCol[1:length(levels(factor(Type$m6Acluster)))]
names(m6aCluCol)=levels(factor(Type$m6Acluster))
ann_colors[["m6Acluster"]]=m6aCluCol

pdf("heatmap.pdf", width=7, height=4.5)
pp=pheatmap(data,
            annotation=Type,
            annotation_colors = ann_colors,
            color = colorRampPalette(c(rep("blue",2), "white", rep("red",2)))(100),
            cluster_cols =F,
            cluster_rows =T,
            scale="row",
            show_colnames=F,
            show_rownames=T,
            fontsize=6,
            fontsize_row=7,
            fontsize_col=6)

dev.off()

data=melt(rt, id.vars=c("m6Acluster"))
colnames(data)=c("m6Acluster", "Gene", "Expression")

p=ggboxplot(data, x="Gene", y="Expression", color = "m6Acluster", 
            ylab="Gene expression",
            xlab="",
            legend.title="m6Acluster",
            palette = m6aCluCol,
            width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=m6Acluster),
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

pdf(file="boxplot.pdf", width=6, height=5)
print(p1)
dev.off()

#Figure 3F
expFile="diffGeneExp.txt"

data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)

pca=prcomp(data, scale=TRUE)
value=predict(pca)
m6Ascore=value[,1]+value[,2]
m6Ascore=as.data.frame(m6Ascore)
scoreOut=rbind(id=colnames(m6Ascore), m6Ascore)
write.table(scoreOut, file="m6Ascore.txt", sep="\t", quote=F, col.names=F)

library(readxl)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
table <- read_excel("F:\\8.小论文\\3.FTO\\生信分析 - 副本\\19. High-low risk\\m6A评分差异分析\\m6AClusterm6AScore.xlsx")
table$ID <- as.factor(table$ID) # Transform values to factor
table$Type <- as.factor(table$Type)

cairo_ps("m6AClusterm6AScore.eps")
my_comparison <- list(c("A", "B"))
ggboxplot(table, x="Type", y="m6Ascore", 
          color = "Type", add = "jitter")+
  ylim(-8, 8)+
  theme(legend.position='right')+
  stat_compare_means(comparisons = my_comparison, method = "wilcox.test",
                     label = "p.signif", label.y=c(7))+
  xlab("")+
  ylab("m6A Score")+
  border(size=0.2)
dev.off()

#Figure 3G-H
library(limma)
library(reshape2)
library(ggpubr)
library(corrplot)

clusterFile="Cluster.txt"
immFile="CIBERSORT-Results.txt"

immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)

group=gsub("(.*)\\_(.*)", "\\2", row.names(immune))
data=immune[group=="treat",,drop=F]
data1=immune[group=="treat",,drop=F]

Cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(Cluster))
rt=cbind(data[sameSample,,drop=F], Cluster[sameSample,"Cluster",drop=F])
rt=rt[order(rt$Cluster, decreasing=F),]
conNum=nrow(rt[rt$Cluster=="A",])
treatNum=nrow(rt[rt$Cluster=="B",])

data=t(rt[,-ncol(rt)])
pdf(file="barplot.pdf", width=14.5, height=8)
col=rainbow(nrow(data), s=0.7, v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data, col=col, xaxt="n", yaxt="n", ylab="Relative Percent", cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
rect(xleft = a1[1]-0.5, ybottom = -0.01, xright = a1[conNum]+0.5, ytop = -0.06,col="green")
text(a1[conNum]/2,-0.035,"A",cex=2)
rect(xleft = a1[conNum]+0.5, ybottom = -0.01, xright =a1[length(a1)]+0.5 , ytop = -0.06,col="red")
text((a1[length(a1)]+a1[conNum])/2,-0.035,"B",cex=2)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3)
dev.off()

data=rt
data=melt(data, id.vars=c("Cluster"))
colnames(data)=c("Cluster", "Immune", "Expression")

group=levels(factor(data$Cluster))
bioCol=c("#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", color="Cluster",
                  xlab="",
                  ylab="Fraction",
                  legend.title="Cluster",
                  add="point",
                  width=0.8,
                  palette=bioCol)+
  rotate_x_text(50)+
  stat_compare_means(aes(group=Cluster),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")

pdf(file="immune.diff.pdf", width=8, height=6)
print(boxplot)
dev.off()