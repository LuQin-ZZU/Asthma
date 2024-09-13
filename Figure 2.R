library(limma)
library(pheatmap)
library(reshape2)
library(ggpubr)

#Figure 2A-B
expFile="m6aGeneExp.txt"

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
exp=data

Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))

sigVec=c()
sigGeneVec=c()
for(i in row.names(data)){
	test=wilcox.test(data[i,] ~ Type)
	pvalue=test$p.value
	Sig=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
	if(pvalue<0.05){
	sigVec=c(sigVec, paste0(i, Sig))
	sigGeneVec=c(sigGeneVec, i)}
}
a=data[sigGeneVec,]
outTab=rbind(ID=colnames(data), data)
write.table(outTab, file="diffGeneExp.txt", sep="\t", quote=F, col.names=F)
row.names(data)=sigVec

names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf("heatmap.pdf", width=7.5, height=4.7)
pheatmap(data,
         annotation=Type,
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

exp=as.data.frame(t(exp))
exp=cbind(exp, Type=Type)
data=melt(exp, id.vars=c("Type"))
colnames(data)=c("Type", "Gene", "Expression")

boxplot(data, x="Gene", y="Expression", color = "Type", 
	     xlab="",
	     ylab="Gene expression",
	     legend.title="Type",
	     palette = c("blue", "red"),
	     width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
	      method="wilcox.test",
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

pdf(file="boxplot.pdf", width=7.5, height=5)
print(p1)
dev.off()

#Figure 2C
library(corrplot)
library(circlize)

inputFile="diffGeneExp.txt"

data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

group=gsub("(.*)\\_(.*)", "\\2", colnames(data))
rt=t(data)

cor1=cor(rt)

col = c(rgb(1,0,0,seq(1,0,length=32)),rgb(0,0,1,seq(0,1,length=32)))
cor1[cor1==1]=0
c1 = ifelse(c(cor1)>=0,rgb(1,0,0,abs(cor1)),rgb(0,0,1,abs(cor1)))
col1 = matrix(c1,nc=ncol(rt))

pdf(file="corrplot.pdf", width=7, height=7)
corrplot(cor1,
         method = "pie",
         order = "hclust",
         type = "upper",
         col=colorRampPalette(c("blue", "white", "red"))(50)
)
dev.off()
