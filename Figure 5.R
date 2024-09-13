#Figure 5A-B
library(rms)
library(rmda)
library(Hmisc)

riskFile="risk.txt"

data=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

group=gsub("(.*?)\\_(.*)", "\\2", row.names(data))
row.names(data)=gsub("(.*?)\\_(.*)", "\\1", row.names(data))
sameSample=row.names(data)
outTab=cbind(data[sameSample,1:(ncol(data)-1),drop=F], data[sameSample,"riskScore",drop=F])
data=data[,2:(ncol(data)-1)]
rt=cbind(as.data.frame(data), Type=group)

ddist=datadist(rt)
options(datadist="ddist")

lrmModel=lrm(Type~ ., data=rt, x=T, y=T)
nomo=nomogram(lrmModel, fun=plogis,
	fun.at=c(0.0001,0.1,0.3,0.6,0.9,0.99),
	lp=F, funlabel="Risk of Disease")

pdf("Nomo.pdf", width=8, height=6)
plot(nomo, cex.axis=0.8)
dev.off()

nomoRisk=predict(lrmModel, type="fitted")
outTab=cbind(outTab, Nomogram=nomoRisk)
outTab=rbind(id=colnames(outTab), outTab)
write.table(outTab, file="nomoRisk.txt", sep="\t", quote=F, col.names=F)

cali=calibrate(lrmModel, method="boot", B=1000)
pdf("Calibration.pdf", width=5.5, height=5.5)
plot(cali,
     xlab="Predicted probability",
     ylab="Actual probability", sub=F)

group=ifelse(group=="con", 0, 1)
cindex=rcorrcens(group~nomoRisk)
se=cindex[,"SD"]/2
c_index=sprintf("%.03f", cindex[,"C"])
c_index.ci_low=sprintf("%.03f", cindex[,"C"]-(se*1.96))
c_index.ci_high=sprintf("%.03f", cindex[,"C"]+(se*1.96))
cindexLabel=paste0(c_index, " (95% CI: ", c_index.ci_low, "-", c_index.ci_high, ")")
text(0.1, 0.85, "C-index:")
text(0.28, 0.78, cindexLabel)
dev.off()

#Figure 5C-D
library(rms)
library(rmda)

inputFile="nomoRisk.txt"

rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)

riskScore=decision_curve(Type ~ riskScore, data=rt, 
                         family = binomial(link ='logit'),
                         thresholds= seq(0,1,by = 0.01),
                         confidence.intervals = 0.95)
Nomogram=decision_curve(Type ~ Nomogram, data=rt, 
                        family = binomial(link ='logit'),
                        thresholds= seq(0,1,by = 0.01),
                        confidence.intervals = 0.95)

pdf(file="DCA.pdf", width=5.5, height=5)
plot_decision_curve(list(riskScore, Nomogram),
                    curve.names=c("riskScore", "Nomogram"),
                    xlab="Threshold probability",
                    cost.benefit.axis=T,
                    confidence.intervals=FALSE,
                    standardize=FALSE)
dev.off()

pdf(file="clinical_impact.pdf", width=6, height=6)
plot_clinical_impact(riskScore,
                     confidence.intervals=T,
                     col = c("red", "blue"))
dev.off()

#Figure 5E
library(readxl)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
table <- read_excel("F:\\8.小论文\\3.FTO\\生信分析 - 副本\\19. High-low risk\\m6A评分差异分析\\riskGroupcore.xlsx")
table$ID <- as.factor(table$ID) # Transform values to factor
table$Type <- as.factor(table$Type)

cairo_ps("riskGroupm6AScore.eps")
my_comparison <- list(c("Low", "High"))
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

#Figure 5F
library(glmnet)
library(pROC)
library(ggsci)

riskFile="nomoRisk.txt"

rt=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
y=rt[,"Type"]

bioCol=pal_simpsons(palette=c("springfield"), alpha=1)(length(2:ncol(rt)))
aucText=c()
k=0
for(x in colnames(rt)[ncol(rt):ncol(rt)]){
  k=k+1

  roc1=roc(y, as.numeric(rt[,x]))
  if(k==1){
    pdf(file="ROC.pdf", width=5.5, height=4.8)
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="")
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }else{
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", add=TRUE)
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }
}

legend("bottomright", aucText, lwd=2, bty="n", col=bioCol, cex=0.7)
dev.off()

#Figure 5G-H
library(limma)               
expFile="geneMatrix.txt"    
conFile="s1.txt"             
treatFile="s2.txt"          


rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

s1=read.table(conFile, header=F, sep="\t", check.names=F)
sampleName1=as.vector(s1[,1])
conData=data[,sampleName1]

s2=read.table(treatFile, header=F, sep="\t", check.names=F)
sampleName2=as.vector(s2[,1])
treatData=data[,sampleName2]

data=cbind(conData, treatData)
data=normalizeBetweenArrays(data)

um=ncol(conData)
treatNum=ncol(treatData)
Type=c(rep("con",conNum),rep("treat",treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData, file="normalize.txt", sep="\t", quote=F, col.names=F)

library(glmnet)
expFile="test.normalize.txt"
geneFile="modelGene.txt"

rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

RT=read.table(geneFile, header=T, sep="\t", check.names=F)
sameGene=intersect(as.vector(geneRT[,1]), row.names(rt))
data=rt[sameGene,]

#???=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
rt=as.data.frame(data)
rt$Type=ifelse(group=="con", 0, 1)

#???glm(Type ~ ., family="binomial", data=rt)
pred=predict(fit, type="response")
outTab=cbind(rt[,c("Type", sameGene)], riskScore=pred)
outTab=cbind(id=row.names(outTab), outTab)
write.table(outTab, file="risk.test.txt", sep="\t", quote=F, row.names=F)

library(glmnet)
library(pROC)
library(ggsci)

riskFile="nomoRisk.txt"

rt=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
y=rt[,"Type"]

bioCol=pal_simpsons(palette=c("springfield"), alpha=1)(length(2:ncol(rt)))
aucText=c()
k=0
for(x in colnames(rt)[ncol(rt):ncol(rt)]){
  k=k+1
  
  roc1=roc(y, as.numeric(rt[,x]))
  if(k==1){
    pdf(file="ROC.pdf", width=5.5, height=4.8)
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="")
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }else{
    plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", add=TRUE)
    aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
  }
}

legend("bottomright", aucText, lwd=2, bty="n", col=bioCol, cex=0.7)
dev.off()