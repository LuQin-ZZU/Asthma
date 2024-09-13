#DEGs between asthma and control
logFoldChange=0.3
adjustP=0.05

library(limma)

rt=read.table("geneMatrix.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

modType=c(rep("normal",20),rep("asthma",88),rep("normal",27),rep("asthma",128))
design <- model.matrix(~0+factor(modType))
colnames(design) <- c("treat","con")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiff=cbind(geneNames=row.names(allDiff), allDiff)
row.names(allDiff)=NULL
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F,row.names=F)

diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & adj.P.Val < adjustP )), ]
write.table(diffSig,file="diff.xls",sep="\t",quote=F,row.names=F)
diffUp <- allDiff[with(allDiff, (logFC>logFoldChange & adj.P.Val < adjustP )), ]
write.table(diffUp,file="up.xls",sep="\t",quote=F,row.names=F)
diffDown <- allDiff[with(allDiff, (logFC<(-logFoldChange) & adj.P.Val < adjustP )), ]
write.table(diffDown,file="down.xls",sep="\t",quote=F,row.names=F)

hmExp=rt[as.vector(diffSig[,1]),]
diffExp=rbind(id=colnames(hmExp),hmExp)
write.table(diffExp,file="diffExp.txt",sep="\t",quote=F,col.names=F)

#Figure 6A
library(GEOquery)
library(limma)
library(ggplot2)
library(ggrepel)
library(ggthemes)

x=read.table("limmaTab.xls",sep="\t",header=T,check.names=F)
rownames(x)=x[,1]

logFCcut <- 0.3
pvalCut <- 0.05
sum(x$adj.P.Val < pvalCut)

deg.adj <- x[x$adj.P.Val <pvalCut & (x$logFC > logFCcut | x$logFC < -logFCcut),]


x$change <- ifelse((x$adj.P.Val < pvalCut & x$logFC > logFCcut), "lightcoral", ifelse((x$P.Value < pvalCut & x$logFC < -logFCcut), "steelblue1","grey30"))
size <- ifelse((x$adj.P.Val < pvalCut & abs(x$logFC) > logFCcut), 2, 1)

xmin <- (range(x$logFC)[1]-(range(x$logFC)[1]+1))
xmax <- (range(x$logFC)[1]+(1-range(x$logFC)[1]))
ymin <- 0
ymax <- 20

p1 <- ggplot(data=x, aes(x=logFC, y=-log10(adj.P.Val), label = geneNames)) +
  geom_point(alpha = 0.6, size=size, colour=x$change, shape = 16) +
  scale_color_manual(values = c("lightgrey", "navy", "red")) + 
  
  labs(x=bquote(~Log[2]~"(fold change)"), y=bquote(~-Log[10]~italic("(adj.P-value)")), title="Control vs. Asthma") + 
  ylim(c(ymin,ymax)) + 
  scale_x_continuous(
    breaks = c(-2, -1, 0, 1, 2),
    labels = c(-2, -1, 0, 1, 2),
    limits = c(-2.5, 2.5)
  ) +
  
  geom_vline(xintercept = logFCcut, color="grey40", linetype="longdash", size=0.5) +
  geom_vline(xintercept = -logFCcut, color="grey40", linetype="longdash", size=0.5) +
  geom_hline(yintercept = -log10(pvalCut), color="grey40", linetype="longdash", size=0.5) +
  
  guides(colour = guide_legend(override.aes = list(shape=16)))+
  
  theme_bw(base_size = 12, base_family = "Times") +
  theme(legend.position="top",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Times", size=8),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="plain", color="black", size=13),
        axis.text.y = element_text(face="plain",  color="black", size=13),
        axis.title.x = element_text(face="bold", color="black", size=15),
        axis.title.y = element_text(face="bold",color="black", size=15))
p1

p1 + geom_text_repel(aes(x = logFC, y = -log10(adj.P.Val), 
                         label = ifelse(logFC > 3, rownames(x),"")),
                     colour="darkred", size = 5, box.padding = unit(0.35, "lines"), 
                     point.padding = unit(0.3, "lines"))

selectedGeneID <- read.table("highlight_gene.txt",header = T)
selectgenes <- x[row.names(x) %in% selectedGeneID$Symbol,]
head(selectgenes)

p2 <- p1 + 
  
  geom_point(data = selectgenes, alpha = 1, size = 2.1, shape = 1, 
             stroke = 1, 
             color = "black") +
  
  geom_text_repel(data = selectgenes, 
                  colour="black", size = 5)
p2

ggsave("volcano_simple.pdf", width=6,height=5)
ggsave("volcano_simple.tiff", width=6,height=5)
dev.off()