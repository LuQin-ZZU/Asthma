#Figure 7E
library(readxl)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
table <- read_excel("FTO.xlsx")
table$DataSet <- as.factor(table$DataSet) # Transform values to factor
table$Type <- as.factor(table$Type)

cairo_ps("FTO.eps", width = 12, height = 10)
ggboxplot(table, x="DataSet", y="Level", 
          color = "Type", add = "jitter")+
  ylim(4.5, 15)+
  theme(legend.position='right')+
  xlab("")+
  ylab("FTO expression")+
  border(size=0.2)
dev.off()

#Figure 7F-G
library(pheatmap)

bioRiskPlot=function(inputFile=null, project=null){
  rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
  rt=rt[order(rt$FTO),] 
  
  riskClass=rt[,"FTOType"]
  lowLength=length(riskClass[riskClass=="Low"])
  highLength=length(riskClass[riskClass=="High"])
  lowMax=max(rt$FTO[riskClass=="Low"])
  line=rt[,"FTO"]
  line[line>10]=10
  pdf(file=paste0(project, ".FTOExpression.pdf"), width=7, height=4)
  plot(line, type="p", pch=20,
       xlab="Patients (increasing FTO expression)",
       ylab="FTO expression",
       col=c(rep("blue",lowLength),rep("red",highLength)) )
  abline(h=lowMax,v=lowLength,lty=2)
  legend("topleft", c("High","Low"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
  dev.off()
  
  color=as.vector(rt$ACQControl)
  color[color==0]="red"
  color[color==1]="blue"
  pdf(file=paste0(project, ".ACQControl.pdf"), width=7, height=4)
  plot(rt$FEV1Reversibility, pch=19,
       xlab="Patients (increasing FTO expression)",
       ylab="FEV1 reversibility (%)",
       col=color)
  legend("topleft", c("Uncontrolled","Controlled"),bty="n",pch=19,col=c("red","blue"),cex=1.2)
  abline(v=lowLength,lty=2)
  dev.off()
}


bioRiskPlot(inputFile="risk.txt", project="all")

#Figure 7H
library(survival)
library(survminer)

bioSurvival=function(inputFile=null, outFile=null){
  
  rt=read.table(inputFile, header=T, sep="\t")
  #?Ƚϸߵͷ????????????죬?õ???????pֵ
  diff=survdiff(Surv(FEV1Reversibility, ACQControl) ~Type,data = rt)
  pValue=1-pchisq(diff$chisq,df=1)
  if(pValue<0.001){
    pValue="p<0.001"
  }else{
    pValue=paste0("p=",sprintf("%.03f",pValue))
  }
  fit <- survfit(Surv(FEV1Reversibility, ACQControl) ~ Type, data = rt)
  
  surPlot=ggsurvplot(fit, 
                     data=rt,
                     conf.int=T,
                     pval=pValue,
                     pval.size=6,
                     legend.title="Risk",
                     legend.labs=c("High", "Low"),
                     xlab="FEV1 reversibility (%)",
                     ylab="ACQ control",
                     break.time.by = 4,
                     palette=c("red", "blue"),
                     risk.table=TRUE,
                     risk.table.title="",
                     risk.table.col = "strata",
                     risk.table.height=.25)
  
  #????ͼ??
  pdf(file=outFile, width = 6.5, height =5.5, onefile = FALSE)
  print(surPlot)
  dev.off()
}

bioSurvival(inputFile="risk.txt", outFile="surv.pdf")

#Figure 7I
library(survival)

bioForest=function(coxFile=null, forestFile=null, forestCol=null){
  rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt$"HR")
  hrLow  <- sprintf("%.3f",rt$"HR.95L")
  hrHigh <- sprintf("%.3f",rt$"HR.95H")
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  pdf(file=forestFile, width=6.6, height=4.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  
  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  
  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
  dev.off()
}

indep=function(riskFile=null, cliFile=null, project=null){
  risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
  cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
  
  sameSample=intersect(row.names(cli),row.names(risk))
  risk=risk[sameSample,]
  cli=cli[sameSample,]
  rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
  
  uniCoxFile=paste0(project,".uniCox.txt")
  uniCoxPdf=paste0(project,".uniCox.pdf")
  uniTab=data.frame()
  for(i in colnames(rt[,3:ncol(rt)])){
    cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
    coxSummary = summary(cox)
    uniTab=rbind(uniTab,
                 cbind(id=i,
                       HR=coxSummary$conf.int[,"exp(coef)"],
                       HR.95L=coxSummary$conf.int[,"lower .95"],
                       HR.95H=coxSummary$conf.int[,"upper .95"],
                       pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
    )
  }
  write.table(uniTab,file=uniCoxFile,sep="\t",row.names=F,quote=F)
  bioForest(coxFile=uniCoxFile, forestFile=uniCoxPdf, forestCol="green")
  
  multiCoxFile=paste0(project,".multiCox.txt")
  multiCoxPdf=paste0(project,".multiCox.pdf")
  uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<1,]
  rt1=rt[,c("futime","fustat",as.vector(uniTab[,"id"]))]
  multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
  multiCoxSum=summary(multiCox)
  multiTab=data.frame()
  multiTab=cbind(
    HR=multiCoxSum$conf.int[,"exp(coef)"],
    HR.95L=multiCoxSum$conf.int[,"lower .95"],
    HR.95H=multiCoxSum$conf.int[,"upper .95"],
    pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
  multiTab=cbind(id=row.names(multiTab),multiTab)
  write.table(multiTab, file=multiCoxFile, sep="\t", row.names=F, quote=F)
  bioForest(coxFile=multiCoxFile, forestFile=multiCoxPdf, forestCol="red")
}

indep(riskFile="risk.txt", cliFile="clinical.txt", project="all"