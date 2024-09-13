#Figure 4A-C
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(xgboost)
library(pROC)

set.seed(123)
inputFile="diffGeneExp.txt"

data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))
data=as.data.frame(data)
data$Type=group

inTrain<-createDataPartition(y=data$Type, p=0.7, list=F)
train<-data[inTrain,]
test<-data[-inTrain,]

control=trainControl(method="repeatedcv", number=9, s100veaedictions=TRUE)
mod_rf = train(Type ~ ., data = train, method='rf', trControl = control)

mod_svm=train(Type ~., data = train, method = "svmRadial", prob.model=TRUE, trControl=control)

mod_xgb=train(Type ~., data = train, method = "xgbDART", trControl=control)

mod_glm=train(Type ~., data = train, method = "glm", family="binomial", trControl=control)

mod_gbm=train(Type ~., data = train, method = "gbm", trControl=control)

mod_lasso=train(Type ~., data = train, method = "glmnet", trControl=control)

p_fun=function(object, newdata){
  predict(object, newdata=newdata, type="prob")[,2]
}
yTest=ifelse(test$Type=="con", 0, 1)

explainer_rf=explain(mod_rf, label = "RF",
                     data = test, y = yTest,
                     predict_function = p_fun,
                     verbose = FALSE)
mp_rf=model_performance(explainer_rf)

explainer_svm=explain(mod_svm, label = "SVM",
                      data = test, y = yTest,
                      predict_function = p_fun,
                      verbose = FALSE)
mp_svm=model_performance(explainer_svm)

explainer_xgb=explain(mod_xgb, label = "XGB",
                      data = test, y = yTest,
                      predict_function = p_fun,
                      verbose = FALSE)
mp_xgb=model_performance(explainer_xgb)

explainer_glm=explain(mod_glm, label = "GLM",
                      data = test, y = yTest,
                      predict_function = p_fun,
                      verbose = FALSE)
mp_glm=model_performance(explainer_glm)

explainer_gbm=explain(mod_gbm, label = "GBM",
                      data = test, y = yTest,
                      predict_function = p_fun,
                      verbose = FALSE)
mp_gbm=model_performance(explainer_gbm)

explainer_lasso=explain(mod_lasso, label = "LASSO",
                        data = test, y = yTest,
                        predict_function = p_fun,
                        verbose = FALSE)
mp_lasso=model_performance(explainer_lasso)

pdf(file="residual.pdf", width=6, height=6)
p1 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm, mp_gbm, mp_lasso)
print(p1)
dev.off()

pdf(file="boxplot.pdf", width=6, height=6)
p2 <- plot(mp_rf, mp_svm, mp_xgb, mp_glm, mp_gbm, mp_lasso, geom = "boxplot")
print(p2)
dev.off()

pred1=predict(mod_rf, newdata=test, type="prob")
pred2=predict(mod_svm, newdata=test, type="prob")
pred3=predict(mod_xgb, newdata=test, type="prob")
pred4=predict(mod_glm, newdata=test, type="prob")
pred5=predict(mod_gbm, newdata=test, type="prob")
pred6=predict(mod_lasso, newdata=test, type="prob")

roc1=roc(yTest, as.numeric(pred1[,2]))
roc2=roc(yTest, as.numeric(pred2[,2]))
roc3=roc(yTest, as.numeric(pred3[,2]))
roc4=roc(yTest, as.numeric(pred4[,2]))
roc5=roc(yTest, as.numeric(pred5[,2]))
roc6=roc(yTest, as.numeric(pred6[,2]))

pdf(file="ROC.pdf", width=5, height=5)
plot(roc1, print.auc=F, legacy.axes=T, main="", col="red")
plot(roc2, print.auc=F, legacy.axes=T, main="", col="blue", add=T)
plot(roc3, print.auc=F, legacy.axes=T, main="", col="green", add=T)
plot(roc4, print.auc=F, legacy.axes=T, main="", col="yellow", add=T)
plot(roc5, print.auc=F, legacy.axes=T, main="", col="orange", add=T)
plot(roc6, print.auc=F, legacy.axes=T, main="", col="purple", add=T)

legend('bottomright',
       c(paste0('RF: ',sprintf("%.03f",roc1$auc)),
         paste0('SVM: ',sprintf("%.03f",roc2$auc)),
         paste0('XGB: ',sprintf("%.03f",roc3$auc)),
         paste0('GLM: ',sprintf("%.03f",roc4$auc)),
         paste0('GBM: ',sprintf("%.03f",roc5$auc)),
         paste0('LASSO: ',sprintf("%.03f",roc6$auc))),
       
       col=c("red","blue","green","yellow","orange","purple"), lwd=2, bty = 'n')
dev.off()

importance_rf<-variable_importance(
  explainer_rf,
  loss_function = loss_root_mean_square
)
importance_svm<-variable_importance(
  explainer_svm,
  loss_function = loss_root_mean_square
)
importance_glm<-variable_importance(
  explainer_glm,
  loss_function = loss_root_mean_square
)
importance_xgb<-variable_importance(
  explainer_xgb,
  loss_function = loss_root_mean_square
)
importance_gbm<-variable_importance(
  explainer_gbm,
  loss_function = loss_root_mean_square
)
importance_lasso<-variable_importance(
  explainer_lasso,
  loss_function = loss_root_mean_square
)

geneNum=6
write.table(importance_rf[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.RF.txt", sep="\t", quote=F, row.names=F)
write.table(importance_svm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.SVM.txt", sep="\t", quote=F, row.names=F)
write.table(importance_xgb[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.XGB.txt", sep="\t", quote=F, row.names=F)
write.table(importance_glm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.GLM.txt", sep="\t", quote=F, row.names=F)
write.table(importance_gbm[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.GBM.txt", sep="\t", quote=F, row.names=F)
write.table(importance_lasso[(ncol(data)-geneNum+2):(ncol(data)+1),], file="importanceGene.LASSO.txt", sep="\t", quote=F, row.names=F)

#RF is the optimal model
riskScores_raw <- predict(mod_rf, newdata = data, type = "prob")
riskScores <- riskScores_raw[,2]
cut_off <- median(riskScores)
breaks <- c(0, cut_off, max(riskScores, na.rm = TRUE))
labels <- c("Low", "High")
riskGroups <- cut(riskScores, breaks = breaks, labels = labels)
risk <- data.frame(RiskScores = riskScores, RiskGroups = riskGroups)
rownames(risk) <- rownames(riskScores_raw)
write.table(risk, file = "risk_data.txt", row.names = TRUE, col.names = TRUE, sep = "\t")

#Figure 4D-E
library(randomForest)
set.seed(123)

inputFile="diffGeneExp.txt"

data=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
data=t(data)
group=gsub("(.*)\\_(.*)", "\\2", row.names(data))

rf=randomForest(as.factor(group)~., data=data, ntree=500)
pdf(file="forest.pdf", width=6, height=6)
plot(rf, main="Random forest", lwd=2)
dev.off()

optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(as.factor(group)~., data=data, ntree=optionTrees)

importance=importance(x=rf2)

pdf(file="geneImportance.pdf", width=6.2, height=5.8)
varImpPlot(rf2, main="")
dev.off()

rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes=names(rfGenes[1:4])
write.table(rfGenes, file="rfGenes.txt", sep="\t", quote=F, col.names=F, row.names=F)

sigExp=t(data[,rfGenes])
sigExpOut=rbind(ID=colnames(sigExp),sigExp)
write.table(sigExpOut, file="rfGeneExp.txt", sep="\t", quote=F, col.names=F)