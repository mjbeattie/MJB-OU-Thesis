# NSDUH Decision Tree Analysis:  Prediction
# Matthew J. Beattie
# September 2018

# This script conducts decision tree analysis of the NSDUH dataset
library(party)          # Includes cforest, the conditional tree random forest algorithm
library(ROCR)           # for graphics
library(caret)          # Statisical library


#********** Version 1 Analysis:  Model heroin use on ALL variables except ones dependent on positive heroin use
# Load dataset
load("RF_data_03to05ADLT.Rda")
fit.df03to05ADLT <- readRDS("03to05ADLT_cf3.rds")

# Create test and training sets
# 20% of the sample size
smp_size <- floor(0.20 * nrow(df03to05ADLT))
set.seed(123)
train_ind <- sample(seq_len(nrow(df03to05ADLT)), size = smp_size)
df03to05ADLT.train <- df03to05ADLT[train_ind, ]
df03to05ADLT.test <- df03to05ADLT[-train_ind, ]

rm(df03to05ADLT)

# Take random sample of test set for predict
smp_size <- floor(0.25 * nrow(df03to05ADLT.test))
set.seed(456)
train_ind <- sample(seq_len(nrow(df03to05ADLT.test)), size = smp_size)
ADLT.test <- df03to05ADLT.test[train_ind, ]

#******************** Run Prediction ********************

# Examine performance against train dataset
pred <- as.data.frame(Predict(fit.df03to05ADLT, newdata=df03to05ADLT.train))
colnames(pred)[1] <- "fitted.values"
pred$prediction <- ifelse(pred$fitted.values>=0.2,1,0)
confusionMatrix(as.factor(pred$prediction), as.factor(df03to05ADLT.train$herflag))

# Examine performance against test dataset
pred <- as.data.frame(Predict(fit.df03to05ADLT, newdata=ADLT.test))
colnames(pred)[1] <- "fitted.values"
pred$prediction <- ifelse(pred$fitted.values>=0.2,1,0)
confusionMatrix(as.factor(pred$prediction), as.factor(ADLT.test$herflag))

pred$actual <- ADLT.test$herflag
saveRDS(pred, "cf3_test_pred.rds")


# Examine AUROC
#fit.df03to05ADLT <- readRDS("df03to05ADLT_cf1.rds")
#pred <- as.data.frame(Predict(fit.df03to05ADLT, df03to05ADLT.test))
#rf.pr = pred[,1]
#rf.pred = prediction(rf.pr, df03to05ADLT.test$herflag)
#rf.perf = performance(rf.pred,"tpr","fpr")
#plot(rf.perf,main="ROC for ADLT 03-05:  Run 1",col=2,lwd=2)
#abline(a=0,b=1,lwd=2,lty=2,col="gray")
#perf <- performance(rf.pred, measure = "auc")
#print(c("AUC: ", perf@y.values))
