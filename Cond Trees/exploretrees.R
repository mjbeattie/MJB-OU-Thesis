# NSDUH Decision Tree Analysis
# Matthew J. Beattie
# September 2018

# This script conducts decision tree analysis of the NSDUH dataset
library(party)          # Includes cforest, the conditional tree random forest algorithm
#library(ggplot2)        # for graphics
library(ROCR)           # for graphics
#library(rattle)      		# fancy tree plot
#library(rpart.plot)			# enhanced tree plots
#library(RColorBrewer)		# color selection for fancy tree plot
library(caret)          # Statisical library
#library(dplyr)          # Used for mutation of data frames


#********** Version 1 Analysis:  Model heroin use on ALL variables except ones dependent on positive heroin use
# Load dataset
fit.df03to05ADLT <- readRDS("03to05ADLT_cf1.rds")
load("RF_data_03to05ADLT.Rda")

# Create test and training sets
# 20% of the sample size
smp_size <- floor(0.20 * nrow(df03to05ADLT))
set.seed(123)
train_ind <- sample(seq_len(nrow(df03to05ADLT)), size = smp_size)
df03to05ADLT.train <- df03to05ADLT[train_ind, ]
df03to05ADLT.test <- df03to05ADLT[-train_ind, ]

rm(df03to05ADLT)                                        # Clear memory

# Split test set into smaller chunks for prediction
ADLT1.test <- df03to05ADLT.test[1:10000,]
ADLT2.test <- df03to05ADLT.test[10001:20000,]
ADLT3.test <- df03to05ADLT.test[20001:30000,]
ADLT4.test <- df03to05ADLT.test[30001:40000,]
ADLT5.test <- df03to05ADLT.test[40001:50000,]
ADLT6.test <- df03to05ADLT.test[50001:60000,]
ADLT7.test <- df03to05ADLT.test[60001:70000,]
ADLT8.test <- df03to05ADLT.test[70001:80000,]
ADLT9.test <- df03to05ADLT.test[80001:89249,]

#******************** Explore random forest output ********************
# Run Conditional Inference Random Forest routine
summary(fit.df03to05ADLT)

varimps <- varimp(fit.df03to05ADLT)                     # Calculate variable importance scores

# Visualize tree
#pt <- prettytree(fit.df03to05ADLT@ensemble[[1]], names(fit.df03to05ADLT@data@get("input"))) 
#nt <- new("BinaryTree") 
#nt@tree <- pt 
#nt@data <- fit.df03to05ADLT@data 
#nt@responses <- fit.df03to05ADLT@responses 
#plot(nt, type="simple")

# Examine performance against train dataset
pred <- as.data.frame(Predict(fit.df03to05ADLT, newdata=df03to05ADLT.train))
colnames(pred)[1] <- "fitted.values"
pred$prediction <- ifelse(pred$fitted.values>=0.2,1,0)
confusionMatrix(as.factor(pred$prediction), as.factor(df03to05ADLT.train$herflag))

# Examine performance against test dataset
pred1 <- as.data.frame(Predict(fit.df03to05ADLT, newdata=ADLT1.test))
pred2 <- as.data.frame(Predict(fit.df03to05ADLT, newdata=ADLT2.test))
pred3 <- as.data.frame(Predict(fit.df03to05ADLT, newdata=ADLT3.test))
pred4 <- as.data.frame(Predict(fit.df03to05ADLT, newdata=ADLT4.test))
pred5 <- as.data.frame(Predict(fit.df03to05ADLT, newdata=ADLT5.test))
pred6 <- as.data.frame(Predict(fit.df03to05ADLT, newdata=ADLT6.test))
pred7 <- as.data.frame(Predict(fit.df03to05ADLT, newdata=ADLT7.test))
pred8 <- as.data.frame(Predict(fit.df03to05ADLT, newdata=ADLT8.test))
pred9 <- as.data.frame(Predict(fit.df03to05ADLT, newdata=ADLT9.test))





colnames(pred)[1] <- "fitted.values"
pred$prediction <- ifelse(pred$fitted.values>=0.2,1,0)
confusionMatrix(as.factor(pred$prediction), as.factor(ADLT1.test$herflag))

# Examine AUROC
fit.df03to05ADLT <- readRDS("df03to05ADLT_cf1.rds")
pred <- as.data.frame(Predict(fit.df03to05ADLT, df03to05ADLT.test))
rf.pr = pred[,1]
rf.pred = prediction(rf.pr, df03to05ADLT.test$herflag)
rf.perf = performance(rf.pred,"tpr","fpr")
plot(rf.perf,main="ROC for ADLT 03-05:  Run 1",col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")
perf <- performance(rf.pred, measure = "auc")
print(c("AUC: ", perf@y.values))