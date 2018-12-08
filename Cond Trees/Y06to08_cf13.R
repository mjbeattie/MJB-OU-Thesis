# NSDUH Random Conditional Inference Forest Tree Analysis
# Matthew J. Beattie
# September 2018

# This script conducts decision tree analysis of the NSDUH dataset
library(party)          # Includes cforest, the conditional tree random forest algorithm
library(ROCR)           # for graphics
library(caret)          # Statisical library


#********** Version 1 Analysis:  Model heroin use on ALL variables except ones dependent on positive heroin use
# Load dataset
#load("RF_data_06to08YTH.Rda")
load("RF_data_06to08YTH.Rda")

# Create test and training sets
# 25% of the sample size
smp_size <- floor(0.35 * nrow(df06to08YTH))
set.seed(123)
train_ind <- sample(seq_len(nrow(df06to08YTH)), size = smp_size)
df06to08YTH.train <- df06to08YTH[train_ind, ]
df06to08YTH.test <- df06to08YTH[-train_ind, ]

rm(df06to08YTH)                                # Clear memory

# Take random sample of train set for predict
frac <- min(1, 20000/nrow(df06to08YTH.train))
smp_size <- floor(frac * nrow(df06to08YTH.train))
set.seed(456)
train_ind <- sample(seq_len(nrow(df06to08YTH.train)), size = smp_size)
YTH.train <- df06to08YTH.train[train_ind, ]

# Take random sample of test set for predict
frac <- min(1, 20000/nrow(df06to08YTH.test))
smp_size <- floor(frac * nrow(df06to08YTH.test))
set.seed(456)
train_ind <- sample(seq_len(nrow(df06to08YTH.test)), size = smp_size)
YTH.test <- df06to08YTH.test[train_ind, ]

rm(df06to08YTH.test)                            # Clear memory


#******************** Run Random Forest ********************
# Run Conditional Inference Random Forest routine
starttime <- Sys.time()
fit.df06to08YTH <- cforest(herflag ~ ., data = df06to08YTH.train, 
                            controls = cforest_unbiased(mtry=10, ntree=300))
endtime <- Sys.time()
cat("Random Forest running time was: ", endtime - starttime)

summary(fit.df06to08YTH)

varimps <- varimp(fit.df06to08YTH)                     # Calculate variable importance scores
saveRDS(fit.df06to08YTH, "06to08YTH_cf13.rds")         # Save fit dataset
write.csv(varimps, file ="06to08YTH_cf13_imp.csv", 
          fileEncoding="UTF-8")                         # Save variable importance scores to CSV file

# Visualize tree
#pt <- prettytree(fit.df06to08YTH@ensemble[[1]], names(fit.df06to08YTH@data@get("input"))) 
#nt <- new("BinaryTree") 
#nt@tree <- pt 
#nt@data <- fit.df06to08YTH@data 
#nt@responses <- fit.df06to08YTH@responses 
#plot(nt, type="simple")

# Examine performance against train dataset
pred <- as.data.frame(Predict(fit.df06to08YTH, newdata=YTH.train))
colnames(pred)[1] <- "fitted.values"
pred$prediction <- ifelse(pred$fitted.values>=0.2,1,0)
confusionMatrix(as.factor(pred$prediction), as.factor(YTH.train$herflag))

# Examine performance against test dataset
pred <- as.data.frame(Predict(fit.df06to08YTH, newdata=YTH.test))
colnames(pred)[1] <- "fitted.values"
pred$prediction <- ifelse(pred$fitted.values>=0.2,1,0)
confusionMatrix(as.factor(pred$prediction), as.factor(YTH.test$herflag))

# Save test predictions for analysis
pred$actual <- YTH.test$herflag
saveRDS(pred, "cf13_test_pred.rds")

# Examine AUROC
#rf.pred = prediction(pred$fitted.values, pred$actual)
#rf.perf = performance(rf.pred,"tpr","fpr")
#plot(rf.perf,main="ROC for YTH 03-05:  Run 1",col=2,lwd=2)
#abline(a=0,b=1,lwd=2,lty=2,col="gray")
#perf <- performance(rf.pred, measure = "auc")
#print(c("AUC: ", perf@y.values))