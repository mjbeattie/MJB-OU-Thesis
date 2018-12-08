# NSDUH Random Conditional Inference Forest Tree Analysis
# Matthew J. Beattie
# September 2018

# This script conducts decision tree analysis of the NSDUH dataset
library(party)          # Includes cforest, the conditional tree random forest algorithm
library(ROCR)           # for graphics
library(caret)          # Statisical library


#********** Version 1 Analysis:  Model heroin use on ALL variables except ones dependent on positive heroin use
# Load dataset
load("RF_data_09to11ADLT.Rda")
#load("RF_data_09to11YTH.Rda")

# Create test and training sets
# 25% of the sample size
smp_size <- floor(0.35 * nrow(df09to11ADLT))
set.seed(123)
train_ind <- sample(seq_len(nrow(df09to11ADLT)), size = smp_size)
df09to11ADLT.train <- df09to11ADLT[train_ind, ]
df09to11ADLT.test <- df09to11ADLT[-train_ind, ]

rm(df09to11ADLT)                                # Clear memory

# Take random sample of train set for predict
frac <- min(1, 20000/nrow(df09to11ADLT.train))
smp_size <- floor(frac * nrow(df09to11ADLT.train))
set.seed(456)
train_ind <- sample(seq_len(nrow(df09to11ADLT.train)), size = smp_size)
ADLT.train <- df09to11ADLT.train[train_ind, ]

# Take random sample of test set for predict
frac <- min(1, 20000/nrow(df09to11ADLT.test))
smp_size <- floor(frac * nrow(df09to11ADLT.test))
set.seed(456)
train_ind <- sample(seq_len(nrow(df09to11ADLT.test)), size = smp_size)
ADLT.test <- df09to11ADLT.test[train_ind, ]

rm(df09to11ADLT.test)                            # Clear memory


#******************** Run Random Forest ********************
# Run Conditional Inference Random Forest routine
starttime <- Sys.time()
fit.df09to11ADLT <- cforest(herflag ~ ., data = df09to11ADLT.train, 
                            controls = cforest_unbiased(mtry=10, ntree=300))
endtime <- Sys.time()
cat("Random Forest running time was: ", endtime - starttime)

summary(fit.df09to11ADLT)

varimps <- varimp(fit.df09to11ADLT)                     # Calculate variable importance scores
saveRDS(fit.df09to11ADLT, "09to11ADLT_cf10.rds")         # Save fit dataset
write.csv(varimps, file ="09to11ADLT_cf10_imp.csv", 
          fileEncoding="UTF-8")                         # Save variable importance scores to CSV file

# Visualize tree
#pt <- prettytree(fit.df09to11ADLT@ensemble[[1]], names(fit.df09to11ADLT@data@get("input"))) 
#nt <- new("BinaryTree") 
#nt@tree <- pt 
#nt@data <- fit.df09to11ADLT@data 
#nt@responses <- fit.df09to11ADLT@responses 
#plot(nt, type="simple")

# Examine performance against train dataset
pred <- as.data.frame(Predict(fit.df09to11ADLT, newdata=ADLT.train))
colnames(pred)[1] <- "fitted.values"
pred$prediction <- ifelse(pred$fitted.values>=0.2,1,0)
confusionMatrix(as.factor(pred$prediction), as.factor(ADLT.train$herflag))

# Examine performance against test dataset
pred <- as.data.frame(Predict(fit.df09to11ADLT, newdata=ADLT.test))
colnames(pred)[1] <- "fitted.values"
pred$prediction <- ifelse(pred$fitted.values>=0.2,1,0)
confusionMatrix(as.factor(pred$prediction), as.factor(ADLT.test$herflag))

# Save test predictions for analysis
pred$actual <- ADLT.test$herflag
saveRDS(pred, "cf10_test_pred.rds")

# Examine AUROC
#rf.pred = prediction(pred$fitted.values, pred$actual)
#rf.perf = performance(rf.pred,"tpr","fpr")
#plot(rf.perf,main="ROC for ADLT 03-05:  Run 1",col=2,lwd=2)
#abline(a=0,b=1,lwd=2,lty=2,col="gray")
#perf <- performance(rf.pred, measure = "auc")
#print(c("AUC: ", perf@y.values))