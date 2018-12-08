# NSDUH Random Conditional Inference Forest Tree Analysis
# Matthew J. Beattie
# September 2018

# This script conducts decision tree analysis of the NSDUH dataset
library(party)          # Includes cforest, the conditional tree random forest algorithm
library(ROCR)           # for graphics
library(caret)          # Statisical library


#********** Version 1 Analysis:  Model heroin use on ALL variables except ones dependent on positive heroin use
# Load dataset
load("RF_data_03to05ADLT.Rda")
#load("RF_data_03to05YTH.Rda")

# Create test and training sets
# 25% of the sample size
smp_size <- floor(0.35 * nrow(df03to05ADLT))
set.seed(123)
train_ind <- sample(seq_len(nrow(df03to05ADLT)), size = smp_size)
df03to05ADLT.train <- df03to05ADLT[train_ind, ]
df03to05ADLT.test <- df03to05ADLT[-train_ind, ]

rm(df03to05ADLT)                                # Clear memory

# Take random sample of train set for predict
frac <- min(1, 20000/nrow(df03to05ADLT.train))
smp_size <- floor(frac * nrow(df03to05ADLT.train))
set.seed(456)
train_ind <- sample(seq_len(nrow(df03to05ADLT.train)), size = smp_size)
ADLT.train <- df03to05ADLT.train[train_ind, ]

# Take random sample of test set for predict
frac <- min(1, 20000/nrow(df03to05ADLT.test))
smp_size <- floor(frac * nrow(df03to05ADLT.test))
set.seed(456)
train_ind <- sample(seq_len(nrow(df03to05ADLT.test)), size = smp_size)
ADLT.test <- df03to05ADLT.test[train_ind, ]

rm(df03to05ADLT.test)                            # Clear memory


#******************** Run Random Forest ********************
# Run Conditional Inference Random Forest routine
starttime <- Sys.time()
fit.df03to05ADLT <- cforest(herflag ~ ., data = df03to05ADLT.train, 
                            controls = cforest_unbiased(mtry=10, ntree=200))
endtime <- Sys.time()
cat("Random Forest running time was: ", endtime - starttime)

summary(fit.df03to05ADLT)

varimps <- varimp(fit.df03to05ADLT)                     # Calculate variable importance scores
saveRDS(fit.df03to05ADLT, "03to05ADLT_cf7.rds")         # Save fit dataset
write.csv(varimps, file ="03to05ADLT_cf7_imp.csv", 
          fileEncoding="UTF-8")                         # Save variable importance scores to CSV file

# Visualize tree
#pt <- prettytree(fit.df03to05ADLT@ensemble[[1]], names(fit.df03to05ADLT@data@get("input"))) 
#nt <- new("BinaryTree") 
#nt@tree <- pt 
#nt@data <- fit.df03to05ADLT@data 
#nt@responses <- fit.df03to05ADLT@responses 
#plot(nt, type="simple")

# Examine performance against train dataset
pred <- as.data.frame(Predict(fit.df03to05ADLT, newdata=ADLT.train))
colnames(pred)[1] <- "fitted.values"
pred$prediction <- ifelse(pred$fitted.values>=0.2,1,0)
confusionMatrix(as.factor(pred$prediction), as.factor(ADLT.train$herflag))

# Examine performance against test dataset
pred <- as.data.frame(Predict(fit.df03to05ADLT, newdata=ADLT.test))
colnames(pred)[1] <- "fitted.values"
pred$prediction <- ifelse(pred$fitted.values>=0.2,1,0)
confusionMatrix(as.factor(pred$prediction), as.factor(ADLT.test$herflag))

# Save test predictions for analysis
pred$actual <- ADLT.test$herflag
saveRDS(pred, "cf7_test_pred.rds")

# Examine AUROC
#rf.pred = prediction(pred$fitted.values, pred$actual)
#rf.perf = performance(rf.pred,"tpr","fpr")
#plot(rf.perf,main="ROC for ADLT 03-05:  Run 1",col=2,lwd=2)
#abline(a=0,b=1,lwd=2,lty=2,col="gray")
#perf <- performance(rf.pred, measure = "auc")
#print(c("AUC: ", perf@y.values))