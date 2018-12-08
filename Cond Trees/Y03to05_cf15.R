# NSDUH Random Conditional Inference Forest Tree Analysis
# Matthew J. Beattie
# September 2018

# This script conducts decision tree analysis of the NSDUH dataset
library(party)          # Includes cforest, the conditional tree random forest algorithm
library(ROCR)           # for graphics
library(caret)          # Statisical library


#********** Version 1 Analysis:  Model heroin use on ALL variables except ones dependent on positive heroin use
# Load dataset
load("RF_data_03to05YTH.Rda")

users <- df03to05YTH[df03to05YTH$herflag==1,]
nonusrs <- df03to05YTH[df03to05YTH$herflag!=1,]

# Take random sample of 19000 from nonusers
frac <- min(1, 19000/nrow(nonusrs))
smp_size <- floor(frac * nrow(nonusrs))
set.seed(456)
train_ind <- sample(seq_len(nrow(nonusrs)), size = smp_size)
nonusrs <- nonusrs[train_ind, ]

# Combine users and nonusers dataframes
dftrain <- rbind(users, nonusrs)

# Clear memory
rm(df03to05YTH, users, nonusrs)

#******************** Run Random Forest ********************
# Run Conditional Inference Random Forest routine
starttime <- Sys.time()
fit.df03to05YTH <- cforest(herflag ~ ., data = dftrain, 
                            controls = cforest_unbiased(mtry=20, ntree=300))
endtime <- Sys.time()
cat("Random Forest running time was: ", endtime - starttime)

summary(fit.df03to05YTH)

varimps <- varimp(fit.df03to05YTH)                     # Calculate variable importance scores
saveRDS(fit.df03to05YTH, "03to05YTH_cf15.rds")         # Save fit dataset
write.csv(varimps, file ="03to05YTH_cf15_imp.csv", 
          fileEncoding="UTF-8")                         # Save variable importance scores to CSV file

# Examine performance against train dataset
pred <- as.data.frame(Predict(fit.df03to05YTH, newdata=dftrain))
colnames(pred)[1] <- "fitted.values"
pred$prediction <- ifelse(pred$fitted.values>=0.2,1,0)
confusionMatrix(as.factor(pred$prediction), as.factor(dftrain$herflag))

# Save predictions for analysis
pred$actual <- dftrain$herflag
saveRDS(pred, "cf15_test_pred.rds")

