# NSDUH Random Conditional Inference Forest Tree Analysis
# Matthew J. Beattie
# September 2018

# This script conducts decision tree analysis of the NSDUH dataset
library(party)          # Includes cforest, the conditional tree random forest algorithm
library(ROCR)           # for graphics
library(caret)          # Statisical library


#********** Version 1 Analysis:  Model heroin use on ALL variables except ones dependent on positive heroin use
# Load dataset
load("RF_data_12to14YTH.Rda")

users <- df12to14YTH[df12to14YTH$herflag==1,]
nonusrs <- df12to14YTH[df12to14YTH$herflag!=1,]

# Take random sample of 19000 from nonusers
frac <- min(1, 19000/nrow(nonusrs))
smp_size <- floor(frac * nrow(nonusrs))
set.seed(456)
train_ind <- sample(seq_len(nrow(nonusrs)), size = smp_size)
nonusrs <- nonusrs[train_ind, ]

# Combine users and nonusers dataframes
dftrain <- rbind(users, nonusrs)

# Clear memory
rm(df12to14YTH, users, nonusrs)

#******************** Run Random Forest ********************
# Run Conditional Inference Random Forest routine
starttime <- Sys.time()
fit.df12to14YTH <- cforest(herflag ~ ., data = dftrain, 
                            controls = cforest_unbiased(mtry=20, ntree=300))
endtime <- Sys.time()
cat("Random Forest running time was: ", endtime - starttime)

summary(fit.df12to14YTH)

varimps <- varimp(fit.df12to14YTH)                     # Calculate variable importance scores
saveRDS(fit.df12to14YTH, "12to14YTH_cf18.rds")         # Save fit dataset
write.csv(varimps, file ="12to14YTH_cf18_imp.csv", 
          fileEncoding="UTF-8")                         # Save variable importance scores to CSV file

# Examine performance against train dataset
pred <- as.data.frame(Predict(fit.df12to14YTH, newdata=dftrain))
colnames(pred)[1] <- "fitted.values"
pred$prediction <- ifelse(pred$fitted.values>=0.2,1,0)
confusionMatrix(as.factor(pred$prediction), as.factor(dftrain$herflag))

# Save predictions for analysis
pred$actual <- dftrain$herflag
saveRDS(pred, "cf18_test_pred.rds")

