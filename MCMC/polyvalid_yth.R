# polyvalid_yth.R
# Matthew J. Beattie
# November 2018

# This script investigates the model output from MCMC on the validation data

# Load dataset
setwd("~/OUThesis/MCMC2")
load("RF_data_12to14YTH.Rda")
newdf <- df12to14YTH
rm(df12to14YTH)

# Create a flag for oxycodone products other than OxyContin
newdf$onlyoxycod <- ifelse((newdf$oxyflag==0 & newdf$OXYCODP2==1), 1, 0)

# Create a flag for hallucinogens other than LSD, PCP, ecstasy.
newdf$otherhal <- ifelse((newdf$halflag==1 & newdf$ecsflag==0 & newdf$lsdflag==0 &
                            newdf$pcpflag==0 & newdf$PSILCY2==0 & newdf$MESC2==0), 1, 0)

# Create a flag for stimulants other than methamphetamine
newdf$otherstm <- ifelse((newdf$cpnstmfg==1 & newdf$cpnmthfg==0), 1, 0)

# Create a flag for pain killers other than morphine
newdf$otherpain <- ifelse((newdf$othanl==1 & newdf$MORPHIN2==0), 1, 0)

# Create a variable for the number of heroin risk drugs used
newdf$polyabuse <- newdf$cocflag + newdf$crkflag + newdf$ecsflag + newdf$onlyoxycod +
                   newdf$oxyflag

# Create a needle use flag
newdf$needleuse <- ifelse(newdf$cocneedl==1 | newdf$mthneedl==1 | newdf$otdgnedl==1,1,0)

# Set perceived risk of heroin use to 2 if regular ok, 1 if try ok, else 0
newdf$herrisk <- ifelse(newdf$grskhreg==0,2,
                        ifelse(newdf$grskhtry==0,1,0))

# Reduce dataset to drug flag variables
smalldf<- subset(newdf, select = c(polyabuse, needleuse, txilalev, rdifher, herflag, herrisk))

# Load parameter values from MCMC analysis
beta0 <- -8.0889
alpha1 <- 1.6043
alpha2 <- 1.0587
alpha4 <- 1.2091
alpha5 <- 0.6689
alpha6 <- 1.4954

# Calculate probability of heroin use for each respondent in validation set
# lmu is the linear combination of parameters and factor values
smalldf$lmu <- beta0 + alpha1*smalldf$polyabuse + alpha2*smalldf$needleuse + alpha4*smalldf$txilalev + 
  alpha5*smalldf$herrisk + alpha6*smalldf$rdifher
# mu represents the probability of heroin use
smalldf$mu <- exp(smalldf$lmu)/(1+exp(smalldf$lmu))

# Extract test set data.  This is the opposite of the random sample used in the model creation.
frac <- 0.70
smp_size <- floor(frac * nrow(smalldf))
set.seed(456)
test_ind <- sample(seq_len(nrow(smalldf)), size = smp_size)
valdf <- smalldf[-test_ind, ]


library(ggplot2)
g = ggplot(data=valdf, aes(x = mu, y=herflag)) + ylab("Used Heroin") + xlab("Probability")
g = g + geom_point(position = position_jitter(height = 0.015, width =0)) +
    ggtitle("Actual Use vs Probability: Youth 2012-2014") + theme_bw()
g

# Calculate AUROC
N = 10000
thr = 0
auc = 0
sens <- vector()
spec <- vector()
tpr <- vector()
fpr <- vector()
for (i in 1:N) {
  thr = thr + 1/N
  tp = sum(valdf$mu > thr & valdf$herflag)
  fn = sum(valdf$mu < thr & valdf$herflag)
  tn = sum(valdf$mu < thr & valdf$herflag==0)
  fp = sum(valdf$mu > thr & valdf$herflag==0)
  sens[i] = tp/(tp + fn)
  spec[i] = tn/(tn + fp)
  tpr[i] = sens[i]
  fpr[i] = 1 - spec[i]
}

# Calculate area of ROC trapezoid (http://stats.stackexchange.com/a/146174/46761)
auc = 0
for (i in 2:N) {
  dfpr = -(fpr[i]-fpr[i-1])
  auc = auc + (tpr[i]+tpr[i-1])/2 * dfpr
}
auclab = paste("AUROC = ", round(auc,4))

df <- data.frame(fpr,tpr)
names(df) <- c("fpr","tpr")

# Plot the ROC for the validation set
g <- ggplot(data=df, aes(x=fpr,y=tpr)) + ylab("True Positive Rate") + xlab("False Positive Rate") +
     ggtitle("ROC for Revised Youth Model 2012-2014") + theme_bw()
g <- g + geom_line(color='red', size=1) + geom_abline(intercept=0,slope=1,linetype="longdash") +
     geom_text(x=.25, y=.75, label=auclab)
g


# Plot the distribution of mu for the validation set
usemean <- paste("Mean =", round(mean(valdf[valdf$herflag==1,]$mu),4))
usemed <- paste("Median =", round(median(valdf[valdf$herflag==1,]$mu),4))
nonusemean <- paste("Mean =", round(mean(valdf[valdf$herflag==0,]$mu),4))
nonusemed <- paste("Median =", round(median(valdf[valdf$herflag==0,]$mu),4))


df2 <- valdf[valdf$herflag==1,]
df3 <- valdf[valdf$herflag==0,]
g2 <- ggplot(data=df2, aes(x=mu)) + geom_density(color='red') + theme_bw() + 
      ggtitle("Distribution of mu for Heroin Users, Youth 2012-2014") +
      theme(plot.title=element_text(size=10)) + geom_text(x=.5, y=2, label=usemean) +
      geom_text(x=.5, y=1.75, label=usemed)
g2

g3 <- ggplot(data=df3, aes(x=mu)) + geom_density(color='blue') + theme_bw() +
      ggtitle("Distribution of mu for Heroin Non-Users, Youth 2012-2014") +
      theme(plot.title=element_text(size=10)) + geom_text(x=.5, y=2000, label=nonusemean) +
      geom_text(x=.5, y=1700, label=nonusemed)
g3


