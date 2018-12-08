# y4flags_50r2.R
# Matthew J. Beattie
# November 2018

# This script conducts MCMC analysis of a 50% sample of the NSDUH dataset
# Years:  2012-2014
# Setup:  
# 1) Pull out rows with no OxyContin data.
# 2) Model: P(heroin) = Bernoulli(mu), logit(mu) = lin(drug flags)
# 3) Priors set based upon assumptions of three categories:  Uncertain, mean=0, Var=0.5;
#    Probably positive, mean=0.6931, Var=0.5;  Probability Strong, mean=1.386, Var=0.5


#********** Setup dataset and samples ************
# Load dataset
start_time <- Sys.time()
setwd("~/MCMC2")
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

# Reduce dataset to drug flag variables
varlist<-c("benzos","cocflag","CODEINE2","cpnmthfg","otherstm","crkflag",
           "DILAUD2","ecsflag","otherhal","lsdflag","MESC2","METHDON2","MORPHIN2",
           "otherpain","oxyflag","pcpflag","PSILCY2","trqflag","onlyoxycod","herflag")
smalldf<-newdf[,c(varlist)]
rm(newdf)

# Remove data for observations with no oxycontin data
smalldf <- smalldf[smalldf$oxyflag != -9,]

# Take random sample
frac <- 0.70
smp_size <- floor(frac * nrow(smalldf))
set.seed(456)
test_ind <- sample(seq_len(nrow(smalldf)), size = smp_size)
smalldf <- smalldf[test_ind, ]

# Prepare data for JAGS by creating list
dataList = list(y = smalldf$herflag, d1 = smalldf$benzos,
                d2 = smalldf$cocflag, d3 = smalldf$CODEINE2, d4 = smalldf$cpnmthfg,
                d5 = smalldf$otherstm, d6 = smalldf$crkflag, d7 = smalldf$DILAUD2, d8 = smalldf$ecsflag,
                d9 = smalldf$otherhal, d10 = smalldf$lsdflag, d11 = smalldf$MESC2,
                d12 = smalldf$METHDON2, d13 = smalldf$MORPHIN2, d14 = smalldf$otherpain, d15 = smalldf$onlyoxycod,
                d16 = smalldf$oxyflag, d17 = smalldf$pcpflag, d18 = smalldf$PSILCY2, d19 = smalldf$trqflag,
                N = length(smalldf$herflag))



#************** RUN JAGS ROUTINES ***************
# Uses runjags to take advantage of multiple CPU cores
library(runjags)
source("DBDA2E-utilities.R")

#Define the model:
modelString = "
model
{
  for (i in 1:N) 
  {
    y[i]  ~ dbern(mu[i])
    logit(mu[i])  <- beta0 +
              omega1*d1[i] + omega2*d2[i] + omega3*d3[i] + omega4*d4[i] + omega5*d5[i] +
              omega6*d6[i] + omega7*d7[i] + omega8*d8[i] + omega9*d9[i] + omega10*d10[i] +
              omega11*d11[i] + omega12*d12[i] + omega13*d13[i] + omega14*d14[i] +
              omega15*d15[i] + omega16*d16[i] + omega17*d17[i] + omega18*d18[i] +
              omega19*d19[i]
   }
    beta0      ~ dnorm(0,0.5)
    omega1     ~ dnorm(0,0.5)
    omega2     ~ dnorm(0.6931,0.5)
    omega3     ~ dnorm(0,0.5)
    omega4     ~ dnorm(0,0.5)
    omega5     ~ dnorm(0,0.5)
    omega6     ~ dnorm(0.6931,0.5)
    omega7     ~ dnorm(0.6931,0.5)
    omega8     ~ dnorm(0,0.5)
    omega9     ~ dnorm(0,0.5)
    omega10    ~ dnorm(0,0.5)
    omega11    ~ dnorm(0,0.5)
    omega12    ~ dnorm(0,0.5)
    omega13    ~ dnorm(0,0.5)
    omega14    ~ dnorm(0,0.5)
    omega15    ~ dnorm(0,0.5)
    omega16    ~ dnorm(1.386,0.5)
    omega17    ~ dnorm(0,0.5)
    omega18    ~ dnorm(0,0.5)
    omega19    ~ dnorm(0,0.5)
}
" # close quote for modelString
writeLines( modelString , con="TEMPmodel4.txt" )

initsList = list(beta0=0, omega1=0, omega2=0, omega3=0, omega4=0, omega5=0, omega6=0,
                 omega7=0, omega8=0, omega9=0, omega10=0, omega11=0, omega12=0, omega13=0,
                 omega14=0, omega15=0, omega16=0, omega17=0, omega18=0, omega19=0)

# Run the chains:
jagsModel = jags.model( file="TEMPmodel4.txt" , 
                        data=dataList , 
                        inits=initsList , 
                        n.chains=3 , 
                        n.adapt=500 )

runJagsOut <- run.jags( method="parallel",
                        model="TEMPmodel4.txt",
                        monitor=c("beta0",
                                  "omega1", "omega2", "omega3", "omega4", "omega5", "omega6",
                                  "omega7", "omega8", "omega9", "omega10", "omega11", "omega12", "omega13",
                                  "omega14", "omega15", "omega16", "omega17", "omega18", "omega19"),
                        data=dataList,
                        inits=initsList,
                        n.chains=3,
                        adapt=500,
                        burnin=500,
                        sample=33340,
                        summarise=FALSE,
                        plots=FALSE)

codaSamples = as.mcmc.list( runJagsOut )
save( codaSamples , file=paste0("y4flags_50r2","coda.Rdata") )

s1 <- summary(codaSamples)
s1

#dicSamples <- extract(runJagsOut, what='dic') 
#save( dicSamples, file=paste0("y4flags_50r2","DIC.Rdata"))

library(ggmcmc)
s = ggs(codaSamples)
save (s, file=paste0("y4flags_50r2", "ggs.Rdata"))

#d = ggs_density(s)

#c = ggs_crosscorrelation(s)

#fileNameRoot = "y4flags_50r2"

# Examine the chains:
# Convergence diagnostics:
#diagMCMC( codaObject=codaSamples , parName="beta0")
#saveGraph( file=paste0(fileNameRoot,"beta0Diag") , type="jpg" )
#diagMCMC( codaObject=codaSamples , parName="omega1")
#saveGraph( file=paste0(fileNameRoot,"alpha1Diag") , type="jpg" )
#diagMCMC( codaObject=codaSamples , parName="omega2" )
#saveGraph( file=paste0(fileNameRoot,"sigmaDiag") , type="jpg" )

# Posterior descriptives:
#openGraph(height=3,width=4)
#par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
#plotPost( codaSamples[,"beta0"] , main="beta0" , xlab=bquote(beta0), cenTend="mean" )
#saveGraph( file=paste0(fileNameRoot,"beta0Post") , type="jpg" )
#openGraph(height=3,width=4)
#par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
#plotPost( codaSamples[,"sigma"] , main="sigma" , xlab=bquote(sigma), cenTend="mean" )
#saveGraph( file=paste0(fileNameRoot,"sigmaPost") , type="jpg" )
#openGraph(height=3,width=4)
#par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
#plotPost( codaSamples[,"beta1"] , main="beta1" , xlab=bquote(beta1), cenTend="mean" )
#saveGraph( file=paste0(fileNameRoot,"beta1Post") , type="jpg" )

end_time <- Sys.time()
end_time - start_time

