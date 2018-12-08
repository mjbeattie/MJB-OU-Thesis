# y3poly_50r3.R
# Matthew J. Beattie
# November 2018

# This script examines top drug, treatment, and polyabuse patterns
# Years:  2009-2011
# Setup:  
# 1) Pull out rows with no OxyContin data.
# 2) Model: P(heroin) = Bernoulli(mu), logit(mu) = lin(drug flags)
# 3) Priors set based upon assumptions of three categories:  Uncertain, mean=0, Var=0.5;
#    Probably positive, mean=0.6931, Var=0.5;  Probability Strong, mean=1.386, Var=0.5
# 4) Examines polyabuse, age of first use, treatment, perceived risk, access, 
#    needle use, and frequency of use


#********** Setup dataset and samples ************
# Load dataset
start_time <- Sys.time()
setwd("~/MCMC2")
load("RF_data_09to11YTH.Rda")
newdf <- df09to11YTH
rm(df09to11YTH)

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
newdf$polyabuse <- newdf$cocflag + newdf$cpnmthfg + newdf$crkflag + newdf$ecsflag +
                   newdf$METHDON2 + newdf$MORPHIN2 + newdf$otherpain + newdf$PSILCY2 +
                   newdf$trqflag

# Create a needle use flag
newdf$needleuse <- ifelse(newdf$cocneedl==1 | newdf$mthneedl==1 | newdf$otdgnedl==1,1,0)

# Create a minimum AFU variable
newdf$minafu <- pmin(newdf$ircocage, newdf$irmthage, newdf$ircrkage, newdf$irecsage,
                     newdf$iranlage, newdf$irhalage, newdf$irtrnage)

# Create a maximum frequency variable
newdf$maxfy <- pmax(newdf$ircocfy, newdf$irmthfy, newdf$ircrkfy, newdf$iranlfy,
                    newdf$irhalfy, newdf$irtrnfy)

# Reduce dataset to drug use variables
varlist<-c("polyabuse","needleuse","minafu","txilalev","grskhtry","grskhreg",
           "rdifher","maxfy","herflag")
smalldf<-newdf
rm(newdf)

# Set missing AFU flag equal to mean for present values
afusamp <- smalldf[smalldf$minafu<991,]
avgafu <- round(mean(afusamp$minafu))
smalldf$minafu <- ifelse(smalldf$minafu==991,avgafu,smalldf$minafu)
rm(afusamp)

# Set missing FY flag equal to mean for present values
fysamp <- smalldf[smalldf$maxfy>0,]
avgfy <- round(mean(fysamp$maxfy))
smalldf$maxfy <- ifelse(smalldf$maxfy==0,avgfy,smalldf$maxfy)
rm(fysamp)

# Set perceived risk of heroin use to 2 if regular ok, 1 if try ok, else 0
smalldf$herrisk <- ifelse(smalldf$grskhreg==0,2,
                          ifelse(smalldf$grskhtry==0,1,0))

# Take random sample
frac <- 0.70
smp_size <- floor(frac * nrow(smalldf))
set.seed(456)
test_ind <- sample(seq_len(nrow(smalldf)), size = smp_size)
smalldf <- smalldf[test_ind, ]

# Prepare data for JAGS by creating list
dataList = list(y=smalldf$herflag, g1=smalldf$polyabuse, g2=smalldf$needleuse,
                g4=smalldf$txilalev, g5=smalldf$herrisk,
                g6=smalldf$rdifher, N=length(smalldf$herflag))



#************** RUN JAGS ROUTINES ***************
# Uses runjags to take advantage of multiple CPU cores
library(runjags)
source("DBDA2E-utilities.R")

#Define the model:
modelString = "
data 
{
  g1bar <- mean(g1)
  g1sd <- sd(g1)
  g5bar <- mean(g5)
  g5sd <- sd(g5)
  for (i in 1:N)
  {
    zg1[i] <- (g1[i]-g1bar)/g1sd
    zg5[i] <- (g5[i]-g5bar)/g5sd
  }
}
model
{
  for (i in 1:N) 
  {
    y[i]  ~ dbern(mu[i])
    logit(mu[i])  <- zeta0 +
              zeta1*zg1[i] + alpha2*g2[i] + alpha4*g4[i] + 
              zeta5*zg5[i] + alpha6*g6[i]
  }
  # Set priors based upon Greenland method
  zeta0      ~ dnorm(0,0.5)
  zeta1      ~ dnorm(0.6931,0.5)
  alpha2     ~ dnorm(1.386,0.5)
  alpha4     ~ dnorm(0,0.5)
  zeta5      ~ dnorm(0,0.5)
  alpha6     ~ dnorm(0.6931,0.5)

  # Transform to original scale
  beta0 <- zeta0 - (zeta1*g1bar/g1sd) - (zeta5*g5bar/g5sd)
  alpha1 <- zeta1/g1sd
  alpha5 <- zeta5/g5sd
}
" # close quote for modelString
writeLines( modelString , con="TEMPmodelP3.txt" )

initsList = list(zeta0=0, zeta1=0, alpha2=0, alpha4=0, zeta5=0, alpha6=0)

# Run the chains:
jagsModel = jags.model( file="TEMPmodelP3.txt" , 
                        data=dataList , 
                        inits=initsList , 
                        n.chains=3 , 
                        n.adapt=500 )

runJagsOut <- run.jags( method="parallel",
                        model="TEMPmodelP3.txt",
                        monitor=c("beta0","alpha1", "alpha2", 
                                  "alpha4", "alpha5", "alpha6"),
                        data=dataList,
                        inits=initsList,
                        n.chains=3,
                        adapt=500,
                        burnin=500,
                        sample=33340,
                        summarise=FALSE,
                        plots=FALSE)

codaSamples = as.mcmc.list( runJagsOut )
save( codaSamples , file=paste0("y3poly_50r3","coda.Rdata") )

s1 <- summary(codaSamples)
s1

#dicSamples <- extract(runJagsOut, what='dic') 
#save( dicSamples, file=paste0("y3poly_50r3","DIC.Rdata"))

library(ggmcmc)
s = ggs(codaSamples)
save (s, file=paste0("y3poly_50r3", "ggs.Rdata"))

#d = ggs_density(s)

#c = ggs_crosscorrelation(s)
#c

#fileNameRoot = "y3poly_50r3"

# Examine the chains:
# Convergence diagnostics:
#diagMCMC( codaObject=codaSamples , parName="beta0")
#saveGraph( file=paste0(fileNameRoot,"beta0Diag") , type="jpg" )
#diagMCMC( codaObject=codaSamples , parName="alpha1")
#saveGraph( file=paste0(fileNameRoot,"alpha1Diag") , type="jpg" )
#diagMCMC( codaObject=codaSamples , parName="alpha3" )
#saveGraph( file=paste0(fileNameRoot,"alpha3Diag") , type="jpg" )

# Posterior descriptives:
#openGraph(height=3,width=4)
#par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
#plotPost( codaSamples[,"beta0"] , main="beta0" , xlab=bquote(beta0), cenTend="mean" )
#saveGraph( file=paste0(fileNameRoot,"beta0Post") , type="jpg" )
#openGraph(height=3,width=4)
#par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
#plotPost( codaSamples[,"alpha1"] , main="alpha1" , xlab=bquote(alpha1), cenTend="mean" )
#saveGraph( file=paste0(fileNameRoot,"alpha1Post") , type="jpg" )
#openGraph(height=3,width=4)
#par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
#plotPost( codaSamples[,"alpha3"] , main="alpha3" , xlab=bquote(alpha3), cenTend="mean" )
#saveGraph( file=paste0(fileNameRoot,"alpha3Post") , type="jpg" )

end_time <- Sys.time()
end_time - start_time

