# NSDUH PCA Analysis
# Matthew J. Beattie
# August 2018

# This script conducts decision tree analysis of the NSDUH dataset
library(rpart)          # for decision tree modeling
library(party)          # for visualizing trees
library(partykit)       # for visualizing trees
library(ggplot2)        # for graphics
library(ROCR)           # for graphics
library(rattle)      		# fancy tree plot
library(rpart.plot)			# enhanced tree plots
library(RColorBrewer)		# color selection for fancy tree plot
library(randomForest)   # Random Forest classification tree method
library(VIM)            # Missing data analysis
library(mice)           # Missing data analysis
library(caret)          # Statisical library
library(dplyr)          # Used for mutation of data frames
library(FactoMineR)     # Robust PCA packages that includes MCA

# Load analysis main dataset
load("NSDUH_analysis_file.Rda")


#********** Perform PCA *****************************
# Create z-scored dataset for PCA
scaledf <- as.data.frame(scale(newdf))

# Create 25% random sample
smp_size <- floor(0.25 * nrow(scaledf))
set.seed(123)
train_ind <- sample(seq_len(nrow(scaledf)), size = smp_size)
scaledf_samp <- scaledf[train_ind, ]

#********** Run 1 **********
pca1<- PCA(scaledf_samp, ncp=10, graph = FALSE)
saveRDS(pca1, "pca1.rds")

# Save eigenvalues to file
output<- capture.output(pca1$eig[1:300,])
cat(output, file="pca1_eig.txt", sep="\n")

# Show correlation to dimensions and save to file
output<- capture.output(dimdesc(pca1, axes=c(1:10)))
cat(output, file="pca1_dim.txt", sep="\n")
write.csv(file="pca1_com.csv", pca1$var$cor)


#*********** Run 2 **********
# Drop CATAG3 age category variable
newdf <- subset(newdf, select = -c(CATAG3))

# Create MHSICKEVR (mental health sick ever) and PHYSICKEVR (physical sick ever)
newdf$MHSICKEVR<- ifelse((newdf$anxdlif==1 | newdf$anxdyr==1 | newdf$deprslif==1 | newdf$deprsyr==1),1,
                         ifelse((newdf$anxdlif==0 & newdf$anxdyr==0 & newdf$deprslif==0 & newdf$deprsyr==0),0,-9))

newdf$PHYSICKEVR<- ifelse((newdf$asmalif==1 | newdf$asmayr==1 | newdf$bronclif==1 | newdf$broncyr==1 | 
                             newdf$cirrlif==1 | newdf$cirryr==1 | newdf$diablif==1 | newdf$diabyr==1 |
                             newdf$hartdlif==1 | newdf$hartdy==1 | newdf$hepatlif==1 | newdf$hepatyr==1 |
                             newdf$hbplif==1 | newdf$hbpyr==1 | newdf$hivlif==1 | newdf$hivyr==1 |
                             newdf$luncalif==1 | newdf$luncayr==1 | newdf$pancrlif==1 | newdf$pancryr==1 |
                             newdf$pneulif==1 | newdf$pneuyr==1 | newdf$stdslif==1 | newdf$stdsyr==1 |
                             newdf$sinuslif==1 | newdf$sinusyr==1 | newdf$slpaplif==1 | newdf$slpapyr==1 |
                             newdf$stroklif==1 | newdf$strokyr==1 | newdf$tinnlif==1 | newdf$tinnyr==1 |
                             newdf$tubrclif==1 | newdf$tubrcyr==1 | newdf$ulcerlif==1 | newdf$ulceryr==1),1,
                          ifelse((newdf$asmalif==0 & newdf$asmayr==0 & newdf$bronclif==0 & newdf$broncyr==0 & 
                                    newdf$cirrlif==0 & newdf$cirryr==0 & newdf$diablif==0 & newdf$diabyr==0 &
                                    newdf$hartdlif==0 & newdf$hartdy==0 & newdf$hepatlif==0 & newdf$hepatyr==0 &
                                    newdf$hbplif==0 & newdf$hbpyr==0 & newdf$hivlif==0 & newdf$hivyr==0 &
                                    newdf$luncalif==0 & newdf$luncayr==0 & newdf$pancrlif==0 & newdf$pancryr==0 &
                                    newdf$pneulif==0 & newdf$pneuyr==0 & newdf$stdslif==0 & newdf$stdsyr==0 &
                                    newdf$sinuslif==0 & newdf$sinusyr==0 & newdf$slpaplif==0 & newdf$slpapyr==0 &
                                    newdf$stroklif==0 & newdf$strokyr==0 & newdf$tinnlif==0 & newdf$tinnyr==0 &
                                    newdf$tubrclif==0 & newdf$tubrcyr==0 & newdf$ulcerlif==0 & newdf$ulceryr==0),0,-9))

# Drop specific physical and mental health variables
newdf<- subset(newdf, select = -c(anxdlif,anxdyr,asmalif,asmayr,bronclif,broncyr,cirrlif,cirryr,
                                      deprslif,deprsyr,diablif,diabyr,hartdlif,hartdyr,hepatlif,
                                      hepatyr,hbplif,hbpyr,hivlif,hivyr,luncalif,luncayr,pancrlif,
                                      pancryr,pneulif,pneuyr,stdslif,stdsyr,slpaplif,slpapyr,
                                      stroklif,strokyr,tinnlif,tinnyr,tubrclif,tubrcyr,ulcerlif,
                                      ulceryr,sinusyr,sinuslif))

dfmiss<-newdf
dfmiss.agg<-aggr(dfmiss)
summary(dfmiss.agg)

# Separate file into YOUTH and ADULT datasets and remove variables unique to the other
dfYouth<-newdf[newdf$YOUTH==1,]
dfYouth<- subset(dfYouth, select = -c(AMHINP2,AMHOUTP3,AMHRX2,AMHTXND2,K6SMXADJ,K6SCMAX,spdyradj,
                                      spdyr,ajamdelt,ajamdeyr,ahltmde,aaltmde,arxmdeyr,EDUCCAT2,
                                      YOUTH,AGE12T17,AGE18T25,AGE26T34,AGE35T49,AGE50PLUS,AMHTXRC3,MHRCOST2,
                                      MHRNBRS2,MHRJOBS2,MHRNCOV2,MHRENUF2,MHRWHER2,MHRCFID2,MHRCMIT2,
                                      MHRNOND2,MHRHAND2,MHRNOHP2,MHRTIME2,MHRFOUT2,MHRTRAN2,MHRSOTH2,
                                      SNMOV5Y2,snysell,snystole,snyattak,snrlgsvc,snrlgimp,snrldcsn))

dfAdult<-newdf[newdf$YOUTH==0,]
dfAdult<- subset(dfAdult, select = -c(MVIN5YR2,schfelt,avggrade,stndscig,stndsmj,stndalc,
                                      parchkhw,talkprob,PRTALK3,YTHACT2,DRPRVME3,ANYEDUC3,
                                      rlgattd,rlgimpt,rlgdcsn,anymhin,anymhout,ymdelt,ymdeyr,
                                      YOUTH,AGE12T17,imother,ifather))

# Create z-scored datasets for PCA
dfYouth_scale <- as.data.frame(scale(dfYouth))
dfAdult_scale <- as.data.frame(scale(dfAdult))

# Create 25% random samples
smp_size <- floor(0.25 * nrow(dfYouth_scale))
set.seed(123)
train_ind <- sample(seq_len(nrow(dfYouth_scale)), size = smp_size)
dfYouth_samp <- dfYouth_scale[train_ind, ]
smp_size <- floor(0.25 * nrow(dfAdult_scale))
train_ind <- sample(seq_len(nrow(dfAdult_scale)), size = smp_size)
dfAdult_samp <- dfAdult_scale[train_ind, ]

# Run PCA for Youth file
pca2_yth<- PCA(dfYouth_samp, ncp=10, graph = FALSE)
saveRDS(pca2_yth, "pca2_yth.rds")

# Save eigenvalues to file
output<- capture.output(pca2_yth$eig[1:300,])
cat(output, file="pca2_yth_eig.txt", sep="\n")

# Show correlation to dimensions and save to file
output<- capture.output(dimdesc(pca2_yth, axes=c(1:10)))
cat(output, file="pca2_yth_dim.txt", sep="\n")
write.csv(file="pca2_yth_com.csv", pca2_yth$var$cor)

# Run PCA for Adult file
pca2_adt<- PCA(dfAdult_samp, ncp=10, graph = FALSE)
saveRDS(pca2_adt, "pca2_adt.rds")

# Save eigenvalues to file
output<- capture.output(pca2_adt$eig[1:300,])
cat(output, file="pca2_adt_eig.txt", sep="\n")

# Show correlation to dimensions and save to file
output<- capture.output(dimdesc(pca2_adt, axes=c(1:10)))
cat(output, file="pca2_adt_dim.txt", sep="\n")
write.csv(file="pca2_adt_com.csv", pca2_adt$var$cor)


#*********** Run 3 ************
# RELOAD NEWDF, DFYOUTH, AND DFADULT FROM RUNS 1 AND 2 FIRST!!!
# Drop reasons for not receiving Mental Health and Illicit Drug/Alcohol Treatment, psychotheratpeutic category
dfYouth<- subset(dfYouth, select = -c(TXLTALC2,TXLTANL2,TXLTCOC2,
                                      TXLTHAL2,TXLTINH2,TXLTMJ2,TXLTSED2,TXLTSTM2,TXLTTRN2,
                                      txrnbusy,txrnfout,txrnhndl,txrnjob,txrnlmcv,txrnnbr,
                                      txrnnhcv,txrnnhlp,txrnnond,txrnnstp,txrnntsp,txrnopen,
                                      txrnstig,txrntype,txrnwher,cpnpsyfg,cpnpsymn,cpnpsyyr,
                                      PSYAGE2,OXYCONT2,totmj,totcoke,totcrack,totdrink,tothall,
                                      totinhal))

dfAdult<- subset(dfAdult, select = -c(MHRCFID2,MHRCMIT2,MHRCOST2,MHRENUF2,MHRFOUT2,MHRHAND2,
                                      MHRJOBS2,MHRNBRS2,MHRNCOV2,MHRNOHP2,MHRNOND2,MHRSOTH2,
                                      MHRTIME2,MHRTRAN2,MHRWHER2,TXLTALC2,TXLTANL2,TXLTCOC2,
                                      TXLTHAL2,TXLTINH2,TXLTMJ2,TXLTSED2,TXLTSTM2,TXLTTRN2,
                                      txrnbusy,txrnfout,txrnhndl,txrnjob,txrnlmcv,txrnnbr,
                                      txrnnhcv,txrnnhlp,txrnnond,txrnnstp,txrnntsp,txrnopen,
                                      txrnstig,txrntype,txrnwher,cpnpsyfg,cpnpsymn,cpnpsyyr,
                                      PSYAGE2,OXYCONT2,totmj,totcoke,totcrack,totdrink,tothall,
                                      totinhal))

# Create z-scored datasets for PCA
dfYouth_scale <- as.data.frame(scale(dfYouth))
dfAdult_scale <- as.data.frame(scale(dfAdult))

# Create 25% random samples
smp_size <- floor(0.25 * nrow(dfYouth_scale))
set.seed(123)
train_ind <- sample(seq_len(nrow(dfYouth_scale)), size = smp_size)
dfYouth_samp <- dfYouth_scale[train_ind, ]
smp_size <- floor(0.25 * nrow(dfAdult_scale))
train_ind <- sample(seq_len(nrow(dfAdult_scale)), size = smp_size)
dfAdult_samp <- dfAdult_scale[train_ind, ]

# Run PCA for Youth file
pca3_yth<- PCA(dfYouth_samp, ncp=10, graph = FALSE)
saveRDS(pca3_yth, "pca3_yth.rds")

# Save eigenvalues to file
output<- capture.output(pca3_yth$eig)
cat(output, file="pca3_yth_eig.txt", sep="\n")

# Show correlation to dimensions and save to file
output<- capture.output(dimdesc(pca3_yth, axes=c(1:10)))
cat(output, file="pca3_yth_dim.txt", sep="\n")
write.csv(file="pca3_yth_com.csv", pca3_yth$var$cor)

# Run PCA for Adult file
pca3_adt<- PCA(dfAdult_samp, ncp=10, graph = FALSE)
saveRDS(pca3_adt, "pca3_adt.rds")

# Save eigenvalues to file
output<- capture.output(pca3_adt$eig[1:300,])
cat(output, file="pca3_adt_eig.txt", sep="\n")

# Show correlation to dimensions and save to file
output<- capture.output(dimdesc(pca3_adt, axes=c(1:10)))
cat(output, file="pca3_adt_dim.txt", sep="\n")
write.csv(file="pca3_adt_com.csv", pca3_adt$var$cor)


# ****** Save final Youth and Adult Datasets for further use *****
saveRDS(dfYouth, "postPCA_yth.rds")
saveRDS(dfAdult, "postPCA_adt.rds")

# ***** Save post-PCA datasets by year for regression tree analysis
# Save dataset into files by period
df03to05YTH<-dfYouth[dfYouth$PERIODNUM==1,]
df06to08YTH<-dfYouth[dfYouth$PERIODNUM==2,]
df09to11YTH<-dfYouth[dfYouth$PERIODNUM==3,]
df12to14YTH<-dfYouth[dfYouth$PERIODNUM==4,]
df03to05ADLT<-dfAdult[dfAdult$PERIODNUM==1,]
df06to08ADLT<-dfAdult[dfAdult$PERIODNUM==2,]
df09to11ADLT<-dfAdult[dfAdult$PERIODNUM==3,]
df12to14ADLT<-dfAdult[dfAdult$PERIODNUM==4,]

# Save datasets
save(df03to05YTH,file="RF_data_03to05YTH.Rda")
save(df06to08YTH,file="RF_data_06to08YTH.Rda")
save(df09to11YTH,file="RF_data_09to11YTH.Rda")
save(df12to14YTH,file="RF_data_12to14YTH.Rda")
save(df03to05ADLT,file="RF_data_03to05ADLT.Rda")
save(df06to08ADLT,file="RF_data_06to08ADLT.Rda")
save(df09to11ADLT,file="RF_data_09to11ADLT.Rda")
save(df12to14ADLT,file="RF_data_12to14ADLT.Rda")

