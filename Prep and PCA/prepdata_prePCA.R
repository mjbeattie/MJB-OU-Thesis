# NSDUH Dataset Prepartion
# Matthew J. Beattie
# September 2018

# This script reads in the 2002-2014 NSDUH public use file and extracts the
# variables required for CART and MCMC analysis.

library(VIM)          # R package to investigate missingness
library(mice)         # R package to investigate missingness

#********** Read in all variables of interest into a new dataframe **********#
alldata<-CONCATENATEDPUF_092116

# Create list of desired variables and save to new small dataframe
varlist<-c("year","ircigrc","ircgrrc","irpiplf","irsltrc","irchwrc","irsnfrc",
           "iralcrc","irmjrc","ircocrc","ircrkrc","irherrc","irhalrc","irlsdrc",
           "irpcprc","irecsrc","irinhrc","iranlrc","iroxyrc","irtrnrc","irstmrc",
           "irmthrc","irsedrc","iralcfy","irmjfy","ircocfy","ircrkfy","irherfy",
           "irhalfy","irinhfy","iranlfy","iroxyfy","irtrnfy","irstmfy","irmthfy",
           "irsedfy","ircigfm","ircgrfm","irchwfm","irsnffm","iralcfm","IRALCD5",
           "irmjfm","ircocfm","ircrkfm","irherfm","irhalfm","irinhfm","ircigage",
           "ircduage","ircgrage","irsltage","irchwage","irsnfage","iralcage","irmjage",
           "ircocage","ircrkage","irherage","irhalage","irlsdage","irpcpage","irecsage",
           "irinhage","iranlage","iroxyage","irtrnage","irstmage","irmthage","irsedage",
           "cigflag","cigyr","cigmon","cgrflag","cgryr","cgrmon","pipflag",
           "pipmon","smkflag","smkyr","smkmon","chwflag","chwyr","chwmon",
           "snfflag","snfyr","snfmon","alcflag","alcyr","alcmon","OXYCONT2",
           "mrjflag","mrjyr","mrjmon","cocflag","cocyr","cocmon","crkflag",
           "crkyr","crkmon","herflag","heryr","hermon","halflag","halyr",
           "halmon","lsdflag","lsdyr","lsdmon","pcpflag","pcpyr","pcpmon",
           "ecsflag","ecsyr","ecsmon","inhflag","inhyr","inhmon","anlflag",
           "anlyr","anlmon","oxyflag","oxyyr","oxymon","trqflag","trqyr",
           "trqmon","cpnstmfg","cpnstmyr","cpnstmmn","cpnmthfg","cpnmthyr","cpnmthmn",
           "sedflag","sedyr","sedmon","PSYAGE2","PSYYFU2","cpnpsyfg","cpnpsyyr",
           "cpnpsymn","bingehvy","PEYOTE2","MESC2","PSILCY2","AMYLNIT2","CLEFLU2",
           "GAS2","GLUE2","ETHER2","SOLVENT2","LGAS2","NITOXID2","SPPAINT2",
           "AEROS2","DARVTYL2","PERCTYL2","VICOLOR2","CODEINE2","DEMEROL2","DILAUD2",
           "FIORICT2","FIORINL2","HYDROCD2","METHDON2","MORPHIN2","PHENCOD2","PROPOXY2",
           "SK65A2","STADOL2","TALACEN2","TALWIN2","TALWINX2","TRAMADL2","ULTRAM2",
           "othanl","OXYCODP2","hydcodop","tramadp","KLONOPI2","XNAXATV2","VALMDIA2",
           "ATARAX2","BUSPAR2","EQUANIL2","FLEXERL2","LIBRIUM2","LIMBTRL2","MEPROB2",
           "MILTOWN2","ROHYPNL2","SERAX2","SOMA2","TRANXEN2","VISTAR2","othtrn",
           "benzos","meprobpd","muscrelx","METHDES2","DIETPIL2","RITMPHE2","CYLERT2",
           "DEXED2","DETAMP2","DIDREX2","ESKAT2","IONAMIN2","MAZANOR2","OBLA2",
           "PLEGINE2","PRELUDN2","SANOREX2","TENUATE2","othstm","amdxphen","mazindol",
           "METHAQ2","NEMBBAR2","RESTTMA2","AMYTAL2","BUTISOL2","CHHYD2","DALMANE2",
           "HALCION2","PHENOBR2","PLACIDY2","TUINAL2","othsed","rtdalhal","anybarb",
           "totdrink","totmj","totcoke",
           "totcrack","tothero","tothall","totinhal","PR11","PR21","TR11",
           "ST11","ST21","SV11","mthneedl","ostneedl","cocneedl","otdgnedl",
           "herneedl","grskpkcg","grskmocc","grskmreg","grskcocc","grskcreg","grskhtry",
           "grskhreg","grskltry","grsklreg","GRSKD4_5","GRSKD5WK","rdifmj","rdifcoc",
           "rdifcrk","rdifher","rdiflsd","appseldg","ndssdnsp","depndalc","depndanl",
           "depndcoc","depndhal","depndher","depndinh","depndmrj","depndsed","depndstm",
           "depndtrn","depndill","abusealc","abuseanl","abusecoc","abusehal","abuseher",
           "abuseinh","abusemrj","abusesed","abusestm","abusetrn","abuseill","booked",
           "NOBOOKY2","mxmjpnlt","anlltsc","TRNLTS2","STMLTS2","MTHLTSC2",
           "SEDLTS2","txilalev","txillalc","TXLTALC2","TXLTMJ2","TXLTCOC2","TXLTHER2",
           "TXLTHAL2","TXLTINH2","TXLTANL2","TXLTTRN2","TXLTSTM2","TXLTSED2","TXPINS2",
           "TXPCARE2","TXPCAID2","TXPPUBP2","TXPSAVE2","TXPFMLY2","TXPCORT2","TXPMILC2",
           "TXPEMPL2","lochosp","locrfin","locrfop","locmhc","loctxer","locdoc",
           "locjail","locshg","txndilal","flndilal","ndefia","txgpilal","txrnnhcv",
           "txrnlmcv","txrnntsp","txrntype","txrnnstp","txrnopen","txrnwher","txrnnbr",
           "txrnjob","txrnnond","txrnhndl","txrnnhlp","txrnbusy","txrnfout","txrnstig",
           "preg","anxdlif","anxdyr","asmalif","asmayr","bronclif","broncyr",
           "cirrlif","cirryr","deprslif","deprsyr","diablif","diabyr","hartdlif",
           "hartdyr","hepatlif","hepatyr","hbplif","hbpyr","hivlif","hivyr",
           "luncalif","luncayr","pancrlif","pancryr","pneulif","pneuyr","stdslif",
           "stdsyr","sinuslif","sinusyr","slpaplif","slpapyr","stroklif","strokyr",
           "tinnlif","tinnyr","tubrclif","tubrcyr","ulcerlif","ulceryr","AMHINP2",
           "AMHOUTP3","AMHRX2","AMHTXRC3","AMHTXND2","MHRCOST2","MHRNBRS2","MHRJOBS2",
           "MHRNCOV2","MHRENUF2","MHRWHER2","MHRCFID2","MHRCMIT2","MHRNOND2","MHRHAND2",
           "MHRNOHP2","MHRTIME2","MHRFOUT2","MHRTRAN2","MHRSOTH2","SNMOV5Y2",
           "snysell","snystole","snyattak","snrlgsvc",
           "snrlgimp","snrldcsn","MVIN5YR2","schfelt","avggrade","stndscig","stndsmj",
           "stndalc","parchkhw","talkprob","PRTALK3","YTHACT2","DRPRVME3","ANYEDUC3",
           "rlgattd","rlgimpt","rlgdcsn","K6SMXADJ","K6SCMAX","spdyradj","spdyr",
           "ajamdelt","ajamdeyr","ahltmde","aaltmde","arxmdeyr","anymhin","anymhout",
           "ymdelt","ymdeyr","IRHH65_2","IRHHSIZ2","IRKI17_2","imother","ifather",
           "IRFAMIN3","IRPINC3","govtprog","irmcdchp","irprvhlt","IRINSUR4","irfamwag",
           "langver","AGE2","LFTSCHA2","MOVESPY2","NOMARR2","schtype","schenrl",
           "service","milstat","sdntftpt","actdever","health","irsex","NEWRACE2",
           "irmarit","EDUCCAT2","empstaty","COUTYP2","MAIIN002","AWTC10_C","MHSAWT_C",
           "SPDWT_C","DEPWT_C","ANALWC1","ANALWC2","ANALWC3","ANALWC4","ANALWC5",
           "ANALWC6","ANALWC7","ANALWC8","ANALWC9","ANALWC10","ANALWC11","ANALWC12",
           "ANALWC13","CATAG3","parol","prob")

smalldf<-alldata[,varlist]
save(smalldf,file="NSDUH_small_file.Rda")

# Extract duplicative imputed drug use variables
newdf<- subset(smalldf, select = -c(PERCTYL2,VICOLOR2,HYDROCD2,TRAMADL2,ULTRAM2,
                                    KLONOPI2,XNAXATV2,VALMDIA2,LIBRIUM2,LIMBTRL2,
                                    ROHYPNL2,SERAX2,TRANXEN2,EQUANIL2,MEPROB2,
                                    MILTOWN2,DIETPIL2,DEXED2,DETAMP2,IONAMIN2,
                                    AMYTAL2,BUTISOL2,PHENOBR2,TUINAL2,MAZANOR2,SANOREX2,
                                    FLEXERL2,SOMA2) )

# Remove 2002 to make evenly sized time periods
newdf<-newdf[!(newdf$year==2002),]

# Create numeric PERIODNUM to represent blocks of years (required for PCA, etc)
newdf$PERIODNUM<-ifelse(newdf$year < 2006, 1,
                     ifelse(newdf$year < 2009, 2,
                            ifelse(newdf$year < 2012, 3, 4)))


# Convert frequency of use flags to zero for none in period or ever
newdf$iralcfy<-ifelse((newdf$iralcfy==991 | newdf$iralcfy==993),0,newdf$iralcfy)
newdf$irmjfy<-ifelse((newdf$irmjfy==991 | newdf$irmjfy==993),0,newdf$irmjfy)
newdf$ircocfy<-ifelse((newdf$ircocfy==991 | newdf$ircocfy==993),0,newdf$ircocfy)
newdf$ircrkfy<-ifelse((newdf$ircrkfy==991 | newdf$ircrkfy==993),0,newdf$ircrkfy)
newdf$irhalfy<-ifelse((newdf$irhalfy==991 | newdf$irhalfy==993),0,newdf$irhalfy)
newdf$irinhfy<-ifelse((newdf$irinhfy==991 | newdf$irinhfy==993),0,newdf$irinhfy)
newdf$iranlfy<-ifelse((newdf$iranlfy==991 | newdf$iranlfy==993),0,newdf$iranlfy)
newdf$iroxyfy<-ifelse((newdf$iroxyfy==991 | newdf$iroxyfy==993),0,newdf$iroxyfy)
newdf$irtrnfy<-ifelse((newdf$irtrnfy==991 | newdf$irtrnfy==993),0,newdf$irtrnfy)
newdf$irstmfy<-ifelse((newdf$irstmfy==991 | newdf$irstmfy==993),0,newdf$irstmfy)
newdf$irmthfy<-ifelse((newdf$irmthfy==991 | newdf$irmthfy==993),0,newdf$irmthfy)
newdf$irsedfy<-ifelse((newdf$irsedfy==991 | newdf$irsedfy==993),0,newdf$irsedfy)
newdf$iralcfm<-ifelse((newdf$iralcfm==91 | newdf$iralcfm==93),0,newdf$iralcfm)
newdf$irmjfm<-ifelse((newdf$irmjfm==91 | newdf$irmjfm==93),0,newdf$irmjfm)
newdf$ircocfm<-ifelse((newdf$ircocfm==91 | newdf$ircocfm==93),0,newdf$ircocfm)
newdf$ircrkfm<-ifelse((newdf$ircrkfm==91 | newdf$ircrkfm==93),0,newdf$ircrkfm)
newdf$irhalfm<-ifelse((newdf$irhalfm==91 | newdf$irhalfm==93),0,newdf$irhalfm)
newdf$irinhfm<-ifelse((newdf$irinhfm==91 | newdf$irinhfm==93),0,newdf$irinhfm)
newdf$ircigfm<-ifelse((newdf$ircigfm==91 | newdf$ircigfm==93),0,newdf$ircigfm)
newdf$ircgrfm<-ifelse((newdf$ircgrfm==91 | newdf$ircgrfm==93),0,newdf$ircgrfm)
newdf$irchwfm<-ifelse((newdf$irchwfm==91 | newdf$irchwfm==93),0,newdf$irchwfm)
newdf$irsnffm<-ifelse((newdf$irsnffm==91 | newdf$irsnffm==93),0,newdf$irsnffm)
newdf$IRALCD5<-ifelse((newdf$IRALCD5==91 | newdf$IRALCD5==93),0,newdf$IRALCD5)

newdf$totdrink<-ifelse((newdf$totdrink==994 | newdf$totdrink==997 | newdf$totdrink==998),0,newdf$totdrink)
newdf$totmj<-ifelse((newdf$totmj==994 | newdf$totmj==997 | newdf$totmj==998),0,newdf$totmj)
newdf$totcoke<-ifelse((newdf$totcoke==994 | newdf$totcoke==997 | newdf$totcoke==998),0,newdf$totcoke)
newdf$totcrack<-ifelse((newdf$totcrack==994 | newdf$totcrack==997 | newdf$totcrack==998),0,newdf$totcrack)
newdf$tothall<-ifelse((newdf$tothall==994 | newdf$tothall==997 | newdf$tothall==998),0,newdf$tothall)
newdf$totinhal<-ifelse((newdf$totinhal==994 | newdf$totinhal==997 | newdf$totinhal==998),0,newdf$totinhal)
newdf$tothero<-ifelse((newdf$tothero==994 | newdf$tothero==997 | newdf$tothero==998),0,newdf$tothero)

newdf$PR11<-ifelse((newdf$PR11==994 | newdf$PR11==997 | newdf$PR11==998),0,newdf$PR11)
newdf$PR21<-ifelse((newdf$PR21==994 | newdf$PR21==997 | newdf$PR21==998),0,newdf$PR21)
newdf$TR11<-ifelse((newdf$TR11==994 | newdf$TR11==997 | newdf$TR11==998),0,newdf$TR11)
newdf$ST11<-ifelse((newdf$ST11==994 | newdf$ST11==997 | newdf$ST11==998),0,newdf$ST11)
newdf$ST21<-ifelse((newdf$ST21==994 | newdf$ST21==997 | newdf$ST21==998),0,newdf$ST21)
newdf$SV11<-ifelse((newdf$SV11==994 | newdf$SV11==997 | newdf$SV11==998),0,newdf$SV11)


# Convert answers for needle use, crime to yes or no only (converts no=2 to no=0, assume skip=no)
newdf$mthneedl<-ifelse((newdf$mthneedl==1 | newdf$mthneedl==3),1,0)
newdf$ostneedl<-ifelse((newdf$ostneedl==1 | newdf$ostneedl==3),1,0)
newdf$cocneedl<-ifelse((newdf$cocneedl==1 | newdf$cocneedl==3),1,0)
newdf$otdgnedl<-ifelse((newdf$otdgnedl==1 | newdf$otdgnedl==3),1,0)
newdf$NOBOOKY2<-ifelse((newdf$NOBOOKY2>3),0,newdf$NOBOOKY2)
newdf$booked<-ifelse((newdf$booked==1 | newdf$booked==3),1,0)
newdf$SNMOV5Y2<-ifelse((newdf$SNMOV5Y2>6),0,newdf$SNMOV5Y2)
newdf$NOMARR2<-ifelse((newdf$NOMARR2>2),0,newdf$NOMARR2)
newdf$irmarit<-ifelse((newdf$irmarit>4),4,newdf$irmarit)
newdf$empstaty<-ifelse((newdf$empstaty>4),4,newdf$empstaty)

# Convert specific crimes so that skips, etc. are no
newdf$snysell<-ifelse(newdf$snysell>5,1,newdf$snysell)
newdf$snystole<-ifelse(newdf$snystole>5,1,newdf$snystole)
newdf$snyattak<-ifelse(newdf$snyattak>5,1,newdf$snyattak)
newdf$snrlgsvc<-ifelse(newdf$snrlgsvc>6,1,newdf$snrlgsvc)
newdf$snrlgimp<-ifelse(newdf$snrlgimp>4,1,newdf$snrlgimp)
newdf$snrldcsn<-ifelse(newdf$snrldcsn>4,1,newdf$snrldcsn)

# Set military status to yes or no flag
newdf$service<-ifelse(newdf$service==1,1,0)

# Set student status to yes or no flag
newdf$student<-ifelse(newdf$sdntftpt<=2,1,0)
newdf$sdntftpt<-NULL

# Define age categories
newdf$YOUTH<-ifelse(newdf$AGE2<=6,1,0)
newdf$AGE12T17<-ifelse(newdf$AGE2<=6,1,0)
newdf$AGE18T25<-ifelse((newdf$AGE2>=7 & newdf$AGE2<=12),1,0)
newdf$AGE26T34<-ifelse((newdf$AGE2>=13 & newdf$AGE2<=14),1,0)
newdf$AGE35T49<-ifelse((newdf$AGE2==15),1,0)
newdf$AGE50PLUS<-ifelse((newdf$AGE2>= 16),1,0)

# Remove recency flags, which are duplicated by past year, etc.
newdf<- subset(newdf, select = -c(ircigrc,ircgrrc,irpiplf,irsltrc,irchwrc,irsnfrc,
                                    iralcrc,irmjrc,ircocrc,ircrkrc,irhalrc,irlsdrc,
                                    irpcprc,irecsrc,irinhrc,iranlrc,iroxyrc,irtrnrc,
                                    irstmrc,irmthrc,irsedrc) )

# Remove miscellaneous variables
newdf<- subset(newdf, select = -c(year,AGE2,PSYYFU2,MOVESPY2,LFTSCHA2,MOVESPY2,schtype,milstat) )

# Remove weighting factor variables
newdf<- subset(newdf, select = -c(AWTC10_C,MHSAWT_C,SPDWT_C,DEPWT_C,ANALWC1,ANALWC2,
                                  ANALWC3,ANALWC4,ANALWC5,ANALWC6,ANALWC7,ANALWC8,
                                  ANALWC9,ANALWC10,ANALWC11,ANALWC12,ANALWC13))

# Remove variables that are defined by heroin use (e.g. age of first heroin use)
newdf<- subset(newdf, select = -c(irherrc,irherfy,irherfm,irherage,heryr,hermon,tothero,
                                  herneedl,depndher,depndill,abuseher,abuseill,
                                  TXLTHER2) )

# Correct missing variables -- Random Forest can't handle missing values 
newdf$grskpkcg[is.na(newdf$grskpkcg)] <- 0
newdf$grskmocc[is.na(newdf$grskmocc)] <- 0
newdf$grskmreg[is.na(newdf$grskmreg)] <- 0
newdf$grskcocc[is.na(newdf$grskcocc)] <- 0
newdf$grskcreg[is.na(newdf$grskcreg)] <- 0
newdf$grskhtry[is.na(newdf$grskhtry)] <- 0
newdf$grskhreg[is.na(newdf$grskhreg)] <- 0
newdf$grskltry[is.na(newdf$grskltry)] <- 0
newdf$grsklreg[is.na(newdf$grsklreg)] <- 0
newdf$GRSKD4_5[is.na(newdf$GRSKD4_5)] <- 0
newdf$GRSKD5WK[is.na(newdf$GRSKD5WK)] <- 0
newdf$rdifmj[is.na(newdf$rdifmj)] <- 0
newdf$rdifcoc[is.na(newdf$rdifcoc)] <- 0
newdf$rdifcrk[is.na(newdf$rdifcrk)] <- 0
newdf$rdifher[is.na(newdf$rdifher)] <- 0
newdf$rdiflsd[is.na(newdf$rdiflsd)] <- 0
newdf$appseldg[is.na(newdf$appseldg)] <- 0
newdf$parol[is.na(newdf$parol)] <- 2
newdf$prob[is.na(newdf$prob)] <- 2
newdf$anlltsc[is.na(newdf$anlltsc)] <- -9
newdf$TRNLTS2[is.na(newdf$TRNLTS2)] <- -9
newdf$STMLTS2[is.na(newdf$STMLTS2)] <- -9
newdf$MTHLTSC2[is.na(newdf$MTHLTSC2)] <- -9
newdf$SEDLTS2[is.na(newdf$SEDLTS2)] <- -9
newdf$txrnnhcv[is.na(newdf$txrnnhcv)] <- -9
newdf$txrnlmcv[is.na(newdf$txrnlmcv)] <- -9
newdf$txrnntsp[is.na(newdf$txrnntsp)] <- -9
newdf$txrntype[is.na(newdf$txrntype)] <- -9
newdf$txrnnstp[is.na(newdf$txrnnstp)] <- -9
newdf$txrnopen[is.na(newdf$txrnopen)] <- -9
newdf$txrnwher[is.na(newdf$txrnwher)] <- -9
newdf$txrnnbr[is.na(newdf$txrnnbr)] <- -9
newdf$txrnjob[is.na(newdf$txrnjob)] <- -9
newdf$txrnnond[is.na(newdf$txrnnond)] <- -9
newdf$txrnhndl[is.na(newdf$txrnhndl)] <- -9
newdf$txrnnhlp[is.na(newdf$txrnnhlp)] <- -9
newdf$txrnbusy[is.na(newdf$txrnbusy)] <- -9
newdf$txrnfout[is.na(newdf$txrnfout)] <- -9
newdf$txrnstig[is.na(newdf$txrnstig)] <- -9
newdf$anxdlif[is.na(newdf$anxdlif)] <- 0
newdf$anxdyr[is.na(newdf$anxdyr)] <- 0
newdf$asmalif[is.na(newdf$asmalif)] <- 0
newdf$asmayr[is.na(newdf$asmayr)] <- 0
newdf$bronclif[is.na(newdf$bronclif)] <- 0
newdf$broncyr[is.na(newdf$broncyr)] <- 0
newdf$cirrlif[is.na(newdf$cirrlif)] <- 0
newdf$cirryr[is.na(newdf$cirryr)] <- 0
newdf$deprslif[is.na(newdf$deprslif)] <- 0
newdf$deprsyr[is.na(newdf$deprsyr)] <- 0
newdf$diablif[is.na(newdf$diablif)] <- 0
newdf$diabyr[is.na(newdf$diabyr)] <- 0
newdf$hartdlif[is.na(newdf$hartdlif)] <- 0
newdf$hartdyr[is.na(newdf$hartdyr)] <- 0
newdf$hepatlif[is.na(newdf$hepatlif)] <- 0
newdf$hepatyr[is.na(newdf$hepatyr)] <- 0
newdf$hbplif[is.na(newdf$hbplif)] <- 0
newdf$hbpyr[is.na(newdf$hbpyr)] <- 0
newdf$hivlif[is.na(newdf$hivlif)] <- 0
newdf$hivyr[is.na(newdf$hivyr)] <- 0
newdf$luncalif[is.na(newdf$luncalif)] <- 0
newdf$luncayr[is.na(newdf$luncayr)] <- 0
newdf$pancrlif[is.na(newdf$pancrlif)] <- 0
newdf$pancryr[is.na(newdf$pancryr)] <- 0
newdf$pneulif[is.na(newdf$pneulif)] <- 0
newdf$pneuyr[is.na(newdf$pneuyr)] <- 0
newdf$stdslif[is.na(newdf$stdslif)] <- 0
newdf$stdsyr[is.na(newdf$stdsyr)] <- 0
newdf$sinuslif[is.na(newdf$sinuslif)] <- 0
newdf$sinusyr[is.na(newdf$sinusyr)] <- 0
newdf$slpaplif[is.na(newdf$slpaplif)] <- 0
newdf$slpapyr[is.na(newdf$slpapyr)] <- 0
newdf$stroklif[is.na(newdf$stroklif)] <- 0
newdf$strokyr[is.na(newdf$strokyr)] <- 0
newdf$tinnlif[is.na(newdf$tinnlif)] <- 0
newdf$tinnyr[is.na(newdf$tinnyr)] <- 0
newdf$tubrclif[is.na(newdf$tubrclif)] <- 0
newdf$tubrcyr[is.na(newdf$tubrcyr)] <- 0
newdf$ulcerlif[is.na(newdf$ulcerlif)] <- 0
newdf$ulceryr[is.na(newdf$ulceryr)] <- 0
newdf$AMHINP2[is.na(newdf$AMHINP2)] <- 2
newdf$AMHOUTP3[is.na(newdf$AMHOUTP3)] <- 2
newdf$AMHRX2[is.na(newdf$AMHRX2)] <- 2
newdf$AMHTXRC3[is.na(newdf$AMHTXRC3)] <- 2
newdf$AMHTXND2[is.na(newdf$AMHTXND2)] <- 2
newdf$MHRCOST2[is.na(newdf$MHRCOST2)] <- -9
newdf$MHRNBRS2[is.na(newdf$MHRNBRS2)] <- -9
newdf$MHRJOBS2[is.na(newdf$MHRJOBS2)] <- -9
newdf$MHRNCOV2[is.na(newdf$MHRNCOV2)] <- -9
newdf$MHRENUF2[is.na(newdf$MHRENUF2)] <- -9
newdf$MHRWHER2[is.na(newdf$MHRWHER2)] <- -9
newdf$MHRCFID2[is.na(newdf$MHRCFID2)] <- -9
newdf$MHRCMIT2[is.na(newdf$MHRCMIT2)] <- -9
newdf$MHRNOND2[is.na(newdf$MHRNOND2)] <- -9
newdf$MHRHAND2[is.na(newdf$MHRHAND2)] <- -9
newdf$MHRNOHP2[is.na(newdf$MHRNOHP2)] <- -9
newdf$MHRTIME2[is.na(newdf$MHRTIME2)] <- -9
newdf$MHRFOUT2[is.na(newdf$MHRFOUT2)] <- -9
newdf$MHRTRAN2[is.na(newdf$MHRTRAN2)] <- -9
newdf$MHRSOTH2[is.na(newdf$MHRSOTH2)] <- -9
newdf$MVIN5YR2[is.na(newdf$MVIN5YR2)] <- -9
newdf$schfelt[is.na(newdf$schfelt)] <- -9
newdf$avggrade[is.na(newdf$avggrade)] <- -9
newdf$stndscig[is.na(newdf$stndscig)] <- -9
newdf$stndsmj[is.na(newdf$stndsmj)] <- -9
newdf$stndalc[is.na(newdf$stndalc)] <- -9
newdf$parchkhw[is.na(newdf$parchkhw)] <- -9
newdf$talkprob[is.na(newdf$talkprob)] <- -9
newdf$PRTALK3[is.na(newdf$PRTALK3)] <- -9
newdf$YTHACT2[is.na(newdf$YTHACT2)] <- -9
newdf$DRPRVME3[is.na(newdf$DRPRVME3)] <- -9
newdf$ANYEDUC3[is.na(newdf$ANYEDUC3)] <- -9
newdf$rlgattd[is.na(newdf$rlgattd)] <- -9
newdf$rlgimpt[is.na(newdf$rlgimpt)] <- -9
newdf$rlgdcsn[is.na(newdf$rlgdcsn)] <- -9
newdf$K6SMXADJ[is.na(newdf$K6SMXADJ)] <- -9
newdf$K6SCMAX[is.na(newdf$K6SCMAX)] <- -9
newdf$spdyradj[is.na(newdf$spdyradj)] <- -9
newdf$spdyr[is.na(newdf$spdyr)] <- -9
newdf$ajamdelt[is.na(newdf$ajamdelt)] <- -9
newdf$ajamdeyr[is.na(newdf$ajamdeyr)] <- -9
newdf$ahltmde[is.na(newdf$ahltmde)] <- -9
newdf$aaltmde[is.na(newdf$aaltmde)] <- -9
newdf$arxmdeyr[is.na(newdf$arxmdeyr)] <- -9
newdf$anymhin[is.na(newdf$anymhin)] <- -9
newdf$anymhout[is.na(newdf$anymhout)] <- -9
newdf$ymdelt[is.na(newdf$ymdelt)] <- -9
newdf$ymdeyr[is.na(newdf$ymdeyr)] <- -9

# Demonstrate missingness of raw values
dfmiss<-newdf[,1:100]
dfmiss.agg<-aggr(dfmiss)
summary(dfmiss.agg)

# Save modified dataset for C&RT
save(newdf,file="NSDUH_analysis_file.Rda")

