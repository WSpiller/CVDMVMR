######################
### Load libraries ###
######################

library(TwoSampleMR)
library(MRInstruments)
library(RadialMR)
library(MVMR)
library(ggplot2)
library(reshape2)

######################
### Univariable MR ###
######################

#Data extraction

data(gwas_catalog)

exposure_data_Men<-extract_instruments("ukb-a-315",clump = T)
exposure_data_BMI<-extract_instruments("ukb-a-248",clump = T)
exposure_data_SBP<-extract_instruments("ukb-a-360",clump = T)
exposure_data_LDL<-extract_instruments("ieu-a-300",clump = T)
exposure_data_HDL<-extract_instruments("ieu-a-299",clump = T)
exposure_data_TG<-extract_instruments("ieu-a-302",clump = T)

outcome_data_MenCHD<- extract_outcome_data(
  snps = exposure_data_Men$SNP, outcomes = "ieu-a-7",proxies=F)

outcome_data_MenMI<- extract_outcome_data(
  snps = exposure_data_Men$SNP, outcomes = "ieu-a-798",proxies=F)

outcome_data_MenILA<- extract_outcome_data(
  snps = exposure_data_Men$SNP, outcomes = "ebi-a-GCST006907",proxies=F)

outcome_data_MenIC<- extract_outcome_data(
  snps = exposure_data_Men$SNP, outcomes = "ebi-a-GCST006910",proxies=F)

outcome_data_MenISV<- extract_outcome_data(
  snps = exposure_data_Men$SNP, outcomes = "ebi-a-GCST006909",proxies=F)

outcome_data_BMICHD<- extract_outcome_data(
  snps = exposure_data_BMI$SNP, outcomes = "ieu-a-7",proxies=F)

outcome_data_BMIMI<- extract_outcome_data(
  snps = exposure_data_BMI$SNP, outcomes = "ieu-a-798",proxies=F)

outcome_data_BMIILA<- extract_outcome_data(
  snps = exposure_data_BMI$SNP, outcomes = "ebi-a-GCST006907",proxies=F)

outcome_data_BMIIC<- extract_outcome_data(
  snps = exposure_data_BMI$SNP, outcomes = "ebi-a-GCST006910",proxies=F)

outcome_data_BMIISV<- extract_outcome_data(
  snps = exposure_data_BMI$SNP, outcomes = "ebi-a-GCST006909",proxies=F)

outcome_data_SBPCHD<- extract_outcome_data(
  snps = exposure_data_SBP$SNP, outcomes = "ieu-a-7",proxies=F)

outcome_data_SBPMI<- extract_outcome_data(
  snps = exposure_data_SBP$SNP, outcomes = "ieu-a-798",proxies=F)

outcome_data_SBPILA<- extract_outcome_data(
  snps = exposure_data_SBP$SNP, outcomes = "ebi-a-GCST006907",proxies=F)

outcome_data_SBPIC<- extract_outcome_data(
  snps = exposure_data_SBP$SNP, outcomes = "ebi-a-GCST006910",proxies=F)

outcome_data_SBPISV<- extract_outcome_data(
  snps = exposure_data_SBP$SNP, outcomes = "ebi-a-GCST006909",proxies=F)

outcome_data_LDLCHD<- extract_outcome_data(
  snps = exposure_data_LDL$SNP, outcomes = "ieu-a-7",proxies=F)

outcome_data_LDLMI<- extract_outcome_data(
  snps = exposure_data_LDL$SNP, outcomes = "ieu-a-798",proxies=F)

outcome_data_LDLILA<- extract_outcome_data(
  snps = exposure_data_LDL$SNP, outcomes = "ebi-a-GCST006907",proxies=F)

outcome_data_LDLIC<- extract_outcome_data(
  snps = exposure_data_LDL$SNP, outcomes = "ebi-a-GCST006910",proxies=F)

outcome_data_LDLISV<- extract_outcome_data(
  snps = exposure_data_LDL$SNP, outcomes = "ebi-a-GCST006909",proxies=F)

outcome_data_HDLCHD<- extract_outcome_data(
  snps = exposure_data_HDL$SNP, outcomes = "ieu-a-7",proxies=F)

outcome_data_HDLMI<- extract_outcome_data(
  snps = exposure_data_HDL$SNP, outcomes = "ieu-a-798",proxies=F)

outcome_data_HDLILA<- extract_outcome_data(
  snps = exposure_data_HDL$SNP, outcomes = "ebi-a-GCST006907",proxies=F)

outcome_data_HDLIC<- extract_outcome_data(
  snps = exposure_data_HDL$SNP, outcomes = "ebi-a-GCST006910",proxies=F)

outcome_data_HDLISV<- extract_outcome_data(
  snps = exposure_data_HDL$SNP, outcomes = "ebi-a-GCST006909",proxies=F)

outcome_data_TGCHD<- extract_outcome_data(
  snps = exposure_data_TG$SNP, outcomes = "ieu-a-7",proxies=F)

outcome_data_TGMI<- extract_outcome_data(
  snps = exposure_data_TG$SNP, outcomes = "ieu-a-798",proxies=F)

outcome_data_TGILA<- extract_outcome_data(
  snps = exposure_data_TG$SNP, outcomes = "ebi-a-GCST006907",proxies=F)

outcome_data_TGIC<- extract_outcome_data(
  snps = exposure_data_TG$SNP, outcomes = "ebi-a-GCST006910",proxies=F)

outcome_data_TGISV<- extract_outcome_data(
  snps = exposure_data_TG$SNP, outcomes = "ebi-a-GCST006909",proxies=F)

MenCHD.dat <- harmonise_data(
  exposure_dat = exposure_data_Men, 
  outcome_dat = outcome_data_MenCHD
)

MenMI.dat <- harmonise_data(
  exposure_dat = exposure_data_Men, 
  outcome_dat = outcome_data_MenMI
)

MenIStkLA.dat <- harmonise_data(
  exposure_dat = exposure_data_Men, 
  outcome_dat = outcome_data_MenILA
)

MenIStkC.dat <- harmonise_data(
  exposure_dat = exposure_data_Men, 
  outcome_dat = outcome_data_MenIC
)

MenIStkSV.dat <- harmonise_data(
  exposure_dat = exposure_data_Men, 
  outcome_dat = outcome_data_MenISV
)

BMICHD.dat <- harmonise_data(
  exposure_dat = exposure_data_BMI, 
  outcome_dat = outcome_data_BMICHD
)

BMIMI.dat <- harmonise_data(
  exposure_dat = exposure_data_BMI, 
  outcome_dat = outcome_data_BMIMI
)

BMIIStkLA.dat <- harmonise_data(
  exposure_dat = exposure_data_BMI, 
  outcome_dat = outcome_data_BMIILA
)

BMIIStkC.dat <- harmonise_data(
  exposure_dat = exposure_data_BMI, 
  outcome_dat = outcome_data_BMIIC
)

BMIIStkSV.dat <- harmonise_data(
  exposure_dat = exposure_data_BMI, 
  outcome_dat = outcome_data_BMIISV
)

SBPCHD.dat <- harmonise_data(
  exposure_dat = exposure_data_SBP, 
  outcome_dat = outcome_data_SBPCHD
)

SBPMI.dat <- harmonise_data(
  exposure_dat = exposure_data_SBP, 
  outcome_dat = outcome_data_SBPMI
)

SBPIStkLA.dat <- harmonise_data(
  exposure_dat = exposure_data_SBP, 
  outcome_dat = outcome_data_SBPILA
)

SBPIStkC.dat <- harmonise_data(
  exposure_dat = exposure_data_SBP, 
  outcome_dat = outcome_data_SBPIC
)

SBPIStkSV.dat <- harmonise_data(
  exposure_dat = exposure_data_SBP, 
  outcome_dat = outcome_data_SBPISV
)

LDLCHD.dat <- harmonise_data(
  exposure_dat = exposure_data_LDL, 
  outcome_dat = outcome_data_LDLCHD
)

LDLMI.dat <- harmonise_data(
  exposure_dat = exposure_data_LDL, 
  outcome_dat = outcome_data_LDLMI
)

LDLIStkLA.dat <- harmonise_data(
  exposure_dat = exposure_data_LDL, 
  outcome_dat = outcome_data_LDLILA
)

LDLIStkC.dat <- harmonise_data(
  exposure_dat = exposure_data_LDL, 
  outcome_dat = outcome_data_LDLIC
)

LDLIStkSV.dat <- harmonise_data(
  exposure_dat = exposure_data_LDL, 
  outcome_dat = outcome_data_LDLISV
)

HDLCHD.dat <- harmonise_data(
  exposure_dat = exposure_data_HDL, 
  outcome_dat = outcome_data_HDLCHD
)

HDLMI.dat <- harmonise_data(
  exposure_dat = exposure_data_HDL, 
  outcome_dat = outcome_data_HDLMI
)

HDLIStkLA.dat <- harmonise_data(
  exposure_dat = exposure_data_HDL, 
  outcome_dat = outcome_data_HDLILA
)

HDLIStkC.dat <- harmonise_data(
  exposure_dat = exposure_data_HDL, 
  outcome_dat = outcome_data_HDLIC
)

HDLIStkSV.dat <- harmonise_data(
  exposure_dat = exposure_data_HDL, 
  outcome_dat = outcome_data_HDLISV
)

TGCHD.dat <- harmonise_data(
  exposure_dat = exposure_data_TG, 
  outcome_dat = outcome_data_TGCHD
)

TGMI.dat <- harmonise_data(
  exposure_dat = exposure_data_TG, 
  outcome_dat = outcome_data_TGMI
)

TGIStkLA.dat <- harmonise_data(
  exposure_dat = exposure_data_TG, 
  outcome_dat = outcome_data_TGILA
)

TGIStkC.dat <- harmonise_data(
  exposure_dat = exposure_data_TG, 
  outcome_dat = outcome_data_TGIC
)

TGIStkSV.dat <- harmonise_data(
  exposure_dat = exposure_data_TG, 
  outcome_dat = outcome_data_TGISV
)

#Results

mr_resultsMenCHD <- generate_odds_ratios(mr(MenCHD.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsMenMI <- generate_odds_ratios(mr(MenMI.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsMenILA <- generate_odds_ratios(mr(MenIStkLA.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsMenIC <- generate_odds_ratios(mr(MenIStkC.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsMenISV <- generate_odds_ratios(mr(MenIStkSV.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

mr_resultsBMICHD <- generate_odds_ratios(mr(BMICHD.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsBMIMI <- generate_odds_ratios(mr(BMIMI.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsBMIILA <- generate_odds_ratios(mr(BMIIStkLA.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsBMIIC <- generate_odds_ratios(mr(BMIIStkC.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsBMIISV <- generate_odds_ratios(mr(BMIIStkSV.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

mr_resultsSBPCHD <- generate_odds_ratios(mr(SBPCHD.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsSBPMI <- generate_odds_ratios(mr(SBPMI.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsSBPILA <- generate_odds_ratios(mr(SBPIStkLA.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsSBPIC <- generate_odds_ratios(mr(SBPIStkC.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsSBPISV <- generate_odds_ratios(mr(SBPIStkSV.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

mr_resultsLDLCHD <- generate_odds_ratios(mr(LDLCHD.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsLDLMI <- generate_odds_ratios(mr(LDLMI.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsLDLILA <- generate_odds_ratios(mr(LDLIStkLA.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsLDLIC <- generate_odds_ratios(mr(LDLIStkC.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsLDLISV <- generate_odds_ratios(mr(LDLIStkSV.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

mr_resultsHDLCHD <- generate_odds_ratios(mr(HDLCHD.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsHDLMI <- generate_odds_ratios(mr(HDLMI.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsHDLILA <- generate_odds_ratios(mr(HDLIStkLA.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsHDLIC <- generate_odds_ratios(mr(HDLIStkC.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsHDLISV <- generate_odds_ratios(mr(HDLIStkSV.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

mr_resultsTGCHD <- generate_odds_ratios(mr(TGCHD.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsTGMI <- generate_odds_ratios(mr(TGMI.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsTGILA <- generate_odds_ratios(mr(TGIStkLA.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsTGIC <- generate_odds_ratios(mr(TGIStkC.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsTGISV <- generate_odds_ratios(mr(TGIStkSV.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

#Mean F-statistics and mean SNP numbers

F1<- MenCHD.dat$beta.exposure^2/MenCHD.dat$se.exposure^2
F2<- MenMI.dat$beta.exposure^2/MenMI.dat$se.exposure^2
F3<- MenIStkLA.dat$beta.exposure^2/MenIStkLA.dat$se.exposure^2
F4<- MenIStkC.dat$beta.exposure^2/MenIStkC.dat$se.exposure^2
F5<- MenIStkSV.dat$beta.exposure^2/MenIStkSV.dat$se.exposure^2
F6<- BMICHD.dat$beta.exposure^2/BMICHD.dat$se.exposure^2
F7<- BMIMI.dat$beta.exposure^2/BMIMI.dat$se.exposure^2
F8<- BMIIStkLA.dat$beta.exposure^2/BMIIStkLA.dat$se.exposure^2
F9<- BMIIStkC.dat$beta.exposure^2/BMIIStkC.dat$se.exposure^2
F10<- BMIIStkSV.dat$beta.exposure^2/BMIIStkSV.dat$se.exposure^2
F11<- SBPCHD.dat$beta.exposure^2/SBPCHD.dat$se.exposure^2
F12<- SBPMI.dat$beta.exposure^2/SBPMI.dat$se.exposure^2
F13<- SBPIStkLA.dat$beta.exposure^2/SBPIStkLA.dat$se.exposure^2
F14<- SBPIStkC.dat$beta.exposure^2/SBPIStkC.dat$se.exposure^2
F15<- SBPIStkSV.dat$beta.exposure^2/SBPIStkSV.dat$se.exposure^2
F16<- LDLCHD.dat$beta.exposure^2/LDLCHD.dat$se.exposure^2
F17<- LDLMI.dat$beta.exposure^2/LDLMI.dat$se.exposure^2
F18<- LDLIStkLA.dat$beta.exposure^2/LDLIStkLA.dat$se.exposure^2
F19<- LDLIStkC.dat$beta.exposure^2/LDLIStkC.dat$se.exposure^2
F20<- LDLIStkSV.dat$beta.exposure^2/LDLIStkSV.dat$se.exposure^2
F21<- HDLCHD.dat$beta.exposure^2/HDLCHD.dat$se.exposure^2
F22<- HDLMI.dat$beta.exposure^2/HDLMI.dat$se.exposure^2
F23<- HDLIStkLA.dat$beta.exposure^2/HDLIStkLA.dat$se.exposure^2
F24<- HDLIStkC.dat$beta.exposure^2/HDLIStkC.dat$se.exposure^2
F25<- HDLIStkSV.dat$beta.exposure^2/HDLIStkSV.dat$se.exposure^2
F26<- TGCHD.dat$beta.exposure^2/TGCHD.dat$se.exposure^2
F27<- TGMI.dat$beta.exposure^2/TGMI.dat$se.exposure^2
F28<- TGIStkLA.dat$beta.exposure^2/TGIStkLA.dat$se.exposure^2
F29<- TGIStkC.dat$beta.exposure^2/TGIStkC.dat$se.exposure^2
F30<- TGIStkSV.dat$beta.exposure^2/TGIStkSV.dat$se.exposure^2
F_stats1<-c(mean(F1),mean(F2),mean(F3),mean(F4),mean(F5),
            mean(F6),mean(F7),mean(F8),mean(F9),mean(F10),
            mean(F11),mean(F12),mean(F13),mean(F14),mean(F15),
            mean(F16),mean(F17),mean(F18),mean(F19),mean(F20),
            mean(F21),mean(F22),mean(F23),mean(F24),mean(F25),
            mean(F26),mean(F27),mean(F28),mean(F29),mean(F30))
min(F_stats1)
max(F_stats1)

SNPnumMen<-mean(mr_resultsMenCHD$nsnp[1],
                mr_resultsMenMI$nsnp[1],
                mr_resultsMenILA$nsnp[1],
                mr_resultsMenIC$nsnp[1],
                mr_resultsMenISV$nsnp[1])
SNPnumBMI<-mean(mr_resultsBMICHD$nsnp[1],
                mr_resultsBMIMI$nsnp[1],
                mr_resultsBMIILA$nsnp[1],
                mr_resultsBMIIC$nsnp[1],
                mr_resultsBMIISV$nsnp[1])
SNPnumSBP<-mean(mr_resultsSBPCHD$nsnp[1],
                mr_resultsSBPMI$nsnp[1],
                mr_resultsSBPILA$nsnp[1],
                mr_resultsSBPIC$nsnp[1],
                mr_resultsSBPISV$nsnp[1])
SNPnumLDL<-mean(mr_resultsLDLCHD$nsnp[1],
                mr_resultsLDLMI$nsnp[1],
                mr_resultsLDLILA$nsnp[1],
                mr_resultsLDLIC$nsnp[1],
                mr_resultsLDLISV$nsnp[1])
SNPnumHDL<-mean(mr_resultsHDLCHD$nsnp[1],
                mr_resultsHDLMI$nsnp[1],
                mr_resultsHDLILA$nsnp[1],
                mr_resultsHDLIC$nsnp[1],
                mr_resultsHDLISV$nsnp[1])
SNPnumTG<-mean(mr_resultsTGCHD$nsnp[1],
                mr_resultsTGMI$nsnp[1],
                mr_resultsTGILA$nsnp[1],
                mr_resultsTGIC$nsnp[1],
                mr_resultsTGISV$nsnp[1])
sum(SNPnumMen,SNPnumBMI,SNPnumSBP,SNPnumLDL,SNPnumHDL,SNPnumTG)

# Multiple testing correction for IVW

#BF_threshold<-0.05/6
#BF_threshold
#mr_resultsHDLILA$pval[1]

# Plotting univariable estimates

res_data<-data.frame(c(mr_resultsMenCHD$or,mr_resultsMenILA$or,mr_resultsMenIC$or,mr_resultsMenISV$or,
                       mr_resultsBMICHD$or,mr_resultsBMIILA$or,mr_resultsBMIIC$or,mr_resultsBMIISV$or,
                       mr_resultsSBPCHD$or,mr_resultsSBPILA$or,mr_resultsSBPIC$or,mr_resultsSBPISV$or,
                       mr_resultsLDLCHD$or,mr_resultsLDLILA$or,mr_resultsLDLIC$or,mr_resultsLDLISV$or,
                       mr_resultsHDLCHD$or,mr_resultsHDLILA$or,mr_resultsHDLIC$or,mr_resultsHDLISV$or,
                       mr_resultsTGCHD$or,mr_resultsTGILA$or,mr_resultsTGIC$or,mr_resultsTGISV$or),
                     c(mr_resultsMenCHD$or_lci95,mr_resultsMenILA$or_lci95,mr_resultsMenIC$or_lci95,mr_resultsMenISV$or_lci95,
                       mr_resultsBMICHD$or_lci95,mr_resultsBMIILA$or_lci95,mr_resultsBMIIC$or_lci95,mr_resultsBMIISV$or_lci95,
                       mr_resultsSBPCHD$or_lci95,mr_resultsSBPILA$or_lci95,mr_resultsSBPIC$or_lci95,mr_resultsSBPISV$or_lci95,
                       mr_resultsLDLCHD$or_lci95,mr_resultsLDLILA$or_lci95,mr_resultsLDLIC$or_lci95,mr_resultsLDLISV$or_lci95,
                       mr_resultsHDLCHD$or_lci95,mr_resultsHDLILA$or_lci95,mr_resultsHDLIC$or_lci95,mr_resultsHDLISV$or_lci95,
                       mr_resultsTGCHD$or_lci95,mr_resultsTGILA$or_lci95,mr_resultsTGIC$or_lci95,mr_resultsTGISV$or_lci95),
                     c(mr_resultsMenCHD$or_uci95,mr_resultsMenILA$or_uci95,mr_resultsMenIC$or_uci95,mr_resultsMenISV$or_uci95,
                       mr_resultsBMICHD$or_uci95,mr_resultsBMIILA$or_uci95,mr_resultsBMIIC$or_uci95,mr_resultsBMIISV$or_uci95,
                       mr_resultsSBPCHD$or_uci95,mr_resultsSBPILA$or_uci95,mr_resultsSBPIC$or_uci95,mr_resultsSBPISV$or_uci95,
                       mr_resultsLDLCHD$or_uci95,mr_resultsLDLILA$or_uci95,mr_resultsLDLIC$or_uci95,mr_resultsLDLISV$or_uci95,
                       mr_resultsHDLCHD$or_uci95,mr_resultsHDLILA$or_uci95,mr_resultsHDLIC$or_uci95,mr_resultsHDLISV$or_uci95,
                       mr_resultsTGCHD$or_uci95,mr_resultsTGILA$or_uci95,mr_resultsTGIC$or_uci95,mr_resultsTGISV$or_uci95),
                     rep(c("IVW","MR-Egger","Weighted Median","Weighted Mode"),24),
                     rep(c(rep("CHD",4),rep("I.Stroke (large-artery)",4),rep("I.Stroke (cardioembolic)",4),rep("I.Stroke (small vessel)",4)),6),
                     c(rep("menarche",16),rep("BMI",16),rep("SBP",16),rep("LDL",16),rep("HDL",16),rep("TG",16)))

names(res_data)<-c("or","orlci_95","oruci_95","method","outcome","exposure")

resdat_ivw<-res_data[seq(1,96,4),]

resdat_ivw<-resdat_ivw[c(rev(seq(4,24,4)),rev(seq(3,24,4)),rev(seq(2,24,4)),rev(seq(1,24,4))),]

resdat_ivw$index<-c(1,1.2,1.4,1.6,1.8,2,
                    3,3.2,3.4,3.6,3.8,4,
                    5,5.2,5.4,5.6,5.8,6,
                    7,7.2,7.4,7.6,7.8,8)

uniplot1<-ggplot(resdat_ivw)+geom_point(aes(x=index, y=or))+
  coord_flip()+geom_errorbar(aes(ymin=orlci_95,ymax=oruci_95,x=index),width=0,orientation = "x")+
  scale_y_continuous(limits = c(-2,4),name = "Effect (odds ratio)")+geom_hline(yintercept=1,linetype = "dashed")+theme_bw()+
  scale_x_continuous(limits = c(1,9),breaks=c(2.5,4.5,6.5,8.5),labels= resdat_ivw$outcome[c(1,7,17,24)],name = "Outcome")+
  scale_color_discrete(name = "Exposure")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png("Figure3_1.png",width=6,height=6,units="in",res=200)
uniplot1
dev.off()

#heterogeneity

hetMenCHD <- mr_heterogeneity(MenCHD.dat)
hetMenILA <- mr_heterogeneity(MenIStkLA.dat)
hetMenIC <- mr_heterogeneity(MenIStkC.dat)
hetMenISV <- mr_heterogeneity(MenIStkSV.dat)

hetBMICHD <- mr_heterogeneity(BMICHD.dat)
hetBMIILA <- mr_heterogeneity(BMIIStkLA.dat)
hetBMIIC <- mr_heterogeneity(BMIIStkC.dat)
hetBMIISV <- mr_heterogeneity(BMIIStkSV.dat)

hetSBPCHD <- mr_heterogeneity(SBPCHD.dat)
hetSBPILA <- mr_heterogeneity(SBPIStkLA.dat)
hetSBPIC <-  mr_heterogeneity(SBPIStkC.dat)
hetSBPISV <- mr_heterogeneity(SBPIStkSV.dat)

hetLDLCHD <- mr_heterogeneity(LDLCHD.dat)
hetLDLILA <- mr_heterogeneity(LDLIStkLA.dat)
hetLDLIC <- mr_heterogeneity(LDLIStkC.dat)
hetLDLISV <- mr_heterogeneity(LDLIStkSV.dat)

hetHDLCHD <- mr_heterogeneity(HDLCHD.dat)
hetHDLILA <- mr_heterogeneity(HDLIStkLA.dat)
hetHDLIC <- mr_heterogeneity(HDLIStkC.dat)
hetHDLISV <- mr_heterogeneity(HDLIStkSV.dat)

hetTGCHD <- mr_heterogeneity(TGCHD.dat)
hetTGILA <- mr_heterogeneity(TGIStkLA.dat)
hetTGIC <- mr_heterogeneity(TGIStkC.dat)
hetTGISV <- mr_heterogeneity(TGIStkSV.dat)

############
### MVMR ###
############

#Extract all exposure data

exp_list<-c("ukb-a-315","ukb-a-248","ukb-a-360","ieu-a-300","ieu-a-299","ieu-a-302")

output<-NULL

for(i in 1:6){
  
  reference<-exp_list[i]
  
  r.exposures<-exp_list[-i]
  
  A<-extract_instruments(reference,clump = F)
  
  for(j in 1:5){
    
    B<-extract_instruments(r.exposures[j],clump = F)
    
    SNPlist<-c(A$SNP,B$SNP)
    SNPlist<-unique(SNPlist)
    
    splitting_function <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
    
    SNP_chunks<-splitting_function(SNPlist,20)
    
    for(k in 1:20){
      refassocs <- extract_outcome_data(
        snps = unlist(SNP_chunks[k]), outcomes = reference,proxies=F)
      
      r.expassocs <- extract_outcome_data(
        snps = unlist(SNP_chunks[k]), outcomes = r.exposures[j],proxies=F)
      
      if(k == 1){
        refcomb1<-refassocs
        r.expcomb1<-r.expassocs
      }else{
        refcomb1<-rbind(refcomb1,refassocs)
        r.expcomb1<-rbind(r.expcomb1,r.expassocs)
      }
    }
    
    SNPlistA<-refcomb1$SNP[refcomb1$SNP %in% r.expcomb1$SNP]
    SNPlistB<-r.expcomb1$SNP[r.expcomb1$SNP %in% refcomb1$SNP]
    SNPlistcomb<-unique(c(SNPlistA,SNPlistB))
    
    SNP_chunks<-splitting_function(SNPlistcomb,20)
    
    for(k in 1:20){
      refassocs <- extract_outcome_data(
        snps = unlist(SNP_chunks[k]), outcomes = reference,proxies=F)
      
      r.expassocs <- extract_outcome_data(
        snps = unlist(SNP_chunks[k]), outcomes = r.exposures[j],proxies=F)
      
      if(k == 1){
        refcomb1<-refassocs
        r.expcomb1<-r.expassocs
      }else{
        refcomb1<-rbind(refcomb1,refassocs)
        r.expcomb1<-rbind(r.expcomb1,r.expassocs)
      }
    }
    
    #Reformatting Menarche data for processing (clumping SNPs etc)
    Clump.datA<-format_data(
      refcomb1,
      type = "exposure",
      effect_allele_col = "effect_allele.outcome",
      beta_col = "beta.outcome",
      se_col = "se.outcome",
      other_allele_col = "other_allele.outcome",
      pval_col = "pval.outcome",
      eaf_col = "eaf.outcome",
      phenotype_col = "outcome",
      header = T
    )
    
    A_clump<-clump_data(Clump.datA, clump_r2 = 0.01)
    
    #Reformatting BMI data for processing (clumping SNPs etc)
    Clump.datB<-format_data(
      r.expcomb1,
      type = "exposure",
      effect_allele_col = "effect_allele.outcome",
      beta_col = "beta.outcome",
      se_col = "se.outcome",
      other_allele_col = "other_allele.outcome",
      pval_col = "pval.outcome",
      eaf_col = "eaf.outcome",
      phenotype_col = "outcome",
      header = T
    )
    
    B_clump<-clump_data(Clump.datB, clump_r2 = 0.01)
    
    if(nrow(A_clump) < nrow(B_clump)){
      
      B_clump<-Clump.datB[Clump.datB$SNP %in% A_clump$SNP,]
      
    }else{
      
      A_clump<-Clump.datA[Clump.datA$SNP %in% B_clump$SNP,]
      
    }
    
    # Obtain CHD outcome associations
    Outcome.data<- extract_outcome_data(
      snps = B_clump$SNP, outcomes = "ieu-a-7",proxies=F)
    
    # Harmonize exposure data
    HarmA_data <- harmonise_data(
      exposure_dat = A_clump, 
      outcome_dat = Outcome.data
    )
    
    HarmB_data <- harmonise_data(
      exposure_dat = B_clump, 
      outcome_dat = Outcome.data
    )
    
    HarmA_data<-HarmA_data[HarmA_data$mr_keep == T,]
    
    HarmB_data<-HarmB_data[HarmB_data$mr_keep == T,]
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmB_data$SNP,]
    
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmA_data$SNP,]
    
    HarmA_data<-HarmA_data[order(HarmA_data$SNP),]
    
    HarmB_data<-HarmB_data[order(HarmB_data$SNP),]
    
    #Create combined dataframe for CHD MVMR
    tBXGs<-matrix(c(HarmA_data$beta.exposure,HarmB_data$beta.exposure),ncol=2)
    
    tseXGs<-matrix(c(HarmA_data$se.exposure,HarmA_data$se.exposure),ncol=2)
    
    tMVMR.dat<-format_mvmr(tBXGs,HarmA_data$beta.outcome,tseXGs,HarmA_data$se.outcome,HarmA_data$SNP)
    
    #Estimate conditional F-statistics for CHD
    print("CF stats for")
    print(reference)
    print("and")
    print(r.exposures[j])
    CF_stat<-strength_mvmr(tMVMR.dat)
    
    if(is.null(output)){
      output<-data.frame(reference,r.exposures[j],CF_stat)
    }else{
      output<-rbind(output,data.frame(reference,r.exposures[j],CF_stat))
    }
    
  }
  
}

write.csv(output,"F-outs.csv",row.names=F)

F_matrix<-matrix(rep(0,36),ncol=6)
F_matrix[1,]<-c(0,output$exposure1[1:5])
F_matrix[2,]<-c(output$exposure1[6],0,output$exposure1[7:10])
F_matrix[3,]<-c(output$exposure1[11:12],0,output$exposure1[13:15])
F_matrix[4,]<-c(output$exposure1[16:18],0,output$exposure1[19:20])
F_matrix[5,]<-c(output$exposure1[21:24],0,output$exposure1[25])
F_matrix[6,]<-c(output$exposure1[26:30],0)

meanF_matrix<-matrix(rep(0,36),ncol=6)

for(i in 1:6){
  for(j in 1:6){
    meanF_matrix[i,j]<-mean(c(F_matrix[i,j],F_matrix[j,i]))
  }
}

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

melted_fdat <- na.omit(melt(get_upper_tri(meanF_matrix)))

Fstat_map <- ggplot(melted_fdat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(0,100), space = "Lab", 
                       name="Conditional\nF-Statistic") +
  scale_x_continuous(breaks = 1:6,labels = c("Menarche","BMI","SBP","LDL","HDL","TG"))+
  scale_y_continuous(breaks = 1:6,labels = c("Menarche","BMI","SBP","LDL","HDL","TG"))+
  xlab("MVMR Exposure 1")+ylab("MVMR Exposure 2")+
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 10, hjust = 1),
        axis.text.y = element_text(angle = 45, vjust = -1.5,
                                   size = 10, hjust = 1),
        panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  coord_fixed()+
  geom_text(aes(Var2, Var1, label = round(value,digits = 2)), color = "black", size = 3)

png("Figure4_1.png",width=6,height=6,units="in",res=200)
Fstat_map
dev.off()

#LDL and TG are strongest conditional pair - assessing additional

exp_list<-c("ukb-a-315","ukb-a-248","ukb-a-360","ieu-a-299")

output2<-NULL

reference1<-"ieu-a-300"

reference2<-"ieu-a-302"

r.exposures<-exp_list

A<-extract_instruments(reference1,clump = F)

B<-extract_instruments(reference2,clump = F)

for(i in 1:4){
  
  C<-extract_instruments(r.exposures[i],clump = F)
  
  SNPlist<-c(A$SNP,B$SNP,C$SNP)
  SNPlist<-unique(SNPlist)
  
  splitting_function <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
  
  SNP_chunks<-splitting_function(SNPlist,20)
  
  for(k in 1:20){
    ref1assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = reference1,proxies=F)
    
    ref2assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = reference2,proxies=F)
    
    r.expassocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = r.exposures[i],proxies=F)
    
    if(k == 1){
      refcomb1<-ref1assocs
      refcomb2<-ref2assocs
      r.expcomb1<-r.expassocs
    }else{
      refcomb1<-rbind(refcomb1,ref1assocs)
      refcomb2<-rbind(refcomb2,ref2assocs)
      r.expcomb1<-rbind(r.expcomb1,r.expassocs)
    }
  }
  
  SNPlistA<-refcomb1$SNP[refcomb1$SNP %in% r.expcomb1$SNP]
  SNPlistB<-r.expcomb1$SNP[r.expcomb1$SNP %in% SNPlistA]
  SNPlistcomb<-unique(c(SNPlistA,SNPlistB))
  
  SNP_chunks<-splitting_function(SNPlistcomb,20)
  
  for(k in 1:20){
    ref1assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = reference1,proxies=F)
    
    ref2assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = reference2,proxies=F)
    
    r.expassocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = r.exposures[i],proxies=F)
    
    if(k == 1){
      refcomb1<-ref1assocs
      refcomb2<-ref2assocs
      r.expcomb1<-r.expassocs
    }else{
      refcomb1<-rbind(refcomb1,ref1assocs)
      refcomb2<-rbind(refcomb2,ref2assocs)
      r.expcomb1<-rbind(r.expcomb1,r.expassocs)
    }
  }
  
  #Reformatting Menarche data for processing (clumping SNPs etc)
  Clump.datA<-format_data(
    refcomb1,
    type = "exposure",
    effect_allele_col = "effect_allele.outcome",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    eaf_col = "eaf.outcome",
    phenotype_col = "outcome",
    header = T
  )
  
  A_clump<-clump_data(Clump.datA, clump_r2 = 0.01)
  
  Clump.datB<-format_data(
    refcomb2,
    type = "exposure",
    effect_allele_col = "effect_allele.outcome",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    eaf_col = "eaf.outcome",
    phenotype_col = "outcome",
    header = T
  )
  
  B_clump<-clump_data(Clump.datB, clump_r2 = 0.01)
  
  Clump.datC<-format_data(
    r.expcomb1,
    type = "exposure",
    effect_allele_col = "effect_allele.outcome",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    eaf_col = "eaf.outcome",
    phenotype_col = "outcome",
    header = T
  )
  
  C_clump<-clump_data(Clump.datC, clump_r2 = 0.01)
  
  if(nrow(A_clump) < nrow(B_clump)){
    
    B_clump<-Clump.datB[Clump.datB$SNP %in% A_clump$SNP,]
    
  }else{
    
    A_clump<-Clump.datA[Clump.datA$SNP %in% B_clump$SNP,]
    
  }
  
  if(nrow(C_clump) < nrow(B_clump)){
    
    A_clump<-Clump.datA[Clump.datA$SNP %in% C_clump$SNP,]
    B_clump<-Clump.datB[Clump.datB$SNP %in% C_clump$SNP,]
    
  }else{
    
    C_clump<-Clump.datC[Clump.datC$SNP %in% A_clump$SNP,]
    
  }
  
  # Obtain CHD outcome associations
  Outcome.data<- extract_outcome_data(
    snps = A_clump$SNP, outcomes = "ieu-a-7",proxies=F)
  
  # Harmonize exposure data
  HarmA_data <- harmonise_data(
    exposure_dat = A_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmB_data <- harmonise_data(
    exposure_dat = B_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmC_data <- harmonise_data(
    exposure_dat = C_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmA_data<-HarmA_data[HarmA_data$mr_keep == T,]
  
  HarmB_data<-HarmB_data[HarmB_data$mr_keep == T,]
  
  HarmC_data<-HarmC_data[HarmC_data$mr_keep == T,]
  
  HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmC_data$SNP,]
  
  HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmC_data$SNP,]
  
  HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmA_data$SNP,]
  
  HarmA_data<-HarmA_data[order(HarmA_data$SNP),]
  
  HarmB_data<-HarmB_data[order(HarmB_data$SNP),]
  
  HarmC_data<-HarmC_data[order(HarmC_data$SNP),]
  
  #Create combined dataframe for CHD MVMR
  tBXGs<-matrix(c(HarmA_data$beta.exposure,HarmB_data$beta.exposure,HarmC_data$beta.exposure),ncol=3)
  
  tseXGs<-matrix(c(HarmA_data$se.exposure,HarmA_data$se.exposure,HarmC_data$se.exposure),ncol=3)
  
  tMVMR.dat<-format_mvmr(tBXGs,HarmA_data$beta.outcome,tseXGs,HarmA_data$se.outcome,HarmA_data$SNP)
  
  #Estimate conditional F-statistics for CHD
  print("CF stats for")
  print(reference1)
  print("and")
  print(reference2)
  print("on")
  print(r.exposures[i])
  CF_stat<-strength_mvmr(tMVMR.dat)
  
  if(is.null(output2)){
    output2<-data.frame(reference1,reference2,r.exposures[i],CF_stat)
  }else{
    output2<-rbind(output2,data.frame(reference1,reference2,r.exposures[i],CF_stat))
  }
  
}

write.csv(output2,"F-out2.csv",row.names=F)

#LDL TG, and HDL are strongest conditional pair - assessing additional exposures

exp_list<-c("ukb-a-315","ukb-a-360","ukb-a-248")

output3<-NULL

reference1<-"ieu-a-300"

reference2<-"ieu-a-302"

reference3<-"ukb-a-299"

r.exposures<-exp_list

A<-extract_instruments(reference1,clump = F)

B<-extract_instruments(reference2,clump = F)

C<-extract_instruments(reference3,clump = F)

for(i in 1:3){
  
  D<-extract_instruments(r.exposures[i],clump = F)
  
  SNPlist<-c(A$SNP,B$SNP,C$SNP,D$SNP)
  SNPlist<-unique(SNPlist)
  
  splitting_function <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
  
  SNP_chunks<-splitting_function(SNPlist,20)
  
  for(k in 1:20){
    ref1assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = reference1,proxies=F)
    
    ref2assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = reference2,proxies=F)
    
    ref3assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = reference3,proxies=F)
    
    r.expassocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = r.exposures[i],proxies=F)
    
    if(k == 1){
      refcomb1<-ref1assocs
      refcomb2<-ref2assocs
      refcomb3<-ref3assocs
      r.expcomb1<-r.expassocs
    }else{
      refcomb1<-rbind(refcomb1,ref1assocs)
      refcomb2<-rbind(refcomb2,ref2assocs)
      refcomb3<-rbind(refcomb3,ref3assocs)
      r.expcomb1<-rbind(r.expcomb1,r.expassocs)
    }
  }
  
  SNPlistA<-refcomb1$SNP[refcomb1$SNP %in% r.expcomb1$SNP]
  SNPlistB<-refcomb2$SNP[refcomb1$SNP %in% SNPlistA]
  SNPlistC<-refcomb3$SNP[refcomb3$SNP %in% SNPlistB]
  SNPlistD<-r.expcomb1$SNP[r.expcomb1$SNP %in% SNPlistC]
  SNPlistcomb<-unique(c(SNPlistA,SNPlistB,SNPlistC,SNPlistD))
  
  SNP_chunks<-splitting_function(SNPlistcomb,20)
  
  for(k in 1:20){
    ref1assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = reference1,proxies=F)
    
    ref2assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = reference2,proxies=F)
    
    ref3assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = reference3,proxies=F)
    
    r.expassocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = r.exposures[i],proxies=F)
    
    if(k == 1){
      refcomb1<-ref1assocs
      refcomb2<-ref2assocs
      refcomb3<-ref3assocs
      r.expcomb1<-r.expassocs
    }else{
      refcomb1<-rbind(refcomb1,ref1assocs)
      refcomb2<-rbind(refcomb2,ref2assocs)
      refcomb3<-rbind(refcomb3,ref3assocs)
      r.expcomb1<-rbind(r.expcomb1,r.expassocs)
    }
  }
  
  #Reformatting Menarche data for processing (clumping SNPs etc)
  Clump.datA<-format_data(
    refcomb1,
    type = "exposure",
    effect_allele_col = "effect_allele.outcome",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    eaf_col = "eaf.outcome",
    phenotype_col = "outcome",
    header = T
  )
  
  A_clump<-clump_data(Clump.datA, clump_r2 = 0.01)
  
  Clump.datB<-format_data(
    refcomb2,
    type = "exposure",
    effect_allele_col = "effect_allele.outcome",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    eaf_col = "eaf.outcome",
    phenotype_col = "outcome",
    header = T
  )
  
  B_clump<-clump_data(Clump.datB, clump_r2 = 0.01)
  
  Clump.datC<-format_data(
    refcomb3,
    type = "exposure",
    effect_allele_col = "effect_allele.outcome",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    eaf_col = "eaf.outcome",
    phenotype_col = "outcome",
    header = T
  )
  
  C_clump<-clump_data(Clump.datC, clump_r2 = 0.01)
  
  Clump.datD<-format_data(
    r.expcomb1,
    type = "exposure",
    effect_allele_col = "effect_allele.outcome",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    eaf_col = "eaf.outcome",
    phenotype_col = "outcome",
    header = T
  )
  
  D_clump<-clump_data(Clump.datD, clump_r2 = 0.01)
  
  if(nrow(A_clump) > nrow(B_clump)){
    
    B_clump<-Clump.datB[Clump.datB$SNP %in% A_clump$SNP,]
    
  }else{
    
    A_clump<-Clump.datA[Clump.datA$SNP %in% B_clump$SNP,]
    
  }
  
  if(nrow(C_clump) > nrow(B_clump)){
    
    A_clump<-Clump.datA[Clump.datA$SNP %in% C_clump$SNP,]
    B_clump<-Clump.datB[Clump.datB$SNP %in% C_clump$SNP,]
    
  }else{
    
    C_clump<-Clump.datC[Clump.datC$SNP %in% A_clump$SNP,]
    
  }
  
  if(nrow(D_clump) > nrow(C_clump)){
    
    A_clump<-Clump.datA[Clump.datA$SNP %in% D_clump$SNP,]
    B_clump<-Clump.datB[Clump.datB$SNP %in% D_clump$SNP,]
    C_clump<-Clump.datC[Clump.datC$SNP %in% D_clump$SNP,]
    
  }else{
    
    D_clump<-Clump.datD[Clump.datD$SNP %in% A_clump$SNP,]
    
  }
  
  # Obtain CHD outcome associations
  Outcome.data<- extract_outcome_data(
    snps = A_clump$SNP, outcomes = "ieu-a-7",proxies=F)
  
  # Harmonize exposure data
  HarmA_data <- harmonise_data(
    exposure_dat = A_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmB_data <- harmonise_data(
    exposure_dat = B_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmC_data <- harmonise_data(
    exposure_dat = C_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmD_data <- harmonise_data(
    exposure_dat = D_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmA_data<-HarmA_data[HarmA_data$mr_keep == T,]
  
  HarmB_data<-HarmB_data[HarmB_data$mr_keep == T,]
  
  HarmC_data<-HarmC_data[HarmC_data$mr_keep == T,]
  
  HarmD_data<-HarmD_data[HarmD_data$mr_keep == T,]
  
  
  if(nrow(HarmA_data) < nrow(HarmB_data)){
    
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmA_data$SNP,]
    
  }else{
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmB_data$SNP,]
    
  }
  
  if(nrow(HarmC_data) < nrow(HarmB_data)){
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmC_data$SNP,]
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmC_data$SNP,]
    
  }else{
    
    HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmA_data$SNP,]
    
  }
  
  if(nrow(HarmD_data) < nrow(HarmC_data)){
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmD_data$SNP,]
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmD_data$SNP,]
    HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmD_data$SNP,]
    
  }else{
    
    HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmA_data$SNP,]
    
  }
  
  HarmA_data<-HarmA_data[order(HarmA_data$SNP),]
  
  HarmB_data<-HarmB_data[order(HarmB_data$SNP),]
  
  HarmC_data<-HarmC_data[order(HarmC_data$SNP),]
  
  HarmD_data<-HarmD_data[order(HarmD_data$SNP),]
  
  #Create combined dataframe for CHD MVMR
  tBXGs<-matrix(c(HarmA_data$beta.exposure,HarmB_data$beta.exposure,HarmC_data$beta.exposure,HarmD_data$beta.exposure),ncol=4)
  
  tseXGs<-matrix(c(HarmA_data$se.exposure,HarmA_data$se.exposure,HarmC_data$se.exposure,HarmD_data$se.exposure),ncol=4)
  
  tMVMR.dat<-format_mvmr(tBXGs,HarmA_data$beta.outcome,tseXGs,HarmA_data$se.outcome,HarmA_data$SNP)
  
  #Estimate conditional F-statistics for CHD
  print("CF stats for")
  print(reference1)
  print("and")
  print(reference2)
  print("and")
  print(reference3)
  print("on")
  print(r.exposures[i])
  CF_stat<-strength_mvmr(tMVMR.dat)
  
  if(is.null(output3)){
    output3<-data.frame(reference1,reference2,reference3,r.exposures[i],CF_stat)
  }else{
    output3<-rbind(output3,data.frame(reference1,reference2,reference3,r.exposures[i],CF_stat))
  }
  
}

write.csv(output3,"F-out3.csv",row.names=F)

#LDL, HDL, TG, and SBP are strongest - assessing additional exposures

exp_list<-c("ukb-a-315","ukb-a-248")

output4<-NULL

reference1<-"ieu-a-300"

reference2<-"ieu-a-302"

reference3<-"ukb-a-360"

reference4<-"ieu-a-299"

r.exposures<-exp_list

A<-extract_instruments(reference1,clump = F)

B<-extract_instruments(reference2,clump = F)

C<-extract_instruments(reference3,clump = F)

D<-extract_instruments(reference4,clump = F)

for(i in 1:2){
  
  E<-extract_instruments(r.exposures[i],clump = F)
  
  SNPlist<-c(A$SNP,B$SNP,C$SNP,D$SNP,E$SNP)
  SNPlist<-unique(SNPlist)
  
  splitting_function <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
  
  SNP_chunks<-splitting_function(SNPlist,20)
  
  for(k in 1:20){
    ref1assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = reference1,proxies=F)
    
    ref2assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = reference2,proxies=F)
    
    ref3assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = reference3,proxies=F)
    
    ref4assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = reference4,proxies=F)
    
    r.expassocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = r.exposures[i],proxies=F)
    
    if(k == 1){
      refcomb1<-ref1assocs
      refcomb2<-ref2assocs
      refcomb3<-ref3assocs
      refcomb4<-ref4assocs
      r.expcomb1<-r.expassocs
    }else{
      refcomb1<-rbind(refcomb1,ref1assocs)
      refcomb2<-rbind(refcomb2,ref2assocs)
      refcomb3<-rbind(refcomb3,ref3assocs)
      refcomb4<-rbind(refcomb4,ref4assocs)
      r.expcomb1<-rbind(r.expcomb1,r.expassocs)
    }
  }
  
  SNPlistA<-refcomb1$SNP[refcomb1$SNP %in% r.expcomb1$SNP]
  SNPlistB<-refcomb2$SNP[refcomb1$SNP %in% SNPlistA]
  SNPlistC<-refcomb3$SNP[refcomb3$SNP %in% SNPlistB]
  SNPlistD<-refcomb4$SNP[refcomb4$SNP %in% SNPlistC]
  SNPlistE<-r.expcomb1$SNP[r.expcomb1$SNP %in% SNPlistD]
  SNPlistcomb<-unique(c(SNPlistA,SNPlistB,SNPlistC,SNPlistD,SNPlistE))
  
  SNP_chunks<-splitting_function(SNPlistcomb,20)
  
  for(k in 1:20){
    ref1assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = reference1,proxies=F)
    
    ref2assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = reference2,proxies=F)
    
    ref3assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = reference3,proxies=F)
    
    ref4assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = reference4,proxies=F)
    
    r.expassocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = r.exposures[i],proxies=F)
    
    if(k == 1){
      refcomb1<-ref1assocs
      refcomb2<-ref2assocs
      refcomb3<-ref3assocs
      refcomb4<-ref4assocs
      r.expcomb1<-r.expassocs
    }else{
      refcomb1<-rbind(refcomb1,ref1assocs)
      refcomb2<-rbind(refcomb2,ref2assocs)
      refcomb3<-rbind(refcomb3,ref3assocs)
      refcomb4<-rbind(refcomb4,ref4assocs)
      r.expcomb1<-rbind(r.expcomb1,r.expassocs)
    }
  }
  
  #Reformatting Menarche data for processing (clumping SNPs etc)
  Clump.datA<-format_data(
    refcomb1,
    type = "exposure",
    effect_allele_col = "effect_allele.outcome",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    eaf_col = "eaf.outcome",
    phenotype_col = "outcome",
    header = T
  )
  
  A_clump<-clump_data(Clump.datA, clump_r2 = 0.01)
  
  Clump.datB<-format_data(
    refcomb2,
    type = "exposure",
    effect_allele_col = "effect_allele.outcome",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    eaf_col = "eaf.outcome",
    phenotype_col = "outcome",
    header = T
  )
  
  B_clump<-clump_data(Clump.datB, clump_r2 = 0.01)
  
  Clump.datC<-format_data(
    refcomb3,
    type = "exposure",
    effect_allele_col = "effect_allele.outcome",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    eaf_col = "eaf.outcome",
    phenotype_col = "outcome",
    header = T
  )
  
  C_clump<-clump_data(Clump.datC, clump_r2 = 0.01)
  
  Clump.datD<-format_data(
    refcomb4,
    type = "exposure",
    effect_allele_col = "effect_allele.outcome",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    eaf_col = "eaf.outcome",
    phenotype_col = "outcome",
    header = T
  )
  
  D_clump<-clump_data(Clump.datD, clump_r2 = 0.01)
  
  Clump.datE<-format_data(
    r.expcomb1,
    type = "exposure",
    effect_allele_col = "effect_allele.outcome",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    eaf_col = "eaf.outcome",
    phenotype_col = "outcome",
    header = T
  )
  
  E_clump<-clump_data(Clump.datE, clump_r2 = 0.01)
  
  if(nrow(A_clump) > nrow(B_clump)){
    
    B_clump<-Clump.datB[Clump.datB$SNP %in% A_clump$SNP,]
    
  }else{
    
    A_clump<-Clump.datA[Clump.datA$SNP %in% B_clump$SNP,]
    
  }
  
  if(nrow(C_clump) > nrow(B_clump)){
    
    A_clump<-Clump.datA[Clump.datA$SNP %in% C_clump$SNP,]
    B_clump<-Clump.datB[Clump.datB$SNP %in% C_clump$SNP,]
    
  }else{
    
    C_clump<-Clump.datC[Clump.datC$SNP %in% A_clump$SNP,]
    
  }
  
  if(nrow(D_clump) > nrow(C_clump)){
    
    A_clump<-Clump.datA[Clump.datA$SNP %in% D_clump$SNP,]
    B_clump<-Clump.datB[Clump.datB$SNP %in% D_clump$SNP,]
    C_clump<-Clump.datC[Clump.datC$SNP %in% D_clump$SNP,]
    
  }else{
    
    D_clump<-Clump.datD[Clump.datD$SNP %in% A_clump$SNP,]
    
  }
  
  if(nrow(E_clump) > nrow(D_clump)){
    
    A_clump<-Clump.datA[Clump.datA$SNP %in% E_clump$SNP,]
    B_clump<-Clump.datB[Clump.datB$SNP %in% E_clump$SNP,]
    C_clump<-Clump.datC[Clump.datC$SNP %in% E_clump$SNP,]
    D_clump<-Clump.datD[Clump.datD$SNP %in% E_clump$SNP,]
    
  }else{
    
    E_clump<-Clump.datE[Clump.datE$SNP %in% A_clump$SNP,]
    
  }
  
  # Obtain CHD outcome associations
  Outcome.data<- extract_outcome_data(
    snps = A_clump$SNP, outcomes = "ieu-a-7",proxies=F)
  
  # Harmonize exposure data
  HarmA_data <- harmonise_data(
    exposure_dat = A_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmB_data <- harmonise_data(
    exposure_dat = B_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmC_data <- harmonise_data(
    exposure_dat = C_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmD_data <- harmonise_data(
    exposure_dat = D_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmE_data <- harmonise_data(
    exposure_dat = D_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmA_data<-HarmA_data[HarmA_data$mr_keep == T,]
  
  HarmB_data<-HarmB_data[HarmB_data$mr_keep == T,]
  
  HarmC_data<-HarmC_data[HarmC_data$mr_keep == T,]
  
  HarmD_data<-HarmD_data[HarmD_data$mr_keep == T,]
  
  HarmE_data<-HarmE_data[HarmE_data$mr_keep == T,]
  
  
  if(nrow(HarmA_data) < nrow(HarmB_data)){
    
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmA_data$SNP,]
    
  }else{
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmB_data$SNP,]
    
  }
  
  if(nrow(HarmC_data) < nrow(HarmB_data)){
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmC_data$SNP,]
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmC_data$SNP,]
    
  }else{
    
    HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmA_data$SNP,]
    
  }
  
  if(nrow(HarmD_data) < nrow(HarmC_data)){
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmD_data$SNP,]
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmD_data$SNP,]
    HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmD_data$SNP,]
    
  }else{
    
    HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmA_data$SNP,]
    
  }
  
  if(nrow(HarmE_data) < nrow(HarmD_data)){
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmE_data$SNP,]
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmE_data$SNP,]
    HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmE_data$SNP,]
    HarmC_data<-HarmD_data[HarmD_data$SNP %in% HarmE_data$SNP,]
    
  }else{
    
    HarmE_data<-HarmE_data[HarmE_data$SNP %in% HarmA_data$SNP,]
    
  }
  
  HarmA_data<-HarmA_data[order(HarmA_data$SNP),]
  
  HarmB_data<-HarmB_data[order(HarmB_data$SNP),]
  
  HarmC_data<-HarmC_data[order(HarmC_data$SNP),]
  
  HarmD_data<-HarmD_data[order(HarmD_data$SNP),]
  
  HarmE_data<-HarmE_data[order(HarmE_data$SNP),]
  
  #Create combined dataframe for CHD MVMR
  tBXGs<-matrix(c(HarmA_data$beta.exposure,HarmB_data$beta.exposure,HarmC_data$beta.exposure,HarmD_data$beta.exposure,HarmE_data$beta.exposure),ncol=5)
  
  tseXGs<-matrix(c(HarmA_data$se.exposure,HarmA_data$se.exposure,HarmC_data$se.exposure,HarmD_data$se.exposure,HarmE_data$se.exposure),ncol=5)
  
  tMVMR.dat<-format_mvmr(tBXGs,HarmA_data$beta.outcome,tseXGs,HarmA_data$se.outcome,HarmA_data$SNP)
  
  #Estimate conditional F-statistics for CHD
  print("CF stats for")
  print(reference1)
  print("and")
  print(reference2)
  print("and")
  print(reference3)
  print("and")
  print(reference4)
  print("on")
  print(r.exposures[i])
  CF_stat<-strength_mvmr(tMVMR.dat)
  
  if(is.null(output4)){
    output4<-data.frame(reference1,reference2,reference3,reference4,r.exposures[i],CF_stat)
  }else{
    output4<-rbind(output4,data.frame(reference1,reference2,reference3,reference4,r.exposures[i],CF_stat))
  }
  
}

write.csv(output4,"F-out4.csv",row.names=F)





# Obtain effect estimates for each outcome

LDL<-"ieu-a-300"

TG<-"ieu-a-302"

SBP<-"ukb-a-360"

HDL<-"ieu-a-299"

A<-extract_instruments(LDL,clump = F)

B<-extract_instruments(HDL,clump = F)

C<-extract_instruments(TG,clump = F)

D<-extract_instruments(SBP,clump = F)

SNPlist<-c(A$SNP,B$SNP,C$SNP,D$SNP)
SNPlist<-unique(SNPlist)
  
splitting_function <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
  
SNP_chunks<-splitting_function(SNPlist,20)
  
for(k in 1:20){
  ref1assocs <- extract_outcome_data(
    snps = unlist(SNP_chunks[k]), outcomes = LDL,proxies=F)
    
  ref2assocs <- extract_outcome_data(
    snps = unlist(SNP_chunks[k]), outcomes = HDL,proxies=F)
    
  ref3assocs <- extract_outcome_data(
    snps = unlist(SNP_chunks[k]), outcomes = TG,proxies=F)
    
  ref4assocs <- extract_outcome_data(
    snps = unlist(SNP_chunks[k]), outcomes = SBP,proxies=F)

    if(k == 1){
      refcomb1<-ref1assocs
      refcomb2<-ref2assocs
      refcomb3<-ref3assocs
      refcomb4<-ref4assocs
    }else{
      refcomb1<-rbind(refcomb1,ref1assocs)
      refcomb2<-rbind(refcomb2,ref2assocs)
      refcomb3<-rbind(refcomb3,ref3assocs)
      refcomb4<-rbind(refcomb4,ref4assocs)
    }
  }
  
  SNPlistA<-refcomb1$SNP[refcomb1$SNP %in% refcomb2$SNP]
  SNPlistB<-refcomb2$SNP[refcomb1$SNP %in% SNPlistA]
  SNPlistC<-refcomb3$SNP[refcomb3$SNP %in% SNPlistB]
  SNPlistD<-refcomb4$SNP[refcomb4$SNP %in% SNPlistC]
  SNPlistcomb<-unique(c(SNPlistA,SNPlistB,SNPlistC,SNPlistD))
  
  SNP_chunks<-splitting_function(SNPlistcomb,20)
  
  for(k in 1:20){
    ref1assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = LDL,proxies=F)
    
    ref2assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = HDL,proxies=F)
    
    ref3assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = TG,proxies=F)
    
    ref4assocs <- extract_outcome_data(
      snps = unlist(SNP_chunks[k]), outcomes = SBP,proxies=F)
    
    if(k == 1){
      refcomb1<-ref1assocs
      refcomb2<-ref2assocs
      refcomb3<-ref3assocs
      refcomb4<-ref4assocs
    }else{
      refcomb1<-rbind(refcomb1,ref1assocs)
      refcomb2<-rbind(refcomb2,ref2assocs)
      refcomb3<-rbind(refcomb3,ref3assocs)
      refcomb4<-rbind(refcomb4,ref4assocs)
    }
  }
  
  #Reformatting Menarche data for processing (clumping SNPs etc)
  Clump.datA<-format_data(
    refcomb1,
    type = "exposure",
    effect_allele_col = "effect_allele.outcome",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    eaf_col = "eaf.outcome",
    phenotype_col = "outcome",
    header = T
  )
  
  A_clump<-clump_data(Clump.datA, clump_r2 = 0.01)
  
  Clump.datB<-format_data(
    refcomb2,
    type = "exposure",
    effect_allele_col = "effect_allele.outcome",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    eaf_col = "eaf.outcome",
    phenotype_col = "outcome",
    header = T
  )
  
  B_clump<-clump_data(Clump.datB, clump_r2 = 0.01)
  
  Clump.datC<-format_data(
    refcomb3,
    type = "exposure",
    effect_allele_col = "effect_allele.outcome",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    eaf_col = "eaf.outcome",
    phenotype_col = "outcome",
    header = T
  )
  
  C_clump<-clump_data(Clump.datC, clump_r2 = 0.01)
  
  Clump.datD<-format_data(
    refcomb4,
    type = "exposure",
    effect_allele_col = "effect_allele.outcome",
    beta_col = "beta.outcome",
    se_col = "se.outcome",
    other_allele_col = "other_allele.outcome",
    pval_col = "pval.outcome",
    eaf_col = "eaf.outcome",
    phenotype_col = "outcome",
    header = T
  )
  
  D_clump<-clump_data(Clump.datD, clump_r2 = 0.01)

  if(nrow(A_clump) > nrow(B_clump)){
    
    B_clump<-Clump.datB[Clump.datB$SNP %in% A_clump$SNP,]
    
  }else{
    
    A_clump<-Clump.datA[Clump.datA$SNP %in% B_clump$SNP,]
    
  }
  
  if(nrow(C_clump) > nrow(B_clump)){
    
    A_clump<-Clump.datA[Clump.datA$SNP %in% C_clump$SNP,]
    B_clump<-Clump.datB[Clump.datB$SNP %in% C_clump$SNP,]
    
  }else{
    
    C_clump<-Clump.datC[Clump.datC$SNP %in% A_clump$SNP,]
    
  }
  
  if(nrow(D_clump) > nrow(C_clump)){
    
    A_clump<-Clump.datA[Clump.datA$SNP %in% D_clump$SNP,]
    B_clump<-Clump.datB[Clump.datB$SNP %in% D_clump$SNP,]
    C_clump<-Clump.datC[Clump.datC$SNP %in% D_clump$SNP,]
    
  }else{
    
    D_clump<-Clump.datD[Clump.datD$SNP %in% A_clump$SNP,]
    
  }
  
  ###########
  ### CHD ###
  ###########
  
  
  # Obtain CHD outcome associations
  Outcome.data<- extract_outcome_data(
    snps = A_clump$SNP, outcomes = "ieu-a-7",proxies=F)
  
  # Harmonize exposure data
  HarmA_data <- harmonise_data(
    exposure_dat = A_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmB_data <- harmonise_data(
    exposure_dat = B_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmC_data <- harmonise_data(
    exposure_dat = C_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmD_data <- harmonise_data(
    exposure_dat = D_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmA_data<-HarmA_data[HarmA_data$mr_keep == T,]
  
  HarmB_data<-HarmB_data[HarmB_data$mr_keep == T,]
  
  HarmC_data<-HarmC_data[HarmC_data$mr_keep == T,]
  
  HarmD_data<-HarmD_data[HarmD_data$mr_keep == T,]
  
  
  if(nrow(HarmA_data) < nrow(HarmB_data)){
    
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmA_data$SNP,]
    
  }else{
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmB_data$SNP,]
    
  }
  
  if(nrow(HarmC_data) < nrow(HarmB_data)){
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmC_data$SNP,]
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmC_data$SNP,]
    
  }else{
    
    HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmA_data$SNP,]
    
  }
  
  if(nrow(HarmD_data) < nrow(HarmC_data)){
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmD_data$SNP,]
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmD_data$SNP,]
    HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmD_data$SNP,]
    
  }else{
    
    HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmC_data$SNP,]
    
  }
  
  if(nrow(HarmD_data) < nrow(HarmC_data)){
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmD_data$SNP,]
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmD_data$SNP,]
    HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmD_data$SNP,]
    
  }else{
    
    HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmC_data$SNP,]
    
  }
  
  HarmA_data<-HarmA_data[order(HarmA_data$SNP),]
  
  HarmB_data<-HarmB_data[order(HarmB_data$SNP),]
  
  HarmC_data<-HarmC_data[order(HarmC_data$SNP),]
  
  HarmD_data<-HarmD_data[order(HarmD_data$SNP),]
  
  #Create combined dataframe for CHD MVMR
  CHDBXGs<-matrix(c(HarmA_data$beta.exposure,HarmB_data$beta.exposure,HarmC_data$beta.exposure,HarmD_data$beta.exposure),ncol=4)
  
  CHDseXGs<-matrix(c(HarmA_data$se.exposure,HarmA_data$se.exposure,HarmC_data$se.exposure,HarmD_data$se.exposure),ncol=4)
  
  CHDMVMR.dat<-format_mvmr(CHDBXGs,HarmA_data$beta.outcome,CHDseXGs,HarmA_data$se.outcome,HarmA_data$SNP)
  
  ###########
  ### ILA ###
  ###########
  
  
  # Obtain CHD outcome associations
  Outcome.data<- extract_outcome_data(
    snps = A_clump$SNP, outcomes = "ebi-a-GCST006907",proxies=F)
  
  # Harmonize exposure data
  HarmA_data <- harmonise_data(
    exposure_dat = A_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmB_data <- harmonise_data(
    exposure_dat = B_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmC_data <- harmonise_data(
    exposure_dat = C_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmD_data <- harmonise_data(
    exposure_dat = D_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmA_data<-HarmA_data[HarmA_data$mr_keep == T,]
  
  HarmB_data<-HarmB_data[HarmB_data$mr_keep == T,]
  
  HarmC_data<-HarmC_data[HarmC_data$mr_keep == T,]
  
  HarmD_data<-HarmD_data[HarmD_data$mr_keep == T,]
  
  
  if(nrow(HarmA_data) < nrow(HarmB_data)){
    
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmA_data$SNP,]
    
  }else{
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmB_data$SNP,]
    
  }
  
  if(nrow(HarmC_data) < nrow(HarmB_data)){
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmC_data$SNP,]
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmC_data$SNP,]
    
  }else{
    
    HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmA_data$SNP,]
    
  }
  
  if(nrow(HarmD_data) < nrow(HarmC_data)){
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmD_data$SNP,]
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmD_data$SNP,]
    HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmD_data$SNP,]
    
  }else{
    
    HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmC_data$SNP,]
    
  }
  
  if(nrow(HarmD_data) < nrow(HarmC_data)){
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmD_data$SNP,]
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmD_data$SNP,]
    HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmD_data$SNP,]
    
  }else{
    
    HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmC_data$SNP,]
    
  }
  
  HarmA_data<-HarmA_data[order(HarmA_data$SNP),]
  
  HarmB_data<-HarmB_data[order(HarmB_data$SNP),]
  
  HarmC_data<-HarmC_data[order(HarmC_data$SNP),]
  
  HarmD_data<-HarmD_data[order(HarmD_data$SNP),]
  
  #Create combined dataframe for CHD MVMR
  ILABXGs<-matrix(c(HarmA_data$beta.exposure,HarmB_data$beta.exposure,HarmC_data$beta.exposure,HarmD_data$beta.exposure),ncol=4)
  
  ILAseXGs<-matrix(c(HarmA_data$se.exposure,HarmA_data$se.exposure,HarmC_data$se.exposure,HarmD_data$se.exposure),ncol=4)
  
  ILAMVMR.dat<-format_mvmr(ILABXGs,HarmA_data$beta.outcome,ILAseXGs,HarmA_data$se.outcome,HarmA_data$SNP)
  
  ##########
  ### IC ###
  ##########
  
  # Obtain CHD outcome associations
  Outcome.data<- extract_outcome_data(
    snps = A_clump$SNP, outcomes = "ebi-a-GCST006910",proxies=F)
  
  # Harmonize exposure data
  HarmA_data <- harmonise_data(
    exposure_dat = A_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmB_data <- harmonise_data(
    exposure_dat = B_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmC_data <- harmonise_data(
    exposure_dat = C_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmD_data <- harmonise_data(
    exposure_dat = D_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmA_data<-HarmA_data[HarmA_data$mr_keep == T,]
  
  HarmB_data<-HarmB_data[HarmB_data$mr_keep == T,]
  
  HarmC_data<-HarmC_data[HarmC_data$mr_keep == T,]
  
  HarmD_data<-HarmD_data[HarmD_data$mr_keep == T,]
  
  
  if(nrow(HarmA_data) < nrow(HarmB_data)){
    
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmA_data$SNP,]
    
  }else{
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmB_data$SNP,]
    
  }
  
  if(nrow(HarmC_data) < nrow(HarmB_data)){
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmC_data$SNP,]
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmC_data$SNP,]
    
  }else{
    
    HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmA_data$SNP,]
    
  }
  
  if(nrow(HarmD_data) < nrow(HarmC_data)){
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmD_data$SNP,]
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmD_data$SNP,]
    HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmD_data$SNP,]
    
  }else{
    
    HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmC_data$SNP,]
    
  }
  
  if(nrow(HarmD_data) < nrow(HarmC_data)){
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmD_data$SNP,]
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmD_data$SNP,]
    HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmD_data$SNP,]
    
  }else{
    
    HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmC_data$SNP,]
    
  }
  
  HarmA_data<-HarmA_data[order(HarmA_data$SNP),]
  
  HarmB_data<-HarmB_data[order(HarmB_data$SNP),]
  
  HarmC_data<-HarmC_data[order(HarmC_data$SNP),]
  
  HarmD_data<-HarmD_data[order(HarmD_data$SNP),]
  
  #Create combined dataframe for CHD MVMR
  ICBXGs<-matrix(c(HarmA_data$beta.exposure,HarmB_data$beta.exposure,HarmC_data$beta.exposure,HarmD_data$beta.exposure),ncol=4)
  
  ICseXGs<-matrix(c(HarmA_data$se.exposure,HarmA_data$se.exposure,HarmC_data$se.exposure,HarmD_data$se.exposure),ncol=4)
  
  ICMVMR.dat<-format_mvmr(ICBXGs,HarmA_data$beta.outcome,ICseXGs,HarmA_data$se.outcome,HarmA_data$SNP)
  
  







  ##########
  ### ISV ###
  ##########
  
  # Obtain ISV outcome associations
  Outcome.data<- extract_outcome_data(
    snps = A_clump$SNP, outcomes = "ebi-a-GCST006909",proxies=F)
  
  # Harmonize exposure data
  HarmA_data <- harmonise_data(
    exposure_dat = A_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmB_data <- harmonise_data(
    exposure_dat = B_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmC_data <- harmonise_data(
    exposure_dat = C_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmD_data <- harmonise_data(
    exposure_dat = D_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmA_data<-HarmA_data[HarmA_data$mr_keep == T,]
  
  HarmB_data<-HarmB_data[HarmB_data$mr_keep == T,]
  
  HarmC_data<-HarmC_data[HarmC_data$mr_keep == T,]
  
  HarmD_data<-HarmD_data[HarmD_data$mr_keep == T,]
  
  
  if(nrow(HarmA_data) < nrow(HarmB_data)){
    
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmA_data$SNP,]
    
  }else{
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmB_data$SNP,]
    
  }
  
  if(nrow(HarmC_data) < nrow(HarmB_data)){
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmC_data$SNP,]
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmC_data$SNP,]
    
  }else{
    
    HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmA_data$SNP,]
    
  }
  
  if(nrow(HarmD_data) < nrow(HarmC_data)){
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmD_data$SNP,]
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmD_data$SNP,]
    HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmD_data$SNP,]
    
  }else{
    
    HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmC_data$SNP,]
    
  }
  
  if(nrow(HarmD_data) < nrow(HarmC_data)){
    
    HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmD_data$SNP,]
    HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmD_data$SNP,]
    HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmD_data$SNP,]
    
  }else{
    
    HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmC_data$SNP,]
    
  }
  
  HarmA_data<-HarmA_data[order(HarmA_data$SNP),]
  
  HarmB_data<-HarmB_data[order(HarmB_data$SNP),]
  
  HarmC_data<-HarmC_data[order(HarmC_data$SNP),]
  
  HarmD_data<-HarmD_data[order(HarmD_data$SNP),]
  
  #Create combined dataframe for CHD MVMR
  ISVBXGs<-matrix(c(HarmA_data$beta.exposure,HarmB_data$beta.exposure,HarmC_data$beta.exposure,HarmD_data$beta.exposure),ncol=4)
  
  ISVseXGs<-matrix(c(HarmA_data$se.exposure,HarmA_data$se.exposure,HarmC_data$se.exposure,HarmD_data$se.exposure),ncol=4)
  
  ISVMVMR.dat<-format_mvmr(ISVBXGs,HarmA_data$beta.outcome,ISVseXGs,HarmA_data$se.outcome,HarmA_data$SNP)

  mvmrresdatCHD<-ivw_mvmr(CHDMVMR.dat)
  mvmrresdatILA<-ivw_mvmr(ILAMVMR.dat)
  mvmrresdatIC<-ivw_mvmr(ICMVMR.dat)
  mvmrresdatISV<-ivw_mvmr(ISVMVMR.dat)
  
  #Plotting multivariable estimates
  
  or<-exp(c(mvmrresdatCHD[1,1],mvmrresdatCHD[2,1],mvmrresdatCHD[3,1],mvmrresdatCHD[4,1],
            mvmrresdatILA[1,1],mvmrresdatILA[2,1],mvmrresdatILA[3,1],mvmrresdatILA[4,1],
            mvmrresdatIC[1,1],mvmrresdatIC[2,1],mvmrresdatIC[3,1],mvmrresdatIC[4,1],
            mvmrresdatISV[1,1],mvmrresdatISV[2,1],mvmrresdatISV[3,1],mvmrresdatISV[4,1]))

  ses<-c(mvmrresdatCHD[1,2],mvmrresdatCHD[2,2],mvmrresdatCHD[3,2],mvmrresdatCHD[4,2],
         mvmrresdatILA[1,2],mvmrresdatILA[2,2],mvmrresdatILA[3,2],mvmrresdatILA[4,2],
         mvmrresdatIC[1,2],mvmrresdatIC[2,2],mvmrresdatIC[3,2],mvmrresdatIC[4,2],
         mvmrresdatISV[1,2],mvmrresdatISV[2,2],mvmrresdatISV[3,2],mvmrresdatISV[4,2])

  lci<-or-1.96*ses
  
  uci<-or+1.96*ses

  exposure<-rep(c("LDL","HDL","TG","SBP"),4)

  outcome<-c(rep("CHD",4),rep("I.Stroke (large-artery)",4),rep("I.Stroke (cardioembolic)",4),rep("I.Stroke (small vessel)",4))
  
  mvmrres_data<-data.frame(or,lci,uci,outcome,exposure)

  names(mvmrres_data)<-c("or","orlci_95","oruci_95","outcome","exposure")
  
  mvmrres_data<-mvmrres_data[c(rev(1:4),rev(5:8),rev(9:12),rev(13:16)),]
  
  mvmrres_data<-mvmrres_data[c(13:16,9:12,5:8,1:4),]
  
  mvmrres_data$index<-c(seq(1,1.75,0.25),
                        seq(2.5,3.25,0.25),
                        seq(4,4.75,0.25),
                        seq(5.5,6.25,0.25))
  
  LDL_CF<-mean(c(strength_mvmr(CHDMVMR.dat)[1,1],strength_mvmr(ILAMVMR.dat)[1,1],
         strength_mvmr(ICMVMR.dat)[1,1],strength_mvmr(ISVMVMR.dat)[1,1]))
  
  HDL_CF<-mean(c(strength_mvmr(CHDMVMR.dat)[1,2],strength_mvmr(ILAMVMR.dat)[1,2],
                 strength_mvmr(ICMVMR.dat)[1,2],strength_mvmr(ISVMVMR.dat)[1,2]))
  
  TG_CF<-mean(c(strength_mvmr(CHDMVMR.dat)[1,3],strength_mvmr(ILAMVMR.dat)[1,3],
                 strength_mvmr(ICMVMR.dat)[1,3],strength_mvmr(ISVMVMR.dat)[1,3]))
  
 SBP_CF<-mean(c(strength_mvmr(CHDMVMR.dat)[1,4],strength_mvmr(ILAMVMR.dat)[1,4],
                 strength_mvmr(ICMVMR.dat)[1,4],strength_mvmr(ISVMVMR.dat)[1,4]))
  
  ######
  
  mvmrplot1<-ggplot(mvmrres_data)+geom_point(aes(x=index, y=or))+
    coord_flip()+geom_errorbar(aes(ymin=orlci_95,ymax=oruci_95,x=index),width=0,orientation = "x")+
    scale_y_continuous(limits = c(-2,4),name = "Effect (odds ratio)")+geom_hline(yintercept=1,linetype = "dashed")+theme_bw()+
    scale_x_continuous(limits = c(1,7.25),breaks=rev(c(1.75,3.25,4.75,6.25)+0.375),labels= mvmrres_data$outcome[c(14,9,5,1)],name = "Outcome")+
    scale_color_discrete(name = "Exposure")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  png("Figure5_1.png",width=6,height=6,units="in",res=200)
  mvmrplot1
  dev.off()

































