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
exposure_data_BMI<-extract_instruments("ukb-b-19953",clump = T)
exposure_data_SBP<-extract_instruments("ukb-a-360",clump = T)
exposure_data_DBP<-extract_instruments("ukb-a-359",clump = T)
exposure_data_LDL<-extract_instruments("ieu-b-110",clump = T)
exposure_data_HDL<-extract_instruments("ieu-b-109",clump = T)
exposure_data_TG<-extract_instruments("ieu-b-111",clump = T)
exposure_data_PA<-extract_instruments("ukb-a-511",clump = T)
exposure_data_DIA<-extract_instruments("ieu-a-26",clump = T)
exposure_data_SMK<-extract_instruments("ukb-b-2134",clump = T)
#Flip smoking to reflect increase in smoking amount
exposure_data_SMK$beta.exposure<-exposure_data_SMK$beta.exposure*-1
exposure_data_SLP<-extract_instruments("ukb-b-4424",clump = T)
exposure_data_ALC<-extract_instruments("ieu-a-1283",clump = T)
exposure_data_KF<-extract_instruments("ukb-a-333",clump = T)

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

outcome_data_DBPCHD<- extract_outcome_data(
  snps = exposure_data_DBP$SNP, outcomes = "ieu-a-7",proxies=F)

outcome_data_DBPMI<- extract_outcome_data(
  snps = exposure_data_DBP$SNP, outcomes = "ieu-a-798",proxies=F)

outcome_data_DBPILA<- extract_outcome_data(
  snps = exposure_data_DBP$SNP, outcomes = "ebi-a-GCST006907",proxies=F)

outcome_data_DBPIC<- extract_outcome_data(
  snps = exposure_data_DBP$SNP, outcomes = "ebi-a-GCST006910",proxies=F)

outcome_data_DBPISV<- extract_outcome_data(
  snps = exposure_data_DBP$SNP, outcomes = "ebi-a-GCST006909",proxies=F)

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

outcome_data_PACHD<- extract_outcome_data(
  snps = exposure_data_PA$SNP, outcomes = "ieu-a-7",proxies=F)

outcome_data_PAMI<- extract_outcome_data(
  snps = exposure_data_PA$SNP, outcomes = "ieu-a-798",proxies=F)

outcome_data_PAILA<- extract_outcome_data(
  snps = exposure_data_PA$SNP, outcomes = "ebi-a-GCST006907",proxies=F)

outcome_data_PAIC<- extract_outcome_data(
  snps = exposure_data_PA$SNP, outcomes = "ebi-a-GCST006910",proxies=F)

outcome_data_PAISV<- extract_outcome_data(
  snps = exposure_data_PA$SNP, outcomes = "ebi-a-GCST006909",proxies=F)

outcome_data_DIACHD<- extract_outcome_data(
  snps = exposure_data_DIA$SNP, outcomes = "ieu-a-7",proxies=F)

outcome_data_DIAMI<- extract_outcome_data(
  snps = exposure_data_DIA$SNP, outcomes = "ieu-a-798",proxies=F)

outcome_data_DIAILA<- extract_outcome_data(
  snps = exposure_data_DIA$SNP, outcomes = "ebi-a-GCST006907",proxies=F)

outcome_data_DIAIC<- extract_outcome_data(
  snps = exposure_data_DIA$SNP, outcomes = "ebi-a-GCST006910",proxies=F)

outcome_data_DIAISV<- extract_outcome_data(
  snps = exposure_data_DIA$SNP, outcomes = "ebi-a-GCST006909",proxies=F)

outcome_data_SMKCHD<- extract_outcome_data(
  snps = exposure_data_SMK$SNP, outcomes = "ieu-a-7",proxies=F)

outcome_data_SMKMI<- extract_outcome_data(
  snps = exposure_data_SMK$SNP, outcomes = "ieu-a-798",proxies=F)

outcome_data_SMKILA<- extract_outcome_data(
  snps = exposure_data_SMK$SNP, outcomes = "ebi-a-GCST006907",proxies=F)

outcome_data_SMKIC<- extract_outcome_data(
  snps = exposure_data_SMK$SNP, outcomes = "ebi-a-GCST006910",proxies=F)

outcome_data_SMKISV<- extract_outcome_data(
  snps = exposure_data_SMK$SNP, outcomes = "ebi-a-GCST006909",proxies=F)

outcome_data_SLPCHD<- extract_outcome_data(
  snps = exposure_data_SLP$SNP, outcomes = "ieu-a-7",proxies=F)

outcome_data_SLPMI<- extract_outcome_data(
  snps = exposure_data_SLP$SNP, outcomes = "ieu-a-798",proxies=F)

outcome_data_SLPILA<- extract_outcome_data(
  snps = exposure_data_SLP$SNP, outcomes = "ebi-a-GCST006907",proxies=F)

outcome_data_SLPIC<- extract_outcome_data(
  snps = exposure_data_SLP$SNP, outcomes = "ebi-a-GCST006910",proxies=F)

outcome_data_SLPISV<- extract_outcome_data(
  snps = exposure_data_SLP$SNP, outcomes = "ebi-a-GCST006909",proxies=F)

outcome_data_ALCCHD<- extract_outcome_data(
  snps = exposure_data_ALC$SNP, outcomes = "ieu-a-7",proxies=F)

outcome_data_ALCMI<- extract_outcome_data(
  snps = exposure_data_ALC$SNP, outcomes = "ieu-a-798",proxies=F)

outcome_data_ALCILA<- extract_outcome_data(
  snps = exposure_data_ALC$SNP, outcomes = "ebi-a-GCST006907",proxies=F)

outcome_data_ALCIC<- extract_outcome_data(
  snps = exposure_data_ALC$SNP, outcomes = "ebi-a-GCST006910",proxies=F)

outcome_data_ALCISV<- extract_outcome_data(
  snps = exposure_data_ALC$SNP, outcomes = "ebi-a-GCST006909",proxies=F)

outcome_data_KFCHD<- extract_outcome_data(
  snps = exposure_data_KF$SNP, outcomes = "ieu-a-7",proxies=F)

outcome_data_KFMI<- extract_outcome_data(
  snps = exposure_data_KF$SNP, outcomes = "ieu-a-798",proxies=F)

outcome_data_KFILA<- extract_outcome_data(
  snps = exposure_data_KF$SNP, outcomes = "ebi-a-GCST006907",proxies=F)

outcome_data_KFIC<- extract_outcome_data(
  snps = exposure_data_KF$SNP, outcomes = "ebi-a-GCST006910",proxies=F)

outcome_data_KFISV<- extract_outcome_data(
  snps = exposure_data_KF$SNP, outcomes = "ebi-a-GCST006909",proxies=F)


##################################################################################


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

DBPCHD.dat <- harmonise_data(
  exposure_dat = exposure_data_DBP, 
  outcome_dat = outcome_data_DBPCHD
)

DBPMI.dat <- harmonise_data(
  exposure_dat = exposure_data_DBP, 
  outcome_dat = outcome_data_DBPMI
)

DBPIStkLA.dat <- harmonise_data(
  exposure_dat = exposure_data_DBP, 
  outcome_dat = outcome_data_DBPILA
)

DBPIStkC.dat <- harmonise_data(
  exposure_dat = exposure_data_DBP, 
  outcome_dat = outcome_data_DBPIC
)

DBPIStkSV.dat <- harmonise_data(
  exposure_dat = exposure_data_DBP, 
  outcome_dat = outcome_data_DBPISV
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

PACHD.dat <- harmonise_data(
  exposure_dat = exposure_data_PA, 
  outcome_dat = outcome_data_PACHD
)

PAMI.dat <- harmonise_data(
  exposure_dat = exposure_data_PA, 
  outcome_dat = outcome_data_PAMI
)

PAIStkLA.dat <- harmonise_data(
  exposure_dat = exposure_data_PA, 
  outcome_dat = outcome_data_PAILA
)

PAIStkC.dat <- harmonise_data(
  exposure_dat = exposure_data_PA, 
  outcome_dat = outcome_data_PAIC
)

PAIStkSV.dat <- harmonise_data(
  exposure_dat = exposure_data_PA, 
  outcome_dat = outcome_data_PAISV
)

DIACHD.dat <- harmonise_data(
  exposure_dat = exposure_data_DIA, 
  outcome_dat = outcome_data_DIACHD
)

DIAMI.dat <- harmonise_data(
  exposure_dat = exposure_data_DIA, 
  outcome_dat = outcome_data_DIAMI
)

DIAIStkLA.dat <- harmonise_data(
  exposure_dat = exposure_data_DIA, 
  outcome_dat = outcome_data_DIAILA
)

DIAIStkC.dat <- harmonise_data(
  exposure_dat = exposure_data_DIA, 
  outcome_dat = outcome_data_DIAIC
)

DIAIStkSV.dat <- harmonise_data(
  exposure_dat = exposure_data_DIA, 
  outcome_dat = outcome_data_DIAISV
)

SMKCHD.dat <- harmonise_data(
  exposure_dat = exposure_data_SMK, 
  outcome_dat = outcome_data_SMKCHD
)

SMKMI.dat <- harmonise_data(
  exposure_dat = exposure_data_SMK, 
  outcome_dat = outcome_data_SMKMI
)

SMKIStkLA.dat <- harmonise_data(
  exposure_dat = exposure_data_SMK, 
  outcome_dat = outcome_data_SMKILA
)

SMKIStkC.dat <- harmonise_data(
  exposure_dat = exposure_data_SMK, 
  outcome_dat = outcome_data_SMKIC
)

SMKIStkSV.dat <- harmonise_data(
  exposure_dat = exposure_data_SMK, 
  outcome_dat = outcome_data_SMKISV
)

SLPCHD.dat <- harmonise_data(
  exposure_dat = exposure_data_SLP, 
  outcome_dat = outcome_data_SLPCHD
)

SLPMI.dat <- harmonise_data(
  exposure_dat = exposure_data_SLP, 
  outcome_dat = outcome_data_SLPMI
)

SLPIStkLA.dat <- harmonise_data(
  exposure_dat = exposure_data_SLP, 
  outcome_dat = outcome_data_SLPILA
)

SLPIStkC.dat <- harmonise_data(
  exposure_dat = exposure_data_SLP, 
  outcome_dat = outcome_data_SLPIC
)

SLPIStkSV.dat <- harmonise_data(
  exposure_dat = exposure_data_SLP, 
  outcome_dat = outcome_data_SLPISV
)

ALCCHD.dat <- harmonise_data(
  exposure_dat = exposure_data_ALC, 
  outcome_dat = outcome_data_ALCCHD
)

ALCMI.dat <- harmonise_data(
  exposure_dat = exposure_data_ALC, 
  outcome_dat = outcome_data_ALCMI
)

ALCIStkLA.dat <- harmonise_data(
  exposure_dat = exposure_data_ALC, 
  outcome_dat = outcome_data_ALCILA
)

ALCIStkC.dat <- harmonise_data(
  exposure_dat = exposure_data_ALC, 
  outcome_dat = outcome_data_ALCIC
)

ALCIStkSV.dat <- harmonise_data(
  exposure_dat = exposure_data_ALC, 
  outcome_dat = outcome_data_ALCISV
)

KFCHD.dat <- harmonise_data(
  exposure_dat = exposure_data_KF, 
  outcome_dat = outcome_data_KFCHD
)

KFMI.dat <- harmonise_data(
  exposure_dat = exposure_data_KF, 
  outcome_dat = outcome_data_KFMI
)

KFIStkLA.dat <- harmonise_data(
  exposure_dat = exposure_data_KF, 
  outcome_dat = outcome_data_KFILA
)

KFIStkC.dat <- harmonise_data(
  exposure_dat = exposure_data_KF, 
  outcome_dat = outcome_data_KFIC
)

KFIStkSV.dat <- harmonise_data(
  exposure_dat = exposure_data_KF, 
  outcome_dat = outcome_data_KFISV
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

mr_resultsDBPCHD <- generate_odds_ratios(mr(DBPCHD.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsDBPMI <- generate_odds_ratios(mr(DBPMI.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsDBPILA <- generate_odds_ratios(mr(DBPIStkLA.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsDBPIC <- generate_odds_ratios(mr(DBPIStkC.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsDBPISV <- generate_odds_ratios(mr(DBPIStkSV.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

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

mr_resultsPACHD <- generate_odds_ratios(mr(PACHD.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsPAMI <- generate_odds_ratios(mr(PAMI.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsPAILA <- generate_odds_ratios(mr(PAIStkLA.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsPAIC <- generate_odds_ratios(mr(PAIStkC.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsPAISV <- generate_odds_ratios(mr(PAIStkSV.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

mr_resultsDIACHD <- generate_odds_ratios(mr(DIACHD.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsDIAMI <- generate_odds_ratios(mr(DIAMI.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsDIAILA <- generate_odds_ratios(mr(DIAIStkLA.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsDIAIC <- generate_odds_ratios(mr(DIAIStkC.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsDIAISV <- generate_odds_ratios(mr(DIAIStkSV.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

mr_resultsSMKCHD <- generate_odds_ratios(mr(SMKCHD.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsSMKMI <- generate_odds_ratios(mr(SMKMI.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsSMKILA <- generate_odds_ratios(mr(SMKIStkLA.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsSMKIC <- generate_odds_ratios(mr(SMKIStkC.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsSMKISV <- generate_odds_ratios(mr(SMKIStkSV.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

mr_resultsSLPCHD <- generate_odds_ratios(mr(SLPCHD.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsSLPMI <- generate_odds_ratios(mr(SLPMI.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsSLPILA <- generate_odds_ratios(mr(SLPIStkLA.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsSLPIC <- generate_odds_ratios(mr(SLPIStkC.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsSLPISV <- generate_odds_ratios(mr(SLPIStkSV.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

mr_resultsALCCHD <- generate_odds_ratios(mr(ALCCHD.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsALCMI <- generate_odds_ratios(mr(ALCMI.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsALCILA <- generate_odds_ratios(mr(ALCIStkLA.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsALCIC <- generate_odds_ratios(mr(ALCIStkC.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsALCISV <- generate_odds_ratios(mr(ALCIStkSV.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

mr_resultsKFCHD <- generate_odds_ratios(mr(KFCHD.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsKFMI <- generate_odds_ratios(mr(KFMI.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsKFILA <- generate_odds_ratios(mr(KFIStkLA.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsKFIC <- generate_odds_ratios(mr(KFIStkC.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))
mr_resultsKFISV <- generate_odds_ratios(mr(KFIStkSV.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

#Mean F-statistics and mean SNP numbers

F1<- MenCHD.dat$beta.exposure^2/MenCHD.dat$se.exposure^2
F2<- MenIStkLA.dat$beta.exposure^2/MenIStkLA.dat$se.exposure^2
F3<- MenIStkC.dat$beta.exposure^2/MenIStkC.dat$se.exposure^2
F4<- MenIStkSV.dat$beta.exposure^2/MenIStkSV.dat$se.exposure^2
F5<- BMICHD.dat$beta.exposure^2/BMICHD.dat$se.exposure^2
F6<- BMIIStkLA.dat$beta.exposure^2/BMIIStkLA.dat$se.exposure^2
F7<- BMIIStkC.dat$beta.exposure^2/BMIIStkC.dat$se.exposure^2
F8<- BMIIStkSV.dat$beta.exposure^2/BMIIStkSV.dat$se.exposure^2
F9<- SBPCHD.dat$beta.exposure^2/SBPCHD.dat$se.exposure^2
F10<- SBPIStkLA.dat$beta.exposure^2/SBPIStkLA.dat$se.exposure^2
F11<- SBPIStkC.dat$beta.exposure^2/SBPIStkC.dat$se.exposure^2
F12<- SBPIStkSV.dat$beta.exposure^2/SBPIStkSV.dat$se.exposure^2
F13<- DBPCHD.dat$beta.exposure^2/SBPCHD.dat$se.exposure^2
F14<- DBPIStkLA.dat$beta.exposure^2/DBPIStkLA.dat$se.exposure^2
F15<- DBPIStkC.dat$beta.exposure^2/DBPIStkC.dat$se.exposure^2
F16<- DBPIStkSV.dat$beta.exposure^2/DBPIStkSV.dat$se.exposure^2
F17<- LDLCHD.dat$beta.exposure^2/LDLCHD.dat$se.exposure^2
F18<- LDLIStkLA.dat$beta.exposure^2/LDLIStkLA.dat$se.exposure^2
F19<- LDLIStkC.dat$beta.exposure^2/LDLIStkC.dat$se.exposure^2
F20<- LDLIStkSV.dat$beta.exposure^2/LDLIStkSV.dat$se.exposure^2
F21<- HDLCHD.dat$beta.exposure^2/HDLCHD.dat$se.exposure^2
F22<- HDLIStkLA.dat$beta.exposure^2/HDLIStkLA.dat$se.exposure^2
F23<- HDLIStkC.dat$beta.exposure^2/HDLIStkC.dat$se.exposure^2
F24<- HDLIStkSV.dat$beta.exposure^2/HDLIStkSV.dat$se.exposure^2
F25<- TGCHD.dat$beta.exposure^2/TGCHD.dat$se.exposure^2
F26<- TGIStkLA.dat$beta.exposure^2/TGIStkLA.dat$se.exposure^2
F27<- TGIStkC.dat$beta.exposure^2/TGIStkC.dat$se.exposure^2
F28<- TGIStkSV.dat$beta.exposure^2/TGIStkSV.dat$se.exposure^2
F29<- PACHD.dat$beta.exposure^2/PACHD.dat$se.exposure^2
F30<- PAIStkLA.dat$beta.exposure^2/PAIStkLA.dat$se.exposure^2
F31<- PAIStkC.dat$beta.exposure^2/PAIStkC.dat$se.exposure^2
F32<- PAIStkSV.dat$beta.exposure^2/PAIStkSV.dat$se.exposure^2
F33<- DIACHD.dat$beta.exposure^2/DIACHD.dat$se.exposure^2
F34<- DIAIStkLA.dat$beta.exposure^2/DIAIStkLA.dat$se.exposure^2
F35<- DIAIStkC.dat$beta.exposure^2/DIAIStkC.dat$se.exposure^2
F36<- DIAIStkSV.dat$beta.exposure^2/DIAIStkSV.dat$se.exposure^2
F37<- SMKCHD.dat$beta.exposure^2/SMKCHD.dat$se.exposure^2
F38<- SMKIStkLA.dat$beta.exposure^2/SMKIStkLA.dat$se.exposure^2
F39<- SMKIStkC.dat$beta.exposure^2/SMKIStkC.dat$se.exposure^2
F40<- SMKIStkSV.dat$beta.exposure^2/SMKIStkSV.dat$se.exposure^2
F41<- SLPCHD.dat$beta.exposure^2/SLPCHD.dat$se.exposure^2
F42<- SLPIStkLA.dat$beta.exposure^2/SLPIStkLA.dat$se.exposure^2
F43<- SLPIStkC.dat$beta.exposure^2/SLPIStkC.dat$se.exposure^2
F44<- SLPIStkSV.dat$beta.exposure^2/SLPIStkSV.dat$se.exposure^2
F45<- ALCCHD.dat$beta.exposure^2/ALCCHD.dat$se.exposure^2
F46<- ALCIStkLA.dat$beta.exposure^2/ALCIStkLA.dat$se.exposure^2
F47<- ALCIStkC.dat$beta.exposure^2/ALCIStkC.dat$se.exposure^2
F48<- ALCIStkSV.dat$beta.exposure^2/ALCIStkSV.dat$se.exposure^2
F49<- KFCHD.dat$beta.exposure^2/KFCHD.dat$se.exposure^2
F50<- KFIStkLA.dat$beta.exposure^2/KFIStkLA.dat$se.exposure^2
F51<- KFIStkC.dat$beta.exposure^2/KFIStkC.dat$se.exposure^2
F52<- KFIStkSV.dat$beta.exposure^2/KFIStkSV.dat$se.exposure^2

F_stats1<-c(mean(F1),mean(F2),mean(F3),mean(F4),mean(F5),
            mean(F6),mean(F7),mean(F8),mean(F9),mean(F10),
            mean(F11),mean(F12),mean(F13),mean(F14),mean(F15),
            mean(F16),mean(F17),mean(F18),mean(F19),mean(F20),
            mean(F21),mean(F22),mean(F23),mean(F24),mean(F25),
            mean(F26),mean(F27),mean(F28),mean(F29),mean(F30),
            mean(F31),mean(F32),mean(F33),mean(F34),mean(F35),
            mean(F36),mean(F37),mean(F38),mean(F39),mean(F40),
            mean(F41),mean(F42),mean(F43),mean(F44),mean(F45),
            mean(F46),mean(F47),mean(F48),mean(F49),mean(F50),
            mean(F51),mean(F52))
min(F_stats1)
max(F_stats1)

SNPnumMen<-mean(mr_resultsMenCHD$nsnp[1],
                mr_resultsMenILA$nsnp[1],
                mr_resultsMenIC$nsnp[1],
                mr_resultsMenISV$nsnp[1])
SNPnumBMI<-mean(mr_resultsBMICHD$nsnp[1],
                mr_resultsBMIILA$nsnp[1],
                mr_resultsBMIIC$nsnp[1],
                mr_resultsBMIISV$nsnp[1])
SNPnumSBP<-mean(mr_resultsSBPCHD$nsnp[1],
                mr_resultsSBPILA$nsnp[1],
                mr_resultsSBPIC$nsnp[1],
                mr_resultsSBPISV$nsnp[1])
SNPnumDBP<-mean(mr_resultsDBPCHD$nsnp[1],
                mr_resultsDBPILA$nsnp[1],
                mr_resultsDBPIC$nsnp[1],
                mr_resultsDBPISV$nsnp[1])
SNPnumLDL<-mean(mr_resultsLDLCHD$nsnp[1],
                mr_resultsLDLILA$nsnp[1],
                mr_resultsLDLIC$nsnp[1],
                mr_resultsLDLISV$nsnp[1])
SNPnumHDL<-mean(mr_resultsHDLCHD$nsnp[1],
                mr_resultsHDLILA$nsnp[1],
                mr_resultsHDLIC$nsnp[1],
                mr_resultsHDLISV$nsnp[1])
SNPnumTG<-mean(mr_resultsTGCHD$nsnp[1],
               mr_resultsTGILA$nsnp[1],
               mr_resultsTGIC$nsnp[1],
               mr_resultsTGISV$nsnp[1])
SNPnumPA<-mean(mr_resultsPACHD$nsnp[1],
               mr_resultsPAILA$nsnp[1],
               mr_resultsPAIC$nsnp[1],
               mr_resultsPAISV$nsnp[1])
SNPnumDIA<-mean(mr_resultsDIACHD$nsnp[1],
                mr_resultsDIAILA$nsnp[1],
                mr_resultsDIAIC$nsnp[1],
                mr_resultsDIAISV$nsnp[1])
SNPnumSMK<-mean(mr_resultsSMKCHD$nsnp[1],
                mr_resultsSMKILA$nsnp[1],
                mr_resultsSMKIC$nsnp[1],
                mr_resultsSMKISV$nsnp[1])
SNPnumSLP<-mean(mr_resultsSLPCHD$nsnp[1],
                mr_resultsSLPILA$nsnp[1],
                mr_resultsSLPIC$nsnp[1],
                mr_resultsSLPISV$nsnp[1])
SNPnumALC<-mean(mr_resultsALCCHD$nsnp[1],
                mr_resultsALCILA$nsnp[1],
                mr_resultsALCIC$nsnp[1],
                mr_resultsALCISV$nsnp[1])
SNPnumKF<-mean(mr_resultsKFCHD$nsnp[1],
               mr_resultsKFILA$nsnp[1],
               mr_resultsKFIC$nsnp[1],
               mr_resultsKFISV$nsnp[1])

sum(SNPnumMen,SNPnumBMI,SNPnumDBP,SNPnumSBP,SNPnumLDL,SNPnumHDL,
    SNPnumTG,SNPnumPA,SNPnumDIA,SNPnumSMK,SNPnumSLP,SNPnumALC,
    SNPnumKF)

# Plotting univariable estimates

res_data<-data.frame(c(mr_resultsMenCHD$or,mr_resultsMenISV$or,mr_resultsMenILA$or,mr_resultsMenIC$or,
                       mr_resultsBMICHD$or,mr_resultsBMIISV$or,mr_resultsBMIILA$or,mr_resultsBMIIC$or,
                       mr_resultsDBPCHD$or,mr_resultsDBPISV$or,mr_resultsDBPILA$or,mr_resultsDBPIC$or,
                       mr_resultsSBPCHD$or,mr_resultsSBPISV$or,mr_resultsSBPILA$or,mr_resultsSBPIC$or,
                       mr_resultsLDLCHD$or,mr_resultsLDLISV$or,mr_resultsLDLILA$or,mr_resultsLDLIC$or,
                       mr_resultsHDLCHD$or,mr_resultsHDLISV$or,mr_resultsHDLILA$or,mr_resultsHDLIC$or,
                       mr_resultsTGCHD$or,mr_resultsTGISV$or,mr_resultsTGILA$or,mr_resultsTGIC$or,
                       mr_resultsPACHD$or,mr_resultsPAISV$or,mr_resultsPAILA$or,mr_resultsPAIC$or,
                       mr_resultsDIACHD$or,mr_resultsDIAISV$or,mr_resultsDIAILA$or,mr_resultsDIAIC$or,
                       mr_resultsSMKCHD$or,mr_resultsSMKISV$or,mr_resultsSMKILA$or,mr_resultsSMKIC$or,
                       mr_resultsSLPCHD$or,mr_resultsSLPISV$or,mr_resultsSLPILA$or,mr_resultsSLPIC$or,
                       mr_resultsALCCHD$or,mr_resultsALCISV$or,mr_resultsALCILA$or,mr_resultsALCIC$or,
                       mr_resultsKFCHD$or,mr_resultsKFISV$or,mr_resultsKFILA$or,mr_resultsKFIC$or),
                     c(mr_resultsMenCHD$or_lci95,mr_resultsMenISV$or_lci95,mr_resultsMenILA$or_lci95,mr_resultsMenIC$or_lci95,
                       mr_resultsBMICHD$or_lci95,mr_resultsBMIISV$or_lci95,mr_resultsBMIILA$or_lci95,mr_resultsBMIIC$or_lci95,
                       mr_resultsDBPCHD$or_lci95,mr_resultsDBPISV$or_lci95,mr_resultsDBPILA$or_lci95,mr_resultsDBPIC$or_lci95,
                       mr_resultsSBPCHD$or_lci95,mr_resultsSBPISV$or_lci95,mr_resultsSBPILA$or_lci95,mr_resultsSBPIC$or_lci95,
                       mr_resultsLDLCHD$or_lci95,mr_resultsLDLISV$or_lci95,mr_resultsLDLILA$or_lci95,mr_resultsLDLIC$or_lci95,
                       mr_resultsHDLCHD$or_lci95,mr_resultsHDLISV$or_lci95,mr_resultsHDLILA$or_lci95,mr_resultsHDLIC$or_lci95,
                       mr_resultsTGCHD$or_lci95,mr_resultsTGISV$or_lci95,mr_resultsTGILA$or_lci95,mr_resultsTGIC$or_lci95,
                       mr_resultsPACHD$or_lci95,mr_resultsPAISV$or_lci95,mr_resultsPAILA$or_lci95,mr_resultsPAIC$or_lci95,
                       mr_resultsDIACHD$or_lci95,mr_resultsDIAISV$or_lci95,mr_resultsDIAILA$or_lci95,mr_resultsDIAIC$or_lci95,
                       mr_resultsSMKCHD$or_lci95,mr_resultsSMKISV$or_lci95,mr_resultsSMKILA$or_lci95,mr_resultsSMKIC$or_lci95,
                       mr_resultsSLPCHD$or_lci95,mr_resultsSLPISV$or_lci95,mr_resultsSLPILA$or_lci95,mr_resultsSLPIC$or_lci95,
                       mr_resultsALCCHD$or_lci95,mr_resultsALCISV$or_lci95,mr_resultsALCILA$or_lci95,mr_resultsALCIC$or_lci95,
                       mr_resultsKFCHD$or_lci95,mr_resultsKFISV$or_lci95,mr_resultsKFILA$or_lci95,mr_resultsKFIC$or_lci95),
                     c(mr_resultsMenCHD$or_uci95,mr_resultsMenISV$or_uci95,mr_resultsMenILA$or_uci95,mr_resultsMenIC$or_uci95,
                       mr_resultsBMICHD$or_uci95,mr_resultsBMIISV$or_uci95,mr_resultsBMIILA$or_uci95,mr_resultsBMIIC$or_uci95,
                       mr_resultsDBPCHD$or_uci95,mr_resultsDBPISV$or_uci95,mr_resultsDBPILA$or_uci95,mr_resultsDBPIC$or_uci95,
                       mr_resultsSBPCHD$or_uci95,mr_resultsSBPISV$or_uci95,mr_resultsSBPILA$or_uci95,mr_resultsSBPIC$or_uci95,
                       mr_resultsLDLCHD$or_uci95,mr_resultsLDLISV$or_uci95,mr_resultsLDLILA$or_uci95,mr_resultsLDLIC$or_uci95,
                       mr_resultsHDLCHD$or_uci95,mr_resultsHDLISV$or_uci95,mr_resultsHDLILA$or_uci95,mr_resultsHDLIC$or_uci95,
                       mr_resultsTGCHD$or_uci95,mr_resultsTGISV$or_uci95,mr_resultsTGILA$or_uci95,mr_resultsTGIC$or_uci95,
                       mr_resultsPACHD$or_uci95,mr_resultsPAISV$or_uci95,mr_resultsPAILA$or_uci95,mr_resultsPAIC$or_uci95,
                       mr_resultsDIACHD$or_uci95,mr_resultsDIAISV$or_uci95,mr_resultsDIAILA$or_uci95,mr_resultsDIAIC$or_uci95,
                       mr_resultsSMKCHD$or_uci95,mr_resultsSMKISV$or_uci95,mr_resultsSMKILA$or_uci95,mr_resultsSMKIC$or_uci95,
                       mr_resultsSLPCHD$or_uci95,mr_resultsSLPISV$or_uci95,mr_resultsSLPILA$or_uci95,mr_resultsSLPIC$or_uci95,
                       mr_resultsALCCHD$or_uci95,mr_resultsALCISV$or_uci95,mr_resultsALCILA$or_uci95,mr_resultsALCIC$or_uci95,
                       mr_resultsKFCHD$or_uci95,mr_resultsKFISV$or_uci95,mr_resultsKFILA$or_uci95,mr_resultsKFIC$or_uci95),
                     rep(c("IVW","MR-Egger","Weighted Median","Weighted Mode"),52),
                     rep(c(rep("CHD",4),rep("I.Stroke (small vessel)",4),rep("I.Stroke (large-artery)",4),rep("I.Stroke (cardioembolic)",4)),13),
                     c(rep("menarche",16),rep("BMI",16),rep("DBP",16),rep("SBP",16),rep("LDL",16),rep("HDL",16),rep("TG",16),
                       rep("PA",16),rep("DIA",16),rep("SMK",16),rep("SLP",16),rep("ALC",16),rep("KF",16)))

names(res_data)<-c("or","orlci_95","oruci_95","method","outcome","exposure")

resdat_ivw<-res_data[seq(1,208,4),]

resdat_ivw<-resdat_ivw[c(rev(seq(4,52,4)),rev(seq(3,52,4)),rev(seq(2,52,4)),rev(seq(1,52,4))),]

resdat_ivw$index<-c(1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.4,
                    4,4.2,4.4,4.6,4.8,5,5.2,5.4,5.6,5.8,6,6.2,6.4,
                    7,7.2,7.4,7.6,7.8,8,8.2,8.4,8.6,8.8,9,9.2,9.4,
                    10,10.2,10.4,10.6,10.8,11,11.2,11.4,11.6,11.8,12,12.2,12.4)

resdat_ivw$oruci_95[15]<-4.5 

uniplot1<-ggplot(resdat_ivw)+geom_point(aes(x=index, y=or))+
  coord_flip()+geom_errorbar(aes(ymin=orlci_95,ymax=oruci_95,x=index),width=0,orientation = "x")+
  scale_y_continuous(limits = c(-2,4.5),name = "Effect (odds ratio)")+geom_hline(yintercept=1,linetype = "dashed")+theme_bw()+
  scale_x_continuous(limits = c(1,13),breaks=c(3.8,6.8,9.8,12.8),labels= resdat_ivw$outcome[c(1,14,28,42)],name = "Outcome")+
  scale_color_discrete(name = "Exposure")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png("Figure3_1.png",width=6,height=10,units="in",res=200)
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

hetDBPCHD <- mr_heterogeneity(DBPCHD.dat)
hetDBPILA <- mr_heterogeneity(DBPIStkLA.dat)
hetDBPIC <-  mr_heterogeneity(DBPIStkC.dat)
hetDBPISV <- mr_heterogeneity(DBPIStkSV.dat)

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

hetPACHD <- mr_heterogeneity(PACHD.dat)
hetPAILA <- mr_heterogeneity(PAIStkLA.dat)
hetPAIC <- mr_heterogeneity(PAIStkC.dat)
hetPAISV <- mr_heterogeneity(PAIStkSV.dat)

hetDIACHD <- mr_heterogeneity(DIACHD.dat)
hetDIAILA <- mr_heterogeneity(DIAIStkLA.dat)
hetDIAIC <- mr_heterogeneity(DIAIStkC.dat)
hetDIAISV <- mr_heterogeneity(DIAIStkSV.dat)

hetSMKCHD <- mr_heterogeneity(SMKCHD.dat)
hetSMKILA <- mr_heterogeneity(SMKIStkLA.dat)
hetSMKIC <- mr_heterogeneity(SMKIStkC.dat)
hetSMKISV <- mr_heterogeneity(SMKIStkSV.dat)

hetSLPCHD <- mr_heterogeneity(SLPCHD.dat)
hetSLPILA <- mr_heterogeneity(SLPIStkLA.dat)
hetSLPIC <- mr_heterogeneity(SLPIStkC.dat)
hetSLPISV <- mr_heterogeneity(SLPIStkSV.dat)

hetALCCHD <- mr_heterogeneity(ALCCHD.dat)
hetALCILA <- mr_heterogeneity(ALCIStkLA.dat)
hetALCIC <- mr_heterogeneity(ALCIStkC.dat)
hetALCISV <- mr_heterogeneity(ALCIStkSV.dat)

hetKFCHD <- mr_heterogeneity(KFCHD.dat)
hetKFILA <- mr_heterogeneity(KFIStkLA.dat)
hetKFIC <- mr_heterogeneity(KFIStkC.dat)
hetKFISV <- mr_heterogeneity(KFIStkSV.dat)

############
### MVMR ###
############

#Extract all exposure data

exp_list<-c("ukb-a-315","ukb-b-19953","ukb-a-359","ukb-a-360","ieu-b-110","ieu-b-109","ieu-b-111","ukb-a-511","ieu-a-26",
            "ukb-b-2134","ukb-b-4424","ieu-a-1283","ukb-a-333")

output<-NULL

for(i in 1:13){
  
  reference<-exp_list[i]
  
  r.exposures<-exp_list[-i]
  
  A<-extract_instruments(reference,clump = F)
  
  for(j in 1:12){
    
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
    
    tseXGs<-matrix(c(HarmA_data$se.exposure,HarmB_data$se.exposure),ncol=2)
    
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

write.csv(output,"F-outs_new.csv",row.names=F)

F_matrix<-matrix(rep(0,169),ncol=13)
F_matrix[1,]<-c(0,output$exposure1[1:12])
F_matrix[2,]<-c(output$exposure1[13],0,output$exposure1[14:24])
F_matrix[3,]<-c(output$exposure1[25:26],0,output$exposure1[27:36])
F_matrix[4,]<-c(output$exposure1[37:39],0,output$exposure1[40:48])
F_matrix[5,]<-c(output$exposure1[49:52],0,output$exposure1[53:60])
F_matrix[6,]<-c(output$exposure1[61:65],0,output$exposure1[66:72])
F_matrix[7,]<-c(output$exposure1[73:78],0,output$exposure1[79:84])
F_matrix[8,]<-c(output$exposure1[85:91],0,output$exposure1[92:96])
F_matrix[9,]<-c(output$exposure1[97:104],0,output$exposure1[105:108])
F_matrix[10,]<-c(output$exposure1[109:117],0,output$exposure1[118:120])
F_matrix[11,]<-c(output$exposure1[121:130],0,output$exposure1[131:132])
F_matrix[12,]<-c(output$exposure1[133:143],0,output$exposure1[144:144])
F_matrix[13,]<-c(output$exposure1[145:156],0)

F_matrix2<-matrix(rep(0,169),ncol=13)
F_matrix2[1,]<-c(0,output$exposure2[1:12])
F_matrix2[2,]<-c(output$exposure2[13],0,output$exposure2[14:24])
F_matrix2[3,]<-c(output$exposure2[25:26],0,output$exposure2[27:36])
F_matrix2[4,]<-c(output$exposure2[37:39],0,output$exposure2[40:48])
F_matrix2[5,]<-c(output$exposure2[49:52],0,output$exposure2[53:60])
F_matrix2[6,]<-c(output$exposure2[61:65],0,output$exposure2[66:72])
F_matrix2[7,]<-c(output$exposure2[73:78],0,output$exposure2[79:84])
F_matrix2[8,]<-c(output$exposure2[85:91],0,output$exposure2[92:96])
F_matrix2[9,]<-c(output$exposure2[97:104],0,output$exposure2[105:108])
F_matrix2[10,]<-c(output$exposure2[109:117],0,output$exposure2[118:120])
F_matrix2[11,]<-c(output$exposure2[121:130],0,output$exposure2[131:132])
F_matrix2[12,]<-c(output$exposure2[133:143],0,output$exposure2[144:144])
F_matrix2[13,]<-c(output$exposure2[145:156],0)

meanF_matrix<-matrix(rep(0,169),ncol=13)

for(i in 1:13){
  for(j in 1:13){
    if(F_matrix[i,j] < 10){
      meanF_matrix[i,j] <- 0
    }
    if(F_matrix2[i,j] < 10){
      meanF_matrix[i,j] <- 0
    }
    if(F_matrix[i,j] > 10 & F_matrix2[i,j] > 10){
      meanF_matrix[i,j]<-mean(c(F_matrix[i,j],F_matrix[j,i]))*abs((1/(F_matrix[i,j]-F_matrix2[i,j])))
      
    }
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
                       midpoint = 0, limit = c(0,10), space = "Lab", 
                       name="Weighted\nconditional\nF-Statistic") +
  scale_x_continuous(breaks = 1:13,labels = c("Menarche","BMI","SBP","DBP","LDL","HDL","TG",
                                              "Physical activity","Type-II diabetes","Smoking",
                                              "Sleep duration","Alcohol consumption", "Creatinine"))+
  scale_y_continuous(breaks = 1:13,labels = c("Menarche","BMI","SBP","DBP","LDL","HDL","TG",
                                              "Physical activity","Type-II diabetes","Smoking",
                                              "Sleep duration","Alcohol consumption", "Creatinine"))+
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

png("Figure4_1.png",width=8,height=8,units="in",res=200)
Fstat_map
dev.off()


#HDL and SBP are strongest conditional pair - assessing additional

exp_list<-c("ukb-a-315","ukb-b-19953","ukb-a-360","ukb-a-359","ieu-b-110","ukb-a-511","ieu-a-26",
            "ukb-b-2134","ukb-b-4424","ieu-a-1283","ukb-a-333")

output2<-NULL

reference1<-"ieu-b-109"

reference2<-"ieu-b-111"

r.exposures<-exp_list

A<-extract_instruments(reference1,clump = F)

B<-extract_instruments(reference2,clump = F)

for(i in 1:11){
  
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
  
  HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmB_data$SNP,]
  
  HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmC_data$SNP,]
  
  HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmA_data$SNP,]
  
  HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmC_data$SNP,]
  
  HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmA_data$SNP,]
  
  HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmB_data$SNP,]
  
  HarmA_data<-HarmA_data[order(HarmA_data$SNP),]
  
  HarmB_data<-HarmB_data[order(HarmB_data$SNP),]
  
  HarmC_data<-HarmC_data[order(HarmC_data$SNP),]
  
  #Create combined dataframe for CHD MVMR
  tBXGs<-matrix(c(HarmA_data$beta.exposure,HarmB_data$beta.exposure,HarmC_data$beta.exposure),ncol=3)
  
  tseXGs<-matrix(c(HarmA_data$se.exposure,HarmB_data$se.exposure,HarmC_data$se.exposure),ncol=3)
  
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

write.csv(output2,"F-out2_new.csv",row.names=F)

#SBP HDL, and LDL are strongest - assessing additional exposures

exp_list<-c("ukb-a-315","ukb-b-19953","ukb-a-360","ukb-a-359","ukb-a-511","ieu-a-26",
            "ukb-b-2134","ukb-b-4424","ieu-a-1283","ukb-a-333")

output3<-NULL

reference1<-"ieu-b-109"

reference2<-"ieu-b-111"

reference3<-"ieu-b-110"

r.exposures<-exp_list

A<-extract_instruments(reference1,clump = F)

B<-extract_instruments(reference2,clump = F)

C<-extract_instruments(reference3,clump = F)

for(i in 1:10){
  
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
  SNPlistB<-refcomb2$SNP[refcomb2$SNP %in% SNPlistA]
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
  
  HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmB_data$SNP,]
  
  HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmC_data$SNP,]
  
  HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmD_data$SNP,]
  
  HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmA_data$SNP,]
  
  HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmC_data$SNP,]
  
  HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmD_data$SNP,]
  
  HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmA_data$SNP,]
  
  HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmB_data$SNP,]
  
  HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmD_data$SNP,]
  
  HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmA_data$SNP,]
  
  HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmB_data$SNP,]
  
  HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmC_data$SNP,]
  
  HarmA_data<-HarmA_data[order(HarmA_data$SNP),]
  
  HarmB_data<-HarmB_data[order(HarmB_data$SNP),]
  
  HarmC_data<-HarmC_data[order(HarmC_data$SNP),]
  
  HarmD_data<-HarmD_data[order(HarmD_data$SNP),]
  
  #Create combined dataframe for CHD MVMR
  tBXGs<-matrix(c(HarmA_data$beta.exposure,HarmB_data$beta.exposure,HarmC_data$beta.exposure,HarmD_data$beta.exposure),ncol=4)
  
  tseXGs<-matrix(c(HarmA_data$se.exposure,HarmB_data$se.exposure,HarmC_data$se.exposure,HarmD_data$se.exposure),ncol=4)
  
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

write.csv(output3,"F-out3_new.csv",row.names=F)

#LDL, HDL, TG, and SBP are strongest - assessing additional exposures

exp_list<-c("ukb-a-315","ukb-b-19953","ukb-a-359","ukb-a-511","ieu-a-26",
            "ukb-b-2134","ukb-b-4424","ieu-a-1283","ukb-a-333")

output4<-NULL

reference1<-"ieu-b-109"

reference2<-"ieu-b-111"

reference3<-"ieu-b-110"

reference4<-"ukb-a-360"

r.exposures<-exp_list

A<-extract_instruments(reference1,clump = F)

B<-extract_instruments(reference2,clump = F)

C<-extract_instruments(reference3,clump = F)

D<-extract_instruments(reference4,clump = F)

for(i in 1:9){
  
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
  SNPlistB<-refcomb2$SNP[refcomb2$SNP %in% SNPlistA]
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
    exposure_dat = E_clump, 
    outcome_dat = Outcome.data
  )
  
  HarmA_data<-HarmA_data[HarmA_data$mr_keep == T,]
  
  HarmB_data<-HarmB_data[HarmB_data$mr_keep == T,]
  
  HarmC_data<-HarmC_data[HarmC_data$mr_keep == T,]
  
  HarmD_data<-HarmD_data[HarmD_data$mr_keep == T,]
  
  HarmE_data<-HarmE_data[HarmE_data$mr_keep == T,]
  
  HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmB_data$SNP,]
  
  HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmC_data$SNP,]
  
  HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmD_data$SNP,]
  
  HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmE_data$SNP,]
  
  HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmA_data$SNP,]
  
  HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmC_data$SNP,]
  
  HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmD_data$SNP,]
  
  HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmE_data$SNP,]
  
  HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmA_data$SNP,]
  
  HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmB_data$SNP,]
  
  HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmD_data$SNP,]
  
  HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmE_data$SNP,]
  
  HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmA_data$SNP,]
  
  HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmB_data$SNP,]
  
  HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmC_data$SNP,]
  
  HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmE_data$SNP,]
  
  HarmE_data<-HarmE_data[HarmE_data$SNP %in% HarmA_data$SNP,]
  
  HarmE_data<-HarmE_data[HarmE_data$SNP %in% HarmB_data$SNP,]
  
  HarmE_data<-HarmE_data[HarmE_data$SNP %in% HarmC_data$SNP,]
  
  HarmE_data<-HarmE_data[HarmE_data$SNP %in% HarmD_data$SNP,]
  
  HarmA_data<-HarmA_data[order(HarmA_data$SNP),]
  
  HarmB_data<-HarmB_data[order(HarmB_data$SNP),]
  
  HarmC_data<-HarmC_data[order(HarmC_data$SNP),]
  
  HarmD_data<-HarmD_data[order(HarmD_data$SNP),]
  
  HarmE_data<-HarmE_data[order(HarmE_data$SNP),]
  
  #Create combined dataframe for CHD MVMR
  tBXGs<-matrix(c(HarmA_data$beta.exposure,HarmB_data$beta.exposure,HarmC_data$beta.exposure,HarmD_data$beta.exposure,HarmE_data$beta.exposure),ncol=5)
  
  tseXGs<-matrix(c(HarmA_data$se.exposure,HarmB_data$se.exposure,HarmC_data$se.exposure,HarmD_data$se.exposure,HarmE_data$se.exposure),ncol=5)
  
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

write.csv(output4,"F-out4_new.csv",row.names=F)


####################
### MVMR Effects ###
####################

# Obtain effect estimates for each outcome

HDL<-"ieu-b-109"

TG<-"ieu-b-111"

LDL<-"ieu-b-110"

SBP<-"ukb-a-360"


# Obtain effect estimates for each outcome

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
SNPlistB<-refcomb2$SNP[refcomb2$SNP %in% SNPlistA]
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

HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmB_data$SNP,]

HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmC_data$SNP,]

HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmD_data$SNP,]

HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmA_data$SNP,]

HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmC_data$SNP,]

HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmD_data$SNP,]

HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmA_data$SNP,]

HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmB_data$SNP,]

HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmD_data$SNP,]

HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmA_data$SNP,]

HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmB_data$SNP,]

HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmC_data$SNP,]

HarmA_data<-HarmA_data[order(HarmA_data$SNP),]

HarmB_data<-HarmB_data[order(HarmB_data$SNP),]

HarmC_data<-HarmC_data[order(HarmC_data$SNP),]

HarmD_data<-HarmD_data[order(HarmD_data$SNP),]

#Create combined dataframe for CHD MVMR
CHDBXGs<-matrix(c(HarmA_data$beta.exposure,HarmB_data$beta.exposure,HarmC_data$beta.exposure,HarmD_data$beta.exposure),ncol=4)

CHDseXGs<-matrix(c(HarmA_data$se.exposure,HarmB_data$se.exposure,HarmC_data$se.exposure,HarmD_data$se.exposure),ncol=4)

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


HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmB_data$SNP,]

HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmC_data$SNP,]

HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmD_data$SNP,]

HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmA_data$SNP,]

HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmC_data$SNP,]

HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmD_data$SNP,]

HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmA_data$SNP,]

HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmB_data$SNP,]

HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmD_data$SNP,]

HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmA_data$SNP,]

HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmB_data$SNP,]

HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmC_data$SNP,]

HarmA_data<-HarmA_data[order(HarmA_data$SNP),]

HarmB_data<-HarmB_data[order(HarmB_data$SNP),]

HarmC_data<-HarmC_data[order(HarmC_data$SNP),]

HarmD_data<-HarmD_data[order(HarmD_data$SNP),]

#Create combined dataframe for CHD MVMR
ILABXGs<-matrix(c(HarmA_data$beta.exposure,HarmB_data$beta.exposure,HarmC_data$beta.exposure,HarmD_data$beta.exposure),ncol=4)

ILAseXGs<-matrix(c(HarmA_data$se.exposure,HarmB_data$se.exposure,HarmC_data$se.exposure,HarmD_data$se.exposure),ncol=4)

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

HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmB_data$SNP,]

HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmC_data$SNP,]

HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmD_data$SNP,]

HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmA_data$SNP,]

HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmC_data$SNP,]

HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmD_data$SNP,]

HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmA_data$SNP,]

HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmB_data$SNP,]

HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmD_data$SNP,]

HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmA_data$SNP,]

HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmB_data$SNP,]

HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmC_data$SNP,]

HarmA_data<-HarmA_data[order(HarmA_data$SNP),]

HarmB_data<-HarmB_data[order(HarmB_data$SNP),]

HarmC_data<-HarmC_data[order(HarmC_data$SNP),]

HarmD_data<-HarmD_data[order(HarmD_data$SNP),]

#Create combined dataframe for CHD MVMR
ICBXGs<-matrix(c(HarmA_data$beta.exposure,HarmB_data$beta.exposure,HarmC_data$beta.exposure,HarmD_data$beta.exposure),ncol=4)

ICseXGs<-matrix(c(HarmA_data$se.exposure,HarmB_data$se.exposure,HarmC_data$se.exposure,HarmD_data$se.exposure),ncol=4)

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


HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmB_data$SNP,]

HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmC_data$SNP,]

HarmA_data<-HarmA_data[HarmA_data$SNP %in% HarmD_data$SNP,]

HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmA_data$SNP,]

HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmC_data$SNP,]

HarmB_data<-HarmB_data[HarmB_data$SNP %in% HarmD_data$SNP,]

HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmA_data$SNP,]

HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmB_data$SNP,]

HarmC_data<-HarmC_data[HarmC_data$SNP %in% HarmD_data$SNP,]

HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmA_data$SNP,]

HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmB_data$SNP,]

HarmD_data<-HarmD_data[HarmD_data$SNP %in% HarmC_data$SNP,]

HarmA_data<-HarmA_data[order(HarmA_data$SNP),]

HarmB_data<-HarmB_data[order(HarmB_data$SNP),]

HarmC_data<-HarmC_data[order(HarmC_data$SNP),]

HarmD_data<-HarmD_data[order(HarmD_data$SNP),]

#Create combined dataframe for CHD MVMR
ISVBXGs<-matrix(c(HarmA_data$beta.exposure,HarmB_data$beta.exposure,HarmC_data$beta.exposure,HarmD_data$beta.exposure),ncol=4)

ISVseXGs<-matrix(c(HarmA_data$se.exposure,HarmB_data$se.exposure,HarmC_data$se.exposure,HarmD_data$se.exposure),ncol=4)

ISVMVMR.dat<-format_mvmr(ISVBXGs,HarmA_data$beta.outcome,ISVseXGs,HarmA_data$se.outcome,HarmA_data$SNP)

mvmrresdatCHD<-ivw_mvmr(CHDMVMR.dat)
mvmrresdatILA<-ivw_mvmr(ILAMVMR.dat)
mvmrresdatIC<-ivw_mvmr(ICMVMR.dat)
mvmrresdatISV<-ivw_mvmr(ISVMVMR.dat)

#Plotting multivariable estimates

or<-exp(c(mvmrresdatCHD[1,1],mvmrresdatCHD[2,1],mvmrresdatCHD[3,1],mvmrresdatCHD[4,1],
          mvmrresdatISV[1,1],mvmrresdatISV[2,1],mvmrresdatISV[3,1],mvmrresdatISV[4,1],
          mvmrresdatILA[1,1],mvmrresdatILA[2,1],mvmrresdatILA[3,1],mvmrresdatILA[4,1],
          mvmrresdatIC[1,1],mvmrresdatIC[2,1],mvmrresdatIC[3,1],mvmrresdatIC[4,1]))

ses<-c(mvmrresdatCHD[1,2],mvmrresdatCHD[2,2],mvmrresdatCHD[3,2],mvmrresdatCHD[4,2],
       mvmrresdatISV[1,2],mvmrresdatISV[2,2],mvmrresdatISV[3,2],mvmrresdatISV[4,2],
       mvmrresdatILA[1,2],mvmrresdatILA[2,2],mvmrresdatILA[3,2],mvmrresdatILA[4,2],
       mvmrresdatIC[1,2],mvmrresdatIC[2,2],mvmrresdatIC[3,2],mvmrresdatIC[4,2])

lci<-or-1.96*ses

uci<-or+1.96*ses

exposure<-rep(c("LDL","HDL","TG","SBP"),4)

outcome<-c(rep("CHD",4),rep("I.Stroke (small vessel)",4),rep("I.Stroke (large-artery)",4),rep("I.Stroke (cardioembolic)",4))

mvmrres_data<-data.frame(or,lci,uci,outcome,exposure)

names(mvmrres_data)<-c("or","orlci_95","oruci_95","outcome","exposure")

mvmrres_data<-mvmrres_data[c(rev(1:4),rev(5:8),rev(9:12),rev(13:16)),]

mvmrres_data<-mvmrres_data[c(13:16,9:12,5:8,1:4),]

mvmrres_data$index<-c(seq(4,16,4),
                      seq(20,32,4),
                      seq(36,48,4),
                      seq(52,64,4))

LDL_CF<-mean(c(strength_mvmr(CHDMVMR.dat)[1,1],strength_mvmr(ILAMVMR.dat)[1,1],
               strength_mvmr(ICMVMR.dat)[1,1],strength_mvmr(ISVMVMR.dat)[1,1]))

HDL_CF<-mean(c(strength_mvmr(CHDMVMR.dat)[1,2],strength_mvmr(ILAMVMR.dat)[1,2],
               strength_mvmr(ICMVMR.dat)[1,2],strength_mvmr(ISVMVMR.dat)[1,2]))

TG_CF<-mean(c(strength_mvmr(CHDMVMR.dat)[1,3],strength_mvmr(ILAMVMR.dat)[1,3],
              strength_mvmr(ICMVMR.dat)[1,3],strength_mvmr(ISVMVMR.dat)[1,3]))

SBP_CF<-mean(c(strength_mvmr(CHDMVMR.dat)[1,4],strength_mvmr(ILAMVMR.dat)[1,4],
               strength_mvmr(ICMVMR.dat)[1,4],strength_mvmr(ISVMVMR.dat)[1,4]))

######

mvmrres_data$method<-rep("MVMR",16)

umr_res<-resdat_ivw[resdat_ivw$exposure %in% c("SBP","TG","HDL","LDL"),]

umr_res$method<-rep("IVW",16)

umr_res<-umr_res[c(4,1,2,3,8,5,6,7,12,9,10,11,16,13,14,15),]

umr_res$index<-c(seq(5,17,4),
                      seq(21,33,4),
                      seq(37,49,4),
                      seq(53,65,4))

mvmrres_data<-rbind(mvmrres_data,umr_res)




mvmrplot1<-ggplot(mvmrres_data)+geom_point(aes(x=index, y=or,shape=method))+
  coord_flip()+geom_errorbar(aes(ymin=orlci_95,ymax=oruci_95,x=index),width=0,orientation = "x")+
  scale_y_continuous(limits = c(-1,4),name = "Effect (odds ratio)")+geom_hline(yintercept=1,linetype = "dashed")+theme_bw()+
  scale_x_continuous(limits = c(2,68),breaks=rev(c(18,34,51,67)),labels= mvmrres_data$outcome[c(14,9,5,1)],name = "Outcome")+
  scale_color_discrete(name = "Exposure")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png("Figure5_1.png",width=6,height=6,units="in",res=200)
mvmrplot1
dev.off()
































