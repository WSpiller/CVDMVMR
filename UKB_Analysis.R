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

SNPS<-unique(c(exposure_data_Men$SNP,
               exposure_data_BMI$SNP,
               exposure_data_SBP$SNP,
               exposure_data_LDL$SNP,
               exposure_data_HDL$SNP,
               exposure_data_TG$SNP,
               exposure_data_PA$SNP,
               exposure_data_DIA$SNP,
               exposure_data_SMK$SNP,
               exposure_data_SLP$SNP,
               exposure_data_ALC$SNP,
               exposure_data_KF$SNP))

dat<-data.frame(SNPS)

names(dat)<-"SNP"

write.csv(dat,"sum_combine2.csv",row.names=F)

IS_incident<-read.csv("Summary_ISdat_v4_1_formatted.csv",header=T)
IS_incident<-na.omit(IS_incident)
mu<-1042/(10420)

##Menarche

outcome_data_Men<- extract_outcome_data(
  snps = exposure_data_Men$SNP, outcomes = "ukb-a-315",proxies=F)

IS_men<-IS_incident[IS_incident$SNP %in% exposure_data_Men$SNP,]

outcome_data_Men<-outcome_data_Men[outcome_data_Men$SNP %in% IS_men$SNP,]
IS_men<-IS_men[IS_men$SNP %in% outcome_data_Men$SNP,]
IS_men<-unique(IS_men)
IS_men<-IS_men[order(IS_men$SNP),]
outcome_data_Men<-outcome_data_Men[order(outcome_data_Men$SNP),]

outcome_data_Men$effect_allele.outcome <-IS_men$ALLELE1
outcome_data_Men$other_allele.outcome <-IS_men$ALLELE0
outcome_data_Men$beta.outcome<-IS_men$BETA/(mu*(1-mu))
outcome_data_Men$se.outcome<-IS_men$SE/(mu*(1-mu))
outcome_data_Men$eaf.outcome<-IS_men$A1FREQ

MenIS.dat <- harmonise_data(
  exposure_dat = exposure_data_Men, 
  outcome_dat = outcome_data_Men
)

mr_resultsMenIS <- generate_odds_ratios(mr(MenIS.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

##BMI

outcome_data_BMI<- extract_outcome_data(
  snps = exposure_data_BMI$SNP, outcomes = "ukb-b-19953",proxies=F)

IS_bmi<-IS_incident[IS_incident$SNP %in% exposure_data_BMI$SNP,]

outcome_data_BMI<-outcome_data_BMI[outcome_data_BMI$SNP %in% IS_bmi$SNP,]
IS_bmi<-IS_bmi[IS_bmi$SNP %in% outcome_data_BMI$SNP,]
IS_bmi<-unique(IS_bmi)
IS_bmi<-IS_bmi[order(IS_bmi$SNP),]
outcome_data_BMI<-outcome_data_BMI[order(outcome_data_BMI$SNP),]

outcome_data_BMI$effect_allele.outcome <-IS_bmi$ALLELE1
outcome_data_BMI$other_allele.outcome <-IS_bmi$ALLELE0
outcome_data_BMI$beta.outcome<-IS_bmi$BETA/(mu*(1-mu))
outcome_data_BMI$se.outcome<-IS_bmi$SE/(mu*(1-mu))
outcome_data_BMI$eaf.outcome<-IS_bmi$A1FREQ

BMIIS.dat <- harmonise_data(
  exposure_dat = exposure_data_BMI, 
  outcome_dat = outcome_data_BMI
)

mr_resultsBMIIS <- generate_odds_ratios(mr(BMIIS.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))






##DBP

outcome_data_DBP<- extract_outcome_data(
  snps = exposure_data_DBP$SNP, outcomes = "ukb-a-359",proxies=F)

IS_dbp<-IS_incident[IS_incident$SNP %in% exposure_data_DBP$SNP,]

outcome_data_DBP<-outcome_data_DBP[outcome_data_DBP$SNP %in% IS_dbp$SNP,]
IS_dbp<-IS_dbp[IS_dbp$SNP %in% outcome_data_DBP$SNP,]
IS_dbp<-unique(IS_dbp)
IS_dbp<-IS_dbp[order(IS_dbp$SNP),]
outcome_data_DBP<-outcome_data_DBP[order(outcome_data_DBP$SNP),]

outcome_data_DBP$effect_allele.outcome <-IS_dbp$ALLELE1
outcome_data_DBP$other_allele.outcome <-IS_dbp$ALLELE0
outcome_data_DBP$beta.outcome<-IS_dbp$BETA/(mu*(1-mu))
outcome_data_DBP$se.outcome<-IS_dbp$SE/(mu*(1-mu))
outcome_data_DBP$eaf.outcome<-IS_dbp$A1FREQ

DBPIS.dat <- harmonise_data(
  exposure_dat = exposure_data_DBP, 
  outcome_dat = outcome_data_DBP
)

mr_resultsDBPIS <- generate_odds_ratios(mr(DBPIS.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))


##SBP

outcome_data_SBP<- extract_outcome_data(
  snps = exposure_data_SBP$SNP, outcomes = "ukb-a-360",proxies=F)

IS_sbp<-IS_incident[IS_incident$SNP %in% exposure_data_SBP$SNP,]

outcome_data_SBP<-outcome_data_SBP[outcome_data_SBP$SNP %in% IS_sbp$SNP,]
IS_sbp<-IS_sbp[IS_sbp$SNP %in% outcome_data_SBP$SNP,]
IS_sbp<-unique(IS_sbp)
IS_sbp<-IS_sbp[order(IS_sbp$SNP),]
outcome_data_SBP<-outcome_data_SBP[order(outcome_data_SBP$SNP),]

outcome_data_SBP$effect_allele.outcome <-IS_sbp$ALLELE1
outcome_data_SBP$other_allele.outcome <-IS_sbp$ALLELE0
outcome_data_SBP$beta.outcome<-IS_sbp$BETA/(mu*(1-mu))
outcome_data_SBP$se.outcome<-IS_sbp$SE/(mu*(1-mu))
outcome_data_SBP$eaf.outcome<-IS_sbp$A1FREQ

SBPIS.dat <- harmonise_data(
  exposure_dat = exposure_data_SBP, 
  outcome_dat = outcome_data_SBP
)

mr_resultsSBPIS <- generate_odds_ratios(mr(SBPIS.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

##LDL

outcome_data_LDL<- extract_outcome_data(
  snps = exposure_data_LDL$SNP, outcomes = "ieu-b-110",proxies=F)

IS_ldl<-IS_incident[IS_incident$SNP %in% exposure_data_LDL$SNP,]

outcome_data_LDL<-outcome_data_LDL[outcome_data_LDL$SNP %in% IS_ldl$SNP,]
IS_ldl<-IS_ldl[IS_ldl$SNP %in% outcome_data_LDL$SNP,]
IS_ldl<-unique(IS_ldl)
IS_ldl<-IS_ldl[order(IS_ldl$SNP),]
outcome_data_LDL<-outcome_data_LDL[order(outcome_data_LDL$SNP),]

outcome_data_LDL$effect_allele.outcome <-IS_ldl$ALLELE1
outcome_data_LDL$other_allele.outcome <-IS_ldl$ALLELE0
outcome_data_LDL$beta.outcome<-IS_ldl$BETA/(mu*(1-mu))
outcome_data_LDL$se.outcome<-IS_ldl$SE/(mu*(1-mu))
outcome_data_LDL$eaf.outcome<-IS_ldl$A1FREQ

LDLIS.dat <- harmonise_data(
  exposure_dat = exposure_data_LDL, 
  outcome_dat = outcome_data_LDL
)

mr_resultsLDLIS <- generate_odds_ratios(mr(LDLIS.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

##HDL

outcome_data_HDL<- extract_outcome_data(
  snps = exposure_data_HDL$SNP, outcomes = "ieu-b-109",proxies=F)

IS_hdl<-IS_incident[IS_incident$SNP %in% exposure_data_HDL$SNP,]

outcome_data_HDL<-outcome_data_HDL[outcome_data_HDL$SNP %in% IS_hdl$SNP,]
IS_hdl<-IS_hdl[IS_hdl$SNP %in% outcome_data_HDL$SNP,]
IS_hdl<-unique(IS_hdl)
IS_hdl<-IS_hdl[order(IS_hdl$SNP),]
outcome_data_HDL<-outcome_data_HDL[order(outcome_data_HDL$SNP),]

outcome_data_HDL$effect_allele.outcome <-IS_hdl$ALLELE1
outcome_data_HDL$other_allele.outcome <-IS_hdl$ALLELE0
outcome_data_HDL$beta.outcome<-IS_hdl$BETA/(mu*(1-mu))
outcome_data_HDL$se.outcome<-IS_hdl$SE/(mu*(1-mu))
outcome_data_HDL$eaf.outcome<-IS_hdl$A1FREQ

HDLIS.dat <- harmonise_data(
  exposure_dat = exposure_data_HDL, 
  outcome_dat = outcome_data_HDL
)

mr_resultsHDLIS <- generate_odds_ratios(mr(HDLIS.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

##TG

outcome_data_TG<- extract_outcome_data(
  snps = exposure_data_TG$SNP, outcomes = "ieu-b-111",proxies=F)

outcome_data_TG<-outcome_data_TG[!duplicated(outcome_data_TG$SNP),]

IS_tg<-IS_incident[IS_incident$SNP %in% outcome_data_TG$SNP,]

IS_tg<-IS_tg[!duplicated(IS_tg$SNP),]

outcome_data_TG<-outcome_data_TG[outcome_data_TG$SNP %in% IS_tg$SNP,]

IS_tg<-IS_tg[IS_tg$SNP %in% outcome_data_TG$SNP,]
IS_tg<-unique(IS_tg)
IS_tg<-IS_tg[order(IS_tg$SNP),]
outcome_data_TG<-outcome_data_TG[order(outcome_data_TG$SNP),]


outcome_data_TG$effect_allele.outcome <-IS_tg$ALLELE1
outcome_data_TG$other_allele.outcome <-IS_tg$ALLELE0
outcome_data_TG$beta.outcome<-IS_tg$BETA/(mu*(1-mu))
outcome_data_TG$se.outcome<-IS_tg$SE/(mu*(1-mu))
outcome_data_TG$eaf.outcome<-IS_tg$A1FREQ

TGIS.dat <- harmonise_data(
  exposure_dat = exposure_data_TG, 
  outcome_dat = outcome_data_TG
)

mr_resultsTGIS <- generate_odds_ratios(mr(TGIS.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

##PA

outcome_data_PA<- extract_outcome_data(
  snps = exposure_data_PA$SNP, outcomes = "ukb-a-511",proxies=F)

outcome_data_PA<-outcome_data_PA[!duplicated(outcome_data_PA$SNP),]

IS_pa<-IS_incident[IS_incident$SNP %in% outcome_data_PA$SNP,]

IS_pa<-IS_pa[!duplicated(IS_pa$SNP),]

outcome_data_PA<-outcome_data_PA[outcome_data_PA$SNP %in% IS_pa$SNP,]

IS_pa<-IS_pa[IS_pa$SNP %in% outcome_data_PA$SNP,]
IS_pa<-unique(IS_pa)
IS_pa<-IS_pa[order(IS_pa$SNP),]
outcome_data_PA<-outcome_data_PA[order(outcome_data_PA$SNP),]


outcome_data_PA$effect_allele.outcome <-IS_pa$ALLELE1
outcome_data_PA$other_allele.outcome <-IS_pa$ALLELE0
outcome_data_PA$beta.outcome<-IS_pa$BETA/(mu*(1-mu))
outcome_data_PA$se.outcome<-IS_pa$SE/(mu*(1-mu))
outcome_data_PA$eaf.outcome<-IS_pa$A1FREQ

PAIS.dat <- harmonise_data(
  exposure_dat = exposure_data_PA, 
  outcome_dat = outcome_data_PA
)

mr_resultsPAIS <- generate_odds_ratios(mr(PAIS.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

##DIA

outcome_data_DIA<- extract_outcome_data(
  snps = exposure_data_DIA$SNP, outcomes = "ieu-a-26",proxies=F)

outcome_data_DIA<-outcome_data_DIA[!duplicated(outcome_data_DIA$SNP),]

IS_dia<-IS_incident[IS_incident$SNP %in% outcome_data_DIA$SNP,]

IS_dia<-IS_dia[!duplicated(IS_dia$SNP),]

outcome_data_DIA<-outcome_data_DIA[outcome_data_DIA$SNP %in% IS_dia$SNP,]

IS_dia<-IS_dia[IS_dia$SNP %in% outcome_data_DIA$SNP,]
IS_dia<-unique(IS_dia)
IS_dia<-IS_dia[order(IS_dia$SNP),]
outcome_data_DIA<-outcome_data_DIA[order(outcome_data_DIA$SNP),]


outcome_data_DIA$effect_allele.outcome <-IS_dia$ALLELE1
outcome_data_DIA$other_allele.outcome <-IS_dia$ALLELE0
outcome_data_DIA$beta.outcome<-IS_dia$BETA/(mu*(1-mu))
outcome_data_DIA$se.outcome<-IS_dia$SE/(mu*(1-mu))
outcome_data_DIA$eaf.outcome<-IS_dia$A1FREQ

DIAIS.dat <- harmonise_data(
  exposure_dat = exposure_data_DIA, 
  outcome_dat = outcome_data_DIA
)

mr_resultsDIAIS <- generate_odds_ratios(mr(DIAIS.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

##SMK

outcome_data_SMK<- extract_outcome_data(
  snps = exposure_data_SMK$SNP, outcomes = "ukb-b-2134",proxies=F)

outcome_data_SMK<-outcome_data_SMK[!duplicated(outcome_data_SMK$SNP),]

IS_smk<-IS_incident[IS_incident$SNP %in% outcome_data_SMK$SNP,]

IS_smk<-IS_smk[!duplicated(IS_smk$SNP),]

outcome_data_SMK<-outcome_data_SMK[outcome_data_SMK$SNP %in% IS_smk$SNP,]

IS_smk<-IS_smk[IS_smk$SNP %in% outcome_data_SMK$SNP,]
IS_smk<-unique(IS_smk)
IS_smk<-IS_smk[order(IS_smk$SNP),]
outcome_data_SMK<-outcome_data_SMK[order(outcome_data_SMK$SNP),]


outcome_data_SMK$effect_allele.outcome <-IS_smk$ALLELE1
outcome_data_SMK$other_allele.outcome <-IS_smk$ALLELE0
outcome_data_SMK$beta.outcome<-IS_smk$BETA/(mu*(1-mu))
outcome_data_SMK$se.outcome<-IS_smk$SE/(mu*(1-mu))
outcome_data_SMK$eaf.outcome<-IS_smk$A1FREQ

SMKIS.dat <- harmonise_data(
  exposure_dat = exposure_data_SMK, 
  outcome_dat = outcome_data_SMK
)

mr_resultsSMKIS <- generate_odds_ratios(mr(SMKIS.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

##SLP

outcome_data_SLP<- extract_outcome_data(
  snps = exposure_data_SLP$SNP, outcomes = "ukb-b-4424",proxies=F)

outcome_data_SLP<-outcome_data_SLP[!duplicated(outcome_data_SLP$SNP),]

IS_slp<-IS_incident[IS_incident$SNP %in% outcome_data_SLP$SNP,]

IS_slp<-IS_slp[!duplicated(IS_slp$SNP),]

outcome_data_SLP<-outcome_data_SLP[outcome_data_SLP$SNP %in% IS_slp$SNP,]

IS_slp<-IS_slp[IS_slp$SNP %in% outcome_data_SLP$SNP,]
IS_slp<-unique(IS_slp)
IS_slp<-IS_slp[order(IS_slp$SNP),]
outcome_data_SLP<-outcome_data_SLP[order(outcome_data_SLP$SNP),]


outcome_data_SLP$effect_allele.outcome <-IS_slp$ALLELE1
outcome_data_SLP$other_allele.outcome <-IS_slp$ALLELE0
outcome_data_SLP$beta.outcome<-IS_slp$BETA/(mu*(1-mu))
outcome_data_SLP$se.outcome<-IS_slp$SE/(mu*(1-mu))
outcome_data_SLP$eaf.outcome<-IS_slp$A1FREQ

SLPIS.dat <- harmonise_data(
  exposure_dat = exposure_data_SLP, 
  outcome_dat = outcome_data_SLP
)

mr_resultsSLPIS <- generate_odds_ratios(mr(SLPIS.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

##ALC

outcome_data_ALC<- extract_outcome_data(
  snps = exposure_data_ALC$SNP, outcomes = "ieu-a-1283",proxies=F)

outcome_data_ALC<-outcome_data_ALC[!duplicated(outcome_data_ALC$SNP),]

IS_alc<-IS_incident[IS_incident$SNP %in% outcome_data_ALC$SNP,]

IS_alc<-IS_alc[!duplicated(IS_alc$SNP),]

outcome_data_ALC<-outcome_data_ALC[outcome_data_ALC$SNP %in% IS_alc$SNP,]

IS_alc<-IS_alc[IS_alc$SNP %in% outcome_data_ALC$SNP,]
IS_alc<-unique(IS_alc)
IS_alc<-IS_alc[order(IS_alc$SNP),]
outcome_data_ALC<-outcome_data_ALC[order(outcome_data_ALC$SNP),]


outcome_data_ALC$effect_allele.outcome <-IS_alc$ALLELE1
outcome_data_ALC$other_allele.outcome <-IS_alc$ALLELE0
outcome_data_ALC$beta.outcome<-IS_alc$BETA/(mu*(1-mu))
outcome_data_ALC$se.outcome<-IS_alc$SE/(mu*(1-mu))
outcome_data_ALC$eaf.outcome<-IS_alc$A1FREQ

ALCIS.dat <- harmonise_data(
  exposure_dat = exposure_data_ALC, 
  outcome_dat = outcome_data_ALC
)

mr_resultsALCIS <- generate_odds_ratios(mr(ALCIS.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))

##KF

outcome_data_KF<- extract_outcome_data(
  snps = exposure_data_KF$SNP, outcomes = "ukb-a-333",proxies=F)

outcome_data_KF<-outcome_data_KF[!duplicated(outcome_data_KF$SNP),]

IS_kf<-IS_incident[IS_incident$SNP %in% outcome_data_KF$SNP,]

IS_kf<-IS_kf[!duplicated(IS_kf$SNP),]

outcome_data_KF<-outcome_data_KF[outcome_data_KF$SNP %in% IS_kf$SNP,]

IS_kf<-IS_kf[IS_kf$SNP %in% outcome_data_KF$SNP,]
IS_kf<-unique(IS_kf)
IS_kf<-IS_kf[order(IS_kf$SNP),]
outcome_data_KF<-outcome_data_KF[order(outcome_data_KF$SNP),]


outcome_data_KF$effect_allele.outcome <-IS_kf$ALLELE1
outcome_data_KF$other_allele.outcome <-IS_kf$ALLELE0
outcome_data_KF$beta.outcome<-IS_kf$BETA/(mu*(1-mu))
outcome_data_KF$se.outcome<-IS_kf$SE/(mu*(1-mu))
outcome_data_KF$eaf.outcome<-IS_kf$A1FREQ

KFIS.dat <- harmonise_data(
  exposure_dat = exposure_data_KF, 
  outcome_dat = outcome_data_KF
)

mr_resultsKFIS <- generate_odds_ratios(mr(KFIS.dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode")))



























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

################
### I_Stroke ###
################

SNPS<-unique(c(D_clump$SNP))

dat<-data.frame(SNPS)

names(dat)<-"SNP"

write.csv(dat,"sum_combine3.csv",row.names=F)

IS_incident<-read.csv("Summary_ISdat_v4_2_formatted.csv",header=T)

IS_incident<-na.omit(IS_incident)
mu<-1042/(10420)

##Menarche

Outcome.data<- extract_outcome_data(
  snps = D_clump$SNP, outcomes = LDL,proxies=F)

IS_MVMR<-IS_incident[IS_incident$SNP %in% D_clump$SNP,]

Outcome.data<-Outcome.data[Outcome.data$SNP %in% IS_MVMR$SNP,]
IS_MVMR<-IS_MVMR[IS_MVMR$SNP %in% Outcome.data$SNP,]
IS_MVMR<-unique(IS_MVMR)
IS_MVMR<-IS_MVMR[order(IS_MVMR$SNP),]
Outcome.data<-Outcome.data[order(Outcome.data$SNP),]

Outcome.data$effect_allele.outcome <-IS_MVMR$ALLELE1
Outcome.data$other_allele.outcome <-IS_MVMR$ALLELE0
Outcome.data$beta.outcome<-IS_MVMR$BETA/(mu*(1-mu))
Outcome.data$se.outcome<-IS_MVMR$SE/(mu*(1-mu))
Outcome.data$eaf.outcome<-IS_MVMR$A1FREQ

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

MVMR_res<-ivw_mvmr(CHDMVMR.dat)

exp(MVMR_res[,1])

exp(MVMR_res[,1]) - 1.96 * MVMR_res[,2]

exp(MVMR_res[,1]) + 1.96 * MVMR_res[,2]


# Plotting univariable estimates

res_data<-data.frame(c(mr_resultsMenIS$or,mr_resultsBMIIS$or,mr_resultsDBPIS$or,mr_resultsSBPIS$or,
                       mr_resultsLDLIS$or,mr_resultsHDLIS$or,mr_resultsTGIS$or,mr_resultsPAIS$or,
                       mr_resultsDIAIS$or,mr_resultsSMKIS$or,mr_resultsSLPIS$or,mr_resultsALCIS$or,
                       mr_resultsKFIS$or),
                     c(mr_resultsMenIS$or_lci95,mr_resultsBMIIS$or_lci95,mr_resultsDBPIS$or_lci95,mr_resultsSBPIS$or_lci95,
                       mr_resultsLDLIS$or_lci95,mr_resultsHDLIS$or_lci95,mr_resultsTGIS$or_lci95,mr_resultsPAIS$or_lci95,
                       mr_resultsDIAIS$or_lci95,mr_resultsSMKIS$or_lci95,mr_resultsSLPIS$or_lci95,mr_resultsALCIS$or_lci95,
                       mr_resultsKFIS$or_lci95),
                     c(mr_resultsMenIS$or_uci95,mr_resultsBMIIS$or_uci95,mr_resultsDBPIS$or_uci95,mr_resultsSBPIS$or_uci95,
                       mr_resultsLDLIS$or_uci95,mr_resultsHDLIS$or_uci95,mr_resultsTGIS$or_uci95,mr_resultsPAIS$or_uci95,
                       mr_resultsDIAIS$or_uci95,mr_resultsSMKIS$or_uci95,mr_resultsSLPIS$or_uci95,mr_resultsALCIS$or_uci95,
                       mr_resultsKFIS$or_uci95),
                     rep(c("IVW","MR-Egger","Weighted Median","Weighted Mode"),13),
                     c(rep("menarche",4),rep("BMI",4),rep("DBP",4),rep("SBP",4),rep("LDL",4),rep("HDL",4),rep("TG",4),
                       rep("PA",4),rep("DIA",4),rep("SMK",4),rep("SLP",4),rep("ALC",4),rep("KF",4)))

names(res_data)<-c("or","orlci_95","oruci_95","method","exposure")

resdat_ivw<-res_data[seq(1,52,4),]

res_mvmr<-data.frame(exp(MVMR_res[,1]),(exp(MVMR_res[,1]) - (1.96*MVMR_res[,2])),(exp(MVMR_res[,1]) + (1.96*MVMR_res[,2])),
                     rep("MVMR",4),c("LDL","HDL","TG","SBP"))

names(res_mvmr)<-names(resdat_ivw)

resdat_ivw<-rbind(resdat_ivw,res_mvmr)

resdat_ivw$index<-1:17

resdat_ivw$index[1:13]<-c(seq(1,13,1))*-1

resdat_ivw$index[14:17]<-c(5.2,6.2,7.2,4.2)*-1


uniplot1<-ggplot(resdat_ivw)+geom_point(aes(x=index, y=or,shape=method))+
  coord_flip()+geom_errorbar(aes(ymin=orlci_95,ymax=oruci_95,x=index),width=0,orientation = "x")+
  scale_y_continuous(limits = c(0,3),name = "Effect (odds ratio)")+geom_hline(yintercept=1,linetype = "dashed")+theme_bw()+
  scale_color_discrete(name = "Exposure")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

png("Figure6_1.png",width=5,height=5,units="in",res=200)
uniplot1
dev.off()

















































































