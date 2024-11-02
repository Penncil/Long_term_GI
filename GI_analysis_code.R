
library(dplyr)
library(readr)
library(tibble)
library(janitor)
library(pracma)
library(stringr)
library(data.table)

# Packages to connect our database
library(srcr)
library(RPostgres)
library(dbplyr)

# Packages for analysis
#library(Cyclops)
#library(CohortMethod)
library(geex)
library(glmnet)
library(survival)

# SMD and random effect meta
library(tableone)
library(cobalt)
library(metafor)

# Plot
library(ggplot2)
library(ggpubr)
#library(wpa)
#library(stringr)

#############################################
## Source functions for vaccine ##
## Updated @ 11-11-2023 to use Poisson regression with offset
## and use integration for conditional Poisson
setwd("~/Project1")
source("~/Project1/site/run.R")
source("~/Project1/site/site_info.R")
source("~/Project1/analysis/source_functions_GI_07112023.R")

setwd("~/Project1/results/GI_11132023")


case<-results_tbl("round2_s10_case_GI_final_dz")
control<-results_tbl("round2_s10_cont_GI_final_dz")

case_data<-case%>%collect_new()
cont_neg_data<-control%>%collect_new()
case_data<-case_data%>%data.frame()
cont_neg_data<-cont_neg_data%>%data.frame()
case_data$treatment=1
cont_neg_data$treatment=0
case_data=case_data[,names(case_data)%in%names(cont_neg_data)]
cont_neg_data=cont_neg_data[,names(cont_neg_data)%in%names(case_data)]
data=rbind(case_data,cont_neg_data)
rm(case)
rm(control)
rm(case_data)
rm(cont_neg_data)


data=clean_names(data)

visits_signs_list=c("abdominal_pain","bloating","constipation","diarrhea",
                    "nausea","vomiting")
visits_disease_list=c("gastroesophageal_reflux_disease")
visits_chronic_disease_list=c("gastroesophageal_reflux_disease","irritable_bowel_syndrome","functional_dyspepsia")


data$visits_signs<- ifelse(rowSums(data%>%select(paste0('visits_',visits_signs_list)))>0,1,0)
data$base_visits_signs <- ifelse(rowSums(data%>%select(paste0('base_visits_',visits_signs_list)))>0,1,0)
data$chronic_visits_signs<- ifelse(rowSums(data%>%
                                             select(paste0('chronic_visits_',visits_signs_list)))>0,1,0)

data$visits_disease<- ifelse(rowSums(data%>%
                                       select(paste0('visits_',visits_disease_list)))>0,1,0)
data$base_visits_disease <- ifelse(rowSums(data%>%
                                             select(paste0('base_visits_',visits_disease_list)))>0,1,0)

data$chronic_visits_disease<- ifelse(rowSums(data%>%
                                               select(paste0('chronic_visits_',visits_chronic_disease_list)))>0,1,0)
data$base_visits_chronic_disease<- ifelse(rowSums(data%>%
                                                    select(paste0('base_visits_',visits_chronic_disease_list)))>0,1,0)


data$visits_all<- ifelse(rowSums(data%>%
                                   select(paste0('visits_',visits_signs_list),paste0('visits_',visits_disease_list)))>0,1,0)
data$base_visits_all<- ifelse(rowSums(data%>%
                                        select(paste0('base_visits_',visits_signs_list),paste0('base_visits_',visits_disease_list)))>0,1,0)
data$chronic_visits_all<- ifelse(rowSums(data%>%
                                           select(paste0('chronic_visits_',visits_signs_list),paste0('chronic_visits_',visits_chronic_disease_list)))>0,1,0)
data$base_visits_chronic_all<-ifelse(rowSums(data%>%
                                               select(paste0('base_visits_',visits_signs_list),paste0('base_visits_',visits_chronic_disease_list)))>0,1,0)




#data is prepared 
#the rest would be different depends on different subgroup
colnames(data) = make_clean_names(colnames(data))

total_list=c(visits_signs_list,visits_chronic_disease_list,
             "signs","disease","all","chronic_disease","chronic_all")

#data structure being used
nstrata = 6
incidence_mat=matrix(NA,nrow=length(total_list),ncol=3)
chronic_incidence_mat=matrix(NA,nrow=length(total_list),ncol=3)
strata_rslt_mat=matrix(NA,nrow=length(total_list),ncol=8)
chronic_strata_rslt_mat=matrix(NA,nrow=length(total_list),ncol=8)



############################
##step 1: propensity score##
############################
for(i in 1:length(total_list)){
  print(i)
  
  selection_statement_1=(paste0("base_visits_",total_list[i],"==0"))
  mydata=data%>%filter(rlang::eval_tidy(rlang::parse_expr(selection_statement_1)))
  
  
  mydata_medical <- mydata%>%select(starts_with('medical_'))
  cond1 <- colSums(mydata_medical)<0.001*nrow(mydata_medical)
  colnames(mydata_medical)[cond1] # low incidence chronic conditions
  mydata <- mydata %>% select(-colnames(mydata_medical)[cond1])
  
  colnames(mydata) = make_clean_names(colnames(mydata))
  
  # extracted variables needed for the analysis
  #select varaibles （need to change
  xvars = c("site","entry_age", "sex_cat", "eth_cat", "obese",
            "n_ed","n_inpatient","n_outpatient","n_drug","total_negative_tests_cat",
            "cohort_entry_month", "pmca_index", colnames(mydata%>%select(starts_with('medical_'))))
  
  
  mydata_ps <- mydata[,colnames(mydata)%in% c(xvars, "treatment")]
  mydata_ps$n_ed[is.na(mydata_ps$n_ed)]<- rep(0)
  mydata_ps$n_ed[mydata_ps$n_ed>2]<- rep("2+")
  
  mydata_ps$n_inpatient[is.na(mydata_ps$n_inpatient)]<- rep(0)
  mydata_ps$n_inpatient[mydata_ps$n_inpatient>2]<- rep("2+")
  
  mydata_ps$n_outpatient[is.na(mydata_ps$n_outpatient)]<- rep(0)
  mydata_ps$n_outpatient[mydata_ps$n_outpatient>2]<- rep("2+")
  
  mydata_ps$n_drug[is.na(mydata_ps$n_drug)]<- rep(0)
  mydata_ps$n_drug[mydata_ps$n_drug>2]<- rep("2+")
  
  mydata_ps <- as.data.frame(mydata_ps)
  for(i_cov in 1:ncol(mydata_ps)){
    if(is.character(mydata_ps[,i_cov]))
      mydata_ps[,i_cov] <- factor(mydata_ps[,i_cov]) 
  }
  
  
  ###################################################################################
  ######          PS model            #######
  ###################################################################################
  
  
  #############    Step 1: Fit PS model   ###########
  form <- as.formula(paste("treatment ~ ", paste(xvars, collapse= "+"))) # model formula
  Xmat <- grab_design_matrix(data = mydata_ps, rhs_formula = form)
  Y <- mydata_ps$treatment
  nfolds = 10
  set.seed(2022)
  foldid = sample(rep(seq(nfolds), length.out = length(Y)))
  Fit_ps_cv <- cv.glmnet(Xmat, Y, alpha = 1, family = "binomial", nfolds = nfolds, foldid = foldid) # cross validation to select lambda
  Fit_ps <- glmnet(Xmat, Y, alpha = 1, family = "binomial", lambda = Fit_ps_cv$lambda.min)
  propensityScore <- predict(Fit_ps, Xmat, type = "response")
  #need to change the save function
  save(propensityScore, file =paste0(total_list[i], "_propensityScore_11132023.rda"))
  gc()
  
}








############################
##step 2: analysis    231-395 ##
############################
for(i in 1:length(total_list)){
  print(i)
  
  selection_statement_1=(paste0("base_visits_",total_list[i],"==0"))
  mydata=data%>%filter(rlang::eval_tidy(rlang::parse_expr(selection_statement_1)))
  
  
  mydata_medical <- mydata%>%select(starts_with('medical_'))
  cond1 <- colSums(mydata_medical)<0.001*nrow(mydata_medical)
  colnames(mydata_medical)[cond1] # low incidence chronic conditions
  mydata <- mydata %>% select(-colnames(mydata_medical)[cond1])
  
  colnames(mydata) = make_clean_names(colnames(mydata))
  
  # extracted variables needed for the analysis
  #select varaibles （need to change
  xvars = c("site","entry_age", "sex_cat", "eth_cat", "obese",
            "n_ed","n_inpatient","n_outpatient","n_drug","total_negative_tests_cat",
            "cohort_entry_month", "pmca_index", colnames(mydata%>%select(starts_with('medical_'))))
  
  
  mydata_ps <- mydata[,colnames(mydata)%in% c(xvars, "treatment")]
  mydata_ps$n_ed[is.na(mydata_ps$n_ed)]<- rep(0)
  mydata_ps$n_ed[mydata_ps$n_ed>2]<- rep("2+")
  
  mydata_ps$n_inpatient[is.na(mydata_ps$n_inpatient)]<- rep(0)
  mydata_ps$n_inpatient[mydata_ps$n_inpatient>2]<- rep("2+")
  
  mydata_ps$n_outpatient[is.na(mydata_ps$n_outpatient)]<- rep(0)
  mydata_ps$n_outpatient[mydata_ps$n_outpatient>2]<- rep("2+")
  
  mydata_ps$n_drug[is.na(mydata_ps$n_drug)]<- rep(0)
  mydata_ps$n_drug[mydata_ps$n_drug>2]<- rep("2+")
  
  mydata_ps <- as.data.frame(mydata_ps)
  for(i_cov in 1:ncol(mydata_ps)){
    if(is.character(mydata_ps[,i_cov]))
      mydata_ps[,i_cov] <- factor(mydata_ps[,i_cov]) 
  }
  
  
  load(paste0(total_list[i], "_propensityScore_11132023.rda"))
  mydata_ps$propensityScore <- propensityScore
  
  
  person_id=mydata$person_id
  mydata_ps=cbind(person_id,mydata_ps)
  #############   Step 2. Evaluate empirical equipoise      #########################
  # Based on each possible value of se, we can calculate the preference scores 
  equipoiseBounds = c(0.3, 0.7)
  equipoise.vec <- c()
  #mydata_ps$treatment=mydata_ps$tr
  # preference score
  mydata_ps_new <- mydata_ps
  mydata_ps_new$propensityScore = mydata_ps_new$propensityScore
  mydata_ps_new$propensityScore[mydata_ps_new$propensityScore > 0.9999999] <- 0.9999999
  mydata_ps_new <- computePreferenceScore(mydata_ps_new)
  equipoise <- mean(mydata_ps_new$preferenceScore >= equipoiseBounds[1] & mydata_ps_new$preferenceScore <= equipoiseBounds[2])
  equipoise.vec <-  equipoise
  
  
  #mydata_ps_new$treatment=mydata_ps_new$tr
  # plot empirical equipoise when se = 1 (or ignore misclassification)
  equipoise.plot <- plotPs_mis(mydata_ps_new, scale = "preference",
                               showCountsLabel = TRUE, showAucLabel = TRUE, showEquiposeLabel = TRUE)
  #need to change the save function
  ggsave(equipoise.plot, file=paste0(total_list[i],"_equipoise_08142023.png"), width=6, height=6)
  
  rm(mydata_ps_new)
  #need to change the save function
  save.equipoise <- data.frame(equipoise = equipoise.vec)
  save(save.equipoise, file =paste0(total_list[i], "_equipoise_08142023.rda"))
  rm(save.equipoise)
  ###################################################################################
  ######        Step 3. Stratification,6 strata              #######
  ###################################################################################
  
  #period
  period=(as.Date("2022-12-03")-mydata$cohort_entry_date)/(365.25)
  period=as.numeric(period)
  mydata_ps$period=period
  
  
  
  if(!total_list[i]%in%c(
    "chronic_disease","chronic_all"
  )){
    mydata_test=mydata_ps%>%dplyr::left_join(
      mydata%>%dplyr::select(person_id,paste0("visits_",total_list[i]),paste0("chronic_visits_",total_list[i])),
      by="person_id"
    )
    rm(mydata)
    rm(mydata_ps)
    
    
    mydata_test_case=mydata_test%>%filter(treatment==1)
    mydata_test_control=mydata_test%>%filter(treatment==0)
    
    incidence_mat[i,1]=total_list[i]
    incidence_mat[i,2]=table(mydata_test_case[,paste0("visits_",total_list[i])])[2]/length(mydata_test_case[,paste0("visits_",total_list[i])])
    incidence_mat[i,3]=table(mydata_test_control[,paste0("visits_",total_list[i])])[2]/length(mydata_test_control[,paste0("visits_",total_list[i])])
    chronic_incidence_mat[i,1]=paste0("chronic_visits_",total_list[i])
    chronic_incidence_mat[i,2]=table(mydata_test_case[,paste0("chronic_visits_",total_list[i])])[2]/length(mydata_test_case[,paste0("chronic_visits_",total_list[i])])
    chronic_incidence_mat[i,3]=table(mydata_test_control[,paste0("chronic_visits_",total_list[i])])[2]/length(mydata_test_control[,paste0("chronic_visits_",total_list[i])])
    rm(mydata_test_case)
    rm(mydata_test_control)
    res_strata=strata_analysis(nstrata,mydata_test,xvars)
    strata_rslt_mat[i,]=res_strata$rslt
    
    balance_table = cbind(res_strata$smd_before, res_strata$smd_after[,2])
    #change names do it afterwards in local computer
    #smd_before=balance_table[,1:2]
    #change_names(smd_before)
    save(balance_table,file =paste0(total_list[i], "balance_table.rda"))
    #adjusted_rslt_mat[i,]=adjusted_analysis(mydata_test,xvars)
    chronic_strata_rslt_mat[i,]=chronic_strata_analysis(nstrata,mydata_test)
    
    
  }else{
    
    tmp=gsub("chronic_","",total_list[i])
    
    mydata_test=mydata_ps%>%dplyr::left_join(
      mydata%>%dplyr::select(person_id,paste0("visits_",tmp),paste0("chronic_visits_",tmp)),
      by="person_id"
    )
    rm(mydata)
    rm(mydata_ps)
    
    
    mydata_test_case=mydata_test%>%filter(treatment==1)
    mydata_test_control=mydata_test%>%filter(treatment==0)
    
    chronic_incidence_mat[i,1]=paste0("chronic_visits_",tmp)
    chronic_incidence_mat[i,2]=table(mydata_test_case[,paste0("chronic_visits_",tmp)])[2]/length(mydata_test_case[,paste0("chronic_visits_",tmp)])
    chronic_incidence_mat[i,3]=table(mydata_test_control[,paste0("chronic_visits_",tmp)])[2]/length(mydata_test_control[,paste0("chronic_visits_",tmp)])
    rm(mydata_test_case)
    rm(mydata_test_control)
    res_strata=strata_analysis(nstrata,mydata_test,xvars)
    strata_rslt_mat[i,]=res_strata$rslt
    
    balance_table = cbind(res_strata$smd_before, res_strata$smd_after[,2])
    #change names do it afterwards in local computer
    #smd_before=balance_table[,1:2]
    #change_names(smd_before)
    save(balance_table,file =paste0(total_list[i], "balance_table.rda"))
    #adjusted_rslt_mat[i,]=adjusted_analysis(mydata_test,xvars)
    chronic_strata_rslt_mat[i,]=chronic_strata_analysis(nstrata,mydata_test)
    
    
  }
  
  
  rm(mydata_test)
  rm(balance_table)
  
  
  
  gc()
  
  
}

write.csv(incidence_mat,"inc_mat.csv")
write.csv(chronic_incidence_mat,"chronic_inc_mat.csv")
write.csv(strata_rslt_mat,"strata_rslt_mat.csv")
write.csv(chronic_strata_rslt_mat,"chronic_strata_mat.csv")
