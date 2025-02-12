#--------------------------------------------------------------------------------------------------------------------------#
# Analyses for "Blood-based biomarkers of Alzheimer’s disease and incident dementia in the community: a 16-year-long study"
# 1. Identifying the optimal cut-off for each biomarker in the training set
# 2. Test the performance of the biomarker in the test set
# R version 4.3.3
#--------------------------------------------------------------------------------------------------------------------------#



# load packages
library(tidyverse)
library(cutpointr)
library(haven)
library(rlang)



# 1. Randomly split data to training and test set ----
datasets_list <- list(
  "dataset_dem10" = dataset_dem10, # "dataset_dem10" include biomarker and 10-year all-cause dementia data
  "dataset_ad10" = dataset_ad10 # "dataset_ad10" include biomarker and 10-year AD data
)

dataset_training_list <- list()
dataset_test_list <- list()

for (dataset_name in names(datasets_list)) {
  data <- datasets_list[[dataset_name]]  
  split_ratio<-0.8
  
  set.seed(202411)
  training_set <- data %>% sample_frac(split_ratio) # Sample for training set
  
  test_set <- data %>% anti_join(training_set, by = "Lopnr") # The remaining sample is the test set
  
  dataset_training_list[[paste0(dataset_name, "_training")]] <- training_set
  dataset_test_list[[paste0(dataset_name, "_test")]] <- test_set
}




# 2. Identify the cut-off value for each biomarker that maximizes Youden's index ----
# ******************************************************************************************
# Maximize Youden's index
# (1) "x=" specifies the biomarker variable
# (2) "class=" specifies the true dementia status
# (3) "method=" specifies how to optimize the metric, e.g., maximizing the metric values
# (4) If "method = oc_manual", a single cutoff value is specified, this is used when 
#     testing the optimal cutoff in the test set.
# (5) "direction=" specify whether biomarker level above or below the cutoff is positive
# ******************************************************************************************


cp_out <- function(training_data, test_data, marker, outcome, direction) {
  # "marker" and "outcome" in the function evaluated as column names
  marker <- enquo(marker)
  outcome <- enquo(outcome)
  
  # Identify the optimal cut-off in the training set
  # The function returns the optimal cut-off and common measures - accuracy, sensitivity, specificity, and AUC
  set.seed(202411)
  cp_training_temp <- cutpointr(
    data = training_data, x = !!marker, class = !!outcome, 
    pos_class = 1, direction = direction, 
    method = maximize_boot_metric, # maxmize metric after bootstrapping
    metric = youden, # Youden's index as the metric
    boot_runs = 5000 # 5000 bootstrap replications
  )
  

  # calculate PPV etc in the training set
  # The function also returns common measures - accuracy, sensitivity, specificity, and AUC
  set.seed(202411)
  cp_training_temp_ppv <- cutpointr(
    data = training_data, x = !!marker, class = !!outcome, 
    pos_class = 1, direction = direction, 
    method = oc_manual, # manually set the cut-off as the optimal cut-off derived from the training set
    cutpoint = cp_training_temp$optimal_cutpoint,
    metric = ppv, #PPV
    boot_runs = 5000
  )
  
  
  # calculate NPV in the training set
  set.seed(202411)
  cp_training_temp_npv <- cutpointr(
    data = training_data, x = !!marker, class = !!outcome, 
    pos_class = 1, direction = direction, 
    method = oc_manual, # manually set the cut-off as the optimal cut-off derived from the training set
    cutpoint = cp_training_temp$optimal_cutpoint,
    metric = npv, #NPV
    boot_runs = 5000
  )
  
  
  # Test the performance of each biomarker in the test set using the cut-off derived from the training set
  # calculate PPV etc in the test set
  # The function also returns common measures - accuracy, sensitivity, specificity, and AUC
  set.seed(202411)
  cp_test_temp_ppv <- cutpointr(
    data = test_data, x = !!marker, class = !!outcome, 
    pos_class = 1, direction = direction, 
    method = oc_manual, # manually set the cut-off as the optimal cut-off derived from the training set
    cutpoint = cp_training_temp$optimal_cutpoint,
    metric = ppv, #PPV
    boot_runs = 5000
  )
  
  
  # calculate NPV in the test set
  set.seed(202411)
  cp_test_temp_npv <- cutpointr(
    data = test_data, x = !!marker, class = !!outcome, 
    pos_class = 1, direction = direction, 
    method = oc_manual, # manually set the cut-off as the optimal cut-off derived from the training set
    cutpoint = cp_training_temp$optimal_cutpoint,
    metric = npv, #NPV
    boot_runs = 5000
  )
  
  
  # estimates and 95% CI in the training set
  # optimal cut-off and 95% CI
  cp<-as.numeric(cp_training_temp$optimal_cutpoint)
  cp_ci<-boot_ci(cp_training_temp,optimal_cutpoint)
  cp_ci<-paste0(format(round(cp,digits = 3),nsmall = 3)," (",
                format(round(as.numeric(cp_ci[1,2]),digits = 3),nsmall = 3),"-",
                format(round(as.numeric(cp_ci[2,2]),digits = 3),nsmall = 3),")")
  
  # PPV and 95% CI
  train_ppv<-as.numeric(cp_training_temp_ppv$ppv)
  train_ppv_ci<-boot_ci(cp_training_temp_ppv,ppv)
  train_ppv_ci<-paste0(format(round(train_ppv*100,digits = 1),nsmall = 1)," (",
                       format(round(train_ppv_ci[1,2]*100,digits = 1),nsmall = 1),"-",
                       format(round(train_ppv_ci[2,2]*100,digits = 1),nsmall = 1),")")
  
  # NPV and 95% CI
  train_npv<-as.numeric(cp_training_temp_npv$npv)
  train_npv_ci<-boot_ci(cp_training_temp_npv,npv)
  train_npv_ci<-paste0(format(round(train_npv*100,digits = 1),nsmall = 1)," (",
                       format(round(train_npv_ci[1,2]*100,digits = 1),nsmall = 1),"-",
                       format(round(train_npv_ci[2,2]*100,digits = 1),nsmall = 1),")")
  
  # Accuracy and 95% CI
  train_acc<-as.numeric(cp_training_temp_ppv$acc)
  train_acc_ci<-boot_ci(cp_training_temp_ppv,acc)
  train_acc_ci<-paste0(format(round(train_acc*100,digits = 1),nsmall = 1)," (",
                       format(round(train_acc_ci[1,2]*100,digits = 1),nsmall = 1),"-",
                       format(round(train_acc_ci[2,2]*100,digits = 1),nsmall = 1),")")
  
  # Sensitivity and 95% CI
  train_sens<-as.numeric(cp_training_temp_ppv$sensitivity)
  train_sens_ci<-boot_ci(cp_training_temp_ppv,sensitivity)
  train_sens_ci<-paste0(format(round(train_sens*100,digits = 1),nsmall = 1)," (",
                       format(round(train_sens_ci[1,2]*100,digits = 1),nsmall = 1),"-",
                       format(round(train_sens_ci[2,2]*100,digits = 1),nsmall = 1),")")
  
  # Specificity and 95% CI
  train_spe<-as.numeric(cp_training_temp_ppv$specificity)
  train_spe_ci<-boot_ci(cp_training_temp_ppv,specificity)
  train_spe_ci<-paste0(format(round(train_spe*100,digits = 1),nsmall = 1)," (",
                       format(round(train_spe_ci[1,2]*100,digits = 1),nsmall = 1),"-",
                       format(round(train_spe_ci[2,2]*100,digits = 1),nsmall = 1),")")
  
  # AUC and 95% CI
  train_auc<-as.numeric(cp_training_temp_ppv$AUC)
  train_auc_ci<-boot_ci(cp_training_temp_ppv,AUC)
  train_auc_ci<-paste0(format(round(train_auc,digits = 3),nsmall = 3)," (",
                       format(round(as.numeric(train_auc_ci[1,2]),digits = 3),nsmall = 3),"-",
                       format(round(as.numeric(train_auc_ci[2,2]),digits = 3),nsmall = 3),")")
  
  # estimates and 95% CI in the test set
  # PPV and 95% CI
  test_ppv<-as.numeric(cp_test_temp_ppv$ppv)
  test_ppv_ci<-boot_ci(cp_test_temp_ppv,ppv)
  test_ppv_ci<-paste0(format(round(test_ppv*100,digits = 1),nsmall = 1)," (",
                       format(round(test_ppv_ci[1,2]*100,digits = 1),nsmall = 1),"-",
                       format(round(test_ppv_ci[2,2]*100,digits = 1),nsmall = 1),")")
  
  # NPV and 95% CI
  test_npv<-as.numeric(cp_test_temp_npv$npv)
  test_npv_ci<-boot_ci(cp_test_temp_npv,npv)
  test_npv_ci<-paste0(format(round(test_npv*100,digits = 1),nsmall = 1)," (",
                      format(round(test_npv_ci[1,2]*100,digits = 1),nsmall = 1),"-",
                      format(round(test_npv_ci[2,2]*100,digits = 1),nsmall = 1),")")
  
  # Accuracy and 95% CI
  test_acc<-as.numeric(cp_test_temp_ppv$acc)
  test_acc_ci<-boot_ci(cp_test_temp_ppv,acc)
  test_acc_ci<-paste0(format(round(test_acc*100,digits = 1),nsmall = 1)," (",
                       format(round(test_acc_ci[1,2]*100,digits = 1),nsmall = 1),"-",
                       format(round(test_acc_ci[2,2]*100,digits = 1),nsmall = 1),")")
  
  # Sensitivity and 95% CI
  test_sens<-as.numeric(cp_test_temp_ppv$sensitivity)
  test_sens_ci<-boot_ci(cp_test_temp_ppv,sensitivity)
  test_sens_ci<-paste0(format(round(test_sens*100,digits = 1),nsmall = 1)," (",
                        format(round(test_sens_ci[1,2]*100,digits = 1),nsmall = 1),"-",
                        format(round(test_sens_ci[2,2]*100,digits = 1),nsmall = 1),")")
  
  # Specificity and 95% CI
  test_spe<-as.numeric(cp_test_temp_ppv$specificity)
  test_spe_ci<-boot_ci(cp_test_temp_ppv,specificity)
  test_spe_ci<-paste0(format(round(test_spe*100,digits = 1),nsmall = 1)," (",
                       format(round(test_spe_ci[1,2]*100,digits = 1),nsmall = 1),"-",
                       format(round(test_spe_ci[2,2]*100,digits = 1),nsmall = 1),")")
  
  # AUC and 95% CI
  test_auc<-as.numeric(cp_test_temp_ppv$AUC)
  test_auc_ci<-boot_ci(cp_test_temp_ppv,AUC)
  test_auc_ci<-paste0(format(round(test_auc,digits = 3),nsmall = 3)," (",
                       format(round(as.numeric(test_auc_ci[1,2]),digits = 3),nsmall = 3),"-",
                       format(round(as.numeric(test_auc_ci[2,2]),digits = 3),nsmall = 3),")")
  
  return(data.frame(cp = cp,
    cp_ci = cp_ci, train_ppv_ci = train_ppv_ci, train_npv_ci = train_npv_ci,
    train_acc_ci = train_acc_ci,train_sens_ci = train_sens_ci, 
    train_spe_ci = train_spe_ci,train_auc_ci = train_auc_ci, 
    test_ppv_ci = test_ppv_ci, test_npv_ci = test_npv_ci,
    test_acc_ci = test_acc_ci, test_sens_ci = test_sens_ci,
    test_spe_ci = test_spe_ci, test_auc_ci = test_auc_ci
  ))
}

empty_frame <- data.frame(data = character(0), cp=numeric(0),cp_ci = character(0), 
                          train_ppv_ci = character(0),train_npv_ci = character(0),
                          train_acc_ci = character(0),train_sens_ci = character(0),
                          train_spe_ci = character(0), train_auc_ci = character(0),
                          test_ppv_ci = character(0),test_npv_ci = character(0),
                          test_acc_ci = character(0),test_sens_ci = character(0), 
                          test_spe_ci = character(0),test_auc_ci = character(0)) 

# 2.1. ab ratio ----
cp_abratio<-empty_frame

# ab ratio for all-cause dementia
for (i in grep("dem10",names(dataset_training_list))){
  cp_abratio<-rbind(cp_abratio,
                    cbind(data=gsub("_training", "", names(dataset_training_list)[i]),
                          cp_out(training_data = dataset_training_list[[i]],
                           test_data = dataset_test_list[[i]],
                           marker = abetaratio, outcome = dem10, 
                           direction = "<=")) %>% 
                      as.data.frame() %>% 
                      set_names(names(empty_frame)))
}


# ab ratio for AD
for (i in grep("ad10",names(dataset_training_list))){
  cp_abratio<-rbind(cp_abratio,
                    cbind(data=gsub("_training", "", names(dataset_training_list)[i]),
                          cp_out(training_data = dataset_training_list[[i]],
                           test_data = dataset_test_list[[i]],
                           marker = abetaratio, outcome = ad10, 
                           direction = "<=")) %>% 
                      as.data.frame() %>% 
                      set_names(names(empty_frame)))
}


# 2.2. ptau181 ----
cp_ptau181<-empty_frame

# ptau181 for all-cause dementia
for (i in grep("dem10",names(dataset_training_list))){
cp_ptau181<-rbind(cp_ptau181,
                 cbind(data=gsub("_training", "", names(dataset_training_list)[i]),
                       cp_out(training_data = dataset_training_list[[i]],
                        test_data = dataset_test_list[[i]],
                        marker = ptau181, outcome = dem10, 
                        direction = ">=")) %>%
                   as.data.frame() %>%
                   set_names(names(empty_frame)))
}


# ptau181 for AD
for (i in grep("ad10",names(dataset_training_list))){
  cp_ptau181<-rbind(cp_ptau181,
                 cbind(data=gsub("_training", "", names(dataset_training_list)[i]),
                  cp_out(training_data = dataset_training_list[[i]],
                        test_data = dataset_test_list[[i]],
                        marker = ptau181, outcome = ad10, 
                        direction = ">=")) %>%
                   as.data.frame() %>%
                   set_names(names(empty_frame)))
}



# 2.3. ptau217 ----
cp_ptau217<-empty_frame

# ptau217 for all-cause dementia
for (i in grep("dem10",names(dataset_training_list))){
  cp_ptau217<-rbind(cp_ptau217,
                    cbind(data=gsub("_training", "", names(dataset_training_list)[i]),
                     cp_out(training_data = dataset_training_list[[i]],
                           test_data = dataset_test_list[[i]],
                           marker = ptau217, outcome = dem10, 
                           direction = ">=")) %>%
                      as.data.frame() %>%
                      set_names(names(empty_frame)))
}


# ptau217 for AD
for (i in grep("ad10",names(dataset_training_list))){
  cp_ptau217<-rbind(cp_ptau217,
                    cbind(data=gsub("_training", "", names(dataset_training_list)[i]),
                          cp_out(training_data = dataset_training_list[[i]],
                           test_data = dataset_test_list[[i]],
                           marker = ptau217, outcome = ad10, 
                           direction = ">=")) %>%
                      as.data.frame() %>%
                      set_names(names(empty_frame)))
}


# 2.4. ttau ----
cp_ttau<-empty_frame

# ttau for all-cause dementia
for (i in grep("dem10",names(dataset_training_list))){
  cp_ttau<-rbind(cp_ttau,
                 cbind(data=gsub("_training", "", names(dataset_training_list)[i]),
                  cp_out(training_data = dataset_training_list[[i]],
                        test_data = dataset_test_list[[i]],
                        marker = ttau, outcome = dem10, 
                        direction = ">=")) %>%
                   as.data.frame() %>%
                   set_names(names(empty_frame)))
}



# ttau for AD
for (i in grep("ad10",names(dataset_training_list))){
  cp_ttau<-rbind(cp_ttau,
                 cbind(data=gsub("_training", "", names(dataset_training_list)[i]),
                  cp_out(training_data = dataset_training_list[[i]],
                        test_data = dataset_test_list[[i]],
                        marker = ttau, outcome = ad10, 
                        direction = ">=")) %>%
                   as.data.frame() %>%
                   set_names(names(empty_frame)))
}




# 2.5. nfl ----
cp_nfl<-empty_frame

# nfl for all-cause dementia
for (i in grep("dem10",names(dataset_training_list))){
  cp_nfl<-rbind(cp_nfl,
                cbind(data=gsub("_training", "", names(dataset_training_list)[i]),
                      cp_out(training_data = dataset_training_list[[i]],
                       test_data = dataset_test_list[[i]],
                       marker = nfl, outcome = dem10, 
                       direction = ">=")) %>% 
                  as.data.frame() %>% 
                  set_names(names(empty_frame)))
}



# nfl for AD
for (i in grep("ad10",names(dataset_training_list))){
  cp_nfl<-rbind(cp_nfl,
                 cbind(data=gsub("_training", "", names(dataset_training_list)[i]),
                       cp_out(training_data = dataset_training_list[[i]],
                        test_data = dataset_test_list[[i]],
                        marker = nfl, outcome = ad10, 
                        direction = ">=")) %>% 
                  as.data.frame() %>% 
                  set_names(names(empty_frame)))
}




# 2.6. gfap ----
cp_gfap<-empty_frame

# gfap for all-cause dementia
for (i in grep("dem10",names(dataset_training_list))){
  cp_gfap<-rbind(cp_gfap,
                cbind(data=gsub("_training", "", names(dataset_training_list)[i]),
                      cp_out(training_data = dataset_training_list[[i]],
                       test_data = dataset_test_list[[i]],
                       marker = gfap, outcome = dem10, 
                       direction = ">=")) %>% 
                  as.data.frame() %>% 
                  set_names(names(empty_frame)))
}



# gfap for AD
for (i in grep("ad10",names(dataset_training_list))){
  cp_gfap<-rbind(cp_gfap,
                cbind(data=gsub("_training", "", names(dataset_training_list)[i]),
                      cp_out(training_data = dataset_training_list[[i]],
                       test_data = dataset_test_list[[i]],
                       marker = gfap, outcome = ad10, 
                       direction = ">=")) %>% 
                  as.data.frame() %>% 
                  set_names(names(empty_frame)))
}



# 3. Combine all tables together ----
table_sum<-cp_abratio %>% 
  mutate(marker="abetaratio") %>% 
  
  rbind(cp_ptau181 %>% 
          mutate(marker="ptau181")) %>% 
  
  rbind(cp_ptau217 %>% 
          mutate(marker="ptau217")) %>% 
  
  rbind(cp_ttau %>% 
          mutate(marker="ttau")) %>% 
  
  rbind(cp_nfl %>% 
          mutate(marker="nfl")) %>% 
  
  rbind(cp_gfap %>% 
          mutate(marker="gfap"))


table_sum<-table_sum %>% 
  mutate(outcome=ifelse(grepl("dem",data),"Dementia","AD")) %>% 
  select(marker,outcome,everything()) %>% 
  select(-data)


# Tidy up the final summary table and edit column names
table_sum_tidy<-table_sum %>%
  select(-cp)

colnames(table_sum_tidy)<-c("Biomarker","Outcome","Cut-off value",
                       "PPV in the training set","NPV in the training set",
                       "Accuracy in the training set","Sensitivity in the training set",
                       "Specificity in the training set","AUC in the training set",
                       "PPV in the test set","NPV in the test set",
                       "Accuracy in the test set","Sensitivity in the test set",
                       "Specificity in the test set","AUC in the test set")


  
