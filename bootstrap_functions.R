################### Bootstrap functions ######################

source("data_pre_processing_functions.R")

# Function of the statistic
f_statistic <- function(df, b=1:nrow(df)){
  
  # allow boot() to select sample
  df <- df[b,] 
  
  # Calculate fold change from baseline by fitting
  # a linear model to the data. Subsequently, change from 
  # baseline (%) is calculated.
  
  df <- df %>% 
    pivot_longer(cols = c(t0, t2, t4, t6), names_to = "time") 
  df$time <- as.numeric(c("t0" = "0", "t2" = "2", "t4" = "4", "t6"="6")[df$time])
  
  # Normalise data 
  overall_mean <- mean(df$value, na.rm = TRUE)
  
  df <- df %>% 
    group_by(ID) %>% 
    mutate(patient_mean = mean(value, na.rm = TRUE)) %>% 
    ungroup() %>% 
    mutate(value = value - patient_mean + overall_mean ) %>% 
    select(-patient_mean) 
  
  # Fit model to logit transformed values
  model <- lm(formula = value~time, data = df)
  
  prediction_baseline_logit <- predict(model, data.frame(time = c(0)))
  prediction_6mo_logit <- predict(model, data.frame(time = c(6)))
  
  # Transform logit to normal (with logistic function)
  prediction_baseline <- logistic(prediction_baseline_logit)
  prediction_6mo <- logistic(prediction_6mo_logit)
  
  fold_change <- prediction_6mo / prediction_baseline
  
  return( unname(fold_change) ) 
}

##############################################################
# Function to perform the bootstrap
bootstrap <- function(df){
  
  # Create bootstrap object
  boot.obj <- boot(data = df, statistic = f_statistic, R = 2000)
  
  return(c(boot.obj$t0, boot.ci(boot.obj)$percent))
}