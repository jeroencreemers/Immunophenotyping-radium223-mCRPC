library(tidyverse)
library(magrittr)

#### Functions ###

# Logit transformation function
logit <- function(x){
  x <- x/100
  log(x/(1-x))
}


# Logistic transformation (transforms logit data back to normal proportions)
logistic <- function(x){
  ( 1/(1+exp(-x)) ) * 100
}


# Function to normalise Radium-223 data
f_normalise <- function(x){
  df %<>% 
    # Calculate overall group mean and add to dataframe
    group_by(cell_type) %>% 
    mutate(overall_mean = mean(value, na.rm = TRUE)) %>%
    
    # Calculate patient mean per cell type
    group_by(ID, cell_type) %>%
    mutate(patient_mean = mean(value, na.rm = TRUE)) %>%
    ungroup() %>%
    
    # Normalize values by patient
    # (first: subtract patient mean per marker from the values, then: add mean of marker)
    mutate(value = value - patient_mean + overall_mean ) %>% 
    select(-c(patient_mean, overall_mean)) %>%
    
    # Calculate group mean per timepoint
    group_by(cell_type, time) %>%
    mutate(group_mean_per_timepoint = mean(value, na.rm = TRUE)) %>% 
    ungroup 
}