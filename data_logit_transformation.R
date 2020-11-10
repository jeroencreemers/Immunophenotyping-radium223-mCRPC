library(tidyverse)
library(magrittr)
library(readxl)

source("data_pre_processing_functions.R")

#### Create a dataframe that contains all cell subset data from the Radium-223 study ####
# Read dataframes
df_a <- read_excel("radium223_data.xlsx", sheet = "panel_a_subset")

df_b <- read_excel("radium223_data.xlsx", sheet = "panel_b_subset") 

df_c <- read_excel("radium223_data.xlsx", sheet = "panel_c_subset") 

df_d <- read_excel("radium223_data.xlsx", sheet = "panel_d_subset")

# Create dataframe (wide format)
df <- list(df_a, df_b, df_c, df_d) %>% 
  reduce(left_join, by = c("study_code", "time")) %>% 
  
  # Create dataframe in long format
  pivot_longer(-c(time, study_code)) %>% 
  rename(cell_type = name, ID = study_code) %>%
  mutate(time = time - 1, 
         ID = as.factor(ID), 
         cell_type = as.factor(cell_type)) %>%
  drop_na() %>% 
  
  # Replace values of 100 with 99 (to be able to tranform data; 1 out of 2789 observations has a value of 100.0)
  mutate(value=replace(value, value==100.0, 99)) %>% 
  
  # Logit transformation
  mutate(value = logit(value))

save(df, file = "data_logit_radium223.RData")
