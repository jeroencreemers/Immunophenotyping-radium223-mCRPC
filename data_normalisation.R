library(tidyverse)
library(magrittr)

load(file = "data_logit_radium223.RData")
source(file = "data_pre_processing_functions.R")

#### Normalisation of dataset of Radium-223 study ####

df <- f_normalise(df)

save(df, file = "data_logit_norm_radium223.RData")







