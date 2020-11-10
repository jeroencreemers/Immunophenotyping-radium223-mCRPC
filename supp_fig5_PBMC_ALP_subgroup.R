### Supplementary figure 5: PBMC availability - ALP response subgroup ###
library(tidyverse)
library(readxl)

# Load data
load("data_logit_norm_radium223.RData")

# Data munging
df_alp <- read_xlsx("radium223_alp_response.xlsx")

df <- full_join(df, df_alp, by = "ID") %>%
  mutate(ID = factor(ID), 
         time = factor(time)) %>%
  select(-group_mean_per_timepoint)

save(df, file = "data_alp_response_subgroup.RData")

# Plot
ggplot(data = df, aes(x = time, y = ID, z = cell_type)) + 
  geom_tile(fill = ifelse(df$alp_response=="no", "white", "grey"), 
            width = 0.8, height = 0.8) +
  scale_x_discrete(name = "Injection number", 
                   labels = c("BL", 2, 4, 6)) +
  scale_y_discrete(name = "Patient") + 
  theme_classic()