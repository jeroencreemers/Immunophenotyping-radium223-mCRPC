## Supplementary figure 2: PBMC Availability - Radium-223 entire population ##
library(dplyr)
library(ggplot2)

load("data_logit_radium223.RData")

df <- df %>% select(ID ,time, cell_type, value) %>% 
  mutate(time = factor(time))

ggplot(data = df, aes(x = time, y = ID, z = cell_type)) + 
  geom_tile(fill = "grey", width = 0.8, height = 0.8) +
  scale_x_discrete(name = "Injection number",
                   labels = c("BL", 2, 4, 6)) +
  scale_y_discrete(name = "Patient") +
  theme_classic()
