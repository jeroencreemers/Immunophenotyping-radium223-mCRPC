#### Figure 1: Absolute lymphocyte & monocyte counts ####
library(tidyverse)
library(readxl)
library(ggbeeswarm)
library(cowplot)
source("data_pre_processing_functions.R")

################# Data munging ##################
# Read dataframe from excel and tranform to long format
df_lymphocyte <- read_excel("abs_lymph_monocyte_count.xlsx", sheet = "lymphocyte_count") %>% 
  pivot_longer(-ID) %>%
  mutate(cell_type = "lymphocyte")

df_monocyte <- read_excel("abs_lymph_monocyte_count.xlsx", sheet = "monocyte_count") %>% 
  pivot_longer(-ID) %>%
  mutate(cell_type = "monocyte")

df <- full_join(df_lymphocyte, df_monocyte) %>% 
  rename(time = name) 

# Data normalisation
df <- f_normalise(df)
df <- df %>% 
  mutate(time = factor(time, ordered = T, c("BL", "1", "2", "3", "4", "5", "6"))) 

#################### Plot #######################
# Note: plots show the normalized absolute cell counts.

# Add layers to ggplot
add_layers <- 
  list(
    geom_quasirandom(size = 0.5, col = "grey50"),
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"), 
          axis.title = element_text(size = 10)), 
    scale_y_continuous(name = expression(Cells~(x10^{textstyle("9")}/L) )  ), 
    scale_x_discrete(name = "Injection cycle"), 
    geom_smooth(aes(as.numeric(time), value), method = "lm"), 
    geom_point(aes(x = time, y = group_mean_per_timepoint), col = "red", size = 0.5, shape = 16)
  )

# Plot lymphocytes
abs_lymph_plot <- df %>% filter(cell_type == "lymphocyte") %>% 
  ggplot(aes(x = time, y = value)) +
  add_layers + 
  ggtitle(label = "Lymphocytes")

#### Summary of regression model (lymphocytes)
summary(lm(value~as.numeric(time), subset(df, cell_type == "lymphocyte")))

# Plot monocytes
abs_monocyte_plot <- df %>% filter(cell_type == "monocyte") %>% 
  ggplot(aes(x = time, y = value)) +
  add_layers + 
  ggtitle(label = "Monocytes") +
  theme(axis.title.y = element_blank())

# Summary of regression model (monocytes)
summary(lm(value~as.numeric(time), subset(df, cell_type == "monocyte")))

plot_grid(abs_lymph_plot, abs_monocyte_plot, nrow = 1, align = 'vh', axis = 'lb')
