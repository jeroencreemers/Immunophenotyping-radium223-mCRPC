### Supplementary figure 3 ###
library(tidyverse)
library(readxl)
library(ggbeeswarm)

source("data_pre_processing_functions.R")

#### Create a dataframe that contains all cell subset dMFIs from the Radium-223 study ####
# Read dataframes
df_a_dmfi <- read_excel("radium223_data.xlsx", sheet = "panel_a_dmfi") 
df_b_dmfi <- read_excel("radium223_data.xlsx", sheet = "panel_b_dmfi") 
df_c_dmfi <- read_excel("radium223_data.xlsx", sheet = "panel_c_dmfi") 
df_d_dmfi <- read_excel("radium223_data.xlsx", sheet = "panel_d_dmfi") 
  
df <- list(df_a_dmfi, df_b_dmfi, df_c_dmfi, df_d_dmfi) %>% 
  reduce(left_join, by = c("study_code", "time")) %>%

# Create dataframe in long format
pivot_longer(-c(time, study_code)) %>% 
  rename(cell_type = name, ID = study_code) %>%
  mutate(time = time - 1, 
         ID = as.factor(ID), 
         cell_type = as.factor(cell_type)) %>%
  drop_na()

save(df, file = "checkpoint_dmfi.RData")

# Normalise dataset
df <- f_normalise(df)

save(df, file = "checkpoint_dmfi_norm.RData")

#### Plot ####
# Set aesthetics
aesthetics <- aes(x = time, y = value)

# Add universal layers to ggplot
add_universal_layers <- 
  list(
    geom_quasirandom(size = 0.5, col = "grey50"),  
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"), 
          axis.title = element_text(size = 10)),
    #scale_y_continuous(labels=function(x)round(exp(x),digits=0)), 
    scale_x_continuous(name = "Time", breaks = c(0,2,4,6), labels = c("BL", 2, 4, 6)),  
    geom_point(aes(x = time, y = group_mean_per_timepoint), col = "red", size = 0.5, shape = 16), 
    geom_smooth(aes(time, value), method = "lm", size = 0.5)
  )

# Add CD4-specific layers
add_cd4_layers <- 
  list(
    add_universal_layers, 
    theme(axis.title = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.line.x = element_blank())
  )

# Add CD8-specific layers
add_cd8_layers <- 
  list(
    add_universal_layers, 
    theme(axis.title = element_blank())
  )

# Individual plots
p_cd4_ctla4 <- df %>% filter(cell_type == "CD4_CTLA-4") %>% 
  ggplot(aesthetics) +
  ggtitle(label = "CTLA-4") +
  add_universal_layers +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.line.x = element_blank()) +
  ylab(label="dMFI")

p_cd4_pd1 <- df %>% filter(cell_type == "CD4_PD-1") %>% 
  ggplot(aesthetics) +
  ggtitle(label = "PD1") +
  add_cd4_layers

p_cd4_pdl1 <- df %>% filter(cell_type == "CD4_PD-L1") %>% 
  ggplot(aesthetics) +
  ggtitle(label = "PD-L1") +
  add_cd4_layers

p_cd4_icos <- df %>% filter(cell_type == "CD4_ICOS") %>% 
  ggplot(aesthetics) +
  ggtitle(label = "ICOS") +
  add_cd4_layers

p_cd4_tim1 <- df %>% filter(cell_type == "CD4_TIM-1") %>% 
  ggplot(aesthetics) +
  ggtitle(label = "TIM-1") +
  add_cd4_layers

p_cd4_tim3 <- df %>% filter(cell_type == "CD4_TIM-3") %>% 
  ggplot(aesthetics) +
  ggtitle(label = "TIM-3") +
  add_cd4_layers

p_cd8_ctla4 <- df %>% filter(cell_type == "CD8_CTLA-4") %>% 
  ggplot(aesthetics) +
  add_universal_layers +
  theme(axis.title.x = element_blank())+
  ylab(label = "dMFI")

p_cd8_pd1 <- df %>% filter(cell_type == "CD8_PD-1") %>% 
  ggplot(aesthetics) +
  add_cd8_layers

p_cd8_pdl1 <- df %>% filter(cell_type == "CD8_PD-L1") %>% 
  ggplot(aesthetics) +
  add_cd8_layers

p_cd8_icos <- df %>% filter(cell_type == "CD8_ICOS") %>% 
  ggplot(aesthetics) +
  add_cd8_layers

p_cd8_tim1 <- df %>% filter(cell_type == "CD8_TIM-1") %>% 
  ggplot(aesthetics) +
  add_cd8_layers

p_cd8_tim3 <- df %>% filter(cell_type == "CD8_TIM-3") %>% 
  ggplot(aesthetics) +
  add_cd8_layers


# Panel plot
cowplot::plot_grid(p_cd4_ctla4, p_cd4_icos, p_cd4_pd1, p_cd4_pdl1, p_cd4_tim1, p_cd4_tim3, 
                    p_cd8_ctla4, p_cd8_icos, p_cd8_pd1, p_cd8_pdl1, p_cd8_tim1, p_cd8_tim3,
                    ncol = 6, align = "vh")
