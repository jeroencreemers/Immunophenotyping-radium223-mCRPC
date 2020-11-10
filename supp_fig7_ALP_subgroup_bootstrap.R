### Supplementary figure 7: bootstrapping ALP responding subset###

library(tidyverse)
library(readxl)
library(boot)

source("bootstrap_functions.R")

# Load logit-transformed dataset
load("data_logit_radium223.RData")

# Filter only ALP responders
df_alp <- read_xlsx("radium223_alp_response.xlsx")

df <- full_join(df, df_alp, by = "ID") %>%
  mutate(ID = factor(ID)) %>%
  filter(alp_response == "yes") %>% 
  select(-alp_response)

###########################################################
# Data munging
df_bootstrap <- df %>% 
  pivot_wider(names_from = time, values_from = value) %>% 
  rename(t0 = "0", t2 = "2", t4 = "4", t6 = "6")  

# Bootstrapping
boot_outcome <- df_bootstrap %>% 
  group_by(cell_type) %>% 
  group_map(~ bootstrap(.x))

# Convert to readable dataframe
boot_outcome <- as.data.frame(do.call(rbind, boot_outcome)) 

# Create vector of cell types 
vec_cell_types <- df %>% 
  pull(cell_type) %>% 
  unique() %>% 
  levels()

# Bootstrap outcome table (contains fold changes!)
boot_outcome <- cbind(vec_cell_types, boot_outcome) %>% 
  rename(cell_population = vec_cell_types, fold_change = V1, conf_int = V2, lower_bound = V5, upper_bound = V6) %>% 
  mutate_if(is.numeric, round, digits=3) %>%
  select(-c(conf_int, V3, V4))

# Create dataframe with relative changes
boot_outcome_rel <- boot_outcome %>% 
  mutate(rel_change = round((fold_change - 1) *100, digits = 1), 
         lower_bound = round((lower_bound - 1) *100, digits = 1),
         upper_bound = round((upper_bound - 1) *100, digits = 1)) %>%
  select(cell_population, rel_change, lower_bound, upper_bound)

### PLOTS ###
# Add layers to ggplot
add_layers <- 
  list(
    geom_hline(yintercept = (c(0, -1, -0.5, -0.1, 0.1, 0.5, 1)+1), 
               color = "grey", 
               linetype = c("dashed", rep("dotted", 6)), 
               alpha=c(1, 0.3, 0.5, 0.8, 0.8, 0.5, 0.3)),
    geom_linerange(),
    coord_flip(ylim = c(0.5, 4)),
    geom_pointrange(size = 0.2),
    scale_x_discrete(name = "Cell population"),
    scale_y_continuous(name = "6-Months change (% of baseline)",
                       trans = "log2",
                       breaks = (c(0, -1, -0.5, -0.1, 0.1, 0.5, 1, 2, 3)+1), 
                       labels = c(0, -100, -50, -10, 10, 50, 100, 200, 300)
    )
  )

# Plot main T Cells
p_main_tcells <- boot_outcome %>% 
  filter(cell_population %in% c("CD3", "CD4", "CD8")) %>%
  ggplot(aes(x = factor(cell_population, levels = c('CD8', 'CD4', 'CD3')), y = fold_change, ymin = lower_bound, ymax = upper_bound)) +
  add_layers +
  ggtitle("Main T cells") +
  theme(axis.title = element_blank(), 
        axis.text.x = element_blank())

# Plot checkpoints on T cells
p_checkpoints <- boot_outcome %>% 
  filter(str_detect(cell_population, 'ICOS|PD-|CTLA|TIM') ) %>%
  ggplot(aes(x= factor(cell_population, levels = c("eMDSC_PD-L1", "M-MDSC_PD-L1", "Tregs_CTLA", 
                                                   "CD8_TIM-3", "CD8_TIM-1", "CD8_PD-L1", "CD8_PD-1", "CD8_ICOS", "CD8_CTLA-4", 
                                                   "CD4_TIM-3", "CD4_TIM-1", "CD4_PD-L1", "CD4_PD-1", "CD4_ICOS", "CD4_CTLA-4")), 
             y = fold_change, ymin = lower_bound, ymax = upper_bound)) +
  add_layers +
  ggtitle("Checkpoints") + 
  theme(axis.title.y = element_blank())

# Plot memory
p_memory <- boot_outcome %>% 
  filter(str_detect(cell_population, 'CD45RO')) %>%
  ggplot(aes(x= factor(cell_population, levels = c("CD8_CD45RO-CCR7-", "CD8_CD45RO+CCR7-", "CD8_CD45RO+CCR7+",
                                                   "CD4_CD45RO-CCR7-", "CD4_CD45RO+CCR7-", "CD4_CD45RO+CCR7+")), 
             y = fold_change, ymin = lower_bound, ymax = upper_bound)) +
  add_layers +
  ggtitle("Memory/effector T cells") +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank()) +
  scale_x_discrete(name = "Cell type", labels = c("RO-CCR7-","RO+CCR7-","RO+CCR7+", 
                                                  "RO-CCR7-","RO+CCR7-","RO+CCR7+"))
# Plot MDSC
p_mdsc <- boot_outcome %>% 
  filter(cell_population %in% c("Tregs", "M-MDSC", "eMDSC")) %>%
  ggplot(aes(x= factor(cell_population, levels = c("eMDSC", "M-MDSC", "Tregs")), 
             y = fold_change, ymin = lower_bound, ymax = upper_bound)) +
  add_layers +
  ggtitle("Inhibitory cells") +
  theme(axis.title.y = element_blank())

# Order plots in panel
plots <- cowplot::align_plots(p_main_tcells, p_checkpoints, align = 'h', axis = 't')
left_col <- cowplot::plot_grid(plots[[1]], p_memory, p_mdsc, ncol =1, align = 'hv', rel_heights = c(4,6,3))
cowplot::plot_grid(left_col, plots[[2]], nrow = 1)