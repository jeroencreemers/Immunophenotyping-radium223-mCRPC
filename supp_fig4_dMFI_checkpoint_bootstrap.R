### Supplementary figure 4: Bootstrapping (delta)MFI ###

library(tidyverse)
library(boot)

load("checkpoint_dmfi.RData")
source("bootstrap_dmfi_functions.R")

# Data munging of
min_value <- abs(min(df$value))+1 # absolute min value +1

df_bootstrap <- df %>% 
  mutate(value = value + min_value) %>% # translation step to obtain only positive values, necessary for plotting on log2-transformed axis
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
    coord_flip(ylim = c(0.5, 2)),
    geom_pointrange(size = 0.2),
    scale_x_discrete(name = "Cell population"),
    scale_y_continuous(name = "6-Months change (% of baseline)",
                       trans = "log2",
                       breaks = (c(0, -1, -0.5, -0.1, 0.1, 0.5, 1, 2)+1), 
                       labels = c(0, -100, -50, -10, 10, 50, 100, 200)
    )
  )

# Plot dMFI checkpoints on T cells
p_checkpoints <- boot_outcome %>% 
  filter(str_detect(cell_population, 'ICOS|PD-|CTLA|TIM') ) %>%
  ggplot(aes(x= factor(cell_population, levels = c("eMDSC_PD-L1", "M-MDSC_PD-L1", "Tregs_CTLA-4", 
                                                   "CD8_TIM-3", "CD8_TIM-1", "CD8_PD-L1", "CD8_PD-1", "CD8_ICOS", "CD8_CTLA-4", 
                                                   "CD4_TIM-3", "CD4_TIM-1", "CD4_PD-L1", "CD4_PD-1", "CD4_ICOS", "CD4_CTLA-4")), 
             y = fold_change, ymin = lower_bound, ymax = upper_bound)) +
  add_layers +
  ggtitle("Checkpoints") + 
  theme(axis.title.y = element_blank())