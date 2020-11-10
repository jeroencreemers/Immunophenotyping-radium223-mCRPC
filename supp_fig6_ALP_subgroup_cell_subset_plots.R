### Supplementary figure 6: Radium-223 subsets in ALP responderds ###

## Note: all plots show normalized values!!

library(ggplot2)
library(dplyr)
library(magrittr)
library(ggbeeswarm)
library(scales)
library(readxl)

source("data_pre_processing_functions.R")

# Load the complete dataset (logit transformed + normalised)
load("data_logit_norm_radium223.RData")

# Filter only ALP responders
df_alp <- read_xlsx("radium223_alp_response.xlsx")

df <- full_join(df, df_alp, by = "ID") %>%
  mutate(ID = factor(ID)) %>%
  filter(alp_response == "yes") %>% 
  select(-alp_response)

# Add layers to ggplot
add_layers <- 
  list(
    geom_quasirandom(size = 0.5, col = "grey50"),  
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"), 
          axis.title = element_text(size = 10)), 
    scale_y_continuous(breaks = logit(c(1,2,5,10, 20, 30, 40, 50, 60, 70, 80, 90)), labels=function(x)round(100*plogis(x),digits=1)), 
    scale_x_continuous(name = "Time", breaks = c(0,2,4,6), labels = c("BL", 2, 4, 6)), 
    geom_point(aes(x = time, y = group_mean_per_timepoint), col = "red", size = 0.5, shape = 16),
    geom_smooth(aes(time, value), method = "lm", size = 0.5)
  )

# Set aesthetic
aesthetics <- aes(x = time, y = value)
#############################################################
#### Panel B ####
#### CD3+ T cells###
p_cd3_tcells <- df %>% filter(cell_type == "CD3") %>% 
  ggplot(aesthetics) +
  ggtitle(label = "CD3") +
  add_layers +
  theme(axis.title.x = element_blank()) 

#### CD4+ T cells###
p_cd4_tcells <- df %>% filter(cell_type == "CD4") %>% 
  ggplot(aesthetics) +
  ggtitle(label = "CD4") +
  add_layers +
  theme(axis.title = element_blank())

#### CD8+ T cells###
p_cd8_tcells <- df %>% filter(cell_type == "CD8") %>% 
  ggplot(aesthetics) +
  ggtitle(label = "CD8") +
  add_layers +
  theme(axis.title = element_blank()) 

#############################################################
#### Panel C ####
remove_title_and_xaxis <- theme(axis.title = element_blank(), 
                                axis.text.x = element_blank(), 
                                axis.line.x = element_blank(), 
                                axis.ticks.x = element_blank())

p_cd4_ctla4_tcells <- df %>% filter(cell_type == "CD4_CTLA-4") %>% 
  ggplot(aesthetics) +
  ggtitle("CTLA-4") + 
  add_layers +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.line.x = element_blank(), 
        axis.ticks.x = element_blank()) +
  ylab(label = expression( CD4^{textstyle("+")}~cells~(textstyle("%")) ))

#### CD4 ICOS ####
p_cd4_icos_tcells <- df %>% filter(cell_type == "CD4_ICOS") %>% 
  ggplot(aesthetics) +
  ggtitle("ICOS") + 
  add_layers +
  remove_title_and_xaxis

#### CD4 PD-1 ####
p_cd4_pd1_tcells <- df %>% filter(cell_type == "CD4_PD-1") %>% 
  ggplot(aesthetics) +
  ggtitle("PD-1") + 
  add_layers +
  remove_title_and_xaxis

#### CD4 PD-L1 ####
p_cd4_pdl1_tcells <- df %>% filter(cell_type == "CD4_PD-L1") %>% 
  ggplot(aesthetics) +
  ggtitle("PD-L1") + 
  add_layers +
  remove_title_and_xaxis

#### CD4 TIM-1 ####
p_cd4_tim1_tcells <- df %>% filter(cell_type == "CD4_TIM-1") %>% 
  ggplot(aesthetics) +
  ggtitle("TIM-1") + 
  add_layers +
  remove_title_and_xaxis

#### CD4 TIM-3 ####
p_cd4_tim3_tcells <- df %>% filter(cell_type == "CD4_TIM-3") %>% 
  ggplot(aesthetics) +
  ggtitle("TIM-3") + 
  add_layers +
  remove_title_and_xaxis

################################################################################
#### Part II ####

#### CD8 CTLA-4 ####
p_cd8_ctla4_tcells <- df %>% filter(cell_type == "CD8_CTLA-4") %>% 
  ggplot(aesthetics) +
  add_layers +
  theme(axis.title.x = element_blank()) +
  ylab(label = expression( CD8^{textstyle("+")}~cells~(textstyle("%")) ))

#### CD8 ICOS ####
p_cd8_icos_tcells <- df %>% filter(cell_type == "CD8_ICOS") %>% 
  ggplot(aesthetics) +
  add_layers +
  theme(axis.title = element_blank())

#### CD8 PD-1 ####
p_cd8_pd1_tcells <- df %>% filter(cell_type == "CD8_PD-1") %>% 
  ggplot(aesthetics) +
  add_layers +
  theme(axis.title = element_blank())

#### CD8 PD-L1 ####
p_cd8_pdl1_tcells <- df %>% filter(cell_type == "CD8_PD-L1") %>% 
  ggplot(aesthetics) +
  add_layers +
  theme(axis.title = element_blank())

#### CD8 TIM-1 ####
p_cd8_tim1_tcells <- df %>% filter(cell_type == "CD8_TIM-1") %>% 
  ggplot(aesthetics) +
  add_layers +
  theme(axis.title = element_blank())

#### CD8 TIM-3 ####
p_cd8_tim3_tcells <- df %>% filter(cell_type == "CD8_TIM-3") %>% 
  ggplot(aesthetics) +
  add_layers +
  theme(axis.title = element_blank())

################################################################################
#### Panel D #### 

#### CD4 CD45RO- CCR7- ####
p_cd4_cd45roneg_ccr7neg <- df %>% filter(cell_type == "CD4_CD45RO-CCR7-") %>% 
  ggplot(aesthetics) +
  add_layers +
  ggtitle(expression( bold (CD45RO^{textstyle("-")}*CCR7^{textstyle("-")} ) ) ) +
  ylab(label = "Cells (%)") +
  theme(axis.title.x = element_blank())

#### CD4 CD45RO+ CCR7- ####
p_cd4_cd45ropos_ccr7neg <- df %>% filter(cell_type == "CD4_CD45RO+CCR7-") %>% 
  ggplot(aesthetics) +
  add_layers +
  ggtitle(expression( bold (CD45RO^{textstyle("+")}*CCR7^{textstyle("-")} ) ) ) +
  theme(axis.title = element_blank())

#### CD4 CD45RO+ CCR7+ ####
p_cd4_cd45ropos_ccr7pos <- df %>% filter(cell_type == "CD4_CD45RO+CCR7+") %>% 
  ggplot(aesthetics) +
  add_layers +
  ggtitle(expression( bold (CD45RO^{textstyle("+")}*CCR7^{textstyle("+")} ) ) ) +
  theme(axis.title = element_blank())

#### CD8 CD45RO- CCR7- ####
p_cd8_cd45roneg_ccr7neg <- df %>% filter(cell_type == "CD8_CD45RO-CCR7-") %>% 
  ggplot(aesthetics) +
  add_layers +
  ggtitle(expression( bold (CD45RO^{textstyle("-")}*CCR7^{textstyle("-")} ) ) ) +
  theme(axis.title = element_blank())

#### CD8 CD45RO+ CCR7- ####
p_cd8_cd45ropos_ccr7neg <- df %>% filter(cell_type == "CD8_CD45RO+CCR7-") %>% 
  ggplot(aesthetics) +
  add_layers +
  ggtitle(expression( bold (CD45RO^{textstyle("+")}*CCR7^{textstyle("-")} ) ) ) +
  theme(axis.title = element_blank())

#### CD8 CD45RO+ CCR7+ ####
p_cd8_cd45ropos_ccr7pos <- df %>% filter(cell_type == "CD8_CD45RO+CCR7+") %>% 
  ggplot(aesthetics) +
  add_layers +
  ggtitle(expression( bold (CD45RO^{textstyle("+")}*CCR7^{textstyle("+")} ) ) ) +
  theme(axis.title = element_blank())

#### Tregs cells###
p_tregs_tcells <- df %>% filter(cell_type == "Tregs") %>% 
  ggplot(aesthetics) +
  ggtitle(label = "Tregs") +
  add_layers +
  theme(axis.title.x = element_blank()) +
  ylab(label = "Cells (%)")

#### Tregs CTLA-4 cells###
p_tregs_ctla_tcells <- df %>% filter(cell_type == "Tregs_CTLA") %>% 
  ggplot(aesthetics) +
  ggtitle(label = "Tregs CTLA-4") +
  add_layers +
  theme(axis.title = element_blank())

#### MDSC ####
p_mdsc <- df %>% filter(cell_type == "M-MDSC") %>% 
  ggplot(aesthetics) +
  ggtitle(label = "M-MDSC") +
  add_layers +
  theme(axis.title = element_blank())

#### MDSC PD-L1 ####
p_mdsc_pdl1 <- df %>% filter(cell_type == "M-MDSC_PD-L1") %>% 
  ggplot(aesthetics) +
  ggtitle(label = "M-MDSC PD-L1") +
  add_layers +
  theme(axis.title = element_blank()) +
  ylab(label = expression( PD-L1^{textstyle("+")}~cells~(textstyle("%")) ) )

#### eMDSC ####
p_emdsc <- df %>% filter(cell_type == "eMDSC") %>% 
  ggplot(aesthetics) +
  ggtitle(label = "e-MDSC") +
  add_layers +
  theme(axis.title = element_blank())

#### eMDSC PD-L1 ####
p_emdsc_pdl1 <- df %>% filter(cell_type == "eMDSC_PD-L1") %>% 
  ggplot(aesthetics) +
  ggtitle(label = "eMDSC PD-L1") +
  add_layers +
  theme(axis.title = element_blank())

# Panel plot
cowplot::plot_grid(p_cd3_tcells, p_cd4_tcells, p_cd8_tcells, NULL, NULL, NULL,  
                    p_cd4_ctla4_tcells, p_cd4_icos_tcells, p_cd4_pd1_tcells, p_cd4_pdl1_tcells, p_cd4_tim1_tcells, p_cd4_tim3_tcells,
                    p_cd8_ctla4_tcells, p_cd8_icos_tcells, p_cd8_pd1_tcells, p_cd8_pdl1_tcells, p_cd8_tim1_tcells, p_cd8_tim3_tcells,
                    p_cd4_cd45roneg_ccr7neg , p_cd4_cd45ropos_ccr7neg, p_cd4_cd45ropos_ccr7pos, p_cd8_cd45roneg_ccr7neg, p_cd8_cd45ropos_ccr7neg, p_cd8_cd45ropos_ccr7pos,
                    p_tregs_tcells, p_tregs_ctla_tcells, p_mdsc, p_mdsc_pdl1, p_emdsc, p_emdsc_pdl1,
                    ncol = 6, align = "vh")