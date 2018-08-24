### 1B - repellent screen:

plan <- drake::drake_plan(
  library(tidyverse),
  library(patchwork),

#filepath = here::here("data/Figure1B_avoidance.csv"),

screen_data = read.csv(drake::file_in(here::here("data/Figure1B_avoidance.csv"))) %>%
  filter(assay %in% c('oct', 'nonanone')) %>%
  format_CIdata() %>% mutate(strain = fct_relevel(strain, 'OP50')),

meanOP50 = screen_data %>%
  dplyr::filter(strain == 'OP50') %>%
  group_by(assay) %>%
  summarise(meanLogit = mean(logit.p)),

# mean_nonanone = dplyr::filter(meanOP50, assay == 'nonanone') %>% dplyr::select(meanLogit) %>% as.numeric(),
# mean_octanol = dplyr::filter(meanOP50, assay == 'oct') %>% dplyr::select(meanLogit) %>% as.numeric(),

nonanone_data = screen_data %>%
  dplyr::filter(assay == 'nonanone') %>%
  mutate(rel.Logit = logit.p - mean(dplyr::filter(., strain == 'OP50')$logit.p)),

octanol_data = screen_data %>%
  dplyr::filter(assay == 'oct') %>%
  mutate(rel.Logit = logit.p - mean(dplyr::filter(., strain == 'OP50')$logit.p)),

plotColors = c('#827E7E','#C93030', '#2F8A34','#484CC7','#B8B1B1', '#E6DDDD', 'gray94'),

plot1 = nonanone_data %>% dplyr::filter(genotype == 'N2') %>%
  plot_plasticityIndex(xvar = strain,
                       plotColors = c('#827E7E','#C93030', '#2F8A34','#484CC7','#B8B1B1', '#E6DDDD', 'gray94')) +
  guides(fill = FALSE),

plot2 =  octanol_data %>% dplyr::filter(genotype == 'N2') %>%
  plot_plasticityIndex(xvar = strain,
                       plotColors = c('#827E7E','#C93030', '#2F8A34','#484CC7','#B8B1B1', '#E6DDDD', 'gray94')),

repellents = plot1 + plot2

)
