# bacterial choice assay
library(tidyverse)
theme_set(theme_my)
library(ProvidenciaChemo)

choice <- read_csv(here::here("data/Bac_choice.csv")) %>%
  mutate(food = fct_relevel(food, c("OP50","JUb39")),
         spacing = fct_relevel(spacing, "near"),
    index = (N_JUb39 - N_OP50) / (N_JUb39 + N_OP50), plateID = factor(seq(1:nrow(.))))

choice %>% ggplot(aes(x = food, y = index)) +
  add.median(index, colour = "red") +
  add.quartiles(index) +
  geom_point(aes(colour = food)) + facet_wrap(~spacing) +
  scale_x_discrete(labels = function(food) str_wrap(food, width = 10)) +
  labs(y = "Providencia preference index") +
  theme(panel.spacing = unit(4, "lines")) +
  guides(colour = FALSE) +
  scale_color_plot(palette = "2-Ps", drop = TRUE)


lme4::glmer(data = choice, cbind(N_JUb39, N_OP50) ~ food*spacing + (1|plateID), family = binomial) %>%
  emmeans::emmeans(pairwise ~ food | spacing)

lm(data = choice, index ~ food*spacing) %>%
  emmeans::emmeans(pairwise ~ food | spacing)

rstanarm::stan_glmer
