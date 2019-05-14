library(tidyverse)

ring <- read.csv(here::here("data","Fruc_ring.csv")) %>%
  mutate(food = fct_relevel(food, "OP50"),
         plate = factor(seq(1:nrow(.))))

ring %>%
  filter(odorant == "none",
         #food %in% c("OP50", "JUb39"),
         genotype == "N2") %>%
  mutate(fraction_out =  Nout/Ntotal) %>%
  ggplot(aes(x = food, y = fraction_out)) +
  ggbeeswarm::geom_quasirandom(aes(color = food)) +
  scale_color_plot(palette = "2-Ps", drop = TRUE) +
  #theme_black() +
  stat_summary(aes(y=fraction_out),
               fun.y = median,
               fun.ymin = median,
               fun.ymax = median,
               geom = "crossbar",
               width = 0.5,
               lwd = 0.35,
               colour = "black") +
  facet_grid(~genotype) +
  guides(colour = FALSE) +
  #add.n(food) +
  add.quartiles(fraction_out)

ring %>%
  filter(odorant == "none",
         #food %in% c("OP50", "JUb39"),
         genotype == "N2") %>%
  lme4::glmer(data = ., cbind(Nout, Ntotal - Nout) ~ food + (1|plate), family = "binomial") %>% summary()


