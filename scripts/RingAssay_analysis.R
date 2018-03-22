library(tidyverse)

ring <- read.csv(here::here("data","Fruc_ring.csv"))

ring %>%
  filter(odorant == "none") %>%
  ggplot(aes(x = food, y = Nout/Ntotal)) +
  geom_point(aes(color = food)) +
  scale_color_manual(values = c("grey", "blue")) +
  theme_black() +
  stat_summary(aes(y=Nout/Ntotal),
               fun.y = mean,
               fun.ymin = mean,
               fun.ymax = mean,
               geom = "crossbar",
               width = 0.5,
               lwd = 0.35,
               colour = "white") +
  facet_grid(~genotype) +
  add.n.categorical(dark)
