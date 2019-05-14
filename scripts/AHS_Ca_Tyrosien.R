# analyzing OP50 JUb39 fructose reponses on NGM plates:
#non-background subtracted:
theme_set(theme_classic())

OP50 <- MF.matR::plotGCaMP_multi("OP50", cue = "180mM_fructose",
                                 neuron = ASH,
                                 backsub = FALSE,
                                 food = OP50,
                                 matlab = FALSE,
                                 show.plots = FALSE,
                                 heatmap_limits = c(-.3, 0, 1),
                                 center_on_pulse = "ON")
JUb39 <- MF.matR::plotGCaMP_multi("JUb39_",
                                  cue = "180mM_fructose",
                                  neuron = ASH,
                                  backsub = FALSE,
                                  food = JUb39,
                                  matlab = FALSE,
                                  show.plots = FALSE,
                                  heatmap_limits = c(-.3, 0, 1),
                                  center_on_pulse = "ON")
JU2KO <- MF.matR::plotGCaMP_multi("JUb392KO",
                                  cue = "180mM_fructose",
                                  neuron = ASH,
                                  backsub = FALSE,
                                  food = JUb39_2KO,
                                  matlab = FALSE,
                                  show.plots = FALSE,
                                  heatmap_limits = c(-.3, 0, 1),
                                  center_on_pulse = "ON")

#trace plot
p1 <- rbind(OP50$data, JUb39$data, JU2KO$data) %>%
  unnest() %>%
  filter(delF < 2) %>%
  mutate(food = fct_relevel(food, "OP50"),
         date = substr(animal, 1,10)) %>%
  ggplot(aes(x = time, y = delF)) +
  geom_line(aes(group = animal, colour = food), alpha = 0.1) +
  geom_smooth(aes(colour = food), method = "loess", span = 0.05) +
  ProvidenciaChemo::scale_color_plot(palette = "grey-blue", drop = TRUE) +
  theme_classic() +
  geom_segment(x = 30, xend = 60, y = 0.6, yend = 0.6) +
  coord_cartesian(ylim = c(-.01, 1)) +
  labs(y = paste0(expression(delta),"F/F"),
       color = "grown on:",
       x = "time (s)")


#mean and SEM trace plot
p2 <- rbind(OP50$data, JUb39$data, JU2KO$data) %>%
  #rbind(OP50$data, JUb39$data) %>%
  unnest() %>%
  filter(delF < 2) %>%
  group_by(food, time) %>%
  summarize(mean = mean(delF),
            n = n(),
            sem_low = mean - sd(delF)/n^0.5,
            sem_high = mean + sd(delF)/n^0.5) %>%
  ungroup() %>%
  mutate(food = fct_relevel(food, "OP50")) %>%
  ggplot(aes(x = time)) +
  geom_line(aes(y = mean, colour = food)) +
  geom_ribbon(aes(ymin = sem_low, ymax = sem_high, fill = food), alpha = 0.2) +
  ProvidenciaChemo::scale_color_plot(palette = "grey-blue", drop = TRUE) +
  ProvidenciaChemo::scale_fill_plot(palette = "grey-blue", drop = TRUE)

#heat maps
heatmaps <- rbind(OP50$data, JUb39$data, JU2KO$data) %>%
  unnest() %>%
  mutate(food = fct_relevel(food, "OP50"),
                date = substr(animal, 1,10)) %>%
  group_by(animal,food) %>%
  ggplot(aes(x = time, y = fct_reorder(animal,
                                       maxD))) +
  geom_tile(aes(fill = delF)) +
  scale_fill_viridis_c(option = "magma",
                       breaks = c(-.05,0,1),
                       labels = as.character(c(-.1,0,1)),
                       limits = c(-.05,1),
                       oob = squish) +
  theme_classic() +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text.y = element_blank()) +
  labs(y = "Animal number") + facet_wrap(~food, scales = 'free_y')

#### histogram of amplitudes ####
rbind(OP50$data, JUb39$data, JU2KO$data) %>%
#rbind(OP50$data, JUb39$data) %>%
  unnest() %>%
  mutate(food = fct_relevel(food, "OP50"),
         date = substr(animal, 1,10)) %>%
  group_by(animal, food, maxD) %>% nest() %>%
  ggplot(aes(x = food, y = maxD)) +
  ggbeeswarm::geom_quasirandom(width = 0.05) +
  add.mean(maxD) +
  add.quartiles(maxD)

  rbind(OP50$data, JUb39$data, JU2KO$data) %>%
    unnest() %>%
    mutate(food = fct_relevel(food, "OP50"),
           date = substr(animal, 1,10)) %>%
    group_by(animal, food, maxD) %>% nest() %>%
    filter(maxD > 0) %>%
    #rstanarm::stan_glm(data = ., maxD ~ food, family = gaussian) %>% plot()
    lm(data = ., maxD ~ food) %>% summary

