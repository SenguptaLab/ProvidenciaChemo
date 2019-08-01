OPmeans <- amine_data %>%
  mutate(response.time = case_when(
    response.time < 1 ~ 1,
    TRUE ~ response.time)) %>%
  group_by(genotype, food, date) %>%
  summarise(meanOP = mean(log(response.time))) %>%
  filter(food == "OP50") %>%
  ungroup() %>%
  select(genotype, date, meanOP)


plot <- full_join(amine_data, OPmeans) %>%
  mutate(
    rel_log = log(response.time) - meanOP,
    food = fct_relevel(food, "OP50")) %>%
  ggplot(aes(x = data_type, y = rel_log)) +
  stat_summary(geom = "bar", fun.y = mean, aes(fill = food), width = 0.75, alpha = 0.5) +
  ggbeeswarm::geom_quasirandom(aes(colour = food), width = 0.2, alpha = 0.75) +
  scale_color_plot("grey-blue", drop = TRUE) +
  scale_fill_plot("grey-blue", drop = TRUE) +
  facet_grid(.~genotype + food,
             labeller = labeller(genotype = label_wrap_gen(10),
                                 food = as_labeller(c("OP50" = "",
                                                      "JUb39" = "")))) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.2) +
  add.n('data_type', y.pos = -1) +
  # stat_pointinterval(aes(y=response.time, x = 1.5),
  #                    data = fitted, fatten_point = 0,
  #                    size_range = c(0.3, 1), colour = "grey") +
  # stat_summary(data = fitted,
  #              aes(y=response.time, x = 1.5),
  #              fun.y = median,
  #              fun.ymin = median,
  #              fun.ymax = median,
  #              geom = "crossbar",
  #              width = 0.05,
  #              lwd = 0.35,
  #              colour = "grey") +
  labs(x = "genotype",
       y = "time to reversal (s)") +
  guides(colour = FALSE) +
  theme(axis.text.x = element_blank())

#emmeans(stan_mod, pairwise ~ (food | genotype))$contrasts %>% c
#oda::as.mcmc() %>% bayesplot::mcmc_intervals(prob = 0.66, prob_outer = 0.95)


plot <- full_join(amine_data, OPmeans) %>%
  mutate(
    rel_log = log(response.time) - meanOP,
    food = fct_relevel(food, "OP50"),
    genotype = fct_relevel(genotype, c("N2", "tdc-1", "tbh-1"))) %>%
  ggplot(aes(x = data_type)) +
  geom_relLatency(fitted = fitted1,
                  fillvar = food,
                  dotvar = food,
                  yvar = rel_log) +
  scale_color_plot("grey-blue", drop = TRUE) +
  scale_fill_plot("grey-blue", drop = TRUE) +
  facet_grid(.~genotype+food) +
  add.n('data_type', y.pos = -1.6) +
  labs(x = "genotype",
       y = "relative reversal latency [log(s)]") +
  guides(colour = FALSE) +
  theme(axis.text.x = element_blank())

gt <- ggplot_gtable(ggplot_build(plot))
gt$widths[c(8,12,16,20)] = 5*gt$widths[c(8,12,16,20)]
gt$widths[c(6,10,14,18,18)] = .3*gt$widths[c(6,10,14,16,18)]
grid::grid.draw(gt)


#tyrosine data
OPmeans <- tyrosine_data %>%
  mutate(response.time = case_when(
    response.time < 1 ~ 1,
    TRUE ~ response.time)) %>%
  group_by(food, date, cond) %>%
  summarise(meanOP = mean(log(response.time))) %>%
  filter(food == "OP50") %>%
  ungroup() %>%
  select(cond, date, meanOP)


plot <- mutate(tyrosine_data, response.time = case_when(
  response.time < 1 ~ 1,
  TRUE ~ response.time)) %>%
  full_join(., OPmeans) %>%
  mutate(
    rel_log = log(response.time) - meanOP,
    food = fct_relevel(food, "OP50")) %>%
  ggplot(aes(x = data_type, y = rel_log)) +
  ggbeeswarm::geom_quasirandom(aes(colour = food), width = 0.2, alpha = 0.75) +
  geom_boxplot(aes(fill = food), alpha = 0.5)+
  scale_color_plot("grey-blue", drop = TRUE) +
  scale_fill_plot("grey-blue", drop = TRUE) +
  facet_grid(.~cond + food,
             labeller = labeller(genotype = label_wrap_gen(10),
                                 food = as_labeller(c("OP50" = "",
                                                      "JUb39" = "")))) +
  add.n('data_type', y.pos = -1) +
  labs(x = "genotype",
       y = "time to reversal (s)") +
  guides(colour = FALSE) +
  add.mean(rel_log, colour = "red") +
  theme(axis.text.x = element_blank())

gt <- ggplot_gtable(ggplot_build(plot))
gt$widths[c(8,12)] = 2.5*gt$widths[c(8,12)]
gt$widths[c(6,10,14)] = .3*gt$widths[c(6,10,14)]
grid::grid.draw(gt)

#tyrosine data
OPmeans <- receptors %>%
  mutate(response.time = case_when(
    response.time < 1 ~ 1,
    TRUE ~ response.time)) %>%
  group_by(food, date, genotype) %>%
  summarise(meanOP = mean(log(response.time))) %>%
  filter(food == "OP50") %>%
  ungroup() %>%
  select(genotype, date, meanOP)


plot <- mutate(receptors, response.time = case_when(
  response.time < 1 ~ 1,
  TRUE ~ response.time)) %>%
  full_join(., OPmeans) %>%
  mutate(
    rel_log = log(response.time) - meanOP,
    food = fct_relevel(food, "OP50")) %>%
  ggplot(aes(x = data_type, y = rel_log)) +
  ggbeeswarm::geom_quasirandom(aes(colour = food), width = 0.2, alpha = 0.75) +
  geom_boxplot(aes(fill = food), alpha = 0.5)+
  scale_color_plot("grey-blue", drop = TRUE) +
  scale_fill_plot("grey-blue", drop = TRUE) +
  facet_grid(.~genotype + food,
             labeller = labeller(genotype = label_wrap_gen(10),
                                 food = as_labeller(c("OP50" = "",
                                                      "JUb39" = "")))) +
  add.n('data_type', y.pos = -1) +
  labs(x = "genotype",
       y = "time to reversal (s)") +
  guides(colour = FALSE) +
  add.mean(rel_log, colour = "red") +
  theme(axis.text.x = element_blank())

gt <- ggplot_gtable(ggplot_build(plot))
gt$widths[c(8,12)] = 2.5*gt$widths[c(8,12)]
gt$widths[c(6,10,14)] = .3*gt$widths[c(6,10,14)]
grid::grid.draw(gt)

