octr1_Ps <- MF.matR::plotGCaMP_multi(FileFilter = "Ps", gentoype = octr1, center_on_pulse = "ON", cue = octanol, food = Ps, show.plots = FALSE)
octr1_OP <- MF.matR::plotGCaMP_multi(FileFilter = "OP50", gentoype = octr1, center_on_pulse = "ON", cue = octanol, food = OP50, show.plots = FALSE)

octr1 <- rbind(octr1_OP$data,octr1_Ps$data) %>% unnest %>% group_by(food, animal, animal_num, maxD) %>%
  nest()

octr1 %>% unnest() %>% group_by(animal, animal_num, food, maxD) %>%
  nest() %>% ggplot(aes(x = food, y = maxD)) + ggbeeswarm::geom_quasirandom(width = 0.2) +
  ProvidenciaChemo::add.mean("maxD") +
  ProvidenciaChemo::add.quartiles("maxD") +
  ProvidenciaChemo::theme_my_ppt +
  ProvidenciaChemo::figure.axes()

octr1 %>% unnest() %>% group_by(animal, animal_num, food, maxD) %>%
  nest() %$% wilcox.test(maxD ~ food)

octr1_OP$plot + ProvidenciaChemo::theme_my_ppt


N2_Ps <- MF.matR::plotGCaMP_multi(genotype = "N2_Ps", FileFilter = "Ps", cue = octanol, center_on_pulse = "ON", food = JUb39)
N2_OP <- MF.matR::plotGCaMP_multi(genotype = "N2_OP",FileFilter = "OP", cue = octanol, center_on_pulse = "ON", food = OP50)

N2 <- rbind(N2_OP$data,N2_Ps$data)
N2 %>% mutate(maxD = case_when(
  maxD < 0 ~ 0,
  TRUE ~ maxD)) %>%
  unnest() %>% group_by(animal, animal_num, food, maxD) %>%
  nest() %>% ggplot(aes(x = fct_relevel(food, "OP50"), y = maxD)) +

  ggbeeswarm::geom_quasirandom(width = 0.1) +
  ProvidenciaChemo::add.mean("maxD") +
  ProvidenciaChemo::add.quartiles("maxD") +
  #ProvidenciaChemo::add.n.categorical(group = "food") +
  ProvidenciaChemo::theme_my_ppt #+
  #ProvidenciaChemo::figure.axes()

#barplot
N2.maxD <- N2 %>%
  mutate(maxD = case_when(
  maxD < 0 ~ 0,
  TRUE ~ maxD)) %>%
  unnest() %>% mutate(
    food = factor(food, levels = c("OP50", "JUb39"))) %>%
  group_by(animal, animal_num, food, maxD) %>%
  nest()

sciplot::bargraph.CI(data = N2.maxD,
                     x.factor = food,
                     response = maxD,
                     ylim = c(0,1.25),
                     col=c("#827E7E","#484CC7"),
                     err.width = 0.05)

N2 %>% unnest() %>% group_by(animal, animal_num, food, maxD) %>%
  ggplot(aes(x = time, y = delF, color = fct_relevel(food, "OP50"))) +
  geom_line(aes(group = animal), alpha = 0.2) +
  geom_smooth(method = "loess", span = 0.1) +
  ProvidenciaChemo::theme_my_ppt +
  ProvidenciaChemo::scale_color_plot("grey-blue", drop = TRUE) +
  guides(color = FALSE)

octr1 %>% unnest() %>% group_by(animal, animal_num, food, maxD) %>%
  ggplot(aes(x = time, y = delF, color = fct_relevel(food, "OP50"))) +
  geom_line(aes(group = animal), alpha = 0.2) +
  geom_smooth(method = "loess", span = 0.1) +
  ProvidenciaChemo::theme_my_ppt +
  ProvidenciaChemo::scale_color_plot("grey-blue", drop = TRUE) +
  guides(color = FALSE)

sciplot::bargraph.CI(data = octr1,
                     x.factor = food,
                     response = maxD,
                     ylim = c(0,1.5),
                     col=c("#827E7E","#484CC7"),
                     err.width = 0.05)

N2 %>% mutate(maxD = case_when(
  maxD < 0 ~ 0,
  TRUE ~ maxD)) %>%
    unnest() %>%
    group_by(animal, animal_num, food, maxD) %>%
  nest() %$% wilcox.test(maxD ~ food)
