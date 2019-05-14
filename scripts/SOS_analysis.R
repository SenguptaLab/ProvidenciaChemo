library(tidyverse)
library(magrittr)
library(ProvidenciaChemo)
library(emmeans)

theme_set(theme_my)

cutoff <- 10 #cutoff time for binary
##### data read ####
SOS <- readr::read_csv(here::here("extdata","SOS_OP50_JUb39.csv")) %>%
  mutate(food = factor(food, levels = c("OP50",
                                        "JUb39",
                                        "BL21",
                                        "BL21_tdcGOF",
                                        "BL21_tdcK>A",
                                        "655",
                                        "655-tdcDel",
                                        "JUb39; tdcDel::cmR",
                                        "JUb39; tdcopDel",
                                        "JUb39; tdcDel::cmR delAADC",
                                        "JUb39; delAADC")),
         genotype = factor(genotype, levels = c("N2",
                                                "tbh-1",
                                                "tdc-1",
                                                "octr-1",
                                                "tyra-2",
                                                "tyra-2; ex[sra6p::tyra-2]",
                                                "tyra-2; octr-1",
                                                "octr-1; ex[sra6p::octr-1]",
                                                "tdc-1; tyra-2",
                                                "tdc-1; tbh-1")),
         bin.time = if_else(response.time < cutoff, 1, 0),
         paralyzed = case_when(response.time == 20 ~ "paralyzed",
                               TRUE ~ "non-paralyzed")) %>%
  group_by(date, genotype, food, cond, cue) %>%
  tidyr::nest() %>% mutate(plateID = factor(1:nrow(.))) %>% tidyr::unnest()



##### plot 1 SOS OP50 v JUb39 ####

SOS_OP_JUb39 <- SOS %>% dplyr::filter(cond == "off_food",
                      genotype %in% c("N2", "tdc-1"),
                      #date == "2019_01_14",
                      response.time < 20,
                      food %in% c("OP50", "JUb39")) #%>%
(p1 <- ggplot(SOS_OP_JUb39, aes(x = food, y = response.time)) +
  ggbeeswarm::geom_quasirandom(aes(colour = food), width = 0.2, alpha = 0.75) +
  scale_color_plot("grey-blue", drop = TRUE) +
  facet_grid(.~genotype) +
  add.median('response.time', colour = "red") +
  add.quartiles('response.time') +
  theme_my +
  labs(x = "food",
      y = "time to reversal (s)") +
  #figure.axes(no.x = FALSE) +
  guides(color = FALSE) +
    facet_wrap(~genotype) +
    add.n(food)) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))

SOS_OP_JUb39 %$%
  t.test(log(response.time) ~ food)

SOS_OP_JUb39 %$% wilcox.test(log(response.time) ~ food) #need to deal with non-independence of worms on each plate

#### plot 2  tdc GOF ####
GOF_data <- SOS %>% dplyr::filter(cond == "off_food",
                                  genotype == "N2",
                                  food %in% c("OP50", "JUb39", "BL21", "BL21_tdcGOF", "BL21_tdcK>A"),
                                  !date %in% c("10_17_17", "2018_11_13"))

GOF_data %>%
  ggplot(aes(x = food, y = log(response.time))) +
  ggbeeswarm::geom_quasirandom(aes(colour = date), width = 0.2) +
  #scale_color_plot("multi-control", drop = TRUE, discrete = TRUE) +
  facet_grid(.~genotype) +
  #add.mean('response.time') +
  #add.quartiles('response.time') +
  #scale_y_continuous(trans = "log") +
  labs(x = "food",
       y = "time to reversal") +
  figure.axes(no.x = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+
  guides(color = FALSE)

p1 %+% GOF_data

# model is lognormal except for threshold 20 values
library(bayesplot)
GOF_mod <- rstanarm::stan_glmer(data = GOF_data,
                              formula = log(response.time) ~ food + (1|date),
                              family = gaussian)
posterior <- as.matrix(GOF_mod)

bayesplot::mcmc_areas(posterior, regex_pars = "food", prob = c(0.8, 0.99)) + ggtitle("posterior intervals", "with medians and 80% intervals")
bayesplot::color_scheme_set("green")
bayesplot::ppc_dens_overlay(y = GOF_mod$y,
                 yrep = rstanarm::posterior_predict(GOF_mod, draws = 50))

GOF_mod %>% rstanarm::posterior_predict(draws = 500) %>%
  ppc_stat_grouped(y = log(GOF_data$response.time),
                   group = GOF_data$food,
                   stat = "median")



tdc_ttest <- SOS %>% dplyr::filter(genotype == "N2",
                                 food %in% c("BL21",
                                             "BL21_tdcGOF")) %$%
  t.test(response.time ~ food)

#### plot 3 tdc, tbh-1, on Tyrosine ####
amine_data <- SOS %>% dplyr::filter(cond == "Tyrosine",
                      cue == "oct_1.0",
                      genotype %in% c("N2","tbh-1", "tdc-1","tdc-1; tbh-1"),
                      #date %in% c("10_24_17", "2018_12_6")) %>%
                      date %in% c("2019_01_07", "2019_01_14", "2019_04_16", "2019_04_27", "2019_04_29"),
                      response.time < 20,
                      food %in% c("OP50", "JUb39")) %>%
  mutate(tdc = case_when(genotype %in% c("N2", "tbh-1") ~ "+",
                         TRUE ~ "mutant"),
         tbh = case_when(genotype %in% c("N2", "tdc-1") ~ "+",
                         TRUE ~ "mutant"))
amine_data %>%
  ggplot(aes(x = food, y = response.time)) +
  ggbeeswarm::geom_quasirandom(aes(colour = food), width = 0.2, alpha = 0.75) +
  scale_color_plot("grey-blue", drop = TRUE) +
  facet_grid(.~genotype) +
  add.median('response.time', colour = "red") +
  add.quartiles('response.time') +
  add.n(food) +
  labs(x = "food",
       y = "time to reversal (s)")

lmer <- amine_data %>%
  lme4::lmer(data = ., log(response.time) ~ food * tdc * tbh + (1|date) + (1|plateID))

lmer2 <- amine_data %>%
  lme4::lmer(data = ., log(response.time) ~ food * genotype + (1|date) + (1|plateID))

emm.lmer <- lmer %>% emmeans(pairwise ~ food  | tdc | tbh)
lmer2 %>% emmeans(pairwise ~ food | genotype)
car::Anova(lmer)

stan_lmer <- amine_data %>%
  rstanarm::stan_glmer(data = ., log(response.time) ~ food * tdc * tbh + (1|date) + (1|plateID),
                       family = gaussian,
                       chains = 6,
                       cores = 6)

emmeans(stan_lmer, pairwise ~ food | tdc * tbh)


#### plot4 tbh-1 and tdc-1 on NGM####
  SOS %>% dplyr::filter(cond == "off_food",
                        genotype %in% c("N2","tbh-1", "tdc-1"),
                        food %in% c("OP50", "JUb39"),
                        date %in% c("6_23_18", "10_24_17", "20180904")) %>%
    ggplot(aes(x = food, y = response.time)) +
    ggbeeswarm::geom_quasirandom(aes(colour = food), width = 0.2) +
    scale_color_plot("grey-blue", drop = TRUE) +
    facet_grid(.~genotype) +
    add.median('response.time', colour = "red") +
    add.quartiles('response.time') +
    theme_my +
    labs(x = "food",
         y = "time to reversal (s)") #+
  figure.axes(no.x = FALSE) #+ theme_black()

SOS.octr <- SOS %>% dplyr::filter(cond == "off_food", date == "10_24_17") %>%
  mutate(cell.norm = response.time,
         food = factor(food, levels = c("OP50", "JUb39")))
lm <- lm(data = SOS.octr, formula = response.time ~ food + genotype)
lm1 <- update(lm, formula = .~. + food*genotype)
lsmeans::lsmeans(lm1, pairwise ~ food | genotype)
glm <- lm(data = SOS.octr, formula = log(response.time) ~ food * genotype)

SOS.tdc <- SOS %>% dplyr::filter(cond == "off_food" & genotype %in% c("N2", "tdc-1"))

###Ps tdc knockout on Tyrosine ####
SOSKO <-  SOS %>% dplyr::filter(cond == "Tyrosine",
                                cue == "oct_1.0",
                      genotype %in% c("N2"),
                      food %in% c("OP50", "JUb39", "JUb39; tdcDel::cmR", "JUb39; tdcDel::cmR delAADC", "JUb39; delAADC"),
                      response.time < 20,
                      #date %in% c("2019_02_13","2019_02_15"),
                      date %in% c("2019_01_07", #need to subset days
                                  "2019_01_14",
                                  "2018_10_23",
                                  "2018_10_29",
                                  "2018_11_13",
                                  "2019_02_13",
                                  "2019_02_15",
                                  "2019_03_06",
                                  "2019_03_19",
                                  '2019_04_17',
                                  '2019_04_23',
                                  '2019_04_27',
                                  '2019_04_29')) %>%
  mutate(food = fct_relevel(food, c("OP50", "JUb39", "JUb39; tdcDel::cmR", "JUb39; delAADC")))%>%
  droplevels() %>%
  mutate(log.response = log(response.time))


# plot
p1 %+% SOSKO +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_color_plot("grey-blue-light", drop = TRUE)

(lmer <- SOSKO %>%
    lme4::lmer(., formula = log10(response.time) ~ food  + (1|date) + (1|plateID)))

  (lmer %>%
  emmeans::ref_grid() %>%
  emmeans::contrast(method = "pairwise", type = "response"))

library(rstanarm)

stan_lm <- SOSKO %>%
  rstanarm::stan_lmer(., formula = response.time ~ food + (1 | date) + (1|plateID))

sjstats::equi_test(stan_lm)

stan_lm %>%
  emmeans::ref_grid() %>%
  emmeans::contrast(method = "pairwise", type = "response") %>%
  #emmeans::emmeans(type = "response")
  #emmeans::contrast(method = "pairwise", type = "response") %>%
  coda::as.mcmc() %>%
  bayesplot::mcmc_areas(prob = 0.80, # 80% intervals
                        prob_outer = 0.95, # 95%
                        point_est = "median") +
  geom_vline(xintercept = c(-4:4), alpha = 0.1)

#bayesplot::color_scheme_set("purple")
# plot Bayes HDI
stan_lm %>% emmeans::ref_grid() %>% coda::as.mcmc() %>%
  # bayesplot::mcmc_areas_ridges(prob = 0.80,
  #                       prob_outer = 0.95,
  #                       point_est = "median") +
  # bayesplot::mcmc_areas(prob = 0.80,
  #                       prob_outer = 0.95,
  #                       point_est = "median") +
  bayesplot::mcmc_intervals(prob = 0.80,
                        prob_outer = 0.95,
                        point_est = "median") +
  # bayesplot::mcmc_areas_ridges_data(prob = 0.80,
  #                              prob_outer = 0.95,
  #                              point_est = "median") +
  scale_x_continuous(limits = c(0,20)) +
  ggbeeswarm::geom_quasirandom(data = SOSKO,
                               aes(x = response.time,
                                   y =as.numeric(forcats::fct_rev(food))-0.2,
                                   colour = food), alpha = 0.5, groupOnX = FALSE, width = 0.1, method = "smiley") +
  scale_color_plot("grey-blue", drop = TRUE)

stan_mcmc_data <- stan_lm %>%  emmeans::ref_grid() %>% coda::as.mcmc() %>%
  bayesplot::mcmc_areas_data(prob = 0.80,
                               prob_outer = 0.95) %>%
  tidyr::separate(parameter, into = c("parameter", "food"), sep = "food ") %>%
  mutate(food = factor(food, levels = levels(SOSKO$food)))

#plot bayes posterior distribution
ggplot(SOSKO) +
  ggbeeswarm::geom_quasirandom(data = SOSKO,
                               aes(x = response.time, y = -0.2, color = food),
                               alpha = 0.4, groupOnX = FALSE, width = 0.1) +
  geom_ribbon(data = stan_mcmc_data,
                aes(ymax = density, x = x, ymin = 0, fill = interval), alpha = 0.5) +

  # geom_rug(data = mutate(SOSKO, food = as.numeric(food)), aes(x = response.time), sides = "b") +
  scale_fill_plot(drop = TRUE) +
  scale_color_plot(drop = TRUE) +
  facet_grid(.~food) +
  coord_flip(ylim = c(-1, 1)) +
  theme_my



#### plot 5 Receptor and rescue ####
SOS %>% dplyr::filter(cond == "Tyrosine",
                      genotype %in% c("N2",
                                      "tyra-2",
                                      "octr-1",
                                      "tyra-2; octr-1",
                                      "octr-1; ex[sra6p::octr-1]",
                                      "tyra-2; ex[sra6p::tyra-2]",
                                      "tdc-1",
                                      "tbh-1"),
                      food %in% c("OP50", "JUb39"), #) %>%
                      #!date %in% c("2018_10_29"),
                      #date %in% c("2018_10_23","2018_11_13"), #weirdness with 10_29 + tdc-1 was contaminated before 11_13
                      date %in% c("2019_01_07", "2019_01_14"),
                      #response.time < 20,
                      !(genotype == "tyra-2" & food == "JUb39; tdcDel::cmR")) %>%
  ggplot(aes(x = food, y = response.time)) +
  ggbeeswarm::geom_quasirandom(aes(colour = food), width = 0.2, alpha = 0.3, method = "smiley") +
  scale_color_plot("grey-blue", drop = TRUE) +
  facet_grid(.~genotype, scales = "free", labeller = label_wrap_gen(width = 2)) +
  add.median('response.time', colour = "red", group = food)  +
  add.quartiles('response.time', group = food) +
  add.n(food) +
  theme_my +
  labs(x = "food",
       y = "time to reversal (s)") +
  theme(axis.title.x = ggplot2::element_blank(),
             axis.text.x = ggplot2::element_blank(),
             axis.text.y = ggplot2::element_text(size = 15),
        panel.spacing.x=unit(1, "lines")) +
  stat_summary(aes(group=date, colour = date), fun.y=median, geom="line")
#figure.axes()
  #scale_x_discrete(labels = c("control", "+Tyrosine")) +
  scale_y_log10()
figure.axes(no.x = FALSE) #+ theme_black()

#### plot octr-1 and rescue ####
SOS %>% dplyr::filter(cond == "Tyrosine",
                      genotype %in% c("N2",
                                      "tyra-2",
                                      "octr-1",
                                      "octr-1; ex[sra6p::octr-1]"),
                      food %in% c("OP50", "JUb39"), #) %>%
                      #date == c("2019_02_15"),
                      #date %in% c("2018_10_23","2018_11_13"), #weirdness with 10_29 + tdc-1 was contaminated before 11_13
                      date %in% c("2019_01_07",
                                  "2019_01_14",
                                  "2019_02_05",
                                  "2019_02_15",
                                  "2019_03_06",
                                  "2019_04_02",
                                  "2019_04_23"),
                      response.time < 20,
                      !(genotype == "tyra-2" & food == "JUb39; tdcDel::cmR")) %>%
  #mutate(response.time = log10(response.time)) %>%
  ggplot(aes(x = food, y = response.time)) +
  ggbeeswarm::geom_quasirandom(aes(colour = food), width = 0.2, alpha = 0.3, method = "smiley") +
  scale_color_plot("grey-blue", drop = TRUE) +
  facet_grid(.~genotype, scales = "free", labeller = label_wrap_gen(width = 2)) +
  add.median('response.time', colour = "red", group = food)  +
  add.quartiles('response.time', group = food) +
  add.n(food, y.pos = -.3) +
  theme_my +
  labs(x = "food",
       y = "time to reversal (s)") +
  theme(axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 15),
        panel.spacing.x=unit(1, "lines")) #+
  #stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.2) +
  #geom_halfeyeh(aes(y = response.time, group = food, fill = food, alpha = 0.2), data = mod.layer)
  #ggridges::geom_density_ridges(aes(y = response.time), data = mod.layer)
  stat_pointinterval(aes(y = response.time, x = food, colour = food),
                     data = mod.layer,
                     .width = c(.99, .95, .75),
                     alpha = 0.75)


mod <- SOS %>% dplyr::filter(cond == "Tyrosine",
                      genotype %in% c("N2",
                                      "octr-1",
                                      "octr-1; ex[sra6p::octr-1]"),
                      food %in% c("OP50", "JUb39"), #) %>%
                      #date == c("2019_02_15"),
                      #date %in% c("2018_10_23","2018_11_13"), #weirdness with 10_29 + tdc-1 was contaminated before 11_13
                      date %in% c("2019_01_07",
                                  "2019_01_14",
                                  "2019_02_05",
                                  "2019_02_15",
                                  "2019_03_06",
                                  "2019_04_02",
                                  "2019_04_23"),
                      response.time < 20,
                      !(genotype == "tyra-2" & food == "JUb39; tdcDel::cmR")) %>%
  #mutate(response.time = log(response.time)) %>%
  #lme4::lmer(data = ., response.time ~ genotype * food + (1 | plateID)) %>% emmeans::emmeans(pairwise ~ food | genotype)
  rstanarm::stan_glmer(data = ., log10(response.time) ~ genotype*food + (1|date) + (1|plateID), family = gaussian, cores = 4, chains = 4)

library(tidybayes)
library(emmeans)
library(bayesplot)
mod.layer <- mod %>% emmeans::emmeans(pairwise ~ food | genotype) %>%
  gather_emmeans_draws() %>%
  filter(.grid == "emmeans") %>%
  mutate(response.time = .value) #%>%

ggplot(mod.layer, aes(x = response.time, y = food)) +
  geom_halfeyeh(aes(group = food, fill = food, alpha = 0.2)) +
  facet_grid(~genotype)

#### tdc KO on NGM v Tyrosine ####
SOS %>% dplyr::filter(cond %in% c("off_food", "Tyrosine"),
                      genotype == "tyra-2",
                      food %in% c("OP50", "JUb39", "JUb39; tdcDel::cmR"),
                      date %in% c("2018_10_23", "2018_10_29")) %>%
  ggplot(aes(x = cond, y = response.time, fill = date)) +
  geom_boxplot(outlier.shape =NA) +
  geom_point(pch = 21, position = position_jitterdodge()) +
  #ggbeeswarm::geom_quasirandom(aes(colour = date), width = 0.2) +
  #scale_color_plot("grey-blue", drop = TRUE) +
  facet_grid(.~food) +
  # add.median('response.time', colour = "red") +
  # add.quartiles('response.time') +
  theme_my +
  labs(x = "food",
       y = "time to reversal (s)") +
  scale_x_discrete(labels = c("control", "+Tyrosine")) #+ theme_black()

#### plot 6 tyra-2 and tdc-1 ####
SOS %>% dplyr::filter(cond %in% c("Tyrosine"),
                      genotype %in% c("N2","tyra-2", "tdc-1", "tyra-2; ex[sra6p::tyra-2]", "tdc-1; tyra-2"),
                      response.time < 20,
                      food %in% c("OP50", "JUb39", "JUb39; tdcDel::cmR", "JUb39; tdcopDel")) %>%
                      #date %in% c("2018_11_27")) %>% #"2018_10_23", "2018_10_29","2018_11_13"
  ggplot(aes(x = food, y = response.time)) +
  #geom_boxplot(outlier.shape =NA) +
  #geom_point(pch = 21) +#, position = position_jitterdodge()) +
  ggbeeswarm::geom_quasirandom(aes(colour = food), width = 0.2) +
  scale_color_plot("grey-blue", drop = TRUE) +
  facet_grid(.~genotype) +
  add.median('response.time', colour = "red") +
  add.quartiles('response.time') +
  theme_my +
  labs(x = "food",
       y = "time to reversal (s)")
#scale_x_discrete(labels = c("control", "+Tyrosine")) #+
#figure.axes(no.x = FALSE) #+ theme_black()

#### plot 7 tdc-1 only ####
SOS %>% dplyr::filter(cond %in% c("off_food"),
                      genotype %in% c("tdc-1", "N2"),
                      food %in% c("OP50", "JUb39", "JUb39; tdcDel::cmR"),
                      date == "2019_04_16") %>%
                      #date %in% c("2018_11_13", "2018_11_27")) %>% #"2018_10_23", "2018_10_29",
  ggplot(aes(x = food, y = response.time)) +
  #geom_boxplot(outlier.shape =NA) +
  #geom_point(pch = 21) +#, position = position_jitterdodge()) +
  ggbeeswarm::geom_quasirandom(aes(colour = food), width = 0.2) +
  scale_color_plot("grey-blue", drop = TRUE) +
  facet_grid(.~genotype) +
  add.median('response.time', colour = "red") +
  add.quartiles('response.time') +
  theme_my +
  labs(x = "food",
       y = "time to reversal (s)")

#### for 30% octanol ####
SOS %>% dplyr::filter(cond %in% c("20min_off_food"),
                      response.time < 20,
                      genotype %in% c("tdc-1", "N2"),
                      food %in% c("OP50", "JUb39", "JUb39; tdcDel::cmR"),
                      date %in% c("2019_03_18", "2019_04_09", "2019_04_13", "04_29_19")) %>%
  ggplot(aes(x = genotype, y = response.time)) +
  #geom_boxplot(outlier.shape =NA) +
  #geom_point(pch = 21) +#, position = position_jitterdodge()) +
  ggbeeswarm::geom_quasirandom(aes(colour = food), width = 0.2) +
  scale_color_plot("grey-blue", drop = TRUE) +
  facet_grid(.~food) +
  add.mean('response.time', colour = "red") +
  add.quartiles('response.time') +
  theme_my +
  labs(x = "food",
       y = "time to reversal (s)")

#brm
SOS %>% dplyr::filter(cond %in% c("20min_off_food"),
                      response.time < 20,
                      genotype %in% c("tdc-1", "N2"),
                      food %in% c("OP50", "JUb39", "JUb39; tdcDel::cmR"),
                      date %in% c("2019_03_18", "2019_04_09", "2019_04_09")) %>%
  lme4::lmer(data = ., response.time ~ genotype*food + (1|date)) %>%
  emmeans::ref_grid() %>% emmeans::emmeans(pairwise ~ genotype | food)
  rstanarm::stan_glmer(data = ., response.time ~ genotype*food + (1|date), family = gaussian) %>% plot()

#### TA-saturation ####
TAGOF <- SOS %>%
  dplyr::filter(
    cond %in% c("Tyr+TA", "Tyrosine"),
    cue == "oct_1.0",
    response.time < 20,
    !date %in% c("2018_10_23", "2018_10_29"),
    genotype == "N2",
    food %in% c("OP50", "JUb39", "JUb39; tdcDel::cmR delAADC")
  ) %>%
  mutate(
    date_group = case_when(
      date %in% c("2019_03_19", "2019_04_17", "2019_04_23", "2019_04_27", "2019_05_06") ~ "TA-dataset",
      TRUE ~ "old_data"
    ),
    cond = fct_relevel(cond, "Tyrosine")
  ) %>%
  filter(
    date_group == "TA-dataset" #& date %in% c("2019_03_19", "2019_04_17")
  )

  TAGOF %>%
  ggplot(aes(x = cond, y = response.time)) +
  ggbeeswarm::geom_quasirandom(width = 0.2, aes(colour = food), alpha = 0.75) +
  scale_color_plot("2-Ps", drop = TRUE) +
  scale_alpha_manual(values = c(0.3, 1)) +
  facet_grid(.~food) +
  add.median('response.time', colour = "red") +
  add.quartiles('response.time') +
  add.n(cond) +
  theme_my +
  guides(colour = FALSE) +
  labs(x = "food",
       y = "time to reversal (s)")

  TAGOF %>% lm(data = ., log(response.time) ~ food*cond + date) %>%
    emmeans(pairwise ~ cond | food)

#### 2-nonanone avoidance ####
SOS %>% dplyr::filter(cond %in% c("Tyrosine"),
                      cue == "nonanone_1.0",
                      response.time < 20,
                      genotype == "N2",
                      food %in% c("OP50", "JUb39")) %>%
  ggplot(aes(x = food, y = response.time)) +
  ggbeeswarm::geom_quasirandom(width = 0.2, aes(colour = food)) +
  scale_color_plot("grey-blue", drop = TRUE) +
  scale_alpha_manual(values = c(0.3, 1)) +
  add.median('response.time', colour = "red") +
  add.quartiles('response.time') +
  add.n(food) +
  theme_my +
  labs(x = "food",
       y = "time to reversal (s)",
       title = "2-nonanone avoidance")


#### both nonanone and oct on Tyrosine together ####
  SOS_control <- SOS %>% dplyr::filter(cond %in% c("Tyrosine"),
                                       cue %in% (c("oct_1.0", "nonanone_1.0")),
                                       response.time < 20,
                                       genotype == "N2",
                                       #date %in% c("2019_03_19", "2019_04_02"),
                                       food %in% c("OP50", "JUb39")) %>%
    mutate(food = fct_relevel(food, "OP50"),
           cue = fct_relevel(cue, "oct_1.0"))


  plot1 <- SOS_control %>%
    ggplot(aes(x = food, y = response.time)) +
    ggbeeswarm::geom_quasirandom(width = 0.2, aes(colour = food)) +
    scale_color_plot("grey-blue", drop = TRUE) +
    scale_alpha_manual(values = c(0.3, 1)) +
    add.median('response.time', colour = "red") +
    add.quartiles('response.time') +
    guides(colour = FALSE) +
    #coord_cartesian(ylim = c(0,20)) +
    facet_grid(.~cue, labeller = as_labeller(c(`nonanone_1.0` = "100% nonanone",
                                             `oct_1.0` = "100% octanol"))) +
    theme(strip.text = element_text(size = 12)) +
    labs(x = "",
         y = "time to reversal (s)",
         title = "Acute avoidance")

  plot1

  lme4::lmer(SOS_control, formula = log(response.time) ~ food*cue + (1|plateID) + (1|date)) %>%
    emmeans(pairwise ~ food | cue)

#### plot the fraction of paralyzed worms by condition: ####
paralyzed_barplot <-  function(df) {
  df %>% ggplot(aes(x=food, fill=fct_relevel(paralyzed, "non-paralyzed"))) +
    geom_bar(position='fill') +
    scale_fill_manual(values=c("#CCCCCC", "#CC6666", "#990000")) +
    theme_my +
    guides(fill=guide_legend(title="proportion paralyzed")) +
  facet_grid(cond~genotype, scales = "free_x")
  }



SOS %>% paralyzed_barplot()

#### tbh-1 ####
dates <- SOS %>% filter(genotype == "tbh-1") %>% droplevels() %>% select(date) %>% unique()

SOS %>% filter(genotype %in% c("N2", "tbh-1"), date %in% dates$date) %>% paralyzed_barplot()

#### octr-1 ####
dates <- SOS %>% filter(genotype == "octr-1") %>% droplevels() %>% select(date) %>% unique()

SOS %>% filter(genotype %in% c("N2", "octr-1"), date %in% dates$date) %>% paralyzed_barplot()







##### temporal response #####




# ######### simulations #######
# ggplot(filter(SOS, genotype == "N2",
#               response.time != 20,
#               food %in% c("OP50","JUb39")),
#        aes(x = food, y = log(response.time))) +
#   geom_point()
# # data excluding 20s looks log-normal, let's check
#
# lm.normal <- lm(filter(SOS, genotype == "N2", response.time != 20), formula = response.time ~ food)
# lm.log <- lm(filter(SOS, genotype == "N2", response.time != 20), formula = log(response.time) ~ food)
# lm.log2 <- lm(filter(SOS, genotype == "N2", response.time != 20), formula = log2(response.time) ~ food)
#
#
# DHARMa::plotSimulatedResiduals(simulationOutput = DHARMa::simulateResiduals(lm.normal))
# DHARMa::plotSimulatedResiduals(simulationOutput = DHARMa::simulateResiduals(lm.log))
# # log looks much better on QQ-plot
# # so simulate for mean (N2, OP50) = 3.3s (estimated from model) = exp(1.2)
# SOS.N2OP <- dplyr::filter(SOS, genotype == "N2", food == "OP50")$response.time
# g1 <- exp(rnorm(75, 1.4, 2))
# g1[g1 > 20] <- 20
# data.frame(group = c(rep("real", length(SOS.N2OP)),
#                      rep("sim", length(g1))),
#            response.time = c(SOS.N2OP, g1)) %>%
#   ggplot(aes(x = group, y = response.time)) + ggbeeswarm::geom_quasirandom(width = 0.2)
# #settings look good
#
# t.sim <- function(df) {
#   genotype2 = df %$% t.test(response.time~group)$p.value
#   model = "t"
#   Fp = NA
#   Chisq.p = NA
#   p.val <- data.frame(cbind(model,genotype2, Fp, Chisq.p))
#   return(p.val)
# }
#
# t.log.sim <- function(df) {
#   genotype2 = df %$% t.test(log(response.time)~group)$p.value
#   model = "t-log"
#   Fp = NA
#   Chisq.p = NA
#   p.val <- data.frame(cbind(model,genotype2, Fp, Chisq.p))
#   return(p.val)
# }
#
# lm.sim <-   function(df) {
#   modsum <- df %>% lm(formula = response.time~group) %>% summary()
#   genotype2 <- as.numeric(modsum$coefficients[,4][2])
#   Fp <- as.numeric(1-pf(modsum$fstatistic[1],modsum$fstatistic[2], modsum$fstatistic[3]))
#   Chisq.p = NA
#   model <- "anova"
#   p.val <- data.frame(cbind(model, genotype2,Fp, Chisq.p))
#   return(p.val)
# }
#
# #data.frame(model = "Rank sum", p.value = wilcox.test(data1, data2)$p.value)
#
# lm.log.sim <-   function(df) {
#   modsum <- df %>% lm(formula = log(response.time)~group) %>% summary()
#   genotype2 <- as.numeric(modsum$coefficients[,4][2])
#   Fp <- as.numeric(1-pf(modsum$fstatistic[1],modsum$fstatistic[2], modsum$fstatistic[3]))
#   Chisq.p = NA
#   model <- "anova_log"
#   p.val <- data.frame(cbind(model, genotype2, Fp, Chisq.p))
#   return(p.val)
# }
#
# wilc.sim <- function(df) {
#   genotype2 = df %$% wilcox.test(response.time~group)$p.value
#   model = "wilcox"
#   Fp = NA
#   Chisq.p = NA
#   p.val <- data.frame(cbind(model,genotype2, Fp, Chisq.p))
#   return(p.val)
# }
#
# simulate_SOS <- function(n1 = 75, n2 = 75, mean1 = 1.4, mean2 = mean1, var1 = 1.1, var2 = var1, ...) {
#   library(magrittr)
#   g1 <- exp(rnorm(75, mean1, var1))
#   g1[g1 > 20] <- 20
#   g2 <- exp(rnorm(75, mean2, var2))
#   g2[g2 > 20] <- 20
#   data <- data.frame(group = c(rep("g1", length(g1)),
#                        rep("g2", length(g2))),
#              response.time = c(g1, g2))
#   t <- t.sim(data)
#   tlog <- t.log.sim(data)
#   lm <- lm.sim(data)
#   lm.log <- lm.log.sim(data)
#   wilc <- wilc.sim(data)
#
#   p.val <- rbind(t, tlog, lm, lm.log, wilc)
#   return(p.val)
#
# }
#
# get_alpha <- function(df,...) {
#   # for t
#   t1 <- df %>% dplyr::filter(model == "t") %>% dplyr::count(genotype2 < 0.05)
#   #t2 <- df %>% dplyr::filter(model == "t") %>% dplyr::count(genotype3 < 0.05)
#
#   tlog <- df %>% dplyr::filter(model == "t-log") %>% dplyr::count(genotype2 < 0.05)
#   # for anova
#   lm1 <- df %>% dplyr::filter(model == "anova") %>% dplyr::count(genotype2 < 0.05 & Fp < 0.05)
#   #lm2 <- df %>% dplyr::filter(model == "anova") %>% dplyr::count(genotype3 < 0.05 & Fp < 0.05)
#
#   # for anova_log
#   lm.log1 <- df %>% dplyr::filter(model == "anova_log") %>% dplyr::count(genotype2 < 0.05 & Fp < 0.05)
#
#   # for Wilcox test
#   wilc <- df %>% dplyr::filter(model == "wilcox") %>% dplyr::count(genotype2 < 0.05)
#
#
#   # for glmm
#   #glmm1 <- df %>% dplyr::filter(model == "glmm") %>% dplyr::count(genotype2 < 0.05 & Chisq.p < 0.05)
#   #glmm2 <- df %>% dplyr::filter(model == "glmm") %>% dplyr::count(genotype3 < 0.05 & Chisq.p < 0.05)
#
#   #for stan
#   #stan1 <- df %>% dplyr::filter(model == "stan") %>% dplyr::count(genotype2 < 0.05)
#   #stan2 <- df %>% dplyr::filter(model == "stan") %>% dplyr::count(genotype3 < 0.05)
#
#   nsim <- nrow(df)/length(levels(df$model))
#   alpha = list(t1 = t1,
#                tlog = tlog,
#                lm1 = lm1,
#                lm.log1 = lm.log1,
#                wilc = wilc)
#
#   # prop_sig <- function(df) {
#   #   # number sig in contingency table = [2,2]
#   #   df[2,2]/nsim
#   # }
#   # alpha = c(lapply(alpha, prop_sig), nsim = nsim)
#   return(alpha)
# }
#
# set.seed(18263)
# sims <- do.call(rbind, replicate(1000, simulate_SOS(), simplify = FALSE)) %>%
#   mutate(genotype2 = as.numeric(as.character(genotype2)),
#          Fp = as.numeric(as.character(Fp)))
#
# sims %>% get_alpha()
# # all are very close to .05
#
# #### use a low n ####
# set.seed(18263)
# sims <- do.call(rbind, replicate(1000, simulate_SOS(n1 = 20, n2 = 20), simplify = FALSE)) %>%
#   mutate(genotype2 = as.numeric(as.character(genotype2)),
#          Fp = as.numeric(as.character(Fp)))
#
# sims %>% get_alpha()
# # all are very close to .05
#
# #### unequal variance ####
# set.seed(18263)
# sims <- do.call(rbind, replicate(1000, simulate_SOS(n1 = 20, n2 = 20, var2 = 1.6), simplify = FALSE)) %>%
#   mutate(genotype2 = as.numeric(as.character(genotype2)),
#          Fp = as.numeric(as.character(Fp)))
#
# sims %>% get_alpha()
# # all are very close to .05
#
#
# #### for power (n = 20), set mean2 = 2.2, exp(2.2) ~ 9s
# set.seed(18264)
# sims <- do.call(rbind, replicate(1000, simulate_SOS(mean2 = 1.7, n1 = 20, n2 = 20), simplify = FALSE)) %>%
#   mutate(genotype2 = as.numeric(as.character(genotype2)),
#          Fp = as.numeric(as.character(Fp)))
#
# sims %>% get_alpha()
#
#
# set.seed(17255)
# sims <- do.call(rbind, replicate(10000, simulate_SOS(g2 = 30), simplify = FALSE)) %>%
#   mutate(genotype2 = as.numeric(as.character(genotype2)),
#          Fp = as.numeric(as.character(Fp)))
#
# sims %>% get_alpha()
# #also pretty conservative
#
#
# # now use a mean closer to 10s ie exp(2.2) ~ 9
# set.seed(23654)
# sims <- do.call(rbind, replicate(10000, simulate_SOS(mean1 = 2.2), simplify = FALSE)) %>%
#   mutate(genotype2 = as.numeric(as.character(genotype2)),
#          Fp = as.numeric(as.character(Fp)))
#
# sims %>% get_alpha()
# # still fine
#
#
#
