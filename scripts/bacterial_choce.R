# bacterial choice assay
library(tidyverse)
library(ProvidenciaChemo)
library(bayesplot)
library(modelr)
library(tidybayes)
theme_set(theme_my)

choice <- read_csv(here::here("data/Bac_choice.csv")) %>%
  filter(spacing == "near",
         is.na(note)) %>%
  mutate(food = fct_relevel(food, c("OP50","JUb39")),
         spacing = fct_relevel(spacing, "near"),
    index = (N_Test - N_OP50) / (N_Test + N_OP50), plateID = factor(seq(1:nrow(.))),
    p = N_Test / (N_Test + N_OP50),
    logit_p = boot::logit(p),
    strain = food,
    data_type = "raw") %>%
  mutate(dataset = case_when(
    date %in% c('2019_07_12', '2019_07_13', '2019_07_22') ~ 'new_data',
    TRUE ~ 'old_data'
  ))

new_data <- filter(choice, date %in% c('2019_07_12', '2019_07_13', '2019_07_22'))

new_data %>% ggplot(aes(x = data_type, y = logit_p)) +
  geom_bardots(fillvar = strain, dotvar = strain) +
  #ggbeeswarm::geom_quasirandom(aes(colour = strain), width = 0.2, alpha = 0.75) +
  facet_grid(~test_bac + strain, scales = "free_x") +
  scale_x_discrete(labels = function(strain) str_wrap(strain, width = 10)) +
  scale_color_plot("grey-blue-green", drop = TRUE) +
  scale_fill_plot("grey-blue-green", drop = TRUE) +
  labs(y = "Providencia preference index") +
  guides(colour = FALSE, fill = FALSE)

choice %>%
  ggplot(aes(x = data_type, y = p)) +
  ggbeeswarm::geom_quasirandom(aes(colour = strain), width = 0.2, alpha = 0.75) +
  facet_grid(~test_bac + strain, scales = "free_x") +
  scale_x_discrete(labels = function(strain) str_wrap(strain, width = 10)) +
  scale_color_plot("grey-blue-green", drop = TRUE) +
  add.mean(p, colour = "red") +
  add.quartiles(p) +
  labs(y = "Providencia preference index") +
    stat_pointinterval(aes(y=.value, x = 1.2),
                       data = fitted, fatten_point = 0,
                       size_range = c(0.3, 1), colour = "grey") +
    stat_summary(data = fitted,
                 aes(y=.value, x = 1.2),
                 fun.y = median,
                 fun.ymin = median,
                 fun.ymax = median,
                 geom = "crossbar",
                 width = 0.05,
                 lwd = 0.35,
                 colour = "grey") +
  guides(colour = FALSE) #+
  #scale_color_plot(palette = "2-Ps", drop = TRUE)

  plot <- choice %>%
    filter(is.na(note)) %>%
    ggplot(aes(x = data_type, y = logit_p)) +
    #stat_summary(geom = "bar", fun.y = mean, aes(fill = strain), width = 0.5, alpha = 0.75) +
    # add.mean(logit_p, colour = "red") +
    # add.quartiles(logit_p) +
    #ggbeeswarm::geom_quasirandom(aes(colour = strain), width = 0.1, alpha = 0.75) +
    geom_bardots(fillvar = strain, dotvar = strain) +
    #stat_summary(geom = "errorbar", fun.data = mean_se, width = 0.3) +
    facet_grid(. ~ test_bac + strain, scales = "free_x") +
    scale_x_discrete(labels = function(strain) str_wrap(strain, width = 10)) +
    labs(y = "Test bacteria preference (log-odds)") +
    #theme(panel.spacing = unit(4, "lines")) +
  #stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  guides(colour = FALSE) +
  scale_color_plot(palette = "grey-blue-green", drop = TRUE) +
    scale_fill_plot(palette = "grey-blue-green", drop = TRUE) +
    stat_pointinterval(aes(y=logit_p, x = 1.4),
                       data = fitted, fatten_point = 0,
                       size_range = c(0.3, 1), colour = "grey") +
    stat_summary(data = fitted,
                 aes(y=logit_p, x = 1.4),
                 fun.y = median,
                 fun.ymin = median,
                 fun.ymax = median,
                 geom = "crossbar",
                 width = 0.05,
                 lwd = 0.35,
                 colour = "grey") +
    guides(fill = FALSE) +
    #figure.axes() +
    add.n(data_type)

choice %>%
  ggplot(aes(y = index, x = N_OP50 + N_Test)) +
  geom_point(aes(colour = factor(date))) + facet_wrap(test_bac~strain, scales = "free_x") +
  geom_smooth(method = "lm")


#mark outliers:
lin_mod <- lm(data = choice, index ~ strain * test_bac)

choice_noOutliers <- flag_outliers(lin_mod, choice, threshold = 4) %>% filter(outlier.status == FALSE)

glmm_mod1 <- lme4::glmer(data = choice %>%
                          mutate(strain = fct_relevel(strain, "JUb39")),
            cbind(N_Test, N_OP50) ~ strain + test_bac + (1|plateID), family = binomial) #%>%
glmm_mod2 <- lme4::glmer(data = choice %>%
                           mutate(strain = fct_relevel(strain, "JUb39")),
                         cbind(N_Test, N_OP50) ~ strain*test_bac + (1|plateID), family = binomial) #%>%

glmm_mod3 <- lme4::glmer(data = filter(choice, test_bac == "JUb39, 2KO") %>%
                           mutate(strain = fct_relevel(strain, "JUb39")),
                         cbind(N_Test, N_OP50) ~ strain + (1|plateID), family = binomial)
  glmm_mod2 %>% emmeans::emmeans(~ strain | test_bac) %>% emmeans::contrast("trt.vs.ctrl")


stan_mod <- rstanarm::stan_glmer(data = choice_noOutliers,
                     cbind(N_Test, N_OP50) ~ strain*test_bac + (1|plateID) + (strain|date),
                     family = binomial,
                     cores = 6,
                     chains = 6,
                     adapt_delta = 0.99)

stan_mod %>% emmeans::emmeans(~ strain | test_bac) %>% emmeans::contrast(method = "pairwise") %>%
  coda::as.mcmc() %>% bayesplot::mcmc_intervals()


fitted <- choice_noOutliers %>%
  data_grid(strain, test_bac) %>%
  add_fitted_draws(stan_mod, re_formula = NA) %>%
  mutate(logit_p = boot::logit(.value), data_type = "fit")

