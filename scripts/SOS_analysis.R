library(tidyverse)

cutoff <- 10 #cutoff time for binary
SOS <- read.csv(here::here("extdata","SOS_OP50_JUb39.csv")) %>%
  mutate(food = factor(food, levels = c("OP50", "JUb39", "BL21", "BL21_tdcGOF","655", "655-tdcDel")),
         genotype = factor(genotype, levels = c("N2", "tdc-1", "octr-1", "tbh-1")),
         bin.time = if_else(response.time < cutoff, 1, 0))


##### plot 1 SOS OP50 v JUb39 ####

SOS %>% dplyr::filter(cond == "off_food", genotype == "N2", food %in% c("OP50", "JUb39"), ) %>%
  ggplot(aes(x = food, y = response.time)) +
  ggbeeswarm::geom_quasirandom(aes(colour = food), width = 0.2) +
  scale_color_plot("grey-blue", drop = TRUE) +
  facet_grid(.~genotype) +
  add.mean('response.time', colour = "red") +
  add.quartiles('response.time') +
  theme_my +
  labs(x = "food",
      y = "time to reversal (s)") +
  #figure.axes(no.x = FALSE) +
  guides(color = FALSE) + facet_wrap(~date)

SOS %>%
  dplyr::filter(cond == "off_food", genotype == "N2", food %in% c("OP50", "JUb39")) %$%
  t.test(response.time ~ food)

SOS %>% dplyr::filter(cond == "off_food", genotype == "N2", food %in% c("OP50", "JUb39")) %$%
  t.test(log(response.time) ~ food)

SOS %>% dplyr::filter(cond == "off_food", genotype == "N2", food %in% c("OP50", "JUb39")) %$%
  wilcox.test(log(response.time) ~ food)

#### plot 2  tdc GOF ####
SOS %>% dplyr::filter(cond == "off_food",
                      genotype == "N2",
                      food %in% c("OP50", "JUb39", "BL21", "BL21_tdcGOF", "655", "655-tdcDel")) %>%
  ggplot(aes(x = food, y = response.time)) +
  ggbeeswarm::geom_quasirandom(aes(colour = food), width = 0.2) +
  scale_color_plot("multi-control", drop = TRUE, discrete = TRUE) +
  facet_grid(.~genotype) +
  add.mean('response.time') +
  add.quartiles('response.time') +
  theme_my_ppt +
  labs(x = "food",
       y = "time to reversal (s)") +
  figure.axes(no.x = FALSE) #+ theme_black()

tdc_ttest <- SOS %>% dplyr::filter(genotype == "N2",
                                 food %in% c("BL21",
                                             "BL21_tdcGOF")) %$%
  t.test(response.time ~ food)

#### plot 3 tdc, octr-1 ####
SOS %>% dplyr::filter(cond == "off_food",
                      genotype %in% c("N2","tdc-1","octr-1"),
                      food %in% c("OP50", "JUb39"),
                      date == "10_24_17") %>%
  ggplot(aes(x = food, y = response.time)) +
  ggbeeswarm::geom_quasirandom(aes(colour = food), width = 0.2) +
  scale_color_plot("grey-blue", drop = TRUE) +
  facet_grid(.~genotype) +
  add.median('response.time', colour = "red") +
  add.quartiles('response.time') +
  theme_my_ppt +
  labs(x = "food",
       y = "time to reversal (s)") #+
  figure.axes(no.x = FALSE) #+ theme_black()


#### plot4 tbh-1 ####
  SOS %>% dplyr::filter(cond == "off_food",
                        genotype %in% c("N2","tbh-1"),
                        food %in% c("OP50", "JUb39"),
                        date == "6_23_18") %>%
    ggplot(aes(x = food, y = response.time)) +
    ggbeeswarm::geom_quasirandom(aes(colour = food), width = 0.2) +
    scale_color_plot("grey-blue", drop = TRUE) +
    facet_grid(.~genotype) +
    add.median('response.time', colour = "red") +
    add.quartiles('response.time') +
    theme_my_ppt +
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

######### simulations #######
ggplot(filter(SOS, genotype == "N2",
              response.time != 20,
              food %in% c("OP50","JUb39")),
       aes(x = food, y = log(response.time))) +
  geom_point()
# data excluding 20s looks log-normal, let's check

lm.normal <- lm(filter(SOS, genotype == "N2", response.time != 20), formula = response.time ~ food)
lm.log <- lm(filter(SOS, genotype == "N2", response.time != 20), formula = log(response.time) ~ food)
lm.log2 <- lm(filter(SOS, genotype == "N2", response.time != 20), formula = log2(response.time) ~ food)


DHARMa::plotSimulatedResiduals(simulationOutput = DHARMa::simulateResiduals(lm.normal))
DHARMa::plotSimulatedResiduals(simulationOutput = DHARMa::simulateResiduals(lm.log))
# log looks much better on QQ-plot
# so simulate for mean (N2, OP50) = 3.3s (estimated from model) = exp(1.2)
SOS.N2OP <- dplyr::filter(SOS, genotype == "N2", food == "OP50")$response.time
g1 <- exp(rnorm(75, 1.4, 2))
g1[g1 > 20] <- 20
data.frame(group = c(rep("real", length(SOS.N2OP)),
                     rep("sim", length(g1))),
           response.time = c(SOS.N2OP, g1)) %>%
  ggplot(aes(x = group, y = response.time)) + ggbeeswarm::geom_quasirandom(width = 0.2)
#settings look good

t.sim <- function(df) {
  genotype2 = df %$% t.test(response.time~group)$p.value
  model = "t"
  Fp = NA
  Chisq.p = NA
  p.val <- data.frame(cbind(model,genotype2, Fp, Chisq.p))
  return(p.val)
}

t.log.sim <- function(df) {
  genotype2 = df %$% t.test(log(response.time)~group)$p.value
  model = "t-log"
  Fp = NA
  Chisq.p = NA
  p.val <- data.frame(cbind(model,genotype2, Fp, Chisq.p))
  return(p.val)
}

lm.sim <-   function(df) {
  modsum <- df %>% lm(formula = response.time~group) %>% summary()
  genotype2 <- as.numeric(modsum$coefficients[,4][2])
  Fp <- as.numeric(1-pf(modsum$fstatistic[1],modsum$fstatistic[2], modsum$fstatistic[3]))
  Chisq.p = NA
  model <- "anova"
  p.val <- data.frame(cbind(model, genotype2,Fp, Chisq.p))
  return(p.val)
}

#data.frame(model = "Rank sum", p.value = wilcox.test(data1, data2)$p.value)

lm.log.sim <-   function(df) {
  modsum <- df %>% lm(formula = log(response.time)~group) %>% summary()
  genotype2 <- as.numeric(modsum$coefficients[,4][2])
  Fp <- as.numeric(1-pf(modsum$fstatistic[1],modsum$fstatistic[2], modsum$fstatistic[3]))
  Chisq.p = NA
  model <- "anova_log"
  p.val <- data.frame(cbind(model, genotype2, Fp, Chisq.p))
  return(p.val)
}

wilc.sim <- function(df) {
  genotype2 = df %$% wilcox.test(response.time~group)$p.value
  model = "wilcox"
  Fp = NA
  Chisq.p = NA
  p.val <- data.frame(cbind(model,genotype2, Fp, Chisq.p))
  return(p.val)
}

simulate_SOS <- function(n1 = 75, n2 = 75, mean1 = 1.4, mean2 = mean1, var1 = 1.1, var2 = var1, ...) {
  library(magrittr)
  g1 <- exp(rnorm(75, mean1, var1))
  g1[g1 > 20] <- 20
  g2 <- exp(rnorm(75, mean2, var2))
  g2[g2 > 20] <- 20
  data <- data.frame(group = c(rep("g1", length(g1)),
                       rep("g2", length(g2))),
             response.time = c(g1, g2))
  t <- t.sim(data)
  tlog <- t.log.sim(data)
  lm <- lm.sim(data)
  lm.log <- lm.log.sim(data)
  wilc <- wilc.sim(data)

  p.val <- rbind(t, tlog, lm, lm.log, wilc)
  return(p.val)

}

get_alpha <- function(df,...) {
  # for t
  t1 <- df %>% dplyr::filter(model == "t") %>% dplyr::count(genotype2 < 0.05)
  #t2 <- df %>% dplyr::filter(model == "t") %>% dplyr::count(genotype3 < 0.05)

  tlog <- df %>% dplyr::filter(model == "t-log") %>% dplyr::count(genotype2 < 0.05)
  # for anova
  lm1 <- df %>% dplyr::filter(model == "anova") %>% dplyr::count(genotype2 < 0.05 & Fp < 0.05)
  #lm2 <- df %>% dplyr::filter(model == "anova") %>% dplyr::count(genotype3 < 0.05 & Fp < 0.05)

  # for anova_log
  lm.log1 <- df %>% dplyr::filter(model == "anova_log") %>% dplyr::count(genotype2 < 0.05 & Fp < 0.05)

  # for Wilcox test
  wilc <- df %>% dplyr::filter(model == "wilcox") %>% dplyr::count(genotype2 < 0.05)


  # for glmm
  #glmm1 <- df %>% dplyr::filter(model == "glmm") %>% dplyr::count(genotype2 < 0.05 & Chisq.p < 0.05)
  #glmm2 <- df %>% dplyr::filter(model == "glmm") %>% dplyr::count(genotype3 < 0.05 & Chisq.p < 0.05)

  #for stan
  #stan1 <- df %>% dplyr::filter(model == "stan") %>% dplyr::count(genotype2 < 0.05)
  #stan2 <- df %>% dplyr::filter(model == "stan") %>% dplyr::count(genotype3 < 0.05)

  nsim <- nrow(df)/length(levels(df$model))
  alpha = list(t1 = t1,
               tlog = tlog,
               lm1 = lm1,
               lm.log1 = lm.log1,
               wilc = wilc)

  # prop_sig <- function(df) {
  #   # number sig in contingency table = [2,2]
  #   df[2,2]/nsim
  # }
  # alpha = c(lapply(alpha, prop_sig), nsim = nsim)
  return(alpha)
}

set.seed(18263)
sims <- do.call(rbind, replicate(1000, simulate_SOS(), simplify = FALSE)) %>%
  mutate(genotype2 = as.numeric(as.character(genotype2)),
         Fp = as.numeric(as.character(Fp)))

sims %>% get_alpha()
# all are very close to .05

#### use a low n ####
set.seed(18263)
sims <- do.call(rbind, replicate(1000, simulate_SOS(n1 = 20, n2 = 20), simplify = FALSE)) %>%
  mutate(genotype2 = as.numeric(as.character(genotype2)),
         Fp = as.numeric(as.character(Fp)))

sims %>% get_alpha()
# all are very close to .05

#### unequal variance ####
set.seed(18263)
sims <- do.call(rbind, replicate(1000, simulate_SOS(n1 = 20, n2 = 20, var2 = 1.6), simplify = FALSE)) %>%
  mutate(genotype2 = as.numeric(as.character(genotype2)),
         Fp = as.numeric(as.character(Fp)))

sims %>% get_alpha()
# all are very close to .05


#### for power (n = 20), set mean2 = 2.2, exp(2.2) ~ 9s
set.seed(18264)
sims <- do.call(rbind, replicate(1000, simulate_SOS(mean2 = 1.7, n1 = 20, n2 = 20), simplify = FALSE)) %>%
  mutate(genotype2 = as.numeric(as.character(genotype2)),
         Fp = as.numeric(as.character(Fp)))

sims %>% get_alpha()


set.seed(17255)
sims <- do.call(rbind, replicate(10000, simulate_SOS(g2 = 30), simplify = FALSE)) %>%
  mutate(genotype2 = as.numeric(as.character(genotype2)),
         Fp = as.numeric(as.character(Fp)))

sims %>% get_alpha()
#also pretty conservative


# now use a mean closer to 10s ie exp(2.2) ~ 9
set.seed(23654)
sims <- do.call(rbind, replicate(10000, simulate_SOS(mean1 = 2.2), simplify = FALSE)) %>%
  mutate(genotype2 = as.numeric(as.character(genotype2)),
         Fp = as.numeric(as.character(Fp)))

sims %>% get_alpha()
# still fine



