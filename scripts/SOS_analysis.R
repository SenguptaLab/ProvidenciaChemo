library(tidyverse)

cutoff <- 10 #cutoff time for binary
SOS <- read.csv(here::here("data","SOS_OP50_JUb39.csv")) %>%
  mutate(cell.norm = response.time,
         food = factor(food, levels = c("OP50", "JUb39")),
         bin.time = if_else(response.time < cutoff, 1, 0))


SOS %>% dplyr::filter(cond == "off_food", date == "10_24_17") %>%
  ggplot(aes(x = food, y = response.time)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.25, aes(fill = food)) +
  scale_fill_manual(values = c("grey", "blue")) +
  facet_grid(.~genotype) + add.median() + add.mean() + theme_my +
  labs(x = "food",
       y = "time to reversal (s)") + theme_black()

SOS %>% dplyr::filter(cond == "off_food", date == "10_24_17", response.time != "20") %>%
  ggplot(aes(x = food, y = response.time)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.25, aes(fill = food)) +
  scale_fill_manual(values = c("grey", "blue")) +
  facet_grid(.~genotype) + add.median() + add.mean() + theme_my +
  labs(x = "food",
       y = "time to reversal (s)") + theme_black()

SOS.octr <- SOS %>% dplyr::filter(cond == "off_food", date == "10_24_17") %>%
  mutate(cell.norm = response.time,
         food = factor(food, levels = c("OP50", "JUb39")))
lm <- lm(data = SOS.octr, formula = response.time ~ food + genotype)
lm1 <- update(lm, formula = .~. + food*genotype)


SOS.tdc <- SOS %>% dplyr::filter(cond == "off_food" & genotype %in% c("N2", "tdc-1"))


SOS.tdc %>%
  ggplot(aes(x = food, y = response.time)) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.3, aes(fill = food)) +
  scale_fill_manual(values = c("grey", "blue")) +
  stat_summary(aes(x=as.numeric(as.factor(food)) + 0.3, y=0),
               fun.data = fun_length, geom = "text", size = 3, colour = "white") +
  facet_grid(.~genotype) + add.mean() + theme_my +
  labs(x = "food",
       y = "time to reversal (s)") + theme_black()

glm <- glm(data = SOS.octr, formula = bin.time ~ food + genotype, family = binomial)
glm1 <- glm(data = SOS.octr, formula = bin.time ~ food * genotype, family = binomial)




