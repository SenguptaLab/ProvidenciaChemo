
library(tidyverse)

Op <- c(19,1,3,3,8,8,2,1,10,5,2,11,9,3,8,4,6,5,3)
Ps <- c(20,4,4,4,5,3,2,4,3,5,2,5,4,20,10,5,4,11,3)
top10tdc <- c(2,2,5,3,3,3,3,5,1,4,4,6,2,6,6,4,9,5,3)
top10Con <- c(2,1,3,2,5,6,4,17,4,4,5,3,2,4,5,4,2,4,6)

BL21tdc0min <- c(1,2,4,1,1,3,2,8,3,3,2,3,2,1,2,10,2,2,6,3)
BL210min <- c(1,1,2,2,1,2,6,2,2,2,2,2,1,3,3,2,2,3,1,2)

BL21tdc20min <- c(1,13,13,19,6,11,14,10,2,6,6,3,5,7,3,14,5,3,6,5)
BL2120min <- c(9,3,3,3,8,5,8,2,1,6,6,4,4,4,2,4,5,3,6,2)





data <- data_frame(
  food = c(
    rep("Op", length(Op)),
    rep("Ps", length(Ps)),
    rep("top10 tdc+", length(top10tdc)),
    rep("top10 Con", length(top10Con)),
    rep("BL21", length(BL210min)),
    rep("BL21 + tdc", length(BL21tdc0min)),
    rep("BL21 20min", length(BL2120min)),
    rep("BL21 + tdc 20min", length(BL21tdc20min))
  ),
  time = c(Op, Ps, top10tdc, top10Con, BL210min, BL21tdc0min, BL2120min,BL21tdc20min),
  oct_conc = c(
    rep("30% On food", length(c(Op, Ps, top10tdc, top10Con))),
    rep("100% 0 min off food", length(c(BL210min, BL21tdc0min))),
    rep("100% 20 min off food", length(c(BL21tdc20min, BL2120min)))
  )
)


data %>% droplevels() %>% ggplot(data = ., aes(x =food, y = time)) +
  geom_jitter(width = 0.2) +
  ProvidenciaChemo::add.mean(time) + theme_classic() +
  facet_grid(.~oct_conc, drop = TRUE, scales = "free")
