############ exponential fitting genotypes ############
#setWD to matlab files#
path<-getwd()
files<-list.files(path, pattern = "reg.mat")

files<-list.files(path, pattern = "N2")
#WT
WT_corrected<-data.frame(sapply(files, FUN=exp.fit.all)) %>% mutate(time = (1:nrow(.))/4) %>%
  melt(id="time", variable.name = "animal", value.name = "signal") %>% mutate(genotype = "WT")

#plot avg
WT_corrected %<>% group_by(animal) %>%
  mutate(lag.sig = lag(signal, n=9),bins = ntile(time,n=13))


setwd("/Volumes/sengupta-lab/sengupta-lab/mikeod/OlfactoryScope/6_13_17/061317_glycerol_ASH_osm5/sum/")
path<-getwd()
files<-list.files(path, pattern = "ASH")

#mutant
files<-list.files(path, pattern = "mut")
mut_corrected<-data.frame(sapply(files, FUN=exp.fit.all)) %>% mutate(time = (1:nrow(.))/4) %>%
  melt(id="time", variable.name = "animal", value.name = "signal") %>% mutate(genotype = "ft7")

#plot avg
mut_corrected %<>% group_by(animal) %>%
  mutate(lag.sig = lag(signal, n=9),bins = ntile(time,n=13))

library(viridis)
rbind(WT_corrected,mut_corrected) %>% group_by(genotype,time) %>%
  summarise(deltaF = mean(signal), se = sd(signal)/sqrt(n())) %>%
  #mutate(mean.delF = mean(lag.sig), sd.delF = sd(lag.sig)) %>%
  ggplot(aes(x=time, y = deltaF)) +
  geom_ribbon(aes(ymin = deltaF - se, ymax = deltaF + se, fill = genotype), alpha = 0.5) +
  geom_line(aes(colour = genotype)) +
  viridis::scale_fill_viridis(discrete = TRUE) +
  scale_colour_manual(values = c("black", "black")) +
  geom_segment(aes(x = 30, y = 1.5, xend = 60, yend = 1.5)) +
  theme_classic()


### ASH soma in OP vs Ps octanol
#OP
# setwd("/Volumes/sengupta-lab/sengupta-lab/mikeod/OlfactoryScope/6_13_17/061317_octanol_ASH_N2OP/sum_ASH_soma/")
# path<-getwd()
# files<-list.files(path, pattern = "OP")
# Soma_corrected_OP<-data.frame(sapply(files, FUN=exp.fit.all.log.lin, skip.time = 10)) %>%
#   mutate(time = (1:nrow(.))/4) %>%
#   melt(id="time", variable.name = "animal", value.name = "signal") %>%
#   mutate(neuron = "ASH", food = "OP", pos = "soma")
setwd("/Volumes/sengupta-lab/sengupta-lab/mikeod/OlfactoryScope/6_13_17/061317_octanol_ASH_N2OP/sum_ASH_process/")
path<-getwd()
files<-list.files(path, pattern = "ASH")
axon_corrected_OP<-data.frame(sapply(files, FUN=exp.fit.all.log.lin, skip.time = 10)) %>%
  mutate(time = (1:nrow(.))/4) %>%
  melt(id="time", variable.name = "animal", value.name = "signal") %>%
  mutate(neuron = "ASH", food = "OP", pos = "axon")
#Ps
# setwd("/Volumes/sengupta-lab/sengupta-lab/mikeod/OlfactoryScope/6_13_17/061317_octanol_ASH_N2Ps/sum_ASH_soma/")
# path<-getwd()
# files<-list.files(path, pattern = "Ps")
# Soma_corrected_Ps<-data.frame(sapply(files, FUN=exp.fit.all.log.lin, skip.time = 10)) %>% mutate(time = (1:nrow(.))/4) %>%
#   melt(id="time", variable.name = "animal", value.name = "signal") %>% mutate(neuron = "ASH", food = "Ps", pos = "soma")
setwd("/Volumes/sengupta-lab/sengupta-lab/mikeod/OlfactoryScope/6_13_17/061317_octanol_ASH_N2Ps/sum_ASH_process/")
path<-getwd()
files<-list.files(path, pattern = "ASH")
axon_corrected_Ps<-data.frame(sapply(files, FUN=exp.fit.all.log.lin, skip.time = 10)) %>%
  mutate(time = (1:nrow(.))/4) %>%
  melt(id="time", variable.name = "animal", value.name = "signal") %>%
  mutate(neuron = "ASH", food = "Ps", pos = "axon")




rbind(axon_corrected_OP,axon_corrected_Ps) %>%
#rbind(Soma_corrected_OP,Soma_corrected_Ps) %>%
group_by(food,time,pos) %>%
  summarise(deltaF = mean(signal), se = sd(signal)/sqrt(n())) %>%
  #mutate(mean.delF = mean(lag.sig), sd.delF = sd(lag.sig)) %>%
  ggplot(aes(x=time, y = deltaF)) +
  geom_ribbon(aes(ymin = deltaF - se, ymax = deltaF + se, fill = food), alpha = 0.2) +
  geom_line(aes(colour = food)) +
  scale_colour_manual(values = c("grey", "blue")) +
  scale_fill_manual(values = c("grey", "blue")) +
  geom_segment(aes(x = 30, y = max(deltaF + 0.05), xend = 60, yend = max(deltaF + 0.05)), colour = "white") +
  facet_grid(pos~.) +
  geom_segment(aes(x = 30, y = max(deltaF + 0.2), xend = 60, yend = max(deltaF + 0.2))) +
  coord_cartesian(xlim=c(20,80)) +
  theme_classic() + theme_black()

#by animal:
rbind(Soma_corrected_OP,axon_corrected_OP) %>%
  ggplot(aes(x=time, y = signal)) +
  geom_line(aes(colour = pos, lty=pos)) +
  viridis::scale_fill_viridis(discrete = TRUE) +
  geom_segment(aes(x = 30, y = 0.4, xend = 60, yend = 0.4)) +
  theme_classic() + facet_wrap(~animal)

rbind(Soma_corrected_Ps,axon_corrected_Ps) %>%
  ggplot(aes(x=time, y = signal)) +
  geom_line(aes(colour = pos, lty=pos)) +
  viridis::scale_fill_viridis(discrete = TRUE) +
  geom_segment(aes(x = 30, y = 0.4, xend = 60, yend = 0.4)) +
  theme_classic() + facet_wrap(~animal)


# ASI ---------------------------------------------------------------------
##### ASI soma glycerol OP vs PS
#OP
#setwd("/Volumes/sengupta-lab/sengupta-lab/mikeod/OlfactoryScope/6_13_17/061317_glycerol_ASH_N2OP/sum_ASI_soma/")
path<-getwd()
files<-list.files(path, pattern = "OP")
ASI_glycerol_OP<-data.frame(sapply(files, FUN=exp.fit.all.log.lin, skip.time = 10)) %>%
  mutate(time = (1:nrow(.))/4) %>%
  melt(id="time", variable.name = "animal", value.name = "signal") %>%
  mutate(neuron = "ASH", food = "OP", pos = "soma")


#PS
#setwd("/Volumes/sengupta-lab/sengupta-lab/mikeod/OlfactoryScope/6_13_17/061317_glycerol_ASH_N2Ps/sum_ASI_soma/")
#path<-getwd()
files<-list.files(path, pattern = "Ps")
ASI_glycerol_Ps<-data.frame(sapply(files, FUN=exp.fit.all.log.lin, skip.time = 10)) %>%
  mutate(time = (1:nrow(.))/4) %>%
  melt(id="time", variable.name = "animal", value.name = "signal") %>%
  mutate(neuron = "ASH", food = "Ps", pos = "soma")

rbind(ASI_glycerol_OP, ASI_glycerol_Ps) %>% group_by(food,time,pos) %>%
  summarise(deltaF = mean(signal), se = sd(signal)/sqrt(n())) %>%
  #mutate(mean.delF = mean(lag.sig), sd.delF = sd(lag.sig)) %>%
  ggplot(aes(x=time, y = deltaF)) +
  geom_ribbon(aes(ymin = deltaF - se, ymax = deltaF + se, fill = food), alpha = 0.5) +
  geom_line(aes(colour = pos, lty=food)) +
  viridis::scale_fill_viridis(discrete = TRUE) +
  scale_colour_manual(values = c("black", "black")) +
  facet_grid(pos~.) +
  geom_segment(aes(x = 30, y = max(deltaF + 0.), xend = 60, yend = max(deltaF + 0.2))) +
  coord_cartesian(xlim=c(20,80)) +
  theme_classic()


# ASH soma glycerol OP vs PS ----------------------------------------------

#OP
setwd("/Volumes/sengupta-lab/sengupta-lab/mikeod/OlfactoryScope/6_13_17/061317_glycerol_ASH_N2OP/sum_ASH_soma/")
path<-getwd()
files<-list.files(path, pattern = "ASH")
Soma_glycerol_OP<-data.frame(sapply(files, FUN=exp.fit.all.log.lin, skip.time = 10)) %>%
  mutate(time = (1:nrow(.))/4) %>%
  melt(id="time", variable.name = "animal", value.name = "signal") %>%
  mutate(neuron = "ASH", food = "OP", pos = "soma")

#PS
setwd("/Volumes/sengupta-lab/sengupta-lab/mikeod/OlfactoryScope/6_13_17/061317_glycerol_ASH_N2Ps/sum_ASH_soma/")
path<-getwd()
files<-list.files(path, pattern = "ASH")
Soma_glycerol_Ps<-data.frame(sapply(files, FUN=exp.fit.all.log.lin, skip.time = 10)) %>%
  mutate(time = (1:nrow(.))/4) %>%
  melt(id="time", variable.name = "animal", value.name = "signal") %>%
  mutate(neuron = "ASH", food = "Ps", pos = "soma")

rbind(Soma_glycerol_OP, Soma_glycerol_Ps) %>% group_by(food,time,pos) %>%
  summarise(deltaF = mean(signal), se = sd(signal)/sqrt(n())) %>%
  #mutate(mean.delF = mean(lag.sig), sd.delF = sd(lag.sig)) %>%
  ggplot(aes(x=time, y = deltaF)) +
  geom_ribbon(aes(ymin = deltaF - se, ymax = deltaF + se, fill = food), alpha = 0.2) +
  geom_line(aes(colour = food)) +
  scale_colour_manual(values = c("grey", "blue")) +
  scale_fill_manual(values = c("grey", "blue")) +
  facet_grid(pos~.) +
  geom_segment(aes(x = 30, y = max(deltaF + 0.2), xend = 60, yend = max(deltaF + 0.2)), colour = "white") +
  coord_cartesian(xlim=c(20,80)) +
  theme_classic() + theme_black()

