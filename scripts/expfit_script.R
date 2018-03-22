library(R.matlab)
library(tidyr)
library(reshape2)
library(dplyr)

exp.fit<-function (data) {
  signal<-data$signal
  time<-(1:length(signal))/4
  df<-data.frame(time, signal)
  rm(signal)
  rm(time)

  #fit to first 4 to 30 sec
  fit1<-lm(data=df[c(10:120,1440:1560),],signal~log(time))

  #fit all data
  fit2<-lm(data=df[-(1:10),],signal~log(time))

  df$fitted1<-predict(fit1,newdata=df)
  df$fitted2<-predict(fit2,newdata=df)
  # correct after fitted values go below zero (~ 20s) could do this in function
  df$corrected<-df$signal
  df$corrected[80:nrow(df)] <- (df$signal[80:nrow(df)] - df$fitted2[80:nrow(df)])
  df$corrected1<-df$signal
  df$corrected1[80:nrow(df)] <- (df$signal[80:nrow(df)] - df$fitted1[80:nrow(df)])

  #plot fits to inspect
  p<- ggplot(df, aes(x=time)) + geom_line(aes(y=signal), colour="black") +
    geom_line(aes(y=fitted1), colour = "red", linetype="dashed") +
    geom_line(aes(y=fitted2), colour = "blue", linetype="dashed") +
    geom_line(aes(y=corrected), colour = "blue") +
    geom_line(aes(y=corrected1), colour = "red") +
    theme_classic()
  print(p)
  return <- df$corrected
}

exp.fit.all.log.lin<-function (filename, skip.time) {
  matfile<-readMat(filename, fixNames = TRUE)
  signal<-matfile$signal
  time<-(1:length(signal))/4
  df<-data.frame(time, signal)
  rm(signal)
  rm(time)

  #fit to first N(skip.time) seconds to 30 sec
  fit1<-lm(data=df[c(skip.time:120,300:360),],signal~log(time)+time) # plus last 15s

  fitted<-predict(fit1,newdata=df)

  df$fitted <- fitted
  # }

  # correct after fitted values go below zero (~ 20s) could do this in function
  df$corrected<-df$signal
  df$corrected[skip.time:nrow(df)] <- (df$signal[skip.time:nrow(df)] - df$fitted[skip.time:nrow(df)])

  #plot fits to inspect
  p<- ggplot(df, aes(x=time)) + geom_line(aes(y=signal), colour="black") +
    geom_line(aes(y=fitted), colour = "red", linetype="dashed") +
    geom_line(aes(y=corrected), colour = "blue") +
    theme_classic()
  print(p)
  return <- df$corrected
}

exp.fit.all.log<-function (filename,skip.time) {
  matfile<-readMat(filename, fixNames = TRUE)
  signal<-matfile$signal
  time<-(1:length(signal))/4
  df<-data.frame(time, signal)
  rm(signal)
  rm(time)

  #fit to first 4 to 30 sec and
  fit1<-lm(data=df[c(skip.time:120,300:360),],signal~log(time)) # plus last 15s
  #fit2<-lm(data=df[-(1:10),],signal~log(time))
  fitted<-predict(fit1,newdata=df)

  # if(fitted[1] < fitted[length(fitted)]) {
  #   df$fitted <- 0 # don't correct if signal increases over time
  # } else {
  # fitted[fitted > 0.001] <- 0
  df$fitted <- fitted
  # }

  # correct after fitted values go below zero (~ 20s) could do this in function
  df$corrected<-df$signal
  df$corrected[skip.time:nrow(df)] <- (df$signal[skip.time:nrow(df)] - df$fitted[skip.time:nrow(df)])

  #plot fits to inspect
  p<- ggplot(df, aes(x=time)) + geom_line(aes(y=signal), colour="black") +
    geom_line(aes(y=fitted), colour = "red", linetype="dashed") +
    geom_line(aes(y=corrected), colour = "blue") +
    theme_classic()
  print(p)
  return <- df$corrected
}

exp.fit.all<-function (filename) {
  matfile<-readMat(filename, fixNames = TRUE)
  signal<-matfile$signal
  time<-(1:length(signal))/4
  df<-data.frame(time, signal)
  rm(signal)
  rm(time)

  #fit to first 4 to 30 sec (frame rate = 4/sec)
  fit1<-lm(data=df[c(10:120,1500:1560),],signal~log(time+40)) # plus last 15s
  #fit1<-lm(data=df[c(10:120),],signal~log(time+40))
  fitted<-predict(fit1,newdata=df)
  fitted[fitted > 0.05] <- 0
  df$fitted<-fitted
  # correct after fitted values go below zero (~ 20s) could do this in function
  df$corrected<-df$signal
  df$corrected <- (df$signal - df$fitted)

  #plot fits to inspect
  p<- ggplot(df, aes(x=time)) + geom_line(aes(y=signal), colour="black") +
    geom_line(aes(y=fitted), colour = "red", linetype="dashed") +
    geom_line(aes(y=corrected), colour = "blue") +
    theme_my
  print(p)
  return <- df$corrected
}


library(RcppRoll)
library(zoo)

# plot by animal, bin
rbind(WTcorrected, tubcorrected) %>% group_by(animal) %>%
  mutate(lag.sig = lag(signal, n=9),bins = ntile(time,n=13)) %>%
  ggplot(aes(x=time, y=lag.sig, colour = animal)) + geom_line() + facet_grid(.~genotype)

#####plot by delta for each curve - need function for this
rbind(WTcorrected, tubcorrected) %>% group_by(animal) %>%
  mutate(lag.sig = lag(signal, n=9),bins = ntile(time,n=13)) %>% filter(bins >2) %>% filter(bins %% 2 == 1) %>%
  group_by(animal, bins) %>% mutate(delta = max(lag.sig) - min(lag.sig)) %>%
  ggplot(aes(x=bins, y=delta, colour = animal)) + geom_point() + facet_grid(.~genotype) + geom_smooth(method = "lm", formula = y ~ I(log(x)))

rbind(WTcorrected, tubcorrected) %>% filter(animal != "X6pulseWT_11_reg.mat") %>%
  ggplot(aes(x=time, y=signal, colour = genotype)) + stat_summary(geom="line", fun.y=mean)

p <- rbind(WTcorrected, tubcorrected) %>%  mutate(smooth = rollmean(signal, 10, na.pad = TRUE)) %>% filter(!animal %in% c("tub130_1.mat", "X6pulseWT_11_reg.mat")) %>%
  ggplot(aes(x=time, y=smooth, colour = genotype)) + stat_summary(geom="line", fun.y=mean)
p

p <- rbind(WTcorrected, tubcorrected) %>%  mutate(smooth = rollmean(signal, 10, na.pad = TRUE)) %>% ggplot(aes(x=time, y=smooth, colour = genotype)) + stat_summary(geom="line", fun.y=mean)


ggplot(aes(x=time, y=signal, colour = "red")) + stat_summary(geom="line", fun.y=mean)

ggplot(WTcorrected_long, aes(x=time, y=signal)) + stat_summary(geom="line", fun.y=mean)
ggplot(tubcorrected_long, aes(x=time, y=signal)) + stat_summary(geom="line", fun.y=mean)


ggplot(WTcorrected_long, aes(x=time, y=signal)) + stat_summary(geom="line", fun.y=median, colour = "red") +
  stat_summary(aes(x=time, y=signal), data=tubcorrected_long, fun.y=median, colour="blue")


#WTcorrected$mean<-mean(WTcorrected[1:nrow(WTcorrected),])
#WTcorrected$mean<-lapply(WTcorrected[,1:length(WTfiles)], FUN=mean)

