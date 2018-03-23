#' exp.fit.all.log
#'
#' Function takes a list of calcium imaging files by filenames in the current working directory (.mat files)
#' and performs an exponential fit to normalize the deltaF/F signal. _This function does not incorporates a linear
#' term to account for gradual rise in signal, which can alter the sign of the single exponential._
#' Will output plots showing orignal and fitted values, as well as corrected values. If an object
#' is assigned, this will be a vector of corrected values
#' @param filename filepaths of .mat files which have a "signal"  and "time" field, or if files are in current
#' working directory, then just list of filenames.
#' @param skip.time number of seconds to skip at the beginning for the exponential fit. N=10 improves the fit.
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' @export
#' @examples data <- exp.fit.all.log(files[1], skip.time = 10)

exp.fit.all.log<-function (filename,skip.time) {
  matfile<-R.matlab::readMat(filename, fixNames = TRUE)
  signal<-matfile$signal
  time<-(1:length(signal))/4
  df<-data.frame(time, signal)
  rm(signal)
  rm(time)

  #fit to first N(skip.time) seconds to 30 sec and
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
