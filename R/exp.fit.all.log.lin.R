#' exp.fit.all.log.lin
#'
#' Function takes a list of calcium imaging files by filenames in the current working directory (.mat files)
#' and performs an exponential fit to normalize the deltaF/F signal. This function incorporates a linear
#' to account for gradual rise in signal, which can alter the single exponential. Will output plots showing orignal
#' and fitted values, as well as corrected values. If an object is assigned, this will be a vector of corrected values
#' @param filename filepaths of .mat files which have a "signal"  and "time" field.
#' @param skip.time number of seconds to skip at the beginning for the exponential fit. N=10 improves the fit.
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' @export
#' @examples data <- exp.fit.all.log.lin(files[1], skip.time = 10)
#'
exp.fit.all.log.lin <- function(filename, skip.time) {
  matfile <- R.matlab::readMat(filename, fixNames = TRUE)
  signal <- matfile$signal
  time <- (1:length(signal)) / 4
  df <- data.frame(time, signal)
  rm(signal)
  rm(time)

  animal_name <- basename(filename)
    #quo_name(enquo(filename))

  # fit to first N(skip.time) seconds to 30 sec
  #fit1 <- lm(data = df[c(skip.time:120, 300:360), ], signal ~ log(time) + time) # plus last 15s
  fit1 <- lm(data = df[c(skip.time:120, 300:360), ], signal ~ log(time) + time) # plus last 15s
  if(fit1$coefficients[2] > 0) {
    fit1 <- lm(data = df[c(skip.time:120, 300:360), ], signal ~ log(time))
  }


  fitted <- predict(fit1, newdata = df)

  df$fitted <- fitted
  # }

  # correct after fitted values go below zero (~ 20s) could do this in function
  df %<>% mutate(corrected = case_when(
    fit1$coefficients[2] > 0 ~ signal, #ignore inverted exp fit
    fitted > 0 ~ signal, #ignor cases which have linear fit > 0
    TRUE ~ signal - fitted
  ))

  # plot fits to inspect
  p <- ggplot(df, aes(x = time, y = signal)) +
    geom_line(colour = "black") +
    geom_line(aes(y = fitted), colour = "red", linetype = "dashed") +
    geom_line(aes(y = corrected), colour = "blue") +
    annotate("text", label = animal_name, y = max(c(df$corrected,df$signal) + 0.1), x = 50) +
    theme_classic()
  print(p)
  delF <- df$corrected
  time <- df$time
  return <- data.frame(delF, time)
}
