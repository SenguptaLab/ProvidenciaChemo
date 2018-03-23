#' plotGCaMP_multi
#'
#' Function is a wrapper for exp.fit.all.log.lin which outputs corrected GCaMP signals as well as plots showing
#' original and corrected signals. Also plots average trace for GCaMP signal. Inputs are matfiles.
#' The script searches recursively for matfiles from a starter file (can be any placeholder file). Need to exclude
#' matfiles that are not GCaMP files either using FileFilter, or by putting only relevant files in the folder.
#' Plan is to update to include stim
#' time etc...
#' @param FileFilter string to search/subset filenames
#' @param genotype label the genotype for these data
#' @param cue label the stimulus cue.
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' @export
#' @examples data <- plotGCaMP_multi(N2, genotype = N2, cue = octanol)
#'
plotGCaMP_multi <- function(FileFilter, genotype, cue, food) {
  library(tidyverse)
  library(magrittr)
  FileFilter <- quo_name(enquo(FileFilter)) # make Filter usable inside other functions
  genotype <- quo_name(enquo(genotype))
  cue <- quo_name(enquo(cue))
  food <- quo_name(enquo(food))

  folderPath <- dirname(file.choose())
  files <- list.files(file.path(folderPath), pattern = "*.mat", recursive = TRUE)
  files <- files[stringr::str_detect(files, pattern = paste0(FileFilter))]
  filenames <- files
  files <- file.path(folderPath, files)
  #df <- data.frame(x = 1, genotype = genotype, cue = cue)

  data <- map(files, ~ exp.fit.all.log.lin(filename = ., skip.time = 10))
  data %<>% data_frame(data = .,
                      animal = filenames,
                      animal_num = seq(from = 1, to = length(filenames)),
                      genotype = genotype,
                      cue = cue,
                      food = food)
  plot <- data %>% unnest() %>%
    ggplot(aes(x = time, y = delF)) +
    geom_line(aes(colour = animal)) +
    guides(color = FALSE) +
    geom_smooth(method = "loess", span = 0.05) +
    theme_classic()
  return(list(data = data, plot = plot))
}
