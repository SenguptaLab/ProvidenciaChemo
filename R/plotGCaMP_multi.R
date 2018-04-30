#' plotGCaMP_multi
#'
#' Function is a wrapper for exp.fit.all.log.lin which outputs corrected GCaMP signals as well as plots showing
#' original and corrected signals. Also plots average trace for GCaMP signal. Inputs are matfiles.
#' The script searches recursively for matfiles from a starter file (can be any placeholder file). Need to exclude
#' matfiles that are not GCaMP files either using FileFilter, or by putting only relevant files in the folder.
#' Requires max_delta helper function
#' time etc...
#' @param FileFilter string to search/subset filenames
#' @param genotype label the genotype for these data
#' @param cue label the stimulus cue.
#' @param food label to food cue
#' @param start begin of stimulus
#' @param end end time of stimulus
#' @param center_on_pulse optional parameter to center delF values by the mean of the stimulus duration
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' @export
#' @examples data <- plotGCaMP_multi(N2, genotype = N2, cue = octanol)
#'
plotGCaMP_multi <- function(FileFilter,
                            genotype,
                            cue = cue,
                            food = OP50,
                            start = 30,
                            end = 60,
                            center_on_pulse = FALSE,
                            ...) {
  library(tidyverse)
  library(magrittr)
  library(patchwork)
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

  if(center_on_pulse %in% c(1,2)) {
    means <- data %>% unnest %>%
      group_by(animal) %>%
      filter(time > start & time < end) %>%
      summarise(mean_pulse_delF = mean(delF))
    data <- full_join(data,means)
    data %<>% unnest %>%
      group_by(animal) %>%
      mutate(delF = delF - mean_pulse_delF) %>%
      nest()
  }

  plot1 <- data %>% unnest() %>%
    ggplot(aes(x = time, y = delF)) +
    geom_line(aes(group = animal), alpha = 0.2) +
    geom_smooth(method = "loess", span = 0.05) +
    theme_classic() +
    geom_segment(aes(x = !!start,
                     y = max(delF) + 0.1*max(delF),
                     xend = !!end,
                     yend =max(delF) + 0.1*max(delF))) +
    annotate(geom = "text",
             label = cue,
             y = max(unnest(data)$delF) + 0.2*max(unnest(data)$delF),
             x = (start + end)/2,
             size = 10) +
    annotate(geom = "text",
             label = genotype,
             y = max(unnest(data)$delF) + 0.2*max(unnest(data)$delF),
             x = 10,
             fontface = "italic",
             size = 10) +
    theme(axis.text = element_text(size = 14))

if(center_on_pulse == 1) {
  plot2 <- data %>%
    unnest() %>%
    group_by(animal_num) %>%
    ggplot(aes(x = time, y = fct_reorder(factor(animal_num), delF, max_delta))) +
    geom_tile(aes(fill = delF)) +
    scale_fill_viridis_c(option = "magma",
                         breaks = c(-1,0,1),
                         labels = c("-1", "0", "1"),
                         limits = c(-1,1),
                         oob =squish) +
    theme_classic() +
    theme(axis.text = element_text(size = 14))
} else {
  if(center_on_pulse == 2) {
  plot2 <- data %>%
    unnest() %>%
    group_by(animal_num) %>%
    ggplot(aes(x = time, y = fct_reorder(factor(animal_num), delF, max_negdelta))) +
    geom_tile(aes(fill = delF)) +
    scale_fill_viridis_c(option = "magma",
                         breaks = c(-1,0,1),
                         labels = c("-1", "0", "1"),
                         limits = c(-1,1),
                         oob =squish) +
    theme_classic() +
    theme(axis.text = element_text(size = 14))
  } else {
  plot2 <- data %>%
    unnest() %>%
    group_by(animal_num) %>%
    ggplot(aes(x = time, y = fct_reorder(factor(animal_num), delF, max_delta))) +
    geom_tile(aes(fill = delF)) +
    scale_fill_viridis_c(option = "magma",
                         breaks = c(-.5,0,1.5),
                         labels = c("-.5", "0.0", "1.5+"),
                         limits = c(-.5,1.5),
                         oob =squish) +
    theme_classic() +
    theme(axis.text = element_text(size = 14))
  }
}

  plots <- plot1 + plot2 + plot_layout(ncol = 1, heights = c(2,1))

  ggsave(plots, filename = file.path(folderPath,paste0(genotype,"_",cue,"plots.png")),
         width = 11, height = 8.5, units = "in")

  return(list(data = data, plot = plots))
}
