#' format_AvoidData_median
#'
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' @export
#' @examples data %>% format_AvoidData_median()

format_AvoidData_median <- function(data, day.correct = "OP50", center.data = FALSE, min_p = 0.005, ...) {
  data %>% mutate(
         strain = fct_relevel(strain, 'OP50'),
         plate = factor(seq(1:nrow(.))),
         nCue = RowA + RowB,
         nControl = RowE + RowF,
         nAll = nCue + nControl,
         CI = (nCue - nControl) / nAll,
         p = nCue  / nAll,
         logit.p = case_when(
            p == 0 ~ min(boot::logit(min_p), boot::logit(1/nAll)),
            #p == 0 ~ boot::logit(1/nAll),
            p == 1 ~ boot::logit(0.995),
            TRUE ~ boot::logit(p)),
         data_type = "raw") %>%
    mutate(plate = interaction(date,assay,plate)) -> data


if(day.correct == "OP50") {
  medians <- data %>%
    filter(strain == "OP50", genotype == "N2") %>%
    group_by(date) %>%
    summarise(medianOP50 = median(logit.p))
}

if(day.correct == "genotype") {
    medians <- data %>%
      filter(strain == "OP50") %>%
      group_by(genotype, date) %>%
      summarise(medianOP50 = median(logit.p))
}

  if(day.correct == "treatment") {
    medians <- data %>%
      filter(strain == "OP50", treatment %in% c("control", "none")) %>%
      group_by(date) %>%
      summarise(medianOP50 = median(logit.p))
  }

  if(day.correct == "treatment_overall") {
    medians <- data %>%
      filter(strain == "OP50", treatment %in% c("control", "none")) %>%
      group_by(genotype) %>%
      summarise(medianOP50 = median(logit.p))
  }

  if(day.correct == "genotype+treatment") {
    medians <- data %>%
      filter(strain == "OP50") %>%
      group_by(genotype, date, treatment) %>%
      summarise(medianOP50 = median(logit.p))
  }


  if(day.correct == FALSE) {
    medians <- data %>%
      filter(strain == "OP50") %>%
      group_by(genotype) %>%
      summarise(medianOP50 = median(logit.p))
  }


data <- full_join(data, medians) %>% mutate(rel.Logit = logit.p - medianOP50) #[,c(1,2,4)])


  # if(day.correct) {
  #   data %>% mutate(rel.Logit = logit.p - medianOP50)
  # } else {
  #   data %>% mutate(rel.Logit = logit.p - median(dplyr::filter(., genotype == "N2" & strain == 'OP50')$logit.p))
  # }

}
