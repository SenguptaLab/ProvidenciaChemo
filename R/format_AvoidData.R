#' format_AvoidData
#'
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' @export
#' @examples data %>% format_AvoidData()

format_AvoidData <- function(data, day.correct = "OP50", center.data = FALSE, ...) {
  data %>% mutate(
         strain = fct_relevel(strain, 'OP50'),
         plate = seq(1:nrow(.)),
         nCue = RowA + RowB,
         nControl = RowE + RowF,
         nAll = nCue + nControl,
         CI = (nCue - nControl) / nAll,
         p = nCue  / nAll,
         logit.p = case_when(
            #p == 0 ~ boot::logit(0.005),
            p == 0 ~ boot::logit(1/nAll),
            p == 1 ~ boot::logit(0.995),
            TRUE ~ boot::logit(p))) %>%
    mutate(plate = interaction(date,assay,plate)) -> data


if(day.correct == "OP50") {
  means <- data %>%
    filter(strain == "OP50", genotype == "N2") %>%
    group_by(date) %>%
    summarise(meanOP50 = mean(logit.p))
}

if(day.correct == "genotype") {
    means <- data %>%
      filter(strain == "OP50") %>%
      group_by(genotype, date) %>%
      summarise(meanOP50 = mean(logit.p))
}

  if(day.correct == "treatment") {
    means <- data %>%
      filter(strain == "OP50", treatment == "live") %>%
      group_by(date) %>%
      summarise(meanOP50 = mean(logit.p))
  }

  if(day.correct == FALSE) {
    means <- data %>%
      filter(strain == "OP50") %>%
      group_by(genotype) %>%
      summarise(meanOP50 = mean(logit.p))
  }


data <- full_join(data, means) %>% mutate(rel.Logit = logit.p - meanOP50) #[,c(1,2,4)])


  # if(day.correct) {
  #   data %>% mutate(rel.Logit = logit.p - meanOP50)
  # } else {
  #   data %>% mutate(rel.Logit = logit.p - mean(dplyr::filter(., genotype == "N2" & strain == 'OP50')$logit.p))
  # }

}
