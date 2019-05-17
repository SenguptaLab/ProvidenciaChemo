#' recenter_fittedValues
#'
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' @export
#' @examples data %>% recenter_fittedValues

recenter_fittedValues <- function(data, model, day.correct = day.correct, transform = "logit", ...) {
  day.correct <- quo_name(enquo(day.correct))

  mod_predictors <- ref_grid(model) %>%
    broom::tidy() %>%
    select(1,2) %>% colnames()

  data %>% group_by_at(vars(mod_predictors)) %>%
    data_grid(!!!syms(mod_predictors)) %>%
    add_fitted_draws(model, re_formula = ~ 0) %>%
    mutate(logit.p = boot::logit(.value)) %>%
    ungroup() -> data

  # if(day.correct == "OP50") {
  #   means <- data %>%
  #     filter(strain == "OP50", genotype == "N2") %>%
  #     group_by(date) %>%
  #     summarise(meanOP50 = mean(logit.p))
  # }
  #
  # if(day.correct == "genotype") {
  #   means <- data %>%
  #     filter(strain == "OP50") %>%
  #     group_by(genotype, date) %>%
  #     summarise(meanOP50 = mean(logit.p))
  # }
  #
  if(day.correct == "treatment") {
    means <- data %>%
      ungroup() %>%
      #group_by_at(vars(mod_predictors[!mod_predictors %in% "treatment"])) %>%
      filter(strain == "OP50", treatment %in% c("control", "none")) %>%
      summarise(meanOP50 = mean(boot::logit(.value)))
  }

  # if(day.correct == "genotype+treatment") {
  #   means <- data %>%
  #     filter(strain == "OP50") %>%
  #     group_by(genotype, date, treatment) %>%
  #     summarise(meanOP50 = mean(logit.p))
  # }
  #
  #
  # if(day.correct == FALSE) {
  #   means <- data %>%
  #     filter(strain == "OP50") %>%
  #     group_by(genotype) %>%
  #     summarise(meanOP50 = mean(logit.p))
  # }


  #data <- full_join(data, means) %>% mutate(rel.Logit = logit.p - meanOP50) #[,c(1,2,4)])

  data <- data %>% mutate(meanOP50 = as.numeric(means), rel.Logit = logit.p - meanOP50)

  # if(day.correct) {
  #   data %>% mutate(rel.Logit = logit.p - meanOP50)
  # } else {
  #   data %>% mutate(rel.Logit = logit.p - mean(dplyr::filter(., genotype == "N2" & strain == 'OP50')$logit.p))
  # }

}
