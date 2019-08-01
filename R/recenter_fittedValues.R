#' recenter_fittedValues
#'
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' @export
#' @examples data %>% recenter_fittedValues

recenter_fittedValues <- function(data,
                                  model,
                                  day.correct = FALSE,
                                  transform = "logit",
                                  BayesFit = "fitted_draws",
                                  ...) {
  day.correct <- quo_name(enquo(day.correct))

  mod_predictors <- ref_grid(model) %>%
    broom::tidy() %>%
    select(-matches("estimate|lower.HPD|upper.HPD")) %>% colnames()

  if(BayesFit == "fitted_draws") {
    #this will not marginalize over group effects
    data %>% group_by_at(vars(mod_predictors)) %>%
    data_grid(!!!syms(mod_predictors)) %>%
    add_fitted_draws(model, re_formula = NA) -> data
    if(transform == "logit") {
      data %>%
        mutate(logit.p = boot::logit(.value)) %>%
        ungroup() -> data
    } else {

    }
    comment(data) <- "fitted_draws"
  }

  if(BayesFit == "HDI") {
    # only for comparing all to a common control
    #sjstats::hdi(stan_glmm, type = "fixed", prob = c(0.66, .95)) %>% tibble()
    mod.data <- emmeans::ref_grid(model) %>%
      emmeans::contrast(method = "trt.vs.ctrl") %>%
      coda::as.mcmc() %>%
      summary(quantiles = c(0.025, 0.17, 0.5, 0.83, 0.975))
    mod.data <- mod.data[[2]] %>%
      data.frame() %>%
      rownames_to_column("contrast")
    if(length(mod_predictors > 1)) {
      mod_predictors <- rev(mod_predictors)
    }
    print(mod_predictors)
    predictors <- data %>% modelr::data_grid(!!!syms(mod_predictors)) %>%
      slice(2:n())
    data <- cbind(mod.data, predictors) %>%
      mutate(interval_type = "HDI",
             data_type = "fit") %>%
      #dplyr::select(1:10) %>%
      rename(lower.2.5 = 2,
             lower.17 = 3,
             median = 4,
             upper.83 = 5,
             upper.97.5 = 6)
  }


  if(BayesFit == "fitted_draws" & transform == "logit") {
    if(day.correct == "OP50_N2") {
    means <- data %>%
      filter(strain == "OP50", genotype == "N2") %>%
      group_by(date) %>%
      summarise(meanOP50 = mean(logit.p))
    }

    if(day.correct == "treatment") {
    data <- mutate(data, genotype = "N2")
    means <- data %>%
      ungroup() %>%
      filter(strain == "OP50", treatment %in% c("control", "none")) %>%
      mutate(genotype = "N2") %>%
      group_by(genotype) %>%
      summarise(meanOP50 = mean(boot::logit(.value)))
    }

    if(day.correct == "OP50_by_genotype") {
      means <- data %>%
        ungroup() %>%
        filter(strain == "OP50") %>%
        group_by(genotype) %>%
        summarise(meanOP50 = mean(boot::logit(.value)))
    }

    if(day.correct == FALSE) {
    data <- mutate(data, genotype = "N2")
    means <- data %>%
      filter(strain == "OP50") %>%
      group_by(genotype) %>%
      summarise(meanOP50 = mean(logit.p))
    }
  }

  if(BayesFit == "fitted_draws" & transform == "logit") {
    data <- full_join(data, means) %>% mutate(rel.Logit = logit.p - meanOP50,
                                              interval_type = "fitted_draws")
  } else {
  return(data)
    }
}
