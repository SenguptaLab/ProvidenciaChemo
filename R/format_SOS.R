#' format_AvoidData
#'
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' @export
#' @examples data %>% format_SOS()

format_SOS <- function(dataset, day_correct = genotype, reorder_levels = NA, ... ) {


  #print(glue::glue("Your are correcting all values relative to {day_correct}"))

 dataset <- dataset %>%
      mutate(response.time = case_when(
        response.time < 1 ~ 1,
        TRUE ~ response.time)) %>%
      group_by({{ day_correct }}, food, date) %>%
      summarise(meanOP = mean(log(response.time))) %>%
      filter(food == "OP50") %>%
      ungroup() %>%
      select({{ day_correct }}, date, meanOP) %>%
      full_join(., dataset) %>%
      mutate(rel_log = log(response.time) - meanOP,
             food = fct_relevel(food, "OP50"))


  if(!is.na(reorder_levels)) {
    column <- rlang::enquo(day_correct)
    column_name <- rlang::quo_name(order)

    dataset <- dataset %>%
      mutate(!!column_name :=
               fct_relevel({{ day_correct }}, {{ reorder_levels }}))
  }

  return(dataset)
}
