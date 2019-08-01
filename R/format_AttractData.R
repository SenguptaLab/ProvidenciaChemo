#' format_AttractData
#'
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' @export
#' @examples data %>% format_AttractData()

format_AttractData <- function(data) {
  data %>% mutate(
         strain = fct_relevel(strain, 'OP50'),
         nAll = nCue + nControl,
         plate = seq(1:nrow(.)),
         #CI = (nCue - nControl) / (nAll + outside),
         CI = (nCue - nControl) / (nAll),
         p = nCue  / nAll,
         logit.p = case_when(
            p == 0 ~ boot::logit(0.005),
            p == 1 ~ boot::logit(0.995),
            TRUE ~ boot::logit(p))) %>%
    dplyr::filter(nAll > 10) %>%
    mutate(rel.Logit = logit.p - mean(dplyr::filter(., strain == 'OP50')$logit.p),
           plate = interaction(date,assay,plate))
}
