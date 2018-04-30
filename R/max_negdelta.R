#' max_negdelta
#'
#' Helper function for reordering levels of calcium response based on the maximum change in reponse
#' aligned to the end of a pulse, in downward responses. Mainly for use in plotGCaMP_multi.R. End value is derived from input
#' in plotGCaMP_multi.R
#'
#' @export
#' @examples ggplot(aes(x = time, y = fct_reorder(animal_num, delF, max_delta)))

max_negdelta <- function(x, end = 60) {
  if(min(x[(end*4 - 10):(end*4 + 20)]) < 0) {
    abs(min(x[(end*4 - 10):(end*4 + 20)]) -
      max(x[(end*4 - 10):(end*4 + 20)]))
  } else {
    abs(min(x[(end*4 - 10):(end*4 + 20)]) +
      max(x[(end*4 -10):(end*4 + 20)]))
  }
}
