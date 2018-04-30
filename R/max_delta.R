#' max_delta
#'
#' Helper function for reordering levels of calcium response based on the maximum change in reponse
#' aligned to the end of a pulse. Mainly for use in plotGCaMP_multi.R. End value is derived from input
#' in plotGCaMP_multi.R
#'
#' @export
#' @examples ggplot(aes(x = time, y = fct_reorder(animal_num, delF, max_delta)))

max_delta <- function(x, end = 60) {
  if(max(x[(end*4 - 10):(end*4 + 20)]) > 0) {
      max(x[(end*4 - 10):(end*4 + 20)]) -
        min(x[(end*4 - 10):(end*4 + 20)])
  } else {
      max(x[(end*4 - 10):(end*4 + 20)]) +
        min(x[(end*4 -10):(end*4 + 20)])
  }
}
