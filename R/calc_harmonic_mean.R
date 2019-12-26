#' Calculate the Harmonic Mean
#'
#' @inheritParams psych::harmonic.mean
#'
#' @return The harmonic mean(s)
#' @export
#'
#' @examples
#'
#' x <- 1:5
#' calc_harmonic_mean(x)
calc_harmonic_mean <- function(x, na.rm = TRUE, zero = TRUE) {
  psych::harmonic.mean(x = x, na.rm = na.rm, zero = zero)
}
