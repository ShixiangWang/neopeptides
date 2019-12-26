#' Calculate the Geometric Mean
#'
#' @inheritParams psych::geometric.mean
#'
#' @return The geometric mean(s)
#' @export
#'
#' @examples
#'
#' x <- 1:5
#' calc_geometric_mean(x)
calc_geometric_mean <- function(x, na.rm = TRUE) {
  psych::geometric.mean(x = x, na.rm = na.rm)
}
