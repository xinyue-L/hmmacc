#' auxillary function for fitting Gaussian distribution
#'
#' @param y the vector of observations for fitting zero-inflated truncated Gaussian.
#' @param weights an optional vector of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector.
#'
#' @keywords internal
#'
#' @return a list of two elements:
#' @return \item{coefficients}{mean.}
#' @return \item{dispersion}{standard deviation.}

Gaussian <- function(y,weights=NULL) {
  if(is.null(weights)) {
    return(list(coefficients=mean(y),dispersion=sd(y)))
  } else {
    m <- sum(y*weights)/sum(weights)
    s <- sqrt(sum(weights*((y-m)^2))/(sum(weights)-1))
    return(list(coefficients=m,
                dispersion=s))
  }
}
