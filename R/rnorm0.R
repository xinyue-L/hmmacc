#' Random Generation for the Zero-Inflated Gaussian Distribution (Approximation)
#'
#' \code{dnorm0} random generation for the zero-inflated Gaussian distribution (approximation).
#'
#' @param n the number of observations.
#' @param alpha the inflated zero component.
#' @param mu the mean of the truncated Gaussian part.
#' @param sigma the standard deviation of the truncated Gaussian part.
#' @param positive whether only positive values are allowed. The default is TRUE.
#'
#' @keywords Zero-Inflated Gaussian
#'
#' @return the vector of generated observations of length \code{n}.
#'
#' @examples
#' rnorm0(n=100,alpha=0.3,mu=5,sigma=1,positive=TRUE)
#'
#' @export

rnorm0 <- function(n,alpha,mu,sd1,positive) {
  binom <- rbinom(n=n,size=1,prob=1-alpha)
  normm <- rnorm(n=n,mean=mu,sd=sd1)
  if(positive==TRUE) {
    normm[normm<0] <- 0
  }
  return(binom*normm)
}
