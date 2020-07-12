#' Random Generation for the Zero-Inflated Truncated Gaussian Distribution
#'
#' \code{rZITGaussian} random generation for the zero-inflated truncated Gaussian distribution.
#'
#' @param n the number of observations.
#' @param alpha the inflated zero component.
#' @param mu the mean of the truncated Gaussian part.
#' @param sigma the standard deviation of the truncated Gaussian part.
#'
#' @keywords Zero-Inflated truncated Gaussian
#'
#' @return the vector of generated observations of length \code{n}.
#'
#' @examples
#' rZITGaussian(n=100,alpha=0.3,mu=5,sigma=1)
#'
#' @export

rZITGaussian <- function(n,alpha,mu,sigma) {
  binom <- rbinom(n=n,size=1,prob=1-alpha)
  normm <- rnorm(n=n,mean=mu,sd=sigma)
  while(any(normm <= 0)) {
    normm[which(normm < 0)] <- rnorm(n=sum(normm < 0),mean=mu,sd=sigma)
  }
  return(binom*normm)
}
