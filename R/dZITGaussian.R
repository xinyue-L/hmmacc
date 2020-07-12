#' Density for Zero-Inflated Truncated Gaussian Distribution
#'
#' \code{dZITGaussian} calculates the density for zero-inflated truncated Gaussian distribution.
#'
#' @param x the vector of observations.
#' @param alpha the inflated zero component.
#' @param mu the mean of the truncated Gaussian part.
#' @param sigma the standard deviation of the truncated Gaussian part.
#' @param log calculates the density on the log scale or not. The default is False.
#'
#' @keywords Zero-Inflated Truncated Gaussian
#'
#' @return the vector of calculated densities for vector \code{x}; of the same length as \code{x}.
#'
#' @examples
#' dZITGaussian(x=c(1,3,5,7),alpha=0.3,mu=5,sigma=1,log=F)
#'
#' @export

dZITGaussian <- Vectorize(function(x,alpha,mu,sigma,log=F) {
  if(log==T) {
    if(x==0) {
      return(log(alpha+(1-alpha)/sigma*dnorm(0,mean=mu,sd=sigma)/pnorm(0,mean=mu,sd=sigma,lower.tail=FALSE)))
    } else {
      return(log(1-alpha)-log(sigma)+dnorm(x,mean=mu,sd=sigma,log=T)-pnorm(0,mean=mu,sd=sigma,lower.tail=FALSE,log.p=T))
    }
  } else {
    if(x==0) {
      return(alpha+(1-alpha)/sigma*dnorm(0,mean=mu,sd=sigma)/pnorm(0,mean=mu,sd=sigma,lower.tail=FALSE))
    } else {
      return((1-alpha)/sigma*dnorm(x,mean=mu,sd=sigma)/pnorm(0,mean=mu,sd=sigma,lower.tail=FALSE))
    }
  }
},vec="x")
