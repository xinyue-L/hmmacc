#' Density for Zero-Inflated Gaussian Distribution
#'
#' \code{dnorm0} calculates the density for zero-inflated Gaussian distribution.
#'
#' @param x the vector of observations.
#' @param alpha the inflated zero component.
#' @param mu the mean of the Gaussian part.
#' @param sigma the standard deviation of the Gaussian part.
#' @param log calculates the density on the log scale or not. The default is False.
#'
#' @keywords Zero-Inflated Gaussian
#'
#' @return the vector of calculated densities for vector \code{x}; of the same length as \code{x}.
#'
#' @examples
#' dnorm0(x=c(1,3,5,7),alpha=0.3,mu=5,sigma=1,log=F)
#'
#' @export

dnorm0 <- Vectorize(function(y,alpha,mu,sd1,log=F) {
  if(y==0) {
    return(ifelse(log==F,
                  alpha+(1-alpha)*pnorm(y,mean=mu,sd=sd1),
                  log(alpha+(1-alpha)*pnorm(y,mean=mu,sd=sd1))))
  } else {
    return(ifelse(log==F,
                  (1-alpha)*dnorm(y,mean=mu,sd=sd1),
                  log(1-alpha)+dnorm(y,mean=mu,sd=sd1,log=TRUE)))
  }
},vec="y")
