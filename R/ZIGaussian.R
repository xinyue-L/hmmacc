#' Fitting Zero-Inflated Gaussian Distribution
#'
#' \code{ZIGaussian} is used to fit zero-inflated Gaussian distribution.
#'
#' @param y the vector of observations for fitting zero-inflated Gaussian.
#' @param weights an optional vector of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector.
#' @param maxstep the maximal number of steps for fitting the Gaussian part and the zero part. The default is 100.
#' @param reltol the relative tolerance of log likelihood in the fitting procedure. The default is 10^-6.
#'
#' @keywords Zero-Inflated Gaussian
#'
#' @return a list of two elements.
#' @return \item{count}{a list of two elements for the Gaussian part: mean and sigma.}
#' @return \item{zero}{the estimated proportion of zeros for the zero inflation part.}
#'
#' @examples
#' data(ob)
#' ZIGaussian(ob,weights=NULL,maxstep=100,reltol=10^-6)
#'
#' @export

ZIGaussian <- function(y,weights=NULL,maxstep=100,reltol=10^-6) {
  invlogit <- function(x) exp(x)/(exp(x)+1)
  ilogit <- function(p) log(p/(1-p))

  ziNorm <- function(parms,wt) {
    mu <- as.vector(X * parms[1:kx] + offsetx)
    phi <- as.vector(invlogit(Z * parms[(kx + 1):(kx + kz)] +
                                offsetz))
    loglik0 <- log(phi + (1-phi)*dnorm(0,mu,sig)/pnorm(0,mu,sig,lower.tail=F))
    loglik1 <- log(1 - phi) - log(pnorm(0,mu,sig,lower.tail=F)) + dnorm(y, mean=mu, sd=sig, log = TRUE)
    loglik <- sum(wt[y0] * loglik0[y0]) + sum(wt[y1] *
                                                loglik1[y1])
    loglik
  }

  if(is.null(weights)) {
    weights.o <- rep(1,length(y))
  } else {
    weights.o <- weights
  }

  y0 <- y <= 0
  y1 <- y > 0
  kx <- kz <- 1
  X <- Z <- rep(1,length(y))
  offsetx <- offsetz <- 0
  sig <- 1
  model_count <- Gaussian(y,weights=weights.o) #fit normal as a proxy to truncated normal
  model_zero <- suppressWarnings(glm(as.integer(y0)~1,family=binomial,weights=weights.o))

  start <- list(count = model_count$coefficients, zero = model_zero$coefficients)

  mui <- rep(model_count$coefficients,length(y))
  sig <- sqrt(model_count$dispersion)
  probi <- model_zero$fitted
  probi <- probi/(probi + (1 - probi) * pnorm(0, mui,sd=sig))
  probi[y1] <- 0
  ll_new <- ziNorm(c(start$count, start$zero),wt=weights.o)
  ll_old <- 2 * ll_new
  step <- 0
  while (abs((ll_old - ll_new)/ll_old) > reltol & step < maxstep) {
    step <- step + 1
    ll_old <- ll_new

    mui <- rep(sum(y*(weights.o*(1-probi)))/sum((weights.o*(1-probi))),length(y))
    sig <- sqrt(sum((weights.o*(1-probi))*((y-mui[1])^2))/(sum(weights.o*(1-probi))-1))

    model_zero <- suppressWarnings(glm(probi~1,family=binomial,weights=weights.o,
                                       start=start$zero))

    probi <- model_zero$fitted
    probi <- probi/(probi + (1 - probi) * pnorm(0, mui,sd=sig))
    probi[y1] <- 0
    start <- list(count = mui[1],
                  sigma = sig,
                  zero = model_zero$coefficients)
    ll_new <- ziNorm(c(start$count, start$zero),wt=weights.o)
    #print(paste("step = ",step,"; loglik = ",ll_new))
  }
  if(step>=maxstep) print(paste("max step = ",maxstep, " reached"))
  return(list(count=list(mean=start$count,sigma=start$sigma),
              zero=invlogit(start$zero)))
}
