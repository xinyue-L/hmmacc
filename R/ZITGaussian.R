#' Fitting Zero-Inflated Truncated Gaussian Distribution
#'
#' \code{ZITGaussian} is used to fit zero-inflated truncated Gaussian distribution.
#' 
#' @import"truncreg"
#' 
#' @param y the vector of observations for fitting zero-inflated truncated Gaussian.
#' @param weights an optional vector of ‘prior weights’ to be used in the fitting process. Should be NULL or a numeric vector.
#' @param maxstep the maximal number of steps for fitting the truncated Gaussian part and the zero part. The default is 100.
#' @param reltol the relative tolerance of log likelihood in the fitting procedure. The default is 10^-6.
#'
#' @keywords Zero-Inflated Truncated Gaussian
#'
#' @return a list of two elements.
#' @return \item{count}{a list of two elements for the truncated Gaussian part: mean and sigma.}
#' @return \item{zero}{the estimated proportion of zeros for the zero inflation part.}
#'
#' @examples
#' ZITGaussian(y,weights=NULL,maxstep=100,reltol=10^-6)
#'
#' @export

ZITGaussian <- function(y,weights=NULL,maxstep=100,reltol=10^-6) {

  invlogit <- function(x) exp(x)/(exp(x)+1)
  ilogit <- function(p) log(p/(1-p))

  ziNorm <- function(parms,wt) {
    ## use global variables: X, Z, y, y0, y1!!!!!!!!!!
    mu <- as.vector(X * parms[1:kx] + offsetx)
    phi <- as.vector(invlogit(Z * parms[(kx + 1):(kx + kz)] +
                                offsetz))
    loglik0 <- log(phi + (1-phi)*dnorm(0,mu,sig)/pnorm(0,mu,sig,lower.tail=FALSE))
    loglik1 <- log(1 - phi) - log(pnorm(0,mu,sig,lower.tail=FALSE)) + dnorm(y, mean=mu, sd=sig, log = TRUE)
    loglik <- sum(wt[y0] * loglik0[y0]) + sum(wt[y1] * loglik1[y1])
    loglik
  }

  if(is.null(weights)) {
    w <- rep(1,length(y))
  } else {
    w <- weights
  }

  wok <- w!=0
  Xm <- rep(1,length(y))
  Xm <- Xm[wok]
  Ym <- y[wok]
  wm <- w[wok]
  Xmfit <- Xm*sqrt(wm) ##for fitting truncreg
  Ymfit <- Ym*sqrt(wm) ##for fitting truncreg
  Xmfit <- model.matrix(Ymfit~Xmfit-1) ##for fitting truncreg

  y <- Ym
  y0 <- y <= 0
  y1 <- y > 0
  kx <- kz <- 1
  X <- Z <- rep(1,length(y))
  offsetx <- offsetz <- 0
  sig <- 1

  model_count <- truncreg:::truncreg.fit(X=Xmfit, y=Ymfit, point=0, direction="left", scaled=F) #fit truncated normal
  model_zero <- glm(as.integer(y0)~1,family=binomial,weights=wm)

  start <- list(count = model_count$coefficients, zero = model_zero$coefficients)

  mui <- rep(model_count$coefficients["Xmfit"],length(y))
  sig <- model_count$coefficients["sigma"]
  probi <- model_zero$fitted
  probi <- probi/(probi + (1 - probi) * dnorm(0, mui,sd=sig))
  probi[y1] <- 0
  ll_new <- ziNorm(c(start$count, start$zero),wt=wm)## use global variables: X, Z, y, y0, y1!!!!!!!!!!
  ll_old <- 2 * ll_new
  step <- 0
  while (abs((ll_old - ll_new)/ll_old) > reltol & step < maxstep) {
    step <- step + 1
    ll_old <- ll_new
    Xmfit <- Xm*sqrt(wm*(1 - probi)) ##for fitting truncreg
    Ymfit <- Ym*sqrt(wm*(1 - probi)) ##for fitting truncreg
    Xmfit <- model.matrix(Ymfit~Xmfit-1) ##for fitting truncreg
    model_count <- truncreg:::truncreg.fit(X=Xmfit, y=Ymfit, point=0, direction="left", scaled=F)
    model_zero <- suppressWarnings(glm(probi~1,family=binomial,weights=wm,
                                       start=start$zero))

    mui <- rep(model_count$coefficients["Xmfit"],length(y))
    sig <- model_count$coefficients["sigma"]
    probi <- model_zero$fitted
    probi <- probi/(probi + (1 - probi) * dnorm(0, mui,sd=sig))
    probi[y1] <- 0
    start <- list(count = model_count$coefficients,
                  zero = model_zero$coefficients)
    ll_new <- ziNorm(c(start$count, start$zero),wt=wm)

  }
  if(step>=maxstep) print(paste("max step = ",maxstep, " reached"))
  return(list(count=list(mean=start$count[-length(start$count)],
                         sigma=start$count["sigma"]),
              zero=invlogit(start$zero)))
}

