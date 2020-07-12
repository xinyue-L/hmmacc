#' The HMM-based Sleep/Wake Identification Algorithm Assuming Gaussian
#'
#' Depending on the device and the calculation method of the summary activity counts, the assumptions for 
#' sleep and wake states may vary. For certain activity count data, there might not be many zeros but
#' small numbers in the sleep state, and in this case, one may consider using Gaussian to approximate the
#' distribution. Make sure to check the empirical distribution and see which assumption is more suitable.
#' 
#' @import"depmixS4"
#' 
#' @param Yob the vector of observations from wearable devices.
#' @param paramInit the initial emission parameters for the sleep state and the wake state.
#' @param tranInit the initial 2 by 2 transition matrix.
#' @param piInit the initial probabilities of length 2: \eqn{\pi_1},\eqn{\pi_2}.
#' @param ntime the vector specifying the start of each independent time series of observations. If not specified, the observations are considered as one single time series.
#' @param maxIter the number of times of trying the fitting procedure. The default is 20. Depending on the initial probabilities assigned, the algorithm may fail to converge. Thus, it is good to try it several times.
#'
#' @keywords Gaussian
#'
#' @return a list. Assume that \code{s} steps in the fitting procedure have been performed:
#' @return \item{param}{the estimated parameters for the sleep state and the wake state.}
#' @return \item{transit}{the estimated transition matrix.}
#' @return \item{init}{the estimated initial probabilities.}
#' @return \item{seq}{the inferred state sequence.}
#' @return \item{iterMax}{logical; if TRUE, the maximal number of iterations have been reached.}

#' @examples
#' data(ob2)
#' 
#' piI <- c(0.9,0.1)  #state 1 and 2
#' paramI <- list(c(3,1),  #state 1
#'                c(6,1))      #state 2
#' tranI <- c(0.9,0.1,0.2,0.8)
#' 
#' re <- HMM_G(Yob=ob2,paramInit=paramI,tranInit=tranI,piInit=piI,ntime=length(ob2),maxIter=20)
#'
#' @export

HMM_G <- function(Yob,paramInit,tranInit,piInit,ntime=NULL,maxIter=20,...) {
  ##depmixS4 for fitting two-state Gaussian
  tmp.df <- data.frame(Y=Yob)
  mod1 <- depmix(response = Y ~ 1, data = tmp.df, nstates = 2, ntimes=ntime,
                 respstart=unlist(paramInit),trstart = tranInit,instart = piInit,...)
  re <- list()
  ##try fitting the model; maximal times: maxIter
  step <- 1; maxtry <- maxIter
  fm1 <- try(fit(mod1))
  while(step <= maxtry & inherits(fm1,"try-error")) {
    step <- step + 1
    fm1 <- try(fit(mod1))
  }
  
  ##write data if the fitting is successful
  if(!inherits(fm1,"try-error")) {
    ##converged or not; converged: "Log likelihood converged to within tol. (relative change)"
    re$itermax <- ifelse(fm1@message=="'maxit' iterations reached in EM without convergence.",TRUE,FALSE) 
    
    ##param & transit
    tmp.param <- c(fm1@response[[1]][[1]]@parameters$coefficients,
                   fm1@response[[1]][[1]]@parameters$sd,
                   fm1@response[[2]][[1]]@parameters$coefficients,
                   fm1@response[[2]][[1]]@parameters$sd)
    tmp.trans <- c(fm1@transition[[1]]@parameters$coefficients,
                   fm1@transition[[2]]@parameters$coefficients)
    tmp.init <- fm1@init
    ##state 1: sedentary; state 2: active
    if(tmp.param[1]>tmp.param[3]) {
      ##switch
      tmp.param <- tmp.param[c(3,4,1,2)]
      tmp.trans <- tmp.trans[c(4,3,2,1)]
      tmp.init <- tmp.init[c(2,1)]
      ##estimated states from Viterbi
      tmp.seq <- ifelse(fm1@posterior$state==1,2,1)
    } else {
      
      ##estimated states from Viterbi
      tmp.seq <- fm1@posterior$state
    }
    re$param <- tmp.param
    re$transit <- tmp.trans
    re$init <- fm1@init
    re$seq <- tmp.seq
  }
  return(re)
}

