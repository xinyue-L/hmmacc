#' The HMM-based Sleep/Wake Identification Algorithm Assuming Zero-Inflated Truncated Gaussian
#'
#' In the sleep/wake identification algorithm, the HMM assumes two states: sleep (state 1)
#' and wake (state 2). State 1 assumes zero-inflated truncated Gaussian, while state 2
#' assumes Gaussian.
#' 
#' It is noted that the estimation of truncated Gaussian using
#' the function \code{truncreg.fit} in the package \code{truncreg} is sometimes unstable.
#' Thus, we recommend using \code{HMM_ZIG}.
#' 
#' @import"truncreg"
#' 
#' @param Yob the vector of observations from wearable devices.
#' @param paramInit the list of initial HMM parameters: list(c(\eqn{\alpha},\eqn{\mu_1},\eqn{\sigma_1}),c(\eqn{\mu_2},\eqn{\sigma_2})).
#' @param tranInit the initial 2 by 2 transition matrix.
#' @param piInit the initial probabilities of length 2: \eqn{\pi_1},\eqn{\pi_2}.
#' @param ntime the vector specifying the start of each independent time series of observations. If not specified, the observations are considered as one single time series.
#' @param maxIter the maximal number of iteration in the fitting procedure. The default is 20.
#' @param reltol the relative tolerance of log likelihood in the fitting procedure. The default is 10^-2.
#'
#' @keywords Zero-Inflated Truncated Gaussian
#'
#' @return a list. Assume that \code{s} steps in the fitting procedure have been performed:
#' @return \item{alpha11}{the estimated zero component \eqn{\alpha_1} for the sleep state (state 1); a vector of length \code{s}.}
#' @return \item{mu1}{the estimated mean \eqn{\mu_1} for the sleep state (state 1); a vector of length \code{s}.}
#' @return \item{sd1}{the estimated standard deviation \eqn{\sigma_1} for the sleep state (state 1); a vector of length \code{s}.}
#' @return \item{mu2}{the estimated mean \eqn{\mu_2} for the wake state (state 2); a vector of length \code{s}.}
#' @return \item{sd2}{the estimated standard deviation \eqn{\sigma_2} for the wake state (state 2); a vector of length \code{s}.}
#' @return \item{loglik}{the estimated log likelihood; a vector of length \code{s}.}
#' @return \item{transit}{the estimated transition matrix; a list of length \code{s}.}
#' @return \item{init}{the estimated initial probabilities; a list of length \code{s}.}
#' @return \item{iterMax}{logical; if TRUE, the maximal number of iterations have been reached.}

#' @examples
#' data(ob)
#' ob <- ob[1:840] #first day observation
#'
#' piI <- c(0.9,0.1)  #state 1 and 2
#' paramI <- list(c(0.4,3,1),  #state 1
#'                c(6,1))      #state 2
#' tranI <- matrix(c(0.9,0.1,0.2,0.8),nrow=2,ncol=2,byrow=T)
#'
#' re <- HMM_ZITG(Yob=ob[1:840],paramInit=paramI,tranInit=tranI,piInit=piI,ntime=c(1),maxIter=30,reltol=10^-6)
#'
#' @export

HMM_ZITG <- function(Yob,paramInit,tranInit,piInit,ntime=NULL,maxIter=20,reltol=10^-2) {
  ####results to be returned
  re <- list(alpha1=numeric(0),
             mu1=numeric(0),
             sd1=numeric(0),
             mu2=numeric(0),
             sd2=numeric(0),
             loglik=numeric(0),
             transit=list(),
             init=list(),
             iterMax=F)

  ## observations and initial state
  if(is.null(ntime)) {
    Yo <- list(Yob)
    init <- list(piInit)
  } else {
    Yo <- init <- list()
    for(i in 1:length(ntime)) {
      init[[i]] <- piInit
      if(i==length(ntime)) {
        ## last sequence of observation
        Yo[[i]] <- Yob[ntime[i]:length(Yob)]
      } else {
        Yo[[i]] <- Yob[ntime[i]:(ntime[i+1]-1)]
      }
    }
  }

  ## initial state
  Ns <- length(paramInit)

  ## transition matrix
  tran <- tranInit
  # state 1
  alpha1 <- paramInit[[1]][1]
  muhat1 <- paramInit[[1]][2]
  sdhat1 <- paramInit[[1]][3]
  # state 2
  muhat2 <- paramInit[[2]][1]
  sdhat2 <- paramInit[[2]][2]
  ## emission probability functions
  emit <- function(o,sig1,mu1,sd1,mu2,sd2,state,dlog=F) {
    return(switch(state,
                  dZITGaussian(o,sig1,mu1,sd1,log=dlog),
                  dnorm(o,mu2,sd2,log=dlog)))
  }

  #####################################################################
  ##
  ## auxillary functions
  ##
  #####################################################################
  invlogit <- function(x) exp(x)/(exp(x)+1)

  ###############################
  ###  forward algorithm
  ###############################
  ##ai_j: i time; j state
  forward <- function(Y,initprob) {
    Nt <- length(Y)
    a <- matrix(nrow=2,ncol=Nt)
    a_aux <- numeric(Nt)
    a[1,1] <- initprob[1]*emit(Y[1],alpha1,muhat1,sdhat1,muhat2,sdhat2,1,dlog=F)
    a[2,1] <- initprob[2]*emit(Y[1],alpha1,muhat1,sdhat1,muhat2,sdhat2,2,dlog=F)

    a_aux[1] <- sum(a[,1])
    a[,1] <- a[,1]/a_aux[1]
    ##go from 2 to Nt
    for (i in 2:Nt) {
      a[,i] <- 0
      ##a_t(k)
      for (k in 1:Ns) {
        ##summation: sum_i^Ns
        for (j in 1:Ns) {
          a[k,i] <- a[k,i] + a[j,i-1]*tran[j,k] ##a_k(j)*aj_k
        }
        a[k,i] <- a[k,i]*emit(Y[i],alpha1,muhat1,sdhat1,muhat2,sdhat2,k,dlog=F) ##b_k(Y_i)
      }
      a_aux[i] <- sum(a[,i])
      a[,i] <- a[,i]/a_aux[i]
    }
    return(list(a=a,loglik=sum(log(a_aux))+log(sum(a[,Nt]))))
  }

  ###############################
  ###  backward algorithm
  ###############################
  ##bi_j: i time; j state
  backward <- function(Y) {
    Nt <- length(Y)
    b <- matrix(nrow=2,ncol=Nt)
    b_aux <- numeric(Nt)
    b[,Nt] <- 1
    b_aux[Nt] <- sum(b[,Nt])
    b[,Nt] <- b[,Nt]/b_aux[Nt]

    ##go from Nt-1 to 1
    for (t in (Nt-1):1) {
      b[,t] <- 0
      ##b_t(i)
      for (i in 1:Ns) {
        ##summbtion: sum_i^Ns
        for (j in 1:Ns) {
          ##\sum_j ai_j*b_j(O_{t+1})*\beta_{t+1}(j)
          b[i,t] <- b[i,t]+tran[i,j]*emit(Y[t+1],alpha1,muhat1,sdhat1,muhat2,sdhat2,j,dlog=F)*b[j,t+1]
        }
      }
      b_aux[t] <- sum(b[,t])
      b[,t] <- b[,t]/b_aux[t]
    }
    return(list(b=b))
  }

  ###############################
  ###  Baum-Welch Algorithm
  ###  for matrix e and matrix g
  ###############################
  ## p{q_i|O,lambda}
  ## g_t(i) <- a[i,t]*b[i,t]/sum
  ## e_t(i,j)=p{q_t=i,q_{t+1}=j,O|\lambda}/p{O|\lambda}
  ## e_t(i,j): state i at t; state j at t+1
  bw <- function(Y,a,b) {
    Nt <- length(Y)

    e <- matrix(nrow=Ns*Ns,ncol=Nt-1)
    eindex <- function(ii,jj,Nss) (ii-1)*Nss+jj
    for (t in 1:(Nt-1)) {
      e[,t] <- 0
      for (i in 1:Ns) {
        for (j in 1:Ns) {
          e[eindex(i,j,Ns),t] <- a[i,t]*tran[i,j]*b[j,t+1]*emit(Y[t+1],alpha1,muhat1,sdhat1,muhat2,sdhat2,j,dlog=F)
        }
      }
      temp <- sum(e[,t])
      e[,t] <- e[,t]/temp
    }

    g2 <- matrix(nrow=Ns,ncol=Nt)
    for (t in 1:Nt) {
      g2[,t] <- 0
      for (i in 1:Ns) {
        g2[i,t] <- a[i,t]*b[i,t]
      }
      temp <- sum(g2[,t])
      g2[,t] <- g2[,t]/temp
    }

    return(list(e=e,g=g2))
  }

  ###############################
  ###  Viterbi Algorithm: on the log scale
  ###############################
  ##di_j: i time; j state
  va <- function(Y,initprob) {
    Nt <- length(Y)
    d <- matrix(nrow=Ns,ncol=Nt)
    d[1,1] <- log(initprob[1])+emit(Y[1],alpha1,muhat1,sdhat1,muhat2,sdhat2,1,dlog=T)
    d[2,1] <- log(initprob[2])+emit(Y[1],alpha1,muhat1,sdhat1,muhat2,sdhat2,2,dlog=T)
    ##dmax: record arg max d_{i-1}(j)*aj_k
    dmax <- matrix(nrow=Ns,ncol=Nt)
    ##go from 2 to Nt
    for (i in 2:Nt) {
      d[,i] <- 0
      ##d_t(k)
      for (k in 1:Ns) {
        ##find max d_{i-1}(j)*aj_k: store as temp
        temp <- numeric(Ns)
        for (j in 1:Ns) {
          temp[j] <- d[j,(i-1)]+log(tran[j,k]) ## d[j,(i-1)] already on log scale
        }
        dmax[k,i] <- which.max(temp)
        d[k,i] <- max(temp)+emit(Y[i],alpha1,muhat1,sdhat1,muhat2,sdhat2,k,dlog=T)
      }
    }

    ##winning states
    ds <- numeric(Nt)
    ds[Nt] <- which.max(d[,Nt])
    for (i in (Nt-1):1) {
      ds[i] <- dmax[ds[i+1],i+1]
    }
    return(ds)
  }

  #max number of iterations
  step <- 1; loglik_ori <- 1; loglik_delta <- 100
  while (step < maxIter & (abs(loglik_delta/loglik_ori) > reltol)) {

    ##forward, backward, and BW algorithm
    re_a <- lapply(1:length(Yo),function(x) forward(Yo[[x]],init[[x]]))
    re_b <- lapply(Yo,backward)
    re_bw <- lapply(1:length(Yo),function(x) bw(Yo[[x]],re_a[[x]]$a,re_b[[x]]$b))

    ##likelihood: p(O|lambda)
    re$loglik[step] <- loglik <- sum(unlist(lapply(re_a,function(x) x$loglik)))
    print(paste("step = ",step,"; loglik = ",round(loglik,3), sep=""))

    ## transition matrix estimation
    tem_e <- rowSums(do.call("cbind",lapply(re_bw,function(x) rowSums(x$e))))
    tem_g <- rowSums(do.call("cbind",lapply(re_bw,function(x) rowSums(x$g[,-ncol(x$g)]))))
    re$transit[[step]] <- tran <- matrix(tem_e/rep(tem_g,each=2),nrow=Ns,byrow=T)

    ## initial probability estimation
    re$init[[step]] <- init <- lapply(re_bw,function(x) x$g[,1])

    ## state 1: zero-inflated Gaussian estimation
    tem_g2 <- do.call("cbind",lapply(re_bw, function(x) x$g))
    b_hat1 <- ZITGaussian(Yob,weights=tem_g2[1,])
    re$alpha1[step] <- alpha1 <- b_hat1$zero
    re$mu1[step] <- muhat1 <- b_hat1$count$mean
    re$sd1[step] <- sdhat1 <- b_hat1$count$sigma

    ## state 2: Gaussian estimation
    b_hat2 <- Gaussian(Yob,weights=tem_g2[2,])
    re$mu2[step] <- muhat2 <- b_hat2$coefficients
    re$sd2[step] <- sdhat2 <- b_hat2$dispersion

    ##log likelihood change
    loglik_delta <- abs(loglik_ori-loglik)
    loglik_ori <- loglik
    step <- step + 1
  }
  if(step>=maxIter) {
    re$iterMax <- T
    print(paste("max EM iteration = ",maxIter, " reached"))
  }

  ## winning state
  re$states <- unlist(lapply(1:length(Yo),function(x) va(Yo[[x]],init[[x]])))

  return(re)
}
