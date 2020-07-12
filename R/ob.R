#' Activity data for one individual of two non-consecutive days
#'
#' A vector of length 1440, giving activity counts on the log scale.It is for illustration
#' purposes of the functions in the package \code{hmmacc} only.
#'
#' @docType data
#'
#' @usage data(ob)
#'
#' @format An object of class \code{vector}.
#'
#' @keywords datasets
#'
#' @examples
#' data(ob)
#'
#' piI <- c(0.9,0.1)  #state 1 and 2
#' paramI <- list(c(0.4,3,1),  #state 1
#'                c(6,1))      #state 2
#' tranI <- matrix(c(0.9,0.1,0.2,0.8),nrow=2,ncol=2,byrow=T)
#'
#' re <- HMM_ZIG(Yob=ob,paramInit=paramI,tranInit=tranI,piInit=piI,ntime=c(1,841),maxIter=30,reltol=10^-6)
#'
#' @seealso \code{\link{HMM_ZIG}}
