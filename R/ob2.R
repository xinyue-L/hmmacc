#' Illustrative data for function HMM_G 
#'
#' A vector of length 440. It is for illustration
#' purposes of the function in the package \code{hmmacc} only.
#'
#' @docType data
#'
#' @usage data(ob2)
#'
#' @format An object of class \code{vector}.
#'
#' @keywords datasets
#'
#' @examples
#' piI <- c(0.9,0.1)  #state 1 and 2
#' paramI <- list(c(3,1),  #state 1
#'                c(6,1))      #state 2
#' tranI <- c(0.9,0.1,0.2,0.8)
#' 
#' re <- HMM_G(Yob=ob2,paramInit=paramI,tranInit=tranI,piInit=piI,ntime=length(ob2),maxIter=30)
#'
#' @seealso \code{\link{HMM_G}}
