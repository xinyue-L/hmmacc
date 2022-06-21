#' The function to prepare data before implementing the HMM-based Sleep/Wake Identification Algorithm directly without GGIR pre-processing, i.e. prepare summary count data. The example format comes from NHANES.
#'
#' In the sleep/wake identification algorithm, the HMM assumes two states: sleep (state 1)
#' and wake (state 2). As the activity levels follow different distributions under different
#' states, HMM can differentiate the two states accordingly.
#'
#' @param filename The name of the file containing summary count data to be reformatted.
#' @param directory The directory in which the file is stored.
#'
#' @keywords sleep variables
#'
#' @return a csv file named "nightsummary.csv" by default that contains sleep information for all individuals.
#'
#' @examples
#'#df_ex <- read.data("example.csv", "C:/myfolder", sep=",", header=T, stringsAsFactors = F)
read.data <- function(filename,directory="",...) {
  df <- read.table(file=paste(directory,filename,sep="/"),...)

  ##the column variable that contains summary counts: count
  colnames(df) <- c("id","time","count")

  ##time format
  df$time <- as.POSIXlt(df$time)

  return(df)
}
