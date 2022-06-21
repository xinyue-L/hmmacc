#' The HMM-based Sleep/Wake Identification Algorithm to provide sleep information: preparation step
#'
#' It uses the shell function in GGIR to read in accelerometer/actigraph data and pre-process using R. The pre-processed data can be further used for sleep/wake identification.
#'
#' @import"GGIR"
#'
#' @param datadir Directory where the accelerometer files are stored.
#' @param outputdir Directory where the output needs to be stored. Folders are created in this directory to store outputs
#' @param f0 File index to start with (default = 1). Index refers to the filenames sorted in increasing order
#' @param f1 File index to finish with (defaults to number of files available)
#' @param overwrite Overwrite previously generated milestone data by this function for this particular dataset. If FALSE then it will skip the previously processed files (default = FALSE).
#'
#' @keywords sleep variables
#'
#' @return No output will be generated. It is a pre-processing step of accelerometer data in preparation for downstream analysis sleep.part3 through sleep.part6.
#'
#' @examples
#'#datadir = "C:/myfolder/datadir"
#'#outputdir = "C:/myfolder/outputdir"
#'#sleep.prep(datadir=datadir,
#'#outputdir=outputdir,
#'#f0=1,f1=0, overwrite=TRUE)
#'
#' @export
sleep.prep <- function(datadir = c(), outputdir = c(), f0=1, f1=0, overwrite=FALSE, ...) {
  g.shell.GGIR(mode=1:2, datadir = datadir, outputdir = outputdir,
               f0=f0, f1=f1, overwrite = overwrite, ...)
}
