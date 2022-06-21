#' The HMM-based Sleep/Wake Identification Algorithm to summarize the sleep information of all individuals and output in a csv file.
#'
#' In the sleep/wake identification algorithm, the HMM assumes two states: sleep (state 1)
#' and wake (state 2). As the activity levels follow different distributions under different
#' states, HMM can differentiate the two states accordingly.
#'
#' @param metadatadir Directory that holds a folder ’meta’.
#' @param outputdir Directory that the output csv file will be saved into.
#'
#' @keywords sleep variables
#'
#' @return a csv file named allsummary.csv that contains sleep information for all individuals.
#'
#' @examples
#'#metadir = "C:/myfolder/meta" # the meta folder containing the results from previous steps
#'#sleep.part6(metadatadir=metadir)
#'
#' @export
sleep.part6 <- function(metadatadir = c(), outputdir = c()) {

  folder5 <- paste(metadatadir, "/meta/ms5.out.new",sep = "")
  ## check results from step 5
  if (!file.exists(folder5)) {
    cat(paste0("Run step 5 first to generate files in the ms5.out.new folder."))
  } else if (length(list.files(folder5))==0) {
    cat(paste0("The ms5.out.new folder is empty without any data."))
  } else {

    ## pool results
    if(is.null(outputdir)) {
      folder6 <- paste(metadatadir, "/meta/ms6.out.new", sep = "")
    } else {
      folder6 <- outputdir
    }

    if (file.exists(folder6)) {
    } else {
      dir.create(file.path(paste(metadatadir, "/meta",
                                 sep = ""), "ms6.out.new"))
    }

    re <- list()
    file5 <- list.files(folder5)
    for(i in 1:length(file5)) {
      re[[i]] <- get(load(file.path(folder5,file5[i])))
    }
    re <- do.call("rbind",re)

    write.table(re, file = paste0(folder6,"/allsummary.csv"), sep=",",quote=FALSE,
              row.names = F)
  }
}
