#' The HMM-based Sleep/Wake Identification Algorithm to output sleep variables
#'
#' In the sleep/wake identification algorithm, the HMM assumes two states: sleep (state 1)
#' and wake (state 2). As the activity levels follow different distributions under different
#' states, HMM can differentiate the two states accordingly.
#'
#' @importFrom"zoo"
#' rollapply
#'
#' @param metadatadir Directory that holds a folder ’meta’.
#' @param datadir Directory where the accelerometer files are stored or list of accelerometer filenames and directories.
#' @param f0 File index to start with (default = 1, and usually working with the default is fine). Index refers to the file names sorted in the increasing order.
#' @param f1 File index to finish with (defaults to the number of files available, and setting to 0 is fine).
#' @param daysleeplength the threshold for calculating daytime sleep. The default value is 0.5 in the unit of hour (= 30 minutes), meaning that daytime sleep episodes that are longer than 30 minutes will be considered.
#' @param nonwearhour Remove the detected nonwear intervals that are shorter than L hours. The default value of \code{L} is 2 (hours).
#' @param overwrite Whether to overwrite data that were generated by this function previously. If overwrite = FALSE then it will not overwrite but skip the previously generated files (default = FALSE).
#' @param do.parallel Whether to use multi-core processing. Parallel processing only works if at least 4 CPU cores are available (default = TRUE).
#'
#'
#' @keywords sleep variables
#'
#' @return a dataset on sleep variables.
#' @return \item{sleeponset}{the evening sleep onset timing for the night sleep, expressed in the unit of hour. For example, value = 21 means 21:00/9pm.}
#' @return \item{wakeup}{the morning wakeup timing for the night sleep, expressed in the unit of hour. For example, value = 8 means 08:00/8am.}
#' @return \item{SptDuration}{the sleep duration, calculated as the time between sleep onset and wakeup, expressed in the unit of hour. For example, value = 10 means 10 hours of nighttime sleep.}
#' @return \item{guider_onset}{HMM crude guider of the sleep onset timing, expressed in the unit of hour. For example, value = 21 means 21:00/9pm.}
#' @return \item{guider_wakeup}{HMM crude guider of the wakeup timing, expressed in the unit of hour. For example, value = 8 means 08:00/8am.}
#' @return \item{guider_SptDuration}{The crude estimate of sleep duration by HMM crude guiders of sleep onset and wakeup.}
#' @return \item{sleeponset_ts}{the sleep onset timing, expressed in the format hh:mm:ss.}
#' @return \item{wakeup_ts}{the wakeup timing, expressed in the format hh:mm:ss.}
#' @return \item{guider_onset_ts}{HMM guider of the sleep onset timing, expressed in the format hh:mm:ss.}
#' @return \item{guider_wakeup_ts}{HMM guider of the wakeup timing, expressed in the format hh:mm:ss.}
#' @return \item{sleep_true}{calculated as the amount of true sleep time in the defined sleep duration period.}
#' @return \item{sleep_efficiency}{calculated as the amount of true sleep time divided by the night sleep duration.}
#' @return \item{waso}{wake after sleep onset, calculated as the amount of wake time in the defined sleep duration period.}
#' @return \item{sleep_daytime_Nmin_more}{daytime sleep time, only counting the period that is longer than N minutes. The default N = 30 minutes.}
#' @return \item{sleep_day_episode_Nmin_more}{the number of daytime sleep episodes, only counting the period that is longer than N minutes. The default N = 30 minutes.}
#' @return \item{sleep_total}{total sleep time – 24 hr assessment of sleep.}
#' @return \item{L5_midpoint}{the midpoint of L5, which is the least 5 active hour period of the 24 hour day. For example, value = 2 means that the midpoint of L5 is 2 hours/2am.}
#' @return \item{M10_midpoint}{the midpoint of M10, which is the most active 10 hour period of the 24 hour day. For example, value = 16 means that the midpoint of M10 is 16 hours/4pm. Limited by each day data, NA will be reported as -999.}
#'
#' @examples
#'#metadir = "C:/myfolder/meta" # the meta folder containing the results from previous steps
#'#sleep.part5(metadatadir=metadir, overwrite=TRUE)
#'
#' @export
sleep.part5 <- function(metadatadir = c(), f0=1, f1=0, datadir = c(),
                        daysleeplength = 0.5, nonwearhour = 2,
                        overwrite = FALSE, do.parallel = TRUE) {

  if (file.exists(paste(metadatadir, sep = ""))) {
  } else {
    dir.create(file.path(metadatadir))
  }
  if (file.exists(paste(metadatadir, "/meta/ms5.out.new",
                        sep = ""))) {
  } else {
    dir.create(file.path(paste(metadatadir, "/meta",
                               sep = ""), "ms5.out.new"))
  }

  fnames = dir(paste(metadatadir, "/meta/ms2.out", sep = ""))
  if (f1 > length(fnames) | f1 == 0)
    f1 = length(fnames)
  if (f0 > length(fnames) | f0 == 0)
    f0 = 1
  ffdone = fdone = dir(paste(metadatadir, "/meta/ms5.out.new",
                             sep = ""))
  if (length(fdone) > 0) {
    for (ij in 1:length(fdone)) {
      tmp = unlist(strsplit(fdone[ij], ".RData"))
      ffdone[ij] = tmp[1]
    }
  }
  else {
    ffdone = c()
  }
  nightsperpage = 7
  if (do.parallel == TRUE) {
    closeAllConnections()
    cores = parallel::detectCores()
    Ncores = cores[1]
    if (Ncores > 3) {
      cl <- parallel::makeCluster(Ncores - 1)
      doParallel::registerDoParallel(cl)
    } else {
      cat(paste0("\nparallel processing not possible because number of available cores (",
                 Ncores, ") < 4"))
      do.parallel = FALSE
    }
  }
  t1 = Sys.time()
  if (do.parallel == TRUE) {
    cat(paste0("\n Busy processing ... see ", metadatadir,
               "/ms3", " for progress\n"))
  }

  GGIRinstalled = is.element("GGIR", installed.packages()[,
                                                          1])
  packages2passon = functions2passon = NULL
  GGIRloaded = "GGIR" %in% .packages()
  if (GGIRloaded) {
    packages2passon = c("GGIR","depmixS4")
    errhand = "pass"
  }
  else {
    functions2passon = c("g.sib.det", "g.detecmidnight",
                         "iso8601chartime2POSIX", "g.sib.plot",
                         "g.sib.sum","depmix")
    errhand = "stop"
  }

  ##nonwear write column names
  nonwear = t(c("ID","filename","day","total_hour_observation","nonwear_hour_observation","nonwear_percentage"))
  write.table(nonwear, file = paste0(metadatadir, "/meta/nonwear.csv"), sep=",",quote=FALSE,
              row.names = F,col.names=F)

  fe_dopar = foreach::`%dopar%`
  fe_do = foreach::`%do%`
  i = 0
  `%myinfix%` = ifelse(do.parallel, fe_dopar, fe_do)
  output_list = foreach::foreach(i = f0:f1, .packages = packages2passon,
                                 .export = functions2passon, .errorhandling = errhand) %myinfix%
    {
      tryCatchResult = tryCatch({
        FI = file.info(paste(metadatadir, "/meta/ms2.out/",
                             fnames[i], sep = ""))
        if (is.na(FI$size) == TRUE)
          FI$size = 0
        if (FI$size == 0 | is.na(FI$size) == TRUE | length(FI$size) ==
            0) {
          cat(paste("ms2.out file ", fnames[i], sep = ""))
          cat("Filename not recognised")
        }
        fname = unlist(strsplit(fnames[i], ".RData"))[1]
        if (length(ffdone) > 0) {
          if (length(which(ffdone == fname)) > 0) {
            skip = 1
          } else {
            skip = 0
          }
        } else {
          skip = 0
        }
        if (overwrite == TRUE)
          skip = 0
        if (skip == 0) {
          cat(paste(" ", i, sep = ""))
          load(paste(metadatadir, "/meta/basic/meta_",
                     fnames[i], sep = ""))
          load(paste(metadatadir, "/meta/ms3.out.new/",
                     fnames[i], sep = ""))
          load(paste(metadatadir, "/meta/ms4.out.new/",
                     fnames[i], sep = ""))

          ##output file
          #-sleep onset (for night time sleep)
          #-wakeup time (for morning wakeup)
          #-sleep duration period (calculated as time between sleep onset and wake up time)
          output <- nightsummary[,c("ID","night","calendar_date","weekday","sleeponset","wakeup","SptDuration",
                                    "guider_onset","guider_wakeup","guider_SptDuration",
                                    "sleeponset_ts","wakeup_ts","guider_onset_ts","guider_wakeup_ts","filename")]

          #-auxillary functions
          time2index <- function(x) {
            if(x>=24) x <- x-24
            h <- floor(x)
            m <- floor((x-h)*60)
            s <- round((x-h-m/60)*3600)
            if(nchar(h)==1) h <- paste0(0,h)
            if(nchar(m)==1) m <- paste0(0,m)
            if(nchar(s)==1) s <- paste0(0,s)
            return(paste(h,m,s,sep=":"))
          }
          time2index <- Vectorize(time2index,vec="x")


          index2time24 <- function(x) {
            timeRecord <- unlist(as.POSIXlt(x, format = "%Y-%m-%dT%H:%M:%S",""))
            sum(as.numeric(timeRecord[c("hour","min", "sec")])/c(1, 60, 3600))
          }
          index2time24 <- Vectorize(index2time24,ve="x")

          index2time24_plus <- function(x) {
            timeRecord <- unlist(as.POSIXlt(x, format = "%Y-%m-%dT%H:%M:%S",""))
            re <- sum(as.numeric(timeRecord[c("hour","min", "sec")])/c(1, 60, 3600))
            if(re<12) re <- re+24
            return(re)
          }
          index2time24_plus <- Vectorize(index2time24_plus,ve="x")

          #-mid-point of sleep (middle of night time sleep duration)
          nightsummary$midpoint_sleep <- time2index((nightsummary$sleeponset+nightsummary$wakeup)/2)

          #-sleep efficiency (calculated as time spent sleep during total sleep time)
          sib.cla.sum$sib.onset.time24 <- index2time24_plus(sib.cla.sum$sib.onset.time)
          sib.cla.sum$sib.end.time24 <- index2time24_plus(sib.cla.sum$sib.end.time)
          output$sleep_true <- sapply(unique(sib.cla.sum$night),function(j) {
            df_night <- sib.cla.sum[sib.cla.sum$night==j,]
            edge <- c(nightsummary$sleeponset[which(nightsummary$night==j)],nightsummary$wakeup[which(nightsummary$night==j)])
            tm_st <- which(df_night$sib.onset.time24>=(edge[1]-0.00001))[1]
            tm_end <- which(df_night$sib.end.time24>=(edge[2]-0.00001))[1]
            sum(df_night$tot.sib.dur.hrs[tm_st:tm_end])
          })
          output$sleep_efficiency <- output$sleep_true/output$SptDuration


          #-WASO/wake after sleep onset (calculated as the amount of wake time in defined sleep duration period)
          output$waso <- output$SptDuration-output$sleep_true


          #-total daytime sleep -  sleep time outside of previously defined sleep duration period
          output$sleep_daytime_15min_more <- sapply(unique(sib.cla.sum$night),function(j) {
            df_night <- sib.cla.sum[sib.cla.sum$night==j,]
            edge <- c(nightsummary$sleeponset[which(nightsummary$night==j)],nightsummary$wakeup[which(nightsummary$night==j)])
            tm_st <- which(df_night$sib.onset.time24>=(edge[1]-0.00001))[1]
            tm_end <- which(df_night$sib.end.time24>=(edge[2]-0.00001))[1]
            df_select <- df_night$tot.sib.dur.hrs[-c(tm_st:tm_end)]
            sum(df_select[which(df_select>=daysleeplength)]) #at least 15 minutes for daytime sleep
          })

          #-total daytime sleep -  sleep time outside of previously defined sleep duration period
          output$sleep_day_episode_15min_more <- sapply(unique(sib.cla.sum$night),function(j) {
            df_night <- sib.cla.sum[sib.cla.sum$night==j,]
            edge <- c(nightsummary$sleeponset[which(nightsummary$night==j)],nightsummary$wakeup[which(nightsummary$night==j)])
            tm_st <- which(df_night$sib.onset.time24>=(edge[1]-0.00001))[1]
            tm_end <- which(df_night$sib.end.time24>=(edge[2]-0.00001))[1]
            df_select <- df_night$tot.sib.dur.hrs[-c(tm_st:tm_end)]
            sum(df_select>=daysleeplength) #at least 15 minutes for daytime sleep
          })

          #-total sleep time – 24 hr assessment of sleep
          output$sleep_total <- output$sleep_daytime+output$SptDuration

          #-sleep latency (for sleep log guided HMM algorithm: time between reported sleep and sleep onset)
          #unavailable right now

          #-L5 – least 5 active hour period of the 24 hour day
          #-M10 – most active 10 hour period of the 24 hour day
          bwidth <- (60/M$windowsizes[1])*5 #precision 5-min is ok
          enmo5 <- data.frame(timestamp=M$metashort$timestamp[(1:(nrow(M$metashort)/(bwidth)))*bwidth-(bwidth-1)],
                              ENMO5=sapply(1:(nrow(M$metashort)/(bwidth)),function(i) sum(M$metashort$ENMO[((i-1)*bwidth+1):(i*bwidth)])),
                              nday=A$night[(1:(nrow(M$metashort)/(bwidth)))*bwidth-(bwidth-1)])
          enmo5$time24 <- index2time24(enmo5$timestamp)

          #first find L5
          roll5 <- zoo::rollapply(enmo5$ENMO5,((60/M$windowsizes[1])*60*5)/bwidth,mean)
          roll5 <- data.frame(time24=enmo5$time24[1:length(roll5)],roll5=roll5,nday=enmo5$nday[1:length(roll5)])
          L5_index <- sapply(nightsummary$night,function(i) {
            tmp_roll <- roll5[roll5$nday==i,]
            which(roll5$nday==i & roll5$time24==tmp_roll$time24[which.min(tmp_roll$roll5)[1]])
          })
          output$L5_midpoint <- roll5$time24[L5_index]+2.5 # to adjust to midpoint
          output$L5_midpoint[which(output$L5_midpoint>=24)] <- output$L5_midpoint[which(output$L5_midpoint>=24)]-24

          #then find M10: between L5; day 1 NA -999
          roll10 <- zoo::rollapply(enmo5$ENMO5,((60/M$windowsizes[1])*60*10)/bwidth,mean)
          roll10<- data.frame(time24=enmo5$time24[1:length(roll10)],roll10=roll10)
          L5_index2 <- L5_index
          if(L5_index2[length(L5_index2)]>length(roll10)) L5_index2[length(L5_index2)] <- length(roll10)
          output$M10_midpoint <- c(-999,sapply(1:(length(L5_index2)-1),function(i) {
            tmp_roll <- roll10[L5_index[i]:L5_index[i+1],]
            tmp_roll$time24[which.max(tmp_roll$roll10)[1]]+5 # to adjust to midpoint
          }))


          #non-wear data output
          hourL <- 60*60/M$windowsizes[1]
          nonwear <- data.frame(ID=output$ID[1],filename=output$filename[1],day=unique(A.ori$night))
          nonwear$total_hour_observation <- sapply(nonwear$day,function(k) sum(A.ori$night==k)/hourL)

          if (length(which(A.ori$invalid == 1)) > 0) {

            ##remove < 2-hour non-wear intervals
            invalid_rle <- rle(A.ori$invalid)
            invalid_rle$values[which(invalid_rle$values==1 & invalid_rle$lengths<=(60/M$windowsizes[1])*60*nonwearhour)] <- 0
            A.ori$invalid <- inverse.rle(invalid_rle)
          }
          nonwear$nonwear_hour_observation <- sapply(nonwear$day,function(k) sum(A.ori$night==k & A.ori$invalid==1)/hourL)
          nonwear$nonwear_percentage <- nonwear$nonwear_hour_observation/nonwear$total_hour_observation
          write.table(nonwear, file = paste0(metadatadir, "/meta/nonwear.csv"), sep=",",quote=FALSE,
                      col.names= F, row.names = F, append = TRUE)

          save(output,
               file = paste(metadatadir, "/meta/ms5.out.new/",
                            fname, ".RData", sep = ""))
        }
      })
      return(tryCatchResult)
    }
  if (do.parallel == TRUE) {
    on.exit(parallel::stopCluster(cl))
    for (oli in 1:length(output_list)) {
      if (is.null(unlist(output_list[oli])) == FALSE) {
        cat(paste0("\nErrors and warnings for ",
                   fnames[oli]))
        print(unlist(output_list[oli]))
      }
    }
  }
}
