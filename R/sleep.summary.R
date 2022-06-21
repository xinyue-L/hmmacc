#' The function to implement the HMM-based Sleep/Wake Identification Algorithm directly without GGIR pre-processing, i.e. using summary count data. The example format comes from NHANES.
#'
#' In the sleep/wake identification algorithm, the HMM assumes two states: sleep (state 1)
#' and wake (state 2). As the activity levels follow different distributions under different
#' states, HMM can differentiate the two states accordingly.
#'
#' @param df The directory that holds a folder ’meta’.
#' @param interval The summary count data interval in seconds, i.e. summary count data every 5 seconds, 30 seconds, 60 seconds (interval=60). The default is interval=60 (60 seconds or 1 minute).
#' @param minhour The minimal number of hours available to be considered for night analysis; otherwise sleep parameters will not be calculated for the night. The default is minhour = 16.
#' @param maxtry The maximal number of trials that the program will run to fit HMM. Randomness is involved in initiation and there can be convergence issues due to data quality. The default is maxtry = 20.
#' @param filename The name of the file to store results. The default is "nightsummary.csv".
#' @param outputdir The directory that the output csv file will be saved into.
#'
#' @keywords sleep variables
#'
#' @return a csv file named "nightsummary.csv" by default that contains sleep information for all individuals.
#'
#' @examples
#'#df_ex <- read.data("example.csv", "C:/myfolder", sep=",", header=T, stringsAsFactors = F)
#'#sleep.summary(df_ex,interval = 60, minhour = 16, maxtry = 20, filename = "example.csv", outputdir = "C:/outdir")
#'
#' @export
sleep.summary <- function(df,interval = 60, minhour = 16, maxtry = 20,
                          filename = "", outputdir = "") {

  ############################################################# zero step
  hours <- df$time$hour
  tmp_12 <- as.numeric(hours==12)
  tmp_12_2 <- rle(tmp_12)
  tmp_12_3 <- cumsum(tmp_12_2$lengths)
  tmp_12_3 <- c(1,tmp_12_3[-length(tmp_12_3)]+1)
  tmp_12_4 <- tmp_12_3[which(tmp_12_2$values==1)]
  tmp_12_5 <- c(tmp_12_4[1]-1,tmp_12_4[-1]-tmp_12_4[-length(tmp_12_4)],
                length(hours)-tmp_12_4[length(tmp_12_4)]+1)
  tmp_12_5 <- tmp_12_5[tmp_12_5!=0]
  night <- unlist(lapply(1:length(tmp_12_5),function(s) rep(s,tmp_12_5[s])))
  rm(tmp_12,tmp_12_2,tmp_12_3,tmp_12_4,tmp_12_5)
  valid_night <- as.numeric(names(table(night))[which(table(night)>=16*60*60/interval)])
  ##count data
  #tmp_log <- log(df$count+1)
  ############################################################# first step

  ##depmixS4 for fitting two-state
  tmp_df <- data.frame(Y=df$count,
                       Ylog=log(df$count+1),
                       hour=hours)
  mod1 <- depmix(response = Y ~ 1, data = tmp_df, nstates = 2,family=poisson(),
                 trstart = c(0.9,0.1,0.1,0.9),instart = c(0.9,0.1))

  ##try fitting the model; maximal 20 times
  step <- 1
  fm1 <- try(fit(mod1))
  while(step <= maxtry & inherits(fm1,"try-error")) {
    step <- step + 1
    fm1 <- try(fit(mod1))
  }

  ##write data if the fitting is successful
  if(!inherits(fm1,"try-error")) {

    ##state 1: sedentary; state 2: active
    if(fm1@response[[1]][[1]]@parameters$coefficients > fm1@response[[2]][[1]]@parameters$coefficients) {

      ##estimated states from Viterbi
      tmp.state <- ifelse(fm1@posterior$state==1,2,1)
    } else {

      ##estimated states from Viterbi
      tmp.state <- fm1@posterior$state
    }

    tmp.ori <- tmp.state
    tmp.rle <- rle(tmp.state)
    while(sum(tmp.rle$lengths<=5)>0) {
      tmp.rle$values[which.min(tmp.rle$lengths)] <- switch(tmp.rle$values[which.min(tmp.rle$lengths)],2,1)
      tmp.rle <- rle(inverse.rle(tmp.rle))
    }
    tmp.state <- inverse.rle(tmp.rle)

  } else {
    ## try Gaussian
    mod1 <- depmix(response = Ylog ~ 1, data = tmp_df, nstates = 2,family=gaussian(),
                   trstart = c(0.9,0.1,0.1,0.9),instart = c(0.9,0.1))

    ##try fitting the model; maximal 20 times
    step <- 1
    fm1 <- try(fit(mod1))
    while(step <= maxtry & inherits(fm1,"try-error")) {
      step <- step + 1
      fm1 <- try(fit(mod1))
    }

    ##write data if the fitting is successful
    if(!inherits(fm1,"try-error")) {

      ##state 1: sedentary; state 2: active
      if(fm1@response[[1]][[1]]@parameters$coefficients > fm1@response[[2]][[1]]@parameters$coefficients) {

        ##estimated states from Viterbi
        tmp.state <- ifelse(fm1@posterior$state==1,2,1)
      } else {

        ##estimated states from Viterbi
        tmp.state <- fm1@posterior$state
      }

      tmp.ori <- tmp.state
      tmp.rle <- rle(tmp.state)
      while(sum(tmp.rle$lengths<=5)>0) {
        tmp.rle$values[which.min(tmp.rle$lengths)] <- switch(tmp.rle$values[which.min(tmp.rle$lengths)],2,1)
        tmp.rle <- rle(inverse.rle(tmp.rle))
      }
      tmp.state <- inverse.rle(tmp.rle)

    }
  }


  ############################################################# second step
  bwidth <- (60/interval)*5*6
  Y_guide <- sapply(1:(nrow(tmp_df)/(bwidth)),
                         function(i) sum(tmp_df$Y[((i-1)*bwidth+1):(i*bwidth)]))


  m2 <- depmix(y~1,ns=2,family=gaussian(),data=data.frame(y=Y_guide))
  fm2 <- fit(m2,em=em.control(maxit=1000))

  if(fm2@response[[1]][[1]]@parameters$coefficients < fm2@response[[2]][[1]]@parameters$coefficients) {
    state <- fm2@posterior$state
  } else {
    state <- ifelse(fm2@posterior$state==2,1,2)
  }

  tmp_rle <- rle(state)
  while(sum(tmp_rle$lengths<=2)>0) {
    tmp_rle$values[which.min(tmp_rle$lengths)] <- switch(tmp_rle$values[which.min(tmp_rle$lengths)],2,1)
    tmp_rle <- rle(inverse.rle(tmp_rle))
  }
  state <- inverse.rle(tmp_rle)
  guide <- rep(c(state,state[length(state)]),each=bwidth)
  guide <- guide[1:nrow(df)]

  df_state <- data.frame(time=df$time,
                         ori=tmp.ori,
                         guide=guide,
                         night=night,
                         state=tmp.state)


  index2time24 <- function(x) {
    timeRecord <- unlist(as.POSIXlt(x, format = "%Y-%m-%dT%H:%M:%S",""))
    sum(as.numeric(timeRecord[c("hour","min", "sec")])/c(1, 60, 3600))
  }
  index2time24 <- Vectorize(index2time24,ve="x")

  ############################################################# third step
  #id	filename	weekday	calendardate	night_number	acc_onset	daysleeper
  #cleaningcode	acc_available	window_length_in_hours	dur_nightsleep_min
  #dur_day_min	dur_night_min	dur_nightandday_min	N_atleast5minwakenight
  #sleep_efficiency	L5TIME	L5VALUE	L5TIME_num	M10TIME	M10VALUE	M10TIME_num
  output = as.data.frame(matrix(0, length(valid_night), 10))
  colnames(output) <- c("ID", "date", "weekday", "night_number",
                        "sleeponset", "wakeup", "sleep_duration", "sleep_efficiency",
                        "L5_midpoint", "M10_midpoint")

  output$ID <- filename
  output$night_number <- valid_night
  output$date <- do.call("c",lapply(valid_night,function(s) as.Date(df_state$time[df_state$night==s][1])))
  output$weekday <- weekdays(output$date)

  if(!inherits(fm1,"try-error")) {

    ##sleep duration and sleep_efficiency
    for(i in 1:nrow(output)) {
      df_tmp <- df_state[df_state$night==output$night_number[i],]
      guide_index <- which(df_tmp$guide==1)
      tmp_st <- rle(df_tmp$state)
      tmp_index <- c(1,cumsum(tmp_st$lengths[-length(tmp_st$lengths)])+1)
      tmp1 <- tmp_index[tmp_st$values==1]
      st <- tmp1[which.min(abs(tmp1-guide_index[1]))]

      tmp_index2 <- cumsum(tmp_st$lengths)
      tmp2 <- tmp_index2[tmp_st$values==1]
      end <- tmp2[which.min(abs(tmp2-guide_index[length(guide_index)]))]

      output$sleeponset[i] <- index2time24(df_tmp$time[st])
      output$wakeup[i] <- index2time24(df_tmp$time[end])

      output$sleep_efficiency[i] <- sum(df_tmp$ori[st:end]==1)/(end-st)

      output$sleep_duration[i] <- (end-st)*interval/60/60

    }

  } else {

    ##sleep duration and sleep_efficiency
    for(i in 1:nrow(output)) {
      df_tmp <- df_state[df_state$night==output$night_number[i],]
      guide_index <- which(df_tmp$guide==1)
      st <- guide_index[1]
      end <- guide_index[length(guide_index)]
      output$sleeponset[i] <- index2time24(df_tmp$time[st])
      output$wakeup[i] <- index2time24(df_tmp$time[end])

      output$sleep_efficiency[i] <- sum(df_tmp$ori[st:end]==1)/(end-st)
      output$sleep_duration[i] <- (end-st)*interval/60/60
    }

  }



  #-L5 – least 5 active hour period of the 24 hour day
  #-M10 – most active 10 hour period of the 24 hour day
  bwidth2 <- (60/interval)*5 #precision 5-min is ok
  enmo5 <- data.frame(time=df$time[(1:(nrow(df)/(bwidth2)))*bwidth2-(bwidth2-1)],
                      ENMO5=sapply(1:(nrow(df)/(bwidth2)),function(i) sum(df$count[((i-1)*bwidth2+1):(i*bwidth2)])),
                      nday=night[(1:(nrow(df)/(bwidth2)))*bwidth2-(bwidth2-1)])
  enmo5$time24 <- index2time24(enmo5$time)

  #first find L5
  roll5 <- zoo::rollapply(enmo5$ENMO5,((60/interval)*60*5)/bwidth2,mean)
  roll5 <- data.frame(time24=enmo5$time24[1:length(roll5)],roll5=roll5,nday=enmo5$nday[1:length(roll5)])
  L5_index <- sapply(valid_night,function(i) {
    tmp_roll <- roll5[roll5$nday==i,]
    which(roll5$nday==i & roll5$time24==tmp_roll$time24[which.min(tmp_roll$roll5)[1]])
  })
  output$L5_midpoint <- roll5$time24[L5_index]+2.5 # to adjust to midpoint
  output$L5_midpoint[which(output$L5_midpoint>=24)] <- output$L5_midpoint[which(output$L5_midpoint>=24)]-24

  #then find M10: between L5; day 1 NA -999
  roll10 <- zoo::rollapply(enmo5$ENMO5,((60/interval)*60*10)/bwidth2,mean)
  roll10<- data.frame(time24=enmo5$time24[1:length(roll10)],roll10=roll10)
  L5_index2 <- L5_index
  if(L5_index2[length(L5_index2)]>length(roll10)) L5_index2[length(L5_index2)] <- length(roll10)
  output$M10_midpoint <- c(-999,sapply(1:(length(L5_index2)-1),function(i) {
    tmp_roll <- roll10[L5_index[i]:L5_index[i+1],]
    tmp_roll$time24[which.max(tmp_roll$roll10)[1]]+5 # to adjust to midpoint
  }))

  if(outputdir=="") outputdir <- getwd()

  output_files <- list.files(outputdir)
  if(filename == "") filename = "nightsummary.csv"
  if(filename %in% output_files) {
    write.table(output, file = paste0(outputdir,"/",filename), sep=",",quote=FALSE,
                row.names = F, col.names = F, append = TRUE)
  } else {
    write.table(output, file = paste0(outputdir,"/",filename), sep=",",quote=FALSE,
                row.names = F)
  }

}









