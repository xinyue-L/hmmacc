#' The HMM-based Sleep/Wake Identification Algorithm to identify sustained inactive periods in the form of RData files.
#'
#' In the sleep/wake identification algorithm, the HMM assumes two states: sleep (state 1)
#' and wake (state 2). State 1 assumes zero-inflated truncated Gaussian, while state 2
#' assumes Gaussian.
#'
#' The input variables are the same as those in GGIR, but part of the functions are modified slightly to facilitate another stream of analysis.
#'
#' @import"GGIR"
#'
#' @param datadir Directory where the accelerometer files are stored or list of accelerometer filenames and directories
#' @param metadatadir Directory that holds a folders 'meta' and inside this a folder 'basic'.
#' @param f0 File index to start with (default = 1). Index refers to the filenames sorted in increasing order
#' @param f1 File index to finish with (defaults to number of files available)
#' @param idloc If value = 1 (default) the code assumes that ID number is stored in the obvious header field. If value = 2 the code uses the character string preceding the character '_' in the filename as the ID number
#' @param loglocation Location of the spreadsheet (csv) with sleep log information. The spreadsheet needs to have the following structure: one column for participant ID, and then followed by alternatingly one column for onset time and one column for waking time. There can be multiple sleeplogs in the same spreadsheet. The first raw of the spreadsheet needs to be filled with column names, it does not matter what these column names are. Timestamps are to be stored without date as in hh:mm:ss. If onset corresponds to lights out or intention to fall asleep, then it is the end-users responsibility to account for this in the interpretation of the results.
#' @param colid Column number in the sleep log spreadsheet in which the participant ID code is stored (default = 1)
#' @param coln1 Column number in the sleep log spreadsheet where the onset of the first night starts
#' @param nnights Number of nights for which sleep log information should be available. It assumes that this is constant within a study. If sleep log information is missing for certain nights then leave these blank
#' @param sleeplogidnum Should the participant identifier as stored in the sleeplog be interpretted as a number (TRUE=default) or a character (FALSE)?
#' @param do.visual If g.part4 is run with do.visual == TRUE then the function will generate a pdf with a visual representation of the overlap between the sleeplog entries and the accelerometer detections. This can be used to visualy verify that the sleeplog entries do not come with obvious mistakes.
#' @param outliers.only Relevant for do.visual == TRUE. Outliers.only == FALSE will visualise all available nights in the data. Outliers.only == TRUE will visualise only for nights with a difference in onset or waking time larger than the variable of argument criterror.
#' @param excludefirstlast If TRUE then the first and last night of the measurement are ignored for the sleep assessment.
#' @param criterror Relevant for do.visual == TRUE and outliers.only == TRUE. criterror specifies the number of minimum number of hours difference between sleep log and accelerometer estimate for the night to be included in the visualisation
#' @param includenightcrit Minimum number of valid hours per night (24 hour window between noon and noon)
#' @param relyonguider If TRUE then sleep onset and waking time are defined based on timestamps derived from the guider. If participants were instructed NOT to wear the accelerometer during waking hours then set to TRUE, in all other scenarios set to FALSE (default).
#' @param relyonsleeplog Now replaced by relyonguider. Values provided to argument relyonsleeplog will be passed on to argument relyonguider to not preserve functionality of old R script.
#' @param def.noc.sleep The time window during which sustained inactivity will be assumed to represent sleep, e.g. def.noc.sleep=c(21,9). This is only used if no sleep log entry is available. If def.noc.sleep is left blank 'def.noc.sleep=c()' then the 12 hour window centred at the least active 5 hours of the 24 hour period will be used instead. Here, L5 is hardcoded and will not change by changing argument winhr in function g.part2. If def.noc.sleep is filled with a single integer, e.g. def.noc.sleep=c(1) then the window will be detected with the method as described in van Hees et al. 2018 Scientific Reports.
#' @param storefolderstructure Store folder structure of the accelerometer data
#' @param overwrite Overwrite previously generated milestone data by this function for this particular dataset. If FALSE then it will skip the previously processed files (default = FALSE).
#' @param desiredtz See g.getmeta in GGIR
#' @param data_cleaning_file Optional path to a csv file you create that holds four columns: ID, day_part5, relyonguider_part4, and night_part4. ID should hold the participant ID. Columns day_part5 and night_part4 allow you to specify which day(s) and night(s) need to be excluded from part 5 and 4, respectively. So, this will be done regardless of whether the rest of GGIR thinks those day(s)/night(s) are valid. Column relyonguider_part4 allows you to specify for which nights part 4 should fully rely on the guider. See also package vignette.
#' @param excludefirst.part4 If TRUE then the first night of the measurement are ignored for the sleep assessment.
#' @param excludelast.part4 If TRUE then the last night of the measurement are ignored for the sleep assessment.
#'
#' @keywords sleep variables
#'
#' @return a RData file containing a data frame.
#' @return \item{nightsummary}{a data frame providing information on sleep but not directly usable in the downstream analysis. Run sleep.part5 for sleep variables and sleep.part6 for generating the corresponding sleep csv files.}
#'
#' @examples
#'
#'#metadir = "C:/myfolder/meta" # the meta folder containing the results from previous steps
#'#sleep.part4(metadatadir=metadir, overwrite=TRUE)
#'
#' @export
sleep.part4 <- function(datadir = c(), metadatadir = c(), f0=1, f1=0,
          idloc = 1, loglocation = c(), colid = 1, coln1 = 2, nnights = 7,
          sleeplogidnum = FALSE, do.visual = FALSE, outliers.only = FALSE,
          excludefirstlast = FALSE, criterror = 1, includenightcrit = 16,
          relyonguider = FALSE, relyonsleeplog = FALSE, def.noc.sleep = 1,
          storefolderstructure = FALSE, overwrite = FALSE, desiredtz = "",
          data_cleaning_file = c(), excludefirst.part4 = FALSE, excludelast.part4 = FALSE)
{
  if (exists("relyonsleeplog") == TRUE & exists("relyonguider") ==
      FALSE)
    relyonguider = relyonsleeplog
  nnpp = 40
  ms3.out = "/meta/ms3.out.new"
  if (file.exists(paste(metadatadir, ms3.out, sep = ""))) {
  } else {
    cat("Warning: First run g.part3 (mode = 3) before running g.part4 (mode = 4)")
  }
  ms4.out = "/meta/ms4.out.new"
  if (file.exists(paste(metadatadir, ms4.out, sep = ""))) {
  } else {
    dir.create(file.path(metadatadir, ms4.out))
  }
  meta.sleep.folder = paste(metadatadir, "/meta/ms3.out.new",
                            sep = "")
  if (length(loglocation) > 0) {
    dolog = TRUE
  } else {
    dolog = FALSE
  }
  if (dolog == TRUE) {
    LL = g.loadlog(loglocation, coln1, colid, nnights, sleeplogidnum)
    sleeplog = LL$sleeplog
    save(sleeplog, file = paste(metadatadir, "/meta/sleeplog.RData",
                                sep = ""))
  }
  fnames = dir(meta.sleep.folder)
  if (f1 > length(fnames))
    f1 = length(fnames)
  if (f0 > length(fnames))
    f0 = 1
  if (f1 == 0 | length(f1) == 0 | f1 > length(fnames))
    f1 = length(fnames)
  cnt = 1
  idlabels = rep(0, nnpp)
  pagei = 1
  cnt67 = 1
  colnamesnightsummary = c("ID", "night", "sleeponset",
                           "wakeup", "SptDuration", "sleepparam",
                           "guider_onset", "guider_wakeup", "guider_SptDuration",
                           "error_onset", "error_wake", "error_dur",
                           "fraction_night_invalid", "SleepDurationInSpt",
                           "duration_sib_wakinghours", "number_sib_sleepperiod",
                           "number_sib_wakinghours", "duration_sib_wakinghours_atleast15min",
                           "sleeponset_ts", "wakeup_ts", "guider_onset_ts",
                           "guider_wakeup_ts", "page", "daysleeper",
                           "weekday", "calendar_date", "filename",
                           "cleaningcode", "sleeplog_used", "acc_available",
                           "guider")
  if (storefolderstructure == TRUE) {
    colnamesnightsummary = c(colnamesnightsummary, "filename_dir",
                             "foldername")
  }
  logdur = rep(0, length(fnames))
  ffdone = c()
  ms4.out = "/meta/ms4.out.new"
  fnames.ms4 = dir(paste(metadatadir, ms4.out, sep = ""))
  fnames.ms4 = sort(fnames.ms4)
  ffdone = fnames.ms4
  fnames = sort(fnames)
  if (storefolderstructure == TRUE) {
    filelist = FALSE
    if (length(datadir) == 1) {
      if (length(unlist(strsplit(datadir, split = "[.](cs|bi|cw)"))) >
          1)
        filelist = TRUE
    }
    else {
      filelist = TRUE
    }
    if (filelist == FALSE) {
      fnamesfull = dir(datadir, recursive = TRUE, pattern = "[.](csv|bin|cwa)")
    }
    else {
      fnamesfull = datadir
    }
    f16 = function(X) {
      out = unlist(strsplit(X, "/"))
      f16 = out[length(out)]
    }
    f17 = function(X) {
      out = unlist(strsplit(X, "/"))
      f17 = out[(length(out) - 1)]
    }
    ffd = ffp = rep("", length(fnamesfull))
    if (length(fnamesfull) > 0) {
      fnamesshort = apply(X = as.matrix(fnamesfull), MARGIN = 1,
                          FUN = f16)
      foldername = apply(X = as.matrix(fnamesfull), MARGIN = 1,
                         FUN = f17)
      for (i in 1:length(fnames)) {
        ff = as.character(unlist(strsplit(fnames[i],
                                          ".RDa"))[1])
        if (length(which(fnamesshort == ff)) > 0) {
          ffd[i] = fnamesfull[which(fnamesshort == ff)]
          ffp[i] = foldername[which(fnamesshort == ff)]
        }
      }
    }
  }
  convertHRsinceprevMN2Clocktime = function(x) {
    if (x > 24)
      x = x - 24
    HR = floor(x)
    MI = floor((x - HR) * 60)
    SE = round(((x - HR) - (MI/60)) * 3600)
    if (SE == 60) {
      MI = MI + 1
      SE = 0
    }
    if (MI == 60) {
      HR = HR + 1
      MI = 0
    }
    if (HR == 24)
      HR = 0
    if (HR < 10)
      HR = paste0("0", HR)
    if (MI < 10)
      MI = paste0("0", MI)
    if (SE < 10)
      SE = paste0("0", SE)
    return(paste0(HR, ":", MI, ":", SE))
  }
  if (length(data_cleaning_file) > 0) {
    DaCleanFile = read.csv(data_cleaning_file)
  }
  for (i in f0:f1) {
    if (overwrite == TRUE) {
      skip = 0
    }
    else {
      skip = 0
      if (length(ffdone) > 0) {
        if (length(which(ffdone == fnames[i])) > 0)
          skip = 1
      }
    }
    if (skip == 0) {
      cat(paste(" ", i, sep = ""))
      if (cnt67 == 1) {
        if (do.visual == TRUE) {
          pdf(file = paste(metadatadir, "/results/visualisation_sleep.pdf",
                           sep = ""), width = 8.27, height = 11.69)
          par(mar = c(4, 5, 1, 2) + 0.1)
          plot(c(0, 0), c(1, 1), xlim = c(12, 36), ylim = c(0,
                                                            nnpp), col = "white", axes = FALSE,
               xlab = "time", ylab = "", main = paste("Page ",
                                                      pagei, sep = ""))
          axis(side = 1, at = 12:36, labels = c(12:24,
                                                1:12), cex.axis = 0.7)
          abline(v = c(18, 24, 30), lwd = 0.2, lty = 2)
          abline(v = c(15, 21, 27, 33), lwd = 0.2, lty = 3,
                 col = "grey")
        }
        cnt67 = 2
      }
      if (storefolderstructure == FALSE) {
        nightsummary = as.data.frame(matrix(0, 0, 31))
      }
      else {
        nightsummary = as.data.frame(matrix(0, 0, 33))
      }
      colnames(nightsummary) = colnamesnightsummary
      sumi = 1
      sptwindow_HDCZA_end = sptwindow_HDCZA_start = L5list = sib.cla.sum = c()
      load(paste(meta.sleep.folder, "/", fnames[i],
                 sep = ""))
      if (nrow(sib.cla.sum) != 0) {
        sib.cla.sum$sib.onset.time = iso8601chartime2POSIX(sib.cla.sum$sib.onset.time,
                                                           tz = desiredtz)
        sib.cla.sum$sib.end.time = iso8601chartime2POSIX(sib.cla.sum$sib.end.time,
                                                         tz = desiredtz)
        if (idloc == 2 | idloc == 5) {
          if (idloc == 2) {
            getCharBeforeUnderscore = function(x) {
              return(as.character(unlist(strsplit(x,
                                                  "_")))[1])
            }
          }
          else {
            getCharBeforeUnderscore = function(x) {
              return(as.character(unlist(strsplit(x,
                                                  " ")))[1])
            }
          }
          accid = apply(as.matrix(as.character(fnames[i])),
                        MARGIN = c(1), FUN = getCharBeforeUnderscore)
          accid_bu = accid
          getLastCharacterValue = function(x) {
            tmp = as.character(unlist(strsplit(x, "")))
            return(tmp[length(tmp)])
          }
          letter = apply(as.matrix(accid), MARGIN = c(1),
                         FUN = getLastCharacterValue)
          for (h in 1:length(accid)) {
            options(warn = -1)
            numletter = as.numeric(letter[h])
            options(warn = 0)
            if (is.na(numletter) == TRUE) {
              accid[h] = as.character(unlist(strsplit(accid[h],
                                                      letter[h]))[1])
            }
          }
          accid = suppressWarnings(as.numeric(accid))
          if (is.na(accid) == TRUE)
            accid = accid_bu
        }
        else {
          newaccid = fnames[i]
          if (length(unlist(strsplit(newaccid, "_"))) >
              1)
            newaccid = unlist(strsplit(newaccid, "_"))[1]
          if (length(unlist(strsplit(newaccid, " "))) >
              1)
            newaccid = unlist(strsplit(newaccid, " "))[1]
          if (length(unlist(strsplit(newaccid, "[.]RDa"))) >
              1)
            newaccid = unlist(strsplit(newaccid, "[.]RDa"))[1]
          if (length(unlist(strsplit(newaccid, "[.]cs"))) >
              1)
            newaccid = unlist(strsplit(newaccid, "[.]cs"))[1]
          accid = newaccid[1]
        }
        if (dolog == TRUE) {
          accid_num = suppressWarnings(as.numeric(accid))
          if (sleeplogidnum == FALSE) {
            wi = which(as.character(sleeplog$ID) == as.character(accid))
            if (length(wi) == 0) {
              wi_alternative = which(sleeplog$ID == accid_num)
              if (length(wi_alternative) > 0) {
                warning("\nArgument sleeplogidnum is set to FALSE, but it seems the identifiers are\n                    stored as numeric values, you may want to consider changing sleeplogidnum to TRUE")
              }
              else {
                warning(paste0("\nSleeplog id is stored as format: ",
                               as.character(sleeplog$ID[1]), ", while\n                           code expects format: ",
                               as.character(accid[1])))
              }
            }
          }
          else {
            wi = which(sleeplog$ID == accid_num)
            if (length(wi) == 0) {
              wi_alternative = which(as.character(sleeplog$ID) ==
                                       as.character(accid))
              if (length(wi_alternative) > 0) {
                warning("\nArgument sleeplogidnum is set to TRUE, but it seems the identifiers are\n                    stored as character values, you may want to consider changing sleeplogidnum to TRUE")
              }
              else {
                if (is.na(accid_num) == TRUE) {
                  warning(paste0("\nSleeplog id is stored as format: ",
                                 as.character(sleeplog$ID[1]), ", while\n                           code expects format: ",
                                 as.character(accid[1])))
                }
              }
            }
          }
        }
        else {
          wi = 1
        }
        if (length(nnights) == 0) {
          nnightlist = 1:max(sib.cla.sum$night)
        }
        else {
          if (max(sib.cla.sum$night) < nnights) {
            nnightlist = 1:max(sib.cla.sum$night)
          }
          else {
            nnightlist = 1:nnights
          }
        }
        if (length(nnightlist) < length(wi))
          nnightlist = nnightlist[1:length(wi)]
        nnights.list = nnightlist
        nnights.list = nnights.list[which(is.na(nnights.list) ==
                                            FALSE & nnights.list != 0)]
        if (excludefirstlast == TRUE & excludelast.part4 ==
            FALSE & excludefirst.part4 == FALSE) {
          if (length(nnights.list) >= 3) {
            nnights.list = nnights.list[2:(length(nnights.list) -
                                             1)]
          }
          else {
            nnights.list = c()
          }
        }
        else if (excludelast.part4 == FALSE & excludefirst.part4 ==
                 TRUE) {
          if (length(nnights.list) >= 2) {
            nnights.list = nnights.list[2:length(nnights.list)]
          }
          else {
            nnights.list = c()
          }
        }
        else if (excludelast.part4 == TRUE & excludefirst.part4 ==
                 FALSE) {
          if (length(nnights.list) >= 2) {
            nnights.list = nnights.list[1:(length(nnights.list) -
                                             1)]
          }
          else {
            nnights.list = c()
          }
        }
        calendar_date = wdayname = rep("", length(nnights.list))
        daysleeper = rep(FALSE, length(nnights.list))
        nightj = 1
        if (dolog == TRUE)
          sleeplog.t = sleeplog[wi, ]
        for (j in nnights.list) {
          if (length(def.noc.sleep) == 0 | length(sptwindow_HDCZA_start) ==
              0) {
            guider = "notavailable"
            if (length(L5list) > 0) {
              defaultSptOnset = L5list[j] - 6
              defaultSptWake = L5list[j] + 6
              guider = "L512"
            }
          }
          else if (length(def.noc.sleep) == 1 | length(loglocation) !=
                   0 & length(sptwindow_HDCZA_start) != 0) {
            defaultSptOnset = sptwindow_HDCZA_start[j]
            defaultSptWake = sptwindow_HDCZA_end[j]
            guider = "HDCZA"
            if (is.na(defaultSptOnset) == TRUE) {
              availableestimate = which(is.na(sptwindow_HDCZA_start) ==
                                          FALSE)
              cleaningcode = 6
              if (length(availableestimate) > 0) {
                defaultSptOnset = mean(sptwindow_HDCZA_start[availableestimate])
              }
              else {
                defaultSptOnset = L5list[j] - 6
                guider = "L512"
              }
            }
            if (is.na(defaultSptWake) == TRUE) {
              availableestimate = which(is.na(sptwindow_HDCZA_end) ==
                                          FALSE)
              cleaningcode = 6
              if (length(availableestimate) > 0) {
                defaultSptWake = mean(sptwindow_HDCZA_end[availableestimate])
              }
              else {
                defaultSptWake = L5list[j] + 6
                guider = "L512"
              }
            }
          }
          else if (length(def.noc.sleep) == 2) {
            defaultSptOnset = def.noc.sleep[1]
            defaultSptWake = def.noc.sleep[2]
            guider = "setwindow"
          }
          if (defaultSptOnset >= 24)
            defaultSptOnset = defaultSptOnset - 24
          if (defaultSptWake >= 24)
            defaultSptWake = defaultSptWake - 24
          defaultdur = defaultSptWake - defaultSptOnset
          sleeplog_used = FALSE
          if (dolog == TRUE) {
            if (is.na(sleeplog[wi[j], 3]) == FALSE) {
              sleeplog_used = TRUE
              cleaningcode = 0
              guider = "sleeplog"
            }
          }
          if (sleeplog_used == FALSE) {
            if (j == nnights.list[1] & dolog == FALSE) {
              sleeplog.t = data.frame(matrix(0, length(nnightlist),
                                             5))
              names(sleeplog.t) = c("ID", "night",
                                    "duration", "sleeponset",
                                    "sleepwake")
            }
            sleeplog.t[j, 1:5] = c(accid, j, defaultdur,
                                   convertHRsinceprevMN2Clocktime(defaultSptOnset),
                                   convertHRsinceprevMN2Clocktime(defaultSptWake))
            cleaningcode = 1
          }
          nightj = nightj + 1
          acc_available = TRUE
          spocum = as.data.frame(matrix(0, 0, 5))
          spocumi = 1
          sleeplog.t2 = sleeplog.t[which(sleeplog.t$night ==
                                           j), ]
          tmp1 = as.character(sleeplog.t2[, which(names(sleeplog.t2) ==
                                                    "sleeponset")])
          tmp2 = unlist(strsplit(tmp1, ":"))
          SptOnset = as.numeric(tmp2[1]) + (as.numeric(tmp2[2])/60) +
            (as.numeric(tmp2[3])/3600)
          tmp4 = as.character(sleeplog.t2[, which(names(sleeplog.t2) ==
                                                    "sleepwake")])
          tmp5 = unlist(strsplit(tmp4, ":"))
          SptWake = as.numeric(tmp5[1]) + (as.numeric(tmp5[2])/60) +
            (as.numeric(tmp5[3])/3600)
          daysleeper[j] = FALSE
          if (is.na(SptOnset) == FALSE & is.na(SptWake) ==
              FALSE & tmp1 != "" & tmp4 != "") {
            doubleDigitClocktime = function(x) {
              x = unlist(strsplit(x, ":"))
              xHR = as.numeric(x[1])
              xMI = as.numeric(x[2])
              xSE = as.numeric(x[3])
              if (xHR < 10)
                xHR = paste("0", xHR, sep = "")
              if (xMI < 10)
                xMI = paste("0", xMI, sep = "")
              if (xSE < 10)
                xSE = paste("0", xSE, sep = "")
              x = paste(xHR, ":", xMI, ":",
                        xSE, sep = "")
              return(x)
            }
            tmp1 = doubleDigitClocktime(tmp1)
            tmp4 = doubleDigitClocktime(tmp4)
            if (SptWake > 12 & SptOnset < 12)
              daysleeper[j] = TRUE
            if (SptWake > 12 & SptOnset > SptWake)
              daysleeper[j] = TRUE
            if (SptOnset < 12)
              SptOnset = SptOnset + 24
            if (SptWake <= 12)
              SptWake = SptWake + 24
            if (SptWake > 12 & SptWake < 18 & daysleeper[j] ==
                TRUE)
              SptWake = SptWake + 24
            if (daysleeper[j] == TRUE) {
              logdur[i] = SptOnset - SptWake
            }
            else {
              logdur[i] = SptWake - SptOnset
            }
          }
          else {
            SptOnset = defaultSptOnset
            SptWake = defaultSptWake + 24
            logdur[i] = SptWake - SptOnset
            cleaningcode = 1
          }
          if (excludefirstlast == FALSE) {
            if (daysleeper[j] == TRUE & j != max(nnights.list)) {
              loaddays = 2
            }
            else {
              loaddays = 1
            }
            if (daysleeper[j] == TRUE & j == max(nnights.list)) {
              daysleeper[j] = FALSE
              loaddays = 1
              if (SptWake > 36)
                SptWake = 36
              logdur[i] = SptWake - SptOnset
            }
          }
          else {
            if (daysleeper[j] == TRUE) {
              loaddays = 2
            }
            else {
              loaddays = 1
            }
          }
          dummyspo = matrix(0, 1, 5)
          dummyspo[1, 1] = 1
          spo_day = c()
          for (loaddaysi in 1:loaddays) {
            qq = sib.cla.sum
            sleepdet = qq[which(qq$night == (j + (loaddaysi -
                                                    1))), ]
            if (nrow(sleepdet) == 0) {
              if (spocumi == 1) {
                spocum = dummyspo
              }
              else {
                spocum = rbind(spocum, dummyspo)
              }
              spocumi = spocumi + 1
              cleaningcode = 3
              acc_available = FALSE
            }
            else {
              acc_available = TRUE
            }
            defs = unique(sleepdet$definition)
            for (k in defs) {
              ki = which(sleepdet$definition == k)
              sleepdet.t = sleepdet[ki, ]
              if (loaddaysi == 1)
                remember_fraction_invalid_day1 = sleepdet.t$fraction.night.invalid[1]
              nsp = length(unique(sleepdet.t$sib.period))
              spo = matrix(0, nsp, 5)
              if (nsp <= 1 & unique(sleepdet.t$sib.period)[1] ==
                  0) {
                spo[1, 1] = 1
                spo[1, 2:4] = 0
                spo[1, 5] = k
                if (daysleeper[j] == TRUE) {
                  tmpCmd = paste("spo_day", k,
                                 "= c()", sep = "")
                  eval(parse(text = tmpCmd))
                }
              }
              else {
                DD = g.create.sp.mat(nsp, spo, sleepdet.t,
                                     daysleep = daysleeper[j])
                if (loaddaysi == 1) {
                  wdayname[j] = DD$wdayname
                  calendar_date[j] = DD$calendar_date
                }
                spo = DD$spo
                reversetime2 = reversetime3 = c()
                if (daysleeper[j] == TRUE) {
                  if (loaddaysi == 1) {
                    w1 = which(spo[, 3] >= 18)
                    if (length(w1) > 0) {
                      spo = as.matrix(spo[w1, ])
                      if (ncol(spo) == 1)
                        spo = t(spo)
                      if (nrow(spo) == 1) {
                        if (spo[1, 2] <= 18)
                          spo[1, 2] = 18
                      }
                      else {
                        spo[which(spo[, 2] <= 18), 2] = 18
                      }
                      tmpCmd = paste("spo_day",
                                     k, "= spo", sep = "")
                      eval(parse(text = tmpCmd))
                    }
                    else {
                      tmpCmd = paste("spo_day",
                                     k, "= c()", sep = "")
                      eval(parse(text = tmpCmd))
                    }
                  }
                  else if (loaddaysi == 2 & length(eval(parse(text = paste0("spo_day",
                                                                            k)))) > 0) {
                    w2 = which(spo[, 2] < 18)
                    if (length(w2) > 0) {
                      spo = as.matrix(spo[w2, ])
                      if (ncol(spo) == 1)
                        spo = t(spo)
                      if (nrow(spo) == 1) {
                        if (spo[1, 3] > 18)
                          spo[1, 3] = 18
                      }
                      else {
                        spo[which(spo[, 3] > 18), 3] = 18
                      }
                      spo[, 2:3] = spo[, 2:3] + 24
                      tmpCmd = paste("spo_day2",
                                     k, "= spo", sep = "")
                      eval(parse(text = tmpCmd))
                    }
                    else {
                      tmpCmd = paste("spo_day2",
                                     k, "= c()", sep = "")
                      eval(parse(text = tmpCmd))
                    }
                    name1 = paste("spo_day", k,
                                  sep = "")
                    name2 = paste("spo_day2", k,
                                  sep = "")
                    tmpCmd = paste("spo = rbind(",
                                   name1, ",", name2, ")",
                                   sep = "")
                    eval(parse(text = tmpCmd))
                  }
                }
                if (daysleeper[j] == TRUE) {
                  if (SptWake < 21 & SptWake > 12 & SptOnset >
                      SptWake)
                    SptWake = SptWake + 24
                }
                relyonguider_thisnight = FALSE
                if (length(data_cleaning_file) > 0) {
                  if (length(which(DaCleanFile$relyonguider_part4 ==
                                   j & DaCleanFile$ID == accid)) > 0) {
                    relyonguider_thisnight = TRUE
                  }
                }
                if (length(which(spo[, 2] < SptWake &
                                 spo[, 3] > SptOnset)) == 0 | relyonguider_thisnight ==
                    TRUE) {
                  cleaningcode = 5
                  newlines = rbind(spo[1, ], spo[1, ])
                  newlines[1, 1:4] = c(nrow(spo) + 1,
                                       SptOnset, SptOnset + 1/60, 1)
                  newlines[2, 1:4] = c(nrow(spo) + 1,
                                       SptWake - 1/60, SptWake, 1)
                  spo = rbind(spo, newlines)
                  spo = spo[order(spo[, 2]), ]
                  spo[, 1] = 1:nrow(spo)
                  relyonguider_thisnight = TRUE
                }
                for (evi in 1:nrow(spo)) {
                  if (spo[evi, 2] < SptWake & spo[evi,
                                                  3] > SptOnset) {
                    spo[evi, 4] = 1
                    if (relyonguider == TRUE | relyonguider_thisnight ==
                        TRUE) {
                      if ((spo[evi, 2] < SptWake & spo[evi,
                                                       3] > SptWake) | (spo[evi, 2] <
                                                                        SptWake & spo[evi, 3] < spo[evi,
                                                                                                    2])) {
                        spo[evi, 3] = SptWake
                      }
                      if ((spo[evi, 2] < SptOnset & spo[evi,
                                                        3] > SptOnset) | (spo[evi, 3] >
                                                                          SptOnset & spo[evi, 3] < spo[evi,
                                                                                                       2])) {
                        spo[evi, 2] = SptOnset
                      }
                    }
                  }
                }
                if (daysleeper[j] == TRUE) {
                  reversetime2 = which(spo[, 2] >= 36)
                  reversetime3 = which(spo[, 3] >= 36)
                  if (length(reversetime2) > 0)
                    spo[reversetime2, 2] = spo[reversetime2,
                                               2] - 24
                  if (length(reversetime3) > 0)
                    spo[reversetime3, 3] = spo[reversetime3,
                                               3] - 24
                }
                spo[, 5] = k
                if (spocumi == 1) {
                  spocum = spo
                }
                else {
                  spocum = rbind(spocum, spo)
                }
                spocumi = spocumi + 1
              }
            }
          }
          if (do.visual == TRUE) {
            if (cnt == (nnpp + 1)) {
              cat(" NEW ")
              pagei = pagei + 1
              if (length(idlabels) < nnpp) {
                idlabels = c(idlabels, rep(" ",
                                           length(idlabels) - nnpp))
              }
              else if (length(idlabels) > nnpp) {
                idlabels = idlabels[1:nnpp]
              }
              axis(side = 2, at = 1:nnpp, labels = idlabels,
                   las = 1, cex.axis = 0.6)
              idlabels = rep(0, nnpp)
              plot(c(0, 0), c(1, 1), xlim = c(12, 36),
                   ylim = c(0, nnpp), col = "white",
                   axes = FALSE, xlab = "time", ylab = "",
                   main = paste("Page ", pagei, sep = ""))
              axis(side = 1, at = 12:36, labels = c(12:24,
                                                    1:12), cex.axis = 0.7)
              abline(v = c(18, 24, 30), lwd = 0.2, lty = 2)
              abline(v = c(15, 21, 27, 33), lwd = 0.2,
                     lty = 3, col = "grey")
              cnt = 1
            }
          }
          if (length(spocum) > 0) {
            if (length(which(spocum[, 5] == "0")) >
                0) {
              spocum = spocum[-which(spocum[, 5] == "0"),
              ]
            }
          }
          if (length(spocum) > 0 & class(spocum)[1] ==
              "matrix" & length(calendar_date) >=
              j) {
            if (nrow(spocum) > 1 & ncol(spocum) >= 5 &
                calendar_date[j] != "") {
              undef = unique(spocum[, 5])
              for (defi in undef) {
                rowswithdefi = which(spocum[, 5] == defi)
                if (length(rowswithdefi) > 1) {
                  spocum.t = spocum[rowswithdefi, ]
                  correct01010pattern = function(x) {
                    x = as.numeric(x)
                    if (length(which(diff(x) == 1)) >
                        1) {
                      minone = which(diff(x) == -1) +
                        1
                      plusone = which(diff(x) == 1)
                      matchingvalue = which(minone %in%
                                              plusone == TRUE)
                      if (length(matchingvalue) > 0)
                        x[minone[matchingvalue]] = 1
                    }
                    return(x)
                  }
                  delta_t1 = diff(as.numeric(spocum.t[,
                                                      3]))
                  spocum.t[, 4] = correct01010pattern(spocum.t[,
                                                               4])
                  nightsummary[sumi, 1] = accid
                  nightsummary[sumi, 2] = j
                  if (is.matrix(spocum.t) == FALSE) {
                    spocum.t = t(as.matrix(spocum.t))
                  }
                  spocum.t = spocum.t[!duplicated(spocum.t),
                  ]
                  if (is.matrix(spocum.t) == FALSE)
                    spocum.t = as.matrix(spocum.t)
                  if (ncol(spocum.t) < 4 & nrow(spocum.t) >
                      3)
                    spocum.t = t(spocum.t)
                  if (length(which(as.numeric(spocum.t[,
                                                       4]) == 1)) > 0) {
                    rtl = which(spocum.t[, 4] == 1)
                    nightsummary[sumi, 3] = spocum.t[rtl[1],
                                                     2]
                    nightsummary[sumi, 4] = spocum.t[rtl[length(rtl)],
                                                     3]
                  }
                  else {
                    cleaningcode = 5
                    nightsummary[sumi, 3] = SptOnset
                    nightsummary[sumi, 4] = SptWake
                  }
                  nightsummary[, 3] = as.numeric(nightsummary[,
                                                              3])
                  nightsummary[, 4] = as.numeric(nightsummary[,
                                                              4])
                  if (nightsummary[sumi, 3] > nightsummary[sumi,
                                                           4] & nightsummary[sumi, 4] < 36 &
                      daysleeper[j] == TRUE) {
                    nightsummary[sumi, 4] = nightsummary[sumi,
                                                         4] + 24
                  }
                  if (nightsummary[sumi, 3] == nightsummary[sumi,
                                                            4] & nightsummary[sumi, 4] == 18) {
                    nightsummary[sumi, 4] = nightsummary[sumi,
                                                         4] + 24
                  }
                  nightsummary[sumi, 5] = nightsummary[sumi,
                                                       4] - nightsummary[sumi, 3]
                  nightsummary[, 5] = as.numeric(nightsummary[,
                                                              5])
                  nightsummary[sumi, 6] = defi
                  if (SptOnset > 36) {
                    nightsummary[sumi, 7] = SptOnset -
                      24
                  }
                  else {
                    nightsummary[sumi, 7] = SptOnset
                  }
                  if (SptWake > 36 & daysleeper[j] ==
                      FALSE) {
                    nightsummary[sumi, 8] = SptWake -
                      24
                  }
                  else {
                    nightsummary[sumi, 8] = SptWake
                  }
                  if (nightsummary[sumi, 7] > nightsummary[sumi,
                                                           8]) {
                    nightsummary[sumi, 9] = abs((36 -
                                                   nightsummary[sumi, 7]) + (nightsummary[sumi,
                                                                                          8] - 12))
                  }
                  else {
                    nightsummary[sumi, 9] = abs(nightsummary[sumi,
                                                             8] - nightsummary[sumi, 7])
                  }
                  nightsummary[sumi, 10] = nightsummary[sumi,
                                                        3] - nightsummary[sumi, 7]
                  nightsummary[sumi, 11] = nightsummary[sumi,
                                                        4] - nightsummary[sumi, 8]
                  if (nightsummary[sumi, 10] > 12)
                    nightsummary[sumi, 10] = -(24 - nightsummary[sumi,
                                                                 10])
                  if (nightsummary[sumi, 10] < -12)
                    nightsummary[sumi, 10] = -(nightsummary[sumi,
                                                            10] + 24)
                  if (nightsummary[sumi, 11] > 12)
                    nightsummary[sumi, 11] = -(24 - nightsummary[sumi,
                                                                 11])
                  if (nightsummary[sumi, 11] < -12)
                    nightsummary[sumi, 11] = -(nightsummary[sumi,
                                                            11] + 24)
                  nightsummary[sumi, 12] = nightsummary[sumi,
                                                        5] - nightsummary[sumi, 9]
                  if (acc_available == TRUE) {
                    nightsummary[sumi, 13] = remember_fraction_invalid_day1
                    if (remember_fraction_invalid_day1 >
                        ((24 - includenightcrit)/24)) {
                      cleaningcode = 2
                    }
                  }
                  else {
                    nightsummary[sumi, 13] = 1
                  }
                  nocs = as.numeric(spocum.t[which(spocum.t[,
                                                            4] == 1), 3]) - as.numeric(spocum.t[which(spocum.t[,
                                                                                                               4] == 1), 2])
                  sibds = as.numeric(spocum.t[which(spocum.t[,
                                                             4] == 0), 3]) - as.numeric(spocum.t[which(spocum.t[,
                                                                                                                4] == 0), 2])
                  negval = which(nocs < 0)
                  if (length(negval) > 0) {
                    kk0 = as.numeric(spocum.t[which(spocum.t[,
                                                             4] == 1), 2])
                    kk1 = as.numeric(spocum.t[which(spocum.t[,
                                                             4] == 1), 3])
                    kk1[negval] = kk1[negval] + 1
                    nocs = kk1 - kk0
                  }
                  if (length(nocs) > 0) {
                    spocum.t.dur.noc = sum(nocs)
                  }
                  else {
                    spocum.t.dur.noc = 0
                  }
                  is_this_a_dst_night_output = is_this_a_dst_night(calendar_date = calendar_date[j],
                                                                   tz = desiredtz)
                  dst_night_or_not = is_this_a_dst_night_output$dst_night_or_not
                  dsthour = is_this_a_dst_night_output$dsthour
                  if (dst_night_or_not == 1) {
                    checkoverlap = spocum.t[which(spocum.t[,
                                                           4] == 1), 2:3]
                    if (length(checkoverlap) > 0 & is.matrix(checkoverlap) ==
                        TRUE) {
                      overlaps = which(checkoverlap[,
                                                    1] <= (dsthour + 24) & checkoverlap[,
                                                                                        2] >= (dsthour + 25))
                    }
                    else {
                      overlaps = c()
                    }
                    if (length(overlaps) > 0) {
                      spocum.t.dur.noc = spocum.t.dur.noc -
                        1
                      nightsummary[sumi, 5] = nightsummary[sumi,
                                                           5] - 1
                      nightsummary[sumi, 9] = nightsummary[sumi,
                                                           9] - 1
                    }
                  }
                  else if (dst_night_or_not == -1) {
                    if (nightsummary[sumi, 3] <= (dsthour +
                                                  24) & nightsummary[sumi, 4] >=
                        (dsthour + 25)) {
                      nightsummary[sumi, 5] = nightsummary[sumi,
                                                           5] + 1
                    }
                    if (nightsummary[sumi, 7] <= (dsthour +
                                                  24) & nightsummary[sumi, 8] >=
                        (dsthour + 25)) {
                      nightsummary[sumi, 9] = nightsummary[sumi,
                                                           9] + 1
                    }
                    correctSptEdgingInDoubleHour = function(nightsummary,
                                                            onsetcol, wakecol, durcol, dsthour,
                                                            delta_t1) {
                      wakeInDoubleHour = nightsummary[,
                                                      wakecol] >= (dsthour + 24) &
                        nightsummary[, wakecol] <= (dsthour +
                                                      25)
                      onsetInDoubleHour = nightsummary[,
                                                       onsetcol] >= (dsthour + 24) &
                        nightsummary[, onsetcol] <= (dsthour +
                                                       25)
                      onsetBeforeDoubleHour = nightsummary[,
                                                           onsetcol] <= (dsthour + 24)
                      wakeAfterDoubleHour = nightsummary[,
                                                         wakecol] >= (dsthour + 25)
                      timeWentBackward = length(which(delta_t1 <
                                                        0)) > 0
                      if (onsetBeforeDoubleHour == TRUE &
                          wakeInDoubleHour == TRUE) {
                        if (timeWentBackward == TRUE) {
                          nightsummary[, durcol] = nightsummary[,
                                                                durcol] + 1
                        }
                        else if (timeWentBackward ==
                                 FALSE) {
                          nightsummary[, durcol] = nightsummary[,
                                                                durcol] + 1
                        }
                      }
                      if (wakeAfterDoubleHour == TRUE &
                          onsetInDoubleHour == TRUE) {
                        if (timeWentBackward == TRUE) {
                          nightsummary[, durcol] = nightsummary[,
                                                                durcol] + 1
                        }
                        else if (timeWentBackward ==
                                 FALSE) {
                          nightsummary[, durcol] = nightsummary[,
                                                                durcol] + 1
                        }
                      }
                      return(nightsummary)
                    }
                    nightsummary[sumi, ] = correctSptEdgingInDoubleHour(nightsummary[sumi,
                    ], onsetcol = 3, wakecol = 4, durcol = 5,
                    dsthour = dsthour, delta_t1 = delta_t1)
                    nightsummary[sumi, ] = correctSptEdgingInDoubleHour(nightsummary[sumi,
                    ], onsetcol = 7, wakecol = 8, durcol = 9,
                    dsthour = dsthour, delta_t1 = delta_t1)
                  }
                  sibds_atleast15min = 0
                  if (length(sibds) > 0) {
                    spocum.t.dur_sibd = sum(sibds)
                    atleast15min = which(sibds >= 1/4)
                    if (length(atleast15min) > 0) {
                      sibds_atleast15min = sibds[atleast15min]
                      spocum.t.dur_sibd_atleast15min = sum(sibds_atleast15min)
                    }
                    else {
                      spocum.t.dur_sibd_atleast15min = 0
                    }
                  }
                  else {
                    spocum.t.dur_sibd = 0
                    spocum.t.dur_sibd_atleast15min = 0
                  }
                  nightsummary[sumi, 14] = spocum.t.dur.noc
                  nightsummary[sumi, 15] = spocum.t.dur_sibd
                  nightsummary[sumi, 16] = length(which(spocum.t[,
                                                                 4] == 1))
                  nightsummary[sumi, 17] = length(which(spocum.t[,
                                                                 4] == 0))
                  nightsummary[sumi, 18] = as.numeric(spocum.t.dur_sibd_atleast15min)
                  acc_onset = nightsummary[sumi, 3]
                  acc_wake = nightsummary[sumi, 4]
                  if (acc_onset > 24)
                    acc_onset = acc_onset - 24
                  if (acc_wake > 24)
                    acc_wake = acc_wake - 24
                  acc_onsetTS = convertHRsinceprevMN2Clocktime(acc_onset)
                  acc_wakeTS = convertHRsinceprevMN2Clocktime(acc_wake)
                  nightsummary[sumi, 19] = acc_onsetTS
                  nightsummary[sumi, 20] = acc_wakeTS
                  nightsummary[sumi, 21] = tmp1
                  nightsummary[sumi, 22] = tmp4
                  nightsummary[sumi, 23] = pagei
                  nightsummary[sumi, 24] = daysleeper[j]
                  nightsummary[sumi, 25] = wdayname[j]
                  nightsummary[sumi, 26] = calendar_date[j]
                  nightsummary[sumi, 27] = fnames[i]
                  if (do.visual == TRUE) {
                    if (defi == undef[1]) {
                      if (outliers.only == TRUE) {
                        if (abs(nightsummary$error_onset[sumi]) >
                            criterror | abs(nightsummary$error_wake[sumi]) >
                            criterror | abs(nightsummary$error_dur[sumi]) >
                            (criterror * 2)) {
                          doplot = TRUE
                          cat(" PLOT ")
                        }
                        else {
                          doplot = FALSE
                        }
                      }
                      else {
                        doplot = TRUE
                      }
                    }
                    if (length(loglocation) > 0) {
                      cleaningcriterion = 1
                    }
                    else {
                      cleaningcriterion = 2
                    }
                    if (doplot == TRUE & cleaningcode <
                        cleaningcriterion) {
                      idlabels[cnt] = paste("ID",
                                            accid, " night", j, sep = "")
                      den = 20
                      defii = which(undef == defi)
                      qtop = ((defii/length(undef)) *
                                0.6) - 0.3
                      qbot = (((defii - 1)/length(undef)) *
                                0.6) - 0.3
                      for (pli in 1:nrow(spocum.t)) {
                        if (spocum.t[pli, 2] > spocum.t[pli,
                                                        3]) {
                          if (pli > 1 & pli < nrow(spocum.t) &
                              abs(as.numeric(spocum.t[pli,
                                                      2]) - as.numeric(spocum.t[pli,
                                                                                3])) < 2) {
                            spocum.t[pli, 2:3] = spocum.t[pli,
                                                          3:2]
                          }
                        }
                        if (spocum.t[pli, 4] == 1) {
                          colb = rainbow(length(undef),
                                         start = 0.7, end = 1)
                        }
                        else {
                          colb = rainbow(length(undef),
                                         start = 0.2, end = 0.4)
                        }
                        if (spocum.t[pli, 2] > spocum.t[pli,
                                                        3]) {
                          rect(xleft = spocum.t[pli,
                                                2], ybottom = (cnt + qbot),
                               xright = 36, ytop = (cnt +
                                                      qtop), col = colb[defii],
                               border = NA)
                          rect(xleft = 12, ybottom = (cnt +
                                                        qbot), xright = spocum.t[pli,
                                                                                 3], ytop = (cnt + qtop),
                               col = colb[defii], border = NA)
                        }
                        else {
                          rect(xleft = spocum.t[pli,
                                                2], ybottom = (cnt + qbot),
                               xright = spocum.t[pli, 3],
                               ytop = (cnt + qtop), col = colb[defii],
                               border = NA)
                        }
                      }
                      SptWaken = SptWake
                      SptOnsetn = SptOnset
                      if (SptWake > 36)
                        SptWaken = SptWake - 24
                      if (SptOnset > 36)
                        SptOnsetn = SptOnset - 24
                      if (defi == undef[length(undef)]) {
                        if (SptOnsetn > SptWaken) {
                          rect(xleft = SptOnsetn, ybottom = (cnt -
                                                               0.3), xright = 36, ytop = (cnt +
                                                                                            0.3), col = "black",
                               border = TRUE, density = den)
                          rect(xleft = 12, ybottom = (cnt -
                                                        0.3), xright = SptWaken,
                               ytop = (cnt + 0.3), col = "black",
                               border = TRUE, density = den)
                        }
                        else {
                          rect(xleft = SptOnsetn, ybottom = (cnt -
                                                               0.3), xright = SptWaken,
                               ytop = (cnt + 0.3), col = "black",
                               border = TRUE, density = den)
                        }
                      }
                    }
                  }
                  nightsummary[sumi, 28] = cleaningcode
                  nightsummary[sumi, 29] = sleeplog_used
                  nightsummary[sumi, 30] = acc_available
                  nightsummary[sumi, 31] = guider
                  if (storefolderstructure == TRUE) {
                    nightsummary[sumi, 32] = ffd[i]
                    nightsummary[sumi, 33] = ffp[i]
                  }
                  sumi = sumi + 1
                }
                if (do.visual == TRUE) {
                  if (cleaningcode < cleaningcriterion &
                      doplot == TRUE) {
                    lines(x = c(12, 36), y = c(cnt, cnt),
                          lwd = 0.2, lty = 2)
                    if (daysleeper[j] == TRUE) {
                      lines(x = c(18, 18), y = c((cnt -
                                                    0.3), (cnt + 0.3)), lwd = 2,
                            lty = 2, col = "black")
                    }
                    cnt = cnt + 1
                  }
                }
              }
            }
          }
        }
        if (length(nnights.list) == 0) {
          nightsummary[sumi, 1:2] = c(accid, 0)
          nightsummary[sumi, 3:26] = NA
          nightsummary[sumi, 27] = fnames[i]
          nightsummary[sumi, 28] = 4
          nightsummary[sumi, 29:31] = c(FALSE, TRUE,
                                        "NA")
          if (storefolderstructure == TRUE) {
            nightsummary[sumi, 32:33] = c(ffd[i], ffp[i])
          }
          sumi = sumi + 1
        }
        save(nightsummary, file = paste(metadatadir,
                                        ms4.out, "/", fnames[i], sep = ""))
      }
    }
  }
  if (cnt67 == 2 & do.visual == TRUE) {
    if (cnt - 1 != (nnpp + 1)) {
      zerolabel = which(idlabels == 0)
      if (length(zerolabel) > 0)
        idlabels[zerolabel] = " "
      axis(side = 2, at = 1:nnpp, labels = idlabels, las = 1,
           cex.axis = 0.5)
    }
    dev.off()
    cnt67 = 1
  }
}
