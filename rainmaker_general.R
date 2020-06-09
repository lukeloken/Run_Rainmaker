# General workflow for Rainmaker
# Last update: 5/15/2020 RBC

# FYI - The Rainmaker workflow follows this path:
# - 1 Aquire and prepare precipitation data for use in `Rainmaker` from either dataRetrieval/NWIS or Aquarius exported rainfall data
# - 2 Optional - Aquire storm start and end times for `RMevents_sample`
# - 3 Determine precipitation event start and end times using `RMevents` or `RMevents_sample`
# - 4 Compute intensities using `RMintensity`
# - 5 Compute erosivity index using `RMerosivity`
# - 6 Compute antecedent rainfall using `RMarf`
# - 7 Output the results to a file

# This script is written for the user to change a few settings at the top, 
# then run the remainder of the script without needing to make any further changes.

#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***
#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***
#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***

#***# SETTINGS/PREFERENCES  #***#***#***#

    #***# RAINFALL DATA #***#***#***#***#
    # Choose if you will use dataRetieval to get data from NWIS, or load a .csv exported from Aquarius (AQ)
    # Opt 1) DataRetrieval
        rain_site <- '424008089050201' # <-- insert station ID
        # Set dates that bound what data to retrieve:
        start_date <- '2015-09-30' # YYYY-MM-DD date/time of study start
        end_date <- '2018-09-30' # YYYY-MM-DD date/time of study end
    
    # Opt 2) Export from Aquarius and read in the .csv
        # Name of the unit-value rainfall data exported from AQ
        rain.filename <- NA #'AQ_precip_405051083391201.csv' # <-- insert your file name here --- MUST BE EXPORTED FROM AQ TO WORK IN THIS SCRIPT

    #***# MANUAL STORM EVENTS  #***#***#***#***#
    # Choose if you will provide storm start and end times that are saved in a .csv
        sampled_events <- 'no' # yes indicates you are providing sample start and end times otherwise say 'no'
    
        sample.filename <- NA #"test_405051083391201.csv" # <-- Insert file name with storm start and end times here.
        # Format headers to be 'storm_start' and 'storm_end'
        # Can fill in with NA if not providing sample times  
        
    #***# DISCHARGE-COMPUTED STORM EVENTS  #***#***#***#***#
    # Choose if you want RainMaker to determine event start and end times based on the discharge record.
    # This option is for sites that go to zero between events, or have a consistent lower discharge threshold.
    # Not recommended for sites with varying baseflow.
        computed_events <- 'yes' # yes indicates you want to run discharge_events function to find sample start and end times otherwise say 'no'
        Q_site <- rain_site # if discharge site is different from rain gage site, otherwise use: <-rain_site    
        Q_ieHr <- 6 # adjust hydrograph interevent period, default is 6 hours.
        Q_thresh <- 0 # discharge greater than this value triggers an event.
        
    #***# FILE PATHS ***#***#***#***#***#
    # INPUT - folder with storm start and end times and/or unit-value rainfall data exported from AQ, if applicable:
        input_path <- NA #"C:/Users/rbcarvin/OneDrive - DOI/Rainmaker" # <-- NA if only using dataRetireval
    # OUTPUT  
        output_path <- 'C:/Users/rbcarvin/OneDrive - DOI/Rainmaker' # <-- where to save output file
        output <- "StormSummary.csv" # <-- insert output filename here
        
    #***# RAINMAKER OPTIONS #***#***#***#***#    
    # Time zone
        site_tz <- "Etc/GMT+6"
        AQ_datetime_format <- "%Y-%m-%d %H:%M"
        Smpl_datetime_format <- "%m/%d/%Y %H:%M"
    # Time between events in hours (interevent)
        ieHr <- 2
    # Amount it must rain to count as an event in tenths of inches (rain threshold) -- note RMevents_sample does not use rainthresh
        rainthresh <- 0.008
    # Antecedent Rainfall in days (ARF.days)
        antecedentDays = c(0.5, 1, 7, 14)
    
#***# SELECT ALL BELOW THIS LINE AND RUN ***#***#***#          
    # keyboard shortcut: Ctrl+Alt+E
        
#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***
#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***
#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***
# Check for Rainmaker updates
devtools::install_github("USGS-R/Rainmaker")
        
# Turn on Rainmaker package
library("Rainmaker")
options(scipen = 999)

# - 1 Aquire data
  # get data from NWIS using dataRetrieval
  if (is.na(rain.filename)) {
          require(dataRetrieval)
          # set readNWISuv parameters
          parameterCd <- "00045"  # Precipitation
          startDate <- as.Date(start_date) - 15 # put a week buffer on the study start date in case storm started prior to first sample date
          endDate <- as.Date(end_date)
          
          # get NWIS data
          message('Pulling precip data from NWIS.')
          start_time <- Sys.time()
          Rain.uv <- readNWISuv(rain_site, parameterCd, tz = site_tz)
          end_time <- Sys.time()
          if (nrow(Rain.uv) > 0){
            message(paste0(nrow(Rain.uv)), ' rows of data pulled from NWIS in ', round(difftime(end_time , start_time, units = 'secs'), 0), ' seconds.')
          } else {
            stop('No precip data pulled from NWIS. Please check inputs to verify correct site number and start and end dates. To debug, see code in "data_processing/2_calc_rain_variables.R"')
          }
          
          # rename columns
          Rain.uv <- renameNWISColumns(Rain.uv)
          
          # rename columns
          names(Rain.uv)[grep('precip_inst$', names(Rain.uv), ignore.case = TRUE)] <- 'rain'
          names(Rain.uv)[grep('dateTime', names(Rain.uv), ignore.case = TRUE)] <- 'pdate'
          
          # print warning if dates of rain do not span dates of study
          if (min(as.Date(Rain.uv$pdate)) > start_date) {
            warning('Data pulled from NWIS does not span the entire study period.')
          }
          
          if (max(as.Date(Rain.uv$pdate)) < end_date) {
            warning('Data pulled from NWIS does not span the entire study period.')
          }
          
        } else {
          
  # OR read in AQ precip file
          Rain.uv <- read.csv(file=paste0(input_path,"/",rain.filename),skip=15,col.names = c("UTC", "dateTime", "rain", "x", "y", "z"))
          Rain.uv[,'pdate'] <- as.POSIXct(Rain.uv[,'dateTime'], tz = site_tz, format = AQ_datetime_format)
          #names(Rain.uv) <- c('pdate', 'rain')
        }
        
# - 2 Sample data
        if(sampled_events=='yes'){
          wq.dat <- read.csv(file=paste0(input_path,"/",sample.filename))
          wq.dat[,'storm_start'] <- as.POSIXct(wq.dat[,'storm_start'],tz= site_tz, format = Smpl_datetime_format)
          wq.dat[,'storm_end'] <- as.POSIXct(wq.dat[,'storm_end'],tz= site_tz, format = Smpl_datetime_format)
        }
        if(computed_events=='yes'){
          # Q dataRetrieval:
          Q_parameterCd <- "00060"  # Discharge
          startDate <- as.Date(start_date) #removes 15 day buffer that was applied to rain data retrieval
          message('Pulling discharge data from NWIS.')
          start_time <- Sys.time()
          Q_raw <- readNWISuv(Q_site, Q_parameterCd, startDate, endDate, tz = site_tz)
          end_time <- Sys.time()
          if (nrow(Q_raw) > 0){
            message(paste0(nrow(Q_raw)), ' rows of data pulled from NWIS in ', round(difftime(end_time , start_time, units = 'secs'), 0), ' seconds.')
          } else {
            stop('No precip data pulled from NWIS. Please check inputs to verify correct site number and start and end dates. To debug, see code in "data_processing/2_calc_rain_variables.R"')
          }
          Q_raw <- renameNWISColumns(Q_raw)
          Q_prep <- Q_raw
          Q_event <- discharge_events(df=Q_prep, ieHr = Q_ieHr, qthresh = Q_thresh, discharge = "Flow_Inst",
                                      time = "dateTime", ieHr_check = TRUE)
        }
# - 3 Find rain depth per event
        if(sampled_events=='yes') {
          events <- RMevents_sample(df = Rain.uv, ieHr = ieHr, rain = 'rain', time = 'pdate', 
                                    dfsamples = wq.dat, bdate = 'storm_start', edate = 'storm_end')
          # extract data from events output
          Rain.event.list <- events$storms
          tipsbystorm <- events$tipsbystorm
        } else if (computed_events=='yes') {
          events <- RMevents_sample(df = Rain.uv, ieHr = ieHr, rain = 'rain', time = 'pdate', 
                                    dfsamples = Q_event, bdate = 'StartDate', edate = 'EndDate')
          # extract data from events output
          Rain.event.list <- events$storms
          tipsbystorm <- events$tipsbystorm
        } else {
          # Determine precipitation event start and end times        
          Rain.events <- RMevents(df=Rain.uv, ieHr=ieHr, rainthresh=rainthresh, rain="rain", time="pdate")
          # extract data from events output
          Rain.event.list <- Rain.events$storms2
          tipsbystorm <- Rain.events$tipsbystorm
        }
          
# - 4 Compute intensities
        StormSummary <- RMintensity(df=tipsbystorm, date="pdate", df.events=Rain.event.list, depth="rain", xmin=c(5,10,15,30,60))

# - 5 Compute erosivity index
        # method 1
        StormSummary <- RMerosivity(df=tipsbystorm, ieHr=ieHr, rain="rain", StormSummary=StormSummary, method=1)
        StormSummary <- rename(StormSummary, 'erosivity_m1' = "erosivity", 'energy_m1' = 'energy')
        # method 2
        StormSummary <- RMerosivity(df= tipsbystorm, ieHr=ieHr, rain="rain", StormSummary=StormSummary, method=2)
        StormSummary <- rename(StormSummary, 'erosivity_m2' = "erosivity", 'energy_m2' = 'energy')
# - 6 Compute antecedent rainfall using `RMarf`
        StormSummary <- RMarf(df = Rain.uv, date = 'pdate', df.events = StormSummary,
                              days = antecedentDays, varnameout = "ARFdays")
# - 7 Output the results to a file
        write.csv(StormSummary,file=output,row.names = FALSE)
#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***#***
