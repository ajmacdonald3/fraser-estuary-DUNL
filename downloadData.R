#
# INSTRUCTIONS:
#
# 1) Make a file called 'pass' in your working directory (no file extension).
# 2) Open this file and put in your password and nothing else.
# 3) Save and close the file.
# 4) Type your Motus login in line 18 of this document
#

downloadMotus_recvData <- function(site.name = '', serno = '', projID = '', dir = getwd()) {

  dir <- ifelse(identical(dir, ''), getwd(), gsub('/$', '', dir))


  if (identical(Sys.getenv('motus_userAuth'), '')) {
    message("Logging in to Motus...")
    if (file.exists('pass')) {
      Sys.setenv(motus_userLogin = 'ajmacdonald3') # PUT YOUR LOGIN NAME HERE (in quotes)
      Sys.setenv(motus_userPassword = readLines('pass'))
    }
    tryCatch({
      motus:::local_auth()
      Sys.setenv(motus_userAuth = TRUE)
      message("Login successful!")
    }, error = function(e) {
      stop("Login failed! '", e$message,"'", call. = F)
    })
  }
  ## Load scripts
  require(motus)
  require(tidyverse)

  # Set session time to GMT
  Sys.setenv(tz = "GMT")

  # Select a folder where data is stored
  #dir <- '../Data/'
  message(paste0(dir, '/receiver-deployments.csv'))

  # Load all receiver deployments
  recv.df <- read.csv(paste0(dir, '/receiver-deployments.csv')) %>%
    mutate(tsStart = as.POSIXct(tsStart, origin = '1970-01-01'),
           tsEnd = as.POSIXct(tsEnd, origin = '1970-01-01'),
           dtStart = as.Date(dtStart),
           dtEnd = as.Date(dtEnd))

  select.recv.df <- recv.df %>%
    filter(receiverID == serno | identical(serno, ''),
           grepl(site.name, deploymentName, T),
           is.na(projID) | recvProjectID == projID | identical(projID, ''))

  if (!identical(serno,'')) {
    recvs <-  c(serno)
  } else {
    recvs <- unique(select.recv.df$receiverID)
  }

  message(paste0(c(site.name, serno, projID), collapse = '-'))
  for (recv in recvs) {

    message(" --- START --- ")

    dbname <- paste0(dir, '/', recv, '.motus')

    message(paste0("File ", dir, '/', recv, ".rds", ifelse(file.exists(paste0(dir, '/', recv, ".rds")), " exists!", " does not exist!")))
   # if(F) {
      if (file.exists(paste0(dir, '/', recv, ".rds")) & file.exists(dbname)) {
        message("Checking if new data needs to be downloaded...")
        dbstatus <- tellme(projRecv = recv, dir = dir)
        downloadData <- rowSums(dbstatus)[[1]] != 0
      } else {
        downloadData <- T
      }
    #}

    #downloadData <- F

    if (downloadData) {

      message("Downloading data...")

      tryCatch({
        motus.sql <- tagme(projRecv = recv, new = !file.exists(dbname), update = T, forceMeta = T, dir = dir)
      }, error = function(x){
        message(paste0('There was an error: ', x))
      })


    } else {

      message("Checking for metaData...")

      tryCatch({
        motus.sql <- tagme(projRecv = recv, new = !file.exists(dbname), update = F, forceMeta = T, dir = dir)
      }, error = function(x){
        message(paste0('There was an error: ', x))
      })
    }

    if (exists("motus.sql")) {
      message("Collapsing 'alltags' virtual table...")
      df <- motus.sql %>% tbl('alltagsGPS') %>% collect() %>% as_tibble()

      message("Getting the Motus runs filters...")
      runs.df <- motus.sql %>% tbl('runs') %>% collect() %>% as_tibble() %>%
        select(runID, motusFilter)

      message("Joining runs filters with alltags data frame...")

      # Save an RDS file (compact, faster loading)
      message("Saving RDS...")
      df %>% left_join(runs.df, by='runID') %>% saveRDS(paste0(dir, '/', recv, ".rds"))
    }

    message("Database is up to date!")

    message(" --- DONE --- ")
  }
}
