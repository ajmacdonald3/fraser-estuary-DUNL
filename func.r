makeReport <- function(receiver.name = '',
                       receiver.deploy.name = '',
                       projectID = NA, # If NA, will look for matching receiver.deploy.names in all deployments
                       printable = F, # Removes hyperlinks, etc.
                       data.dir = 'E:/Data/',
                       save.dir = 'E:/Data/Reports/'
                       ) {
  require(lubridate)
  require(tidyverse)
  receiver.name <- receiver.name
  receiver.deploy.name <- ifelse(is.na(receiver.deploy.name), receiver.name, receiver.deploy.name)
  language <- 'EN' # English = EN, Francais = FR, Espaniol = ES


  pdf.fileName <- gsub(' ', '-', paste0(receiver.name, ' Motus Station Report.pdf'), fixed = T)

  pdf.exists <- file.exists(paste0(save.dir, pdf.fileName))

  if (!pdf.exists |
      (pdf.exists
    &
    Sys.time() - file.info(paste0(save.dir, pdf.fileName))$mtime > days(1))) {

    message(paste0('Starting report for ', receiver.name))


  # for pdf reports
     rmarkdown::render(input = "stationReport_v1.0.rmd",
             output_format = "pdf_document",
             output_file = pdf.fileName,
             output_dir = save.dir)
  } else {
    message(paste0('A report for ', receiver.name, ' has already been made in the past 24 hours'))
  }
}
