##################################################################
####
#### MOTUS STATION REPORT
#### Lucas Berrigan
#### 2 November 2020
####
##################################################################
####
#### ABOUT
####
#### This script uses a function to produce a station summary for one location.
####  - It will find all receiver serial numbers for a single 'station' based on a common name and project ID.
####  - It automatically download data from all receivers during that stations history.
####  - An option for printable reports removes track hyperlinks (keep false if not printing)
####  - To change title and authorship, edit the first two lines of the 'stationReport_v#.#.rmd' document.
####  - Report bugs via GitHub
####
##################################################################
####

# Load the function
source('func.r')

# Run the function
makeReport(
            receiver.name = 'Boundary Bay Field', # The name you want to have displayed on the report
            receiver.deploy.name = 'Boundary Bay Field', # A character string that is common among all deployment names
            projectID = 273, # If NA, will look for matching receiver.deploy.names in all deployments
            printable = F, # Removes hyperlinks, etc.
            data.dir = 'station-reports/data/', # Where the Motus data will be stored
            save.dir = 'station-reports/reports/' # Where the reports will be stored
           )

# If you have multiple reports to make, save the variables in a data frame and
