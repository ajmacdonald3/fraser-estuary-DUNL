# load packages
library(motus)
library(motusData)
library(DBI)
library(RSQLite)
library(tidyverse)
library(dbplyr)
library(lubridate)
library(mapdata)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(readxl)

# set timezone
Sys.setenv(TZ = "GMT")

# source functions
source("./analysis/useful_functions.R")

# set up basic things for maps (need rgdal package installed) 
theme_set(theme_bw())
world <- ne_countries(scale = "medium", returnclass = "sf") # need rgeos package installed
#lakes <- ne_download(scale = "medium", type = 'lakes', category = 'physical',
#                     returnclass = "sf", destdir = "./map-data/lakes") # only need this first time downloading
lakes <- ne_load(type = "lakes", scale = "medium", category = 'physical',
                 returnclass = "sf",
                 destdir = paste0(getwd(), "./map-data/lakes")) # use this if already downloaded shapefiles

# download/update project database and receiver databases (change to new = FALSE after first time downloading)
tagme(projRecv = 349, new = FALSE, update = TRUE, forceMeta = TRUE, dir = "sql-databases/")

# file path to sql database files
path <- "./sql-databases/"

#### pull out metadata for all project 47 tags #### 

# in the motus species list REKN = 4670
# double check species codes
sqldb <- tagme(projRecv = 349, new = FALSE, update = FALSE, dir = path)

sp.list <- c(4820) # dunlin

species <- tbl(sqldb, "species") %>%
  filter(id %in% sp.list) %>%
  collect() %>%
  as.data.frame()

# get tag metadata
tagmet <- tbl(sqldb, "tagDeps")
  
tags.meta <- tagmet %>%
#  select(deployID, tagID, projectID, tsStart, tsEnd, deferSec, speciesID,
#         markerNumber, latitude, longitude, fullID, comments) %>%
  filter(speciesID %in% sp.list) %>%
  collect() %>%
  as.data.frame() %>% # once in df format, can format dates with lubridate
  mutate(tsStart = as_datetime(tsStart, tz = "UTC", origin = "1970-01-01"),
         tsEnd = as_datetime(tsEnd, tz = "UTC", origin = "1970-01-01")) %>% 
  distinct() %>%
  filter(!is.na(longitude)) %>%
  group_by(tagID) %>%
  mutate(n = n()) %>%
  ungroup()

# get ambiguous tag info
allambigs <- tbl(sqldb, "allambigs")
  
ambig <- allambigs %>%
           collect() %>%
           as.data.frame()

# check if any tags interested in have ambiguous detections
tags.ambigs <- ambig %>% 
  group_by(ambigID) %>%
  mutate(inds = ifelse(any(motusTagID %in% tags.meta$tagID), "yes","no")) %>%
  filter(inds == "yes") %>%
  ungroup() %>%
  rename(tagID = motusTagID)

# get tag list (includes tags interested in and relevant ambiguous tags)
tag.list <- tags.meta %>%
  select(tagID) %>%
  bind_rows(., tags.ambigs %>% select(tagID)) %>%
  distinct() %>% pull()

# alltags dataset (~ 4 s to run 1-Mar-2021)
system.time({

alltag <- tbl(sqldb, "alltags")

tags.detect <- alltag %>%
             filter(tagProjID == 349 & motusTagID %in% tag.list) %>%
#             select(hitID, runID, batchID, ts, sig, port, noise, freqsd, motusTagID,
#                    ambigID, runLen, tagProjID, tagDeployID, tagDeployStart, tagDeployEnd,
#                    tagDepLat, tagDepLon, deviceID, recvDeployID, recv,
#                    speciesSci, markerNumber, mfgID) %>%
             collect() %>%
             as.data.frame() %>%
             mutate(ts = as_datetime(ts, tz = "UTC", origin = "1970-01-01"),
                    tagDeployStart = as_datetime(tagDeployStart, tz = "UTC", origin = "1970-01-01"),
                    tagDeployEnd = as_datetime(tagDeployEnd, tz = "UTC", origin = "1970-01-01"))
  
})

# get receiver metadata
# get tag metadata
recvmet <- tbl(sqldb, "recvDeps")
  
recv.meta <- recvmet %>%
           collect() %>%
           as.data.frame() %>% # once in df format, can format dates with lubridate
           mutate(tsStart = as_datetime(tsStart, tz = "UTC", origin = "1970-01-01"),
                  tsEnd = as_datetime(tsEnd, tz = "UTC", origin = "1970-01-01"))
  


recv.meta <- recv.meta %>%
  rename(recvDeployID = deployID, recvProjectID = projectID) %>% 
  distinct() %>%
  filter(!is.na(longitude)) %>%
  group_by(deviceID) %>%
  mutate(n = n()) %>%
  ungroup()

recvs <- tags.detect %>%
  select(receiverID = recv, recvDeployID) %>%
  filter(!is.na(recvDeployID)) %>%
  distinct() %>% 
  left_join(., recv.meta) %>%
  select(receiverID, recvDeployID, recvProjectID, siteName,
         latitude, longitude, receiverType, tsStart, tsEnd, isMobile) %>%
  filter(!is.na(longitude))

# join recv metadata with detection data
tags.detect <-
  left_join(tags.detect, recvs %>%
              select(recv = receiverID, recvDeployID, receiverType, recvName = siteName,
                     recvLat = latitude, recvLon = longitude, isMobile, recvProjID = recvProjectID) %>%
              distinct())

# remove detections that occur outside of deployment time window or if no recv location
# add some useful variables
tags.detect <- tags.detect %>% 
  mutate(date = as.Date(ts),
         year = as.numeric(format(ts,'%Y')),
         month = as.numeric(format(ts,'%m'))) %>%
  filter(ts >= tagDeployStart & ts <= tagDeployEnd,
         !is.na(recvLat))

# determine if any ambiguous detections associated with tags we are interested in
tags.ambigs <- ambig %>% 
  group_by(ambigID) %>%
  mutate(inds = ifelse(any(motusTagID %in% tag.list), "yes","no")) %>%
  filter(inds == "yes") %>%
  distinct() %>%
  ungroup() %>%
  rename(tagID = motusTagID)

# 4 (of 82) ambiguous tags (5%) - remove from dataset, can examine more closely later if desired
# tag list (tags interested in with ambiguous tags removed)
tag.list <- tag.list %>%
  as.data.frame() %>%
  filter(!. %in% tags.ambigs$tagID) %>%
  pull()

tag.list <- tags.meta %>% 
  select(tagID) %>% 
  distinct() %>% pull()

# filter detections to just birds interested in
tags.detect <- tags.detect %>%
  filter(motusTagID %in% tag.list)

# get motus filter results and add to dataframe
# get runs info
runs <- tbl(sqldb, "runs")
  
runs <- runs %>%
           select(runID, motusFilter) %>%
           collect() %>%
           as.data.frame()

# keep only runs associated with detection data
runs <- runs %>% 
  filter(runID %in% tags.detect$runID) %>%
  distinct()

# join runs with detection data
tags.detect <- left_join(tags.detect, runs)

# update filter 
# detections with freqsd > 0.1 (SG recvs only) --> these are generally false
tags.detect <- tags.detect %>% 
  mutate(motusFilter = ifelse(freqsd > 0.1, 0, motusFilter))

# filter out birds that have only detections that are likely false
tags.detect <- tags.detect %>%
  group_by(motusTagID) %>%
  mutate(af = sum(motusFilter)) %>%
  filter(af > 0) %>%
  ungroup() %>%
  select(-af)

# update list of tags after basic filter
tag.list <- tags.detect %>% 
  select(motusTagID) %>% 
  distinct() %>% pull()

# filter tag metadata to only tags remaining after basic filter
tags.meta <- tags.meta %>% 
  filter(tagID %in% tag.list)

# filter detections to only tags remaining after basic filter
tags.detect <- tags.detect %>% 
  filter(tagDeployID %in% tags.meta$deployID)

# save data
saveRDS(tags.detect, file = "processed-data/DUNL-tag-detections.rds")
write.csv(tags.detect, file = "processed-data/DUNL-tag-detections.csv", row.names = FALSE)

saveRDS(tags.meta, file = "processed-data/DUNL-tag-dep-metadata.rds")
write.csv(tags.meta, file = "processed-data/DUNL-tag-dep-metadata.csv", row.names = FALSE)

# load data ####
tags.detpro <- readRDS("processed-data/DUNL-tag-detections.rds")
names(tags.detpro)

# 1. check for tags deployed more than once
dtag <- tags.detpro %>% select(motusTagID, tagDeployID) %>%
  distinct() %>%
  group_by(motusTagID) %>%
  mutate(n = n()) %>%
  ungroup()
any(dtag$n > 1)
# filter out old SESA deployment
# tags.detpro <- tags.detpro %>% 
#   filter(!tagDeployID == 23721)

# 2. examine all of the tags for additional false positives
str(tags.detpro)
summary(tags.detpro)

# calculate percent of detection classified as 'good' by the filter for each tag-recv-day combination
tags.detpro <- tags.detpro %>%
  group_by(date, recvName, motusTagID) %>%
  mutate(n = n(),
         ntrue = sum(motusFilter),
         ptrue = ntrue/n) %>%
  ungroup()

# list of unique motusTagIDs
dunl.tagids <- tags.detpro %>% select(motusTagID) %>% arrange(motusTagID) %>% distinct() %>% pull()

# 3. filter out false positives
#### round 1 of filtering ####
# remove detections where less than half considered 'good' by filter for each tag-recv-day combination
dunl.tags.filt <- tags.detpro %>%
  mutate(ftemp = ifelse(ptrue > 0.50, 1, 0),
         motusFilter = ftemp) %>%
  filter(ftemp == 1)

# list of tags - 44 tags
dunl.tags.filt.list <- dunl.tags.filt %>%
  select(motusTagID) %>%
  arrange(motusTagID) %>%
  distinct() %>%
  pull()

# calculate site transitions
# check for birds with weird movements (e.g. too fast, 'wrong direction')
transit.check <- dunl.tags.filt %>%
  do(add_deploy_loc(.)) %>%
  do(site_transit_min(.))

transit.check.suspect <- transit.check %>% 
  filter(suspect.transit == "suspect")
# none suspect

# classify transitions as connected, or not_connected
tc <- transit.check %>%
  mutate(state = ifelse(dist.min == 0 | rate >= 5, "connected", "not_connected")) %>%
  filter(state == "connected") %>%
  select(motusTagID, ts.x, lat.x, lon.x, lat.y, lon.y)

# map filtered detections and save maps for each bird/tag as a png
# set limits to map based on locations of detections
# simplify detection data for plotting
dunl.tags.path.filt <- fun.getpath(dunl.tags.filt)

for (i in 1:length(dunl.tags.filt.list)){
  
  bird <- dunl.tags.path.filt %>% filter(motusTagID == dunl.tags.filt.list[i])
  
  pbird <- tc %>% filter(motusTagID == dunl.tags.filt.list[i]) %>%
    distinct()
  
  xmin <- min(bird$recvLon, bird$tagDepLon, na.rm = TRUE) - 2
  xmax <- max(bird$recvLon, bird$tagDepLon, na.rm = TRUE) + 2
  ymin <- min(bird$recvLat, bird$tagDepLat, na.rm = TRUE) - 1
  ymax <- max(bird$recvLat, bird$tagDepLat, na.rm = TRUE) + 1
  
  png(filename = paste0("figures/diagnostic-plots/round-1/",
                        bird$motusTagID,"-","map.png"),
      width=8, height=8, units="in", res=600)
  
  print(ggplot(data = world) +
          geom_sf() +
          geom_sf(data = lakes, fill = "white") +
          coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
          geom_point(data = bird, aes(recvLon, recvLat, colour = month), alpha = 0.5, size = 5) +
          geom_segment(data = pbird, aes(x = lon.x, y = lat.x, xend = lon.y, yend = lat.y),
                       linetype = "dashed",
                       arrow = arrow(angle = 18, type = "closed", length = unit(0.1, "inches"))) +
          geom_point(data = bird,
                     aes(tagDepLon, tagDepLat), colour = "red",
                     shape = 8, size = 4) +
          scale_color_gradientn(colors = rev(viridis_pal()(12)), limits = c(1,12)) +
          ggtitle(paste0("motusTagID: ", bird$motusTagID[1])))
  dev.off()
  
}
rm(i)

# make sigs strength plots
for (i in 1:length(dunl.tags.filt.list)){
  
  png(filename = paste0("figures/diagnostic-plots/round-1/", "sigstrength-",
                        dunl.tags.filt.list[i],".png"),
      width=11, height=7.5, units="in", res=600)
  print(plotTagSig_mod(dunl.tags.filt, motusTagID = dunl.tags.filt.list[i]))
  dev.off()
}
rm(i)

# write csv to keep track of issues spotted via visual inspection
write.csv(dunl.tags.filt %>% select(motusTagID) %>% distinct() %>% arrange(motusTagID),
          file = "./processed-data/filtering_notes_tags_round1.csv", row.names=FALSE)

# remove false positives
dunl.tags.filt2 <- dunl.tags.filt %>% 
  filter(!recvName == "Grève de Tadoussac") %>% 
  filter(!recvName == "Clear Lake") %>% 
  filter(!recvName == "Breakwater") %>% 
  filter(!recvName == "Old Cut")

# write rds file of detection data after the basic filter has been run
saveRDS(dunl.tags.filt2, file = "./processed-data/dunl-tag-detections-clean.rds")

# plot dunl map
species_clean <- dunl.tags.filt2 %>% 
  select(speciesSci) %>% 
  distinct() %>% pull()

region <- ne_states(country = c("Canada", "United States of America"), returnclass = "sf") %>% 
  filter(name %in% c("British Columbia", "Yukon", "Northwest Territories", "Washington", "Alaska"))

for (i in 1:length(species_clean)){
  
  species <- dunl.tags.filt2 %>% filter(speciesSci == species_clean[i]) #%>% 
    #filter(!recvName == "Cordova")
  
  xmin <- min(species$recvLon, species$tagDepLon, na.rm = TRUE) - 1
  xmax <- max(species$recvLon, species$tagDepLon, na.rm = TRUE) + 1
  ymin <- min(species$recvLat, species$tagDepLat, na.rm = TRUE) - 0.5
  ymax <- max(species$recvLat, species$tagDepLat, na.rm = TRUE) + 0.5
  
  png(filename = paste0("figures/",
                        species$speciesSci,"-","map2.png"),
      width=8, height=8, units="in", res=600)
  
  print(ggplot() +
          geom_sf(data = region, colour = NA) +
          geom_sf(data = lakes, fill = "white", colour = NA) +
          coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
          geom_path(data = species,
                    aes(x = recvLon, y = recvLat, group = motusTagID),
                    colour = "black", size = 0.5, alpha = 0.3) +
          geom_point(data = species,
                     aes(x = recvLon, y = recvLat, group = motusTagID, colour = month),
                     size = 3, alpha = 0.3) +
          geom_point(data = species,
                     aes(tagDepLon, tagDepLat), colour = "red",
                     shape = 17, size = 1) +
          scale_color_gradientn(colors = rev(viridis_pal()(12)), limits = c(1,12)) +
          ggtitle(species$speciesSci[1]))
  
  dev.off()
}
rm(i)
