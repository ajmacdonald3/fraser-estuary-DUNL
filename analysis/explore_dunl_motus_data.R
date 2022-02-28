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
library(bcmaps)

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

theme_set(theme_bw())

# load data
dunl.data <- readRDS("./processed-data/dunl-tag-detections-clean.rds")

# get path and divide by banding location
dunl.path <- fun.getpath(dunl.data)

dunl.path <- dunl.path %>% 
  mutate(tagDepLoc = ifelse(tagDepLon > (-123.02), "Boundary Bay", "Roberts Bank")) # lat 49.07

dunl.path.bb <- dunl.path %>% 
  filter(tagDepLoc == "Boundary Bay")

dunl.path.rb <- dunl.path %>% 
  filter(tagDepLoc == "Roberts Bank")

# set up map
region <- ne_states(country = c("Canada", "United States of America"), returnclass = "sf") %>% 
  filter(name %in% c("British Columbia", "Washington"))

region <- ne_countries(scale = "large", country = c("Canada", "United States of America"),
                       returnclass = "sf")

xmin <- min(dunl.path$recvLon, dunl.path$tagDepLon, na.rm = TRUE) - 1
xmax <- max(dunl.path$recvLon, dunl.path$tagDepLon, na.rm = TRUE) + 1
ymin <- min(dunl.path$recvLat, dunl.path$tagDepLat, na.rm = TRUE) - 0.5
ymax <- max(dunl.path$recvLat, dunl.path$tagDepLat, na.rm = TRUE) + 0.5

# paths of all connections, direct flights or not  
png(filename = paste0("figures/dunl-detection-paths-map3.png"),
    width=8, height=8, units="in", res=600)

ggplot(data = region) +
  geom_sf(colour = NA) +
  geom_sf(data = lakes, fill = "white", colour = NA) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  geom_path(data = dunl.path.bb, aes(recvLon, recvLat, group = motusTagID),
            colour = "#29AF7FFF", size = 0.5) +
  geom_path(data = dunl.path.rb, aes(recvLon, recvLat, group = motusTagID),
            colour = "#482677FF", size = 0.5)

dev.off()

# summarize transitions
dunl.data <- dunl.data %>% 
  mutate(tagDepLoc = ifelse(tagDepLon > (-123.02), "Boundary Bay", "Roberts Bank"))

dunl.data.bb <- dunl.data %>% 
  filter(tagDepLoc == "Boundary Bay")

dunl.data.rb <- dunl.data %>% 
  filter(tagDepLoc == "Roberts Bank")

transit.bb <- dunl.data.bb %>%
  do(add_deploy_loc(.)) %>%
  do(site_transit_min(.))

transit.rb <- dunl.data.rb %>%
  do(add_deploy_loc(.)) %>%
  do(site_transit_min(.))

trans.sum.bb <- transit.bb %>% 
  group_by(recvName.x, recvName.y, lat.x, lon.x, lat.y, lon.y) %>% 
  summarize(n.trans = n()) %>% 
  ungroup() %>% 
  mutate(tagDepLoc = "Boundary Bay")

trans.sum.rb <- transit.rb %>% 
  group_by(recvName.x, recvName.y, lat.x, lon.x, lat.y, lon.y) %>% 
  summarize(n.trans = n()) %>% 
  ungroup() %>% 
  mutate(tagDepLoc = "Roberts Bank")

trans.sum <- bind_rows(trans.sum.bb, trans.sum.rb)

xmin <- min(trans.sum$lon.x, trans.sum$lon.y, na.rm = TRUE) - 0.25
xmax <- max(trans.sum$lon.x, trans.sum$lon.y, na.rm = TRUE) + 0.25
ymin <- min(trans.sum$lat.x, trans.sum$lat.y, na.rm = TRUE) - 0.25
ymax <- max(trans.sum$lat.x, trans.sum$lat.y, na.rm = TRUE) + 0.25

png(filename = paste0("figures/dunl-transitions3.png"),
    width=8, height=6, units="in", res=600)

ggplot(data = region) +
  geom_sf() +
  geom_sf(data = lakes, fill = "white", colour = NA) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  geom_segment(data = trans.sum, aes(x = lon.x, y = lat.x, xend = lon.y, yend = lat.y,
                                     colour = tagDepLoc, size = n.trans), alpha = 0.7) +
  facet_wrap(. ~ tagDepLoc) +
  scale_size(range = c(0.2, 3)) +
  theme(legend.position="none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"))

dev.off()

# summarize by detection per receiver per deploy loc
detect.sum.day <- dunl.data %>% 
  mutate(doy = yday(ts)) %>% 
  select(motusTagID, tagDepLoc, recvName, recvLat, recvLon, doy) %>%
  distinct() %>% 
  group_by(tagDepLoc, recvName, recvLat, recvLon) %>% 
  summarize(n.det = n()) %>% 
  filter(!recvName == "Cordova") %>% 
  ungroup()

png(filename = paste0("figures/dunl-transitions3.png"),
    width=8, height=6, units="in", res=600)

ggplot(data = region) +
  geom_sf() +
  geom_sf(data = lakes, fill = "white", colour = NA) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  geom_point(data = detect.sum, aes(x = recvLon, y = recvLat,
                                    colour = tagDepLoc, size = n.det), alpha = 0.7) +
  facet_wrap(. ~ tagDepLoc) +
  #scale_size(range = c(0.2, 3)) +
  theme(legend.position="none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"))

dev.off()

#bc maps
bc <- bc_bound_hres()
bc_neigh <- bc_neighbours()
bc_water <- watercourses_15M()

# convert detections to sf object
detect.sum.sf <- st_as_sf(detect.sum.day %>%
                            filter(!recvName %in% c("Breckenridge Bluff", "Sidney Island")),
                          coords = c("recvLon", "recvLat"), crs = 4326)
detect.sum.sf <- st_transform(detect.sum.sf, crs = 3005)

xmin <- st_bbox(detect.sum.sf)[1] - 10000
xmax <- st_bbox(detect.sum.sf)[3] + 5000
ymin <- st_bbox(detect.sum.sf)[2] - 5000
ymax <- st_bbox(detect.sum.sf)[4] + 10000


png(filename = paste0("figures/dunl-days-detected3.png"),
    width=8, height=6, units="in", res=600)

ggplot() +
  #geom_sf(data = bc_neigh %>% filter(name == "Washington"), colour = NA) +
  geom_sf(data = bc, colour = NA) +
  #geom_sf(data = lakes, fill = "white", colour = NA) +
  geom_sf(data = detect.sum.sf,
          aes(colour = tagDepLoc, size = n.det), alpha = 0.7) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  facet_wrap(tagDepLoc ~ .) +
  scale_size(range = c(1, 10), name = "Dunlin detections") +
  scale_color_manual(values = c("#29AF7FFF", "#482677FF"), guide = "none") +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold", size = 14))

dev.off()

# number of days each bird detected
# summarize by detection per receiver per deploy loc
days.detect <- dunl.data %>% 
  mutate(doy = yday(ts)) %>% 
  select(motusTagID, doy) %>%
  distinct() %>% 
  group_by(motusTagID) %>% 
  summarize(n.days = n()) %>% 
  ungroup()

days.recv.detect <- dunl.data %>% 
  mutate(doy = yday(ts)) %>% 
  select(motusTagID, recvName, doy) %>%
  distinct() %>% 
  group_by(motusTagID, recvName) %>% 
  summarize(n.days = n()) %>% 
  ungroup()

detect.recv <- dunl.data %>% 
  group_by(motusTagID) %>% 
  summarize(nRecv = length(unique(recvName))) %>% 
  ungroup()

# plot paths
# within each group repeat each point 
# then slice the first and last out and 
# add a variable called linegroup, which provides grouping for start and endpoints of each line
dunl.path <- dunl.path %>% 
  filter(!recvName %in% c("Cordova", "Breckenridge Bluff", "Sidney Island"))

dunl.path.sf <- dunl.path %>%
  group_by(motusTagID) %>%
  slice(rep(1:n(), each = 2)) %>%
  slice(-c(1, n())) %>%
  mutate(linegroup = lapply(1:(n()/2), function(x) rep(x, 2)) %>% unlist) %>% 
  ungroup()

# create linestring sf object by summarizing the points, 
# grab the last survivor and direction value of each group (i.e. the 'endpoint' value)
dunl.path.sf <- st_as_sf(dunl.path.sf, coords = c("recvLon", "recvLat"), crs = 4326) %>%
  group_by(motusTagID, linegroup) %>%
  summarise(tagDepLoc = last(tagDepLoc), do_union = FALSE) %>%
  st_cast("LINESTRING")

dunl.path.sf <- st_as_sf(dunl.path.sf, coords = c("recvLon", "recvLat"), crs = 4326)
dunl.path.sf <- st_transform(dunl.path.sf, crs = 3005)

recv.locs <- dunl.path %>% 
  select(recvName, recvLat, recvLon) %>% 
  filter(!recvName == "Cordova") %>% 
  distinct() %>% 
  st_as_sf(coords = c("recvLon", "recvLat"), crs = 4326) %>% 
  st_transform(crs = 3005)

xmin <- st_bbox(dunl.path.sf)[1] - 10000
xmax <- st_bbox(dunl.path.sf)[3] + 5000
ymin <- st_bbox(dunl.path.sf)[2] - 5000
ymax <- st_bbox(dunl.path.sf)[4] + 10000

png(filename = paste0("figures/dunl-tracks3.png"),
    width=8, height=8, units="in", res=600)

ggplot() +
  geom_sf(data = bc, colour = NA) +
  #geom_sf(data = bc_water) +
  geom_sf(data = dunl.path.sf, aes(group = motusTagID,
                                     colour = tagDepLoc), size = 0.5, alpha = 0.3) +
  geom_sf(data = recv.locs, colour = "black", size = 2) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  facet_wrap(tagDepLoc ~ .) +
  scale_color_manual(values = c("#29AF7FFF", "#482677FF"), guide = "none") +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold", size = 14))

dev.off()

# individuals
detect.sf <- dunl.data %>% 
  mutate(tagDepLoc = ifelse(tagDepLat > 49.07, "Boundary Bay", "Roberts Bank")) %>% 
  dplyr::select(motusTagID, recvName, recvLat, recvLon, tagDepLoc) %>% 
  st_as_sf(coords = c("recvLon", "recvLat"), crs = 4326) %>% 
  st_transform(crs = 3005)

png(filename = paste0("figures/ind-dunl-tracks.png"),
    width=8, height=6, units="in", res=600)

ggplot() +
  geom_sf(data = bc, colour = NA) +
  #geom_sf(data = bc_water) +
  geom_sf(data = dunl.path.sf, aes(group = motusTagID,
                                   colour = tagDepLoc), size = 0.5) +
  geom_sf(data = detect.sf, aes(colour = tagDepLoc), size = 1, alpha = 0.1) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  facet_wrap(motusTagID ~ ., ncol = 8) +
  scale_color_manual(values = c("#29AF7FFF", "#482677FF")) +
  theme(legend.position = c(0.75, 0.05),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold", size = 6)) +
  guides(color = guide_legend(title = "Tagging\nlocation", title.position = "left"))

dev.off()

# colour by date
detect.sf <- dunl.data %>% 
  mutate(tagDepLoc = ifelse(tagDepLat > 49.07, "Boundary Bay", "Roberts Bank")) %>% 
  dplyr::select(motusTagID, ts, recvName, recvLat, recvLon, tagDepLoc) %>%
  mutate(doy = yday(ts)) %>% 
  select(-ts) %>% 
  distinct() %>% 
  st_as_sf(coords = c("recvLon", "recvLat"), crs = 4326) %>% 
  st_transform(crs = 3005)

png(filename = paste0("figures/date-dunl-tracks.png"),
    width=8, height=6, units="in", res=600)

ggplot() +
  geom_sf(data = bc, colour = NA) +
  #geom_sf(data = bc_water) +
  #geom_sf(data = dunl.path.sf, aes(group = motusTagID,
  #                                 colour = tagDepLoc), size = 0.5) +
  geom_sf(data = detect.sf, aes(color = doy), size = 1) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) +
  facet_wrap(motusTagID ~ ., ncol = 8) +
  scale_color_viridis_c(option = "magma") +
  theme(legend.position = c(0.75, 0.05),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold", size = 6)) +
  guides(color = guide_legend(title = "Date", title.position = "left",
                              direction = "horizontal"))

dev.off()

