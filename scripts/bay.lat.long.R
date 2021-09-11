# get mean lat/long for each bay

library(tidyverse)
library(mgcv)
library(chron)

# this is the version Ben Laurel shared with me on 9/11/20
dat <- read.csv("data/Kodiak.seine.data.2.csv")

head(dat)

dat$lat.dec <- dat$Lat.deg + dat$Lat.min/60
dat$lon.dec <- dat$Long.deg + dat$Long.min/60

ben.dat <- dat %>%
  group_by(Region) %>%
  summarise(lat=as.numeric(mean(lat.dec, na.rm=T)), lon=mean(lon.dec, na.rm=T))

# now the WGOA sites
dat <- read.csv("data/2018 2020 site.csv")

head(dat)

our.dat <- dat %>%
  group_by(Bay) %>%
  summarise(lat=as.numeric(mean(lat, na.rm=T)), lon=mean(long, na.rm=T))

names(ben.dat)[1] <- "Bay"

bay.dat <- rbind(ben.dat, our.dat)

write.csv(bay.dat, "data/bay_lat_long.csv")
