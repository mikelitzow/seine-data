## estimate temperature for each site / bay / julian day / year combination

library(tidyverse)
library(mgcv)
library(chron)

theme_set(theme_bw())
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# load data

# first, long-term Kodiak sites

# this is the version Ben Laurel shared with me on 9/11/20
dat <- read.csv("data/Kodiak.seine.data.2.csv")

head(dat)

names(dat)[1] <- "year"

dat$Date <- dates(as.character(dat$Date))

dat$julian <- lubridate::yday(dat$Date)

ben.dat <- dat %>%
  select(year, julian, Site.Name, Region, Temperature..measured.)

names(ben.dat)[3:4] <- c("site", "bay")
names(ben.dat)[5] <- "temperature"

change <- ben.dat$site == "laminaria #2"
ben.dat$site[change] <- "Laminaria #2"

## and the WGOA sites
wgoa <- read.csv("data/2018 2020 site.csv")

head(wgoa)

wgoa <- wgoa %>%
  select(Date, Site, Bay, Temp.C)

# calculate Julian day
wgoa$Date <- dates(as.character(wgoa$Date))
wgoa$julian <- lubridate::yday(wgoa$Date)
wgoa$year <- years(wgoa$Date)

head(wgoa)


names(wgoa)[2:4] <- c("site", "bay", "temperature")

wgoa <- wgoa %>%
  select(year, julian, site, bay, temperature)

# check for overlap - are 2020 Cook's sites entered twice?
# remove 2020 Cook's from wgoa as these are already entered
View(filter(wgoa, year==2020, bay=="Cooks"))
View(filter(ben.dat, year==2020, bay=="Cook Bay"))

# yes, there are repeats there!
drop <- wgoa$year==2020 & wgoa$bay=="Cooks"
wgoa <- wgoa[!drop,]

temp.dat <- rbind(ben.dat, wgoa)
temp.dat$year <- as.numeric(temp.dat$year)

# drop rows with no temp
temp.dat <- na.omit(temp.dat)

write.csv(temp.dat, "./data/seine_temperatures.csv")
