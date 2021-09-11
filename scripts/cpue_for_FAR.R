## produce cod & pollock cpue and ssb for pre-settlement FAR analysis

library(tidyverse)
library(chron)

dat <- read.csv("./output/wgoa.cod.poll.cpue.csv")

# put into form previously used in brms modeling

dat <- dat %>%
  select(-age.1, -station, -temp.c) %>% # restricting to age 0!
  pivot_wider(names_from = species, values_from = age.0) 


# rename to match earlier conventions
names(dat)[5:6] <- c("cod", "pollock")

dat <- dat %>%
  select(year, julian, site, bay, cod, pollock)

# load long-term data from Kodiak sites
# this is the version Ben Laurel shared with me on 9/11/20
long.dat <- read.csv("./data/Kodiak.seine.data.2.csv")

head(long.dat)

names(long.dat)[1] <- "year"

long.dat$Date <- dates(as.character(long.dat$Date))

long.dat$julian <- lubridate::yday(long.dat$Date)

long.dat$cod <- long.dat$X..Pacific.cod
long.dat$pollock <- long.dat$X..Pollock

head(long.dat)

long.dat <- long.dat %>%
  select(year, julian, Site.Name, Region, cod, pollock)

head(long.dat)
names(long.dat)[3:4] <- c("site", "bay")

unique(long.dat$bay)
unique(long.dat$site)


# now make some exploratory plots!

ggplot(filter(long.dat, bay=="Cook Bay"), aes(log(cod+1))) + 
  geom_histogram(bins=10) + 
  facet_grid(year~site)

# this looks fine - sites are represented throughout, no early sites that aren't represented in other years

ggplot(filter(long.dat, bay=="Anton Larson Bay"), aes(log(cod+1))) + 
  geom_histogram(bins=10) + 
  facet_grid(year~site)

# looks like laminaria #2 should be Laminaria #2!

change <- long.dat$site == "laminaria #2"
long.dat$site[change] <- "Laminaria #2"

# double-check there are no xtra sites in 2006
View(filter(long.dat, year==2006)) #looks good 

head(dat); head(long.dat)

# check for 2020 Cooks Bay sites entered twice
View(filter(dat, bay=="Cooks"))
View(filter(long.dat, year==2020, bay=="Cook Bay"))

# yes, there are repeats - remove Cooks from dat
dat <- dat %>%
  filter(bay != "Cooks")

# combine 
dat <- rbind(dat, long.dat)


# add ssb
poll <- read.csv("./data/pollock.ssb.csv")
dat$pollock.ssb <- poll$ssb[match(dat$year, poll$year)]

cod <- read.csv("./data/cod.ssb.csv")
dat$cod.ssb <- cod$ssb[match(dat$year, cod$year)]

# some plots
ggplot(dat, aes(log(cod+1))) + 
  geom_histogram(bins=10) + 
  facet_grid(bay~year, scale="free_y")

# looks good...

write.csv(dat, "./output/cpue.for.FAR.csv", row.names = F)
