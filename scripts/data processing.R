library(tidyverse)
library(mgcv)
library(chron)

theme_set(theme_bw())
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# test for a threshold temperature effect on cod CPUE

# load data

# first, long-term Kodiak sites

# dat <- read.csv(data/Kodiak.seine.data.csv) # this is the first version I was using, is missing some data!

# this is the version Ben Laurel shared with me on 9/11/20
dat <- read.csv("data/Kodiak.seine.data.2.csv")

head(dat)

names(dat)[1] <- year

dat$Date <- dates(as.character(dat$Date))

dat$julian <- lubridate::yday(dat$Date)

dat$cod <- dat$X..Pacific.cod
dat$pollock <- dat$X..Pollock

head(dat)

# save an indexing version for length analysis, below...
index.dat <- dat %>%
  select(year, julian, Site.Name, Region, cod, pollock, Temperature..measured.)

names(index.dat)[3:4] <- c(site, bay)
names(index.dat)[7] <- temperature

change <- index.dat$site == laminaria #2
index.dat$site[change] <- Laminaria #2

write.csv(index.dat, "data/index long-term cpue sites.csv")
#############################################################

dat <- dat %>%
  select(year, julian, Site.Name, Region, cod, pollock)

head(dat)
names(dat)[3:4] <- c(site, bay)

unique(dat$bay)
unique(dat$site)
# now make some exploratory plots!

ggplot(filter(dat, bay==Cook Bay), aes(log(cod+1))) + 
  geom_histogram(bins=10) + 
  facet_grid(year~site)
# this looks fine - sites are represented throughout, no early sites that aren't represented in other years

ggplot(filter(dat, bay==Anton Larson Bay), aes(log(cod+1))) + 
  geom_histogram(bins=10) + 
  facet_grid(year~site)

# looks like laminaria #2 should be Laminaria #2!

change <- dat$site == laminaria #2
dat$site[change] <- Laminaria #2

# double-check there are no xtra sites in 2006

View(filter(dat, year==2006))

# now load WGOA data
d2 <- read.csv(data/2018 2019 catch.csv)

# combine with complete 2020 catch data!
d3 <- read.csv(data/all_cpue_2020.csv)

# add a measured column
d3$Measured <- NA

unique(d3$species) # looks good, no different names for cod/pollock

#  combine with d2
d3 <- d3 %>%
  select(Station, species, Measured, CPUE)

names(d3) <- c(Station, Species, Measured, Total.catch)

# combine
d2 <- rbind(d2, d3)

# check for repeat spp names
unique(d2$Species)[order(unique(d2$Species))]

# pollock and cod are clean - no different names!
cod <- d2 %>%
  filter(Species==Pacific cod)
cod <- data.frame(Station=cod$Station, cod=cod$Total.catch)

pollock <- d2 %>%
  filter(Species==walleye pollock)
pollock <- data.frame(Station=pollock$Station, pollock=pollock$Total.catch)

wgoa.dat <- data.frame(Station=unique(d2$Station))

wgoa.dat <- left_join(wgoa.dat, cod)
wgoa.dat <- left_join(wgoa.dat, pollock)

# # log transform non-0 catches
# wgoa.dat$cod <- log(wgoa.dat$cod)
# wgoa.dat$pollock <- log(wgoa.dat$pollock)

# change NAs to 0
change <- is.na(wgoa.dat)

wgoa.dat[change] <- 0

# now need to add year, julian day, site, and bay!
d4 <- read.csv("data/2018 2020 site.csv")

head(d4)

# retain only the keeper sets

## NB! changing 115 to cpue==no!
d4$use.for.CPUE[d4$Station==115] <- "no"

d4 <- d4 %>%
  filter(use.for.CPUE=="yes")

d4 <- d4 %>%
  select(Date, Station, Site, Bay, Temp.C)

# calculate Julian day
d4$Date <- dates(as.character(d4$Date))
d4$julian <- lubridate::yday(d4$Date)
d4$year <- years(d4$Date)

head(d4)

d4 <- d4 %>%
  select(-Date)

names(d4)[1] <- "station"
d4$station <- as.factor(d4$station)

d4 <- left_join(d4, wgoa.dat)

names(d4)[2:3] <- c("site", "bay")

d4 <- d4 %>%
  select(year, julian, site, bay, cod, pollock)

# check for overlap - are 2020 Cook's sites entered twice?
# remove 2020 Cook's from wgoa as these are already entered
View(filter(d4, year==2020, bay=="Cooks"))
View(filter(dat, year==2020, bay=="Cook Bay"))

# yes, there are repeats there!
drop <- d4$year==2020 & d4$bay==Cooks
d4 <- d4[!drop,]

dat <- rbind(dat, d4)
dat$year <- as.numeric(dat$year)

# now add SST (placeholder for temp!)
sst <- read.csv(data/goa.winter.sst.csv)

names(sst)[2] <- sst

dat <- left_join(dat, sst)

# remove Kujulik
dat <- dat %>%
  filter(bay != Kujulik)

dat$year <- as.factor(dat$year)

hist(dat$julian)

# remove late samples
dat <- dat %>%
  filter(julian < 255)

unique(dat$bay)

unique(dat$site)
unique(dat$sst)

# now add cod SSB
ssb <-read.csv(data/cod.ssb.csv) # estimated female SSB fom 2019 assessment report
dat$cod.ssb <- ssb$ssb[match(dat$year, ssb$year)]

# and pollock ssb!
p.ssb <-read.csv(data/pollock.ssb.csv) # estimated female SSB fom 2019 assessment report
dat$pollock.ssb <- p.ssb$ssb[match(dat$year, p.ssb$year)]

write.csv(dat, data/cpue.data.csv)

#######################

# now length data

dat <- read.csv("data/Age-0 gadid sizes for Mike _2006-20.csv")

head(dat)

dat[is.na(dat$Region),]
unique(dat$Year) # every year except 2014 and 2019

names(dat)[2] <- "year"

unique(dat$Region) # no repeat names for Cook / Anton's!

dat$Date <- dates(as.character(dat$Date))

dat$julian <- lubridate::yday(dat$Date)

unique(dat$Species)

unique(dat$Site.Name)
names(dat)[15] <- "size"
head(dat)

# lotta odd names there!
# plot to check

ggplot(filter(dat, Species=="pcod"), aes(size)) +
  geom_histogram() +
  facet_grid(year~Site.Name)

# that's a mess!

cod <- dat %>%
  filter(Species=="pcod") %>%
  group_by(Site.Name, year) %>%
  summarise(count=n())

cod

ggplot(cod, aes(year, count)) +
  theme_bw() +
  facet_wrap(~Site.Name, scales="free_y") +
  geom_bar(stat="identity")

# ok...first thing, there's a clear problem with different site names in different years

# make them uniform
# using the same names that are used in the cpue data
change <- grep("Back", dat$Site.Name)
length(change)
dat$Site.Name[change] <- "Back Bay eelgrass"

change <- grep("Bare", dat$Site.Name)
dat$Site.Name[change]
dat$Site.Name[change] <- "Bare North"

change <- grep("Cobble", dat$Site.Name)
dat$Site.Name[change]
dat$Site.Name[change] <- "Cobble point"

change <- dat$Site.Name %in% c("Laminaria point", "Laminaria Point")
sum(change)
dat$Site.Name[change] <- "Laminaria Point"

change <- dat$Site.Name %in% c("Eelgrass point", "Eelgrass Point")
sum(change)
dat$Site.Name[change] <- "Eelgrass point"

change <- dat$Site.Name %in% c("Laminaria #2", "Llaminaria #2")
sum(change)
dat$Site.Name[change] <- "Laminaria #2"

change <- dat$Site.Name %in% c("Wooden Boat", "Wooden boat")
sum(change)
dat$Site.Name[change] <- "Wooden boat"

change <- dat$Site.Name %in% c("Middle Cove", "Middle cove")
sum(change)
dat$Site.Name[change] <- "Middle cove"

# and check
cod <- dat %>%
  filter(Species=="pcod") %>%
  group_by(Site.Name, year) %>%
  summarise(count=n())

ggplot(cod, aes(year, count)) +
  theme_bw() +
  facet_wrap(~Site.Name, scales="free_y") +
  geom_bar(stat="identity")


# limit to the sites used in CPUE analysis
cpue <- read.csv("data/cpue.data.csv", row.names = 1)

cpue <- cpue %>%
  filter(bay %in% c("Cook Bay", "Anton Larson Bay"))

unique(cpue$site)

# and limit length data to the same sites
dat <- dat %>%
  filter(Site.Name %in% unique(cpue$site))

head(dat)

dat <- dat %>%
  select(year, julian, Site.Name, Region, Species, size)

head(dat)
names(dat)[3:5] <- c("site", "bay", "species")

unique(dat$bay)
unique(dat$site)

# now make some exploratory plots!

ggplot(filter(dat, bay=="Cook Bay" & species=="pcod"), aes(size)) + 
  geom_histogram(bins=15) + 
  facet_grid(year~site)

# this looks fine - sites are represented throughout, no early sites that aren't represented in other years

ggplot(filter(dat, bay=="Anton Larson Bay" & species=="pcod"), aes(size)) + 
  geom_histogram(bins=15) + 
  facet_grid(year~site)

# now we need to bring in temperature and cpue!

# create an index column that can be used for matching!
dat$index <- apply(dat[,1:3], 1, paste, collapse = "-")

# and an indexing cpue data frame!
index.cpue <- read.csv("data/index long-term cpue sites.csv", row.names = 1)

check <- data.frame(length.site=unique(dat$site)[order(unique(dat$site))],
                    cpue.site=unique(index.cpue$site)[order(unique(index.cpue$site))])
check
# looks great!

index.cpue$index <- apply(index.cpue[,1:3], 1, paste, collapse = "-")

# keep only the data we need (temp and CPUE)
use.dat <- index.cpue %>%
  select(cod, pollock, temperature, index)

names(use.dat)[1:2] <- c("cod.cpue", "pollock.cpue")

dat <- left_join(dat, use.dat)

# now load WGOA length data
wgoa.len <- read.csv("data/gadid_len.csv")

wgoa.len[is.na(wgoa.len$Station), ]

wgoa.len$Species <- as.character(wgoa.len$Species)
  
# change species names
change <- wgoa.len$Species=="Pacific cod"
wgoa.len$Species[change] <- "pcod"

change <- wgoa.len$Species=="walleye pollock"
wgoa.len$Species[change] <- "poll"

# restrict to age-0
head(wgoa.len)

hist(wgoa.len$L)

# load wgoa catch data
d2 <- read.csv("data/2018 2019 catch.csv")

# combine with complete 2020 catch data!
d3 <- read.csv("data/all_cpue_2020.csv")

# add a measured column
d3$Measured <- NA

unique(d3$species) # looks good, no different names for cod/pollock

#  combine with d2
d3 <- d3 %>%
  select(Station, species, Measured, CPUE)

names(d3) <- c("Station", "Species", "Measured", "Total.catch")

# combine
d2 <- rbind(d2, d3)

# check for repeat spp names
unique(d2$Species)[order(unique(d2$Species))]

# pollock and cod are clean - no different names!
cod <- d2 %>%
  filter(Species=="Pacific cod")
cod <- data.frame(Station=cod$Station, cod=cod$Total.catch)

pollock <- d2 %>%
  filter(Species=="walleye pollock")
pollock <- data.frame(Station=pollock$Station, pollock=pollock$Total.catch)

wgoa.dat <- data.frame(Station=unique(d2$Station))

wgoa.dat <- left_join(wgoa.dat, cod)
wgoa.dat <- left_join(wgoa.dat, pollock)

# change NAs to 0
change <- is.na(wgoa.dat)

wgoa.dat[change] <- 0

# now need to add year, julian day, site, bay, and temperature!
d4 <- read.csv("data/2018 2020 site.csv")

head(d4)

# retain only the keeper sets

## NB! changing 115 to cpue==no!
d4$use.for.CPUE[d4$Station==115] <- "no"

d4 <- d4 %>%
  filter(use.for.CPUE=="yes")

d4 <- d4 %>%
  select("Date", "Station", "Site", "Bay", "Temp.C")

# calculate Julian day
d4$Date <- dates(as.character(d4$Date))
d4$julian <- lubridate::yday(d4$Date)
d4$year <- years(d4$Date)


head(d4)

d4[is.na(d4$Bay),]

d4 <- d4 %>%
  select(-Date)

d4 <- left_join(d4, wgoa.dat)

names(d4)[2:3] <- c("site", "bay")

# check the two wgoa df
View(wgoa.len)
View(d4)

sum(is.na(wgoa.len))
sum(is.na(d4)) # site 117 and 228 are bad so show up as NAs!


wgoa.len <- left_join(wgoa.len, d4)

wgoa.len[is.na(wgoa.len$bay),]

# clean up and combine the two TS

dat <- dat %>%
  select(-index)

names(wgoa.len)

wgoa.len <- wgoa.len %>%
  select(-Station)

names(wgoa.len) <- c("species", "length", "site", "bay", "temperature",  "julian",  "year", "cod.cpue", "pollock.cpue")

names(dat)[6] <- "length"

wgoa.len <- wgoa.len %>%
  select("year", "julian", "site", "bay", "species", "length", "cod.cpue", "pollock.cpue", "temperature")

length.dat <- rbind(dat, wgoa.len)

#####################
# NB! DROPPING NAs
length.dat <- na.omit(length.dat)

# and check for repeat bay names!
unique(length.dat$bay)

change <- length.dat$bay=="Cooks"

length.dat$bay[change] <- "Cook Bay"

# and drop Kujulik!
length.dat <- length.dat %>%
  filter(bay != "Kujulik")

write.csv(length.dat, "data/gadid length data.csv")

length.dat <- read.csv("data/gadid length data.csv", row.names = 1)

cod.length <- length.dat %>%
  filter(species=="pcod")

pollock.length <- length.dat %>%
  filter(species=="poll")

write.csv(cod.length, "data/cod length data.csv")
write.csv(pollock.length, "data/pollock length data.csv")
# # compare 
# 
# check.dat <- na.omit(length.dat) %>%
#   group_by(year, bay) %>%
#   summarize(count=n())
# 
# View(check.dat)
# 
# 
# cc <- length.dat %>%
#   filter(year ==2018 & species=="poll" & bay == "Cooks")
# 
# cc
# 
# site <- read.csv("data/2018 2020 site.csv")
# 
# head(site)
# View(site)

## updating with 2021 WGOA length data
length.2021 <- read.csv("./data/length_2021.csv")

# clean up - select only P. cod and pollock
unique(length.2021$Species)

length.2021 <- length.2021 %>%
  filter(Species %in% c("Pacific cod", "walleye pollock"))

# change species code to match previous data
change <- length.2021$Species == "Pacific cod"
length.2021$Species[change] <- "pcod"

change <- length.2021$Species == "walleye pollock"
length.2021$Species[change] <- "poll"

cod.length <- read.csv("./data/cod length data.csv", row.names = 1)

head(cod.length)

pollock.length <- read.csv("./data/pollock length data.csv", row.names = 1)

head(pollock.length)

old.dat <- rbind(cod.length, pollock.length)

# check site and bay names
unique(length.2021$Site)
unique(old.dat$site)

# shorten "Mitro"
change <- grep("Mitro", length.2021$Site)  
length.2021$Site[change] <- paste(str_sub(length.2021$Site[change], 1,3), str_sub(length.2021$Site[change], 6,7), sep = "")

unique(length.2021$Bay)
unique(old.dat$bay)

# clean up
change <- length.2021$Bay == "Japanese Bay"
length.2021$Bay[change] <- "Japanese"

change <- length.2021$Bay == "Kaiugnak Bay"
length.2021$Bay[change] <- "Kaiugnak"

change <- length.2021$Bay == "Pt Wrangell"
length.2021$Bay[change] <- "Port Wrangell"

change <- length.2021$Bay == "Rodmans Reach"
length.2021$Bay[change] <- "Rodman Reach"

length.2021$Date <- dates(length.2021$Date)
length.2021$julian <- lubridate::yday(length.2021$Date)

length.2021$year <- 2021

# now need to link in with temp and cpue data
site <- read.csv("./data/site_2021.csv")
cpue <- read.csv("./data/cpue_2021.csv")

head(site) 
head(cpue)

# simplify and clean up to join with length
unique(cpue$species)

cpue <- cpue %>%
  filter(species %in% c("Pacific cod", "walleye pollock")) %>%
  select(Station, species, CPUE) %>%
  pivot_wider(names_from = species, values_from = CPUE) %>%
  rename(cod.cpue = `Pacific cod`, pollock.cpue = `walleye pollock`)

# change NA to 0
change <- is.na(cpue)
cpue[change] <- 0

# join to length.2021
length.2021 <- left_join(length.2021, cpue)

# and now temps
site <- site %>%
  select(Station, Temp.C) %>%
  rename(temperature = Temp.C) 

site$Station <- as.numeric(site$Station)

length.2021 <- left_join(length.2021, site)

head(length.2021)
head(old.dat)

# clean up names and put together
length.2021 <- length.2021 %>%
  select(year, julian, Site, Bay, Species, Length, cod.cpue, pollock.cpue, temperature) %>%
  rename(site = Site, bay = Bay, species = Species, length = Length)

new.dat <- rbind(old.dat, length.2021)

unique(new.dat$site)
unique(new.dat$bay)

# and save!
write.csv(new.dat, "./data/cod.pollock.lengths.2006.2021.csv", row.names = F)

########################
# cod condition data
dat <- read.csv("data/Cond_18_20.csv", row.names = 1)

head(dat)

dat <- dat %>%
  filter(age==0) %>%
  select(-age, -location, -date, - month, -year)

names(dat)[1] <- "station"
names(dat)[2] <- "length"

nrow(dat)

# load cpue data
# now load WGOA data
d2 <- read.csv("data/2018 2019 catch.csv")

# combine with complete 2020 catch data!
d3 <- read.csv("data/all_cpue_2020.csv")

# add a measured column
d3$Measured <- NA

unique(d3$species) # looks good, no different names for cod/pollock

#  combine with d2
d3 <- d3 %>%
  select(Station, species, Measured, CPUE)

names(d3) <- c("Station", "Species", "Measured", "Total.catch")

# combine
d2 <- rbind(d2, d3)

# check for repeat spp names
unique(d2$Species)[order(unique(d2$Species))]

# pollock and cod are clean - no different names!
cod <- d2 %>%
  filter(Species=="Pacific cod")
cod <- data.frame(Station=cod$Station, cod=cod$Total.catch)

wgoa.dat <- data.frame(Station=unique(d2$Station))

wgoa.dat <- left_join(wgoa.dat, cod)

names(wgoa.dat) <- c("station", "cod.cpue")
wgoa.dat$station <- as.factor(wgoa.dat$station)
head(dat)

dat <- left_join(dat, wgoa.dat)

# now need to add year, julian day, site, and bay!
d4 <- read.csv("data/2018 2020 site.csv")

head(d4)

# retain only the keeper sets

## NB! changing 115 to cpue==no!
d4$use.for.CPUE[d4$Station==115] <- "no"

d4 <- d4 %>%
  filter(use.for.CPUE=="yes")

d4 <- d4 %>%
  select(Date, Station, Site, Bay, Temp.C)

# calculate Julian day
d4$Date <- dates(as.character(d4$Date))
d4$julian <- lubridate::yday(d4$Date)
d4$year <- years(d4$Date)

head(d4)

d4 <- d4 %>%
  select(-Date)

names(d4)[1] <- "station"
d4$station <- as.factor(d4$station)

all.dat <- left_join(dat, d4)

names(all.dat)[12:14] <- c("site", "bay", "temperature")

str(all.dat)
all.dat$year <- as.numeric(as.character(all.dat$year))

# note that there are two data entry errors

fix <- all.dat$station==203 & all.dat$HSI.dry > 600
all.dat[fix,]
all.dat[fix, "HSI.dry"] <- 7.57

fix <- all.dat$station==260 & all.dat$HSI.dry > 70
all.dat[fix,]
all.dat[fix, "HSI.dry"] <- 7.79

ggplot(all.dat, aes(HSI.dry)) +
  geom_histogram(bins = 100) +
  theme_bw()

write.csv(all.dat, "data/condition_data.csv")
