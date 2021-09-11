library(tidyverse)
library(mgcv)
library(chron)

theme_set(theme_bw())
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# load WGOA data
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
cod <- data.frame(Station=cod$Station, 
                  species="Pacific cod",
                  total=cod$Total.catch,
                  measured=cod$Measured)

pollock <- d2 %>%
  filter(Species=="walleye pollock")
pollock <- data.frame(Station=pollock$Station, 
                      species="walleye pollock",
                      total=pollock$Total.catch, 
                      measured=pollock$Measured)

temp1 <- data.frame(Station=unique(d2$Station))

temp1 <- left_join(temp1, cod)

# fill in NA species as cod
change <- is.na(temp1$species)
temp1$species[change] <- "Pacific cod"

# now pollock
temp2 <- data.frame(Station=unique(d2$Station))
temp2 <- left_join(temp2, pollock)

# fill in NA species as cod
change <- is.na(temp2$species)
temp2$species[change] <- "walleye pollock"

# combine
wgoa.dat <- rbind(temp1, temp2)

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

names(d4)[1] <- names(wgoa.dat)[1] <- "station"

d4 <- left_join(d4, wgoa.dat)

names(d4)[2:4] <- c("site", "bay", "temp.c")

# remove Kujulik
d4 <- d4 %>%
  filter(bay != "Kujulik")


hist(d4$julian)

str(d4)

# now load WGOA length data
wgoa.len <- read.csv("data/gadid_len.csv")

wgoa.len[is.na(wgoa.len$Station), ]

wgoa.len$Species <- as.character(wgoa.len$Species)
  
# restrict to age-0
head(wgoa.len)

hist(wgoa.len$Length, breaks=50)

# get count of fish >= 150 mm
age.1 <- wgoa.len %>%
  filter(Length >= 150) %>%
  group_by(Station, Species) %>%
  summarise(age.1=n())

names(age.1)[1:2] <- c("station", "species") 

d4 <- left_join(d4, age.1)

# replace age.1 NA with 0
change <- is.na(d4$age.1)
d4$age.1[change] <- 0

d4$age.0 <- d4$total-d4$age.1

d4 <- d4 %>%
  select(-total, -measured)

# that's our cpue by age!

# exploratory plots!

## Read in data --------------------------------------------

data <- d4 %>%
  filter(species == "Pacific cod")

data$date <- as.Date(data$julian,
                         origin = paste0(data$year, "-01-01"))

## Explore data --------------------------------------------

## Check distributions

plot(data$age.1)
hist(data$age.1, breaks = 100) ## lots of zeros
tab <- table(data$age.1)
plot(tab)


g <- ggplot(data) +
  aes(x = date, y = age.1, color = site) +
  geom_point() +
  facet_wrap( ~ bay) +
  theme(legend.position = "none")
print(g)

## all looks swell!
write.csv(d4, "./output/wgoa.cod.poll.cpue.csv", row.names = F)
