# add 2021 data to models of cod and pollock abundance

library(tidyverse)

# load cod/pollock data for 2006-2020
d1 <- read.csv("./data/cpue.data.csv")
head(d1)

# load 2021 cpue data
d2 <- read.csv("./data/cpue_2021.csv")
head(d2)
str(d2)

# change years that were incorrectly entered as 2020 - enter as 2021
d2$year <- 2021

# change date to julian
d2$julian <- lubridate::yday(as.POSIXct(paste(d2$year, d2$month, d2$day, sep = "-")))

# examine spp. names 
unique(d2$species)

# total fish caught
sum(d2$CPUE)


# need to remove age-1 cod
d3 <- read.csv("./data/length_2021.csv")
head(d3)
unique(d3$Species)
hist(filter(d3, Species == "Pacific cod")$Length)

# clear break with age 1 >> 150mm

age.1 <- d3 %>%
  filter(Species == "Pacific cod", Length > 150) %>%
  group_by(Site) %>%
  summarize(cod.age.1 = n())
age.1

# clean up names to match d1
names(d2)[5:6] <- c("site", "bay")

# restrict d2 to cod and pollock
d2 <- d2 %>%
  filter(species %in% c("Pacific cod", "walleye pollock")) %>%
  select(Station, year, bay, site, julian, species, CPUE) %>%
  pivot_wider(names_from = species, values_from = CPUE)

# remove age-1 cod
names(age.1)[1] <- "site"

unique(d2$site)
unique(age.1$site)

d2 <- left_join(d2, age.1)

# replace NAs with 0
change <- is.na(d2)
d2[change] <- 0

# clean up again
d2 <- d2 %>%
  mutate(cod.age.0 = `Pacific cod` - cod.age.1) %>%
  mutate(pollock.age.0 = `walleye pollock`) %>%
  select(-`Pacific cod`, -cod.age.1, -`walleye pollock`, -year, -julian, -bay, -site)

# now we need to join to site data to account for sets with no cod or pollock caught
d4 <- read.csv("./data/site_2021.csv")

# change years that were incorrectly entered as 2020 - enter as 2021
d4$year <- 2021

# change date to julian
d4$julian <- lubridate::yday(as.POSIXct(paste(d4$year, d4$month, d4$day, sep = "-")))

# clean up 
d4 <- d4 %>%
  filter(use.for.CPUE == "yes") %>%
  select(Station, year, Bay, Site, julian) %>%
  mutate(Station = as.integer(as.character(Station)))

d4 <- left_join(d4, d2)

# replace NA with 0
change <- is.na(d4)
d4[change] <- 0

# finally, clean up d1/d4 and combine
d1 <- d1 %>%
  select(-temperature)

d4 <- d4 %>%
  select(-Station)

names(d4)[2:3] <- c("bay", "site")

dat <- rbind(d1, d4)

# hot dog

## add 2021 Cooks / Anton's data------------------------------------------------
d5 <- read.csv("./data/Kodiak gadid CPUE 2006-2021.csv")

head(d5)

# clean up to combine with dat
# change to julian day
d5$julian <- lubridate::yday(as.POSIXct(d5$Date, format = "%m/%d/%Y"))

d5 <- d5 %>% 
  filter(Year == 2021) %>%
  select(Year, Region, Site.Name, julian, Pacific.cod, Pollock)

# reset names
names(d5) <- names(dat)

# check for uniform names
unique(d5$bay)
unique(dat$bay)


check.sites <- data.frame(dat = c(str_sort(site.dat), NA),
                          d5 = str_sort(site.d5))

# aha! d5 has "Middle cove" and "Middle Cove"
change <- d5$site == "Middle cove"
d5$site[change] <- "Middle Cove"

# check that worked
site.dat <- unique(filter(dat, bay %in% c("Cook Bay", "Anton Larson Bay"))$site)
site.d5 <-  unique(d5$site)
check.sites <- data.frame(dat = str_sort(site.dat),
                          d5 = str_sort(site.d5))
check.sites

# hot dog

# combine
dat <- rbind(dat, d5)

# check
g <- ggplot(dat) +
  aes(x = year, y = cod.age.0, color = site) +
  geom_point() +
  facet_wrap( ~ bay) +
  theme(legend.position = "none")
print(g)

# some messy names

# first, remove Caton Harbor and NE Harbor (Sanak)
dat <- dat %>%
  filter(bay != "Caton Harbor", bay != "NE Harbor")

# clean up names!

change <- dat$bay == "Cook Bay"
dat$bay[change] <- "Cook"

change <- dat$bay == "Anton Larson Bay"
dat$bay[change] <- "Anton Larson"

# see what's going on with Agripina
unique(dat$bay) # the famous trailing space!

change <- dat$bay == "Agripina "
dat$bay[change] <- "Agripina"

change <- dat$bay == "Japanese Bay"
dat$bay[change] <- "Japanese"

change <- dat$bay == "Kaiugnak Bay"
dat$bay[change] <- "Kaiugnak"

change <- dat$bay == "Pt Wrangell"
dat$bay[change] <- "Port Wrangell"

change <- dat$bay == "Rodmans Reach"
dat$bay[change] <- "Rodman Reach"

# check again
g <- ggplot(dat) +
  aes(x = year, y = cod.age.0, color = site) +
  geom_point() +
  facet_wrap( ~ bay) +
  theme(legend.position = "none")
print(g)

# bays are good - now check sites
str_sort(unique(dat$site))

# looks like some repeat Rodman Reach names
View(filter(dat, bay == "Rodman Reach"))
# yes - the site names are different in 2021

change <- dat$site == "Rod-1"
dat$site[change] <- "RR-1"

change <- dat$site == "Rod-2"
dat$site[change] <- "RR-2"

change <- dat$site == "Rod-4"
dat$site[change] <- "RR-4"

change <- dat$site == "Rod-5"
dat$site[change] <- "RR-5"

change <- dat$site == "Rod-6"
dat$site[change] <- "RR-6"

# check again
str_sort(unique(dat$site))

filter(dat, site == "Rod-1")

# save
write.csv(dat, "./data/age.0_cod_pollock_seine_cpue.csv")