# add 2021 data to models of cod and pollock abundance

# updated at the end of the script to update with 2022 data
# updated at the end of the script to update with 2023 data
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

### add 2022 data ----------------------------

# load 2022 wGOA cpue data
d6 <- read.csv("./data/cpue2022.csv")

# change date to julian
d6$julian <- lubridate::yday(as.POSIXct(paste(d6$year, d6$month, d6$day, sep = "-")))

# examine spp. names 
unique(d6$species)

# need to check for age-1 cod
d7 <- read.csv("./data/length2022.csv")
head(d7)
unique(d7$species)
hist(filter(d7, species == "Pacific cod")$length)
# none!

hist(filter(d7, species == "walleye pollock")$length)
# one - need to remove anything < 100 mm

age.1 <- d7 %>%
  filter(species == "walleye pollock", length > 100) %>%
  group_by(Site) %>%
  summarize(pollock.age.1 = n())
age.1

# clean up names to match d1
names(d6)[6:7] <- c("site", "bay")

# restrict d6 to cod and pollock
d6 <- d6 %>%
  filter(species %in% c("Pacific cod", "walleye pollock")) %>%
  select(Station, year, bay, site, julian, species, CPUE) %>%
  pivot_wider(names_from = species, values_from = CPUE)

# remove age-1 pollock
names(age.1)[1] <- "site"

unique(d6$site)
unique(age.1$site)

d6 <- left_join(d6, age.1)

# replace NAs with 0
change <- is.na(d6)
d6[change] <- 0

# clean up again
d6 <- d6 %>%
  mutate(cod.age.0 = `Pacific cod`) %>%
  mutate(pollock.age.0 = `walleye pollock` - pollock.age.1) %>%
  select(-`Pacific cod`, -pollock.age.1, -`walleye pollock`, -year, -julian, -bay, -site)

# now we need to join to site data to account for sets with no cod or pollock caught
d8 <- read.csv("./data/site2022.csv")


# change date to julian
d8$julian <- lubridate::yday(as.POSIXct(d8$Date, format = "%m/%d/%Y"))

# clean up 
d8 <- d8 %>%
  mutate(year = 2022) %>%
  filter(use.for.CPUE == "yes") %>%
  select(Station, year, Bay, Site, julian) %>%
  mutate(Station = as.integer(as.character(Station)))

d8 <- left_join(d8, d6)

# replace NA with 0
change <- is.na(d8)
d8[change] <- 0

# remove Chief Cove and May samples
d8 <- d8 %>%
  filter(Bay != "Chief Cove",
         julian > 151)

# clean up 
d8 <- d8 %>%
  select(-Station)

names(d8)[2:3] <- c("bay", "site")

dat <- rbind(dat, d8)

## add 2022 Cook / Anton's data-------------------
d9 <- read.csv("./data/Kodiak 2022 seine data - gadid and salmonid.csv")

head(d9)

# clean up to combine with dat
# change to julian day
d9$julian <- lubridate::yday(as.POSIXct(d9$Date, format = "%m/%d/%Y"))

d9 <- d9 %>%
  arrange(desc(X.1..Pacific.cod))
head(d9)

d9 <- d9 %>% 
  filter(Year == 2022) %>%
  select(Year, Region, Site.Name, julian, X..Pacific.cod, X..Pollock) %>%
  rename(Pacific.cod = X..Pacific.cod,
         Pollock = X..Pollock)

# reset names
names(d9) <- names(dat)

unique(d9$bay)
unique(dat$bay)

change <- d9$bay == "Cook Bay"
d9$bay[change] <- "Cook"


change <- d9$bay == "Anton Larson Bay"
d9$bay[change] <- "Anton Larsen"

# change Anton Larsen spelling in dat
change <- dat$bay == "Anton Larson"
dat$bay[change] <- "Anton Larsen"


unique(d9$site)
unique(dat$site)

change <- d9$site == "laminaria #2"
d9$site[change] <- "Laminaria #2"

change <- d9$site == "Middle cove"
d9$site[change] <- "Middle Cove"

change <- dat$site == "Mitro-1"
dat$site[change] <- "Mit-1"

change <- dat$site == "Mitro-2"
dat$site[change] <- "Mit-2"

change <- dat$site == "Mitro-3"
dat$site[change] <- "Mit-3"

change <- dat$site == "Mitro-4"
dat$site[change] <- "Mit-4"

change <- dat$site == "Mitro-5"
dat$site[change] <- "Mit-5"

change <- dat$site == "Mitro-6"
dat$site[change] <- "Mit-6"

change <- dat$bay == "Kiluida"
dat$bay[change] <- "Kiliuda"

change <- dat$bay == "Pt Wrangell"
dat$bay[change] <- "Port Wrangell"

# check that worked
site.dat <- unique(filter(dat, bay %in% c("Cook", "Anton Larsen"))$site)
site.d9 <-  unique(d9$site)
check.sites <- data.frame(dat = str_sort(site.dat),
                          d9 = str_sort(site.d9))

# combine
dat <- rbind(dat, d9)

# check
g <- ggplot(dat) +
  aes(x = year, y = cod.age.0, color = site) +
  geom_point() +
  facet_wrap( ~ bay) +
  theme(legend.position = "none")
print(g)

### add 2023 data ----------------------------

# load 2023 wGOA cpue data
d10 <- read.csv("./data/cpue2023.csv")

# change date to julian
d10$julian <- lubridate::yday(as.POSIXct(d10$date, format = "%m/%d/%Y"))

# examine spp. names 
unique(d10$species)

# need to check for age-1 cod
d11 <- read.csv("./data/length2023.csv")
head(d11)
unique(d11$species)
hist(filter(d11, species == "Pacific cod")$length)

# clear break with age 1 >= 139mm
#now isolate age-1

age.1 <- d11 %>%
  filter(species == "Pacific cod", length > 120) %>%
  group_by(Site) %>%
  summarize(cod.age.1 = n())
age.1

# clean up column names to match d1
names(d11)[3:4] <- c("site", "bay")
names(d10)[6:7] <- c("site", "bay")

# restrict cpue file, d10, to cod and pollock
d10 <- d10 %>%
  filter(species %in% c("Pacific cod", "walleye pollock")) %>%
  select(Station, year, bay, site, julian, species, CPUE) %>%
  pivot_wider(names_from = species, values_from = CPUE)

# remove age-1 cod
names(age.1)[1] <- "site"

unique(d10$site)
unique(age.1$site)

d10 <- left_join(d10, age.1)

# replace NAs with 0
change <- is.na(d10)
d10[change] <- 0

##check if age-1 pollock
hist(filter(d11, species == "walleye pollock")$length)

#need to remove pollock > 100mm
#first isolate age-1 pollock

age.1p <- d11 %>%
  filter(species == "walleye pollock", length > 100) %>%
  group_by(site) %>%
  summarize(pollock.age.1 = n())
age.1p

# remove age-1 pollock
names(age.1p)[1] <- "site"

unique(d10$site)
unique(age.1p$site)

d10 <- left_join(d10, age.1p)

# replace NAs with 0
change <- is.na(d10)
d10[change] <- 0

# clean up again
d10 <- d10 %>%
  mutate(cod.age.0 = `Pacific cod` - cod.age.1) %>%
  mutate(pollock.age.0 = `walleye pollock` - pollock.age.1) 
head(d10)
##here I removed code that selected to remove year, cod.age.0, bay, etc
##not sure if this will mess up data joining or not

# now we need to join to site data to account for sets with no cod or pollock caught
d12 <- read.csv("./data/site2023.csv")

# change date to julian
d12$julian <- lubridate::yday(as.POSIXct(d12$Date, format = "%m/%d/%Y"))

# clean up 
d12 <- d12 %>%
  mutate(year = 2023) %>%
  filter(use.for.CPUE == "yes") %>%
  select(Station, year, Bay, Site, julian) %>%
  mutate(Station = as.integer(as.character(Station)))

d12 <- left_join(d12, d10)

# replace NA with 0
change <- is.na(d12)
d12[change] <- 0

# remove Chief Cove and May samples
d12 <- d12 %>%
  filter(Bay != "Chief Cove",
         julian > 151)

# clean up 
d12 <- d12 %>%
  select(-Station)

names(d12)[2:3] <- c("bay", "site")
#want to make d12 have only 6 columns so matches 'dat'
d12a <- d12 %>%
  select(year, bay, site, julian, cod.age.0, pollock.age.0)
head(d12a)

dat <- rbind(dat, d12a)
#hot dog!

## add 2023 Cook / Anton's data-------------------
d13 <- read.csv("./data/Kodiak_2023_seine_data_gadid_and_pinks.csv")

head(d13)

# clean up to combine with dat
# change to julian day
d13$julian <- lubridate::yday(as.POSIXct(d13$Date, format = "%m/%d/%Y"))

d13 <- d13 %>%
  select(year, bay, site, julian, cod.age.0, pollock.age.0)
head(d13)

# reset names
names(d13) <- names(dat)

#check for repeat names
unique(d13$bay)
unique(dat$bay)
unique(d13$site)
unique(dat$site)

#none found, so don't need anything like next 2 lines
#change <- d9$bay == "Cook Bay"
#d9$bay[change] <- "Cook"

# check that worked
site.dat <- unique(filter(dat, bay %in% c("Cook", "Anton Larsen"))$site)
site.d13 <-  unique(d13$site)
check.sites <- data.frame(dat = str_sort(site.dat),
                          d13 = str_sort(site.d13))

# combine
dat <- rbind(dat, d13)

# check
g <- ggplot(dat) +
  aes(x = year, y = cod.age.0, color = site) +
  geom_point() +
  facet_wrap( ~ bay) +
  theme(legend.position = "none")
print(g)

# save
write.csv(dat, "./data/age.0_cod_pollock_seine_cpue.csv")
