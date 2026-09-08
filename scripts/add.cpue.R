# add 2021 data to models of cod and pollock abundance

# updated at the end of the script to update with 2022 data
# updated at the end of the script to update with 2023 data
# updated at end of script (around line 580) with 2024 data
# updated at end of script with 2025 data
# updated at end of script with 2026 data
library(tidyverse)

# View() only exists in RStudio; make it a no-op when this file is sourced
# non-interactively (e.g. Rscript) so the script runs end-to-end either way.
if (!interactive()) View <- function(x, ...) invisible(x)

# load cod/pollock data for 2006-2020
d1 <- read.csv("./data/cpue.data.csv")
head(d1)

### QA/QC corrections to the 2006-2020 base data ------------------------------
#
# Three errors found by scripts/qaqc_seine_cpue.R are corrected here rather than
# by hand-editing data/cpue.data.csv, so the raw file stays as received and the
# reasoning travels with the code.
#
# They have to run at this point, before the temperature column is dropped
# further down, because temperature is the only thing in the base file that
# tells two sets at one site on one day apart -- the base file records no time.
# The times themselves come from "Kodiak gadid CPUE 2006-2021.csv", which does.

kod.times <- read.csv("./data/Kodiak gadid CPUE 2006-2021.csv", check.names = FALSE)
kod.times$julian <- lubridate::yday(as.POSIXct(kod.times$Date, format = "%m/%d/%Y"))

## (a) 2009 Cook Bay: all 8 sites were seined twice on julian 196, 5-120 min
##     apart. The second pass is re-dated to julian 197 so the two passes are
##     not treated as one sampling day.
cook09 <- kod.times %>%
  filter(Year == 2009, Region == "Cook Bay", julian == 196) %>%
  transmute(site.key = tolower(`Site Name`),
            temp.key = round(Temperature, 1),
            set.time = as.integer(Time))

# the later set at each site
later09 <- cook09 %>%
  group_by(site.key) %>%
  slice_max(set.time, n = 1) %>%
  ungroup()

change <- with(d1, year == 2009 & bay == "Cook Bay" & julian == 196 &
                 paste(tolower(site), round(temperature, 1)) %in%
                 paste(later09$site.key, later09$temp.key))
stopifnot(sum(change) == 8)   # one per site
d1$julian[change] <- 197

## (b) 2010 Cook Bay: Laminaria East appears twice on julian 239 while Middle
##     Cove is missing entirely. The 12:15 set (the one with no temperature) is
##     a mistaken site entry -- it was Middle Cove. Relabelling it restores the
##     full 8-site day. "Middle Cove" is the spelling used in every other year.
lam10 <- which(d1$year == 2010 & d1$bay == "Cook Bay" &
                 d1$julian == 239 & d1$site == "Laminaria East")
mistaken <- lam10[is.na(d1$temperature[lam10])]   # 10:32 set recorded 11.2 C
stopifnot(length(mistaken) == 1)
d1$site[mistaken] <- "Middle Cove"
stopifnot(sum(d1$year == 2010 & d1$bay == "Cook Bay" & d1$julian == 239) == 8,
          !any(duplicated(d1[d1$year == 2010 & d1$bay == "Cook Bay" &
                               d1$julian == 239, "site"])))

## (c) 2012 Anton Larson Bay: Eelgrass patches is recorded twice on julian 195
##     with an identical catch (258 cod, 21 pollock) and no temperature on the
##     second row. An accidental repeat -- drop the second.
dup12 <- which(d1$year == 2012 & d1$bay == "Anton Larson Bay" &
                 d1$site == "Eelgrass patches" & d1$julian == 195)
stopifnot(length(dup12) == 2, is.na(d1$temperature[dup12[2]]))
d1 <- d1[-dup12[2], ]

## (d) 2023 Anton Larsen, Eelgrass patches on julian 198: two sets with clearly
##     different catches but no time recorded anywhere (the 2023 Kodiak file is
##     the only year in that series without a Time column). Left as-is pending
##     clarification from the data provider.

rm(kod.times, cook09, later09, lam10, mistaken, dup12)

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


# NOTE: site.dat / site.d5 used to be referenced here before they were defined,
# which errored on a fresh top-to-bottom run. Defined first now.
site.dat <- unique(filter(dat, bay %in% c("Cook Bay", "Anton Larson Bay"))$site)
site.d5 <-  unique(d5$site)

# small helper: line two name vectors up side by side even when they differ in length
compare.names <- function(a, b, a.name = "dat", b.name = "new") {
  n <- max(length(a), length(b))
  out <- data.frame(str_sort(c(a, rep(NA, n - length(a)))),
                    str_sort(c(b, rep(NA, n - length(b)))))
  names(out) <- c(a.name, b.name)
  out
}

check.sites <- compare.names(site.dat, site.d5, "dat", "d5")

# aha! d5 has "Middle cove" and "Middle Cove"
change <- d5$site == "Middle cove"
d5$site[change] <- "Middle Cove"

# check that worked
site.dat <- unique(filter(dat, bay %in% c("Cook Bay", "Anton Larson Bay"))$site)
site.d5 <-  unique(d5$site)
check.sites <- compare.names(site.dat, site.d5, "dat", "d5")
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
dat$bay[change] <- "Anton Larsen"

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
check.sites <- compare.names(site.dat, site.d9, "dat", "d9")

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

# QA/QC: station number 531 is used twice in this file -- Jap-1 on 7/30 and
# Jap-5 on 7/31. Station numbers are otherwise one-per-set in 2023. The Jap-5
# rows belong to station 529. This does not change the time series (those two
# rows are whitespotted greenling and snake prickleback, no gadids), but it
# stops one station number pointing at two different sites.
change <- d10$Station == 531 & d10$Site == "Jap-5"
stopifnot(sum(change) == 2)
d10$Station[change] <- 529

# change date to julian
d10$julian <- lubridate::yday(as.POSIXct(d10$date, format = "%m/%d/%Y"))

# examine spp. names 
unique(d10$species)

# need to check for age-1 cod
d11 <- read.csv("./data/length2023.csv")
head(d11)
unique(d11$species)
hist(filter(d11, species == "Pacific cod")$length)

#make this into a nice histogram for ESASS talk
cod23len <- filter(d11, species == "Pacific cod")
head(cod23len)
count(cod23len)

library(patchwork)
library(tidyverse)
library(lubridate)
library(ggplot2)

ggplot(cod23len, aes(length, fill =(species))) +
  geom_histogram(binwidth = 10) +
  xlab("Total Length (mm)")+
  xlim(0,250) +
  ylab("Count") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Pacific cod from beach seines in 2023 (n = 1029)")


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
unique(d12$bay)
unique(d12$site)
head(d12)
#there are 2 columns with name bay and 2 columns with name site
#need to rename columns 5 and 6
names(d12)[5:6] <- c("Bay", "Site")
head(d12)

#want to make d12 have only 6 columns so matches 'dat'
d12a <- d12 %>%
  select(year, bay, site, julian, cod.age.0, pollock.age.0)
head(d12a)
unique(d12$site)

dat <- rbind(dat, d12a)
head(dat)
#hot dog!

## add 2023 Cook / Anton's data-------------------
d13 <- read.csv("./data/Kodiak_2023_seine_data_gadid_and_pinks.csv")

# QA/QC: Eelgrass patches is recorded twice on 7/17/2023 while Eelgrass point is
# missing from that day entirely. The set that caught 3 age-0 and 5 age-1 cod
# was Eelgrass point. Relabelling it restores the 4/4 split between the two
# sites that holds in every other year (2023 currently reads 5/3).
#
# NOTE ON SPELLING: the series uses "Eelgrass point" with a lower-case p in all
# 21 years (81 rows). "Eelgrass Point" would create a second, separate site --
# the same split that "Japanese " caused. Lower case is used here deliberately.
change <- d13$Date == "7/17/2023" & d13$site == "Eelgrass patches" &
  d13$cod.age.0 == 3 & d13$cod.age.1 == 5
stopifnot(sum(change) == 1)
d13$site[change] <- "Eelgrass point"

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
check.sites <- compare.names(site.dat, site.d13, "dat", "d13")

# combine
dat <- rbind(dat, d13)
View(dat)
#data look awesome

# check
g <- ggplot(dat) +
  aes(x = year, y = cod.age.0, color = site) +
  geom_point() +
  facet_wrap( ~ bay) +
  theme(legend.position = "none")
print(g)

### add 2024 data ----------------------------

# load 2024 wGOA cpue data
d14 <- read.csv("./data/cpue2024.csv")

# change date to julian
d14$julian <- lubridate::yday(as.POSIXct(d14$date, format = "%m/%d/%Y"))

# examine spp. names 
unique(d14$species)

# need to check for age-1 cod
d15 <- read.csv("./data/length2024.csv")
head(d15)
unique(d15$species)
hist(filter(d15, species == "Pacific cod")$length)

#make this into a nice histogram 
cod24len <- filter(d15, species == "Pacific cod")
head(cod24len)
count(cod24len)


ggplot(cod24len, aes(length, fill =(species))) +
  geom_histogram(binwidth = 10) +
  xlab("Total Length (mm)")+
  xlim(0,250) +
  ylab("Count") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Pacific cod from beach seines in 2024 (n = 910)")


# clear break with age 1 >= 139mm
#now isolate age-1

age.1 <- d15 %>%
  filter(species == "Pacific cod", length > 120) %>%
  group_by(Site) %>%
  summarize(cod.age.1 = n())
age.1

# clean up column names to match d1
names(d15)[3:4] <- c("site", "bay")
names(d14)[6:7] <- c("site", "bay")

# restrict cpue file, d14, to cod and pollock
d14 <- d14 %>%
  filter(species %in% c("Pacific cod", "walleye pollock")) %>%
  select(Station, year, bay, site, julian, species, CPUE) %>%
  pivot_wider(names_from = species, values_from = CPUE)

# remove age-1 cod
names(age.1)[1] <- "site"

unique(d14$site)
unique(age.1$site)

d14 <- left_join(d14, age.1)

# replace NAs with 0
change <- is.na(d14)
d14[change] <- 0

##check if age-1 pollock
hist(filter(d15, species == "walleye pollock")$length)

poll24len <- filter(d15, species == "walleye pollock")
head(poll24len)
count(poll24len)
ggplot(poll24len, aes(length, fill =(species))) +
  geom_histogram(binwidth = 10) +
  xlab("Total Length (mm)")+
  xlim(0,300) +
  ylab("Count") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Walleye pollock from beach seines in 2024 (n = 243)")

#need to remove pollock > 100mm
#first isolate age-1 pollock

age.1p <- d15 %>%
  filter(species == "walleye pollock", length > 100) %>%
  group_by(site) %>%
  summarize(pollock.age.1 = n())
age.1p

# remove age-1 pollock
names(age.1p)[1] <- "site"

unique(d14$site)
unique(age.1p$site)

d14 <- left_join(d14, age.1p)

# replace NAs with 0
change <- is.na(d14)
d14[change] <- 0

# clean up again
d14 <- d14 %>%
  mutate(cod.age.0 = `Pacific cod` - cod.age.1) %>%
  mutate(pollock.age.0 = `walleye pollock` - pollock.age.1) 
head(d14)

# now we need to join to site data to account for sets with no cod or pollock caught
d16 <- read.csv("./data/site2024.csv")

# QA/QC: station 589 (Fox-1, 7/23/2024, 12:10) was a bad tow. It is the only
# retained station in any year with no records at all in its cpue file, so the
# left_join below would enter it as a false 0/0 rather than a missing set.
# Marked "no" here instead of in the raw file, matching the yes/no vocabulary
# site2024.csv already uses.
change <- d16$Station == 589
stopifnot(sum(change) == 1)
d16$use.for.CPUE[change] <- "no"

# change date to julian
d16$julian <- lubridate::yday(as.POSIXct(d16$Date, format = "%m/%d/%Y"))

# clean up 
d16 <- d16 %>%
  mutate(year = 2024) %>%
  filter(use.for.CPUE == "yes") %>%
  select(Station, year, Bay, Site, julian) %>%
  mutate(Station = as.integer(as.character(Station)))

d16 <- left_join(d16, d14)

# replace NA with 0
change <- is.na(d16)
d16[change] <- 0

# remove Chief Cove and May samples
d16 <- d16 %>%
  filter(Bay != "Chief Cove",
         julian > 151)

# clean up 
d16 <- d16 %>%
  select(-Station)
head(d16)

#there are 2 columns with name bay and 2 columns with name site
#need to rename columns 5 and 6
names(d16)[5:6]<- c("BAY", "Site")
head(d16)

names(d16)[2:3] <- c("bay", "site")
unique(d16$bay)
unique(d16$site)
head(d16)

#want to make d16 have only 6 columns so matches 'dat'
d16a <- d16 %>%
  select(year, bay, site, julian, cod.age.0, pollock.age.0)
head(d16a)
unique(d16$site)

dat <- rbind(dat, d16a)
head(dat)
#hot dog!

## add 2024 Cook / Anton's data-------------------
d17 <- read.csv("./data/Kodiak_2024_seine_data_gadid.csv")

head(d17)

# clean up to combine with dat
# change to julian day
d17$julian <- lubridate::yday(as.POSIXct(d17$Date, format = "%m/%d/%Y"))

d17 <- d17 %>%
  select(year, bay, site, julian, cod.age.0, pollock.age.0)
head(d17)

# reset names
names(d17) <- names(dat)

#check for repeat names that are misspelled or dissimmilar
unique(d17$bay)
unique(dat$bay)
unique(d17$site)
unique(dat$site)

#none found, so don't need anything like next 2 lines
#change <- d9$bay == "Cook Bay"
#d9$bay[change] <- "Cook"

# check that worked
site.dat <- unique(filter(dat, bay %in% c("Cook", "Anton Larsen"))$site)
site.d17 <-  unique(d17$site)
check.sites <- compare.names(site.dat, site.d17, "dat", "d17")

# combine
dat <- rbind(dat, d17)
View(dat)
#data look awesome!

# check
g <- ggplot(dat) +
  aes(x = year, y = cod.age.0, color = site) +
  geom_point() +
  facet_wrap( ~ bay) +
  theme(legend.position = "none")
print(g)

# save
write.csv(dat, "./data/age.0_cod_pollock_seine_cpue.csv")

###

### NOTE: a verbatim duplicate of the 2024 block (old lines 792-993) was removed here.
### It re-ran the entire 2024 section a second time, so sourcing this file top-to-bottom
### appended 2024 twice. The saved CSV was never affected because the script was run
### interactively, section by section. Removed so the file can be sourced end-to-end.

#### add 2025 data -------------


# load 2024 wGOA cpue data
d18 <- read.csv("./data/cpue2025.csv")

# change date to julian
d18$julian <- lubridate::yday(as.POSIXct(d18$date, format = "%m/%d/%Y"))

# examine spp. names 
unique(d18$species)

# need to check for age-1 cod
d19 <- read.csv("./data/length2025.csv")
head(d19)
unique(d19$species)
hist(filter(d19, species == "Pacific cod")$length)

#make this into a nice histogram 
cod25len <- filter(d19, species == "Pacific cod")
head(cod25len)
count(cod25len)


ggplot(cod24len, aes(length, fill =(species))) +
  geom_histogram(binwidth = 10) +
  xlab("Total Length (mm)")+
  xlim(0,250) +
  ylab("Count") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Pacific cod from beach seines in 2024 (n = 910)")


# clear break with age 1 >= 139mm
#now isolate age-1

age.1 <- d19 %>%
  filter(species == "Pacific cod", length > 120) %>%
  group_by(Site) %>%
  summarize(cod.age.1 = n())
age.1

# clean up column names to match d1
names(d19)[3:4] <- c("site", "bay")
names(d18)[6:7] <- c("site", "bay")

# restrict cpue file, d18, to cod and pollock
d18 <- d18 %>%
  filter(species %in% c("Pacific cod", "walleye pollock")) %>%
  select(Station, year, bay, site, julian, species, CPUE) %>%
  pivot_wider(names_from = species, values_from = CPUE)

# remove age-1 cod
names(age.1)[1] <- "site"

unique(d18$site)
unique(age.1$site)

d18 <- left_join(d18, age.1)

# replace NAs with 0
change <- is.na(d18)
d18[change] <- 0

##check if age-1 pollock
hist(filter(d19, species == "walleye pollock")$length)

poll25len <- filter(d19, species == "walleye pollock")
head(poll25len)
count(poll25len)
ggplot(poll25len, aes(length, fill =(species))) +
  geom_histogram(binwidth = 10) +
  xlab("Total Length (mm)")+
  xlim(0,300) +
  ylab("Count") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Walleye pollock from beach seines in 2025 (n = 253)")

#need to remove pollock > 100mm
#first isolate age-1 pollock

age.1p <- d19 %>%
  filter(species == "walleye pollock", length > 100) %>%
  group_by(site) %>%
  summarize(pollock.age.1 = n())
age.1p

# remove age-1 pollock
names(age.1p)[1] <- "site"

unique(d18$site)
unique(age.1p$site)

d18 <- left_join(d18, age.1p)

# replace NAs with 0
change <- is.na(d18)
d18[change] <- 0

# clean up again
d18 <- d18 %>%
  mutate(cod.age.0 = `Pacific cod` - cod.age.1) %>%
  mutate(pollock.age.0 = `walleye pollock` - pollock.age.1) 
head(d18)

# now we need to join to site data to account for sets with no cod or pollock caught
d20 <- read.csv("./data/site2025.csv")

# change date to julian
d20$julian <- lubridate::yday(as.POSIXct(d20$Date, format = "%m/%d/%Y"))

# clean up 
d20 <- d20 %>%
  mutate(year = 2025) %>%
  filter(use.for.CPUE == "yes") %>%
  select(Station, year, Bay, Site, julian) %>%
  mutate(Station = as.integer(as.character(Station)))

d20 <- left_join(d20, d18)

# replace NA with 0
change <- is.na(d20)
d20[change] <- 0

# check for May samples
range(d20$julian) # none!

# clean up 
d20 <- d20 %>%
  select(-Station)
head(d20)

#there are 2 columns with name bay and 2 columns with name site
#need to rename columns 5 and 6
names(d20)[5:6]<- c("BAY", "Site")
head(d20)

names(d20)[2:3] <- c("bay", "site")
unique(d20$bay)
unique(d20$site)
head(d20)

#want to make d20 have only 6 columns so matches 'dat'
d20a <- d20 %>%
  select(year, bay, site, julian, cod.age.0, pollock.age.0)
head(d20a)
unique(d20$site)

dat <- rbind(dat, d20a)
head(dat)

#hot dog!

## add 2025 Cook / Anton's data-------------------
d21 <- read.csv("./data/Kodiak_2025_seine_data_gadid.csv")

head(d21)

# clean up to combine with dat
# change to julian day
d21$julian <- lubridate::yday(as.POSIXct(d21$Date, format = "%m/%d/%Y"))

d21 <- d21 %>%
  select(year, bay, site, julian, cod.age.0, pollock.age.0)
head(d21)

# reset names
names(d21) <- names(dat)

#check for repeat names that are misspelled or dissimmilar
unique(d21$bay)
unique(dat$bay)
unique(d21$site)
unique(dat$site)

#none found, so don't need anything like next 2 lines
#change <- d9$bay == "Cook Bay"
#d9$bay[change] <- "Cook"

# check that worked
site.dat <- unique(filter(dat, bay %in% c("Cook", "Anton Larsen"))$site)
site.d21 <-  unique(d21$site)
check.sites <- compare.names(site.dat, site.d21, "dat", "d21")

# combine
dat <- rbind(dat, d21)
View(dat)
#data look awesome!

# check
g <- ggplot(dat) +
  aes(x = year, y = cod.age.0, color = site) +
  geom_point() +
  facet_wrap( ~ bay) +
  theme(legend.position = "none")
print(g)

# save
write.csv(dat, "./data/age.0_cod_pollock_seine_cpue.csv")

#### add 2026 data -------------

# Three things differ from earlier years and are handled explicitly below:
#  1. All four 2026 files carry blank trailing rows (26 in cpue, 2682 in length,
#     20 in site) that must be dropped before anything else.
#  2. site2026.csv reintroduces the trailing space in "Japanese " and contains a
#     transposition typo, "Kiul-4" for station 663. Because the final bay/site
#     columns are taken from the SITE file, both defects would land in the time
#     series -- exactly how "Japanese " got into 2024 and "Kiul-4" into 2025.
#  3. Kiliuda sites are spelled "Kilu-*" in the site/cpue files but "Kil-*" in
#     the length file. The age-1 join used to be keyed on site NAME, so a
#     mismatch there fails silently. It is keyed on Station below instead.

# load 2026 wGOA cpue data (dropping the blank trailing rows)
d22 <- read.csv("./data/cpue2026.csv") %>%
  filter(species != "", !is.na(Station))

# change date to julian
d22$julian <- lubridate::yday(as.POSIXct(d22$date, format = "%m/%d/%Y"))

# examine spp. names
unique(d22$species)

# need to check for age-1 cod
d23 <- read.csv("./data/length2026.csv") %>%
  filter(species != "", !is.na(Station))
head(d23)
unique(d23$species)
hist(filter(d23, species == "Pacific cod")$length)

cod26len <- filter(d23, species == "Pacific cod")
count(cod26len)

ggplot(cod26len, aes(length, fill = (species))) +
  geom_histogram(binwidth = 10) +
  xlab("Total Length (mm)") +
  xlim(0, 250) +
  ylab("Count") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Pacific cod from beach seines in 2026 (n = 1054)")

# age-0 run 26-80 mm, then nothing at all until 190 mm, so the > 120 cut used
# in 2023-2025 sits in a wide empty gap and is unambiguous here.
# keyed on Station, not site name -- see note 3 above
age.1 <- d23 %>%
  filter(species == "Pacific cod", length > 120) %>%
  group_by(Station) %>%
  summarize(cod.age.1 = n())
age.1

# clean up column names to match d1
names(d23)[3:4] <- c("site", "bay")
names(d22)[6:7] <- c("site", "bay")

# restrict cpue file, d22, to cod and pollock
d22 <- d22 %>%
  filter(species %in% c("Pacific cod", "walleye pollock")) %>%
  select(Station, year, bay, site, julian, species, CPUE) %>%
  pivot_wider(names_from = species, values_from = CPUE)

d22 <- left_join(d22, age.1, by = "Station")

# replace NAs with 0
change <- is.na(d22)
d22[change] <- 0

##check if age-1 pollock
hist(filter(d23, species == "walleye pollock")$length)

poll26len <- filter(d23, species == "walleye pollock")
count(poll26len)
ggplot(poll26len, aes(length, fill = (species))) +
  geom_histogram(binwidth = 10) +
  xlab("Total Length (mm)") +
  xlim(0, 300) +
  ylab("Count") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(title = "Walleye pollock from beach seines in 2026 (n = 525)")

# same > 100 mm cut as 2022-2025. The 80-90 mm bin holds 4 fish with empty bins
# on both sides (70-80 and 90-100); they are kept as fast-growing age-0, which
# is what the > 100 cut has always done.
age.1p <- d23 %>%
  filter(species == "walleye pollock", length > 100) %>%
  group_by(Station) %>%
  summarize(pollock.age.1 = n())
age.1p

d22 <- left_join(d22, age.1p, by = "Station")

# replace NAs with 0
change <- is.na(d22)
d22[change] <- 0

# clean up again
d22 <- d22 %>%
  mutate(cod.age.0 = `Pacific cod` - cod.age.1) %>%
  mutate(pollock.age.0 = `walleye pollock` - pollock.age.1)
head(d22)

# now we need to join to site data to account for sets with no cod or pollock caught
d24 <- read.csv("./data/site2026.csv") %>%
  filter(!is.na(Station))

# change date to julian
d24$julian <- lubridate::yday(as.POSIXct(d24$Date, format = "%m/%d/%Y"))

# fix the two defects in the site file before they reach the time series
d24$Bay <- trimws(d24$Bay)    # "Japanese " -> "Japanese"
d24$Site <- trimws(d24$Site)
change <- d24$Site == "Kiul-4"
d24$Site[change] <- "Kilu-4"

# clean up
d24 <- d24 %>%
  mutate(year = 2026) %>%
  filter(use.for.CPUE == "yes") %>%
  select(Station, year, Bay, Site, julian) %>%
  mutate(Station = as.integer(as.character(Station)))

d24 <- left_join(d24, d22)

# replace NA with 0
change <- is.na(d24)
d24[change] <- 0

# check for Chief Cove and May samples
range(d24$julian) # 167-203, so no May samples
d24 <- d24 %>%
  filter(Bay != "Chief Cove",
         julian > 151)

# clean up
d24 <- d24 %>%
  select(-Station)
head(d24)

#there are 2 columns with name bay and 2 columns with name site
#need to rename columns 5 and 6
names(d24)[5:6] <- c("BAY", "Site")
head(d24)

names(d24)[2:3] <- c("bay", "site")
unique(d24$bay)
unique(d24$site)
head(d24)

#want to make d24 have only 6 columns so matches 'dat'
d24a <- d24 %>%
  select(year, bay, site, julian, cod.age.0, pollock.age.0)
head(d24a)

dat <- rbind(dat, d24a)
head(dat)

#hot dog!

## add 2026 Cook / Anton's data-------------------
d25 <- read.csv("./data/Kodiak_2026_seine_data_gadid.csv")

head(d25)

# clean up to combine with dat
# change to julian day
d25$julian <- lubridate::yday(as.POSIXct(d25$Date, format = "%m/%d/%Y"))

d25 <- d25 %>%
  select(year, bay, site, julian, cod.age.0, pollock.age.0)
head(d25)

# reset names
names(d25) <- names(dat)

#check for repeat names that are misspelled or dissimmilar
unique(d25$bay)
unique(dat$bay)
unique(d25$site)
unique(dat$site)

# check that worked
site.dat <- unique(filter(dat, bay %in% c("Cook", "Anton Larsen"))$site)
site.d25 <-  unique(d25$site)
check.sites <- compare.names(site.dat, site.d25, "dat", "d25")

# combine
dat <- rbind(dat, d25)
View(dat)

## repair the same two defects where they already entered the series ----------
# "Japanese " (2024, 7 rows) and "Kiul-4" (2025, 1 row) came in from earlier
# site files by the route described at the top of this section. Trimming the
# whole series also inoculates it against any future recurrence.
dat$bay <- trimws(dat$bay)
dat$site <- trimws(dat$site)

change <- dat$site == "Kiul-4"
dat$site[change] <- "Kilu-4"

# both should now be gone
stopifnot(!any(dat$bay != trimws(dat$bay)))
stopifnot(!any(dat$site == "Kiul-4"))

# no age-0 count should ever be negative (would mean over-subtracted age-1)
stopifnot(all(dat$cod.age.0 >= 0), all(dat$pollock.age.0 >= 0))

## drop bays / sites with no replication --------------------------------------
# BAYS. All four requested bays are already excluded, by three different routes,
# so nothing new is needed here -- only the assertion that they stay out:
#   Canton Harbor  spelled "Caton Harbor" in the data. 2021 only (stations
#                  317-319, sites Sanak-1/2/3). Dropped near the top of this
#                  script by filter(bay != "Caton Harbor").
#   Sanak          the Sanak sites are split over two bay labels; Sanak-4/5
#                  (stations 320-321, 2021 only) sit under bay "NE Harbor",
#                  dropped by the same filter.
#   Chief Cove     2022 only (10 sets: 5 on 22 May, 5 on 2 July). Dropped by
#                  the filter(Bay != "Chief Cove") in each year's block.
#   Kujulik        never sampled. It appears only in bay_lat_long.csv and
#                  seine_temperatures.csv, never in a cpue/site/length file,
#                  so it has never entered the series.
stopifnot(!any(dat$bay %in% c("Caton Harbor", "Canton Harbor", "NE Harbor",
                              "Sanak", "Chief Cove", "Kujulik")))

# SITES. Every site in the series sampled in fewer than three years is dropped as
# unreplicated -- eight sites, 10 sets. Seven are 2018-only; BB-1 is 2018-2019.
# Ugak-1 was seined twice in 2018 (julian 186 and 240), the rest once per year.
# Note on Kilu-7: it was seined again in 2022 (station 368), but on 24 May, so
# the julian > 151 rule already keeps that set out of the series.
# After this the least-replicated remaining site spans four years.
drop.sites <- c("Ugak-1", "Sand-1", "Kilu-9", "Kai-5", "BB-5",
                "Kilu-7", "Kilu-8", "BB-1")

# if a later year ever re-occupies one of these, this fires and that site should
# be taken off the list rather than silently dropped
stopifnot(all(dat$year[dat$site %in% drop.sites] %in% c(2018, 2019)))

change <- dat$site %in% drop.sites
stopifnot(sum(change) == 10)
dat <- dat[!change, ]

# nothing sampled in fewer than three years should survive
reps <- tapply(dat$year, dat$site, function(x) length(unique(x)))
stopifnot(all(reps >= 3))

# NOT dropped: station 546. It is the 2024 wGOA set at Cook / Eelgrass North
# (6 June, julian 158). Eelgrass North is sampled in all 21 years, and the
# early-June wGOA visit to Cook Bay recurs in 2022, 2023, 2025 and 2026, so
# neither the site nor the visit is a 2024 one-off.
stopifnot(n_distinct(dat$year[dat$bay == "Cook" & dat$site == "Eelgrass North"]) > 1)

# check
g <- ggplot(dat) +
  aes(x = year, y = cod.age.0, color = site) +
  geom_point() +
  facet_wrap( ~ bay) +
  theme(legend.position = "none")
print(g)

# save -- eol matches the CRLF line endings the file has always used
write.csv(dat, "./data/age.0_cod_pollock_seine_cpue.csv", eol = "\r\n")

