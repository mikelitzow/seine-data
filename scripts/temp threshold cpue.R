library(tidyverse)
library(mgcv)
library(chron)

theme_set(theme_bw())
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# test for a threshold temperature effect on cod CPUE

# load data

# first, long-term Kodiak sites

dat <- read.csv("data/Kodiak.seine.data.csv")

head(dat)

names(dat)[1] <- "year"

dat$Date <- dates(as.character(dat$Date))

dat$julian <- lubridate::yday(dat$Date)
# calculate Julian day

# dat <- dat %>% 
#   separate(Date, c("month", "day", "drop"), "/") %>%
#   select(-drop)
# 
# 
# dat$month <- as.numeric(dat$month)
# dat$day <- as.numeric(dat$day)
# 
# # dat$year <- ifelse(dat$year<2000, dat$year+2000, dat$year)
# dat$julian <- NA
# for(i in 1:nrow(dat)){
#   
# dat$julian[i] <- julian(dat$month[i], dat$day[i], dat$year[i], origin. = c(1,1,dat$Year[i]))}
# 

dat$cod <- log(dat$X..Pacific.cod+1)
dat$pollock <- log(dat$X..Pollock+1)

head(dat)

dat <- dat %>%
  select(year, julian, Site.Name, Region, cod, pollock)

head(dat)
names(dat)[3:4] <- c("site", "bay")

# now load WGOA data
d2 <- read.csv("data/2018 2019 catch.csv")

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

# log transform
wgoa.dat$cod <- log(wgoa.dat$cod+1)
wgoa.dat$pollock <- log(wgoa.dat$pollock+1)

# now need to add year, juliabn day, site, and bay!
d3 <- read.csv("data/2018 2019 site.csv")

head(d3)

# retain only the keeper sets

## NB! changing 115 to cpue=="no"!
d3$use.for.CPUE[d3$Station==115] <- "no"

d3 <- d3 %>%
  filter(use.for.CPUE=="yes")

d3 <- d3 %>%
  select(Date, Station, Site, Bay)

# calculate Julian day
d3$Date <- dates(as.character(d3$Date))
d3$julian <- lubridate::yday(d3$Date)
d3$year <- years(d3$Date)

# d3 <- d3 %>% 
#   separate(Date, c("month", "day", "year"), "/")
# 
# d3$month <- as.numeric(d3$month)
# d3$day <- as.numeric(d3$day)
# d3$year <- as.numeric(d3$year)+2000
# 
# d3$julian <- julian(d3$month, d3$day, d3$year, origin. = c(1,1,d3$year))

head(d3)

d3 <- d3 %>%
  select(-Date)

d3 <- left_join(d3, wgoa.dat)

names(d3)[2:3] <- c("site", "bay")

# tidyverse not working here for some reason, don't have interet to figure out why!
# d3 <- d3 %>%
#   select(year, julian, site, bay, cod, pollock)

d3 <- d3[,c(5,4,2,3,6,7)]

dat <- rbind(dat, d3)
dat$year <- as.numeric(dat$year)
# now add SST (placeholder for temp!)
sst <- read.csv("data/goa.winter.sst.csv")

names(sst)[2] <- "sst"

dat <- left_join(dat, sst)

# now we can analyze...

# quick check for outliers

# (need to add that when we're not trying to round a cape!)

# remove Kujulik
dat <- dat %>%
  filter(bay != "Kujulik")

dat$year <- as.factor(dat$year)

# first model - site nested in bay, smooth sst and julian date terms...
mod1 <- gam(cod ~  bay/site  + s(sst, bs='cr', k=5) + s(julian, bs='cr', k=5), data=dat)
# Model has more coefficients than data

# might be a better idea to fit a linear era / threshold term to sst instead of a smooth

# let's look at the time series of temp

plot <- dat %>%
  group_by(year) %>%
  summarise(sst=mean(sst))

plot$year <- as.numeric(as.character(plot$year))
ggplot(plot, aes(year, sst)) +
  theme_bw() +
  geom_line() +
  geom_point()

# add a factor indicating whether the temp is above/below threshold
dat$temp <- ifelse(dat$sst < 6.5, "cold", "hot")

mod2 <- gam(cod ~  bay/site  + temp + s(julian, bs='cr', k=5), data=dat)

# Model has more coefficients than data

# try gamm
library(gamm4)
mod3 <- gamm4(cod ~  bay  + temp + s(julian, bs='cr', k=5), random = ~(site|bay), data=dat)
# Error: number of observations (=968) <= number of random effects (=1425) for term (site | bay); 
# the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable

mod4 <- gamm4(cod ~  bay  + temp + s(julian, bs='cr', k=5), random = ~(site), data=dat)
# fixed-effect model matrix is rank deficient so dropping 15 columns / coefficients
summary(mod4$gam) # R2 = 0.03!

dat$cod <- dat$cod+0.1
mod5 <- gamm4(cod ~  bay  + temp + s(julian, bs='cr', k=5), random = ~(site), data=dat, family="inverse.gaussian")

# fixed-effect model matrix is rank deficient so dropping 15 columns / coefficients
# Error in (function (fr, X, reTrms, family, nAGQ = 1L, verbose = 0L, maxit = 100L,  : 
#                       (maxstephalfit) PIRLS step-halvings failed to reduce deviance in pwrssUpdate
#                     In addition: There were 11 warnings (use warnings() to see them)

mod5 <- gam(cod ~  bay  + temp + s(julian, bs='cr', k=5), data=dat, family="inverse.gaussian")
# Warning messages:
#   1: In sqrt(eta) : NaNs produced
# 2: Step size truncated due to divergence 

ggplot(dat, aes(julian, cod)) +
  theme_bw() +
  geom_point() +
  geom_smooth(method='gam')
  
# so it looks like we have a few late-date sets...might pull those
new.dat <- dat %>%
  filter(julian < 255)

mod6 <- gam(cod ~  bay  + temp + s(julian, bs='cr', k=5), data=new.dat, family="inverse.gaussian")
# Warning messages:
#   1: In sqrt(eta) : NaNs produced
# 2: Step size truncated due to divergence 

mod7 <- gam(cod ~  bay  + temp + s(julian, bs='cr', k=5), data=new.dat, family="gamma")
# Error: $ operator is invalid for atomic vectors

mod8 <- gam(cod ~  bay  + temp + s(julian, bs='cr', k=5), data=new.dat, family=gamma())
# Error in gamma() : 0 arguments passed to 'gamma' which requires 1

mod9 <- gam(cod ~  bay  + temp + s(julian, bs='cr', k=5), data=new.dat, inverse.gaussian(link = "1/mu^2"))
# Warning messages:
#   1: In sqrt(eta) : NaNs produced
# 2: Step size truncated due to divergence 

mod10 <- gamm4(cod ~  bay  + temp + s(julian, bs='cr', k=5), random = ~(site), data=new.dat)
summary(mod10$gam)
# R2=0.02!!

plot.gam(mod10$gam)
plot.gam(mod1, resid=T, pch=1)

mod2 <- gam(cod ~  + s(sst, k=5) + s(bay, year, bs="re") + s(julian, year, bs="re"), data=dat)
summary(mod2)

plot.gam(mod2)
plot.gam(mod2, resid=T, pch=1)


cod.sst <- tapply(dat$cod, dat$sst, mean)

plot(names(cod.sst), cod.sst)
