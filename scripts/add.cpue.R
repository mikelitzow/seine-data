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

# now fit model to update time series
library(dplyr)
library(plyr)
library(tidyverse)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
source("./scripts/stan_utils.R")


## prepare data --------------------------------------------
dat$cod <- dat$cod.age.0
dat$bay_fac <- as.factor(dat$bay)
dat$year_fac <- as.factor(dat$year)
dat$site_fac <- as.factor(dat$site)

## brms: setup ---------------------------------------------

## Define model formula
time.series_formula <-  bf(cod ~ year_fac + s(julian, k = 4) + (1 | bay_fac/site_fac),
                       zi ~ year_fac + s(julian, k = 4) + (1 | bay_fac/site_fac))

## Set model distributions
zinb <- zero_inflated_negbinomial(link = "log", link_shape = "log", link_zi = "logit")

## Set priors
priors_zinb <- c(set_prior("normal(0, 3)", class = "b"),
                 set_prior("normal(0, 3)", class = "Intercept"),
                 set_prior("student_t(3, 0, 3)", class = "sd"),
                 set_prior("student_t(3, 0, 3)", class = "sds"),
                 set_prior("gamma(0.01, 0.01)", class = "shape"),
                 set_prior("normal(0, 3)", class = "b", dpar = "zi"),
                 set_prior("logistic(0, 1)", class = "Intercept", dpar = "zi"),
                 set_prior("student_t(3, 0, 3)", class = "sd", dpar = "zi"),
                 set_prior("student_t(3, 0, 3)", class = "sds", dpar = "zi"))


## fit: zero-inflated --------------------------------------
cod_time.series_zinb <- brm(time.series_formula,
                    data = dat,
                    prior = priors_zinb,
                    family = zinb,
                    cores = 4, chains = 4, iter = 4000,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.999, max_treedepth = 10))
cod_time.series_zinb  <- add_criterion(cod_time.series_zinb, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_time.series_zinb, file = "output/cod_time.series_zinb.rds")

cod_time.series_zinb <- readRDS("./output/cod_time.series_zinb.rds")
check_hmc_diagnostics(cod_time.series_zinb$fit)
neff_lowest(cod_time.series_zinb$fit)
rhat_highest(cod_time.series_zinb$fit)
summary(cod_time.series_zinb)
bayes_R2(cod_time.series_zinb)
plot(cod_time.series_zinb$criteria$loo, "k")
plot(conditional_smooths(cod_time.series_zinb), ask = FALSE)
y <- cod.data$cod
yrep_cod_time.series_zinb  <- fitted(cod_time.series_zinb, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_time.series_zinb[sample(nrow(yrep_cod_time.series_zinb), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_time.series_zinb")
pdf("./figs/trace_cod_time.series_zinb.pdf", width = 6, height = 4)
trace_plot(cod_time.series_zinb$fit)
dev.off()

## Predicted effects ---------------------------------------

## 95% CI
ce1s_1 <- conditional_effects(cod_far_zinb, effect = "far_fac", re_formula = NA,
                              probs = c(0.025, 0.975))  

plot <- ce1s_1$far_fac %>%
  select(far_fac, estimate__, lower__, upper__)

plot$far_fac <- reorder(plot$far_fac, desc(plot$far_fac))

fig.2b <- ggplot(plot, aes(far_fac, estimate__)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Fish / set") +
  xlab("FAR") +
  scale_x_discrete(labels=c(expression("<0.98"), expression("">=0.98))) +
  scale_y_continuous(breaks=c(1,5,10,50,100,150)) +
  coord_trans(y = "pseudo_log") + 
  theme_bw()

print(fig.2b)