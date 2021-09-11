## estimate pink salmon abundance by year for WGOA sampling

library(tidyverse)
library(chron)

#  load WGOA data
d2 <- read.csv("data/2018 2019 catch.csv")

# combine with complete 2020 catch data!
d3 <- read.csv("data/all_cpue_2020.csv")

# add a measured column
d3$Measured <- NA

unique(d2$Species)
unique(d3$species) # looks good

#  combine with d2
d3 <- d3 %>%
  select(Station, species, Measured, CPUE)

names(d3) <- c("Station", "Species", "Measured", "Total.catch")

names(d2)

# combine
dat <- rbind(d2, d3)

dat <- dat %>%
  filter(Species=="pink salmon")

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

names(d4)[1] <- "Station"
d4$Station <- as.factor(d4$Station)

dat$Station <- as.factor(dat$Station)

d4 <- left_join(d4, dat)

names(d4)[1:3] <- c("station", "site", "bay")
names(d4)[9] <- "pink"

d4 <- d4 %>%
  select(year, julian, site, bay, pink)

# drop Cooks (only sampled once!)
d4 <- d4 %>% 
  filter(bay != "Cooks") %>%
  filter(bay != "Kujulik")
d4[d4$bay == "Cooks",]

# set regional codes
unique(d4$bay)
d4$region <- ifelse(d4$bay %in% c("Ugak", "Kiluida", "Kaiugnak", "Japanese", "Rodman Reach"),
                    "East side",
                    "Peninsula")
pink.dat <- d4

change <- is.na(pink.dat$pink)
pink.dat$pink[change] <- 0

pink.dat$log.pink <- log(pink.dat$pink+1, 10)

summary <- pink.dat %>%
  group_by(year, region) %>%
  summarise(mean=mean(log.pink), se=sd(log.pink)/sqrt(n()), n=n())

## now estimate annual abundance

## Model up seine recruitment estimates for each year
## and compare with stock assessment model esimated recruitment

library(ggplot2)
library(plyr)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
source("./scripts/stan_utils.R")


## set up data --------------------------------------------
pink.dat$bay_fac <- as.factor(pink.dat$bay)
pink.dat$year_fac <- as.factor(pink.dat$year)
pink.dat$region_fac <- as.factor(pink.dat$region)
change <- is.na(pink.dat$pink)
pink.dat$pink[change] <- 0

pink.dat <- pink.dat %>%
  filter(region == "East side")
## brms: setup ---------------------------------------------

## Define model formulas
pink_formula <-  bf(pink ~ julian + year_fac,
                      zi ~ julian + year_fac)


## Set model distributions
zinb <- zero_inflated_negbinomial(link = "log", link_shape = "log", link_zi = "logit")

## fit: zero-inflated --------------------------------------
pink_zinb <- brm(pink_formula,
                   data = pink.dat,
                   family = zinb,
                   cores = 4, chains = 4, iter = 4000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.99, max_treedepth = 10))
saveRDS(pink_zinb, file = "output/pink_zinb.rds")
pink_zinb  <- add_criterion(pink_zinb,c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(pink_zinb, file = "output/pink_zinb.rds")

pink_zinb <- readRDS("./output/pink_zinb.rds")
check_hmc_diagnostics(pink_zinb$fit)
neff_lowest(pink_zinb$fit)
rhat_highest(pink_zinb$fit)
summary(pink_zinb)
bayes_R2(pink_zinb)
plot(pink_zinb$criteria$loo, "k")

plot(pink_zinb, ask = FALSE)
y <- pink.dat$pink
yrep_pink_zinb  <- fitted(pink_zinb, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_pink_zinb[sample(nrow(yrep_pink_zinb), 25), ]) +
  xlim(0, 500) +
  ggtitle("pink_zinb")
pdf("./figs/trace_pink_zinb.pdf", width = 6, height = 4)
trace_plot(pink_zinb$fit)
dev.off()

## Predicted effects ---------------------------------------

## Year predictions ##

## 95% CI
ce1s_1 <- conditional_effects(pink_zinb, probs = c(0.025, 0.975))
mod.95 <- ce1s_1$year_fac %>%
  select(year_fac, estimate__, lower__, upper__)
names(mod.95)[3:4] <- c("ymin.95", "ymax.95")
mod.95$year <- as.numeric(as.character(mod.95$year_fac))

theme_set(theme_bw())

g1 <- ggplot(mod.95) +
  aes(x = year, y = estimate__) +
  geom_ribbon(aes(ymin = ymin.95, ymax = ymax.95), fill = "grey90", alpha=0.8) +
  # geom_ribbon(aes(ymin = ymin.90, ymax = ymax.90), fill = "grey85") +
  # geom_ribbon(aes(ymin = ymin.80, ymax = ymax.80), fill = "grey80") +
  geom_line(size = 0.5, color = "red3") +
  geom_hline(yintercept = 0.9, lty=2) +
  theme(axis.title.x = element_blank()) +
  ylab("CPUE") +
  scale_x_continuous(breaks=seq(1980, 2040, 10)) 

print(g1)
