## Estimate gear temperature by year / day / bay / site
library(ggplot2)
library(plyr)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
library(tidyverse)
source("./scripts/stan_utils.R")


## Read in data --------------------------------------------
dat <- read.csv("./data/seine_temperatures_2006_2021.csv")

dat$bay_fac <- as.factor(dat$bay)
dat$year_fac <- as.factor(dat$year)
dat$site_fac <- as.factor(dat$site)

## Check distributions
plot(dat$temperature)
hist(dat$temperature, breaks = 50) 

g <- ggplot(dat) +
  aes(x = julian, y = temperature) +
  geom_point()
print(g)


## brms: setup ---------------------------------------------

## Define model formulas
temp_formula_1 <-  bf(temperature ~ year_fac + bay_fac)

temp_formula_2 <-  bf(temperature ~ year_fac + bay_fac/site)

temp_formula_3 <-  bf(temperature ~ year_fac + bay_fac + s(julian, k=5))

temp_formula_4 <-  bf(temperature ~ year_fac + bay_fac + s(julian, by=bay_fac, k=3))


## Set model distribution
Gamma <- Gamma(link = "log")

## fit: brms --------------------------------------
temp1_gauss <- brm(temp_formula_1,
                data = dat,
                cores = 4, chains = 4, iter = 4000,
                save_pars = save_pars(all = TRUE),
                control = list(adapt_delta = 0.999, max_treedepth = 10))
temp1_gauss  <- add_criterion(temp1_gauss, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(temp1_gauss, file = "output/temp1_gauss.rds")

temp1_gauss <- readRDS("./output/temp1_gauss.rds")
check_hmc_diagnostics(temp1_gauss$fit)
neff_lowest(temp1_gauss$fit)
rhat_highest(temp1_gauss$fit)
summary(temp1_gauss)
bayes_R2(temp1_gauss)
y <- dat$temperature
yrep_temp1_gauss  <- fitted(temp1_gauss, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_temp1_gauss[sample(nrow(yrep_temp1_gauss), 25), ]) +
  xlim(0, 500) +
  ggtitle("temp1_gauss")


# temp2_gauss <- brm(temp_formula_2,
#                    data = dat,
#                    cores = 4, chains = 4, iter = 3000,
#                    save_pars = save_pars(all = TRUE),
#                    control = list(adapt_delta = 0.9999, max_treedepth = 15))
# temp2_gauss  <- add_criterion(temp2_gauss, c("loo", "bayes_R2"), moment_match = TRUE)
# saveRDS(temp2_gauss, file = "output/temp2_gauss.rds")
# 
# temp2_gauss <- readRDS("./output/temp2_gauss.rds")
# check_hmc_diagnostics(temp2_gauss$fit)
# neff_lowest(temp2_gauss$fit)
# rhat_highest(temp2_gauss$fit)
# summary(temp2_gauss)
# bayes_R2(temp2_gauss)
# y <- dat$temperature
# yrep_temp2_gauss  <- fitted(temp2_gauss, scale = "response", summary = FALSE)
# ppc_dens_overlay(y = y, yrep = yrep_temp2_gauss[sample(nrow(yrep_temp2_gauss), 25), ]) +
#   xlim(0, 500) +
#   ggtitle("temp2_gauss")


temp3_gauss <- brm(temp_formula_3,
                   data = dat,
                   cores = 4, chains = 4, iter = 4000,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 12))
temp3_gauss  <- add_criterion(temp3_gauss, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(temp3_gauss, file = "output/temp3_gauss.rds")

temp3_gauss <- readRDS("./output/temp3_gauss.rds")
check_hmc_diagnostics(temp3_gauss$fit)
neff_lowest(temp3_gauss$fit)
rhat_highest(temp3_gauss$fit)
summary(temp3_gauss)
bayes_R2(temp3_gauss)
y <- as.vector(na.omit(dat$temperature))
yrep_temp3_gauss  <- fitted(temp3_gauss, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_temp3_gauss[sample(nrow(yrep_temp3_gauss), 25), ])  +
  ggtitle("temp3_gauss")

# this model doesn't fit well!
# temp4_gauss <- brm(temp_formula_4,
#                    data = dat,
#                    cores = 4, chains = 4, iter = 6000,
#                    save_pars = save_pars(all = TRUE),
#                    control = list(adapt_delta = 0.999, max_treedepth = 12))
# temp4_gauss  <- add_criterion(temp4_gauss, c("loo", "bayes_R2"), moment_match = TRUE)
# saveRDS(temp4_gauss, file = "output/temp4_gauss.rds")
# 
# temp4_gauss <- readRDS("./output/temp4_gauss.rds")
# check_hmc_diagnostics(temp4_gauss$fit)
# neff_lowest(temp4_gauss$fit)
# rhat_highest(temp4_gauss$fit)
# summary(temp4_gauss)
# bayes_R2(temp4_gauss)
# y <- dat$temperature
# yrep_temp4_gauss  <- fitted(temp4_gauss, scale = "response", summary = FALSE)
# ppc_dens_overlay(y = y, yrep = yrep_temp4_gauss[sample(nrow(yrep_temp4_gauss), 25), ]) +
#   xlim(0, 500) +
#   ggtitle("temp4_gauss")

loo(temp1_gauss, temp3_gauss)


## Predicted effects ---------------------------------------

## Year predictions ##

## 95% CI
ce1s_1 <- conditional_effects(temp3_gauss, probs = c(0.025, 0.975))
temp.95 <- ce1s_1$year_fac %>%
  select(year_fac, estimate__, lower__, upper__)
names(temp.95)[3:4] <- c("ymin.95", "ymax.95")

## 90% CI
ce1s_2 <- conditional_effects(temp3_gauss, probs = c(0.05, 0.95))
temp.90 <- ce1s_2$year_fac %>%
  select(year_fac, estimate__, lower__, upper__)
names(temp.90)[3:4] <- c("ymin.90", "ymax.90")

## 80% CI
ce1s_3 <- conditional_effects(temp3_gauss, probs = c(0.1, 0.9))
temp.80 <- ce1s_3$year_fac %>%
  select(year_fac, estimate__, lower__, upper__)
names(temp.80)[3:4] <- c("ymin.80", "ymax.80")


pred.mod <- left_join(temp.95, temp.90)
pred.mod <- left_join(pred.mod, temp.80)
pred.mod$year <- as.numeric(as.character(pred.mod$year_fac))

theme_set(theme_bw())

g1 <- ggplot(pred.mod) +
  aes(x = year, y = estimate__) +
  geom_errorbar(aes(ymin = ymin.95, ymax = ymax.95)) +
  geom_point() +
  geom_line(size = 0.5) +
  theme(axis.title.x = element_blank()) +
  ylab("Degrees C") +
  scale_x_continuous(breaks=seq(1980, 2040, 10)) 

print(g1)

# ggsave("./figs/year_predicted_effect_mod_far.png", width = 4.5, height = 2)

# save table version
names(temp.95) <- c("year", "temp", "LCI_95%", "UCI_95%")

write.csv(temp.95, "./output/seine_temp_estimates_2006-2021.csv", row.names = F)
