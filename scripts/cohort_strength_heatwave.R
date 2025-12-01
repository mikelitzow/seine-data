# fit brms models to update time series
library(plyr)
library(dplyr)
library(tidyverse)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)

source("./scripts/stan_utils.R")

theme_set(theme_bw())

dat <- read.csv("./data/age.0_cod_pollock_seine_cpue.csv", row.names = 1)

## prepare cod data --------------------------------------------
dat$cod <- dat$cod.age.0
dat$bay_fac <- as.factor(dat$bay)
dat$year_fac <- as.factor(dat$year)
dat$site_fac <- as.factor(dat$site)

# define heatwave years following Nick Strait
dat <- dat %>%
  mutate(heatwave_fac = 
           as.factor(case_when(
             year < 2015 ~ "before",
             year %in% c(2015:2017, 2019) ~ "during",
             year %in% c(2018, 2020) ~ "adjacent", 
             year > 2020 ~ "after"
           )))

# and filter for specified bays
unique(dat$bay_fac)

keep_bays <- c("Anton Larsen",
               "Cook",
               "Kaiugnak",
               "Agripina",
               "Falmouth",
               "Baralof",
               "Sand Point")

dat_reduced <- dat %>%
  filter(bay %in% keep_bays)

# check
unique(dat_reduced$bay_fac)

# check for site-bay-year combinations
check <- dat_reduced %>%
  group_by(year_fac, bay_fac, site_fac) %>%
  summarise(n = n())

View(check)

# and define region codes used in manuscript
dat_reduced <- dat_reduced %>%
  mutate(region_fac = 
           as.factor(case_when(
             bay_fac == "Agripina" ~ "Agripina",
             bay_fac %in% c("Anton Larsen", "Cook") ~ "Eastern Kodiak",
             bay_fac == "Kaiugnak" ~ "Kaiugnak", 
             bay_fac %in% c("Falmouth", "Baralof", "Sand Point") ~ "Shumagin Islands"
           )))

# check for bay-region-year combinations
check <- dat_reduced %>%
  group_by(year_fac, region_fac, bay_fac) %>%
  summarise(n = n())

View(check) # looks good

## cod brms: setup ---------------------------------------------

## Begin with a global model including heatwave status and region
## Excluding year as a group-level effect because it is conflated with heatwave status

## Define model formula
time.series_formula <-  bf(cod ~ heatwave_fac + region_fac + s(julian, k = 4) + (1 | bay_fac/site_fac),
                           zi ~ heatwave_fac + region_fac + s(julian, k = 4) + (1 | bay_fac/site_fac))

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


## cod fit: zero-inflated --------------------------------------
global_mhw_region_zinb <- brm(time.series_formula,
                            data = dat_reduced,
                            prior = priors_zinb,
                            family = zinb,
                            cores = 4, chains = 4, iter = 4000,
                            save_pars = save_pars(all = TRUE),
                            control = list(adapt_delta = 0.999, max_treedepth = 11))

#cod_time.series_zinb  <- add_criterion(cod_time.series_zinb, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(global_mhw_region_zinb, file = "output/global_mhw_region_zinb.rds")

global_mhw_region_zinb <- readRDS("./output/global_mhw_region_zinb.rds")
check_hmc_diagnostics(global_mhw_region_zinb$fit)
neff_lowest(global_mhw_region_zinb$fit)
rhat_highest(global_mhw_region_zinb$fit)
summary(global_mhw_region_zinb)
bayes_R2(global_mhw_region_zinb)


y <- dat_reduced$cod
yrep_global_mhw_region_zinb  <- fitted(global_mhw_region_zinb, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_global_mhw_region_zinb[sample(nrow(yrep_global_mhw_region_zinb), 25), ]) +
  xlim(0, 500) +
  ggtitle("global_mhw_region_zinb")

trace_plot(global_mhw_region_zinb$fit)


## Cod predicted effects ---------------------------------------
## region estimates
## 95% CI
ce1s_1 <- conditional_effects(global_mhw_region_zinb, effect = "region_fac", re_formula = NA,
                              probs = c(0.025, 0.975))  

plot.cod <- ce1s_1$region_fac %>%
  select(region_fac, estimate__, se__, lower__, upper__)

# reorder heatwave status for plotting
plot.cod <- plot.cod %>%
  mutate(order = case_when(
  region_fac == "Shumagin Islands" ~ 1,
  region_fac == "Agripina" ~ 2,
  region_fac == "Kaiugnak" ~ 3,
  region_fac == "Eastern Kodiak" ~ 4),
  region_fac = reorder(region_fac, order))

ggplot(plot.cod, aes(region_fac, estimate__)) +
  geom_col(color = "black", fill = "grey") +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Age-0 cod / set") +
  scale_y_continuous(breaks=c(1,5,10,50,100,200,300), minor_breaks = NULL) +
  coord_trans(y = "pseudo_log") +
  xlab("Heatwave status")

ggsave("./figs/seine_cod_age0_abundance_region_global_model.png", width = 4, height = 6, units = 'in')

## 95% CI
ce1s_1 <- conditional_effects(global_mhw_region_zinb, effect = "heatwave_fac", re_formula = NA,
                              probs = c(0.025, 0.975))  

plot.cod <- ce1s_1$heatwave_fac %>%
  select(heatwave_fac, estimate__, se__, lower__, upper__)

# reorder heatwave status for plotting
plot.cod <- plot.cod %>%
  mutate(order = case_when(
    heatwave_fac == "before" ~ 1,
    heatwave_fac == "during" ~ 2,
    heatwave_fac == "adjacent" ~ 3,
    heatwave_fac == "after" ~ 4),
    heatwave_fac = reorder(heatwave_fac, order))

ggplot(plot.cod, aes(heatwave_fac, estimate__)) +
  geom_col(color = "black", fill = "grey") +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Age-0 cod / set") +
  scale_y_continuous(breaks=c(1,5,10,50,100,200,300), minor_breaks = NULL) +
  coord_trans(y = "pseudo_log") +
  xlab("Heatwave status")

ggsave("./figs/seine_cod_age0_abundance_heatwave_global_model.png", width = 4, height = 6, units = 'in')


# 
# 
# # round, rename columns, and save
# plot.cod <- plot.cod[,1:5]
# plot.cod[,2:5] <- round(plot.cod[,2:5], 2)
# names(plot.cod) <- c("heatwave status", "cod_per_set", "cod_se", "cod_95percent_LCI", "cod_95percent_UCI")
# write.csv(plot.cod, "./output/seine_cod_age0_abundance_heatwave_restricted_bays.csv", row.names = F)

### model heatwave effects for Eastern Kodiak only

## limit data to Cook and Anton Larsen

dat_temp <- dat_reduced %>%
  filter(region_fac == "Eastern Kodiak")


## Define model formula
time.series_formula <-  bf(cod ~ heatwave_fac + s(julian, k = 4) + (1 | bay_fac/site_fac),
                           zi ~ heatwave_fac + s(julian, k = 4) + (1 | bay_fac/site_fac))

## cod fit: zero-inflated --------------------------------------
mhw_zinb <- brm(time.series_formula,
                            data = dat_temp,
                            prior = priors_zinb,
                            family = zinb,
                            cores = 4, chains = 4, iter = 2500,
                            save_pars = save_pars(all = TRUE),
                            control = list(adapt_delta = 0.9999, max_treedepth = 14))
#cod_time.series_zinb  <- add_criterion(cod_time.series_zinb, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(mhw_zinb, file = "output/cod_heatwave_zinb_cook_anton.rds")

mhw_zinb <- readRDS("./output/cod_heatwave_zinb_cook_anton.rds")
check_hmc_diagnostics(mhw_zinb$fit)
neff_lowest(mhw_zinb$fit)
rhat_highest(mhw_zinb$fit)
summary(mhw_zinb)
bayes_R2(mhw_zinb)


y <- dat_temp$cod
yrep_mhw_zinb  <- fitted(mhw_zinb, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_mhw_zinb[sample(nrow(yrep_mhw_zinb), 25), ]) +
  xlim(0, 500) +
  ggtitle("mhw_zinb")


trace_plot(mhw_zinb$fit)

## Cod predicted effects ---------------------------------------

## 95% CI
ce1s_1 <- conditional_effects(mhw_zinb, effect = "heatwave_fac", re_formula = NA,
                              probs = c(0.025, 0.975))  

plot.cod <- ce1s_1$heatwave_fac %>%
  select(heatwave_fac, estimate__, se__, lower__, upper__)

# reorder heatwave status for plotting
plot.cod <- plot.cod %>%
  mutate(order = case_when(
    heatwave_fac == "before" ~ 1,
    heatwave_fac == "during" ~ 2,
    heatwave_fac == "adjacent" ~ 3,
    heatwave_fac == "after" ~ 4),
    heatwave_fac = reorder(heatwave_fac, order))

ggplot(plot.cod, aes(heatwave_fac, estimate__)) +
  geom_col(color = "black", fill = "grey") +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Age-0 cod / set") +
  scale_y_continuous(breaks=c(1,5,10,50,100,200,300), minor_breaks = NULL) +
  coord_trans(y = "pseudo_log") +
  xlab("Heatwave status")
## does not reduce uncertainty of estimates

### model heatwave effects for entire set of bays to see if we get more precision in estimates!
## Define model formula
time.series_formula <-  bf(cod ~ heatwave_fac + s(julian, k = 4) + (1 | bay_fac/site_fac) + (1 year_fac),
                           zi ~ heatwave_fac + s(julian, k = 4) + (1 | bay_fac/site_fac)+ (1 year_fac))
## cod fit: zero-inflated --------------------------------------
cod_time.series_zinb <- brm(time.series_formula,
                            data = dat,
                            prior = priors_zinb,
                            family = zinb,
                            cores = 4, chains = 4, iter = 4000,
                            save_pars = save_pars(all = TRUE),
                            control = list(adapt_delta = 0.999, max_treedepth = 11))
#cod_time.series_zinb  <- add_criterion(cod_time.series_zinb, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_time.series_zinb, file = "output/cod_heatwave_zinb_all_bays.rds")

cod_time.series_zinb <- readRDS("./output/cod_heatwave_zinb_all_bays.rds")
check_hmc_diagnostics(cod_time.series_zinb$fit)
neff_lowest(cod_time.series_zinb$fit)
rhat_highest(cod_time.series_zinb$fit)
summary(cod_time.series_zinb)
bayes_R2(cod_time.series_zinb)


y <- dat$cod
yrep_cod_time.series_zinb  <- fitted(cod_time.series_zinb, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_time.series_zinb[sample(nrow(yrep_cod_time.series_zinb), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_time.series_zinb")

trace_plot(cod_time.series_zinb$fit)

## Cod predicted effects ---------------------------------------

## 95% CI
ce1s_1 <- conditional_effects(cod_time.series_zinb, effect = "heatwave_fac", re_formula = NA,
                              probs = c(0.025, 0.975))  

plot.cod <- ce1s_1$heatwave_fac %>%
  select(heatwave_fac, estimate__, se__, lower__, upper__)

# reorder heatwave status for plotting
plot.cod <- plot.cod %>%
  mutate(order = case_when(
    heatwave_fac == "before" ~ 1,
    heatwave_fac == "during" ~ 2,
    heatwave_fac == "adjacent" ~ 3,
    heatwave_fac == "after" ~ 4),
    heatwave_fac = reorder(heatwave_fac, order))

ggplot(plot.cod, aes(heatwave_fac, estimate__)) +
  geom_col(color = "black", fill = "grey") +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Age-0 cod / set") +
  scale_y_continuous(breaks=c(1,5,10,50,100,200,300), minor_breaks = NULL) +
  coord_trans(y = "pseudo_log") +
  xlab("Heatwave status")

ggsave("./figs/seine_cod_age0_abundance_heatwave_all_bays.png", width = 4, height = 6, units = 'in')

## cod brms: setup ---------------------------------------------

## And now regional model for 2018-2023 - including year as a group-level effect

## Limit data to desired years
dat_temp <- dat_reduced %>%
  filter(year %in% 2018:2023)


## Define model formula
time.series_formula <-  bf(cod ~ region_fac + s(julian, k = 4) + (1 | bay_fac/site_fac) + (1 | year_fac),
                           zi ~ region_fac + s(julian, k = 4) + (1 | bay_fac/site_fac) + (1 | year_fac))

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


## cod fit: zero-inflated --------------------------------------
region_zinb <- brm(time.series_formula,
                              data = dat_temp,
                              prior = priors_zinb,
                              family = zinb,
                              cores = 4, chains = 4, iter = 2500,
                              save_pars = save_pars(all = TRUE),
                              control = list(adapt_delta = 0.9999, max_treedepth = 11))

#cod_time.series_zinb  <- add_criterion(cod_time.series_zinb, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(region_zinb, file = "output/region_zinb.rds")

region_zinb <- readRDS("./output/region_zinb.rds")
check_hmc_diagnostics(region_zinb$fit)
neff_lowest(region_zinb$fit)
rhat_highest(region_zinb$fit)
summary(region_zinb)
bayes_R2(region_zinb)


y <- dat_temp$cod
yrep_region_zinb  <- fitted(region_zinb, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = region_zinb[sample(nrow(region_zinb), 25), ]) +
  xlim(0, 500) +
  ggtitle("region_zinb")

trace_plot(region_zinb$fit)


## Cod predicted effects ---------------------------------------
## region estimates
## 95% CI
ce1s_1 <- conditional_effects(region_zinb, effect = "region_fac", re_formula = NA,
                              probs = c(0.025, 0.975))  

plot.cod <- ce1s_1$region_fac %>%
  select(region_fac, estimate__, se__, lower__, upper__)

# reorder heatwave status for plotting
plot.cod <- plot.cod %>%
  mutate(order = case_when(
    region_fac == "Shumagin Islands" ~ 1,
    region_fac == "Agripina" ~ 2,
    region_fac == "Kaiugnak" ~ 3,
    region_fac == "Eastern Kodiak" ~ 4),
    region_fac = reorder(region_fac, order))

ggplot(plot.cod, aes(region_fac, estimate__)) +
  geom_col(color = "black", fill = "grey") +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Age-0 cod / set") +
  scale_y_continuous(breaks=c(1,5,10,50,100,200,300), minor_breaks = NULL) +
  coord_trans(y = "pseudo_log") +
  xlab("Heatwave status")

ggsave("./figs/seine_cod_age0_abundance_region_model.png", width = 4, height = 6, units = 'in')


# round, rename columns, and save
plot.cod <- plot.cod[,1:5]
plot.cod[,2:5] <- round(plot.cod[,2:5], 2)
names(plot.cod) <- c("region", "cod_per_set", "cod_se", "cod_95percent_LCI", "cod_95percent_UCI")
write.csv(plot.cod, "./output/seine_cod_age0_abundance_region.csv", row.names = F)
