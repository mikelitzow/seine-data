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


## cod brms: setup ---------------------------------------------

## Define model formula
time.series_formula <-  bf(cod ~ heatwave_fac + s(julian, k = 4) + (1 | bay_fac/site_fac) + (1 | year_fac),
                           zi ~ heatwave_fac + s(julian, k = 4) + (1 | bay_fac/site_fac)+ (1 | year_fac))

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
cod_time.series_zinb <- brm(time.series_formula,
                            data = dat_reduced,
                            prior = priors_zinb,
                            family = zinb,
                            cores = 4, chains = 4, iter = 4000,
                            save_pars = save_pars(all = TRUE),
                            control = list(adapt_delta = 0.999, max_treedepth = 11))
#cod_time.series_zinb  <- add_criterion(cod_time.series_zinb, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(cod_time.series_zinb, file = "output/cod_heatwave_zinb.rds")

cod_time.series_zinb <- readRDS("./output/cod_heatwave_zinb.rds")
check_hmc_diagnostics(cod_time.series_zinb$fit)
neff_lowest(cod_time.series_zinb$fit)
rhat_highest(cod_time.series_zinb$fit)
summary(cod_time.series_zinb)
bayes_R2(cod_time.series_zinb)


y <- dat_reduced$cod
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

ggsave("./figs/seine_cod_age0_abundance_heatwave_restricted_bays.png", width = 4, height = 6, units = 'in')

# round, rename columns, and save
plot.cod <- plot.cod[,1:5]
plot.cod[,2:5] <- round(plot.cod[,2:5], 2)
names(plot.cod) <- c("heatwave status", "cod_per_set", "cod_se", "cod_95percent_LCI", "cod_95percent_UCI")
write.csv(plot.cod, "./output/seine_cod_age0_abundance_heatwave_restricted_bays.csv", row.names = F)

### model heatwave effects for entire set of bays to see if we get more precision in estimates!

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

cod_time.series_zinb <- readRDS("./output/cod_heatwave_zinb.rds")
check_hmc_diagnostics(cod_time.series_zinb$fit)
neff_lowest(cod_time.series_zinb$fit)
rhat_highest(cod_time.series_zinb$fit)
summary(cod_time.series_zinb)
bayes_R2(cod_time.series_zinb)


y <- dat_reduced$cod
yrep_cod_time.series_zinb  <- fitted(cod_time.series_zinb, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_time.series_zinb[sample(nrow(yrep_cod_time.series_zinb), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_time.series_zinb")

trace_plot(cod_time.series_zinb$fit)

