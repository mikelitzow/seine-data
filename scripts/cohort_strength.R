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

# summarize
summ <- dat %>%
  group_by(year) %>%
  summarise(n_bays = length(unique(bay_fac)),
            n_sets = n())

summ

## cod brms: setup ---------------------------------------------

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


## cod fit: zero-inflated --------------------------------------
cod_time.series_zinb <- brm(time.series_formula,
                            data = dat,
                            prior = priors_zinb,
                            family = zinb,
                            cores = 4, chains = 4, iter = 4000,
                            save_pars = save_pars(all = TRUE),
                            control = list(adapt_delta = 0.999, max_treedepth = 11))
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
y <- dat$cod
yrep_cod_time.series_zinb  <- fitted(cod_time.series_zinb, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_time.series_zinb[sample(nrow(yrep_cod_time.series_zinb), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_time.series_zinb")
pdf("./figs/trace_cod_time.series_zinb.pdf", width = 6, height = 4)
trace_plot(cod_time.series_zinb$fit)
dev.off()

## Cod predicted effects ---------------------------------------

## 95% CI
ce1s_1 <- conditional_effects(cod_time.series_zinb, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))  

plot.cod <- ce1s_1$year_fac %>%
  select(year_fac, estimate__, se__, lower__, upper__)


ggplot(plot.cod, aes(year_fac, estimate__)) +
  geom_col(color = "black", fill = "grey") +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Age-0 cod / set") +
  scale_y_continuous(breaks=c(1,5,10,50,100,200,300), minor_breaks = NULL) +
  coord_trans(y = "pseudo_log") +
  theme(axis.title.x = element_blank())

ggsave("./figs/seine_cod_age0_abundance_estimates.png", width = 6, height = 4, units = 'in')

# round, rename columns, and save
plot.cod[,2:5] <- round(plot.cod[,2:5], 2)
names(plot.cod) <- c("year", "cod_per_set", "cod_se", "cod_95percent_LCI", "cod_95percent_UCI")
write.csv(plot.cod, "./output/seine_cod_age0_abundance_estimates.csv", row.names = F)

## prepare pollock data --------------------------------------------
dat$pollock <- dat$pollock.age.0

# restrict to long-term sites and AK Peninsula bays with high proportion of positive catches
levels(dat$bay_fac)
keep <- c("Agripina", "Anton Larsen", "Balboa", "Cook", "Mitrofania", "Port Wrangell") 
dat <- dat %>%
  filter(bay_fac %in% keep)

# check we have all the bays!
unique(dat$bay)

# summarize
summ <- dat %>%
  group_by(year) %>%
  dplyr::summarise(n_sets = n(), n_bays = length(unique(bay_fac))
            )

summ


## pollock brms: setup ---------------------------------------------

## Define model formula
time.series_formula <-  bf(pollock ~ year_fac + s(julian, k = 4) + (1 | bay_fac/site_fac),
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


## pollock fit: zero-inflated --------------------------------------
pollock_time.series_zinb <- brm(time.series_formula,
                            data = dat,
                            prior = priors_zinb,
                            family = zinb,
                            cores = 4, chains = 4, iter = 4000,
                            save_pars = save_pars(all = TRUE),
                            control = list(adapt_delta = 0.999, max_treedepth = 11))
#pollock_time.series_zinb  <- add_criterion(pollock_time.series_zinb, "bayes_R2")
saveRDS(pollock_time.series_zinb, file = "output/pollock_time.series_zinb.rds")

pollock_time.series_zinb <- readRDS("./output/pollock_time.series_zinb.rds")
check_hmc_diagnostics(pollock_time.series_zinb$fit)
neff_lowest(pollock_time.series_zinb$fit)
rhat_highest(pollock_time.series_zinb$fit)
summary(pollock_time.series_zinb)
bayes_R2(pollock_time.series_zinb)

plot(conditional_smooths(pollock_time.series_zinb), ask = FALSE)
y <- dat$pollock
yrep_pollock_time.series_zinb  <- fitted(pollock_time.series_zinb, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_pollock_time.series_zinb[sample(nrow(yrep_pollock_time.series_zinb), 25), ]) +
  xlim(0, 500) +
  ggtitle("pollock_time.series_zinb")
pdf("./figs/trace_pollock_time.series_zinb.pdf", width = 6, height = 4)
trace_plot(pollock_time.series_zinb$fit)
dev.off()

## pollock predicted effects ---------------------------------------

## 95% CI
ce1s_1 <- conditional_effects(pollock_time.series_zinb, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))  

plot.pollock <- ce1s_1$year_fac %>%
  select(year_fac, estimate__, se__, lower__, upper__)


ggplot(plot.pollock, aes(year_fac, estimate__)) +
  geom_col(color = "black", fill = "grey") +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Age-0 pollock / set") +
  scale_y_continuous(breaks=c(0,1,5,10,50,100,200,300), minor_breaks = NULL) +
  coord_trans(y = "pseudo_log") +
  theme(axis.title.x = element_blank())

ggsave("./figs/seine_pollock_age0_abundance_estimates.png", width = 6, height = 4, units = 'in')

# round, rename columns, and save
plot.pollock[,2:5] <- round(plot.pollock[,2:5], 2)
names(plot.pollock) <- c("year", "pollock_per_set", "pollock_se", "pollock_95percent_LCI", "pollock_95percent_UCI")
write.csv(plot.pollock, "./output/seine_pollock_age0_abundance_estimates.csv", row.names = F)

# get time series mean and quantiles
summ <- plot.pollock %>%
  mutate(log_cpue = log(pollock_per_set)) %>%
  select(year, log_cpue)

ggplot(summ, aes(log_cpue)) +
  geom_histogram(bins = 6, fill = "grey", color = "black")
  
mean(summ$log_cpue)
