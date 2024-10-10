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

distinct(dat,bay_fac)
#15 bays for full model

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
#cod_time.series_zinb  <- add_criterion(cod_time.series_zinb, c("loo", "bayes_R2"), moment_match = TRUE)
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

################################ test if only cook and ALB
dat <- read.csv("./data/age.0_cod_pollock_seine_cpue.csv", row.names = 1)

## prepare cod data --------------------------------------------
dat$cod <- dat$cod.age.0
dat$bay_fac <- as.factor(dat$bay)
dat$year_fac <- as.factor(dat$year)
dat$site_fac <- as.factor(dat$site)

## prepare data with only cook and ALB --------------------------------------------
#15 bays for full model, filter so only Cook and ALB
dat <- filter(dat, bay == "Cook" | bay == "Anton Larsen" )
distinct(dat,bay_fac)
#this is only 2 bays

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
cod_time.series_zin_test <- brm(time.series_formula,
                            data = dat,
                            prior = priors_zinb,
                            family = zinb,
                            cores = 4, chains = 4, iter = 4000,
                            save_pars = save_pars(all = TRUE),
                            control = list(adapt_delta = 0.999, max_treedepth = 11))

saveRDS(cod_time.series_zinb_test, file = "output/cod_time.series_zinb_test.rds")

cod_time.series_zinb_test <- readRDS("./output/cod_time.series_zinb_test.rds")
check_hmc_diagnostics(cod_time.series_zinb_test$fit)
neff_lowest(cod_time.series_zinb_test$fit)
rhat_highest(cod_time.series_zinb_test$fit)
summary(cod_time.series_zinb_test)
bayes_R2(cod_time.series_zinb_test)
plot(cod_time.series_zinb_test$criteria$loo, "k")
plot(conditional_smooths(cod_time.series_zinb_test), ask = FALSE)
y <- dat$cod
yrep_cod_time.series_zinb_test  <- fitted(cod_time.series_zinb_test, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_time.series_zinb_test[sample(nrow(yrep_cod_time.series_zinb_test), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_time.series_zinb_test")
pdf("./figs/trace_cod_time.series_zinb_test.pdf", width = 6, height = 4)
trace_plot(cod_time.series_zinb_test$fit)
dev.off()

## Cod predicted effects ---------------------------------------

## 95% CI
ce1s_1 <- conditional_effects(cod_time.series_zinb_test, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))  

plot.cod_test <- ce1s_1$year_fac %>%
  select(year_fac, estimate__, se__, lower__, upper__)


ggplot(plot.cod_test, aes(year_fac, estimate__)) +
  geom_col(color = "black", fill = "grey") +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Age-0 cod / set") +
  scale_y_continuous(breaks=c(1,5,10,50,100,200,300), minor_breaks = NULL) +
  coord_trans(y = "pseudo_log") +
  theme(axis.title.x = element_blank())

ggsave("./figs/seine_cod_age0_abundance_estimates_test.png", width = 6, height = 4, units = 'in')

# round, rename columns, and save
plot.cod_test[,2:5] <- round(plot.cod_test[,2:5], 2)
names(plot.cod_test) <- c("year", "cod_per_set", "cod_se", "cod_95percent_LCI", "cod_95percent_UCI")
write.csv(plot.cod_test, "./output/seine_cod_age0_abundance_estimates_test.csv", row.names = F)
