# fit brms model to estimate age-1 abundance 
# analysis only includes bays with > 1 site catching age-1 cod in at least 1 year

library(dplyr)
library(plyr)
library(tidyverse)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
source("./scripts/stan_utils.R")
theme_set(theme_bw())

dat <- read.csv("./data/juv_cod_pollock_seine_cpue.csv", row.names = 1)

# check data
str(dat)

unique(dat$bay)
unique(dat$site)

## prepare data --------------------------------------------
dat$cod <- dat$cod.age.1
dat$bay_fac <- as.factor(dat$bay)
dat$year_fac <- as.factor(dat$year)
dat$site_fac <- as.factor(dat$site)

# summarize
summ <- dat %>%
  group_by(year) %>%
  dplyr::summarise(n_bays = length(unique(bay_fac)),
            n_sets = n())

summ

# check proportion of zeros by bays / sites
check <- dat %>%
  group_by(bay, site) %>%
  dplyr::summarise(prop_0 = 1-sum(cod>0)/n())

check
View(check) ##95 sites combos with age-1 present

# restrict to positive sites (at least one age-1 cod caught over time)
keep <- check %>%
  filter(prop_0 < 1)

nrow(keep) # 34 sites kept
View(keep)
# restrict dat to these sites 
dat <- dat %>%
  filter(site %in% keep$site)

# check proportion of zeros by bays / sites
check <- dat %>%
  group_by(bay, site) %>%
  dplyr::summarise(prop_0 = 1-sum(cod>0)/n())

check

# how many sites retained per bay?
check2 <- dat %>%
  group_by(bay) %>%
  dplyr::summarise(n_sites = length(unique(site)))

check2

# limit to 11 bays with > 1 site catching age 1 cod

dat <- dat %>%
  filter(bay %in% check2$bay[check2$n_sites > 1])

# plot to check
plot_check <- dat %>%
  group_by(bay, year) %>%
  dplyr::summarise(cpue = mean(cod.age.1))

ggplot(plot_check, aes(year, cpue)) +
  geom_col() +
  facet_wrap(~bay)

## cod brms: setup ---------------------------------------------

## Define model formula
time.series_formula <-  bf(cod ~ year_fac + s(julian, k = 4) + (1 | bay_fac),
                           zi ~ year_fac + s(julian, k = 4) + (1 | bay_fac))

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
saveRDS(cod_time.series_zinb, file = "output/age1_cod_time.series_zinb.rds")

cod_time.series_zinb <- readRDS("./output/age1_cod_time.series_zinb.rds")
check_hmc_diagnostics(cod_time.series_zinb$fit)
neff_lowest(cod_time.series_zinb$fit)
rhat_highest(cod_time.series_zinb$fit)
summary(cod_time.series_zinb)
bayes_R2(cod_time.series_zinb)
plot(cod_time.series_zinb$criteria$loo, "k")
plot(conditional_smooths(cod_time.series_zinb), ask = FALSE)
y <- dat$cod.age.1
yrep_cod_time.series_zinb  <- fitted(cod_time.series_zinb, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_time.series_zinb[sample(nrow(yrep_cod_time.series_zinb), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_time.series_zinb") +
  scale_x_continuous(trans = "pseudo_log")

trace_plot(cod_time.series_zinb$fit)


## Cod predicted effects ---------------------------------------

## 95% CI
ce1s_1 <- conditional_effects(cod_time.series_zinb, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))  

plot.cod <- ce1s_1$year_fac %>%
  select(year_fac, estimate__, se__, lower__, upper__)


ggplot(plot.cod, aes(year_fac, estimate__)) +
  geom_col(color = "black", fill = "grey") +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Age-1 cod / set") +
  scale_y_continuous(breaks=c(1,5,10,50,100,200,300), minor_breaks = NULL) +
  coord_trans(y = "pseudo_log") +
  theme(axis.title.x = element_blank())

ggsave("./figs/seine_cod_age1_abundance_estimates.png", width = 6, height = 4, units = 'in')

# round, rename columns, and save
plot.cod[,2:5] <- round(plot.cod[,2:5], 2)
names(plot.cod) <- c("year", "cod_per_set", "cod_se", "cod_95percent_LCI", "cod_95percent_UCI")
write.csv(plot.cod, "./output/seine_cod_age1_abundance_estimates.csv", row.names = F)
