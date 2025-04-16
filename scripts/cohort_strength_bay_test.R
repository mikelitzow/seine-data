# fit brms models to update time series
library(plyr)
library(dplyr)
library(tidyverse)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
source("./scripts/stan_utils.R")

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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

## prepare data with only cook and ALB --------------------------------------------
#15 bays for full model, filter so only Cook and ALB
dat2 <- filter(dat, bay == "Cook" | bay == "Anton Larsen" )
distinct(dat2,bay_fac)
#this is only 2 bays

## prepare cod data --------------------------------------------
dat2$cod2 <- dat2$cod.age.0
dat2$bay_fac <- as.factor(dat2$bay)
dat2$year_fac <- as.factor(dat2$year)
dat2$site_fac <- as.factor(dat2$site)

## cod brms: setup ---------------------------------------------

## Define model formula - simpler model without site factor
time.series_formula <-  bf(cod2 ~ year_fac + s(julian, k = 4) + (1 | bay_fac),
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

cod_time.series_zinb_test <- brm(time.series_formula,
                            data = dat2,
                            prior = priors_zinb,
                            family = zinb,
                            cores = 4, chains = 4, iter = 3000,
                            save_pars = save_pars(all = TRUE),
                            control = list(adapt_delta = 0.9999, max_treedepth = 12))

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
  scale_y_continuous(breaks=c(1,5,10,50,100,200,500,1000,1500), minor_breaks = NULL) +
  coord_trans(y = "pseudo_log") +
  theme(axis.title.x = element_blank())

ggsave("./figs/seine_cod_age0_abundance_estimates_Anton_Cook.png", width = 6, height = 4, units = 'in')

# round, rename columns, and save
plot.cod_test[,2:5] <- round(plot.cod_test[,2:5], 2)
names(plot.cod_test) <- c("year", "cod_per_set", "cod_se", "cod_95percent_LCI", "cod_95percent_UCI")
write.csv(plot.cod_test, "./output/seine_cod_age0_abundance_estimates_Cook_Anton.csv", row.names = F)

# second test - full set of Bays with no site_fac term---------------------

## Define model formula - simpler model without site factor
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

cod_time.series_zinb_test2 <- brm(time.series_formula,
                                 data = dat,
                                 prior = priors_zinb,
                                 family = zinb,
                                 cores = 4, chains = 4, iter = 3000,
                                 save_pars = save_pars(all = TRUE),
                                 control = list(adapt_delta = 0.999, max_treedepth = 12))

saveRDS(cod_time.series_zinb_test2, file = "output/cod_time.series_zinb_test2.rds")

cod_time.series_zinb_test2 <- readRDS("./output/cod_time.series_zinb_test2.rds")
check_hmc_diagnostics(cod_time.series_zinb_test2$fit)
neff_lowest(cod_time.series_zinb_test2$fit)
rhat_highest(cod_time.series_zinb_test2$fit)
summary(cod_time.series_zinb_test2)
bayes_R2(cod_time.series_zinb_test2)
plot(cod_time.series_zinb_test2$criteria$loo, "k")
plot(conditional_smooths(cod_time.series_zinb_test2), ask = FALSE)
y <- dat$cod
yrep_cod_time.series_zinb_test2  <- fitted(cod_time.series_zinb_test2, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_time.series_zinb_test2[sample(nrow(yrep_cod_time.series_zinb_test2), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_time.series_zinb_test2")
pdf("./figs/trace_cod_time.series_zinb_test2.pdf", width = 6, height = 4)
trace_plot(cod_time.series_zinb_test2$fit)
dev.off()

## Cod predicted effects ---------------------------------------

## 95% CI
ce1s_1 <- conditional_effects(cod_time.series_zinb_test2, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))  

plot.cod_test2 <- ce1s_1$year_fac %>%
  select(year_fac, estimate__, se__, lower__, upper__)


ggplot(plot.cod_test2, aes(year_fac, estimate__)) +
  geom_col(color = "black", fill = "grey") +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Age-0 cod / set") +
  scale_y_continuous(breaks=c(1,5,10,50,100,200,300), minor_breaks = NULL) +
  coord_trans(y = "pseudo_log") +
  theme(axis.title.x = element_blank())

ggsave("./figs/seine_cod_age0_abundance_estimates_all_bays.png", width = 6, height = 4, units = 'in')

# round, rename columns, and save
plot.cod_test2[,2:5] <- round(plot.cod_test2[,2:5], 2)
names(plot.cod_test2) <- c("year", "cod_per_set", "cod_se", "cod_95percent_LCI", "cod_95percent_UCI")
write.csv(plot.cod_test2, "./output/seine_cod_age0_abundance_estimates_all_bays.csv", row.names = F)


## compare the two models ------
test2 <- read.csv("./output/seine_cod_age0_abundance_estimates_all_bays.csv") %>%
  rename(all = cod_per_set,
         all_LCI = cod_95percent_LCI,
         all_UCI = cod_95percent_UCI,
         all_se = cod_se) 

test1 <- read.csv("./output/seine_cod_age0_abundance_estimates_Cook_Anton.csv") %>%
  rename(ant_cook = cod_per_set,
         ant_cook_LCI = cod_95percent_LCI,
         ant_cook_UCI = cod_95percent_UCI,
         ant_cook_se = cod_se) 

plot_both <- left_join(test1, test2)


ggplot(plot_both, aes(ant_cook, all)) +
  geom_point() +
  ylab("All bays cod / set") +
  xlab("Anton's Cook cod / set") +
  geom_abline(slope = 1) +
  geom_smooth(method = "gam", se = F)

ggplot(plot_both, aes(ant_cook_se, all_se)) +
  geom_point() +
  ylab("All bays SE") +
  xlab("Anton's Cook SE") +
  geom_abline(slope = 1) +
  geom_smooth(method = "gam", se = F)


ggplot(filter(plot_both, year %in% 2018:2024), aes(ant_cook, all)) +
  geom_point() +
  ylab("All bays cod / set") +
  xlab("Anton's Cook cod / set") +
  geom_errorbar(aes(ymin = all_LCI, ymax = all_UCI)) +
  geom_errorbarh(aes(xmin = ant_cook_LCI, xmax = ant_cook_UCI))

ggplot(filter(plot_both, year %in% 2018:2024), aes(all, ant_cook)) +
  geom_point() +
  xlab("All bays cod / set") +
  ylab("Anton's Cook cod / set") +
  geom_errorbarh(aes(xmin = all_LCI, xmax = all_UCI)) +
  geom_errorbar(aes(ymin = ant_cook_LCI, ymax = ant_cook_UCI))



ggplot(filter(plot_both, year %in% 2018:2024), aes(all, ant_cook)) +
  geom_point() +
  xlab("All bays cod / set") +
  ylab("Anton's Cook cod / set") +
  geom_errorbarh(aes(xmin = all_LCI, xmax = all_UCI)) +
  geom_errorbar(aes(ymin = ant_cook_LCI, ymax = ant_cook_UCI))

# side by side
test2 <- read.csv("./output/seine_cod_age0_abundance_estimates_all_bays.csv") %>%
  mutate(group = "Full_survey") 

test1 <- read.csv("./output/seine_cod_age0_abundance_estimates_Cook_Anton.csv") %>%
  mutate(group = "Reduced_survey") 

plot_stack <- rbind(test1, test2) %>%
  rename(CPUE = cod_per_set,
         Survey = group,
         LCI = cod_95percent_LCI,
         UCI = cod_95percent_UCI)

ggplot(filter(plot_stack, year %in% 2018:2024), aes(year, CPUE, fill = Survey)) +
  geom_col(position = position_dodge(width = 0.9)) +
  scale_x_continuous(breaks = 2018:2024) +
  scale_fill_manual(values = cb[c(2,6)])

ggsave("./figs/full_reduced_survey_point_estimates.png", width = 6, height = 3.5, units = "in")

ggplot(filter(plot_stack, year %in% 2018:2024), aes(year, CPUE, fill = Survey)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = LCI, ymax = UCI), position = position_dodge(width = 0.9), width = 0.3) +
  scale_x_continuous(breaks = 2018:2024) +
  scale_fill_manual(values = cb[c(2,6)])

ggsave("./figs/full_reduced_survey_point_estimates_.png", width = 6, height = 3.5, units = "in")

test1 <- read.csv("./output/seine_cod_age0_abundance_estimates_Cook_Anton.csv") %>%
  mutate(group = "All")

test2 <- read.csv("./output/seine_cod_age0_abundance_estimates_all_bays.csv") %>%
  mutate(group = "Anton_Cook")

plot_both <- rbind(test1, test2)

ggplot(plot_both, aes(ant_cook, all)) +
  geom_point() +
  ylab("All bays cod / set") +
  xlab("Anton's Cook cod / set") +
  geom_abline(slope = 1) +
  geom_smooth(method = "gam", se = F)



# third test - reduced set of Bays with no site_fac term---------------------

dat3 <- dat %>%
  filter(!bay_fac %in% c("Agripina", "Fox", "Mitrofania", "Port Wrangell", "Rodman Reach", "Ugak"))

## Define model formula - simpler model without site factor
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

cod_time.series_zinb_test3 <- brm(time.series_formula,
                                  data = dat3,
                                  prior = priors_zinb,
                                  family = zinb,
                                  cores = 4, chains = 4, iter = 3000,
                                  save_pars = save_pars(all = TRUE),
                                  control = list(adapt_delta = 0.999, max_treedepth = 12))

saveRDS(cod_time.series_zinb_test3, file = "output/cod_time.series_zinb_test3.rds")

cod_time.series_zinb_test3 <- readRDS("./output/cod_time.series_zinb_test3.rds")
check_hmc_diagnostics(cod_time.series_zinb_test3$fit)
neff_lowest(cod_time.series_zinb_test3$fit)
rhat_highest(cod_time.series_zinb_test3$fit)
summary(cod_time.series_zinb_test3)
bayes_R2(cod_time.series_zinb_test3)
plot(cod_time.series_zinb_test3$criteria$loo, "k")
plot(conditional_smooths(cod_time.series_zinb_test3), ask = FALSE)
y <- dat$cod
yrep_cod_time.series_zinb_test3  <- fitted(cod_time.series_zinb_test3, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod_time.series_zinb_test3[sample(nrow(yrep_cod_time.series_zinb_test3), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod_time.series_zinb_test3")
pdf("./figs/trace_cod_time.series_zinb_test3.pdf", width = 6, height = 4)
trace_plot(cod_time.series_zinb_test3$fit)
dev.off()

## Cod predicted effects ---------------------------------------

## 95% CI
ce1s_1 <- conditional_effects(cod_time.series_zinb_test3, effect = "year_fac", re_formula = NA,
                              probs = c(0.025, 0.975))  

plot.cod_test3 <- ce1s_1$year_fac %>%
  select(year_fac, estimate__, se__, lower__, upper__)


ggplot(plot.cod_test3, aes(year_fac, estimate__)) +
  geom_col(color = "black", fill = "grey") +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  ylab("Age-0 cod / set") +
  scale_y_continuous(breaks=c(1,5,10,50,100,200,300,500), minor_breaks = NULL) +
  coord_trans(y = "pseudo_log") +
  theme(axis.title.x = element_blank())

ggsave("./figs/seine_cod_age0_abundance_estimates_resstricted_bays.png", width = 6, height = 4, units = 'in')

# round, rename columns, and save
plot.cod_test3[,2:5] <- round(plot.cod_test3[,2:5], 2)
names(plot.cod_test3) <- c("year", "cod_per_set", "cod_se", "cod_95percent_LCI", "cod_95percent_UCI")
write.csv(plot.cod_test3, "./output/seine_cod_age0_abundance_estimates_restricted_bays.csv", row.names = F)

## compare the three models ------
test2 <- read.csv("./output/seine_cod_age0_abundance_estimates_all_bays.csv") %>%
  rename(all = cod_per_set,
         all_LCI = cod_95percent_LCI,
         all_UCI = cod_95percent_UCI,
         all_se = cod_se) 

test1 <- read.csv("./output/seine_cod_age0_abundance_estimates_Cook_Anton.csv") %>%
  rename(ant_cook = cod_per_set,
         ant_cook_LCI = cod_95percent_LCI,
         ant_cook_UCI = cod_95percent_UCI,
         ant_cook_se = cod_se) 

test3 <- read.csv("./output/seine_cod_age0_abundance_estimates_restricted_bays.csv") %>%
  rename(restricted = cod_per_set,
         restr_LCI = cod_95percent_LCI,
         restr_UCI = cod_95percent_UCI,
         restr_se = cod_se) 



plot_all <- left_join(test1, test2) %>%
  left_join(.,test3)


plot_ts <- plot_all %>%
  select(year, ant_cook, all, restricted) %>%
  pivot_longer(cols = -year)

ggplot(plot_ts, aes(year, value, color = name)) +
  geom_point() +
  geom_line() +
  ylab("Age-0 cod / set") +
  # scale_y_continuous(breaks=c(1,5,10,50,100,200,300,500), minor_breaks = NULL) +
  # coord_trans(y = "pseudo_log") +
  theme(axis.title.x = element_blank()) +
  scale_color_manual(values = cb[c(2,4,6)])

plot_se <- plot_all %>%
  select(year, ant_cook_se, all_se, restr_se) %>%
  pivot_longer(cols = -year)

ggplot(plot_se, aes(year, value, color = name)) +
  geom_point() +
  geom_line() +
  ylab("Standard error") +
  # scale_y_continuous(breaks=c(1,5,10,50,100,200,300,500), minor_breaks = NULL) +
  # coord_trans(y = "pseudo_log") +
  theme(axis.title.x = element_blank()) +
  scale_color_manual(values = cb[c(2,4,6)])

plot_CI <- plot_all %>%
  mutate(ant_cook_CI = ant_cook_UCI - ant_cook_LCI,
         all_CI = all_UCI - all_LCI,
         restr_CI = restr_UCI - restr_LCI) %>%
  select(year, ant_cook_CI, all_CI, restr_CI) %>%
  pivot_longer(cols = -year)

ggplot(plot_CI, aes(year, value, color = name)) +
  geom_point() +
  geom_line() +
  ylab("Credible interval range") +
  # scale_y_continuous(breaks=c(1,5,10,50,100,200,300,500), minor_breaks = NULL) +
  # coord_trans(y = "pseudo_log") +
  theme(axis.title.x = element_blank()) +
  scale_color_manual(values = cb[c(2,4,6)])

