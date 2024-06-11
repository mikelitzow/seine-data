library(tidyverse)
library(brms)
library(rstan)

source("./scripts/stan_utils.R")

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

theme_set(theme_bw())
###what follows is age-0 from seine and model age-3###

seine <- read.csv("./output/seine_cod_age0_abundance_estimates.csv") %>%
  mutate(seine = as.vector(scale(log(cod_per_set)))) %>%
  select(year, seine)

model <- read.csv("./data/2023_Pcod_SAFE_recruitment.csv") %>%
  mutate(model = as.vector(scale(log(age0)))) %>%
  select(-age0)

names(model)[1] <- "year"

data = left_join(seine, model)

ggplot(data, aes(seine, model)) +
  geom_text(aes(label = year))

## fit brms model ------------------------

formula <- bf(model ~ s(seine, k = 4))

seine_model_brm <- brm(formula,
                      data = data[data$year <= 2020,],
                      cores = 4, chains = 4, iter = 2000,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(seine_model_brm, file = "output/seine_model_brm.rds")

seine_model_brm <- readRDS("./output/seine_model_brm.rds")
check_hmc_diagnostics(seine_model_brm$fit)
neff_lowest(seine_model_brm$fit)
rhat_highest(seine_model_brm$fit)
summary(seine_model_brm)
bayes_R2(seine_model_brm)
y <- data$model[data$year <= 2020]
yrep_seine_model_brm  <- fitted(seine_model_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_seine_model_brm[sample(nrow(yrep_seine_model_brm), 25), ]) +
  xlim(-6, 6) +
  ggtitle("seine_model_brm")
trace_plot(seine_model_brm$fit)

## plot
ce1s_1 <- conditional_effects(seine_model_brm, effect = "seine", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(seine_model_brm, effect = "seine", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(seine_model_brm, effect = "seine", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$seine
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$seine[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$seine[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$seine[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$seine[["lower__"]]

model.plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Seine abundance", y = "Model recruitment") +
  geom_text(data = data[data$year <= 2019,], aes(seine, model, label = year), size=3) 

print(model.plot)

ggsave("./figs/seine_model_regression.png", width = 4, height = 3, units = 'in')

## plot residuals by year-------------------------------------------------------
seine_model_brm <- readRDS("./output/seine_model_brm.rds")

prediction_residuals <- as.data.frame(residuals(seine_model_brm)) %>%
  mutate(year = 2006:2020)


ggplot(prediction_residuals, aes(year, Estimate)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), lty = 0, fill = "grey", alpha = 0.5) +
  geom_hline(yintercept = 0, color = "dark grey") +
  geom_line(color = "red3") +
  geom_point(color = "red3") +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), lty = 0, fill = "grey", alpha = 0.5) +
  scale_x_continuous(breaks = 2006:2019,
                     minor_breaks = NULL) +
  labs(x = "Year class",
       y = "Posterior residuals")

ggsave("./figs/seine_model_residuals_2020.png", width = 6, height = 4, units = 'in')


## fit brms model to restricted set of years ------------------------

formula <- bf(model ~ s(seine, k = 4))

seine_model_brm_2 <- brm(formula,
                         data = data[data$year <= 2016,],
                         cores = 4, chains = 4, iter = 2000,
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(seine_model_brm_2, file = "output/seine_model_brm_2006_2016.rds")

seine_model_brm_2 <- readRDS("./output/seine_model_brm_2006_2016.rds")
check_hmc_diagnostics(seine_model_brm_2$fit)
neff_lowest(seine_model_brm_2$fit)
rhat_highest(seine_model_brm_2$fit)
summary(seine_model_brm_2)
bayes_R2(seine_model_brm_2)
y <- data$model[data$year <= 2019]
yrep_seine_model_brm_2  <- fitted(seine_model_brm_2, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_seine_model_brm_2[sample(nrow(yrep_seine_model_brm_2), 25), ]) +
  xlim(-6, 6) +
  ggtitle("seine_model_brm_2")
trace_plot(seine_model_brm_2$fit)

## plot
ce1s_1 <- conditional_effects(seine_model_brm_2, effect = "seine", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(seine_model_brm_2, effect = "seine", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(seine_model_brm_2, effect = "seine", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$seine
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$seine[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$seine[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$seine[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$seine[["lower__"]]

model.plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Seine abundance", y = "Model recruitment") +
  geom_text(data = data[data$year <= 2016,], aes(seine, model, label = year), size=3) 

print(model.plot)

ggsave("./figs/seine_model_regression_2006_2016.png", width = 4, height = 3, units = 'in')

##plot for ESASS talk 5/2024
model.plot.esass <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  scale_fill_manual(values = cb[c(2,6)]) + 
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "papayawhip") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "papayawhip") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "papayawhip") +
  geom_line(size = 1, color = "navajowhite1") +
  expand_limits(y = c(-4,2))+
  labs(x = "Seine estimate", y = "Assessment model estimate") +
  geom_text(data = data[data$year <= 2016,], aes(seine, model, label = year), size=4) 

print(model.plot.esass)
ggsave("./figs/seine_ESASSmodel_regression_2006_2016.png", width = 4, height = 3, units = 'in')
##done wtih ESASS figure


## fit brms model with era effect ------------------------

data <- data %>%
  mutate(era = if_else(year <= 2016, "2006-2016", "2017-2020"))

formula <- bf(model ~ seine*era)

seine_model_brm_3 <- brm(formula,
                         data = data[data$year <= 2020,],
                         cores = 4, chains = 4, iter = 2000,
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(seine_model_brm_3, file = "output/seine_model_brm_era.rds")

seine_model_brm_3 <- readRDS("./output/seine_model_brm_era.rds")
check_hmc_diagnostics(seine_model_brm_3$fit)
neff_lowest(seine_model_brm_3$fit)
rhat_highest(seine_model_brm_3$fit)
summary(seine_model_brm_3)
bayes_R2(seine_model_brm_3)
y <- data$model[data$year <= 2020]
yrep_seine_model_brm_3  <- fitted(seine_model_brm_3, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_seine_model_brm_3[sample(nrow(yrep_seine_model_brm_3), 25), ]) +
  xlim(-6, 6) +
  ggtitle("seine_model_brm_3")
trace_plot(seine_model_brm_3$fit)

## plot


seine_effect <- conditional_effects(seine_model_brm_3, effect = "seine:era", re_formula = NA,
                              probs = c(0.025, 0.975)) 


ggplot(seine_effect$`seine:era`, aes(effect1__, estimate__, color = era, fill = era)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, color = NA) +
  scale_color_manual(values = cb[c(2,6)]) +
  scale_fill_manual(values = cb[c(2,6)]) + 
  geom_text(data = data[data$year <= 2020,], aes( x = seine, y = model, label = year)) +
  labs(x = "Seine estimate",
       y = "Assessment model estimate")
  

ggsave("./figs/seine_model_regression_by_era.png", width = 6.25, height = 3.5, units = 'in')


## fit model with one seine slope, different era intercepts

formula <- bf(model ~ era + seine)

seine_model_brm_4 <- brm(formula,
                         data = data[data$year <= 2020,],
                         cores = 4, chains = 4, iter = 2000,
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(seine_model_brm_4, file = "output/seine_model_brm_era_intercept_one_slope.rds")

seine_model_brm_4 <- readRDS("./output/seine_model_brm_era_intercept_one_slope.rds")
check_hmc_diagnostics(seine_model_brm_4$fit)
neff_lowest(seine_model_brm_4$fit)
rhat_highest(seine_model_brm_4$fit)
summary(seine_model_brm_4)
bayes_R2(seine_model_brm_4)
y <- data$model[data$year <= 2020]
yrep_seine_model_brm_4  <- fitted(seine_model_brm_4, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_seine_model_brm_4[sample(nrow(yrep_seine_model_brm_4), 25), ]) +
  xlim(-6, 6) +
  ggtitle("seine_model_brm_4")
trace_plot(seine_age1_model_brm$fit)

compare <- loo(seine_model_brm_3, seine_model_brm_4, moment_match = T)

conditions <- make_conditions(seine_model_brm_4, vars = "era")

era_effect <- conditional_effects(seine_model_brm_4, effect = "seine", re_formula = NA, conditions = conditions,
                                    probs = c(0.025, 0.975)) 

data$era <- if_else(data$year <= 2016, "2006-2016", "2017-2020")

ggplot(era_effect$seine, aes(effect1__, estimate__, color = era, fill = era)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, color = NA) +
  scale_color_manual(values = cb[c(2,6)]) +
  scale_fill_manual(values = cb[c(2,6)]) + 
  geom_text(data = data[data$year <= 2020,], aes( x = seine, y = model, label = year)) +
  labs(x = "Seine estimate",
       y = "Assessment model estimate")


ggsave("./figs/seine_model_regression_era_intercepts_one_slope.png", width = 6.25, height = 3.5, units = 'in')


## just for kicks, fit through 2022
## fit model with one seine slope, different era intercepts
data$era <- if_else(data$year <= 2016, "2006-2016", "2017-2022")
formula <- bf(model ~ era + seine)

seine_model_brm_5 <- brm(formula,
                         data = data[data$year <= 2022,],
                         cores = 4, chains = 4, iter = 2000,
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(seine_model_brm_5, file = "output/seine_model_brm_era_intercept_one_slope_through_2022.rds")

seine_model_brm_5 <- readRDS("./output/seine_model_brm_era_intercept_one_slope_through_2022.rds")
check_hmc_diagnostics(seine_model_brm_5$fit)
neff_lowest(seine_model_brm_5$fit)
rhat_highest(seine_model_brm_5$fit)
summary(seine_model_brm_5)
bayes_R2(seine_model_brm_5)
y <- data$model[data$year <= 2020]
yrep_seine_model_brm_5  <- fitted(seine_model_brm_5, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_seine_model_brm_5[sample(nrow(yrep_seine_model_brm_5), 25), ]) +
  xlim(-6, 6) +
  ggtitle("seine_model_brm_5")
trace_plot(seine_model_brm_5$fit)

conditions <- make_conditions(seine_model_brm_5, vars = "era")

era_effect <- conditional_effects(seine_model_brm_5, effect = "seine", re_formula = NA, conditions = conditions,
                                  probs = c(0.025, 0.975)) 

data$era <- if_else(data$year <= 2016, "2006-2016", "2017-2022")

ggplot(era_effect$seine, aes(effect1__, estimate__, color = era, fill = era)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, color = NA) +
  scale_color_manual(values = cb[c(2,6)]) +
  scale_fill_manual(values = cb[c(2,6)]) + 
  geom_text(data = data[data$year <= 2022,], aes( x = seine, y = model, label = year)) +
  labs(x = "Seine estimate",
       y = "Assessment model estimate")


ggsave("./figs/seine_model_regression_era_intercepts_one_slope_through_2022.png", width = 6.25, height = 3.5, units = 'in')


##################
###what follows is age-0 from seine and age-1 from seine###
#need to lag the age-1 data so lagged it in excel and made .csv with Cohort column
#Cohort = cohort year
#year 2023 removed because do not have age-1 from cohort 2023 (not sampled yet)
##thus, this is for years 2006 - 2022 only

juv_seine <- read.csv("./output/seine_cod_age01_abundance_estimates.csv") %>%
  mutate(seine = as.vector(scale(log(cod0_per_set)))) %>%
  select(Cohort, seine)

names(juv_seine)[2] <- "age0"
# View(juv_seine)

juv_seine1 <- read.csv("./output/seine_cod_age01_abundance_estimates.csv") %>%
  mutate(seine1 = as.vector(scale(log(cod1_per_set +0.01)))) %>%
  select(Cohort, seine1)

names(juv_seine1)[2] <- "age1"
# View(juv_seine1)

juv_data = left_join(juv_seine, juv_seine1)

# View(juv_data) #looks good

ggplot(juv_data[juv_data$Cohort <= 2022,], aes(age0, age1)) +
  geom_text(aes(label = Cohort))

## fit brms model ------------------------

formula <- bf(age1 ~ s(age0, k = 4))

seine_model_brm <- brm(formula,
                       data = juv_data[juv_data$Cohort <= 2022,],
                       cores = 4, chains = 4, iter = 2000,
                       save_pars = save_pars(all = TRUE),
                       control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(seine_model_brm, file = "output/juv_seine_model_brm.rds")

seine_model_brm <- readRDS("./output/juv_seine_model_brm.rds")
check_hmc_diagnostics(seine_model_brm$fit)
neff_lowest(seine_model_brm$fit)
rhat_highest(seine_model_brm$fit)
summary(seine_model_brm)
bayes_R2(seine_model_brm)
y <- data$model[data$Cohort <= 2022]
yrep_seine_model_brm  <- fitted(seine_model_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_seine_model_brm[sample(nrow(yrep_seine_model_brm), 25), ]) +
  xlim(-6, 6) +
  ggtitle("seine_model_brm")
trace_plot(seine_model_brm$fit)

## plot
ce1s_1 <- conditional_effects(seine_model_brm, effect = "age0", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(seine_model_brm, effect = "age0", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(seine_model_brm, effect = "age0", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$age0
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$seine[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$seine[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$seine[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$seine[["lower__"]]

model.plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  #geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  #geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Seine abundance age-0", y = "Seine abundance age-1") +
  geom_text(data = juv_data[juv_data$Cohort <= 2022,], aes(age0, age1, label = Cohort), size=3) 

print(model.plot)
## seems like this plot driven by high 2020 and 2022 years. 
#To compare with age-0 to age-3 want to run again with years <=2020

ggsave("./figs/seine_0_1_model_regression2022.png", width = 4, height = 3, units = 'in')

#To compare with age-0 to age-3 era model, want to run age-0 to age-1 again with years <=2020
## fit brms model ------------------------

formula <- bf(age1 ~ s(age0, k = 4))

seine_model_brm <- brm(formula,
                       data = juv_data[juv_data$Cohort <= 2020,],
                       cores = 4, chains = 4, iter = 2000,
                       save_pars = save_pars(all = TRUE),
                       control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(seine_model_brm, file = "output/juv2020_seine_model_brm.rds")

seine_model_brm <- readRDS("./output/juv2020_seine_model_brm.rds")
check_hmc_diagnostics(seine_model_brm$fit)
neff_lowest(seine_model_brm$fit)
rhat_highest(seine_model_brm$fit)
summary(seine_model_brm)
bayes_R2(seine_model_brm)
y <- data$model[data$Cohort <= 2020]
yrep_seine_model_brm  <- fitted(seine_model_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_seine_model_brm[sample(nrow(yrep_seine_model_brm), 25), ]) +
  xlim(-6, 6) +
  ggtitle("seine_model_brm")
trace_plot(seine_model_brm$fit)

## plot
ce1s_1 <- conditional_effects(seine_model_brm, effect = "age0", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(seine_model_brm, effect = "age0", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(seine_model_brm, effect = "age0", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$age0
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$seine[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$seine[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$seine[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$seine[["lower__"]]

model.plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  #geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  #geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Seine abundance age-0", y = "Seine abundance age-1") +
  geom_text(data = juv_data[juv_data$Cohort <= 2020,], aes(age0, age1, label = Cohort), size=3) 

print(model.plot)
## interesting. now look at model output to age-1

ggsave("./figs/seine_0_1_model_regression2020.png", width = 4, height = 3, units = 'in')

###########################################
###plot age-1 vs model output
###what follows is age-1 from seine and model. 
###But I am confused and need help here: the model is age-0 and not age-3
# to compare age-1 to model, do we need the age-3 numbers from the model?
#need to ask Mike to help me lag model output

seine <- read.csv("./output/seine_cod_age1_abundance_estimates.csv") %>%
  mutate(seine = as.vector(scale(log(cod_per_set + 0.01)))) %>%
  select(year, seine)

seine

model <- read.csv("./data/2023_Pcod_SAFE_recruitment.csv") %>%
  mutate(model = as.vector(scale(log(age0)))) %>%
  select(-age0)

model

names(model)[1] <- "year"

# lag seine data
seine_1 <- seine %>%
  mutate(year = year-1)


data = left_join(seine_1, model)

ggplot(data[data$year <= 2022,], aes(seine, model)) +
  geom_text(aes(label = year))

## fit brms model ------------------------

data <- data %>%
  mutate(era = if_else(year <= 2016, "2005-2016", "2017-2022"))

formula <- bf(model ~ seine*era)

seine_age1_model_brm <- brm(formula,
                         data = data,
                         cores = 4, chains = 4, iter = 2000,
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(seine_age1_model_brm, file = "output/seine_age1_model_brm_era.rds")

seine_age1_model_brm <- readRDS("./output/seine_age1_model_brm_era.rds")
check_hmc_diagnostics(seine_age1_model_brm$fit)
neff_lowest(seine_age1_model_brm$fit)
rhat_highest(seine_age1_model_brm$fit)
summary(seine_age1_model_brm)
bayes_R2(seine_age1_model_brm)
y <- data$model[data$year <= 2022]
yrep_seine_age1_model_brm  <- fitted(seine_age1_model_brm, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_seine_age1_model_brm[sample(nrow(yrep_seine_age1_model_brm), 25), ]) +
  xlim(-6, 6) +
  ggtitle("seine_age1_model_brm")
trace_plot(seine_age1_model_brm$fit)

conditions <- make_conditions(seine_age1_model_brm, vars = "era")

era_effect <- conditional_effects(seine_age1_model_brm, effect = "seine", re_formula = NA, conditions = conditions,
                                  probs = c(0.025, 0.975)) 

data$era <- if_else(data$year <= 2016, "2005-2016", "2017-2022")

ggplot(era_effect$seine, aes(effect1__, estimate__, color = era, fill = era)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, color = NA) +
  scale_color_manual(values = cb[c(2,6)]) +
  scale_fill_manual(values = cb[c(2,6)]) + 
  geom_text(data = data, aes( x = seine, y = model, label = year)) +
  labs(x = "Seine estimate (from age-1 catch)",
       y = "Assessment model estimate")


ggsave("./figs/seine_age1_model_regression_era_.png", width = 6.25, height = 3.5, units = 'in')

# compare with a model fitting a single relationship between age-1 seine catches and assessment models

formula <- bf(model ~ seine)

seine_age1_model_brm_2 <- brm(formula,
                            data = data,
                            cores = 4, chains = 4, iter = 2000,
                            save_pars = save_pars(all = TRUE),
                            control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(seine_age1_model_brm_2, file = "output/seine_age1_model_brm_no_era.rds")

seine_age1_model_brm_2 <- readRDS("./output/seine_age1_model_brm_no_era.rds")
check_hmc_diagnostics(seine_age1_model_brm_2$fit)
neff_lowest(seine_age1_model_brm_2$fit)
rhat_highest(seine_age1_model_brm_2$fit)
summary(seine_age1_model_brm_2)
bayes_R2(seine_age1_model_brm_2)

# compare
compare <- loo(seine_age1_model_brm, seine_age1_model_brm_2, moment_match = T)

compare # era model is much better

#### fit era models to evaluate change point

output <- data.frame()

for(i in 2009:2017){ # fit with these years as the last year in era 1 
 # i <- 2013
data$era <- as.factor(if_else(data$year <= i, 1, 2))

mod <- lm(model ~ seine*era, data = data)

output <- rbind(output,
                data.frame(age = "age_1",
                           breakpoint = i + 0.5,
                           AICc = MuMIn::AICc(mod)))
  
}

ggplot(output, aes(breakpoint, AICc)) +
  geom_point() +
  geom_line()

# check best breakpoint (2013.5)
data$era <- as.factor(if_else(data$year <= 2013, 1, 2))

ggplot(data, aes(seine, model, color = era)) +
  geom_text(data = data, aes( x = seine, y = model, label = year))

# compare with stationary model AICc 
stationary_output <- data.frame()

mod <- lm(model ~ seine, data = data)

stationary_output <- rbind(stationary_output,
                           data.frame(age = "age_1",
                                      AICc = MuMIn::AICc(mod)))


# now add age 0 seine comparison with assessment model
seine <- read.csv("./output/seine_cod_age0_abundance_estimates.csv") %>%
  mutate(seine = as.vector(scale(log(cod_per_set)))) %>%
  select(year, seine)

model <- read.csv("./data/2023_Pcod_SAFE_recruitment.csv") %>%
  mutate(model = as.vector(scale(log(age0)))) %>%
  select(-age0)

names(model)[1] <- "year"

data = left_join(seine, model) %>%
  filter(year <= 2022)

for(i in 2009:2017){ # fit with these years as the last year in era 1 
  # i <- 2013
  data$era <- as.factor(if_else(data$year <= i, 1, 2))
  
  mod <- lm(model ~ seine*era, data = data)
  
  output <- rbind(output,
                  data.frame(age = "age_0",
                             breakpoint = i + 0.5,
                             AICc = MuMIn::AICc(mod)))
  
}

# and stationary model
mod <- lm(model ~ seine, data = data)

stationary_output <- rbind(stationary_output,
                           data.frame(age = "age_0",
                                      AICc = MuMIn::AICc(mod)))
##plot only age-0 for ESASS talk

output0 <- output %>%
  filter(age == "age_0")

ggplot(output0, aes(breakpoint, AICc, color = age)) +
  geom_point() +
  geom_line() +
  geom_hline(data = stationary_output, aes(yintercept = AICc), lty = 2) +
  scale_color_manual(values = cb[c(3,4)]) +
  labs(x = "Breakpoint")
ggsave("./figs/breakpoint_age0_AIC_plot.png", width = 6, height = 4, units = 'in')
##
ggplot(output, aes(breakpoint, AICc, color = age)) +
  geom_point() +
  geom_line() +
  geom_hline(data = stationary_output, aes(yintercept = AICc, color = age), lty = 2) +
  scale_color_manual(values = cb[c(3,4)]) +
  labs(x = "Breakpoint")

ggsave("./figs/breakpoint_AIC_plot.png", width = 6, height = 4, units = 'in')


## fit brms model to restricted set of years ------------------------
##because what I really want to know is if new control from 2016 - 2020
formula <- bf(model ~ s(seine, k = 4))

seine_model_brm_16 <- brm(formula,
                         data = data[data$year >= 2017,],
                         cores = 4, chains = 4, iter = 2000,
                         save_pars = save_pars(all = TRUE),
                         control = list(adapt_delta = 0.999, max_treedepth = 10))

saveRDS(seine_model_brm_16, file = "output/seine_model_brm_2016_2022.rds")

seine_model_brm_16 <- readRDS("./output/seine_model_brm_2016_2022.rds")
check_hmc_diagnostics(seine_model_brm_16$fit)
neff_lowest(seine_model_brm_16$fit)
rhat_highest(seine_model_brm_16$fit)
summary(seine_model_brm_16)
bayes_R2(seine_model_brm_16)
y <- data$model[data$year >= 2017]
yrep_seine_model_brm_16  <- fitted(seine_model_brm_16, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_seine_model_brm_16[sample(nrow(yrep_seine_model_brm_16), 25), ]) +
  xlim(-6, 6) +
  ggtitle("seine_model_brm_16")
trace_plot(seine_model_brm_16$fit)

## plot
ce1s_1 <- conditional_effects(seine_model_brm_16, effect = "seine", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(seine_model_brm_16, effect = "seine", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(seine_model_brm_16, effect = "seine", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$seine
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$seine[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$seine[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$seine[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$seine[["lower__"]]

model_16.plot <- ggplot(dat_ce) +
  aes(x = effect1__, y = estimate__) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "grey90") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "grey85") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "grey80") +
  geom_line(size = 1, color = "red3") +
  labs(x = "Seine abundance age-1", y = "Model recruitment") +
  geom_text(data = data[data$year >= 2017,], aes(seine, model, label = year), size=3) 

print(model_16.plot)

ggsave("./figs/seine_model_regression_2017_2022.png", width = 4, height = 3, units = 'in')
