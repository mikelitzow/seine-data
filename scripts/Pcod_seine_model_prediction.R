library(tidyverse)
library(brms)
library(rstan)

source("./scripts/stan_utils.R")

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

theme_set(theme_bw())

seine <- read.csv("./output/seine_cod_age0_abundance_estimates.csv") %>%
  mutate(seine = as.vector(scale(log(cod_per_set)))) %>%
  select(year, seine)

model <- read.csv("./data/2023_Pcod_SAFE_recruitment.csv") %>%
  mutate(model = as.vector(scale(log(age0)))) %>%
  select(-age0)

names(model)[1] <- "year"

data = left_join(seine, model)

ggplot(data[data$year <= 2020,], aes(seine, model)) +
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

ggsave("./figs/seine_model_residuals.png", width = 6, height = 4, units = 'in')


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
trace_plot(seine_model_brm_3$fit)

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
