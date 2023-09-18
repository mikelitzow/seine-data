library(tidyverse)
library(brms)

theme_set(theme_bw())

seine <- read.csv("./output/seine_cod_age0_abundance_estimates.csv") %>%
  mutate(seine = as.vector(scale(log(cod_per_set)))) %>%
  select(year, seine)

model <- read.csv("./data/2022_Pcod_SAFE_recruitment.csv") %>%
  mutate(model = as.vector(scale(log(age0)))) %>%
  select(-age0)

data = left_join(seine, model)

ggplot(data[data$year <= 2019,], aes(seine, model)) +
  geom_text(aes(label = year))

## fit brms model ------------------------

formula <- bf(model ~ s(seine, k = 4))

seine_model_brm <- brm(formula,
                      data = data[data$year <= 2019,],
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
y <- data$model[data$year <= 2019]
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
