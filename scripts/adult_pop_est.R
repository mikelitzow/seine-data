library(tidyverse)
library(brms)
library(rstan)

source("./scripts/stan_utils.R")

cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

theme_set(theme_bw())
###what follows is population estimate by age of Pacific cod from bottom trawl survey###
##data from Pete Hulson. I combined number for all ages and I am only using ages 1-5

getwd()

dat <- read.csv("./data/twl_cod_pop_est.csv")

str(dat)

dat$age <-as.character(dat$age)

ggplot(dat, aes(year, tot_num)) +
  geom_point()+
  facet_wrap(~age)+
    ylab("Total number all sexes combined by age") 

ggplot(dat, aes(age, tot_num)) +
      geom_point()+
      facet_wrap(~year) +
    ylab("Total number all sexes combined") 
##but our seine data begin in 2006, so omit early years

dat_restricted <- filter(dat, year >= 2005)

ggplot(dat_restricted, aes(age, tot_num)) +
  geom_point()+
  facet_wrap(~year) +
  ylab("Total number all sexes combined") 
#I don't think this plot tells us anything about assumed mortality in late juveniles
#it is just the number by age in each year. I want to follow year class to look at mortality

dat2 <- read.csv("./data/twl_cod_pop_est_cohort.csv")
str(dat2)
dat2_restr <- filter(dat2, trawl_year >= 2007)

ggplot(dat2_restr, aes(age, total_num)) +
  geom_point()+
  facet_wrap(~trawl_year)+
  title(main = "Trawl estimates by age for cohort year", xlab = "Age", ylab = "Total number all sexes combined")
##I think this is a mess. not sure what it is telling me. THe year is the yearall messed up.
q()
