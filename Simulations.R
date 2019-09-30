# LIBRARY IMPORTS ==================================================================================
library(brms)
library(bridgesampling)
library(tidyverse)
library(scales)
library(beepr)
library(future)
library(future.apply)
plan(multiprocess, workers = 6)

# GENERATE DATA ====================================================================================
# Define simulation-specific variables
n_labs <- 20
infants_by_lab <- rnorm(n_labs, 35, 20) %>% round() %>% pmax.int(16) # Mean and SD from MB1, min 16
n_infants <- sum(infants_by_lab)
age_min <- 165
age_max <- 320
pr_helper <- 0.64     # From Margoni & Surian (estimated true effect size)
pr_helper_sd <- 0.03  # From Margoni & Surian (confidence interval)

# Generate repeat datasets
n_sims <- 1 # Running multiple sims leads to brm errors with save_all_pars
data_sims <- tibble(lab_id = mapply(rep,
                                    paste0("lab",1:n_labs),
                                    infants_by_lab) %>%
                      unlist() %>%
                      as_factor(),
                    age_days = sample(age_min:age_max,
                                      n_infants,
                                      replace = T),
                    z_age_days = scale(age_days),
                    chose_helper = rbinom(n_infants, 1,
                                          rnorm(n_infants,
                                                mean = pr_helper,
                                                sd = pr_helper_sd)))
# Define informative priors for Bayesian models (sampling fails otherwise)
# Same Intercept whether raw or scaled age
# Same mildly informative prior on age_days whether raw or scaled age 
priors.full <- c(set_prior("normal(.5753641, .1)", # From Margoni & Surian, through logit
                           class = "Intercept"),
                 set_prior("normal(0, .5)",
                           class = "b"),
                 set_prior("student_t(3, 0, 2)",
                           class = "sd"))
priors.nointercept <- c(set_prior("normal(0, .5)",
                                  class = "b"),
                        set_prior("student_t(3, 0, 2)",
                                  class = "sd"))

# RUN STATISTICS: RAW AGE ==========================================================================
save_path <- "simulations_results/raw_age/"
# Run bayesian models and bridge-sample
run_models <- F
if(run_models){
  ## Full model
  ### Run models
  mod.raw_age.full <- brm(chose_helper ~ age_days +
                            (age_days | lab_id),
                          data = data_sims,
                          family = bernoulli,
                          prior = priors.full,
                          iter = 5000,
                          save_all_pars = T,
                          control = list(adapt_delta = .9999,
                                         max_treedepth = 20))
  ### Bridge-sample posteriors
  bridge.raw_age.full <- bridge_sampler(mod.raw_age.full, silent = T)
  ## No Intercept
  ### Run models
  mod.raw_age.nointercept <- brm(chose_helper ~ 0 + age_days +
                                   (0 + age_days | lab_id),
                                 data = data_sims,
                                 family = bernoulli,
                                 prior = priors.nointercept,
                                 iter = 5000,
                                 save_all_pars = T,
                                 control = list(adapt_delta = .9999,
                                                max_treedepth = 20))
  ### Bridge-sample posteriors
  bridge.raw_age.nointercept <- bridge_sampler(mod.raw_age.nointercept, silent = T)
  ## Save all
  saveRDS(mod.raw_age.full, paste0(save_path, "model_full.rds"))
  saveRDS(mod.raw_age.nointercept, paste0(save_path, "model_nointercept.rds"))
  saveRDS(bridge.raw_age.full, paste0(save_path, "bridge_full.rds"))
  saveRDS(bridge.raw_age.nointercept, paste0(save_path, "bridge_nointercept.rds"))
}else{
  ## Read all
  mod.raw_age.full <- readRDS(paste0(save_path, "model_full.rds"))
  mod.raw_age.nointercept <- readRDS(paste0(save_path, "model_nointercept.rds"))
  bridge.raw_age.full <- readRDS(paste0(save_path, "bridge_full.rds"))
  bridge.raw_age.nointercept <- readRDS(paste0(save_path, "bridge_nointercept.rds"))
}
beep(8)

# Get Bayes factors
bf.raw_age <- bayes_factor(bridge.raw_age.full,
                           bridge.raw_age.nointercept)

# RUN STATISTICS: SCALED AGE =======================================================================
save_path <- "simulations_results/scaled_age/"
# Run bayesian models and bridge-sample
run_models <- F
if(run_models){
  ## Full model
  ### Run models
  mod.scaled_age.full <- brm(chose_helper ~ z_age_days +
                               (z_age_days | lab_id),
                             data = data_sims,
                             family = bernoulli,
                             prior = priors.full,
                             iter = 5000,
                             save_all_pars = T,
                             control = list(adapt_delta = .9999,
                                            max_treedepth = 20))
  ### Bridge-sample posteriors
  bridge.scaled_age.full <- bridge_sampler(mod.scaled_age.full, silent = T)
  ## No Intercept
  ### Run models
  mod.scaled_age.nointercept <- brm(chose_helper ~ 0 + z_age_days +
                                      (0 + z_age_days | lab_id),
                                    data = data_sims,
                                    family = bernoulli,
                                    prior = priors.nointercept,
                                    iter = 5000,
                                    save_all_pars = T,
                                    control = list(adapt_delta = .9999,
                                                   max_treedepth = 20))
  ### Bridge-sample posteriors
  bridge.scaled_age.nointercept <- bridge_sampler(mod.scaled_age.nointercept, silent = T)
  ## Save all
  saveRDS(mod.scaled_age.full, paste0(save_path, "model_full.rds"))
  saveRDS(mod.scaled_age.nointercept, paste0(save_path, "model_nointercept.rds"))
  saveRDS(bridge.scaled_age.full, paste0(save_path, "bridge_full.rds"))
  saveRDS(bridge.scaled_age.nointercept, paste0(save_path, "bridge_nointercept.rds"))
}else{
  ## Read all
  mod.scaled_age.full <- readRDS(paste0(save_path, "model_full.rds"))
  mod.scaled_age.nointercept <- readRDS(paste0(save_path, "model_nointercept.rds"))
  bridge.scaled_age.full <- readRDS(paste0(save_path, "bridge_full.rds"))
  bridge.scaled_age.nointercept <- readRDS(paste0(save_path, "bridge_nointercept.rds"))
}
beep(8)

# Get Bayes factors
bf.scaled_age <- bayes_factor(bridge.scaled_age.full,
                              bridge.scaled_age.nointercept)
