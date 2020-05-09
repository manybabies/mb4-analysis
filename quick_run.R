t <- proc.time()
library(tidyverse);library(scales);library(beepr)
library(rstan);library(brms);library(bridgesampling);library(lme4)
library(future);library(future.apply);library(furrr)
plan(multiprocess, workers = 8) # Adapt to the number of cores you want to use
source("StatTools.R")
source("geom_flat_violin.R")

# Define simulation-specific variables ===========
n_labs <- 15
infants_by_lab <- 32
n_infants <- n_labs * infants_by_lab
age_min <- 165
age_max <- 320
pr_helper_mean <- 0.64  # From Margoni & Surian 2018 (estimated true effect size)
diff_pr_helper <- pr_helper_mean - .5
pr_helper_sd <- 0.1     # From Margoni & Surian 2018 (CI*sqrt(mean(N))/Z(95))

# Generate dataset ===============================
data_sim <- tibble(lab_id = rep(paste0("lab",1:n_labs),
                                each = infants_by_lab),
                   age_days = sample(age_min:age_max, n_infants,
                                     replace = T),
                   z_age_days = scale(age_days),
                   pr_helper_rn = rnorm(n_infants, pr_helper_mean, pr_helper_sd),
                   pr_helper = pmax(0, pmin(1, pr_helper_rn)),
                   chose_helper = rbinom(n_infants, 1, pr_helper))

# Set priors =====================================
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
priors.noage <- c(set_prior("normal(.5753641, .1)", # From Margoni & Surian, through logit
                           class = "Intercept"),
                 set_prior("student_t(3, 0, 2)",
                           class = "sd"))

# Raw age ==========================================================================================
save_path <- "simulation_results/raw_age/"
# Run simulations ================================
# Run bayesian models, bridge-sample, then
# save bf, parameter estimates, and HPDIs, to csv
# Running the models takes several hours
n_sims <- 8
run_models <- T
if(run_models){
  print("Initialise raw_age")
  ## Initialise filenames for estimates and bf
  filename.est <- paste0(save_path, "estimates_full.csv")
  filename.bf <- paste0(save_path, "bayes_factors.csv")
  ## Initialise simulation rng seeds
  if(file.exists(filename.est)){
    n_previous <- filename.est %>%
      read_csv(col_types = cols_only(sim_id = col_integer())) %>%
      max()
  }else{n_previous <- 0}
  ## Initialise estimates column names if file doesn't exist
  if(!file.exists(filename.est)){
    estimates <- tibble(`Parameter` = character(),
                        `Estimate` = double(),
                        `Est. Error` = double(),
                        `lower` = double(),
                        `upper` = double(),
                        `95% CI` = character(),
                        `sim_id` = character(),
                        .rows = 0)
    write_csv(estimates, filename.est)
  }
  ## Initialise models
  model.raw_age.full <- brm(chose_helper ~ age_days + (age_days | lab_id),
	                          data = data_sim, family = bernoulli,
	                          prior = priors.full, iter = 10000,
	                          future = T, save_all_pars = T,
	                          control = list(adapt_delta = .9999,
	                                         max_treedepth = 20))
  model.raw_age.nointercept <- brm(chose_helper ~ 0 + age_days + (0 + age_days | lab_id),
	                                 data = data_sim, family = bernoulli,
	                                 prior = priors.nointercept, iter = 10000,
	                                 future = T, save_all_pars = T,
	                                 control = list(adapt_delta = .9999,
	                                                max_treedepth = 20))
  ## Define single simulation function
  sim.raw_age <- function(seed){
    set.seed(seed)
    ### Generate dataset
    df <- tibble(lab_id = rep(paste0("lab",1:n_labs),
                              each = infants_by_lab),
                 age_days = sample(age_min:age_max, n_infants,
                                   replace = T),
                 z_age_days = scale(age_days),
                 pr_helper_rn = rnorm(n_infants, pr_helper_mean, pr_helper_sd),
                 pr_helper = pmax(0, pmin(1, pr_helper_rn)),
                 chose_helper = rbinom(n_infants, 1, pr_helper))
    ### Run models
    m.full <- update(model.raw_age.full, newdata = df, seed = seed)
    m.nointercept <- update(model.raw_age.nointercept, newdata = df, seed = seed)
    ### Get and save estimates and HPDIs for full model
    estimates.brm_fixef(m.full, prob = .95, digits = 7) %>%
      mutate(sim_id = seed) %>%
      write_csv(filename.est, append = T)
    ### Bridge sample
    bridge.full <- bridge_sampler(m.full, silent = T)
    bridge.nointercept <- bridge_sampler(m.nointercept, silent = T)
    ### Bayes factor
    bf <- bayes_factor(bridge.full, bridge.nointercept)$bf
    return(bf)
  }
  ## Get new Bayes factors
  print("Start raw_age")
  bf <- future_map_dbl(n_previous + 1:n_sims,
                       possibly(sim.raw_age, NA),
                       .progress = T)
  bf %>% enframe(name = NULL, value = "b.value") %>%
    write_csv(filename.bf, append = file.exists(filename.bf))
}

# Scaled age =======================================================================================
save_path <- "simulation_results/scaled_age/"
# Run simulations ================================
# Run bayesian models, bridge-sample, then
# save bf, parameter estimates, and HPDIs, to csv
# Running the models takes several hours
n_sims <- 8
run_models <- T
if(run_models){
  print("Initialise scaled_age")
  ## Initialise filenames for estimates and bf
  filename.est <- paste0(save_path, "estimates_full.csv")
  filename.bf <- paste0(save_path, "bayes_factors.csv")
  ## Initialise simulation rng seeds
  if(file.exists(filename.est)){
    n_previous <- filename.est %>%
      read_csv(col_types = cols_only(sim_id = col_integer())) %>%
      max()
  }else{n_previous <- 0}
  ## Initialise estimates column names if file doesn't exist
  if(!file.exists(filename.est)){
    estimates <- tibble(`Parameter` = character(),
                        `Estimate` = double(),
                        `Est. Error` = double(),
                        `lower` = double(),
                        `upper` = double(),
                        `95% CI` = character(),
                        `sim_id` = character(),
                        .rows = 0)
    write_csv(estimates, filename.est)
  }
  ## Initialise models
  model.scaled_age.full <- brm(chose_helper ~ z_age_days + (z_age_days | lab_id),
                          data = data_sim, family = bernoulli,
                          prior = priors.full, iter = 10000,
                          future = T, save_all_pars = T,
                          control = list(adapt_delta = .9999,
                                         max_treedepth = 20))
  model.scaled_age.nointercept <- brm(chose_helper ~ 0 + z_age_days + (0 + z_age_days | lab_id),
                          data = data_sim, family = bernoulli,
                          prior = priors.nointercept, iter = 10000,
                          future = T, save_all_pars = T,
                          control = list(adapt_delta = .9999,
                                         max_treedepth = 20))
  ## Define single simulation function
  sim.scaled_age <- function(seed){
    set.seed(seed)
    ### Generate dataset
    df <- tibble(lab_id = rep(paste0("lab",1:n_labs),
                              each = infants_by_lab),
                 age_days = sample(age_min:age_max, n_infants,
                                   replace = T),
                 z_age_days = scale(age_days),
                 pr_helper_rn = rnorm(n_infants, pr_helper_mean, pr_helper_sd),
                 pr_helper = pmax(0, pmin(1, pr_helper_rn)),
                 chose_helper = rbinom(n_infants, 1, pr_helper))
    ### Run models
    m.full <- update(model.scaled_age.full, newdata = df, seed = seed)
    m.nointercept <- update(model.scaled_age.nointercept, newdata = df, seed = seed)
    ### Get and save estimates and HPDIs for full model
    estimates.brm_fixef(m.full, prob = .95, digits = 7) %>%
      mutate(sim_id = seed) %>%
      write_csv(filename.est, append = T)
    ### Bridge sample
    bridge.full <- bridge_sampler(m.full, silent = T)
    bridge.nointercept <- bridge_sampler(m.nointercept, silent = T)
    ### Bayes factor
    bf <- bayes_factor(bridge.full, bridge.nointercept)$bf
    return(bf)
  }
  ## Get new Bayes factors
  print("Start scaled_age")
  bf <- future_map_dbl(n_previous + 1:n_sims,
                       possibly(sim.scaled_age, NA),
                       .progress = T)
  bf %>% enframe(name = NULL, value = "b.value") %>%
    write_csv(filename.bf, append = file.exists(filename.bf))
}

# Null age =========================================================================================
save_path <- "simulation_results/null_age/"
# Run simulations ================================
# Run bayesian models, bridge-sample, then
# save bf, parameter estimates, and HPDIs, to csv
# Running the models takes several hours
n_sims <- 8
run_models <- T
if(run_models){
  print("Initialise null_age")
  ## Initialise filenames for estimates and bf
  filename.est <- paste0(save_path, "estimates_full.csv")
  filename.bf <- paste0(save_path, "bayes_factors.csv")
  ## Initialise simulation rng seeds
  if(file.exists(filename.est)){
    n_previous <- filename.est %>%
      read_csv(col_types = cols_only(sim_id = col_integer())) %>%
      max()
  }else{n_previous <- 0}
  ## Initialise models (full model already initialised)
  model.scaled_age.noage <- brm(chose_helper ~ 1 + (1 | lab_id),
                                data = data_sim, family = bernoulli,
                                prior = priors.noage, iter = 10000,
                                future = T, save_all_pars = T,
                                control = list(adapt_delta = .9999,
                                               max_treedepth = 20))
  ## Define single simulation function
  sim.null_age <- function(seed){
    set.seed(seed)
    ### Generate dataset
    df <- tibble(lab_id = rep(paste0("lab",1:n_labs),
                              each = infants_by_lab),
                 z_age_days = sample(age_min:age_max, n_infants,
                                   replace = T),
                 z_z_age_days = scale(z_age_days),
                 pr_helper_rn = rnorm(n_infants, pr_helper_mean, pr_helper_sd),
                 pr_helper = pmax(0, pmin(1, pr_helper_rn)),
                 chose_helper = rbinom(n_infants, 1, pr_helper))
    ### Run models
    m.full <- update(model.scaled_age.full, newdata = df, seed = seed)
    m.noage <- update(model.scaled_age.noage, newdata = df, seed = seed)
    ### Bridge sample
    bridge.full <- bridge_sampler(m.full, silent = T)
    bridge.noage <- bridge_sampler(m.noage, silent = T)
    ### Bayes factor
    bf <- bayes_factor(bridge.full, bridge.noage)$bf
    return(bf)
  }
  ## Get new Bayes factors
  print("Start null_age")
  bf <- future_map_dbl(n_previous + 1:n_sims,
                       possibly(sim.null_age, NA),
                       .progress = T)
  bf %>% enframe(name = NULL, value = "b.value") %>%
    write_csv(filename.bf, append = file.exists(filename.bf))
}

# True null ========================================================================================
save_path <- "simulation_results/true_null/"
# Run simulations ================================
# Run bayesian models, bridge-sample, then
# save bf, parameter estimates, and HPDIs, to csv
# Running the models takes several hours
n_sims <- 8
run_models <- T
if(run_models){
  print("Initialise true_null")
  ## Initialise filenames for estimates and bf
  filename.est <- paste0(save_path, "estimates_full.csv")
  filename.bf <- paste0(save_path, "bayes_factors.csv")
  ## Initialise simulation rng seeds
  if(file.exists(filename.est)){
    n_previous <- filename.est %>%
      read_csv(col_types = cols_only(sim_id = col_integer())) %>%
      max()
  }else{n_previous <- 0}
  ## Initialise estimates column names if file doesn't exist
  if(!file.exists(filename.est)){
    estimates <- tibble(`Parameter` = character(),
                        `Estimate` = double(),
                        `Est. Error` = double(),
                        `lower` = double(),
                        `upper` = double(),
                        `95% CI` = character(),
                        `sim_id` = character(),
                        .rows = 0)
    write_csv(estimates, filename.est)
  }
  ## Define single simulation function
  sim.true_null <- function(seed){
    set.seed(seed)
    ### Generate dataset (null)
    df <- tibble(lab_id = rep(paste0("lab",1:n_labs),
                              each = infants_by_lab),
                 z_age_days = sample(age_min:age_max, n_infants,
                                   replace = T),
                 z_z_age_days = scale(z_age_days),
                 pr_helper_rn = rnorm(n_infants, .5, pr_helper_sd),
                 pr_helper = pmax(0, pmin(1, pr_helper_rn)),
                 chose_helper = rbinom(n_infants, 1, pr_helper))
    ### Run models (same as scaled_age sims, but null data above
    m.full <- update(model.scaled_age.full, newdata = df, seed = seed)
    m.nointercept <- update(model.scaled_age.nointercept, newdata = df, seed = seed)
    ### Get and save estimates and HPDIs for full model
    estimates.brm_fixef(m.full, prob = .95, digits = 7) %>%
      mutate(sim_id = seed) %>%
      write_csv(filename.est, append = T)
    ### Bridge sample
    bridge.full <- bridge_sampler(m.full, silent = T)
    bridge.nointercept <- bridge_sampler(m.nointercept, silent = T)
    ### Bayes factor
    bf <- bayes_factor(bridge.full, bridge.nointercept)$bf
    return(bf)
  }
  ## Get new Bayes factors
  print("Start true_null")
  bf <- future_map_dbl(n_previous + 1:n_sims,
                       possibly(sim.true_null, NA),
                       .progress = T)
  bf %>% enframe(name = NULL, value = "b.value") %>%
    write_csv(filename.bf, append = file.exists(filename.bf))
}
t <- proc.time() - t
print(t)
beep(8)
