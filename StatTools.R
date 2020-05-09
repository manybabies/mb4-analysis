library(brms)
library(coda)

# BRM FIXEF ESTIMATES
# Return estimates and HPDIs for the fixed effects from a brms model,
# for the given prob (default 0.89)
estimates.brm_fixef <- function(model, prob = 0.89, digits = 2){
  points <- model %>%
    summary() %>%
    {round(.$fixed[,1:2], digits = digits)} %>%
    as_tibble(rownames = "Parameter") %>%
    rename("Est. Error" = Est.Error)
  CI <- paste0(100*prob, "% CI")
  intervals <- model %>%
    as.mcmc(pars = "^b_", combine_chains = T) %>%
    HPDinterval(prob = prob) %>%
    round(digits = digits) %>%
    as_tibble(rownames = "Parameter") %>%
    mutate(Parameter = str_remove(Parameter, "^b_"),
           lower_str = paste0("[", lower),
           upper_str = paste0(upper, "]")) %>%
    unite(!!CI, lower_str, upper_str, sep = ", ")
  return(left_join(points, intervals))
}