

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(nnet)

# Data --------------------------------------------------------------------

dat <- read_csv('../data/processed/data_for_nnet.csv')

# Constants ---------------------------------------------------------------

COUNTRIES = unique(dat$country)

# Number of bootstrap replicates
N_BOOT = 5e3

# fit model ---------------------------------------------------------------

# Set up vector of observed times
X <- as.integer(colnames(dat)[3:ncol(dat)])

# Check that there's nothing weird with the col names
stopifnot(typeof(X) == 'integer')
stopifnot(min(X) == 1)
stopifnot(max(X) == length(X))
stopifnot(!any(is.na(X)))

# Set X as a data frame
X <- data.frame("t" = X)

# Set up tibble to hold summarized results
full_results <- tibble(lineage = character(),
                       t = integer(),
                       country = character(),
                       p_hat_MLE_mean = double(),
                       p_hat_MLE_ub = double(),
                       p_hat_MLE_lb = double())

# Iterate through the countries to fit the multinomial MLE model to as well as 
# bootstrap confidence intervals
for(COUNTRY in COUNTRIES){
  
  # Y values are the matrix of observed counts
  Y <- dat %>% 
    filter(country == COUNTRY) %>% 
    select(-lineage, -country) %>% 
    as.matrix()
  
  colnames(Y) <- NULL
  
  # nnet::multinom requires that there be at least one observed sequence on a 
  # given day and behaves poorly if a lineage is never observed. Correct for 
  # this by removing lineages and days that are never observed.
  
  # First, save the days and lineages that are being fit to
  vars_analyzed <- unique(dat$lineage)[rowSums(Y) > 0]
  t_analyzed = colSums(Y) > 0
  
  # Remove days and lineages with no observations
  Y <- Y[rowSums(Y) > 0,]
  X <- as.matrix(X)
  X <- X[colSums(Y) > 0,]
  Y <- Y[,colSums(Y) > 0]
  
  # Fit the multinomial regression
  fit <- nnet::multinom(t(Y) ~ 1 + X)
  
  # Check that model has converged
  stopifnot(fit$convergence == 0)
  
  # Set up the array to hold the predictions
  preds <- tibble(cat = integer(),
                  t = integer(),
                  rep = integer(),
                  pred_p = double())
  
  # Set up stuff that will be re-used a bunch of times in the bootstrap fits
  base_case <-  tibble(intercept = 1, 
                       r = 0)
  
  lineage_time_grid <- expand_grid(cat = 1:nrow(Y),
                                   t = 1:length(t_analyzed)) %>% 
    left_join(tibble(lineage = unique(dat$lineage)) %>% 
                mutate(cat = row_number()))
  
  
  # Fit the bootstrap replicates
  for (i in 1:N_BOOT){
    
    # Sample with replacement the indices for the bootstrap fit
    boot_indices <- sample(1:length(X), size = length(X), replace = T)
    
    # Pull out the X and Y values for the bootstrap fit
    X_boot <- X[boot_indices]
    Y_boot <- Y[,boot_indices]
    
    # Fit the model, silencing printed output so that it doesn't overwhelm
    # the console. Only works on *nix & will throw an error on windows
    capture.output(fit_boot <- nnet::multinom(t(Y_boot) ~ 1 + X_boot), 
                   file = '/dev/null')
    
    # Save the output of the bootstrap replicate
    mean.rep <- base_case %>% 
      rbind(tibble(intercept = summary(fit_boot)$coefficients[,1],
                   r = summary(fit_boot)$coefficients[,2])) %>% 
      mutate(cat = row_number()) %>% 
      full_join(lineage_time_grid,
                by = "cat") %>% 
      mutate(pred_linear = intercept + r * t) %>% 
      group_by(t) %>% 
      mutate(pred_p = exp(pred_linear) / sum(exp(pred_linear))) %>% 
      ungroup() %>% 
      select(cat, t, pred_p) %>% 
      mutate(rep = i)
    
    # Append the replicate's predictions to the tibble of predictions
    preds <- rbind(preds, mean.rep)
    
    # Print the percentage completed if a multiple of 100
    if (i %% 100 == 0){
        
        print((i / N_BOOT) * 100)
    }
    
  }
  
  # Generate the full model fit predictions
  base_fit <- base_case %>% 
      rbind(tibble(intercept = summary(fit)$coefficients[,1],
                   r = summary(fit)$coefficients[,2])) %>% 
      mutate(cat = row_number()) %>% 
      full_join(lineage_time_grid,
                by = c("cat")) %>% 
      mutate(pred_linear = intercept + r * t) %>% 
      group_by(t) %>% 
      mutate(p_hat_MLE_mean = exp(pred_linear) / sum(exp(pred_linear))) %>% 
      ungroup() %>% 
      select(cat, t, p_hat_MLE_mean)
  
  # Save the bootstrap CIs appending the mean of the full data fit
  boot_se <- preds %>% 
      group_by(cat, t) %>% 
      summarize(p_hat_MLE_lb = quantile(pred_p, 0.025),
                p_hat_MLE_ub = quantile(pred_p, 0.975),
                .groups = 'drop') %>% 
      left_join(tibble(lineage = vars_analyzed,
                       cat = 1:length(vars_analyzed)),
                by = 'cat') %>% 
      left_join(base_fit,
                by = c('cat', 't')) %>% 
    mutate(country = COUNTRY) %>% 
    select(-cat)
  
  full_results <- rbind(full_results, boot_se)
  
}


# Save output -------------------------------------------------------------

# KAITLYN: This won't be good enough because we need to change this by
# date, but I'm not sure what the right structure is.
write_csv(full_results, paste0('../data/mle/bootstrap_results.csv'))


# Plot output

#dat %>% 
#    pivot_longer(colnames(dat)[3]:colnames(dat)[ncol(dat)],
#                 names_to = 't',
#                 values_to = 'n') %>% 
#    mutate(t = as.integer(t)) %>% 
#  filter(country == 'Brazil',
#         lineage %in% vars_analyzed) %>% 
#  group_by(t) %>% 
#  mutate(p_obs = n / sum(n)) %>% 
#  ungroup() %>%
#  left_join(boot_se) %>% 
#  ggplot()+
#  geom_point(aes(t, p_obs, color = lineage))+
#  geom_line(aes(t, pred_p, color = lineage,  group = lineage))+
#  geom_ribbon(aes(x = t,
#                  ymin = .lower,
#                  ymax = .upper,
#                  fill = lineage),
#              alpha = 0.2)+
#    facet_wrap(~lineage, scale = 'free_y')
#