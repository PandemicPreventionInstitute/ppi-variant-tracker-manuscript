# Author: Kaitlyn Johnson & Zack Susswein
# Purpose of this script is to load in the results from cmdstan from the 
# multicountry model and output datasets for visualization of the results

# Uses the time stamped data from July 1st, 2022

rm(list = ls())
USE_CASE <- 'local'
setwd("~/Documents/variant-tracker/ppi-variant-tracker/code")

# Libraries ---------------------------------------------------------------


library(cmdstanr) # command stan
library(tidybayes) # wrangle command stan objects
library(tidyverse) # data wrangling
library(tibble) # data wrangling
library(janitor) # column naming
library(countrycode) # country codes
library(lubridate) # date times
library(readxl) # excel import
library(zoo) # calculate rolling averages
library(R.utils) # R utilities
library(stringr) # to parse strings in R
library(dplyr) # data wrangling
library(readr) # read_csv
library(lme4) # logistic regression


# Load original dataset ---------------------------------------------------

# Path to the processed GISAID metadata file produced by pre_processing.R
DATA_PATH <- '../data/processed/lineage_t_2022-07-01.csv'
setwd("~/Documents/variant-tracker/ppi-variant-tracker/code")

# Numbers of draws from the posterior distribution to sample for Rt estimation
N_SAMPLED_DRAWS = 100
# Generation interval
MEAN_GI <- 5.8
TRANS_ADV_MULTIPLIER <- 7 # corresponds to the weekly transmission advantage 

# Load data ---------------------------------------------------------------

df.orig <- read_csv(DATA_PATH)


FIRST_DATE<-as.character(min(df.orig$collection_date))
LAST_COLLECTION_DATE <- as.character(max(df.orig$collection_date[df.orig$tot_seq>0]))
LAST_DATE_IN_DATASET <- as.character(max(df.orig$collection_date))


# Wrangle data ------------------------------------------------------------


df <- df.orig %>% 
  select(lineage, country, collection_date, n_seq, t) %>% 
  group_by(lineage, country, collection_date, t) %>% 
  summarize(n = sum(n_seq)) %>% 
  ungroup() %>% 
  droplevels() %>% 
  left_join(df.orig %>% select(lineage, country, collection_date,
                             mid_week_date, mid_week_p_lineage, mid_week_p_lineage_se)) # add these for global plotting
    


# Generate constants for data munging -------------------------------------

last_twenty_timepoints <- df %>% 
    arrange(collection_date) %>%
    pull(collection_date) %>%
    unique() %>% 
    tail(n=20)

most_prevalent_lineage <- df %>% 
    #filter(collection_date %in% last_twenty_timepoints) %>% 
    group_by(lineage) %>% 
    summarize(n = sum(n)) %>% 
    filter(n == max(n)) %>% 
    pull(lineage)

# in case there is more than one most prevalent lineage
most_prevalent_lineage <- most_prevalent_lineage[[1]]
print(most_prevalent_lineage)

# Map to Stan indices -----------------------------------------------------

# Map lineage names to indices
lineages <- df %>% 
  # Leave in the last day in the dataset as a dummy to allow pivoting
  filter(t == max(t)) %>% 
  select(lineage, t) %>% 
  unique() %>%
  # Pivot df to allow specifying first column by id
  pivot_wider(names_from = lineage, values_from = t) %>% 
  relocate((!!sym(most_prevalent_lineage)), .before = everything()) %>% 
  # Return dataframe to wide
  pivot_longer(cols = everything(), 
               names_to = 'lineage',
               values_to = 't') %>% 
  # Add in lineage id and comparison id (which is NA for the base case)
  mutate(lineage_id = row_number(),
         comparison_id = lineage_id - 1L,
         comparison_id = if_else(comparison_id < 1, 
                                 NA_integer_, 
                                 comparison_id),
         relation = paste0(lineage, '_rel_to_', most_prevalent_lineage)) %>% 
  # Remove dummy
  select(-t)


# Map country names to indices
countries <- df %>% 
  select(country) %>% 
  unique() %>% 
  mutate(country_id = row_number())



# Load fitted model -----------------------------------------------------------
# Using the time stamped dataset from July 1st for the manuscript
fit <- cmdstanr::as_cmdstan_fit(c('../data/output/multicountry_output/output_2022_07_011.csv', 
                                  '../data/output/multicountry_output/output_2022_07_012.csv',
                                  '../data/output/multicountry_output/output_2022_07_013.csv',
                                  '../data/output/multicountry_output/output_2022_07_014.csv'))

# Get tidy output for p_hat and Y_tilde -----------------------------------

all_draws<-fit$draws(format = "df") # get all the variables and all the draws

p_hats <- select(all_draws, contains('p_hat')) # subset to p_hats only, with country, t, and lineage_ids
sampled_p_hats <- p_hats %>% sample_n(N_SAMPLED_DRAWS) %>% # sample from the draws
    mutate(draw = row_number()) %>% # relabel the draws
    relocate(draw, .before = everything()) %>% # bring to first column
    pivot_longer(-draw) %>% # pivot
    mutate(indices = substring(name, 7, (str_length(name)-1))) %>% 
    separate(indices,
             into = c("country_id", "t", "lineage_id"),
             sep = ",") %>% 
    mutate(country_id = as.integer(country_id),
           t = as.integer(t),
           lineage_id = as.integer(lineage_id)) %>% 
    left_join(countries, by = "country_id") %>% 
    left_join(lineages, by = "lineage_id") %>% 
    rename(p_hat = value) %>% 
    select(- name, -comparison_id, -relation)

Y_tilde<- select(all_draws, contains('Y_tilde')) # subset to p_hats only, with country, t, and lineage_ids
sampled_Y_tildes <- Y_tilde %>% sample_n(N_SAMPLED_DRAWS) %>% # sample from the draws
    mutate(draw = row_number()) %>% # relabel the draws
    relocate(draw, .before = everything()) %>% # bring to first column
    pivot_longer(-draw) %>% # pivot
    mutate(indices = substring(name, 9, (str_length(name)-1))) %>% 
    separate(indices,
             into = c("country_id", "t", "lineage_id"),
             sep = ",") %>% 
    mutate(country_id = as.integer(country_id),
           t = as.integer(t),
           lineage_id = as.integer(lineage_id)) %>% 
    left_join(countries, by = "country_id") %>% 
    left_join(lineages, by = "lineage_id") %>% 
    rename(Y_tilde = value) %>% 
    select(- name, -comparison_id, -relation)



# Join the p_hat and Y_tilde together
tidy_fitted <- sampled_p_hats %>% left_join(sampled_Y_tildes) 




# Map model output times to dates
t <- seq(from = min(tidy_fitted$t), to = max(tidy_fitted$t), by = 1)
date <- seq(from = ymd(FIRST_DATE), to = ymd(FIRST_DATE) + days(length(t)-1), by = "days")
dates <- tibble(date = date,
                t = as.integer(t))
tidy_fitted <- tidy_fitted %>% left_join(dates, by = "t")


# Save summarized output --------------------------------------------------

# Dataframe containing all the draws 
out_df <- tidy_fitted %>% left_join(df) %>% 
    ungroup() %>% 
    group_by(t, country, draw) %>% 
    mutate(N = sum(n)) %>% 
    ungroup() %>% 
    mutate(obs_prev = n/N) %>%  
    ungroup() %>% # This differs for every draw
    # option to add the 7d average of predicted observed prevalence here
    select( country, lineage, t, date, n,  N, p_hat,
            draw, obs_prev)

out_df$obs_prev[out_df$N ==0] <- NA


# Check that prevalences all add to 1, ignoring days where there are no sequences (NAs)
p_test <-out_df %>% group_by(t, country, draw) %>%
    summarise(p_obs_test = sum(obs_prev),
              p_pred_test = sum(p_hat, na.rm = T),
              #model_obs_p = sum(pred_obs_prev),
              N = max(N))
stopifnot('All prevalences dont add to 1' = (round(max(p_test$p_obs_test, na.rm = T),4) ==1 & round(min(p_test$p_obs_test, na.rm = T),4) ==1))
stopifnot('All prevalences dont add to 1' = (round(max(p_test$p_pred_test),4) ==1 & round(min(p_test$p_pred_test),4) ==1))            

out_df <- out_df %>% select(country, lineage, t, date, n, N, p_hat,
                            draw, obs_prev) %>%  # contains daily (pred_obs_prev) and 7d average (pred_7d_obs_prev) of model predicted prevalence
                            mutate(last_collection_date = ymd(LAST_COLLECTION_DATE),
                                   reference_date = ymd(LAST_DATE_IN_DATASET),
                                   multicountry_period = ifelse(date <= LAST_COLLECTION_DATE,
                                                                'multicountry calibration', 'multicountry nowcast'))


# Make a dataframe that aggregate over draws with summary stats(best for plotting!)
pred_p <- out_df %>% group_by(date, country, lineage) %>% 
    summarise( p_hat_mean = mean(p_hat),
               p_hat_med = quantile(p_hat, probs = 0.5),
               p_hat_lb = quantile(p_hat, probs = 0.025),
               p_hat_ub = quantile(p_hat, probs = 0.975)) %>% distinct() %>% 
    ungroup() %>% 
    rename(collection_date = date)

# Add a test that ensures that at every time point the p_hat means for all variants sum to 1. 
test2 <- pred_p %>% group_by(collection_date, country) %>%
    summarise(summed_prev = sum(p_hat_mean))
stopifnot('For each draw and time point, lineage prevalences dont all add to 1' = (round(min(test2$summed_prev, na.rm = T), 4) == 1) &
              (round(max(test2$summed_prev, na.rm = T), 4) == 1))


# Save a cleaned country-lineage-time-summary-stats dataframe ----------------------------------------
clean_global_df<- pred_p %>% left_join(df %>% left_join(countries) %>% left_join(lineages)) %>% 
    ungroup() %>% 
    group_by(t, country) %>% 
    mutate(N= sum(n),
           collection_date = ymd(collection_date),
           mid_week_date = ymd(mid_week_date)) %>% 
    ungroup() %>% 
    mutate(last_collection_date = ymd(LAST_COLLECTION_DATE),
           reference_date = ymd(LAST_DATE_IN_DATASET),
           multicountry_period = ifelse(collection_date <= ymd(LAST_COLLECTION_DATE), 
                                        'multicountry calibration', 'multicountry nowcast'))




country_test<-c('Israel', 'Portugal', 'United States', 'Kenya', 'Bangladesh', 'Brazil', 'South Africa', 'Ghana', 'Algeria') 
# line plot
clean_global_df %>% filter(country =='Israel') %>% 
    ggplot()+
    geom_point(aes(mid_week_date, mid_week_p_lineage, color = lineage), size = 1)+
    geom_linerange(aes(x=mid_week_date,ymin = mid_week_p_lineage - mid_week_p_lineage_se, ymax = mid_week_p_lineage + mid_week_p_lineage_se, color = lineage))+
    geom_line(aes(x = collection_date, y = p_hat_mean, group = lineage, color = lineage))+
    geom_ribbon(aes(collection_date, ymin = p_hat_lb , ymax = p_hat_ub , group = lineage, fill = lineage), alpha = .2)+
    facet_wrap(~country)+
    theme_bw()+ xlab('Date') + ylab('Lineage prevalence')
# area plot
clean_global_df %>% filter(country %in% country_test) %>% 
    ggplot()+
    geom_area(aes(x = collection_date, y = p_hat_mean, group = lineage, fill = lineage))+
    facet_wrap(~country)+
    theme_bw()+ xlab('Date') + ylab('Lineage prevalence')

# Correlation matrix! 
Omega_hat <- fit %>% 
  spread_draws(Omega_hat[comparison_id_1, comparison_id_2]) %>% 
  mean_qi()  %>% 
  left_join(lineages %>% 
              rename(lineage_1 = lineage,
                     comparison_id_1 = comparison_id) %>% 
              select(-lineage_id, -relation)) %>% 
    left_join(lineages %>% 
              rename(lineage_2 = lineage,
                     comparison_id_2 = comparison_id) %>% 
              select(-lineage_id, -relation)) %>% 
    rename(Omega_hat_mean = Omega_hat,
           Omega_hat_ub = .upper,
           Omega_hat_lb = .lower)

# Posterior fitness advantage estimates for each variant country with all draws!
draws <- sample(unique(r_raw$.draw), N_SAMPLED_DRAWS)
r_raw <- fit %>% 
    spread_draws(r_hat[comparison_id, country_id])

draws <- sample(unique(r_raw$.draw), N_SAMPLED_DRAWS)

r_summary <- r_raw %>% 
    filter(.draw %in% draws) %>% 
    group_by(comparison_id, country_id) %>% 
    mutate(
        transmission_advantage = exp(r_hat*TRANS_ADV_MULTIPLIER) -1) %>% 
    median_qi(r_hat, transmission_advantage) %>% 
    left_join(countries) %>% 
    left_join(lineages %>% 
                  select(-lineage_id)) %>% 
    rename(r_median = r_hat,
           r_ub = r_hat.upper,
           r_lb = r_hat.lower,
           trans_adv_median = transmission_advantage,
           trans_adv_ub = transmission_advantage.upper,
           trans_adv_lb = transmission_advantage.lower) %>% 
    select(-.width, -.point, -.interval) %>% 
    ungroup() %>% 
    mutate(last_collection_date = ymd(LAST_COLLECTION_DATE),
           reference_date = ymd(LAST_DATE_IN_DATASET)) %>% 
    select(-comparison_id)

# This is the equivalent of r_out from the single country model (the one that concatenates all countries individual estimates)


# Global estimates------------------------------------------------------------
# fitness advantage of each variant across all countries
# summary
mu_hat <- fit %>% 
  spread_draws(mu_hat[comparison_id]) %>% 
  median_qi() %>% 
  left_join(lineages %>% 
              select(-lineage_id)) %>% 
    rename(mu_hat_median = mu_hat,
           mu_hat_ub = .upper,
           mu_hat_lb = .lower) %>% 
    select(-.width, -.point, -.interval)
# with draws
mu_all <- fit %>% 
    spread_draws(mu_hat[comparison_id]) %>% 
    left_join(lineages %>% 
                  select(-lineage_id)) %>% 
    rename(draw = .draw) %>% 
    mutate(transmission_advantage = exp(mu_hat*TRANS_ADV_MULTIPLIER) -1)
    #filter(draw %in% draws) 

mu_all %>% ggplot() + geom_density(aes(x = mu_hat,color = lineage)) + theme_bw()
mu_all %>% ggplot() + geom_density(aes(x = transmission_advantage, color = lineage)) + theme_bw()+
    geom_vline(aes(xintercept = 0), linetype = "dashed") +
    xlab(paste0('Transmission advantage over ', most_prevalent_lineage)) + ylab('Density') +
    ggtitle('Global posterior transmission advantage estimates')

# add summaries of fitness advantage to mu_hat
mu_hat <- mu_hat %>% left_join( mu_all %>% group_by (lineage) %>%  
                                    summarize( 
                                        transmission_adv_med = quantile(transmission_advantage, probs = 0.5),
                                        transmission_adv_lb = quantile(transmission_advantage, probs = 0.025),
                                        transmission_adv_ub = quantile(transmission_advantage, probs = 0.975))) %>% 
    mutate(last_collection_date = LAST_COLLECTION_DATE,
           reference_date = LAST_DATE_IN_DATASET)

# Add the global advantage to the country-specific advatnage summary df
r_combined<-r_summary %>% left_join(mu_hat %>% select(-last_collection_date, -reference_date, -comparison_id), by = c("lineage", "relation")) %>% 
    rename(global_trans_adv_median = transmission_adv_med,
           global_trans_adv_lb = transmission_adv_lb,
           global_trans_adv_ub = transmission_adv_ub)


# estimated y-intercept (initial prevalence) for each variant country
b0_hat <- fit %>% 
  spread_draws(b0[comparison_id, country_id]) %>% 
  mean_qi() %>% 
  left_join(countries) %>% 
  left_join(lineages %>% 
              select(-lineage_id)) %>% 
    rename(b0_hat_mean = b0,
           b0_hat_ub = .upper,
           b0_hat_lb = .lower) %>% 
    select(-.width, -.point, -.interval)

# Save output -------------------------------------------------------------

if (USE_CASE == 'local'){
  # Save the Y_tilde and p_hat draws in multicountry_model_pred_p.csv
  write.csv(out_df, '../data/processed/multicountry_model_pred_p_2022-07-01.csv', row.names = F)
  # Save the global fitness advantages estimates with all the draws
  write.csv(mu_all, '../data/processed/multicountry_mu_distrib_2022-07-01.csv', row.names = F)
  
  # Save the summarized quantities in individual csv files
  write.csv(r_combined, '../data/processed/r_summary_2022-07-01_main.csv', row.names = F)
  write.csv(clean_global_df, '../data/processed/clean_global_df_2022-07-01.csv', row.names = F)
  write.csv(Omega_hat, '../data/processed/Omega_hat_2022-07-01.csv', row.names = F)
  write.csv(mu_hat, '../data/processed/mu_hat_2022-07-01.csv', row.names = F)
  write.csv(b0_hat, '../data/processed/b0_hat_2022-07-01.csv', row.names = F)
  
}