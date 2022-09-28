# Author: Kaitlyn Johnson
# Modified by Zack Susswein
# Date initiated: 06-22-2022

# This script loads in the the multicountry model outputs from the runs on the 
# different reference datasets and pivots and summarizes them 


rm(list = ls())
USE_CASE = Sys.getenv("USE_CASE")
if(USE_CASE == ""){
    USE_CASE<-'local'
}
if (USE_CASE == 'local'){
    setwd("~/Documents/variant-tracker/ppi-variant-tracker/historical_validation_code")
}

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
library(DescTools) # Concordance correlation coefficient calculation


# Constants ---------------------------------------------------
# Numbers of draws from the posterior distribution to sample for Rt estimation
N_SAMPLED_DRAWS = 100
# Generation interval
MEAN_GI <- 5.8
TRANS_ADV_MULTIPLIER <- 7 # corresponds to the weekly transmission advantage 

if (USE_CASE == 'local'){
    REFERENCE_DATA_PATH <- '../data/processed/validation_data/reference_data_used.csv'
    VARIANT_T_PATH <- '../data/processed/variant_t_2022-07-01.csv'
}

# Load metadata on reference data ---------------------------------------------
reference_data_df <- read_csv(REFERENCE_DATA_PATH)
variant_t <- read_csv(VARIANT_T_PATH) 
variants_were_tracking <- unique(variant_t$variant)
NUM_DAYS_NOWCAST<-21
MLE_COUNTRIES <-c('Bangladesh', 'Brazil','Denmark', 'India', 'South Africa', 
                  'United Kingdom', 'United States', 'Senegal', 'Portugal') 


df_last.orig<-read_csv(paste0('../data/processed/validation_data/lineage_t_for_comp_', 
                              as.character(reference_data_df$reference_date[nrow(reference_data_df)]), '.csv'))

LAST_COLLECTION_DATE <- as.character(max(df_last.orig$collection_date[df_last.orig$tot_seq>0]))

df_last <- df_last.orig %>% 
    filter(collection_date <= ymd(LAST_COLLECTION_DATE)) %>% 
    select(lineage, country, collection_date, n_seq, t) %>% 
    group_by(lineage, country, collection_date, t) %>% 
    summarize(n = sum(n_seq)) %>% 
    ungroup() %>% 
    droplevels() %>% 
    rename(n_last = n,
           t_last = t) %>% 
    select(lineage, country, collection_date, n_last, t_last)








# Loop through all reference dates and load lineage data and model output -------


for (i in 1:nrow(reference_data_df)-1){
    THIS_REF_DATE<- as.character(reference_data_df$reference_date[i])
    if (USE_CASE == 'local'){
        DATA_PATH <- paste0('../data/processed/validation_data/lineage_t_', THIS_REF_DATE, '.csv')
    }
    df.orig <- read_csv(DATA_PATH)


    FIRST_DATE<-as.character(min(df.orig$collection_date))
    LAST_COLLECTION_DATE <- as.character(max(df.orig$collection_date[df.orig$tot_seq>0]))
    LAST_DATE_IN_DATASET <- as.character(max(df.orig$collection_date))


    # Wrangle data ------------------------------------------------------------


    df <- df.orig %>% 
      select(lineage, country, collection_date, n_seq, t) %>% 
      #filter(country %in% MLE_COUNTRIES) %>% 
      group_by(lineage, country, collection_date, t) %>% 
      summarize(n = sum(n_seq)) %>% 
      ungroup() %>% 
      droplevels() %>% 
      left_join(df.orig %>% select(lineage, country, collection_date,
                                   period, p_lineage_week, p_lineage_week_se, p_lineage, p_lineage_se,
                                 mid_week_date, mid_week_p_lineage, mid_week_p_lineage_se)) 
   
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

fit <- cmdstanr::as_cmdstan_fit(c(paste0('../data/output/multicountry_output/validation/output_', THIS_REF_DATE, '_1.csv'), 
                                  paste0('../data/output/multicountry_output/validation/output_', THIS_REF_DATE,'_2.csv'),
                                  paste0('../data/output/multicountry_output/validation/output_', THIS_REF_DATE, '_3.csv'),
                                  paste0('../data/output/multicountry_output/validation/output_', THIS_REF_DATE, '_4.csv')
            ))
    
   
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
    

    
    # Need to add in the future number of sequences on each day and of each lineage
    df_last_i <- df_last %>% filter(collection_date <= ymd(LAST_COLLECTION_DATE) + days(NUM_DAYS_NOWCAST),
                                    collection_date >=ymd(FIRST_DATE)) %>% 
        select(lineage, country, collection_date, n_last, t_last) %>% 
        mutate(t_last = t_last - min(t_last) + 1)  # Get the number of sequences of that lineage
    # observed on that day from 7-01
        
    df_comb<- df_last_i %>% left_join(df, by = c("country", "collection_date", "lineage", "t_last" = "t"))
    
    out_df <- tidy_fitted %>% 
        rename(collection_date = date) %>% 
        full_join(df_comb, by = c("country", "collection_date", "lineage")) %>% 
        group_by(t, collection_date, draw, country) %>% 
        mutate(N_last = sum(n_last, na.rm = T)) %>% 
        ungroup() %>% 
        mutate(pred_prev = Y_tilde/N_last) # gut check to make sure these look right! 
    
    

    
    
    # Calculate Brier score-----------------------------------------------------
    
    # lineage counts by country
    lin_country_counts <- out_df %>% group_by(country, lineage) %>% 
        summarise(tot_seq_that_lin = sum(n, na.rm = T)/N_SAMPLED_DRAWS)
    # could consider only evaluating based on lineages we're tracking i.e.
    # BA.4, BA.5, BA.2.12.1, BA.1, BA.2 
    
    lin_counts <- out_df %>% 
        #left_join(lin_country_counts, by = c('country', 'lineage'))%>%
        mutate(lineage_collapsed = ifelse( lineage %in% variants_were_tracking, lineage, 'other')) %>% 
        #mutate(lineage_collapsed = ifelse(tot_seq_that_lin>0, lineage, 'other')) %>% 
                group_by(t, collection_date, draw, country, lineage_collapsed) %>% 
        summarise(
            p_hat = sum(p_hat),
            Y_tilde = sum(Y_tilde), 
            n_lineages = sum(n_last)) %>% # This should be from the 07-02 dataset
        ungroup() %>% 
        group_by(country, t, draw) %>% 
        mutate(N_tot_seq_day = sum(n_lineages, na.rm = T),
               sum_p_hat = sum(p_hat)) %>% # check it adds to 1
        ungroup() %>% 
        select(country, lineage_collapsed, t, collection_date, draw, p_hat, n_lineages, Y_tilde, N_tot_seq_day) %>% 
        mutate(
            multicountry_period =case_when( 
                collection_date <= ymd(LAST_COLLECTION_DATE) ~ 'multicountry calibration',
                collection_date > ymd(LAST_COLLECTION_DATE) ~ 'multicountry forecast')) %>% 
        group_by(country, draw) %>% 
        mutate(N_tot_seq = sum(n_lineages, na.rm = T)) %>% 
        ungroup() %>% 
        group_by(country, multicountry_period, draw) %>% 
        mutate(N_tot_seq_period = sum(n_lineages, na.rm = T)) %>% 
        ungroup()
    lin_counts$n_lineages[is.na(lin_counts$n_lineages)]<-0
    lin_counts$Y_tilde[lin_counts$N_tot_seq_day==0]<-0
    
    
    
    # Brier score by country for both calibration and forecast period
    BS_all <- lin_counts %>% 
        group_by(draw, country,  lineage_collapsed) %>% 
        summarise(
            sum_over_t = sum( n_lineages - 2*p_hat*n_lineages + N_tot_seq_day*(p_hat^2)), # using sufficient stat for for (p_hat - Y)^2
            N_tot_seq = max(N_tot_seq),
        ) %>% ungroup() %>% 
        group_by(country, draw) %>% 
        summarise( Brier_score = 1/max(N_tot_seq)*sum(sum_over_t)) %>% 
        group_by(country) %>% 
        summarise(BS_med_MC = quantile(Brier_score, probs = 0.5, na.rm = T),
                  BS_mean_MC = mean(Brier_score, na.rm = T),
                  BS_lb_MC = quantile(Brier_score, probs = 0.025, na.rm = T),
                  BS_ub_MC = quantile(Brier_score, probs = 0.975, na.rm = T))
    
    BS_all %>% filter(country %in% MLE_COUNTRIES) %>% ggplot() + 
        geom_point(aes(x = country, y = BS_med_MC)) +
        geom_linerange(aes(x = country, ymin = BS_lb_MC, ymax = BS_ub_MC))
    
    # Brier score by country broken down by calibration vs forecast period 
    BS_by_period <- lin_counts %>% 
        group_by(draw, country, lineage_collapsed, multicountry_period) %>% 
        summarise(
            sum_over_t = sum( n_lineages - 2*p_hat*n_lineages + N_tot_seq_day*(p_hat^2), na.rm = T),
            N_tot_seq_period = max(N_tot_seq_period)
        ) %>% ungroup() %>% 
        group_by(country, draw, multicountry_period) %>% 
        summarise(Brier_score_period = 1/max(N_tot_seq_period)*sum(sum_over_t, na.rm = T)) %>% 
        group_by(country, multicountry_period) %>% 
        summarise(BS_med_period_MC = quantile(Brier_score_period, probs = 0.5, na.rm = T),
                  BS_mean_period_MC = mean(Brier_score_period, na.rm = T),
                  BS_lb_period_MC = quantile(Brier_score_period, probs = 0.025, na.rm = T),
                  BS_ub_period_MC = quantile(Brier_score_period, probs = 0.975, na.rm = T))

    BS<- BS_by_period %>% left_join(BS_all, by = 'country') %>% 
        mutate(reference_date = THIS_REF_DATE)
    
    
    # CCC by country 
    CCC_all <- lin_counts %>% 
        group_by(draw, country) %>% 
        summarise(
            CCC_by_country = CCC(n_lineages, Y_tilde)$rho.c$est
        ) %>% 
        group_by(country) %>% 
        summarise(
            CCC_country_med_MC = quantile(CCC_by_country, probs = 0.5, na.rm = T),
            CCC_country_mean_MC = mean(CCC_by_country, na.rm = T),
            CCC_country_lb_MC = quantile(CCC_by_country, probs = 0.025, na.rm = T),
            CCC_country_ub_MC = quantile(CCC_by_country, probs = 0.975, na.rm = T)
        )
    
    
    # CCC by country and variant 
    CCC_country_lin <- lin_counts %>% 
        group_by(draw, country, lineage_collapsed) %>% 
        summarise(
            CCC_by_country_lineage = CCC(n_lineages, Y_tilde)$rho.c$est
        ) %>% 
        group_by(country, lineage_collapsed) %>% 
        summarise(
            CCC_country_lin_med_MC = quantile(CCC_by_country_lineage, probs = 0.5, na.rm = T),
            CCC_country_lin_mean_MC = mean(CCC_by_country_lineage, na.rm = T),
            CCC_country_lin_lb_MC = quantile(CCC_by_country_lineage, probs = 0.025, na.rm = T),
            CCC_country_lin_ub_MC = quantile(CCC_by_country_lineage, probs = 0.975, na.rm = T)
        )
    
   
    # CCC by country, period, and variant
    CCC_lin_period <- lin_counts %>% 
        group_by(draw, country, lineage_collapsed, multicountry_period) %>% 
        summarise(
            CCC_by_country_lineage_period = CCC(n_lineages, Y_tilde)$rho.c$est
        ) %>% ungroup() %>% 
        group_by(country, lineage_collapsed, multicountry_period) %>% 
        summarise(
            CCC_country_lin_period_med_MC = quantile(CCC_by_country_lineage_period, probs = 0.5, na.rm = T),
            CCC_country_lin_period_mean_MC = mean(CCC_by_country_lineage_period, na.rm = T),
            CCC_country_lin_period_lb_MC = quantile(CCC_by_country_lineage_period, probs = 0.025, na.rm = T),
            CCC_country_lin_period_ub_MC = quantile(CCC_by_country_lineage_period, probs = 0.975, na.rm =T)
        )
    
    
    # CCC by country and period 
    CCC_period <- lin_counts %>% 
        group_by(draw, country, multicountry_period) %>% 
        summarise(
            CCC_by_country_period = CCC(n_lineages, Y_tilde)$rho.c$est
        ) %>% ungroup() %>% 
        group_by(country, multicountry_period) %>% 
        summarise(
            CCC_country_period_med_MC = quantile(CCC_by_country_period, probs = 0.5, na.rm = T),
            CCC_country_period_mean_MC = mean(CCC_by_country_period, na.rm = T),
            CCC_country_period_lb_MC = quantile(CCC_by_country_period, probs = 0.025, na.rm = T),
            CCC_country_period_ub_MC = quantile(CCC_by_country_period, probs = 0.975, na.rm = T)
        )
    
    # Combine the country level metrics
    country_metrics <- BS %>% left_join(CCC_all, by = "country") %>% 
        left_join(CCC_period, by = c("country", "multicountry_period")) %>% 
        relocate(reference_date, .before = multicountry_period)
    
    #Combine the country variant level CCC metrics
    country_lineage_metrics <- CCC_lin_period %>% 
        left_join(CCC_country_lin, by = c('country', 'lineage_collapsed')) %>% 
        mutate (reference_date = THIS_REF_DATE) %>% 
        rename(lineage= lineage_collapsed) %>% 
        relocate(reference_date, .before = multicountry_period)
    
    
    # Posterior fitness advantage estimates for each variant country with all draws!
    r_raw <- fit %>% 
        spread_draws(r_hat[comparison_id, country_id])
        
    draws <- sample(unique(r_raw$.draw), N_SAMPLED_DRAWS)
    
    # Save the country level r distribs here 
    r_distrib <- r_raw %>% filter(.draw %in% draws) %>% 
        group_by(comparison_id, country_id) %>% 
        mutate(
            transmission_advantage = exp(r_hat*TRANS_ADV_MULTIPLIER) -1) %>% 
        ungroup() %>% 
        left_join(countries) %>% 
        left_join(lineages %>% 
                      select(-lineage_id)) %>% 
        select( -country_id, -comparison_id, -.chain, -.iteration) %>% 
        rename(draw = .draw) %>% 
        relocate(country, .before = everything()) %>% 
        mutate(reference_date = THIS_REF_DATE)
    
    # save the summary stats 
    r_summary<- r_raw %>% filter(.draw %in% draws) %>% 
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
               last_submission_date = ymd(LAST_DATE_IN_DATASET),
               reference_date = ymd(THIS_REF_DATE)) %>% 
        select(-comparison_id)
        
    


# Global estimates------------------------------------------------------------

    # fitness advantage of each variant across all countries
    # summary
    mu_hat <- fit %>% 
      spread_draws(mu_hat[comparison_id]) %>% 
        filter(.draw %in% draws) %>% 
      median_qi() %>% 
      left_join(lineages %>% 
                  select(-lineage_id)) %>% 
        rename(mu_hat_median = mu_hat,
               mu_hat_ub = .upper,
               mu_hat_lb = .lower) %>% 
        select(-.width, -.point, -.interval)
    # get all draws of the global transmission advantage 
    mu_all <- fit %>% 
        spread_draws(mu_hat[comparison_id]) %>% 
        filter(.draw %in% draws) %>% 
        left_join(lineages %>% 
                      select(-lineage_id)) %>% 
        rename(draw = .draw) %>% 
        mutate(global_transmission_advantage = exp(mu_hat*TRANS_ADV_MULTIPLIER) -1) %>% 
        ungroup() %>% 
        select(-.chain, -.iteration, -comparison_id)


    # add summaries of fitness advantage to mu_hat
    mu_hat <- mu_hat %>% left_join( mu_all %>% group_by (lineage) %>%  
                                        summarize( 
                                            global_trans_adv_median = quantile(global_transmission_advantage, probs = 0.5),
                                            global_trans_adv_lb = quantile(global_transmission_advantage, probs = 0.025),
                                            global_trans_adv_ub = quantile(global_transmission_advantage, probs = 0.975))) %>% 
        mutate(last_collection_date = ymd(LAST_COLLECTION_DATE),
               last_submission_date = ymd(LAST_DATE_IN_DATASET),
               reference_date = ymd(THIS_REF_DATE))

    # Add the global advantage to the country-specific advatnage summary df
    r_combined<-r_summary %>% left_join(mu_hat %>% select(-last_collection_date, - reference_date, -last_submission_date), by = c("lineage", "relation")) 
    
    # Add the draws of the global advantage to the country specific (r_distrib)
    r_distrib <- r_distrib %>% left_join (mu_all, by = c("lineage", "relation", "draw"))
    mu_all <- mu_all %>% mutate(reference_date = THIS_REF_DATE)
    

# Save output -------------------------------------------------------------

    if (USE_CASE == 'local'){
      # Save the summarized quantities in individual csv files
      write.csv(country_metrics, paste0('../data/processed/validation_data/country_metrics_', THIS_REF_DATE, '.csv'), row.names = F)
      write.csv(country_lineage_metrics, paste0('../data/processed/validation_data/country_lineage_metrics_', THIS_REF_DATE, '.csv'), row.names = F)
      write.csv(r_combined, paste0('../data/processed/validation_data/r_summary_', THIS_REF_DATE, '.csv'), row.names = F)
      write.csv(r_distrib, paste0('../data/processed/validation_data/r_distrib_', THIS_REF_DATE, '.csv'),  row.names = F)
      write.csv(mu_all, paste0('../data/processed/validation_data/mu_distrib_', THIS_REF_DATE, '.csv'),  row.names = F)
      write.csv(mu_hat, paste0('../data/processed/validation_data/mu_hat_', THIS_REF_DATE, '.csv'), row.names = F)
    }
} # end loop over all reference dates 
