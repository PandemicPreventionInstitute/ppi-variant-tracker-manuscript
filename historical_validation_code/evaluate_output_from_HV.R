# Author: Kaitlyn Johnson
# Modified from Zack Susswein
# Date initiated: 06-23-2022

# This script loads in the results from all of the reference datasets for the 
# multicountry and the MLE models


rm(list = ls())
USE_CASE = Sys.getenv("USE_CASE")
if(USE_CASE == ""){
    USE_CASE<-'local'
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
library(gridExtra)
library(grid)


# Constants ---------------------------------------------------
if (USE_CASE == 'local'){
    REFERENCE_DATA_PATH <- '../data/processed/validation_data/reference_data_used.csv'
    MOST_RECENT_R_SUMMARY_PATH <- '../data/processed/validation_data/r_summary_2022-07-01.csv'
    MOST_RECENT_CLEAN_GLOBAL_DF_PATH <- '../data/processed/validation_data/clean_global_df_2022-07-01.csv'
    MOST_RECENT_MU_HAT_PATH <- '../data/processed/validation_data/mu_hat_2022-07-01.csv'
    MLE_FIT_PATH <- '../data/processed/validation_data/MLE_t.csv' # contains subset of countries, all data
    MLE_REGRESSORS_PATH <- '../data/processed/validation_data/MLE_regressors.csv'
    COUNTRIES_FIT_TO_MLE <- '../data/processed/validation_data/countries_fit_to_MLE.csv'
    MLE_METRICS_PATH <- '../data/processed/validation_data/country_metrics_MLE.csv'
    MLE_LINEAGE_METRICS_PATH <- '../data/processed/validation_data/country_lineage_metrics_MLE.csv'
    #setwd("~/Documents/variant-tracker/ppi-variant-tracker/historical_validation_code")
} 

# Load reference dates & most recent data --------------------------------------
reference_data_df <- read_csv(REFERENCE_DATA_PATH)
r_summary_recent <- read_csv(MOST_RECENT_R_SUMMARY_PATH)
clean_global_df_recent <- read.csv(MOST_RECENT_CLEAN_GLOBAL_DF_PATH)
mu_hat_recent <- read_csv(MOST_RECENT_MU_HAT_PATH)
mle_t <- read_csv(MLE_FIT_PATH)
mle_regressors <- read_csv(MLE_REGRESSORS_PATH)
mle_countries <- read.csv(COUNTRIES_FIT_TO_MLE)
country_list<-mle_countries$x
mle_metrics <- read_csv(MLE_METRICS_PATH) %>% 
    relocate(country, .before = everything())
mle_lin_metrics <- read_csv(MLE_LINEAGE_METRICS_PATH) %>% 
    relocate(country, .before = everything())

# Make new column names for the reference data
clean_global_df_recent <- clean_global_df_recent %>% 
    rename( # rename a bunch of columns so we can join it to the old dataset
        mid_week_date_recent = mid_week_date,
        mid_week_prev_recent = mid_week_p_lineage,
        mid_week_prev_se_recent = mid_week_p_lineage_se,
        tot_seq_recent = N, # daily numbers of sequences
        n_recent = n,
        p_lineage_week_recent = p_lineage_week, 
        p_lineage_week_se_recent = p_lineage_week_se,
        p_lineage_recent = p_lineage,
        p_lineage_se_recent = p_lineage_se) %>% 
    filter(multicountry_period == 'multicountry calibration') %>% 
    select(country, lineage, collection_date,
        mid_week_date_recent,mid_week_prev_recent, mid_week_prev_se_recent,
        p_lineage_week_recent,p_lineage_week_se_recent, p_lineage_recent,
        p_lineage_se_recent)

r_summary_recent <- r_summary_recent %>% 
  rename(trans_adv_median_recent = trans_adv_median,
         global_trans_adv_median_recent = global_trans_adv_median) 
  


r_summary <-c()
clean_global_df <- c()
mu_hat <-c()
country_metrics <-c()
country_lineage_metrics <- c()
mu_distrib <- c()
r_distrib <- c()
# remove the last dataset because we will use this for comparison
reference_data_df <- reference_data_df %>% filter(reference_date != ymd("2022-07-01")) 

# Loop through the outputs from the reference datasets-----------
for (i in 1:nrow(reference_data_df)){
    THIS_REF_DATE<- as.character(reference_data_df$reference_date[i])
    
    # Load in the datasets from that reference date for both model types
    if (USE_CASE == 'local'){
        this_r_summary <- read_csv(paste0('../data/processed/validation_data/r_summary_', THIS_REF_DATE, '.csv'))  
        this_clean_global_df <- read.csv(paste0('../data/processed/validation_data/clean_global_df_', THIS_REF_DATE, '.csv')) 
        this_mu_hat <- read_csv(paste0('../data/processed/validation_data/mu_hat_', THIS_REF_DATE, '.csv'))
        this_country_metrics <- read_csv(paste0('../data/processed/validation_data/country_metrics_', THIS_REF_DATE, '.csv'))
        this_country_lineage_metrics <- read_csv(paste0('../data/processed/validation_data/country_lineage_metrics_', THIS_REF_DATE, '.csv')) %>% 
            filter(country %in% country_list)
        this_mu_distrib <- read_csv(paste0('../data/processed/validation_data/mu_distrib_', THIS_REF_DATE, '.csv'))
        this_r_distrib <- read_csv(paste0('../data/processed/validation_data/r_distrib_', THIS_REF_DATE, '.csv')) %>% 
            filter(country %in% country_list)
    }

    

        # Left join the MLE estimated prevalence over time 
        this_mle_t <- mle_t %>% filter(reference_date == THIS_REF_DATE) %>% select(-reference_date)
        this_clean_global_df <- left_join(this_clean_global_df, this_mle_t, by = c("country", "lineage", "t"))
        
        # Left join the Brier score and CCC metrics 
        this_mle_metrics <- mle_metrics %>% filter(reference_date == THIS_REF_DATE) %>% select(-reference_date)
        this_country_metrics <- left_join(this_country_metrics, this_mle_metrics, by = c("country", "multicountry_period"))
        this_mle_lineage_metrics <- mle_lin_metrics %>% filter(reference_date == THIS_REF_DATE) %>% select(-reference_date)
        this_lineage_metrics <- left_join(this_country_lineage_metrics, this_mle_lineage_metrics, by = c("country", "multicountry_period", "lineage"))
        
        
        
        # Left join the current data to the historical model for direct comparison (gets rid of lots of the model)
        this_clean_global_df <- clean_global_df_recent %>% select(country, lineage, collection_date, mid_week_date_recent,
                                                                  mid_week_prev_recent, mid_week_prev_se_recent,
                                                                  p_lineage_recent) %>% 
            left_join(this_clean_global_df,  by = c("lineage", "country", "collection_date")) %>% 
            mutate(collection_date = ymd(collection_date),
                   reference_date = ymd(THIS_REF_DATE),
                   last_collection_date = ymd(this_clean_global_df$last_collection_date[1]))
        
        this_r_summary <- r_summary_recent %>% select(country, lineage, trans_adv_median_recent, global_trans_adv_median_recent) %>% 
          left_join(this_r_summary, by = c("country", "lineage")) %>% 
          mutate(reference_date = ymd(THIS_REF_DATE), 
                 num = i)

       
        
        # Concatenate the reference dataset outputs
        r_summary <- bind_rows(r_summary, this_r_summary)
        clean_global_df <- bind_rows(clean_global_df, this_clean_global_df)
        mu_hat <- bind_rows(mu_hat, this_mu_hat)
        country_metrics <- bind_rows(country_metrics, this_country_metrics)
        country_lineage_metrics <- bind_rows(country_lineage_metrics, this_lineage_metrics)
        mu_distrib <- bind_rows(mu_distrib, this_mu_distrib)
        r_distrib <- bind_rows(r_distrib, this_r_distrib)


}



# Join MLE fitness advantage estimates by country (again there will be gaps)

r_summary <- r_summary %>% left_join(mle_regressors, by = c("reference_date", "lineage", "country")) %>% 
  select(- transmission_adv_med, - transmission_adv_lb, - transmission_adv_ub) %>% 
  mutate(dif_trans_adv_MC = trans_adv_median - trans_adv_median_recent,
         dif_global_trans_adv_MC = trans_adv_median - global_trans_adv_median,
         dif_trans_adv_MLE = trans_adv_MLE_mean - trans_adv_median_recent,
         dif_global_trans_adv_MLE = trans_adv_MLE_mean -global_trans_adv_median,
         pct_dif_trans_adv_MC = 100*(trans_adv_median - trans_adv_median_recent)/ trans_adv_median_recent,
         pct_dif_global_trans_adv_MC = 100*(trans_adv_median - global_trans_adv_median)/ global_trans_adv_median,
         pct_dif_trans_adv_MLE = 100*(trans_adv_MLE_mean - trans_adv_median_recent)/trans_adv_MLE_mean,
         pct_dif_global_trans_adv_MLE = 100*(trans_adv_MLE_mean -global_trans_adv_median)/global_trans_adv_median)

# Compute difference between MLE and multicountry metrics
country_early <- country_metrics %>% 
  filter(reference_date == reference_data_df$reference_date[1]) %>% 
    mutate(BS_improvement = BS_med_MLE - BS_med_MC,
           CCC_improvement = CCC_country_med_MC - CCC_country_med_MLE,
           BS_period_improvement = BS_med_period_MLE - BS_med_period_MC,
           CCC_period_improvement = CCC_country_period_med_MC - CCC_country_period_med_MLE) %>% 
    unique() %>% 
    arrange(desc(BS_improvement)) %>% 
  relocate(BS_improvement, .before = multicountry_period)

# Find the number of sequences per country lineage and join it to the r_summary data
n_seq_country_lineage<- clean_global_df %>% group_by(country, lineage, reference_date) %>%
    summarise(n_seq_lineage_country = sum(n, na.rm = T))
# Join with r_summary
r_summary <- r_summary %>% left_join(n_seq_country_lineage, by = c('reference_date', 'country', 'lineage'))

# Save all the files 
if (USE_CASE == 'local'){
  write.csv(r_summary, '../data/output/validation/r_summary.csv', row.names = F)
  write.csv(clean_global_df, '../data/output/validation/clean_global_df.csv', row.names = F)
  write.csv(mu_hat, '../data/output/validation/mu_hat.csv', row.names = F)
  write.csv(mu_distrib, '../data/output/validation/mu_distrib.csv', row.names = F)
  write.csv(r_distrib, '../data/output/validation/r_distrib.csv', row.names = F)
  write.csv(country_metrics, '../data/output/validation/country_metrics.csv', row.names = F)
  write.csv(country_lineage_metrics, '../data/output/validation/country_lineage_metrics.csv', row.names = F)
}







