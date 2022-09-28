# Author: Kaitlyn Johnson
# Modified from Zack Susswein
# Date initiated: 06-22-2022

# This script loads in the stansummary prevalence estimates and fitness 
# advantages from each reference dataset


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


# Constants ---------------------------------------------------
# Generation interval
#MEAN_GI <- 5.8
#TRANS_ADV_MULTIPLIER <-7 # corresponds to weekly transmission advantage 
if (USE_CASE == 'local'){
    REFERENCE_DATA_PATH <- '../data/processed/validation_data/reference_data_used.csv'
}

# Load metadata on reference data ---------------------------------------------
reference_data_df <- read_csv(REFERENCE_DATA_PATH)


# Loop through all reference dates and load lineage data and model output -------

for (i in 1:(nrow(reference_data_df)-1)){
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
    
    p_hat <- read.csv(paste0('../data/output/multicountry_output/validation/p_hat_', THIS_REF_DATE,'.csv'))
    Y_tilde <-read.csv(paste0('../data/output/multicountry_output/validation/Y_tilde_', THIS_REF_DATE, '.csv'))
  
    # Map model output times to dates
    t <- seq(from = min(p_hat$t), to = max(p_hat$t), by = 1)
    collection_date <- seq(from = ymd(FIRST_DATE), to = ymd(FIRST_DATE) + days(length(t)-1), by = "days")
    dates <- tibble(collection_date = collection_date,
                    t = as.integer(t))
    
    clean_global_df <- p_hat %>% left_join(dates, by = "t") %>% left_join(countries, by = "country_id") %>% 
        left_join(lineages, by = "lineage_id") %>% 
        rename(
            p_hat_mean = mean,
            p_hat_med = median,
            p_hat_lb = lb,
            p_hat_ub = ub,
            p_hat_lower_25th = lower_25th,
            p_hat_upper_75th = upper_75th
        ) %>% select(-name, -std_dev) %>% 
        left_join(Y_tilde, by = c("lineage_id", "country_id", "t")) %>% 
        rename(
            Y_tilde_mean = mean,
            Y_tilde_med = median,
            Y_tilde_lb = lb,
            Y_tilde_ub = ub,
            Y_tilde_lower_25th = lower_25th,
            Y_tilde_upper_75th = upper_75th
        ) %>% 
        select(country, lineage, t, collection_date, p_hat_mean, p_hat_med, 
               p_hat_lb, p_hat_ub, p_hat_lower_25th, p_hat_upper_75th, 
               Y_tilde_mean, Y_tilde_med, Y_tilde_lb, Y_tilde_ub,
               Y_tilde_lower_25th, Y_tilde_upper_75th) %>% 
        left_join(df, by = c("country", "lineage", "collection_date", "t")) %>% 
        group_by(t, country) %>% 
        mutate(N = sum(n)) %>% 
        # Y_tilde = ifelse( N>0, Y_tilde, 0)) %>% 
        ungroup() %>% 
        mutate(obs_prev = n/N,
               last_collection_date = ymd(LAST_COLLECTION_DATE),
               last_submission_date = ymd(LAST_DATE_IN_DATASET),
               reference_date = ymd(THIS_REF_DATE)) %>% 
        group_by(country) %>% 
        mutate(
            multicountry_period =case_when( 
                    collection_date <= ymd(LAST_COLLECTION_DATE) ~ 'multicountry calibration',
                    collection_date > ymd(LAST_COLLECTION_DATE) ~ 'multicountry forecast')) %>% 
        ungroup() 
   
     clean_global_df$obs_prev[clean_global_df$N==0]<-NA
     


    if (USE_CASE == 'local'){
      write.csv(clean_global_df, paste0('../data/processed/validation_data/clean_global_df_',THIS_REF_DATE, '.csv'), row.names = F)
      
    }else{
        write.csv(clean_global_df, paste0('/mnt/data/processed/clean_global_df_',THIS_REF_DATE, '.csv'), row.names = F)
    }
} # end loop over all reference dates 
