# Author: Kaitlyn Johnson
# Date initiated: 06-22-2022'

# This script runs the cmdstan model on each of the reference date datasets and
# generates the timestamped output files. It then uses stansummary to process the model
# output to get median, mean, and 95% CI for estimated variant prevalence 

rm(list = ls())
library(tidyverse)
library(lubridate)
library(tidybayes)
library(cmdstanr)
library(janitor)




REFERENCE_DATA_PATH <- 'data/processed/validation_data/reference_data_used.csv'
setwd("~/Documents/ppi-variant-tracker-manuscript") # might need to change
reference_data_df <- read_csv(REFERENCE_DATA_PATH)


# Fit each reference dataset to the model
for (i in 1:(nrow(reference_data_df)-1)){
    THIS_REF_DATE <- as.character(reference_data_df$reference_date[i])
    system(paste0("bash run_cmdstan_HV.sh ", THIS_REF_DATE))
}

# process the model output (note this is duplicative as we also
# process from the raw full distribution output in the process_results_from_cmdstan_HV.R)
for (i in 1:(nrow(reference_data_df)-1)){
    THIS_REF_DATE <- as.character(reference_data_df$reference_date[i])
    system(paste0("bash process_output_HV.sh ", THIS_REF_DATE))
    # Gets saved as `validation/processed_output_{REFERENCE_DATE}.csv`
}


# Get model estimated variant prevalence (p_hat) and sequence counts of each lienage (Y_tilde)
for (i in 1:(nrow(reference_data_df)-1)){
    THIS_REF_DATE <- as.character(reference_data_df$reference_date[i])
    processed_output<- read.csv(paste0('data/output/multicountry_output/validation/processed_output_', THIS_REF_DATE, '.csv')) %>% 
        clean_names() %>% 
        rename(median = x50,
               lower_25th = x25,
               upper_75th = x75,
               lb = x2_5,
               ub = x97_5) %>% 
        select(name, mean, std_dev, median, lb, ub, lower_25th, upper_75th)
    
    
    # estimated variant prevalence by lineage-country-timepoint
    p_hat <- processed_output %>% filter( grepl("p_hat", name, fixed = TRUE)) %>% 
        mutate(indices = substring(name, 7, (str_length(name)-1))) %>% 
        separate(indices,
                 into = c("country_id", "t", "lineage_id"),
                 sep = ",")
    
    # model predicted observed prevalence by lineage-country-timepoint
    Y_tilde <- processed_output %>% filter(grepl("Y_tilde", name, fixed = TRUE)) %>% 
        mutate(indices = substring(name,9, (str_length(name)-1))) %>% 
        separate(indices,
                 into = c("country_id", "t", "lineage_id"),
                 sep = ",")

    # Save with date appended for loading into processing script 
    write.csv(p_hat, paste0('data/output/multicountry_output/validation/p_hat_',THIS_REF_DATE, '.csv'), row.names = F)
    write.csv(Y_tilde, paste0('data/output/multicountry_output/validation/Y_tilde_', THIS_REF_DATE, '.csv'), row.names = F)
 
    
    
}

