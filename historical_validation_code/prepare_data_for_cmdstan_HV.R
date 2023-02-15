# Author: Kaitlyn Johnson
# Modified from Zack Susswein

# Date initiated: 06-22-2022

# This script loads in the lineage_t datasets from different reference dates and 
# prepares them for running in cmdstan by generating a json data object

rm(list = ls())

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(tidybayes)
library(cmdstanr)

# Constants ---------------------------------------------------------------
if (USE_CASE == 'local'){
    REFERENCE_DATA_PATH <- '../data/processed/validation_data/reference_data_used.csv'
    setwd("~/Documents/ppi-variant-tracker-manuscript/code") # May need to be adjusted for local path 
}

# Get a subset of countries that have enough sequencing in the time period for comparing 
MLE_COUNTRIES1 <- read_csv('../data/processed/validation_data/lineage_t_2022-04-30.csv') %>% 
    group_by(country) %>% 
    summarise(tot_seq = sum(tot_seq)) %>% 
    filter(tot_seq >=3000) %>% pull(country)
MLE_COUNTRIES2 <- read_csv('../data/processed/validation_data/lineage_t_2022-07-01.csv') %>% 
    group_by(country) %>% 
    summarise(tot_seq = sum(tot_seq)) %>% 
    filter(tot_seq >=3000) %>% pull(country)
MLE_COUNTRIES_FLAG<- MLE_COUNTRIES2 %in% MLE_COUNTRIES1
MLE_COUNTRIES<-MLE_COUNTRIES2[MLE_COUNTRIES_FLAG]


# Load metadata on reference data ---------------------------------------------
reference_data_df <- read_csv(REFERENCE_DATA_PATH)
reference_data_df$most_prevalent_lineage <- rep(NA, nrow(reference_data_df)) # pre-allocate


# Get the last reference data number of sequences for observed predictions------ 

# This has no cut off based on duration, within the loop we will set the time for
# each dataset we are comparing
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







# Loop through all reference datasets and generate time-stamped json object------ 


for (i in 1:nrow(reference_data_df)){
    THIS_REF_DATE<- as.character(reference_data_df$reference_date[i])
    if (USE_CASE == 'local'){
        DATA_PATH <- paste0('../data/processed/validation_data/lineage_t_', THIS_REF_DATE, '.csv')
    }
    df.orig <- read_csv(DATA_PATH)

    # Run model starting this many days in the past
    DURATION <- length(unique(df.orig$t))
    print(paste0('Data from ', THIS_REF_DATE, ' was run for ', DURATION, ' days'))
    
    
    FIRST_DATE<-as.character(min(df.orig$collection_date))
    LAST_COLLECTION_DATE <- as.character(max(df.orig$collection_date[df.orig$tot_seq>0]))
    LAST_DATE_IN_DATASET <- as.character(max(df.orig$collection_date))
    
    
    # Set each model to run 21 days in the future
    NUM_DAYS_NOWCAST <- 21
    
    # Wrangle data ------------------------------------------------------------
    
    df <- df.orig %>% 
        filter(collection_date <= ymd(LAST_COLLECTION_DATE)) %>% 
        select(lineage, country, collection_date, n_seq, t) %>% 
        group_by(lineage, country, collection_date, t) %>% 
        summarize(n = sum(n_seq)) %>% 
        ungroup() %>% 
        droplevels() %>% 
        mutate(t = t - min(t) + 1)
    
    # Generate constants for data munging -------------------------------------
    
    # This assumes that the previous twenty days are observed in the dataset
    last_twenty_timepoints <- df %>% 
        arrange(collection_date) %>%
        pull(collection_date) %>%
        unique() %>% 
        tail(n=20)
    
    most_prevalent_lineage <- df %>% 
        filter(collection_date %in% last_twenty_timepoints) %>% 
        group_by(lineage) %>% 
        summarize(n = sum(n)) %>% 
        filter(n == max(n)) %>% 
        pull(lineage)
    
    # in case there is more than one most prevalent lineage
    most_prevalent_lineage<- most_prevalent_lineage[[1]]
    print(most_prevalent_lineage)
    
    reference_data_df$most_prevalent_lineage[i]<-most_prevalent_lineage
    
    
    # Preprocess data for Stan ------------------------------------------------
    
    # Number of timepoints observed
    N <- length(unique(df$t))
    
    # Number of categories the multinomial could produce
    ncat <- length(unique(df$lineage))
    
    # number of countries
    K <- length(unique(df$country))
    
    # Construct a tibble mapping numbers to lineage names
    # Make the most prevalent lineage have id = 1, so that it can be easily moved
    # to the top of of the df when sorted in the next step
    lineage_map <- tibble(lineage = df %>% 
                              filter(lineage != most_prevalent_lineage) %>%
                              pull(lineage) %>%
                              unique(),
                          lin_id = 2:ncat) %>% 
        bind_rows(tibble(lineage = most_prevalent_lineage,
                         lin_id = 1))
    
    # Construct a tibble mapping numbers to country names
    # This ensures that countries are consistently sorted into the same order
    country_map <- tibble(country = df %>% 
                              pull(country) %>%
                              unique(),
                          country_id = 1:K)
    
    # An N x (ncat x lineage) matrix of observed draws from each category at each time point
    Y <- df %>% 
        full_join(lineage_map) %>%
        arrange(lin_id) %>% 
        mutate(n = as.integer(n)) %>% 
        group_by(lineage, country) %>% 
        arrange(collection_date) %>% 
        ungroup() %>% 
        select(-collection_date) %>% 
        pivot_wider(id_cols = c(lineage, country),
                    names_from = t,
                    values_from = n) %>% 
        left_join(lineage_map) %>% 
        relocate(lin_id, .before = everything()) %>% 
        select(-lineage, -country) %>% 
        as.matrix()
    mode(Y) <- "integer"
    
    # Convert from N x (ncat x lineage) matrix to N x ncat x lineage array
    Y <- simplify2array(by(Y, Y[,1], as.matrix))[1:K, 2:(N+1), 1:ncat]
    
    # rows = countries, columns = timepoints, slices = lineage
    stopifnot(all(dim(Y) == c(K, N, ncat)))
    
    colnames(Y) <- NULL
    
    # validity testing
    stopifnot(nrow(Y) == K)
    stopifnot(ncol(Y) == N)
    stopifnot(typeof(Y) == "integer")
    
    # Total number of seqs observed at each timepoint
    # number of countries rows x number of timepoints columns
    trials <- df %>% 
        full_join(country_map) %>% 
        group_by(t, country) %>% 
        arrange(country_id) %>% 
        summarize(N = sum(n)) %>% 
        ungroup() %>% 
        arrange(t) %>% 
        mutate(N = as.integer(N)) %>% 
        pivot_wider(names_from = t,
                    values_from = N) %>% 
        select(-country) %>% 
        as.matrix()
    colnames(trials) <- NULL
    
    # Total number of sequences observed at each time point in the future dataset
    df_last_i <- df_last %>% filter(collection_date <= ymd(LAST_COLLECTION_DATE) + days(NUM_DAYS_NOWCAST),
                                    collection_date >=ymd(FIRST_DATE))
    
    # This is used to get the number of sequences we now know to have been collected 
    # on each of the days, to see what the model would've predicted their lineages as
    future_trials <- df %>% 
        full_join(df_last_i, by = c("country", "collection_date", "lineage")) %>% # add in the n_last
        left_join(country_map) %>% 
        filter(!is.na(country_id)) %>%  # get the right country id
        select(-t) %>% 
        mutate(t = t_last - min(t_last) + 1) %>% 
        group_by(t, country) %>% 
        arrange(country_id) %>% 
        summarise(N= sum(n_last)) %>% 
        ungroup() %>% 
        arrange(t) %>% 
        mutate(N= as.integer(N)) %>% 
        pivot_wider(names_from = t,
                    values_from = N) %>% 
        select(-country) %>% 
        as.matrix()
    colnames(future_trials) <- NULL

        
    # validity testing
    stopifnot(dim(trials) == c(K, N))
    stopifnot(typeof(trials) == "integer")
    stopifnot(typeof(future_trials) == "integer")
    #stopifnot(dim(future_trials) ==c(K, N+NUM_DAYS_NOWCAST))
    
    # should the likelihood be ignored?
    PRIOR_ONLY = 0
    
    N_future = ncol(future_trials)
    # Package data for Stan ---------------------------------------------------
    
    data <- list(N = N,
                 N_future = N_future,
                 K = K,
                 ncat = ncat,
                 Y = Y,
                 trials = trials,
                 prior_only = PRIOR_ONLY,
                 num_days_nowcast = NUM_DAYS_NOWCAST,
                 fut_trials = future_trials)

    # Write data for cmdstan --------------------------------------------------
    if (USE_CASE == 'local'){
        cmdstanr::write_stan_json(data, paste0('../data/processed/validation_data/data_for_cmdstan_', THIS_REF_DATE, '.json'))
    }
    
    # Package data for nnet for MLE single country model------------------------
    data_nnet <- df %>% 
        filter(country %in% MLE_COUNTRIES) %>% 
        full_join(lineage_map) %>%
        arrange(lin_id) %>% 
        mutate(n = as.integer(n)) %>% 
        group_by(lineage, country) %>% 
        arrange(collection_date) %>% 
        ungroup() %>% 
        select(-collection_date) %>% 
        pivot_wider(id_cols = c(lineage, country),
                    names_from = t,
                    values_from = n) %>% 
        left_join(lineage_map) %>% 
        relocate(lin_id, .before = everything()) %>% 
        select(-lin_id)
    
    write.csv(data_nnet, paste0('../data/processed/validation_data/data_for_nnet_', THIS_REF_DATE, '.csv'), row.names = F)
    
    
    
    
    
    
} # end loop through reference dates

# Keep track of the countries that you fit the MLE to 
MLE_countries<- data.frame(MLE_COUNTRIES)
write.csv(MLE_COUNTRIES,'../data/processed/validation_data/countries_fit_to_MLE.csv', row.names = F )
stopifnot(' Reference lineage is not the same for all datasets'= length(unique(reference_data_df$most_prevalent_lineage))==1)

