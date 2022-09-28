
here::i_am('historical_validation_code/preprocess_flu_data.R')

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(readxl)
library(here)
library(countrycode)
library(rgeos)
library(rworldmap)

# Load data ---------------------------------------------------------------

raw <- read_xls('data/flu/gisaid_epiflu_isolates.xls')
raw2 <- read_xls('data/flu/gisaid_epiflu_isolates (1).xls')


# Function to do all the actual work on each hemisphere -------------------


preprocess_data <- function(df, hemisphere){
  
  df <- df %>% 
    full_join(tibble(
      t = 1:max(df$t),
      collection_date = seq.Date(from = ymd('2021-7-1'),
                                 by = 'day',
                                 length.out = 366)
    )) %>% 
    complete(t, country, subtype, fill = list(n = 0)) %>% 
    filter(!is.na(subtype), !is.na(country), !is.na(t), !is.na(n), country != '') %>% 
    select(-collection_date) %>% 
    left_join(tibble(
      t = 1:max(df$t),
      collection_date = seq.Date(from = ymd('2021-7-1'),
                                 by = 'day',
                                 length.out = 366)))
  
  # Today's date as a string
  # this will be the "day run" (1 day above max collection date, assuming 1 
  # country collected and submitted yesterday)
  
  FIRST_DATE<-as.character(min(df$collection_date))
  LAST_COLLECTION_DATE <- as.character(max(df$collection_date))
  LAST_DATE_IN_DATASET <- as.character(max(df$collection_date))
  
  TODAY_DATE<-"2022-07-01"
  
  # Nowcast up until the last submission date of sequence data
  NUM_DAYS_NOWCAST <- as.integer(ymd("2022-07-19") - ymd(LAST_COLLECTION_DATE))
  
  # Generate constants for data munging -------------------------------------
  
  most_prevalent_subtype <- df %>% 
    #filter(collection_date %in% last_twenty_timepoints) %>% 
    group_by(subtype) %>% 
    summarize(n = sum(n)) %>% 
    filter(n == max(n)) %>% 
    pull(subtype)
  
  # in case there is more than one most prevalent lineage
  most_prevalent_subtype <- most_prevalent_subtype[[1]]
  print(most_prevalent_subtype)
  
  # Preprocess data for Stan ------------------------------------------------
  
  # Number of timepoints observed
  N <- length(unique(df$t))
  
  # Number of categories the multinomial could produce
  ncat <- length(unique(df$subtype))
  
  # number of countries
  K <- length(unique(df$country))
  
  # Construct a tibble mapping numbers to lineage names
  # Make the most prevalent lineage have id = 1, so that it can be easily moved
  # to the top of of the df when sorted in the next step
  subtype_map <- tibble(subtype = df %>% 
                          filter(subtype != most_prevalent_subtype) %>%
                          pull(subtype) %>%
                          unique(),
                        strain_id = 2:ncat) %>% 
    bind_rows(tibble(subtype = most_prevalent_subtype,
                     strain_id = 1))
  
  # Construct a tibble mapping numbers to country names
  # This ensures that countries are consistently sorted into the same order
  country_map <- tibble(country = df %>% 
                          pull(country) %>%
                          unique(),
                        country_id = 1:K)
  
  # An N x (ncat x lineage) matrix of observed draws from each category at each time point
  Y <- df %>% 
    full_join(subtype_map) %>%
    arrange(strain_id) %>% 
    mutate(n = as.integer(n)) %>% 
    group_by(subtype, country) %>% 
    arrange(collection_date) %>% 
    ungroup() %>% 
    select(-collection_date) %>% 
    pivot_wider(id_cols = c(subtype, country),
                names_from = t,
                values_from = n) %>% 
    left_join(subtype_map) %>% 
    relocate(strain_id, .before = everything()) %>% 
    select(-subtype, -country) %>% 
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
  
  # validity testing
  stopifnot(dim(trials) == c(K, N))
  stopifnot(typeof(trials) == "integer")
  
  # should the likelihood be ignored?
  PRIOR_ONLY = 0
  
  # Package data for Stan ---------------------------------------------------
  
  data <- list(N = N,
               K = K,
               ncat = ncat,
               Y = Y,
               trials = trials,
               prior_only = PRIOR_ONLY,
               num_days_nowcast = NUM_DAYS_NOWCAST)
  
  # Compile model -----------------------------------------------------------
  
  cmdstanr::write_stan_json(data, paste0('data/processed/data_for_cmdstan_flu_', hemisphere, '.json'))
  
  write_csv(subtype_map, paste0('data/flu/subtype_map_', hemisphere, '.csv'))
  write_csv(country_map, paste0('data/flu/country_map_', hemisphere, '.csv'))
  
}


# Get centroids -----------------------------------------------------------

# get world map
wmap <- getMap(resolution="low")

# get centroids
centroids <- gCentroid(wmap, byid=TRUE)

# get a data.frame with centroids
hemisphere <- as_tibble(centroids, rownames = 'country.name') %>% 
  mutate(iso3c = countrycode(sourcevar = country.name,
                             origin = 'country.name',
                             destination = 'iso3c'),
         iso3c = if_else(country.name == 'Kosovo', 
                         'XKX',
                         iso3c)) %>% 
  mutate(hemisphere = if_else(y > 0, 'Northern', 'Southern')) %>% 
  select(iso3c, hemisphere)

# Load the data and pslit by hemisphere -----------------------------------

df <- raw %>% 
  janitor::clean_names() %>% 
  select(location, subtype, collection_date) %>% 
  mutate(country = str_split_fixed(location, '/', 3)[,2],
         country = trimws(country),
         collection_date = ymd(collection_date)) %>% 
  group_by(subtype, collection_date, country) %>% 
  summarize(n = n()) %>% 
  ungroup()

df2 <- raw2 %>% 
  janitor::clean_names() %>% 
  select(location, subtype, collection_date) %>% 
  mutate(country = str_split_fixed(location, '/', 3)[,2],
         country = trimws(country),
         collection_date = ymd(collection_date)) %>% 
  group_by(subtype, collection_date, country) %>% 
  summarize(n = n()) %>% 
  ungroup()

df <- rbind(df, df2) %>% 
  filter(!is.na(country), !is.na(subtype), !is.na(collection_date), country != '') %>% 
  mutate(t = as.integer(collection_date - min(collection_date)) + 1,
         iso3c = countrycode(sourcevar = country,
                             origin = 'country.name',
                             destination = 'iso3c'),
         iso3c = if_else(country == 'Kosovo', 
                         'XKX',
                         iso3c)) %>% 
  left_join(hemisphere)


# Save combined flu metadata ----------------------------------------------

write_csv(df, 'data/flu/combined.csv')

# Run the preprocessing ---------------------------------------------------


northern <- df %>% filter(hemisphere == 'Northern')
southern <- df %>% filter(hemisphere == 'Southern')

preprocess_data(northern, "northern")
preprocess_data(southern, "southern")









         