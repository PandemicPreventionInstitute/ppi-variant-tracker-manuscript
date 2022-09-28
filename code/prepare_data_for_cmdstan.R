# This script formats the data for cmdstan to be run in the multicountry model


USE_CASE <- 'local'
# Libraries ---------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(tidybayes)
library(cmdstanr)

# Constants ---------------------------------------------------------------

# Path to the processed GISAID metadata file produced by pre_processing.R
DATA_PATH <- '../data/processed/lineage_t_2022-07-01.csv'

# Path to the Stan model code
MODEL_PATH <- 'stancode/multivariate_variant_multinomial_ncp.stan'

# If the total number of times this lineage is observed is below this number,
# collapse them to `unassigned`.
OBSERVATION_THRESHOLD <- 50 # This is a global threshold

# Load data ---------------------------------------------------------------

df.orig <- read_csv(DATA_PATH)

# Run model starting this many days in the past
DURATION <- length(unique(df.orig$t))

# Today's date as a string
# this will be the "day run" (1 day above max collection date, assuming 1 
# country collected and submitted yesterday)

FIRST_DATE<-as.character(min(df.orig$collection_date))
LAST_COLLECTION_DATE <- as.character(max(df.orig$collection_date[df.orig$tot_seq>0]))
LAST_DATE_IN_DATASET <- as.character(max(df.orig$collection_date))

TODAY_DATE<-"2022-07-01"

# Nowcast up until the last submission date of sequence data
NUM_DAYS_NOWCAST <- as.integer(ymd("2022-07-19") - ymd(LAST_COLLECTION_DATE))

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
    #filter(collection_date %in% last_twenty_timepoints) %>% 
    group_by(lineage) %>% 
    summarize(n = sum(n)) %>% 
    filter(n == max(n)) %>% 
    pull(lineage)

# in case there is more than one most prevalent lineage
most_prevalent_lineage <- most_prevalent_lineage[[1]]
print(most_prevalent_lineage)

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

cmdstanr::write_stan_json(data, '../data/processed/data_for_cmdstan_2022_07_01.json')



# Run the cmdstan model using bash-----------------------------------------

setwd("~/Documents/variant-tracker/ppi-variant-tracker")
system("bash run_cmdstan_paper.sh")
