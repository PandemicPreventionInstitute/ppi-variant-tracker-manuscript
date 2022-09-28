
here::i_am('historical_validation_code/preprocess_monkeypox_data.R')


# Libraries ---------------------------------------------------------------

library(tidyverse)
library(lubridate)

# Read in data ------------------------------------------------------------

df.orig <- read_delim('data/mpxv/nextstrain_monkeypox_hmpxv1_metadata.tsv',
                      delim = '\t') %>% 
  mutate(date = str_split_fixed(date, ' ', 2)[,1],
         date = ymd(date)) %>% 
  filter(date >= ymd('2022-1-1')) %>% 
  group_by(date, clade_membership, country) %>% 
  summarize(n = n()) %>% 
  ungroup()

# Today's date as a string
# this will be the "day run" (1 day above max collection date, assuming 1 
# country collected and submitted yesterday)

FIRST_DATE<-as.character(min(df.orig$date))

# Nowcast up until the last submission date of sequence data
NUM_DAYS_NOWCAST <- 10

# Wrangle data ------------------------------------------------------------

df <- df.orig %>% 
  mutate(t = as.integer(date - min(date)) + 1) %>% 
  mutate(t = t - min(t) + 1)

# Generate constants for data munging -------------------------------------

most_prevalent_lineage <- df %>% 
  #filter(collection_date %in% last_twenty_timepoints) %>% 
  group_by(clade_membership) %>% 
  summarize(n = sum(n)) %>% 
  filter(n == max(n)) %>% 
  pull(clade_membership)

# in case there is more than one most prevalent lineage
most_prevalent_lineage <- most_prevalent_lineage[[1]]
print(most_prevalent_lineage)

# Preprocess data for Stan ------------------------------------------------

# Number of timepoints observed
N <- max(df$t)

# Number of categories the multinomial could produce
ncat <- length(unique(df$clade_membership))

# number of countries
K <- length(unique(df$country))

# Construct a tibble mapping numbers to lineage names
# Make the most prevalent lineage have id = 1, so that it can be easily moved
# to the top of of the df when sorted in the next step
clade_map <- tibble(clade_membership = df %>% 
                        filter(clade_membership != most_prevalent_lineage) %>%
                        pull(clade_membership) %>%
                        unique(),
                      clade_id = 2:ncat) %>% 
  bind_rows(tibble(clade_membership = most_prevalent_lineage,
                   clade_id = 1))

# Construct a tibble mapping numbers to country names
# This ensures that countries are consistently sorted into the same order
country_map <- tibble(country = df %>% 
                        pull(country) %>%
                        unique(),
                      country_id = 1:K)

# An N x (ncat x lineage) matrix of observed draws from each category at each time point
Y <- df %>% 
  full_join(tibble(t = 1:max(df$t),
                   country = NA)) %>% 
  complete(t, country, clade_membership, fill = list(n = 0)) %>% 
  full_join(clade_map) %>%
  arrange(clade_id) %>% 
  mutate(n = as.integer(n)) %>% 
  filter(!is.na(country), !is.na(clade_membership)) %>% 
  group_by(clade_membership, country) %>% 
  arrange(t) %>% 
  ungroup() %>% 
  select(-date) %>% 
  pivot_wider(id_cols = c(clade_membership, country),
              names_from = t,
              values_from = n,
              values_fill = 0) %>% 
  left_join(clade_map) %>% 
  relocate(clade_id, .before = everything()) %>% 
  select(-clade_membership, -country) %>% 
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
  full_join(tibble(t = 1:max(df$t),
                   country = NA)) %>% 
  complete(t, country, clade_membership, fill = list(n = 0)) %>% 
  full_join(clade_map) %>%
  arrange(clade_id) %>% 
  mutate(n = as.integer(n)) %>% 
  filter(!is.na(country), !is.na(clade_membership)) %>%  
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

cmdstanr::write_stan_json(data, 'data/processed/data_for_cmdstan_mpxv.json')

write_csv(clade_map, 'data/mpxv/clade_map.csv')
write_csv(country_map, 'data/mpxv/country_map.csv')

