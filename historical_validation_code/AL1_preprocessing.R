# This file takes GISAID metadata and munges it into shape to run the variant
# model for select countries at the AL1 level.

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(cmdstanr)

# Data --------------------------------------------------------------------

raw <- read_delim('data/2022-07-01_metadata.csv',
                 delim = ',')

# Globals -----------------------------------------------------------------

countries <- c('Brazil',
               'Argentina',
               'Paraguay')

lineages <- c('BA.1',
              'BA.2',
              'BA.2.12.1',
              'BA.4',
              'BA.5')

earliest_date <- ymd('2022-01-01')

NUM_DAYS_NOWCAST <- 20

# Subset data -------------------------------------------------------------

clean <- raw %>% 
  filter(collection_date >= earliest_date) %>% 
  # Pull out countries of interest
  mutate(country = str_split_fixed(location, '/', n = 3)[,2],
         country = str_trim(country)) %>% 
  filter(country %in% countries) %>% 
  # Pull out AL1, discarding any sub-AL1 info
  # and collapse lineages to buckets of interest
  mutate(AL1 = str_split_fixed(location, '/', n = 4)[,3],
         AL1 = str_trim(AL1),
         lineage = if_else(pango_lineage %in% lineages,
                           pango_lineage,
                           'other')) %>% 
  filter(!(AL1 == "" & country != 'Paraguay')) %>% 
  mutate(AL1 = if_else(country == 'Paraguay', 'Paraguay', AL1)) %>% 
  mutate(AL1 = str_to_lower(AL1)) %>% 
  mutate(AL1 = case_when(
    AL1 == "ciudad de buenos aires" ~ "caba",
    AL1 == "distrito federal" ~ "federal district",
    AL1 == "san juan nepomuceno" ~ "caazapa",
    AL1 == "cuidad del este" ~ " alto parana",
    AL1 == "fernando de la mora" ~ "central",
    AL1 == "san lorenzo" ~ "central",
    AL1 == "confluencia" ~ "neuquen",
    AL1 == "altos" ~ "cordillera",
    AL1 == "lambare" ~ "central",
    AL1 == "luque" ~ "central",
    AL1 == "karapa`i" ~ "amambay",
    AL1 == 'capiata' ~ "central",
    AL1 == "itacurubi del rosario" ~ "san pedro",
    AL1 == "villeta" ~ "central",
    AL1 == "nemby" ~ "central",
    AL1 == "itagua" ~ "central",
    AL1 == "limpio" ~ "central",
    AL1 == "villarrica" ~ "guaira",
    AL1 == "mariano roque alonso" ~ "central",
    AL1 == "ita" ~ "central",
    AL1 == "coronel oviedo" ~ "caaguazu",
    AL1 == "villa elisa" ~ "central",
    AL1 == "quiindy" ~ "paraguarí",
    AL1 == "curuguaty" ~ "canindeyu",
    AL1 == "ypane" ~ "central",
    AL1 == "aregua" ~ "central",
    AL1 == "piribebuy" ~ "cordillera",
    AL1 == "santa rosa" ~ "misiones",
    AL1 == "eusebio ayala" ~ "cordillera",
    AL1 == "san berdnardino" ~ "cordillera",
    AL1 == "3 de mayo" ~ "caazapa",
    AL1 == "caacupe" ~ "cordillera",
    AL1 == "caiibary" ~ "san pedro",
    AL1 == "san juan nepomuceno" ~ "caazapa",
    AL1 == "pedro juan caballero" ~ "amambay",
    AL1 == "avellaneda"   ~ "buenos aires",
    AL1 == "bahia blanca" ~ "buenos aires",
    AL1 == "ciudad autonoma de buenos aires" ~ "caba",
    AL1 == 'la matanza' ~ "buenos aires",
    AL1 == 'la plata' ~ "buenos aires",
    AL1 == 'lanus' ~ "buenos aires",
    AL1 == 'lomas de zamora' ~ "buenos aires",
    AL1 == 'longchamps' ~ "buenos aires",
    AL1 == "mar del plata" ~ "buenos aires",
    AL1 == "mercedes" ~ "buenos aires",
    AL1 == "merlo" ~ "buenos aires",
    AL1 == "montegrande" ~ "buenos aires",
    AL1 == "olivos" ~ "buenos aires",
    AL1 == "presidente peron" ~ "buenos aires",
    AL1 == "quilmes" ~ "buenos aires",
    AL1 == "ramos mejia"  ~ "buenos aires",
    AL1 == "rio tercero" ~ "cordoba",
    AL1 == "san antonio de areco" ~ "buenos aires",
    AL1 == "san isidro" ~ "buenos aires",
    AL1 == "san miguel" ~ "buenos aires",
    AL1 == "santa teresita" ~ "buenos aires",
    AL1 == "suipacha" ~ "buenos aires",
    AL1 == "tandil"  ~ "buenos aires",
    AL1 == "villa fiorito"  ~ "buenos aires",
    TRUE ~ as.character(AL1)
  )) %>% 
  group_by(collection_date, country, AL1, lineage) %>% 
  # Reduce to counts
  summarize(n = n()) %>% 
  ungroup() %>% 
  rename(date = collection_date) %>% 
  droplevels()


t_map <- clean %>%
  select(date) %>% 
  unique() %>% 
  mutate(t = row_number())

df <- clean %>% 
  select(date, AL1, lineage) %>% 
  unique() %>% 
  complete(date, AL1, lineage) %>% 
  # add back in countries
  left_join(clean %>% 
              select(AL1, country) %>% 
              unique()) %>% 
  left_join(clean) %>% 
  # make day count t, assumes that all days are within the same year
  mutate(t = yday(date), 
         t = t - min(t) + 1) %>%
  mutate(n = if_else(!is.na(n), n, as.integer(0)))

# Generate constants for data munging -------------------------------------

# This assumes that the previous twenty days are observed in the dataset
last_twenty_timepoints <- df %>% 
  arrange(date) %>%
  pull(date) %>%
  unique() %>% 
  tail(n=20)

most_prevalent_linage <- df %>% 
  #filter(date %in% last_twenty_timepoints) %>% 
  group_by(lineage) %>% 
  summarize(n = sum(n)) %>% 
  filter(n == max(n)) %>% 
  pull(lineage)

# in case there is more than one most prevalent lineage
most_prevalent_linage <- most_prevalent_linage[[1]]
print(most_prevalent_linage)

# Preprocess data for Stan ------------------------------------------------

# Number of timepoints observed
N <- length(unique(df$t))

# Number of categories the multinomial could produce
ncat <- length(unique(df$lineage))

# number of countries
K <- length(unique(df$AL1))

# Construct a tibble mapping numbers to lineage names
# Make the most prevalent lineage have id = 1, so that it can be easily moved
# to the top of of the df when sorted in the next step
lineage_map <- tibble(lineage = df %>% 
                        filter(lineage != most_prevalent_linage) %>%
                        pull(lineage) %>%
                        unique(),
                      lin_id = 2:ncat) %>% 
  bind_rows(tibble(lineage = most_prevalent_linage,
                   lin_id = 1))

# Construct a tibble mapping numbers to AL1 names
# This ensures that countries are consistently sorted into the same order
AL1_map <- tibble(AL1 = df %>% 
                        pull(AL1) %>%
                        unique(),
                      AL1_id = 1:K)

# An N x (ncat x lineage) matrix of observed draws from each category at each time point
Y <- df %>% 
  full_join(lineage_map) %>%
  arrange(lin_id) %>% 
  mutate(n = as.integer(n)) %>% 
  group_by(lineage, AL1) %>% 
  arrange(date) %>% 
  ungroup() %>% 
  select(-date) %>% 
  pivot_wider(id_cols = c(lineage, AL1),
              names_from = t,
              values_from = n) %>% 
  left_join(lineage_map) %>% 
  relocate(lin_id, .before = everything()) %>% 
  select(-lineage, -AL1) %>% 
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
stopifnot(!is.na(sum(Y)))

# Total number of seqs observed at each timepoint
trials <- df %>% 
  full_join(AL1_map) %>% 
  group_by(t, AL1) %>% 
  arrange(AL1_id) %>% 
  summarize(N = sum(n)) %>% 
  ungroup() %>% 
  arrange(t) %>% 
  mutate(N = as.integer(N)) %>% 
  pivot_wider(names_from = t,
              values_from = N) %>% 
  select(-AL1) %>% 
  as.matrix()
colnames(trials) <- NULL

# validity testing
stopifnot(dim(trials) == c(K, N))
stopifnot(typeof(trials) == "integer")
stopifnot(!is.na(sum(trials)))

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

cmdstanr::write_stan_json(data, 'data/processed/AL1_data_for_cmdstan.json')


