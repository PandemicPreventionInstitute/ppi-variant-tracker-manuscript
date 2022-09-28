
here::i_am('historical_validation_code/postprocess_flu.R')

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(cmdstanr)

# Data --------------------------------------------------------------------

subtype_map_northern <- read_csv('data/flu/subtype_map_northern.csv',
                        col_types = 'ci')
country_map_northern <- read_csv('data/flu/country_map_northern.csv')

northern <- read_table('data/output/flu/northern/stansummary.txt',
                 skip = 15,
                 col_names =  c(
                   "name",
                   "mean",
                   "mcse",
                   "sd",
                   "q5",
                   "q50",
                   "q95",
                   "n_eff",
                   "n_eff/s",
                   "rhat"
                 )) %>% 
  mutate(type = "northern")

subtype_map_southern <- read_csv('data/flu/subtype_map_southern.csv',
                                 col_types = 'ci')
country_map_southern <- read_csv('data/flu/country_map_southern.csv')


southern <- read_table('data/output/flu/southern/stansummary.txt',
                       skip = 15,
                       col_names =  c(
                         "name",
                         "mean",
                         "mcse",
                         "sd",
                         "q5",
                         "q50",
                         "q95",
                         "n_eff",
                         "n_eff/s",
                         "rhat"
                       )) %>% 
  mutate(type = "southern")

pal_hemi <- dutchmasters_pal(palette = 'anatomy')(7)

hemi_colors = c('Northern Hemisphere' = pal_hemi[4],
                   'Southern Hemisphere' = pal_hemi[5])

pal_clade <- dutchmasters_pal(palette = "view_of_Delft")(12)

clade_colors


# Data for panel A --------------------------------------------------------

country_map_southern %>% 
  filter(country %in% c('Australia', 'Brazil', 'South Africa', 'Congo, the Democatic Republic of')) %>% 
  pull(country_id)


variables <- expand_grid(country = c(2, 4, 6, 14),
       t = c(1:366),
       clade = c(1:4)) %>% 
  mutate(variable = paste0("p_hat[", country, ',', t, ',', clade, "]"))

southern_draws_raw <- cmdstanr::read_cmdstan_csv(files = list.files('data/output/flu/southern', 
                                                                pattern = "*.csv",
                                                                full.names = T),
                                             variables = variables$variable,
                                             sampler_diagnostics = '',
                                             format = "draws_matrix")

southern_draws <- as_tibble(southern_draws_raw$post_warmup_draws) %>% 
  mutate(draw = row_number()) %>%
  filter(draw %in% sample(1:2000, 100)) %>% 
  pivot_longer(cols = -c(draw)) %>% 
  separate(name, c('name', NA, 'country_id', 't', 'strain_id')) %>% 
  mutate(value = as.double(value))  %>% 
<<<<<<< HEAD
  left_join(subtype_map_southern %>% 
=======
  left_join(subtype_map %>% 
>>>>>>> da466650e1c3330863efc602bbd44d54a3045301
              mutate(strain_id = as.character(strain_id))) %>% 
  left_join(country_map_southern %>% 
              mutate(country_id = as.character(country_id))) %>% 
  left_join(tibble(t = as.character(1:366),
                   date = seq.Date(from = ymd('2021-07-1'),
                                   by = 'day',
                                   length.out = 366))) %>% 
  mutate(t = as.integer(t)) %>% 
  filter(country  %in% c('Australia', 'Brazil', 'South Africa', 'Congo, the Democatic Republic of'))

write_csv(southern_draws, 'data/output/flu/southern/draws_for_fig.csv')


# Data for panel B --------------------------------------------------------

country_map_northern %>% 
  filter(country %in% c('United States', 'India', 'Ghana', 'Panama')) %>% 
  pull(country_id)


variables <- expand_grid(country = c(29, 35, 64, 89),
                         t = c(1:366),
                         clade = c(1:9)) %>% 
  mutate(variable = paste0("p_hat[", country, ',', t, ',', clade, "]"))

northern_draws_raw <- cmdstanr::read_cmdstan_csv(files = list.files('data/output/flu/northern', 
                                                                    pattern = "*.csv",
                                                                    full.names = T),
                                                 variables = variables$variable,
                                                 sampler_diagnostics = '',
                                                 format = "draws_matrix")

northern_draws <- as_tibble(northern_draws_raw$post_warmup_draws) %>% 
  mutate(draw = row_number()) %>%
  filter(draw %in% sample(1:2000, 100)) %>% 
  pivot_longer(cols = -c(draw)) %>% 
  separate(name, c('name', NA, 'country_id', 't', 'strain_id')) %>% 
<<<<<<< HEAD
  mutate(value = as.double(value)) %>% 
  left_join(subtype_map_northern %>% 
              mutate(strain_id = as.character(strain_id))) %>% 
  left_join(country_map_northern %>% 
              mutate(country_id = as.character(country_id))) %>% 
  left_join(tibble(t = as.character(1:366),
                   date = seq.Date(from = ymd('2021-07-1'),
                                   by = 'day',
                                   length.out = 366))) %>% 
  mutate(subtype = if_else(subtype %in% c("A / H1N1",
                                          "B", 
                                          "A / H3",
                                          "A / H3N2"),
                           subtype,
                           "Other")) %>% 
  group_by(draw, country, date, subtype) %>% 
  summarize(value = sum(value))
=======
  mutate(value = as.double(value))
>>>>>>> da466650e1c3330863efc602bbd44d54a3045301

write_csv(northern_draws, 'data/output/flu/northern/draws_for_fig.csv')


# Data for panel C --------------------------------------------------------

northern %>% 
  filter(str_detect(name, 'r_hat')) %>% 
  separate(name, c('name', NA, 'strain_id', "country_id")) %>% 
  left_join(subtype_map_northern %>% 
              mutate(strain_id = as.character(strain_id - 1))) %>% 
  left_join(country_map_northern %>% 
              mutate(country_id = as.character(country_id))) %>% 
  full_join(
    southern %>% 
      filter(str_detect(name, 'r_hat')) %>% 
      separate(name, c('name', NA, 'strain_id', "country_id")) %>% 
      left_join(subtype_map_southern %>% 
                  mutate(strain_id = as.character(strain_id - 1))) %>% 
      left_join(country_map_southern %>% 
                  mutate(country_id = as.character(country_id)))) %>% 
  mutate(type = if_else(type == 'northern',
                        'Northern Hemisphere',
                        'Southern Hemisphere'),
         mean = exp(mean *  7) - 1) %>% 
  write_csv('data/output/flu/r_hats.csv')


# Data for panel D --------------------------------------------------------

northern %>% 
  filter(str_detect(name, 'mu_hat')) %>% 
  separate(name, c('name', NA, 'strain_id')) %>% 
  left_join(subtype_map_northern %>% 
              mutate(strain_id = as.character(strain_id - 1))) %>% 
  full_join(southern %>% 
              filter(str_detect(name, 'mu_hat')) %>% 
              separate(name, c('name', NA, 'strain_id')) %>% 
              left_join(subtype_map_southern %>% 
                          mutate(strain_id = as.character(strain_id - 1)))) %>% 
  mutate(type = if_else(type == 'northern',
                        'Northern Hemisphere',
                        'Southern Hemisphere'),
         q5 = exp(q5 * 7) - 1,
         q95 = exp(q95 * 7) - 1,
         mean = exp(mean *  7) - 1) %>% 
  write_csv('data/output/flu/mu_hats.csv')
