
here::i_am('historical_validation_code/monkeypox_postprocessing.R')

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(lubridate)

# data --------------------------------------------------------------------

clade_map <- read_csv('data/mpxv/clade_map.csv',
                        col_types = 'ci')
country_map <- read_csv('data/mpxv/country_map.csv')

df <- read_table('data/output/mpxv/stansummary.txt',
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
              ))

# Generate figs -----------------------------------------------------------

df %>% 
  filter(str_detect(name, 'mu_hat')) %>% 
  separate(name, c('name', NA, 'clade_id')) %>% 
  left_join(clade_map %>% 
              mutate(clade_id = as.character(clade_id - 1))) %>% 
  ggplot()+
  geom_point(aes(clade_membership,  mean))+
  geom_linerange(aes(clade_membership, ymin = q5, ymax = q95))


df %>% 
  filter(str_detect(name, 'r_hat')) %>% 
  separate(name, c('name', NA, 'clade_id', "country_id")) %>% 
  left_join(clade_map %>% 
              mutate(clade_id = as.character(clade_id - 1))) %>% 
  left_join(country_map %>% 
              mutate(country_id = as.character(country_id))) %>% 
  ggplot()+
  geom_point(aes(country, mean))+
  geom_linerange(aes(country, ymin = q5, ymax = q95))+
  facet_wrap(~clade_membership)


df %>% 
  filter(str_detect(name, 'p_hat')) %>% 
  separate(name, c('name', NA, 'country_id', 't', 'clade_id')) %>% 
  left_join(clade_map %>% 
              mutate(clade_id = as.character(clade_id))) %>% 
  left_join(country_map %>% 
              mutate(country_id = as.character(country_id))) %>% 
  left_join(tibble(t = as.character(1:125),
                       date = seq.Date(from = ymd('2022-04-27'),
                                                  by = 'day',
                                                  length.out = 125))) %>% 
  mutate(t = as.integer(t)) %>% 
  ggplot()+
  geom_line(aes(date, mean, group = clade_membership, color = clade_membership))+
  geom_ribbon(aes(x = date, ymin = q5, ymax = q95, group = clade_membership, fill = clade_membership), alpha = 0.2)+
  facet_wrap(~country)


df %>% 
  filter(str_detect(name, 'Omega_hat')) %>% 
  separate(name, c('name', NA, 'clade_id_1', 'clade_id_2')) %>% 
  left_join(clade_map %>% 
              mutate(clade_id_1 = as.character(clade_id -1)) %>% 
              rename(clade_1 = clade_membership)) %>%
  left_join(clade_map %>% 
              mutate(clade_id_2 = as.character(clade_id -1)) %>% 
              rename(clade_2 = clade_membership) %>% 
              select(-clade_id)) %>%
  filter(clade_2 > clade_1) %>% 
  ggplot(aes(clade_1, clade_2, fill = mean))+
  geom_tile()+
  labs(x = element_blank(),
       y = element_blank())

df %>% 
  filter(str_detect(name, 'Omega_hat')) %>% 
  separate(name, c('name', NA, 'clade_id_1', 'clade_id_2')) %>% 
  left_join(clade_map %>% 
              mutate(clade_id_1 = as.character(clade_id -1)) %>% 
              rename(clade_1 = clade_membership)) %>%
  left_join(clade_map %>% 
              mutate(clade_id_2 = as.character(clade_id -1)) %>% 
              rename(clade_2 = clade_membership) %>% 
              select(-clade_id)) %>%
  filter(clade_1 != clade_2) %>% 
  mutate(clade_id_1 = as.integer(clade_id_1),
         clade_id_2 = as.integer(clade_id_2)) %>% 
  filter(clade_id_1 < clade_id_2) %>% 
  ggplot(aes(interaction(clade_1, clade_2), mean))+
  geom_point()+
  geom_linerange(aes(interaction(clade_1, clade_2), ymin = q5, ymax = q95))



df %>% 
  filter(str_detect(name, 'r_hat')) %>% 
  separate(name, c('name', NA, 'clade_id', "country_id")) %>% 
  left_join(clade_map %>% 
              mutate(clade_id = as.character(clade_id - 1))) %>% 
  left_join(country_map %>% 
              mutate(country_id = as.character(country_id))) %>% 
  mutate(ISOALPHA = countrycode(sourcevar = country,
                                origin = 'country.name',
                                destination = 'iso3c')) %>% 
  left_join(immunity) %>% 
  left_join(df %>% 
              filter(str_detect(name, 'mu_hat')) %>% 
              separate(name, c('name', NA, 'clade_id')) %>% 
              left_join(clade_map %>% 
                          mutate(clade_id = as.character(clade_id - 1))) %>% 
              rename(mu_hat = mean) %>% 
              select(mu_hat, clade_membership)) %>% 
  mutate(deviation = mean - mu_hat) %>% 
  ggplot(aes(mean_susceptibility, deviation))+
  geom_point()+
  facet_wrap(~clade_membership)+
  geom_smooth(method = 'lm')
  


immunity <- read_csv('https://raw.githubusercontent.com/bansallab/mpx_landscape/main/estimates/world-by-country.csv')


immunity

