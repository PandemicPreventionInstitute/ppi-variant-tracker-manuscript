

library(tidyverse)
library(sf)
library(data.table)
library(geojsonsf)
library(here)
library(patchwork)

# Argentina ---------------------------------------------------------------

arg <- geojson_sf('../data/geojson/argentina.geojson') %>% 
  st_simplify(dTolerance = 1000) %>% 
  st_transform(crs = 4326)

st_crs(arg) <- NA

# Paraguay ----------------------------------------------------------------

par <- st_read(dsn = "~/Downloads/pry_adm_dgeec_2020_shp/pry_admbnda_adm0_DGEEC_2020.shp",
               as_tibble = T,
               type = 6) %>% 
  st_simplify(dTolerance = 1000) %>% 
  mutate(ADM1 = 'paraguay') %>% 
  rename(ADM0 = ADM0_ES) %>% 
  select(ADM0, ADM1, geometry)

st_crs(par) <- NA

# Brazil ------------------------------------------------------------------

bra <- st_read(dsn = "~/Downloads/bra_adm_ibge_2020_shp/bra_admbnda_adm1_ibge_2020.shp",
               as_tibble = TRUE,
               type = 6) %>% 
  #st_simplify(dTolerance = 1000) %>% 
  rename(ADM0 = ADM0_PT,
         ADM1 = ADM1_PT) %>% 
  select(ADM0, ADM1, geometry)

# Modeled output ----------------------------------------------------------

df <- read_csv('data/processed/AL1_summary.csv')

prep <- df %>% 
  filter(str_detect(name, 'r_hat')) %>% 
  rename(.mean = Mean,
         .lower = `5%`,
         .median = `50%`,
         .upper = `95%`) %>% 
  filter(!str_detect(name, 'L_Omega')) %>% 
  separate(name, c('name', NA, 'lin_id', 'AL1_id')) %>%
  mutate(lin_id = as.numeric(lin_id) +1,
         .mean = as.numeric(.mean),
         AL1_id = as.integer(AL1_id)) %>% 
  left_join(lineage_map) %>% 
  left_join(AL1_map) %>% 
  left_join(clean %>% 
              select(country, AL1) %>% 
              unique()) %>% 
  arrange(country) 


x <- df %>% 
  filter(str_detect(name, "p_hat")) %>% 
  separate(name, c('name', NA, 'AL1', 't', 'lin_id')) %>% 
  mutate(t = as.integer(t),
         lin_id = as.integer(lin_id),
         AL1_id = as.integer(AL1) ) %>% 
  select(-AL1) %>% 
  left_join(AL1_map) %>% 
  filter(lin_id == 6, t %in% c(120, 150, 180))


df %>% 
  filter(str_detect(name, "p_hat")) %>% 
  separate(name, c('name', NA, 'AL1', 't', 'lin_id')) %>% 
  mutate(t = as.integer(t),
         lin_id = as.integer(lin_id),
         AL1 = as.character(AL1)) %>% 
  janitor::clean_names() %>% 
  rename(AL1 = al1) %>% 
  mutate(mean = as.numeric(mean)) %>% 
  filter(lin_id == 6) %>% 
  ggplot(aes(t, mean, group = AL1))+
  geom_line()


x %>% 
  filter(t == 120, lin_id == 6) %>% 
  select(geometry, AL1, t, mean) %>% 
  write_csv('data/geojson/BA5_first_timepoint.csv')

x %>% 
  filter(t == 150, lin_id == 6) %>% 
  select(geometry, AL1, t, mean) %>% 
  write_csv('data/geojson/BA5_middle_timepoint.csv')

x %>% 
  filter(t == 180, lin_id == 6) %>% 
  select(geometry, AL1, t, mean) %>% 
  write_csv('data/geojson/BA5_late_timepoint.csv')


# Set up map --------------------------------------------------------------

al1 <- rbind(arg, bra) %>% 
  rbind(par) %>% 
  mutate(AL1 = str_to_lower(ADM1),
         AL1 = stringi::stri_trans_general(AL1, "Latin-ASCII")) %>% 
  left_join(x %>% 
              filter(lin_id == 6) %>% 
              select(AL1, Mean, t)) %>% 
  mutate(.mean = as.double(Mean))

prep

p1 <- al1 %>% 
  filter(t == 120) %>% 
  ggplot(aes(geometry = geometry, fill = .mean))+
  geom_sf()+
  scale_fill_gradientn(limits = c(0,1),
                       colours=c("navyblue", "darkmagenta", "darkorange1"))

p2 <- al1 %>% 
  filter(t == 150) %>% 
  ggplot(aes(geometry = geometry, fill = .mean))+
  geom_sf()+
  scale_fill_gradientn(limits = c(0,1),
                       colours=c("navyblue", "darkmagenta", "darkorange1"))

p3 <- al1 %>% 
  filter(t == 180) %>% 
  ggplot(aes(geometry = geometry, fill = .mean))+
  geom_sf()+
  scale_fill_gradientn(limits = c(0,1),
                       colours=c("navyblue", "darkmagenta", "darkorange1"))

p1 + p2 + p3
