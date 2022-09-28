# This file takes the summarized output from the AL1 Stan run and turns it into
# figure summarizing results

here::i_am("historical_validation_code/AL1_postprocessing.R")

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(jsonlite)
library(here)
library(sf)
library(geojsonsf)
library(patchwork)
library(wesanderson)

# Globals -----------------------------------------------------------------

low = "#ddecff"
mid = "#9cc8ff"
high = "#00295D"
na_color = "#8E98A4"
space = '           '


# Argentina ---------------------------------------------------------------

arg <- geojson_sf('data/geojson/argentina.geojson') %>% 
  st_simplify(dTolerance = 1000) %>% 
  st_transform(crs = 4326) %>% 
  rename(AL1 = name) %>% 
  mutate(ADM0 = 'Argentina',
         ADM1 = str_to_lower(AL1)) %>% 
  select(ADM0, ADM1, geometry)

st_crs(arg) <- NA

# Paraguay ----------------------------------------------------------------

par <- st_read(dsn = "data/geojson/pry_admbnda_adm0_DGEEC_2020.shp",
               as_tibble = T,
               type = 6) %>% 
  st_simplify(dTolerance = 1000) %>% 
  mutate(ADM1 = 'paraguay') %>% 
  rename(ADM0 = ADM0_ES) %>% 
  select(ADM0, ADM1, geometry)

st_crs(par) <- NA

# Brazil ------------------------------------------------------------------

bra <- st_read(dsn = "data/geojson/bra_admbnda_adm1_ibge_2020.shp",
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
  filter(lin_id == 6, t %in% c(152, 166, 189)) %>% 
  mutate(x = if_else(t == min(t), 1, 
                     if_else(t == max(t), 
                                      3, 2)))


# Set up map --------------------------------------------------------------


al1 <- rbind(arg, bra) %>% 
  rbind(par) %>% 
  mutate(AL1 = str_to_lower(ADM1),
         AL1 = stringi::stri_trans_general(AL1, "Latin-ASCII")) %>% 
  left_join(x %>% 
              filter(lin_id == 6) %>% 
              select(AL1, Mean, t, x)) %>% 
  mutate(.mean = as.double(Mean))

cols <- c(colorRampPalette(c("#e7f0fa", "#c9e2f6", "#95cbee", "#0099dc", "#4ab04a", "#ffd73e"))(10),
                              colorRampPalette(c("#eec73a", "#e29421", "#e29421", "#f05336","#ce472e"), bias=2)(90))

p1 <- al1 %>% 
  filter(t == 152 | is.na(x)) %>% 
  ggplot(aes(geometry = geometry, fill = .mean))+
  geom_sf(show.legend = F)+
  scale_fill_gradient2(low = low,
                       mid = mid,
                       high =  high,
                       limits = c(0,1), na.value = na_color)+
  theme_void()+
  labs(tag = 'A',
       title  = paste0(space, '2022-06-01'))

p2 <- al1 %>% 
  filter(t == 166 | is.na(x)) %>% 
  ggplot(aes(geometry = geometry, fill = .mean))+
  geom_sf(show.legend = F)+
  scale_fill_gradient2(low = low,
                       mid = mid,
                      high =  high,
                      limits = c(0,1), na.value = na_color)+
  theme_void()+
  labs(title = paste0(space, '2022-06-16'))

p3 <- al1 %>% 
  filter(t == 189 | is.na(t)) %>% 
  #mutate(.mean = if_else(.mean > .75, .75, .mean)) %>% 
  ggplot(aes(geometry = geometry, fill = .mean))+
  geom_sf(show.legend = F)+
  scale_fill_gradient2(low = low,
                      mid = mid,
                      high =  high,
                      limits = c(0,1), na.value = na_color)+
  theme_void()+
  labs(title = paste0(space, '2022-07-01'))

p4 <- df %>% 
  filter(str_detect(name, "p_hat")) %>% 
  separate(name, c('name', NA, 'AL1', 't', 'lin_id')) %>% 
  mutate(t = as.integer(t),
         lin_id = as.integer(lin_id),
         AL1 = as.character(AL1)) %>% 
  janitor::clean_names() %>% 
  rename(AL1 = al1) %>% 
  mutate(mean = as.numeric(mean)) %>% 
  filter(lin_id == 6) %>% 
  left_join(tibble(t = 1:189,
                   date =  seq.Date(from = ymd('2022-01-01'), by  = 'day', length.out = 189))) %>% 
  left_join(lineage_map) %>%
  filter(t > 99) %>% 
  ggplot(aes(date, mean, group = AL1, color = mean))+
  geom_vline(aes(xintercept = ymd('2022-06-01')), alpha = .25)+
  geom_vline(aes(xintercept = ymd('2022-06-16')), alpha = .25)+
  geom_vline(aes(xintercept = ymd('2022-07-01')), alpha = .25)+
  geom_line()+
  theme_bw()+
  scale_color_gradient2(low = low,
                        mid= mid,
                       high =  high,
                       limits = c(0,1), na.value = "#de2d26")+
  labs(tag = "B",
       color = 'Estimated proportion of BA.5',
       y = 'Estimated proportion of BA.5',
       x = element_blank(),
       color = 'Estimated proportion of BA.5')+
  theme(legend.position='bottom')


p5 <- df %>% 
  filter(str_detect(name, "p_hat")) %>% 
  separate(name, c('name', NA, 'AL1', 't', 'lin_id')) %>% 
  mutate(t = as.integer(t),
         lin_id = as.integer(lin_id),
         AL1_id = as.integer(AL1)) %>% 
  select(-AL1) %>% 
  janitor::clean_names() %>% 
  mutate(mean = as.numeric(mean)) %>% 
  filter(t > 100) %>% 
  left_join(tibble(t = 1:169,
                   date =  seq.Date(from = ymd('2022-01-01'), by  = 'day', length.out = 169))) %>% 
  left_join(lineage_map) %>%
  rename(AL1_id = al1_id) %>% 
  left_join(AL1_map) %>% 
  left_join(clean %>% 
              select(country, AL1) %>% 
              unique()) %>% 
  filter(AL1 %in% c('paraguay', 'caba', 'rio de janeiro')) %>% 
  ggplot(aes(x = t, y = mean)) +
  stat_horizon(aes(x = date, y = mean), 
               show.legend = F, 
               reverse = F,
               horizonscale = 6, 
               origin = 'min')+
  scale_fill_hcl(palette = "peach", reverse = T)+
  theme_few()+
  facet_grid(lineage~AL1)+
  theme(
    panel.spacing.y=unit(0, "lines"),
    strip.text.y = element_text(angle = 0),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.border = element_blank(),
    axis.text.x =  element_text(angle = -45, hjust=-.5, vjust = 2,size = 0.7*axis_text_size)
  )+
  labs(x = element_blank(),
       tag = "C")
  

design <- "
  12355
  44455
"

fig <- p1 + p2 + p3 + p4 + p5 + plot_layout(design = design)

ggsave('data/output/figures/figure_5.png', 
       plot = fig,
       width = 15,
       height = 8)

# Data --------------------------------------------------------------------

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



supp_fig <- df %>% 
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
  arrange(country) %>% 
  mutate(country = ordered(country)) %>% 
  filter(lin_id == 6) %>% 
  mutate(AL1 = factor(AL1, levels = unique(prep$AL1), ordered = T)) %>% 
  ggplot(aes(x = AL1, y = exp(.mean * 7) - 1, color = country))+
  geom_point()+
  geom_linerange(aes(x = AL1, ymin = exp(.lower * 7) - 1, ymax = exp(.upper * 7) - 1, color = country))+
  facet_wrap(~lineage)+
  #geom_hline(aes(yintercept = 0), alpha = .5)+
  coord_flip()+
  theme_bw()+
  labs(tag = 'C',
       y = 'Posterior of local BA.5 fitness advantage relative to BA.2 (median with 95% Cl)',
       x = element_blank(),
       color = element_blank())+
  theme(legend.position = 'bottom')



