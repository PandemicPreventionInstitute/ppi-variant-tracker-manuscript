# This file generates the proposed Figure 3 for the variant tracker paper

here::i_am('code/figures.R')

# Libraries ---------------------------------------------------------------
rm(list = ls())

library(tidyverse)
library(ggplot2)
library(scales)
library(ggdist)
library(patchwork)
library(dutchmasters)
library(ggside)
library(ggridges)
library(here)
library(lubridate)
library(ggpubr)
library(latex2exp)
library(wesanderson)
library(rnaturalearth)
library(rnaturalearthdata)
library(countrycode)


# Figure 1 ----------------------------------------------------------------


# Globals -----------------------------------------------------------------

low = "#ddecff"
mid = "#9cc8ff"
high = "#00295D"
na_color = "#8E98A4"
space = '           '

pal <- dutchmasters_pal(palette = 'milkmaid')(13)
pal2 <-wes_palette("Zissou1",21, type = 'continuous')
pal_capacity <- dutchmasters_pal(palette = 'pearl_earring')(13)



# Specify lineages --------------------------------------------------------

lineages <- c('BA.1', 
              'BA.2.12.1', 
              'BA.4',
              'BA.5')

countries<- c('Bangladesh', 
              'Brazil', 
              'Denmark', 
              'India', 
              'Portugal', 
              'Senegal', 
              'South Africa', 
              'United Kingdom', 
              'United States')

# Generate color palette --------------------------------------------------

pal <- dutchmasters_pal(palette = 'milkmaid')(13)
pal2 <-wes_palette("Zissou1",21, type = 'continuous')

lineage_colors <- c('BA.1' = pal[10],
                    'BA.2.12.1' = pal[4],
                    'BA.4' = pal[13],
                    'BA.5' = pal[12],
                    'BA.2' = pal[11])

capacity_colors = c('Top 10%' = pal_capacity[3],
                    'Bottom 90%' = pal_capacity[4])

model_colors <-c('Multicountry' = pal[12],
                 #'MLE' = 'seagreen',
                 'Single country' = 'seagreen')
model_colors_muted <-c('Multicountry' = pal[2],
                       #'MLE' = 'darkseagreen',
                       'Single country' = 'darkseagreen')

# Shared figure settings --------------------------------------------------

y_intercept_alpha = .2
axis_text_size = 10
axis_title_size = 15
tag_size = 22
strip_text_size = 12


# Data --------------------------------------------------------------------

world <- ne_countries(scale = "medium", returnclass = "sf")

world %>% 
  select(adm0_a3, name, geometry) %>% 
  rename(iso3c = adm0_a3)

df <- read_csv('../data/raw/2022-07-01_metadata.csv')

# Prep data for fig -------------------------------------------------------

processed <- df %>% 
  janitor::clean_names() %>% 
  filter(pango_lineage == 'BA.5',
         collection_date > ymd('2022-3-1')) %>% 
  mutate(iso3c = countrycode(sourcevar = country,
                             origin = 'country.name',
                             destination = 'iso3c'),
         iso3c = case_when(country == 'Bonaire' ~ 'NLD',
                           country == 'Saint Martin' ~ 'NLD',
                           TRUE ~ iso3c)) %>% 
  group_by(iso3c, collection_date) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  full_join(tibble(iso3c = NA,
                   collection_date = seq.Date(from = ymd('2022-3-1'),
                                              to = ymd('2022-7-1'),
                                              by = 'day'))) %>% 
  complete(iso3c, collection_date, fill = list(n = 0)) %>% 
  filter(!is.na(collection_date), !is.na(iso3c)) %>% 
  group_by(iso3c) %>% 
  arrange(collection_date) %>% 
  mutate(cum_n = cumsum(n)) %>% 
  ungroup() %>% 
  mutate(country.name = countrycode(sourcevar = iso3c,
                                    destination = 'country.name',
                                    origin = 'iso3c')) %>% 
  full_join(world %>% 
              select(adm0_a3, name, geometry) %>% 
              rename(iso3c = adm0_a3)
  ) %>% 
  mutate(cum_n = if_else(is.na(cum_n), as.integer(0), cum_n))


# Generate panel A --------------------------------------------------------

#p1 <- processed %>% 
#  filter(collection_date == ymd('2022-5-1') | is.na(collection_date), name != 'Antarctica') %>% 
#  ggplot(aes(geometry = geometry,
#             fill = log10(cum_n)))+
#  geom_sf( show.legend = F)+
#  scale_fill_gradient2(low = low,
#                       mid = mid,
#                       high =  high,
#                       na.value = na_color,
#                       limits = c(0, 4.1))+
#  theme_void()+
#  labs(tag = 'A',
#       title  = paste0(space, '2022-05-01'))
#
#
#p2 <- processed %>% 
#  filter(collection_date == ymd('2022-6-1') | is.na(collection_date), name != 'Antarctica') %>% 
#  ggplot(aes(geometry = geometry,
#             fill = log10(cum_n)))+
#  geom_sf(show.legend = F)+
#  scale_fill_gradient2(low = low,
#                       mid = mid,
#                       high =  high,
#                       na.value = na_color,
#                       limits = c(0, 4.1))+
#  theme_void()+
#  labs(title  = paste0(space, '2022-06-01'))
#
#
#p3 <- processed %>% 
#  filter(collection_date == ymd('2022-7-1') | is.na(collection_date), name != 'Antarctica') %>% 
#  ggplot(aes(geometry = geometry,
#             fill = log10(cum_n)))+
#  geom_sf()+
#  scale_fill_gradient2(low = low,
#                       mid = mid,
#                       high =  high,
#                       na.value = na_color,
#                       limits = c(0, 4.1))+
#  theme_void()+
#  labs(title  = paste0(space, '2022-07-01'))+
#  theme(legend.position = 'bottom')
#

inset <- df %>% 
  group_by(country) %>% 
  summarize(n = n()) %>%
  mutate(p = n/sum(n)) %>%
  arrange(desc(p)) %>% 
  filter(p > .0001) %>% 
  mutate(r = row_number(),
         top10 = if_else(r <= 21, 'Top 10%', 'Bottom 90%'),
         p = cumsum(p)) %>% 
  mutate(country = as.factor(country), 
         country = reorder(country, p, sum),) %>% 
  ggplot(aes(country, p, fill = top10))+
  geom_col(show.legend = F)+
  geom_hline(aes(yintercept = 1), alpha = 0.2)+
  geom_hline(aes(yintercept = 0.903), alpha = 1, linetype = 2)+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 9))+
  scale_fill_manual(values = capacity_colors)+
  labs(x = 'Country',
       y = 'Cumulative\nproportion')

p1 <- df %>% 
  janitor::clean_names() %>% 
  mutate(iso3c = countrycode(sourcevar = country,
                             origin = 'country.name',
                             destination = 'iso3c'),
         iso3c = case_when(country == 'Bonaire' ~ 'NLD',
                           country == 'Saint Martin' ~ 'NLD',
                           TRUE ~ iso3c)) %>% 
  group_by(iso3c) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  mutate(log10 = log10(n)) %>% 
  full_join(world %>% 
              select(adm0_a3, name, geometry) %>% 
              rename(iso3c = adm0_a3)) %>% 
  filter(name != 'Antarctica') %>% 
  ggplot(aes(geometry = geometry,
             fill = log10))+
  geom_sf(size = 0.1)+
  theme_void()+
  labs(tag = 'A',
       fill = TeX(r"(Cumulative number of sequences ($log_{10}$))"))+
  scale_fill_gradient(low = 'white',
                      high = "#A65141")+
  theme(legend.position = 'bottom',
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        #panel.background = element_rect(fill = "lightblue"),
        axis.line.y = element_blank(), 
        axis.line.x = element_blank())+
  theme(
    plot.tag = element_text(size = tag_size),
    legend.position= "top",
    legend.text = element_text(size=15),
    legend.title = (element_text(size = 15)))+
  annotation_custom(
    ggplotGrob(inset), 
    xmin = -200, xmax = -100, ymin = -65, ymax = 15
  )



p5_dat <- df %>% 
  janitor::clean_names() %>% 
  filter(!is.na(pango_lineage),
         !is.na(collection_date)) %>% 
  filter((collection_date > ymd('2021-8-1') & pango_lineage == 'BA.1') | pango_lineage != 'BA.1') %>% 
  group_by(pango_lineage, collection_date) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  group_by(pango_lineage) %>% 
  arrange(collection_date) %>% 
  mutate(cum_n = cumsum(n)) %>%
  filter(max(cum_n) > 5000) %>% 
  ungroup() %>% 
  filter(cum_n >= 100, cum_n <= 10000) %>% 
  group_by(pango_lineage) %>% 
  mutate(t_max = as.integer(max(collection_date) - min(collection_date)),
         t = as.integer(collection_date - min(collection_date)) + 1) %>% 
  ungroup() %>%   
  full_join(tibble(pango_lineage = NA,
                   t = 1:520)) %>% 
  complete(pango_lineage, t, fill = list(n = 0)) %>% 
 
   filter(!is.na(pango_lineage)) %>% 
  filter(cum_n > 0)

p2_inset <- p5_dat %>% 
  filter(cum_n > 500) %>% 
  group_by(pango_lineage) %>% 
  filter(collection_date == min(collection_date)) %>% 
  ggplot()+
  geom_histogram(aes(t), bins = 10, color = 'black', alpha = 0.3)+
  geom_segment(aes(x = t, xend = t, y = 0, yend = Inf, color = pango_lineage, group = pango_lineage) , 
               size = 1.3,
               show.legend = F,
               position = position_dodge(width = .2),
               data = p5_dat %>% 
                 filter(cum_n > 500) %>% 
                 group_by(pango_lineage) %>% 
                 filter(collection_date == min(collection_date)) %>% 
                 filter(pango_lineage %in%
                          c('BA.4')))+
  geom_segment(aes(x = t, xend = t, y = 0, yend = Inf, color = pango_lineage) , 
               size = 1.3,
               linetype = 2,
               show.legend = F,,
             data = p5_dat %>% 
               filter(cum_n > 500) %>% 
               group_by(pango_lineage) %>% 
               filter(collection_date == min(collection_date)) %>% 
               filter(pango_lineage %in%
                        c('BA.1', 'BA.2', 'BA.4', 'BA.5', 'BA.2.12.1')))+
  #ylim(0, 25)+
  theme_bw()+
  scale_x_log10()+
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8))+
  scale_color_manual(values = lineage_colors)+
  labs(x = 'Time to 500 sequences',
       y = 'Lineage')

p2 <- ggplot()+
  geom_line(aes(t, 
                cum_n, 
                group = pango_lineage), 
            alpha = .1,
            show.legend = F,
            data = p5_dat %>% 
              filter(!(pango_lineage %in% 
                         c('BA.1', 'BA.2', 'BA.4', 'BA.5', 'BA.2.12.1'))))+
  geom_line(aes(t, 
                cum_n, 
                group = pango_lineage,
                color = pango_lineage),
            alpha = 1, 
            show.legend = T,
            size = 1.4,
            data = p5_dat %>% 
              filter(pango_lineage %in% 
                       c('BA.1', 'BA.2', 'BA.4', 'BA.5', 'BA.2.12.1')))+
  geom_hline(aes(yintercept = 500), alpha = .5)+
  scale_color_manual(values = lineage_colors)+
  theme_bw()+
  labs(y = TeX(r"(Cumulative number of sequences)"),
       x = 'Days from emergence',
       tag = 'B',
       color = element_blank())+
  scale_y_log10()+
  scale_x_log10()+
  theme(
    axis.text = element_text(size = axis_text_size),
    axis.title = element_text(size = axis_title_size),
    plot.tag = element_text(size = tag_size),
    strip.text = element_text(size = strip_text_size),
    legend.position= "top",
    strip.background = element_blank(), #  element_rect(colour = "black", fill = NA)
    legend.text = element_text(size=15))

p3 <- df %>% 
  janitor::clean_names() %>% 
  filter(pango_lineage == 'BA.5') %>% 
  group_by(country, collection_date) %>% 
  summarize(n = n(),
            .groups = 'drop_last') %>% 
  arrange(collection_date) %>% 
  mutate(cum_n = cumsum(n),
         max_cum_n = max(cumsum(n))) %>% 
  ungroup() %>% 
  filter(max_cum_n > 10) %>% 
  complete(country, collection_date, fill = list(n=0)) %>% 
  group_by(country) %>% 
  arrange(collection_date) %>% 
  mutate(cum_n = cumsum(n),
         istop = if_else(country %in% c("USA", "United Kingdom", "Germany", "Denmark", "Canada", "France", 
                                        "Japan", "India", "Sweden", "Brazil", "Switzerland", "Spain", 
                                        "Austria", "Italy", "Belgium", "Netherlands", "Australia", "Israel", 
                                        "Turkey", "Poland", "Ireland"),
                         "Top 10%",
                         "Bottom 90%")) %>% 
  filter(collection_date %in% c( ymd('2022-4-15'),
                                 ymd('2022-5-15'),
                                 ymd('2022-6-15'))) %>% 
  ungroup() %>% 
  mutate(collection_date = as.factor(collection_date),
         country = reorder(country, -cum_n, min)) %>% 
  ggplot(aes(country, cum_n, fill = istop))+
  geom_col()+
  facet_wrap(~collection_date, scales = 'free_x')+
  scale_y_log10()+
  theme_bw()+
  theme(legend.position = 'bottom',
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  scale_fill_brewer(palette='Set1', direction = -1)+
  labs(y = 'Cumulative BA.5 sequences',
       x = 'Country')


p4 <- df %>% 
  janitor::clean_names() %>% 
  filter(pango_lineage %in% c('BA.5', 'BA.1', 'BA.2', 'BA.4')) %>% 
  group_by(country, collection_date, pango_lineage) %>% 
  summarize(n = n()) %>% 
  group_by(country, pango_lineage) %>% 
  arrange(collection_date) %>% 
  mutate(cum_n = cumsum(n),
         max_cum_n = max(cumsum(n))) %>% 
  ungroup() %>% 
  filter(max_cum_n > 10) %>% 
  complete(country, collection_date, pango_lineage, fill = list(n=0)) %>% 
  group_by(country, pango_lineage) %>% 
  arrange(collection_date) %>% 
  mutate(cum_n = cumsum(n),
         istop = if_else(country %in% c("USA", "United Kingdom", "Germany", "Denmark", "Canada", "France", 
                                        "Japan", "India", "Sweden", "Brazil", "Switzerland", "Spain", 
                                        "Austria", "Italy", "Belgium", "Netherlands", "Australia", "Israel", 
                                        "Turkey", "Poland", "Ireland"),
                         "Top 10%",
                         "Bottom 90%")) %>% 
  #filter(collection_date > ymd('2022-4-1')) %>% 
  ungroup() %>% 
  mutate(country = as.factor(country), 
         country = reorder(country, -cum_n, sum)) %>% 
  group_by(collection_date, istop, pango_lineage) %>% 
  summarize(cum_n = sum(cum_n)) %>% 
  ungroup() %>% 
  group_by(pango_lineage, collection_date) %>% 
  filter(sum(cum_n) > 200, sum(cum_n) < 100000) %>% 
  ggplot(aes(collection_date, cum_n, fill=istop))+
  geom_col(position = 'fill', show.legend = T)+
  labs(y = 'Proportion of\nsequences',
       x = element_blank(),
       tag = 'C',
       fill = 'Sequencing Capacity:')+
  scale_fill_manual(values = capacity_colors)+
  theme_bw()+
  facet_wrap(~pango_lineage, scale = 'free_x', nrow = 1)+
  theme(
    axis.text = element_text(size = axis_text_size),
    axis.title = element_text(size = axis_title_size),
    plot.tag = element_text(size = tag_size),
    strip.text = element_text(size = strip_text_size),
    legend.position= "bottom",
    strip.background = element_blank(), #  element_rect(colour = "black", fill = NA)
    legend.text = element_text(size=15),
    legend.title = element_text(size = 15))
  
  
  
# Top 10% of sequencing countries
#p3 <- df %>% 
#  janitor::clean_names() %>% 
#  filter(pango_lineage %in%  c('BA.1', 'BA.2', 'BA.4', 'BA.5')) %>% 
#  mutate(istop = if_else(country %in% c("USA", "United Kingdom", "Germany", "Denmark", "Canada", "France", 
#                                "Japan", "India", "Sweden", "Brazil", "Switzerland", "Spain", 
#                                "Austria", "Italy", "Belgium", "Netherlands", "Australia", "Israel", 
#                                "Turkey", "Poland", "Ireland"),
#                         "Top 10%",
#                         "Bottom 90%"),
#         istop = factor(istop, c("Top 10%",
#                                 "Bottom 90%"))) %>% 
#  filter(!is.na(pango_lineage),
#         !is.na(collection_date),
#         !is.na(istop)) %>% 
#  group_by(pango_lineage, collection_date, istop) %>% 
#  summarize(n = n()) %>% 
#  ungroup() %>% 
#  mutate(t_max = as.integer(max(collection_date) - min(collection_date)),
#         t = as.integer(collection_date - min(collection_date)) + 1) %>% 
#  ungroup() %>%   
#  full_join(tibble(pango_lineage = NA,
#                   t = 1:1000)) %>% 
#  complete(pango_lineage, t, istop, fill = list(n = 0)) %>% 
#  filter(!is.na(t),
#         !is.na(istop)) %>% 
#  group_by(pango_lineage, istop) %>%
#  arrange(t) %>% 
#  mutate(cum_n = cumsum(n)) %>%
#  ungroup() %>% 
#  filter(!is.na(pango_lineage)) %>% 
#  group_by(t, pango_lineage) %>%
#  filter(sum(cum_n) < 10000, sum(cum_n) > 500) %>%
#  ungroup() %>% 
#  group_by(pango_lineage, t) %>% 
#  filter(sum(cum_n) > 1) %>% 
#  ungroup() %>% 
#  group_by(pango_lineage) %>% 
#  mutate(t = t - min(t) + 1) %>% 
#  ggplot(aes(t, cum_n, fill = istop))+
#  geom_col(position = 'stack')+
#  facet_wrap(~pango_lineage, scale = 'free_x', nrow = 1)+
#  scale_fill_brewer(palette = 'Set1')+
#  theme_bw()+
#  labs(x = 'Days since 500th sequence collected',
#       y = 'Cumulative number\nof sequences',
#       fill = 'Sequencing capacity',
#       tag = 'C')+
#    theme(
#      axis.text = element_text(size = axis_text_size),
#      axis.title = element_text(size = axis_title_size),
#      plot.tag = element_text(size = tag_size),
#      legend.position= "bottom",
#      legend.text = element_text(size=15))+
#  theme(
#    axis.text = element_text(size = axis_text_size),
#    axis.title = element_text(size = axis_title_size),
#    plot.tag = element_text(size = tag_size),
#    strip.text = element_text(size = strip_text_size),
#    strip.background = element_blank(), #  element_rect(colour = "black", fill = NA)
#    legend.text = element_text(size=15),
#    legend.title = element_text(size = 15))
  


# Top 10% of sequencing countries
p3 <- df %>% 
  janitor::clean_names() %>% 
  filter(pango_lineage %in%  c('BA.1', 'BA.2', 'BA.4', 'BA.5', 'BA.2.12.1')) %>% 
  mutate(istop = if_else(country %in% c("USA", "United Kingdom", "Germany", "Denmark", "Canada", "France", 
                                "Japan", "India", "Sweden", "Brazil", "Switzerland", "Spain", 
                                "Austria", "Italy", "Belgium", "Netherlands", "Australia", "Israel", 
                                "Turkey", "Poland", "Ireland"),
                         "Top 10%",
                         "Bottom 90%"),
         istop = factor(istop, c("Top 10%",
                                 "Bottom 90%"))) %>% 
  filter(!is.na(pango_lineage),
         !is.na(collection_date),
         !is.na(istop)) %>% 
  group_by(pango_lineage, collection_date, istop) %>% 
  summarize(n = n()) %>% 
  ungroup() %>% 
  mutate(t_max = as.integer(max(collection_date) - min(collection_date)),
         t = as.integer(collection_date - min(collection_date)) + 1) %>% 
  ungroup() %>%   
  full_join(tibble(pango_lineage = NA,
                   t = 1:1000)) %>% 
  complete(pango_lineage, t, istop, fill = list(n = 0)) %>% 
  filter(!is.na(t),
         !is.na(istop)) %>% 
  group_by(pango_lineage, istop) %>%
  arrange(t) %>% 
  mutate(cum_n = cumsum(n)) %>%
  ungroup() %>% 
  filter(!is.na(pango_lineage)) %>% 
  group_by(t, pango_lineage) %>%
  filter(sum(cum_n) < 10000, sum(cum_n) > 500) %>%
  ungroup() %>% 
  group_by(pango_lineage, t) %>% 
  filter(sum(cum_n) > 1) %>% 
  ungroup() %>% 
  group_by(pango_lineage) %>% 
  mutate(t = t - min(t) + 1) %>% 
  ggplot(aes(t, cum_n, fill = istop))+
  geom_col(position = 'stack')+
  facet_wrap(~pango_lineage, scale = 'free', nrow = 1)+
  scale_fill_brewer(palette = 'Set1')+
  theme_bw()+
  labs(x = 'Days since 500th sequence collected',
       y = 'Cumulative number\nof sequences',
       fill = 'Sequencing capacity',
       tag = 'C')

design <- "
 11122
 11122
 11122
 33333
"

fig <- p1 + p2 + p4 + plot_layout(design = design)

ggsave('../data/output/figures/figure_1.tiff', 
       plot = fig,
       width = 15,
       height = 8)
  
  
# Figure 3 ----------------------------------------------------------------


lineage_colors <- c('BA.1' = pal[10],
                    'BA.2.12.1' = pal[4],
                    'BA.4' = pal[13],
                    'BA.5' = pal[12],
                    'BA.2' = pal[11],
                    'other' = pal[6])


# Load data ---------------------------------------------------------------

#setwd(here::here("code"))
clean_global_df <- read_csv('../data/processed/validation_data/clean_global_df_2022-07-01.csv') %>% 
    filter(collection_date <= ymd("2022-07-01"))
r_summary <- read_csv('../data/processed/validation_data/r_summary_2022-07-01.csv')
mu_all <- read.csv('../data/processed/multicountry_mu_distrib_2022-07-01.csv')
variant_t <-read_csv('../data/processed/variant_t_2022-07-01.csv')
variant_draws <-read_csv('../data/processed/multicountry_model_pred_p_2022-07-01.csv') %>% 
    filter(date <= ymd("2022-07-01"))



# Test out that the model results make sense
clean_global_df %>% filter(country %in% countries) %>% 
    ggplot() + geom_point(aes(x =mid_week_date, y = mid_week_p_lineage, color = lineage))+
    geom_linerange(aes(x = mid_week_date, ymin = mid_week_p_lineage-mid_week_p_lineage_se,
                       ymax = mid_week_p_lineage+mid_week_p_lineage_se, color = lineage))+
    geom_line(aes(x = collection_date, y = p_hat_med, color = lineage)) +
    geom_ribbon(aes(x = collection_date, ymin = p_hat_lb, ymax = p_hat_ub, fill = lineage), alpha = 0.4) +
    facet_wrap(~country)

# Generate Figure panel A -------------------------------------------------

variant_summarized <- variant_draws %>% 
  mutate(variant = if_else(lineage %in% c(lineages, 'BA.2', 'other'), 
                           lineage,
                           'other')) %>% 
  group_by(country, date, draw, variant) %>% 
  summarise(p_est = sum(p_hat), # combines lineages were not tracking into other
            multicountry_period = max(multicountry_period)) %>% # keeps multicountry period variable
  ungroup() %>% 
    filter(country %in% countries) %>% 
    left_join(variant_t %>% filter(country %in% countries),
              by = c("variant", "country", "date" = "collection_date"))

most_obs_lineage <- clean_global_df %>% 
    filter(country %in% countries) %>% 
    ungroup() %>% 
    filter(n > 0, 
           !is.na(n)) %>% 
    group_by(country, collection_date) %>% 
    filter(n == max(n)) %>% 
    mutate(n_remaining = n()) %>% 
    mutate(lineage = if_else(n_remaining == 1,
                             lineage,
                             NA_character_)) %>% 
    fill(lineage, .direction = 'updown') %>% 
    ungroup()  %>% 
    select(collection_date, country, lineage) %>% 
    rename(most_obs_lineage = lineage,
           date = collection_date) %>% 
    mutate(most_obs_lineage = if_else(most_obs_lineage %in% c(lineages, 'BA.2'),
                                      most_obs_lineage,
                                      'other'))

cases_calib<- variant_summarized %>% 
  filter(multicountry_period == 'multicountry calibration') %>% 
    left_join(most_obs_lineage)
cases_nowcast<- variant_summarized %>% 
  filter(multicountry_period == 'multicountry nowcast')

p1 <- ggplot()+
    geom_line(mapping = aes(date, p_est, group = interaction(draw, variant), color = variant),
              data = cases_nowcast, 
              alpha = .05,
              show.legend = T)+
    geom_line(mapping = aes(date, p_est, group = interaction(draw, variant), color = variant), 
              data = cases_calib, 
              alpha = .1,
              show.legend = F)+
    facet_wrap(~country)+
    scale_color_manual(values = lineage_colors)+
    scale_fill_manual(values = lineage_colors)+
    guides(colour = guide_legend(override.aes = list(alpha = 1,
                                                     size = 4),
                                 nrow = 1))+
    theme_bw()+
    theme(ggside.panel.scale = .3,
          axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          strip.text = element_text(size = strip_text_size),
          strip.background = element_blank(),
          legend.position= "bottom",
          legend.text = element_text(size=15))+
    labs(y = 'Estimated variant proportion',
         color = element_blank(),
         x = element_blank(),
         fill = element_blank(),
         tag = 'A')+
    geom_xsidecol(aes(x = date, y = n, fill = most_obs_lineage),
                  data = variant_t %>% 
                    filter(country %in% countries) %>% 
                    group_by(collection_date, country) %>% 
                    summarize(n = log10(sum(n_seq))) %>%
                    rename(date = collection_date) %>% 
                    left_join(most_obs_lineage) %>% 
                    ungroup() %>% 
                    mutate(most_obs_lineage = factor(most_obs_lineage)),
                  show.legend = F,
                  alpha = 1)+
    scale_xsidey_continuous(position = 'right',
                            minor_breaks = NULL,
                            n.breaks = 3,
                            guide = guide_axis(TeX(r"(Number of collected sequences (log$_{10}$))")))

# Generate panel B --------------------------------------------------------

p2 <- r_summary %>% 
    filter(lineage %in% lineages,
           country %in% countries) %>% 
    mutate(country = factor(country, ordered = T)) %>% 
    mutate(country = reorder(country, desc(country))) %>% 
    ggplot()+
    geom_vline(aes(xintercept = 0), 
               alpha = y_intercept_alpha)+
    geom_line(aes(x = global_trans_adv_median,
                  y =  country,
                  group = lineage,
                  color = lineage), 
              alpha = 1,

              show.legend = F)+
    geom_ribbon(aes(y = country, 
                    xmin = global_trans_adv_lb, 
                    xmax = global_trans_adv_ub,
                    group = lineage,
                    fill = lineage), 
                alpha = .25,
                show.legend = F)+
    geom_point(aes(y = country, 
                   x = trans_adv_median,
                   color = lineage),
               show.legend = F)+
    geom_linerange(aes(y = country, 
                       xmin = trans_adv_lb, 
                       xmax = trans_adv_ub, 
                       color = lineage),
                   show.legend = F)+
    facet_wrap(~lineage)+
    scale_color_manual(values = lineage_colors)+
    scale_fill_manual(values = lineage_colors)+
    theme_bw()+
    theme(axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          strip.background = element_blank(), # element_rect(colour = "black", fill = NA)
          strip.text = element_text(size = strip_text_size))+ 
    labs(x = 'Posterior of country-specific fitness advantage\nrelative to BA.2 (median with 95% CI)',
         y = element_blank(),
         tag = 'B')

# Generate panel C --------------------------------------------------------

p3 <- mu_all %>% 
    filter(lineage %in% lineages) %>% 
    unique() %>% 
    mutate(lineage = factor(lineage, ordered = T)) %>% 
    mutate(lineage = reorder(lineage, desc(lineage))) %>% 
    ggplot(aes(y = lineage, x = transmission_advantage))+
    geom_vline(aes(xintercept = 0),
               alpha = y_intercept_alpha)+
    stat_halfeye(slab_alpha = 1,
                 aes(fill = lineage),
                 show.legend = F)+
    theme_bw()+
    theme(legend.position=c(.2,.4),
          axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    labs(y = element_blank(), 
         x = 'Posterior distribution of global mean weekly\nfitness advantage relative to BA.2 (95% CI)',
         tag = 'C',
         fill = element_blank())+
    scale_fill_manual(values = lineage_colors)

# Put together panels -----------------------------------------------------

design <- "
  112
  113
"

fig <- p1 + p2 + p3 + plot_layout(design = design)

# Save fig ----------------------------------------------------------------

ggsave('../data/output/figures/figure_3.pdf', 
       plot = fig,
       width = 15,
       height = 8)


# Load in data for figure 4 ---------------------------------------------------

r_summary<- read_csv('../data/output/validation/r_summary.csv')
clean_global_df <- read_csv('../data/output/validation/clean_global_df.csv')
mu_hat <- read_csv('../data/output/validation/mu_hat.csv')
mu_distrib <- read_csv('../data/output/validation/mu_distrib.csv')
r_distrib <- read_csv('../data/output/validation/r_distrib.csv')
country_metrics <- read_csv('../data/output/validation/country_metrics.csv')
country_lineage_metrics <- read_csv('../data/output/validation/country_lineage_metrics.csv')

# Clean the data for visualizations -------------------------------------------
# Use this to decide on which countries to visualize
country_early <- country_metrics %>% 
    filter(reference_date == ymd('2022-04-30')) %>% 
    mutate(BS_improvement = BS_med_MLE - BS_med_MC,
           CCC_improvement = CCC_country_med_MC - CCC_country_med_MLE,
           BS_period_improvement = BS_med_period_MLE - BS_med_period_MC,
           CCC_period_improvement = CCC_country_period_med_MC - CCC_country_period_med_MLE) %>% 
    unique() %>% 
    arrange(desc(BS_improvement)) %>% 
    relocate(BS_improvement, .before = multicountry_period)

BA5_early <- country_lineage_metrics %>% 
    filter(lineage == 'BA.5', 
           reference_date == ymd('2022-04-30')) %>% 
    mutate(CCC_improvement = CCC_country_lin_med_MC - CCC_country_lin_med_MLE,
           CCC_period_improvement = CCC_country_lin_period_med_MC - CCC_country_lin_period_med_MLE) %>% 
    arrange(desc(CCC_improvement)) %>% 
    relocate(CCC_improvement, .before = multicountry_period)

# Make a clean table of metrics for plotting
cleaned_metrics_MC <- country_metrics %>% 
    select(country, reference_date, multicountry_period,
           BS_med_MC, BS_mean_MC, BS_lb_MC, BS_ub_MC,
           BS_med_period_MC, BS_mean_period_MC, BS_lb_period_MC, BS_ub_period_MC) %>% 
    rename(BS_med = BS_med_MC,
           BS_mean = BS_mean_MC,
           BS_lb = BS_lb_MC,
           BS_ub = BS_ub_MC,
           BS_med_period = BS_med_period_MC,
           BS_mean_period = BS_mean_period_MC,
           BS_lb_period = BS_lb_period_MC,
           BS_ub_period = BS_ub_period_MC,) %>% 
    mutate( model = 'Multicountry')
cleaned_metrics_MLE <- country_metrics %>% 
    select(country, reference_date, multicountry_period,
           BS_med_MLE, BS_mean_MLE, BS_lb_MLE, BS_ub_MLE,
           BS_med_period_MLE, BS_mean_period_MLE, BS_lb_period_MLE, BS_ub_period_MLE) %>% 
    rename(BS_med = BS_med_MLE,
           BS_mean = BS_mean_MLE,
           BS_lb = BS_lb_MLE,
           BS_ub = BS_ub_MLE,
           BS_med_period = BS_med_period_MLE,
           BS_mean_period = BS_mean_period_MLE,
           BS_lb_period = BS_lb_period_MLE,
           BS_ub_period = BS_ub_period_MLE) %>% 
    mutate( model = 'Single country')
cleaned_metrics <- bind_rows(cleaned_metrics_MC, cleaned_metrics_MLE)
cleaned_metrics <- cleaned_metrics %>% left_join(cleaned_metrics %>% group_by(reference_date) %>% 
                                                     summarise(n = n()) %>% 
                                                     mutate(num = row_number())) 



# Make figure of global fitness advantage over time
mu_distrib <- mu_distrib %>% mutate(ref_date_num = as.numeric(reference_date))
r_distrib <- r_distrib %>% mutate(ref_date_num = as.numeric(reference_date))

r_comb<- r_distrib %>%
    filter(country %in% countries) %>% 
    left_join(r_summary %>% 
                                     select(country, lineage, reference_date,
                                            trans_adv_median, 
                                            trans_adv_lb, 
                                            trans_adv_ub,
                                            trans_adv_MLE_mean,
                                            trans_adv_MLE_lb,
                                            trans_adv_MLE_ub,
                                            #trans_adv_normal_approx_lb,
                                            #trans_adv_normal_approx_ub,
                                            global_trans_adv_median,
                                            global_trans_adv_lb,
                                            global_trans_adv_ub,),
                                 by = c("reference_date", "country", "lineage")) %>% 
    left_join(r_summary %>% group_by(reference_date) %>% 
                  summarise(n= n()) %>% 
                  mutate(num = row_number()))

# add total number of sequences by reference date and country 
r_summary <- r_summary %>%  left_join(clean_global_df %>% 
                                     group_by(country, reference_date, lineage) %>% 
                                     summarise(n_seq = sum(n, na.rm = T)),
                                               by = c("reference_date", "country", "lineage"))

mu_comb <- mu_distrib %>% left_join(mu_hat %>% select(lineage, reference_date, global_trans_adv_median,
                                                      global_trans_adv_lb,
                                                      global_trans_adv_ub),
                                    by = c("reference_date", "lineage")) %>% 
    left_join(mu_hat %>% group_by(reference_date) %>% 
                  summarise(n = n()) %>% 
                  mutate(num = row_number()))

# Global BA.5 posterior
fig4global <- mu_comb %>% 
    filter(lineage == 'BA.5') %>% 
    ggplot(aes( x = as.factor(num), y = global_transmission_advantage))+
    stat_halfeye(slab_alpha = 0.7,
                 #interval_color = 'gray',
                 #point_color = 'gray',
                 aes(fill = lineage),
                 show.legend = F)+
    theme_bw()+
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    scale_fill_manual(values = lineage_colors)+
    theme(axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2,size = 0.7*axis_text_size),
          axis.title.y = element_text(size = axis_title_size),
          axis.title.x = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    labs(y = paste0('Global weekly fitness advantage of BA.5 relative to BA.2'),
         x = 'Reference date',
         fill = '',
         tag = 'B')
fig4global


# Portugal BA.5-------------------------------------------------

this_country<- 'Portugal'

get_single_country_metrics <- function(this_country){
    global_seq<- clean_global_df %>% group_by(collection_date, reference_date) %>% 
        summarise(seq = log10(sum(N, na.rm = T))) %>% 
        ungroup() %>% 
        mutate(method = 'Multicountry')
    country_df<-clean_global_df %>% filter(country == this_country, lineage == "BA.5")
    country_seq <- country_df  %>% 
        mutate(seq = log10(N)) %>% 
        select(collection_date, reference_date, seq) %>% 
        distinct() %>% 
        mutate(method = 'Single country')
    
    # Create data to make the top panel
    seq_numbers <- bind_rows(global_seq, country_seq)
    seq_numbers$method <- factor(seq_numbers$method, levels = c("Multicountry", "Single country"))
    seq_summary <- seq_numbers %>% filter(method == 'Single country') 
    ref_dates<-unique(seq_summary$reference_date)
    df<- c()
    for (i in 1:length(ref_dates)){
        subset<- seq_summary %>% filter(reference_date == ref_dates[i])
        reference_date <- ref_dates[i]
        last_date_w_obs_seq <- max(subset$collection_date[subset$seq>0], na.rm = T)
        dfi <- data.frame(reference_date, last_date_w_obs_seq)
        df<-bind_rows(df, dfi)
    }
    
    
    
    single_country_df <- country_df %>% select( collection_date,mid_week_date_recent, mid_week_prev_recent, mid_week_prev_se_recent,
                                                mid_week_date, mid_week_p_lineage, mid_week_p_lineage_se, multicountry_period,
                                                reference_date, p_hat_MLE_mean, p_hat_MLE_lb, p_hat_MLE_ub) %>% 
        rename(p_hat_mean = p_hat_MLE_mean,
               p_hat_lb = p_hat_MLE_lb,
               p_hat_ub = p_hat_MLE_ub) %>% 
        left_join(df, by = 'reference_date') %>% 
        mutate(method = 'Single country',
               period = ifelse(collection_date <= last_date_w_obs_seq, 'calibration', 'forecast'))
    multicountry_df<- country_df %>% select(collection_date,mid_week_date_recent, mid_week_prev_recent, mid_week_prev_se_recent,
                                            mid_week_date, mid_week_p_lineage, mid_week_p_lineage_se, multicountry_period,
                                            reference_date, p_hat_mean, p_hat_lb, p_hat_ub) %>% 
        left_join(df, by = 'reference_date') %>% 
        mutate(method = 'Multicountry',
               period = ifelse(multicountry_period == 'multicountry calibration', 'calibration', 'forecast'))
    combined_df <- bind_rows(single_country_df, multicountry_df)
    
    combined_df$method <- factor(combined_df$method, levels = c("Multicountry", "Single country"))
    
    return(combined_df)
}
combined_df <- get_single_country_metrics(this_country)
calib_df<- combined_df %>% filter(period == 'calibration' )
forecast_df <- combined_df %>% filter(period == 'forecast')

ggplot() + geom_col(data = seq_numbers, aes(x = collection_date, y = seq, fill = method))+ facet_grid(method~reference_date, scales = 'free')
fig4a <- ggplot() +
    geom_line(data = calib_df, aes(x = ymd(collection_date), y = p_hat_mean, color = method), size =1) +
    geom_ribbon(data = calib_df, aes(x = ymd(collection_date), ymin = p_hat_lb, ymax = p_hat_ub, fill = method), alpha = 0.4, show.legend = F)+
    geom_line(data = forecast_df, aes(x = ymd(collection_date), y = p_hat_mean, color = method), size =0.2, show.legend = F) +
    geom_ribbon(data = forecast_df, aes(x = ymd(collection_date), ymin = p_hat_lb, ymax = p_hat_ub, fill = method), alpha = 0.2, show.legend = F)+
    geom_point(data = country_df, aes(x = ymd(mid_week_date_recent), y = mid_week_prev_recent), size = 1) +
    geom_linerange(data = country_df, aes(x = ymd(mid_week_date_recent), ymin = mid_week_prev_recent -mid_week_prev_se_recent, 
                                          ymax = mid_week_prev_recent + mid_week_prev_se_recent)) +
    geom_xsidecol(data = seq_numbers, aes(x = ymd(collection_date), y = seq, fill = method),
                  show.legend = F, alpha = 1)+
    scale_xsidey_continuous(position = 'right',
                            minor_breaks = NULL,
                            n.breaks = 3,
                            guide = guide_axis(TeX(r"(Number of collected sequences (log$_{10}$))")))+
    facet_grid(method~reference_date)+
    theme_bw()+ 
    guides(colour = guide_legend(override.aes = list(alpha = 1,
                                                     size = 4),
                                 nrow = 1))+
    theme(ggside.panel.scale = .3,
          axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle =- 45, hjust=.6, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          strip.text = element_text(size = strip_text_size),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size=15),
          legend.position = "bottom",
          strip.background =element_rect(fill="white"))+
    scale_color_manual(values = model_colors)+
    scale_fill_manual(values = model_colors)+
    scale_x_date(date_labels = "%B %d")+
    #coord_cartesian(ylim = c(0,1))+
    labs(x = '',
         y = 'BA.5 prevalence',
         fill = '',
         color = 'Method',
         tag = 'A')+
    ggtitle(paste0(this_country, ' BA.5'))
fig4a



this_country_metrics <-cleaned_metrics %>% filter(country == this_country)
dodge<- position_dodge(width = 0.5)
fig4c<- ggplot()+
    geom_point(data = this_country_metrics, aes(x = as.factor(num), y = BS_mean, 
                                         color = model), size = 3.5, position = dodge, show.legend = F) +
    geom_linerange(data = this_country_metrics, aes(x = as.factor(num), ymin = BS_lb,
                                            ymax = BS_ub, color = model),size = 1.5, position = dodge, show.legend = F) +
    theme_bw()+ 
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    #scale_y_sqrt()+
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2, size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          legend.position = "none")+
    scale_color_manual(values = model_colors)+
    #scale_x_date(date_labels = "%B %d", limits = c(ymd('2022-04-25'), ymd('2022-07-15')))+
    labs(x = 'Reference date',
         y = 'Brier score',
         fill = '',
         color = 'Model',
         tag = 'C') + 
    ggtitle(paste0(this_country, ' BA.5'))
fig4c

# heatmaps!

r_summary <- r_summary %>% mutate( dif_trans_capped_MC = ifelse(dif_trans_adv_MC >4, 4, dif_trans_adv_MC),
                                   dif_trans_capped_MLE = ifelse(dif_trans_adv_MLE >4, 4, dif_trans_adv_MLE))
fig4addtop<- r_summary %>% filter(lineage %in% lineages, country %in% countries) %>% 
    ggplot() + geom_tile(aes(x = as.factor(num), y = lineage, fill = dif_trans_capped_MC)) +
    facet_wrap(~country, ncol = 3) +
    labs(x = 'Reference date',
         y = 'Lineage',
         fill = 'change in\nfitness\nadvantage',
         tag = 'D')+
    #scale_fill_continuous(type = "viridis")+
    scale_fill_gradientn(colours = pal2,
                         na.value = "lightgray", guide = "colourbar", aesthetics = "fill",
                         #breaks=c(-0.5, 0,1, 2,3, 4), labels=c(-1, 0, 1,2,3, 4),
                         limits=c(-0.5, 4))+
    # scale_colour_gradient2(low = pal2[1], mid = pal2[10], high = pal2[21], midpoint = 0,
    #                        na.value = "gray", guide = "colourbar", aesthetics = "fill",
    #                        breaks=c(-0.5, 0,1, 2, 3, 4),labels=c(-0.5, 0,1,2,3,4),
    #                        limits=c(-1, 4))+
    theme_bw()+
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    # coord_fixed()+
    theme(aspect.ratio = 1,
          axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=0, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          strip.background =element_rect(fill="white"))+
    ggtitle('Multicountry model')
fig4addtop
fig4addbottom<- r_summary %>% filter(lineage %in% lineages, country %in% countries) %>% 
    ggplot() + geom_tile(aes(x = as.factor(num), y = lineage, fill = dif_trans_capped_MLE), show.legend = F) +
    facet_wrap(~country, ncol = 3)+
    labs(x = 'Reference date',
         y = 'Lineage',
         fill = '',
         tag = 'E')+
    scale_fill_gradientn(colours = pal2,
                         na.value = "lightgray", guide = "colourbar", aesthetics = "fill",
                         #breaks=c(-0.5, 0,1,3, 2, 4), labels=c(-0.5, 0, 1,3, 2, 4),
                         limits=c(-0.5, 4))+
    theme_bw()+
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    theme(aspect.ratio = 1,
          axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=0, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          #axis.text.y = element_blank(),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          strip.background =element_rect(fill="white"))+
    ggtitle('Single country model')
fig4addbottom
heatmaps <- ggarrange(fig4addtop, fig4addbottom, ncol = 2, common.legend = T, legend = 'bottom', align = 'h')

# Try combining everything together------------------------------
design <- "
  111122
  334455
"

fig <- fig4a + fig4global + fig4c + fig4addtop + fig4addbottom+
    plot_layout(design = design)
fig

# Save main fig ----------------------------------------------------------------

ggsave('../data/output/figures/figure_4_Portugal_BA5.pdf', 
       plot = fig,
       width = 21,
       height = 15)










# Make a figure with panels A and C from a few different countries 
list_of_countries <- c('Denmark', 'Czechia', 'Brazil', 'Spain')
this_country <-list_of_countries[1]
combined_df <- get_single_country_metrics(this_country)
calib_df<- combined_df %>% filter(period == 'calibration' )
forecast_df <- combined_df %>% filter(period == 'forecast')
this_country_metrics <-cleaned_metrics %>% filter(country == this_country)


timeseries <- ggplot() +
    geom_line(data = calib_df, aes(x = ymd(collection_date), y = p_hat_mean, color = method), size =1) +
    geom_ribbon(data = calib_df, aes(x = ymd(collection_date), ymin = p_hat_lb, ymax = p_hat_ub, fill = method), alpha = 0.4, show.legend = F)+
    geom_line(data = forecast_df, aes(x = ymd(collection_date), y = p_hat_mean, color = method), size =0.2, show.legend = F) +
    geom_ribbon(data = forecast_df, aes(x = ymd(collection_date), ymin = p_hat_lb, ymax = p_hat_ub, fill = method), alpha = 0.2, show.legend = F)+
    geom_point(data = country_df, aes(x = ymd(mid_week_date_recent), y = mid_week_prev_recent), size = 1) +
    geom_linerange(data = country_df, aes(x = ymd(mid_week_date_recent), ymin = mid_week_prev_recent -mid_week_prev_se_recent, 
                                          ymax = mid_week_prev_recent + mid_week_prev_se_recent)) +
    geom_xsidecol(data = seq_numbers, aes(x = ymd(collection_date), y = seq, fill = method),
                  show.legend = F, alpha = 1)+
    scale_xsidey_continuous(position = 'right',
                            minor_breaks = NULL,
                            n.breaks = 3,
                            guide = guide_axis(TeX(r"(Number of collected sequences (log$_{10}$))")))+
    facet_grid(method~reference_date)+
    theme_bw()+ 
    guides(colour = guide_legend(override.aes = list(alpha = 1,
                                                     size = 4),
                                 nrow = 1))+
    theme(ggside.panel.scale = .3,
          axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle =- 45, hjust=.6, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          strip.text = element_text(size = strip_text_size),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size=15),
          legend.position = "bottom",
          strip.background =element_rect(fill="white"))+
    scale_color_manual(values = model_colors)+
    scale_fill_manual(values = model_colors)+
    scale_x_date(date_labels = "%B %d")+
    #coord_cartesian(ylim = c(0,1))+
    labs(x = '',
         y = 'BA.5 prevalence',
         fill = '',
         color = 'Method',
         tag = 'A')+
    ggtitle(paste0(this_country, ' BA.5'))
timeseries


Brier_score_comparison<- ggplot()+
    geom_point(data = this_country_metrics, aes(x = as.factor(num), y = BS_mean, 
                                                color = model), size = 3.5, position = dodge, show.legend = F) +
    geom_linerange(data = this_country_metrics, aes(x = as.factor(num), ymin = BS_lb,
                                                    ymax = BS_ub, color = model),size = 1.5, position = dodge, show.legend = F) +
    theme_bw()+ 
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    #scale_y_sqrt()+
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2, size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          legend.position = "none")+
    scale_color_manual(values = model_colors)+
    #scale_x_date(date_labels = "%B %d", limits = c(ymd('2022-04-25'), ymd('2022-07-15')))+
    labs(x = 'Reference date',
         y = 'Brier score',
         fill = '',
         color = 'Model',
         tag = 'B') + 
    ggtitle(paste0(this_country, ' BA.5'))
Brier_score_comparison 


this_country <-list_of_countries[2]
combined_df <- get_single_country_metrics(this_country)
calib_df<- combined_df %>% filter(period == 'calibration' )
forecast_df <- combined_df %>% filter(period == 'forecast')
this_country_metrics <-cleaned_metrics %>% filter(country == this_country)


timeseries2 <- ggplot() +
    geom_line(data = calib_df, aes(x = ymd(collection_date), y = p_hat_mean, color = method), size =1) +
    geom_ribbon(data = calib_df, aes(x = ymd(collection_date), ymin = p_hat_lb, ymax = p_hat_ub, fill = method), alpha = 0.4, show.legend = F)+
    geom_line(data = forecast_df, aes(x = ymd(collection_date), y = p_hat_mean, color = method), size =0.2, show.legend = F) +
    geom_ribbon(data = forecast_df, aes(x = ymd(collection_date), ymin = p_hat_lb, ymax = p_hat_ub, fill = method), alpha = 0.2, show.legend = F)+
    geom_point(data = country_df, aes(x = ymd(mid_week_date_recent), y = mid_week_prev_recent), size = 1) +
    geom_linerange(data = country_df, aes(x = ymd(mid_week_date_recent), ymin = mid_week_prev_recent -mid_week_prev_se_recent, 
                                          ymax = mid_week_prev_recent + mid_week_prev_se_recent)) +
    geom_xsidecol(data = seq_numbers, aes(x = ymd(collection_date), y = seq, fill = method),
                  show.legend = F, alpha = 1)+
    scale_xsidey_continuous(position = 'right',
                            minor_breaks = NULL,
                            n.breaks = 3,
                            guide = guide_axis(TeX(r"(Number of collected sequences (log$_{10}$))")))+
    facet_grid(method~reference_date)+
    theme_bw()+ 
    guides(colour = guide_legend(override.aes = list(alpha = 1,
                                                     size = 4),
                                 nrow = 1))+
    theme(ggside.panel.scale = .3,
          axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle =- 45, hjust=.6, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          strip.text = element_text(size = strip_text_size),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size=15),
          legend.position = "bottom",
          strip.background =element_rect(fill="white"))+
    scale_color_manual(values = model_colors)+
    scale_fill_manual(values = model_colors)+
    scale_x_date(date_labels = "%B %d")+
    #coord_cartesian(ylim = c(0,1))+
    labs(x = '',
         y = 'BA.5 prevalence',
         fill = '',
         color = 'Method',
         tag = 'C')+
    ggtitle(paste0(this_country, ' BA.5'))



Brier_score_comparison2<- ggplot()+
    geom_point(data = this_country_metrics, aes(x = as.factor(num), y = BS_mean, 
                                                color = model), size = 3.5, position = dodge, show.legend = F) +
    geom_linerange(data = this_country_metrics, aes(x = as.factor(num), ymin = BS_lb,
                                                    ymax = BS_ub, color = model),size = 1.5, position = dodge, show.legend = F) +
    theme_bw()+ 
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    #scale_y_sqrt()+
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2, size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          legend.position = "none")+
    scale_color_manual(values = model_colors)+
    #scale_x_date(date_labels = "%B %d", limits = c(ymd('2022-04-25'), ymd('2022-07-15')))+
    labs(x = 'Reference date',
         y = 'Brier score',
         fill = '',
         color = 'Model',
         tag = 'D') + 
    ggtitle(paste0(this_country, ' BA.5'))


this_country <-list_of_countries[3]
combined_df <- get_single_country_metrics(this_country)
calib_df<- combined_df %>% filter(period == 'calibration' )
forecast_df <- combined_df %>% filter(period == 'forecast')
this_country_metrics <-cleaned_metrics %>% filter(country == this_country)


timeseries3 <- ggplot() +
    geom_line(data = calib_df, aes(x = ymd(collection_date), y = p_hat_mean, color = method), size =1) +
    geom_ribbon(data = calib_df, aes(x = ymd(collection_date), ymin = p_hat_lb, ymax = p_hat_ub, fill = method), alpha = 0.4, show.legend = F)+
    geom_line(data = forecast_df, aes(x = ymd(collection_date), y = p_hat_mean, color = method), size =0.2, show.legend = F) +
    geom_ribbon(data = forecast_df, aes(x = ymd(collection_date), ymin = p_hat_lb, ymax = p_hat_ub, fill = method), alpha = 0.2, show.legend = F)+
    geom_point(data = country_df, aes(x = ymd(mid_week_date_recent), y = mid_week_prev_recent), size = 1) +
    geom_linerange(data = country_df, aes(x = ymd(mid_week_date_recent), ymin = mid_week_prev_recent -mid_week_prev_se_recent, 
                                          ymax = mid_week_prev_recent + mid_week_prev_se_recent)) +
    geom_xsidecol(data = seq_numbers, aes(x = ymd(collection_date), y = seq, fill = method),
                  show.legend = F, alpha = 1)+
    scale_xsidey_continuous(position = 'right',
                            minor_breaks = NULL,
                            n.breaks = 3,
                            guide = guide_axis(TeX(r"(Number of collected sequences (log$_{10}$))")))+
    facet_grid(method~reference_date)+
    theme_bw()+ 
    guides(colour = guide_legend(override.aes = list(alpha = 1,
                                                     size = 4),
                                 nrow = 1))+
    theme(ggside.panel.scale = .3,
          axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle =- 45, hjust=.6, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          strip.text = element_text(size = strip_text_size),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size=15),
          legend.position = "bottom",
          strip.background =element_rect(fill="white"))+
    scale_color_manual(values = model_colors)+
    scale_fill_manual(values = model_colors)+
    scale_x_date(date_labels = "%B %d")+
    #coord_cartesian(ylim = c(0,1))+
    labs(x = '',
         y = 'BA.5 prevalence',
         fill = '',
         color = 'Method',
         tag = 'E')+
    ggtitle(paste0(this_country, ' BA.5'))



Brier_score_comparison3<- ggplot()+
    geom_point(data = this_country_metrics, aes(x = as.factor(num), y = BS_mean, 
                                                color = model), size = 3.5, position = dodge, show.legend = F) +
    geom_linerange(data = this_country_metrics, aes(x = as.factor(num), ymin = BS_lb,
                                                    ymax = BS_ub, color = model),size = 1.5, position = dodge, show.legend = F) +
    theme_bw()+ 
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    #scale_y_sqrt()+
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2, size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          legend.position = "none")+
    scale_color_manual(values = model_colors)+
    #scale_x_date(date_labels = "%B %d", limits = c(ymd('2022-04-25'), ymd('2022-07-15')))+
    labs(x = 'Reference date',
         y = 'Brier score',
         fill = '',
         color = 'Model',
         tag = 'F') + 
    ggtitle(paste0(this_country, ' BA.5'))


this_country <-list_of_countries[4]
combined_df <- get_single_country_metrics(this_country)
calib_df<- combined_df %>% filter(period == 'calibration' )
forecast_df <- combined_df %>% filter(period == 'forecast')
this_country_metrics <-cleaned_metrics %>% filter(country == this_country)


timeseries4 <- ggplot() +
    geom_line(data = calib_df, aes(x = ymd(collection_date), y = p_hat_mean, color = method), size =1) +
    geom_ribbon(data = calib_df, aes(x = ymd(collection_date), ymin = p_hat_lb, ymax = p_hat_ub, fill = method), alpha = 0.4, show.legend = F)+
    geom_line(data = forecast_df, aes(x = ymd(collection_date), y = p_hat_mean, color = method), size =0.2, show.legend = F) +
    geom_ribbon(data = forecast_df, aes(x = ymd(collection_date), ymin = p_hat_lb, ymax = p_hat_ub, fill = method), alpha = 0.2, show.legend = F)+
    geom_point(data = country_df, aes(x = ymd(mid_week_date_recent), y = mid_week_prev_recent), size = 1) +
    geom_linerange(data = country_df, aes(x = ymd(mid_week_date_recent), ymin = mid_week_prev_recent -mid_week_prev_se_recent, 
                                          ymax = mid_week_prev_recent + mid_week_prev_se_recent)) +
    geom_xsidecol(data = seq_numbers, aes(x = ymd(collection_date), y = seq, fill = method),
                  show.legend = F, alpha = 1)+
    scale_xsidey_continuous(position = 'right',
                            minor_breaks = NULL,
                            n.breaks = 3,
                            guide = guide_axis(TeX(r"(Number of collected sequences (log$_{10}$))")))+
    facet_grid(method~reference_date)+
    theme_bw()+ 
    guides(colour = guide_legend(override.aes = list(alpha = 1,
                                                     size = 4),
                                 nrow = 1))+
    theme(ggside.panel.scale = .3,
          axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle =- 45, hjust=.6, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          strip.text = element_text(size = strip_text_size),
          legend.title = element_text(size = 15), 
          legend.text = element_text(size=15),
          legend.position = "bottom",
          strip.background =element_rect(fill="white"))+
    scale_color_manual(values = model_colors)+
    scale_fill_manual(values = model_colors)+
    scale_x_date(date_labels = "%B %d")+
    #coord_cartesian(ylim = c(0,1))+
    labs(x = '',
         y = 'BA.5 prevalence',
         fill = '',
         color = 'Method',
         tag = 'G')+
    ggtitle(paste0(this_country, ' BA.5'))



Brier_score_comparison4<- ggplot()+
    geom_point(data = this_country_metrics, aes(x = as.factor(num), y = BS_mean, 
                                                color = model), size = 3.5, position = dodge, show.legend = F) +
    geom_linerange(data = this_country_metrics, aes(x = as.factor(num), ymin = BS_lb,
                                                    ymax = BS_ub, color = model),size = 1.5, position = dodge, show.legend = F) +
    theme_bw()+ 
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    #scale_y_sqrt()+
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2, size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          legend.position = "none")+
    scale_color_manual(values = model_colors)+
    #scale_x_date(date_labels = "%B %d", limits = c(ymd('2022-04-25'), ymd('2022-07-15')))+
    labs(x = 'Reference date',
         y = 'Brier score',
         fill = '',
         color = 'Model',
         tag = 'H') + 
    ggtitle(paste0(this_country, ' BA.5'))


design <- "
  111122
  333344
  555566
  777788
"

fig <- timeseries + Brier_score_comparison +
    timeseries2 + Brier_score_comparison2+
    timeseries3 + Brier_score_comparison3 + 
    timeseries4 + Brier_score_comparison4+
    plot_layout(design = design)
fig

ggsave('../data/output/figures/multiple_countries_retro_eval.pdf', 
       plot = fig,
       width = 15,
       height = 21)




# Make a figure with sequence counts by country-lineage-timepoint
table<- r_summary %>% filter(lineage %in% lineages, country %in% countries) %>% 
    select(country, lineage, n_seq_lineage_country, r_MLE_mean) %>% filter(is.na(r_MLE_mean))

# Quick plot on number of sequences of a lineage in a country to compare
fig_seq_nums<- r_summary %>% filter(lineage %in% lineages, country %in% countries) %>% 
    ggplot() + geom_tile(aes(x = as.factor(num), y = lineage, fill = log10(n_seq_lineage_country))) +
    facet_wrap(~country, ncol = 3)+
    labs(x = 'Reference date',
         y = 'Lineage',
         fill = 'log10(N)')+
    # scale_fill_gradientn(colours = pal2,
    #                      na.value = "lightgray", guide = "colourbar", aesthetics = "fill",
    #                      breaks=c(0, 5,100,1000,2000), labels=c(0,5,100,1000,2000))+
    theme_bw()+
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    theme(aspect.ratio = 1,
          axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=0, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          #axis.text.y = element_blank(),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = 0.5*tag_size),
          strip.background =element_rect(fill="white"))+
    ggtitle('Number of sequences of each variant in each country over time')
fig_seq_nums

ggsave('../data/output/figures/Supp_figure_seq_counts.pdf', 
       plot = fig_seq_nums,
       width = 8,
       height = 8)











# Supplementary figure comparing Portugal BA.5 estimates multicountry vs MLE

sfig4MCPortugalBA5 <- r_comb %>% 
    filter(country == this_country, lineage == 'BA.5') %>% 
    ggplot(aes(x =as.factor(num), y = transmission_advantage), fill = model_colors[1], color = model_colors[1], show.legend = F) + 
    stat_halfeye(slab_alpha = 0.5, point_size = 1, fill = model_colors[1], color = model_colors[1])+
    #geom_density_ridges(scale = 1, fill = model_colors[1], alpha = 0.5, show.legend = F)+
    #geom_point(aes(x =trans_adv_median, y = reference_date), color = model_colors[1], show.legend = F)+
    #geom_errorbarh(aes(xmin = trans_adv_lb, xmax = trans_adv_ub, y = reference_date), width = 0.3, color = model_colors[1], show.legend = F)+
    #geom_point(aes(x =trans_adv_MLE_mean, y = reference_date), color = model_colors[2], show.legend = F)+
    #geom_errorbarh(aes(xmin = trans_adv_MLE_lb, xmax = trans_adv_MLE_ub, y = reference_date), width = 0.5, color = model_colors[2], show.legend = F)+
    theme_bw()+
    ylim(c(0.5, 8.5))+
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    #scale_x_log10()+
    #coord_cartesian(xlim = c(0, 5))+
    #xlim(c(0.5,8.5))+
    #scale_color_manual(values = model_colors)+
    #scale_fill_manual(values = model_colors)+
    theme(axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2,size = 0.7*axis_text_size),
          axis.title.y = element_text(size = 0.7*axis_title_size),
          axis.title.x = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    labs(y = paste0('Weekly fitness advantage of BA.5 relative to BA.2'),
         x = 'Reference date',
         fill = '',
         color = '',
         tag = 'A') + 
    ggtitle('Portugal BA.5 Multicountry model')
sfig4MCPortugalBA5

sfig4MLEPortugalBA5 <- r_comb %>% 
    filter(country == this_country, lineage == 'BA.5') %>%
    ggplot() + 
    geom_point(aes(y =trans_adv_MLE_mean, x = as.factor(num)), color = model_colors[2], show.legend = F)+
    geom_linerange(aes( x = as.factor(num), ymin = trans_adv_MLE_lb, ymax = trans_adv_MLE_ub), width = 0.5, color = model_colors[2], show.legend = F)+
    #geom_point(aes(x =trans_adv_median, y = reference_date), color = model_colors[1], show.legend = F)+
    #geom_errorbarh(aes(xmin = trans_adv_lb, xmax = trans_adv_ub, y = reference_date), width = 0.3, color = model_colors[1], show.legend = F)+
    #geom_point(aes(x =trans_adv_MLE_mean, y = reference_date), color = model_colors[2], show.legend = F)+
    #geom_errorbarh(aes(xmin = trans_adv_MLE_lb, xmax = trans_adv_MLE_ub, y = reference_date), width = 0.5, color = model_colors[2], show.legend = F)+
    theme_bw()+
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    ylim(c(0.5, 8.5))+
    #scale_y_log10()+
    #coord_cartesian(ylim = c(0.5, 3.7))+
    scale_color_manual(values = model_colors)+
    scale_fill_manual(values = model_colors)+
    theme(axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2,size = 0.7*axis_text_size),
          axis.title.y = element_text(size = 0.7*axis_title_size),
          axis.title.x = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    labs(x = 'Reference date',
         y = 'Weekly fitness advantage of BA.5 relative to BA.2',
         fill = '',
         tag = 'B') +
    ggtitle('Portugal BA.5 single country model')
sfig4MLEPortugalBA5

sfig4 <- ggarrange(sfig4MCPortugalBA5, sfig4MLEPortugalBA5, ncol = 2, nrow = 1, align = 'v')
sfig4

ggsave('../data/output/figures/SFIG2_BS.pdf', 
       plot = sfig4,
       width = 12,
       height = 6)

# Same thing but with MLE using normal approx
sfigMLEPortugalBA5_normapprox <- r_comb %>% 
    filter(country == this_country, lineage == 'BA.5') %>%
    ggplot() + 
    geom_point(aes(y =trans_adv_MLE_mean, x = as.factor(num)), color = model_colors[2], show.legend = F)+
    geom_linerange(aes( x = as.factor(num), ymin = trans_adv_normal_approx_lb, ymax = trans_adv_normal_approx_ub), width = 0.5, color = model_colors[2], show.legend = F)+
    #geom_point(aes(x =trans_adv_median, y = reference_date), color = model_colors[1], show.legend = F)+
    #geom_errorbarh(aes(xmin = trans_adv_lb, xmax = trans_adv_ub, y = reference_date), width = 0.3, color = model_colors[1], show.legend = F)+
    #geom_point(aes(x =trans_adv_MLE_mean, y = reference_date), color = model_colors[2], show.legend = F)+
    #geom_errorbarh(aes(xmin = trans_adv_MLE_lb, xmax = trans_adv_MLE_ub, y = reference_date), width = 0.5, color = model_colors[2], show.legend = F)+
    theme_bw()+
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    ylim(c(0.5, 5.2))+
    #scale_y_log10()+
    #coord_cartesian(ylim = c(0.5, 3.7))+
    scale_color_manual(values = model_colors)+
    scale_fill_manual(values = model_colors)+
    theme(axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2,size = 0.7*axis_text_size),
          axis.title.y = element_text(size = 0.7*axis_title_size),
          axis.title.x = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    labs(x = 'Reference date',
         y = 'Weekly fitness advantage of BA.5 relative to BA.2',
         fill = '',
         tag = 'B') +
    ggtitle('Single country model')
sfigMLEPortugalBA5_normapprox


sfigMCPortugalBA5_normapprox <- r_comb %>% 
    filter(country == this_country, lineage == 'BA.5') %>% 
    ggplot(aes(x =as.factor(num), y = transmission_advantage), fill = model_colors[1], color = model_colors[1], show.legend = F) + 
    stat_halfeye(slab_alpha = 0.5, point_size = 1, fill = model_colors[1], color = model_colors[1])+
    #geom_density_ridges(scale = 1, fill = model_colors[1], alpha = 0.5, show.legend = F)+
    #geom_point(aes(x =trans_adv_median, y = reference_date), color = model_colors[1], show.legend = F)+
    #geom_errorbarh(aes(xmin = trans_adv_lb, xmax = trans_adv_ub, y = reference_date), width = 0.3, color = model_colors[1], show.legend = F)+
    #geom_point(aes(x =trans_adv_MLE_mean, y = reference_date), color = model_colors[2], show.legend = F)+
    #geom_errorbarh(aes(xmin = trans_adv_MLE_lb, xmax = trans_adv_MLE_ub, y = reference_date), width = 0.5, color = model_colors[2], show.legend = F)+
    theme_bw()+
    ylim(c(0.5, 5.2))+
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    #scale_x_log10()+
    #coord_cartesian(xlim = c(0, 5))+
    #xlim(c(0.5,8.5))+
    #scale_color_manual(values = model_colors)+
    #scale_fill_manual(values = model_colors)+
    theme(axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2,size = 0.7*axis_text_size),
          axis.title.y = element_text(size = 0.7*axis_title_size),
          axis.title.x = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    labs(y = paste0('Weekly fitness advantage of BA.5 relative to BA.2'),
         #y = 'Reference date',
         x = 'Reference date',
         fill = '',
         color = '',
         tag = 'A') + 
    ggtitle('Multicountry model')
sfigMCPortugalBA5_normapprox

sfig4_norm_approx <- ggarrange(sfigMCPortugalBA5_normapprox, sfigMLEPortugalBA5_normapprox,  ncol = 2, nrow = 1, align = 'v')
sfig4_norm_approx


ggsave('../data/output/figures/SFIG2.png', 
       plot = sfig4_norm_approx,
       width = 12,
       height = 6)




# Supplementary Figure: Brier scores for calibration and forecast periods

dodge<- position_dodge(width = 0.2)

calib_metrics<- cleaned_metrics %>% filter(country == this_country, multicountry_period == 'multicountry calibration')
forecast_metrics<- cleaned_metrics %>% filter(country == this_country, multicountry_period == 'multicountry forecast') 
sfigBScalib<- ggplot()+
    geom_point(data = calib_metrics, aes(x = as.factor(num), y = BS_mean_period, 
                                         color = model), position = dodge, show.legend = F) +
    geom_linerange(data = calib_metrics, aes(x = as.factor(num), ymin = BS_lb_period,
                                             ymax = BS_ub_period, color = model), position = dodge, show.legend = F) +
    theme_bw()+ 
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    #scale_y_sqrt()+
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2, size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          legend.position = "none")+
    scale_color_manual(values = model_colors)+
    #scale_x_date(date_labels = "%B %d", limits = c(ymd('2022-04-25'), ymd('2022-07-15')))+
    labs(x = 'Reference date',
         y = 'Brier score',
         fill = '',
         color = 'Model',
         tag = 'C') + 
    ggtitle('Portugal calibration period')
sfigBScalib

sfigBSforecast<- ggplot()+
    geom_point(data = forecast_metrics, aes(x = as.factor(num), y = BS_mean_period, 
                                            color = model), position = dodge, show.legend = F) +
    geom_linerange(data = forecast_metrics, aes(x = as.factor(num), ymin = BS_lb_period,
                                                ymax = BS_ub_period, color = model), position = dodge, show.legend = F) +
    theme_bw()+ 
    #scale_y_sqrt()+
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2, size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          legend.position = 'none')+
    scale_color_manual(values = model_colors_muted)+
    # scale_x_date(date_labels = "%B %d", limits = c(ymd('2022-04-25'), ymd('2022-07-15')))+
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    labs(x = 'Reference date',
         y = 'Brier score',
         fill = '',
         color = 'Model',
         tag = 'A')+
    ggtitle('Portugal forecast period')
sfigBSforecast
sfigBSperiod<- ggarrange(sfigBScalib, sfigBSforecast, ncol = 2, nrow = 1,
                    align = "v")
sfigBSperiod

# Not currently a part of the manuscript but we could add in 
ggsave('../data/output/figures/sfigure_Portugal_BA5_BS.png', 
       plot = sfigBSperiod,
       width = 12,
       height = 6)

# Supplementary Figure: Heatmaps of r for multicountry and MLE
sfigrheatmapMC<- r_summary %>% filter(lineage %in% lineages, country %in% countries) %>% 
    mutate(time = rank(reference_date)) %>% 
    ggplot() + geom_tile(aes(x = as.factor(num), y = lineage, fill = trans_adv_median)) +
    facet_wrap(~country, ncol = 3) +
    labs(x = 'Reference date',
         y = 'Lineage',
         fill = 'transmission\nadvantage\nestimate',
         tag = 'A')+
    #scale_fill_continuous(type = "viridis")+
    # scale_fill_gradientn(colours = pal2,
    #                      na.value = "gray", guide = "colourbar", aesthetics = "fill",
    #                      breaks=c(-0.5, 0,1,2,3, 4),labels=c(-0.5, 0,1,2,3, 4),
    #                      limits=c(-1, 4))+
    scale_colour_gradient2(low = "red", mid = "white", high = "blue", midpoint = 1,
                           na.value = "gray", guide = "colourbar", aesthetics = "fill",
                           breaks=c(-0.5, 0,1, 2, 3, 4),labels=c(-0.5, 0,1,2,3,4),
                           limits=c(-1, 4))+
    theme_bw()+
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    # coord_fixed()+
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=0, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          strip.background =element_rect(fill="white"))+
    ggtitle('Multicountry model')
sfigrheatmapMC
sfigrheatmapMLE<- r_summary %>% filter(lineage %in% lineages, country %in% countries) %>% 
    mutate(time = rank(reference_date)) %>% 
    ggplot() + geom_tile(aes(x = as.factor(num), y = lineage, fill = trans_adv_MLE_mean)) +
    facet_wrap(~country, ncol = 3)+
    labs(x = 'Reference date',
         y = 'Lineage',
         fill = 'transmission\nadvantage\nestimate',
         tag = 'B')+
    scale_colour_gradient2(low = "red", mid = "white", high = "blue", midpoint = 1,
                           na.value = "gray", guide = "colourbar", aesthetics = "fill",
                           breaks=c(-0.5, 0,1, 2, 3, 4),labels=c(-0.5, 0,1,2,3,4),
                           limits=c(-1, 4))+
    theme_bw()+
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=0, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          strip.background =element_rect(fill="white"))+
    ggtitle('Single country model')
sfigrheatmapMLE
design <- "
  12
"
fig <- sfigrheatmapMC + sfigrheatmapMLE +plot_layout(design = design)

ggsave('../data/output/figures/SFIG3.pdf', 
       plot = fig,
       width = 10,
       height = 5)


# Supplementary Figure: Calibration data vs model for MLE vs MC

sfig4a <- ggplot() +
    geom_line(data = calib_df, aes(x = ymd(collection_date), y = p_hat_mean, color = method), size =1) +
    geom_ribbon(data = calib_df, aes(x = ymd(collection_date), ymin = p_hat_lb, ymax = p_hat_ub, fill = method), alpha = 0.4, show.legend = F)+
    geom_line(data = forecast_df, aes(x = ymd(collection_date), y = p_hat_mean, color = method), size =0.2, show.legend = F) +
    geom_ribbon(data = forecast_df, aes(x = ymd(collection_date), ymin = p_hat_lb, ymax = p_hat_ub, fill = method), alpha = 0.2, show.legend = F)+
    geom_point(data = country_df, aes(x = ymd(mid_week_date), y = mid_week_p_lineage), size = 1, color = 'purple') +
    geom_linerange(data = country_df, aes(x = ymd(mid_week_date), ymin = mid_week_p_lineage -mid_week_p_lineage_se, 
                                          ymax = mid_week_p_lineage + mid_week_p_lineage_se)) +
    facet_grid(method~reference_date)+
    theme_bw()+ 
    guides(colour = guide_legend(override.aes = list(alpha = 1,
                                                     size = 4),
                                 nrow = 1))+
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle =- 45, hjust=.6, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          legend.position = "bottom",
          strip.background =element_rect(fill="white"))+
    scale_color_manual(values = model_colors)+
    scale_fill_manual(values = model_colors)+
    scale_x_date(date_labels = "%B %d")+
    coord_cartesian(ylim = c(0,1))+
    labs(x = '',
         y = 'BA.5 prevalence',
         fill = '',
         color = 'Method',
         tag = '')+
    ggtitle(paste0(this_country, ' BA.5 model compared to calibration data'))
sfig4a


# Save fig ----------------------------------------------------------------

ggsave('../data/output/figures/SFIG1.pdf', 
       plot = sfig4a,
       width = 15,
       height = 6)














# Test out a heat map ----------------------------------------------------

fig4addtop<- r_summary %>% filter(lineage %in% lineages, country %in% countries) %>% 
    mutate(time = rank(reference_date)) %>% 
    ggplot() + geom_tile(aes(x = as.factor(num), y = lineage, fill = trans_adv_median)) +
    facet_wrap(~country, ncol = 3) +
    labs(x = 'Reference date',
         y = 'Lineage',
         fill = 'transmission\nadvantage\nestimate',
         tag = 'G')+
    #scale_fill_continuous(type = "viridis")+
    scale_fill_gradientn(colours = pal2,
                            na.value = "gray", guide = "colourbar", aesthetics = "fill",
                           breaks=c(-0.5, 0,1,2,3, 4),labels=c(-0.5, 0,1,2,3, 4),
                         limits=c(-1, 4))+
    theme_bw()+
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
   # coord_fixed()+
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=0, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          strip.background =element_rect(fill="white"))+
    ggtitle('Multicountry model')
fig4addtop
fig4addbottom<- r_summary %>% filter(lineage %in% lineages, country %in% countries) %>% 
    mutate(time = rank(reference_date)) %>% 
    ggplot() + geom_tile(aes(x = as.factor(num), y = lineage, fill = trans_adv_MLE_mean)) +
    facet_wrap(~country, ncol = 3)+
    labs(x = 'Reference date',
         y = 'Lineage',
         fill = 'transmission\nadvantage\nestimate',
         tag = 'H')+
    scale_colour_gradient2(low = "red", mid = "white", high = "blue", midpoint = 1,
                           na.value = "gray", guide = "colourbar", aesthetics = "fill",
                           breaks=c(-0.5, 0,1, 2, 3, 4),labels=c(-0.5, 0,1,2,3,4),
                           limits=c(-1, 4))+
    theme_bw()+
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=0, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          strip.background =element_rect(fill="white"))+
    ggtitle('MLE model')
fig4addbottom
design <- "
  1111
  2222
"
fig <- fig4addtop + fig4addbottom +plot_layout(design = design)

ggsave('../data/output/figures/fig4_heatmap.png', 
       plot = fig,
       width = 10,
       height = 15)


# Try it with the difference from final
r_summary <- r_summary %>% mutate( dif_trans_capped_MC = ifelse(dif_trans_adv_MC >4, 4, dif_trans_adv_MC),
                                   dif_trans_capped_MLE = ifelse(dif_trans_adv_MLE >4, 4, dif_trans_adv_MLE))
fig4addtop<- r_summary %>% filter(lineage %in% lineages, country %in% countries) %>% 
    ggplot() + geom_tile(aes(x = as.factor(num), y = lineage, fill = dif_trans_capped_MC)) +
    facet_wrap(~country, ncol = 3) +
    labs(x = 'Reference date',
         y = 'Lineage',
         fill = 'change in\nfitness\nadvantage',
         tag = 'F')+
    #scale_fill_continuous(type = "viridis")+
    scale_fill_gradientn(colours = pal2,
                           na.value = "gray", guide = "colourbar", aesthetics = "fill",
                           #breaks=c(-0.5, 0,1, 2,3, 4), labels=c(-1, 0, 1,2,3, 4),
                           limits=c(-0.5, 4))+
    # scale_colour_gradient2(low = pal2[1], mid = pal2[10], high = pal2[21], midpoint = 0,
    #                        na.value = "gray", guide = "colourbar", aesthetics = "fill",
    #                        breaks=c(-0.5, 0,1, 2, 3, 4),labels=c(-0.5, 0,1,2,3,4),
    #                        limits=c(-1, 4))+
    theme_bw()+
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    # coord_fixed()+
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=0, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          strip.background =element_rect(fill="white"))+
    ggtitle('Multicountry model')
fig4addtop
fig4addbottom<- r_summary %>% filter(lineage %in% lineages, country %in% countries) %>% 
    ggplot() + geom_tile(aes(x = as.factor(num), y = lineage, fill = dif_trans_capped_MLE)) +
    facet_wrap(~country, ncol = 3)+
    labs(x = 'Reference date',
         y = 'Lineage',
         fill = 'change in\nfitness\nadvantage',
         tag = 'G')+
    scale_fill_gradientn(colours = pal2,
                           na.value = "gray", guide = "colourbar", aesthetics = "fill",
                           #breaks=c(-0.5, 0,1,3, 2, 4), labels=c(-0.5, 0, 1,3, 2, 4),
                           limits=c(-0.5, 4))+
    theme_bw()+
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=0, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          strip.background =element_rect(fill="white"))+
    ggtitle('MLE model')
fig4addbottom
design <- "
  1111
  2222
"
fig <- fig4addtop + fig4addbottom +plot_layout(design = design)

ggsave('../data/output/figures/fig4_heatmap_dif.png', 
       plot = fig,
       width = 10,
       height = 15)

# Try combining everything together------------------------------
design <- "
  111122
  111122
  334455
"

fig <- fig4a + fig4bc + fig4dejoint + fig4addtop + fig4addbottom+
    plot_layout(design = design)
fig

# Save fig ----------------------------------------------------------------

ggsave('../data/output/figures/figure_4_Portugal_BA5_big.png', 
       plot = fig,
       width = 21,
       height = 20)

# Remove C and D --------------------------------------------------------


fig4addtop<- r_summary %>% filter(lineage %in% lineages, country %in% countries) %>% 
    ggplot() + geom_tile(aes(x = as.factor(num), y = lineage, fill = dif_trans_capped_MC)) +
    facet_wrap(~country, ncol = 3) +
    labs(x = 'Reference date',
         y = 'Lineage',
         fill = 'change in\nfitness\nadvantage',
         tag = 'D')+
    #scale_fill_continuous(type = "viridis")+
    scale_fill_gradientn(colours = pal2,
                         na.value = "gray", guide = "colourbar", aesthetics = "fill",
                         #breaks=c(-0.5, 0,1, 2,3, 4), labels=c(-1, 0, 1,2,3, 4),
                         limits=c(-0.5, 4))+
    # scale_colour_gradient2(low = pal2[1], mid = pal2[10], high = pal2[21], midpoint = 0,
    #                        na.value = "gray", guide = "colourbar", aesthetics = "fill",
    #                        breaks=c(-0.5, 0,1, 2, 3, 4),labels=c(-0.5, 0,1,2,3,4),
    #                        limits=c(-1, 4))+
    theme_bw()+
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    # coord_fixed()+
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=0, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          strip.background =element_rect(fill="white"))+
    ggtitle('Multicountry model')
fig4addtop
fig4addbottom<- r_summary %>% filter(lineage %in% lineages, country %in% countries) %>% 
    ggplot() + geom_tile(aes(x = as.factor(num), y = lineage, fill = dif_trans_capped_MLE)) +
    facet_wrap(~country, ncol = 3)+
    labs(x = 'Reference date',
         y = 'Lineage',
         fill = 'change in\nfitness\nadvantage',
         tag = 'E')+
    scale_fill_gradientn(colours = pal2,
                         na.value = "gray", guide = "colourbar", aesthetics = "fill",
                         #breaks=c(-0.5, 0,1,3, 2, 4), labels=c(-0.5, 0, 1,3, 2, 4),
                         limits=c(-0.5, 4))+
    theme_bw()+
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=0, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size),
          strip.background =element_rect(fill="white"))+
    ggtitle('MLE model')
fig4addbottom

fig4global <- mu_comb %>% 
    filter(lineage == 'BA.5') %>% 
    ggplot(aes( x = as.factor(num), y = global_transmission_advantage))+
    stat_halfeye(slab_alpha = 0.7,
                 #interval_color = 'gray',
                 #point_color = 'gray',
                 aes(fill = lineage),
                 show.legend = F)+
    theme_bw()+
    scale_x_discrete(labels = c("1" = "April 30", "2" = "May 16", "3" = "May 27", "4" = "June 04", "5" = "June 27"))+
    scale_fill_manual(values = lineage_colors)+
    theme(axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2,size = 0.7*axis_text_size),
          axis.title.y = element_text(size = 0.7*axis_title_size),
          axis.title.x = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    labs(y = paste0('Global weekly fitness advantage of BA.5 relative to BA.2'),
         x = 'Reference date',
         fill = '',
         tag = 'B')
fig4global



design <- "
  112
  345
  
"

fig <- fig4a +  fig4global + fig4efjointc + fig4addtop + fig4addbottom +
    plot_layout(design = design)
fig

ggsave('../data/output/figures/figure_4_Portugal_BA5_w_global.png', 
       plot = fig,
       width = 21,
       height = 15)




# Spain BA.5 option --------------------------------------------------------------------


this_country<- 'Spain'
country_df<-clean_global_df %>% filter(country == this_country, lineage == "BA.5")
calib_df<- country_df %>% filter(multicountry_period == 'multicountry calibration')
forecast_df <- country_df %>% filter(multicountry_period == 'multicountry forecast')
fig4a <- ggplot() +
    geom_line(data = calib_df, aes(x = ymd(collection_date), y = p_hat_mean), color = model_colors[1], size =1, show.legend = F) +
    geom_ribbon(data = calib_df, aes(x = ymd(collection_date), ymin = p_hat_lb, ymax = p_hat_ub), fill = model_colors[1], alpha = 0.6, show.legend = F)+
    geom_line(data = forecast_df, aes(x = ymd(collection_date), y = p_hat_mean), color = model_colors[1], size =0.2, show.legend = F) +
    geom_ribbon(data = forecast_df, aes(x = ymd(collection_date), ymin = p_hat_lb, ymax = p_hat_ub), fill = model_colors[1], alpha = 0.2, show.legend = F)+
    geom_point(data = country_df, aes(x = ymd(mid_week_date_recent), y = mid_week_prev_recent), size = 0.5) +
    geom_linerange(data = country_df, aes(x = ymd(mid_week_date_recent), ymin = mid_week_prev_recent -mid_week_prev_se_recent, 
                                          ymax = mid_week_prev_recent + mid_week_prev_se_recent)) +
    facet_wrap(~reference_date, ncol = 5)+
    theme_bw()+ 
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle =- 45, hjust=.6, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    scale_color_manual(values = lineage_colors)+
    scale_fill_manual(values = lineage_colors)+
    scale_x_date(date_labels = "%B %d")+
    coord_cartesian(ylim = c(0,1))+
    labs(x = '',
         y = 'BA.5 prevalence',
         fill = '',
         tag = 'A')+
    ggtitle(paste0(this_country, ' BA.5 multicountry model'))
fig4a

fig4b <- ggplot() + 
    geom_line(data = calib_df, aes(x = ymd(collection_date), y = p_hat_MLE_mean),color = model_colors[2], size =1, show.legend = F) +
    geom_ribbon(data = calib_df, aes(x = ymd(collection_date), ymin = p_hat_MLE_lb, ymax = p_hat_MLE_ub), fill = model_colors[2], alpha = 0.6, show.legend = F)+
    geom_line(data = forecast_df, aes(x = ymd(collection_date), y = p_hat_MLE_mean), color = model_colors[2], size =0.2, show.legend = F) +
    geom_ribbon(data = forecast_df, aes(x = ymd(collection_date), ymin = p_hat_MLE_lb, ymax = p_hat_MLE_ub), fill = model_colors[2], alpha = 0.2, show.legend = F)+
    geom_point(data = country_df, aes(x = ymd(mid_week_date_recent), y = mid_week_prev_recent), size = 0.5) +
    geom_linerange(data = country_df, aes(x = ymd(mid_week_date_recent), ymin = mid_week_prev_recent -mid_week_prev_se_recent, 
                                          ymax = mid_week_prev_recent + mid_week_prev_se_recent)) +
    facet_wrap(~reference_date, ncol = 5)+
    theme_bw()+ 
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=.6, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    scale_color_manual(values = lineage_colors)+
    scale_fill_manual(values = lineage_colors)+
    scale_x_date(date_labels = "%B %d")+
    coord_cartesian(ylim = c(0,1))+
    labs(x = '',
         y = 'BA.5 prevalence',
         fill = '',
         tag = 'B')+
    ggtitle(paste0(this_country, ' BA.5 MLE model'))
fig4b

cleaned_metrics_MC <- country_metrics %>% 
    select(country, reference_date, multicountry_period,
           BS_med_MC, BS_mean_MC, BS_lb_MC, BS_ub_MC) %>% 
    rename(BS_med = BS_med_MC,
           BS_mean = BS_mean_MC,
           BS_lb = BS_lb_MC,
           BS_ub = BS_ub_MC) %>% 
    mutate( model = 'multicountry')
cleaned_metrics_MLE <- country_metrics %>% 
    select(country, reference_date, multicountry_period,
           BS_med_MLE, BS_mean_MLE, BS_lb_MLE, BS_ub_MLE) %>% 
    rename(BS_med = BS_med_MLE,
           BS_mean = BS_mean_MLE,
           BS_lb = BS_lb_MLE,
           BS_ub = BS_ub_MLE) %>% 
    mutate( model = 'MLE')
cleaned_metrics <- bind_rows(cleaned_metrics_MC, cleaned_metrics_MLE)

dodge<- position_dodge(width = 2)
fig4d<- cleaned_metrics %>% filter(country == this_country) %>% 
    ggplot() + geom_point(aes(x = reference_date, y = BS_mean, color = model), position = dodge) +
    geom_errorbar(aes(x = reference_date, ymin = BS_lb, ymax = BS_ub, color = model), position = dodge, width = 0.8) +
    theme_bw()+ 
    scale_y_sqrt()+
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2, size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    scale_color_manual(values = model_colors)+
    scale_x_date(date_labels = "%B %d")+
    labs(x = 'Reference date',
         y = 'Brier score',
         fill = '',
         color = 'Model',
         tag = 'D')
fig4d


fig4c <- r_comb %>% 
    filter(country == this_country, lineage == 'BA.5') %>% 
    ggplot(aes(x = transmission_advantage, y = reference_date, 
               group = ref_date_num)) + 
    coord_flip()+
    geom_density_ridges(scale = 1, fill = model_colors[1], alpha = 0.5, show.legend = F)+
    geom_point(aes(x =trans_adv_median, y = reference_date), color = model_colors[1], show.legend = F)+
    #geom_errorbarh(aes(xmin = trans_adv_lb, xmax = trans_adv_ub, y = reference_date), width = 0.3, color = model_colors[1], show.legend = F)+
    geom_point(aes(x =trans_adv_MLE_mean, y = reference_date), color = model_colors[2], show.legend = F)+
    geom_errorbarh(aes(xmin = trans_adv_MLE_lb, xmax = trans_adv_MLE_ub, y = reference_date), width = 0.5, color = model_colors[2], show.legend = F)+
    theme_bw()+
    scale_y_date(date_labels = "%B %d")+
    scale_x_log10()+
    #coord_cartesian(xlim =c(0,100))+
    scale_color_manual(values = model_colors)+
    scale_fill_manual(values = model_colors)+
    theme(axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2,size = 0.7*axis_text_size),
          axis.title.y = element_text(size = 0.7*axis_title_size),
          axis.title.x = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    labs(x = paste0('Weekly fitness advantage of BA.5 relative to BA.2'),
         y = 'Reference date',
         fill = '',
         tag = 'C')
fig4c


# Put together panels 

design <- "
  11113
  22224
"

fig <- fig4a + fig4b + fig4c + fig4d + plot_layout(design = design)
fig

# Save fig 

ggsave('../data/output/figures/figure_4_Spain_BA5.png', 
       plot = fig,
       width = 15,
       height = 8)


# Portugal BA.2.12.1-------------------------------------------------

this_country<- 'Portugal'
country_df<-clean_global_df %>% filter(country == this_country, lineage == "BA.2.12.1")
calib_df<- country_df %>% filter(multicountry_period == 'multicountry calibration')
forecast_df <- country_df %>% filter(multicountry_period == 'multicountry forecast')
fig4a <- ggplot() +
    geom_line(data = calib_df, aes(x = ymd(collection_date), y = p_hat_mean), color = model_colors[1], size =1, show.legend = F) +
    geom_ribbon(data = calib_df, aes(x = ymd(collection_date), ymin = p_hat_lb, ymax = p_hat_ub), fill = model_colors[1], alpha = 0.6, show.legend = F)+
    geom_line(data = forecast_df, aes(x = ymd(collection_date), y = p_hat_mean), color = model_colors[1], size =0.2, show.legend = F) +
    geom_ribbon(data = forecast_df, aes(x = ymd(collection_date), ymin = p_hat_lb, ymax = p_hat_ub), fill = model_colors[1], alpha = 0.2, show.legend = F)+
    geom_point(data = country_df, aes(x = ymd(mid_week_date_recent), y = mid_week_prev_recent), size = 0.5) +
    geom_linerange(data = country_df, aes(x = ymd(mid_week_date_recent), ymin = mid_week_prev_recent -mid_week_prev_se_recent, 
                                          ymax = mid_week_prev_recent + mid_week_prev_se_recent)) +
    facet_wrap(~reference_date, ncol = 5)+
    theme_bw()+ 
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle =- 45, hjust=.6, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    scale_color_manual(values = lineage_colors)+
    scale_fill_manual(values = lineage_colors)+
    scale_x_date(date_labels = "%B %d")+
    coord_cartesian(ylim = c(0,1))+
    labs(x = '',
         y = 'BA.2.12.1 prevalence',
         fill = '',
         tag = 'A')+
    ggtitle(paste0(this_country, ' BA.2.12.1 multicountry model'))
fig4a

fig4b <- ggplot() + 
    geom_line(data = calib_df, aes(x = ymd(collection_date), y = p_hat_MLE_mean),color = model_colors[2], size =1, show.legend = F) +
    geom_ribbon(data = calib_df, aes(x = ymd(collection_date), ymin = p_hat_MLE_lb, ymax = p_hat_MLE_ub), fill = model_colors[2], alpha = 0.6, show.legend = F)+
    geom_line(data = forecast_df, aes(x = ymd(collection_date), y = p_hat_MLE_mean), color = model_colors[2], size =0.2, show.legend = F) +
    geom_ribbon(data = forecast_df, aes(x = ymd(collection_date), ymin = p_hat_MLE_lb, ymax = p_hat_MLE_ub), fill = model_colors[2], alpha = 0.2, show.legend = F)+
    geom_point(data = country_df, aes(x = ymd(mid_week_date_recent), y = mid_week_prev_recent), size = 0.5) +
    geom_linerange(data = country_df, aes(x = ymd(mid_week_date_recent), ymin = mid_week_prev_recent -mid_week_prev_se_recent, 
                                          ymax = mid_week_prev_recent + mid_week_prev_se_recent)) +
    facet_wrap(~reference_date, ncol = 5)+
    theme_bw()+ 
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=.6, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    scale_color_manual(values = lineage_colors)+
    scale_fill_manual(values = lineage_colors)+
    scale_x_date(date_labels = "%B %d")+
    coord_cartesian(ylim = c(0,1))+
    labs(x = '',
         y = 'BA.2.12.1 prevalence',
         fill = '',
         tag = 'B')+
    ggtitle(paste0(this_country, ' BA.2.12.1 MLE model'))
fig4b

cleaned_metrics_MC <- country_metrics %>% 
    select(country, reference_date, multicountry_period,
           BS_med_MC, BS_mean_MC, BS_lb_MC, BS_ub_MC) %>% 
    rename(BS_med = BS_med_MC,
           BS_mean = BS_mean_MC,
           BS_lb = BS_lb_MC,
           BS_ub = BS_ub_MC) %>% 
    mutate( model = 'multicountry')
cleaned_metrics_MLE <- country_metrics %>% 
    select(country, reference_date, multicountry_period,
           BS_med_MLE, BS_mean_MLE, BS_lb_MLE, BS_ub_MLE) %>% 
    rename(BS_med = BS_med_MLE,
           BS_mean = BS_mean_MLE,
           BS_lb = BS_lb_MLE,
           BS_ub = BS_ub_MLE) %>% 
    mutate( model = 'MLE')
cleaned_metrics <- bind_rows(cleaned_metrics_MC, cleaned_metrics_MLE)

dodge<- position_dodge(width = 2)
fig4d<- cleaned_metrics %>% filter(country == this_country) %>% 
    ggplot() + geom_point(aes(x = reference_date, y = BS_mean, color = model), position = dodge) +
    geom_errorbar(aes(x = reference_date, ymin = BS_lb, ymax = BS_ub, color = model), position = dodge, width = 0.8) +
    theme_bw()+ 
    scale_y_sqrt()+
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2, size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    scale_color_manual(values = model_colors)+
    scale_x_date(date_labels = "%B %d")+
    labs(x = 'Reference date',
         y = 'Brier score',
         fill = '',
         color = 'Model',
         tag = 'D')
fig4d


fig4c <- r_comb %>% 
    filter(country == this_country, lineage == 'BA.2.12.1') %>% 
    ggplot(aes(x = transmission_advantage, y = reference_date, 
               group = ref_date_num)) + 
    coord_flip()+
    geom_density_ridges(scale = 1, fill = model_colors[1], alpha = 0.5, show.legend = F)+
    geom_point(aes(x =trans_adv_median, y = reference_date), color = model_colors[1], show.legend = F)+
    #geom_errorbarh(aes(xmin = trans_adv_lb, xmax = trans_adv_ub, y = reference_date), width = 0.3, color = model_colors[1], show.legend = F)+
    geom_point(aes(x =trans_adv_MLE_mean, y = reference_date), color = model_colors[2], show.legend = F)+
    geom_errorbarh(aes(xmin = trans_adv_MLE_lb, xmax = trans_adv_MLE_ub, y = reference_date), width = 0.5, color = model_colors[2], show.legend = F)+
    theme_bw()+
    scale_y_date(date_labels = "%B %d")+
    scale_x_log10()+
    scale_color_manual(values = model_colors)+
    scale_fill_manual(values = model_colors)+
    theme(axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2,size = 0.7*axis_text_size),
          axis.title.y = element_text(size = 0.7*axis_title_size),
          axis.title.x = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    labs(x = paste0('Weekly fitness advantage of BA.2.12.1 relative to BA.2'),
         y = 'Reference date',
         fill = '',
         tag = 'C')
fig4c


# Put together panels -----------------------------------------------------

design <- "
  11113
  22224
"

fig <- fig4a + fig4b + fig4c + fig4d + plot_layout(design = design)
fig
# Save fig ----------------------------------------------------------------

ggsave('../data/output/figures/figure_4_Portugal_BA2121.png', 
       plot = fig,
       width = 15,
       height = 8)







# US BA.2.12.1-------------------------------------------------

this_country<- 'United States'
country_df<-clean_global_df %>% filter(country == this_country, lineage == "BA.2.12.1")
calib_df<- country_df %>% filter(multicountry_period == 'multicountry calibration')
forecast_df <- country_df %>% filter(multicountry_period == 'multicountry forecast')
fig4a <- ggplot() +
    geom_line(data = calib_df, aes(x = ymd(collection_date), y = p_hat_mean), color = model_colors[1], size =1, show.legend = F) +
    geom_ribbon(data = calib_df, aes(x = ymd(collection_date), ymin = p_hat_lb, ymax = p_hat_ub), fill = model_colors[1], alpha = 0.6, show.legend = F)+
    geom_line(data = forecast_df, aes(x = ymd(collection_date), y = p_hat_mean), color = model_colors[1], size =0.2, show.legend = F) +
    geom_ribbon(data = forecast_df, aes(x = ymd(collection_date), ymin = p_hat_lb, ymax = p_hat_ub), fill = model_colors[1], alpha = 0.2, show.legend = F)+
    geom_point(data = country_df, aes(x = ymd(mid_week_date_recent), y = mid_week_prev_recent), size = 0.5) +
    geom_linerange(data = country_df, aes(x = ymd(mid_week_date_recent), ymin = mid_week_prev_recent -mid_week_prev_se_recent, 
                                          ymax = mid_week_prev_recent + mid_week_prev_se_recent)) +
    facet_wrap(~reference_date, ncol = 5)+
    theme_bw()+ 
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle =- 45, hjust=.6, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    scale_color_manual(values = lineage_colors)+
    scale_fill_manual(values = lineage_colors)+
    scale_x_date(date_labels = "%B %d")+
    coord_cartesian(ylim = c(0,1))+
    labs(x = '',
         y = 'BA.2.12.1 prevalence',
         fill = '',
         tag = 'A')+
    ggtitle(paste0(this_country, ' BA.2.12.1 multicountry model'))
fig4a

fig4b <- ggplot() + 
    geom_line(data = calib_df, aes(x = ymd(collection_date), y = p_hat_MLE_mean),color = model_colors[2], size =1, show.legend = F) +
    geom_ribbon(data = calib_df, aes(x = ymd(collection_date), ymin = p_hat_MLE_lb, ymax = p_hat_MLE_ub), fill = model_colors[2], alpha = 0.6, show.legend = F)+
    geom_line(data = forecast_df, aes(x = ymd(collection_date), y = p_hat_MLE_mean), color = model_colors[2], size =0.2, show.legend = F) +
    geom_ribbon(data = forecast_df, aes(x = ymd(collection_date), ymin = p_hat_MLE_lb, ymax = p_hat_MLE_ub), fill = model_colors[2], alpha = 0.2, show.legend = F)+
    geom_point(data = country_df, aes(x = ymd(mid_week_date_recent), y = mid_week_prev_recent), size = 0.5) +
    geom_linerange(data = country_df, aes(x = ymd(mid_week_date_recent), ymin = mid_week_prev_recent -mid_week_prev_se_recent, 
                                          ymax = mid_week_prev_recent + mid_week_prev_se_recent)) +
    facet_wrap(~reference_date, ncol = 5)+
    theme_bw()+ 
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=.6, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    scale_color_manual(values = lineage_colors)+
    scale_fill_manual(values = lineage_colors)+
    scale_x_date(date_labels = "%B %d")+
    coord_cartesian(ylim = c(0,1))+
    labs(x = '',
         y = 'BA.2.12.1 prevalence',
         fill = '',
         tag = 'B')+
    ggtitle(paste0(this_country, ' BA.2.12.1 MLE model'))
fig4b

cleaned_metrics_MC <- country_metrics %>% 
    select(country, reference_date, multicountry_period,
           BS_med_MC, BS_mean_MC, BS_lb_MC, BS_ub_MC) %>% 
    rename(BS_med = BS_med_MC,
           BS_mean = BS_mean_MC,
           BS_lb = BS_lb_MC,
           BS_ub = BS_ub_MC) %>% 
    mutate( model = 'multicountry')
cleaned_metrics_MLE <- country_metrics %>% 
    select(country, reference_date, multicountry_period,
           BS_med_MLE, BS_mean_MLE, BS_lb_MLE, BS_ub_MLE) %>% 
    rename(BS_med = BS_med_MLE,
           BS_mean = BS_mean_MLE,
           BS_lb = BS_lb_MLE,
           BS_ub = BS_ub_MLE) %>% 
    mutate( model = 'MLE')
cleaned_metrics <- bind_rows(cleaned_metrics_MC, cleaned_metrics_MLE)

dodge<- position_dodge(width = 2)
fig4d<- cleaned_metrics %>% filter(country == this_country) %>% 
    ggplot() + geom_point(aes(x = reference_date, y = BS_mean, color = model), position = dodge) +
    geom_errorbar(aes(x = reference_date, ymin = BS_lb, ymax = BS_ub, color = model), position = dodge, width = 0.8) +
    theme_bw()+ 
    scale_y_sqrt()+
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2, size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    scale_color_manual(values = model_colors)+
    scale_x_date(date_labels = "%B %d")+
    labs(x = 'Reference date',
         y = 'Brier score',
         fill = '',
         color = 'Model',
         tag = 'D')
fig4d


fig4c <- r_comb %>% 
    filter(country == this_country, lineage == 'BA.2.12.1') %>% 
    ggplot(aes(x = transmission_advantage, y = reference_date, 
               group = ref_date_num)) + 
    coord_flip()+
    geom_density_ridges(scale = 1, fill = model_colors[1], alpha = 0.5, show.legend = F)+
    geom_point(aes(x =trans_adv_median, y = reference_date), color = model_colors[1], show.legend = F)+
    #geom_errorbarh(aes(xmin = trans_adv_lb, xmax = trans_adv_ub, y = reference_date), width = 0.3, color = model_colors[1], show.legend = F)+
    geom_point(aes(x =trans_adv_MLE_mean, y = reference_date), color = model_colors[2], show.legend = F)+
    geom_errorbarh(aes(xmin = trans_adv_MLE_lb, xmax = trans_adv_MLE_ub, y = reference_date), width = 0.5, color = model_colors[2], show.legend = F)+
    theme_bw()+
    scale_y_date(date_labels = "%B %d")+
    scale_x_log10()+
    scale_color_manual(values = model_colors)+
    scale_fill_manual(values = model_colors)+
    theme(axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2,size = 0.7*axis_text_size),
          axis.title.y = element_text(size = 0.7*axis_title_size),
          axis.title.x = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    labs(x = paste0('Weekly fitness advantage of BA.2.12.1 relative to BA.2'),
         y = 'Reference date',
         fill = '',
         tag = 'C')
fig4c


# Put together panels -----------------------------------------------------

design <- "
  11113
  22224
"

fig <- fig4a + fig4b + fig4c + fig4d + plot_layout(design = design)
fig
# Save fig ----------------------------------------------------------------

ggsave('../data/output/figures/figure_4_US_BA2121.png', 
       plot = fig,
       width = 15,
       height = 8)


# Switzerland BA.5-------------------------------------------------

this_country<- 'Switzerland'
country_df<-clean_global_df %>% filter(country == this_country, lineage == "BA.5")
calib_df<- country_df %>% filter(multicountry_period == 'multicountry calibration')
forecast_df <- country_df %>% filter(multicountry_period == 'multicountry forecast')
fig4a <- ggplot() +
    geom_line(data = calib_df, aes(x = ymd(collection_date), y = p_hat_mean), color = model_colors[1], size =1, show.legend = F) +
    geom_ribbon(data = calib_df, aes(x = ymd(collection_date), ymin = p_hat_lb, ymax = p_hat_ub), fill = model_colors[1], alpha = 0.6, show.legend = F)+
    geom_line(data = forecast_df, aes(x = ymd(collection_date), y = p_hat_mean), color = model_colors[1], size =0.2, show.legend = F) +
    geom_ribbon(data = forecast_df, aes(x = ymd(collection_date), ymin = p_hat_lb, ymax = p_hat_ub), fill = model_colors[1], alpha = 0.2, show.legend = F)+
    geom_point(data = country_df, aes(x = ymd(mid_week_date_recent), y = mid_week_prev_recent), size = 0.5) +
    geom_linerange(data = country_df, aes(x = ymd(mid_week_date_recent), ymin = mid_week_prev_recent -mid_week_prev_se_recent, 
                                          ymax = mid_week_prev_recent + mid_week_prev_se_recent)) +
    facet_wrap(~reference_date, ncol = 5)+
    theme_bw()+ 
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle =- 45, hjust=.6, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    scale_color_manual(values = lineage_colors)+
    scale_fill_manual(values = lineage_colors)+
    scale_x_date(date_labels = "%B %d")+
    coord_cartesian(ylim = c(0,1))+
    labs(x = '',
         y = 'BA.5 prevalence',
         fill = '',
         tag = 'A')+
    ggtitle(paste0(this_country, ' BA.5 multicountry model'))
fig4a

fig4b <- ggplot() + 
    geom_line(data = calib_df, aes(x = ymd(collection_date), y = p_hat_MLE_mean),color = model_colors[2], size =1, show.legend = F) +
    geom_ribbon(data = calib_df, aes(x = ymd(collection_date), ymin = p_hat_MLE_lb, ymax = p_hat_MLE_ub), fill = model_colors[2], alpha = 0.6, show.legend = F)+
    geom_line(data = forecast_df, aes(x = ymd(collection_date), y = p_hat_MLE_mean), color = model_colors[2], size =0.2, show.legend = F) +
    geom_ribbon(data = forecast_df, aes(x = ymd(collection_date), ymin = p_hat_MLE_lb, ymax = p_hat_MLE_ub), fill = model_colors[2], alpha = 0.2, show.legend = F)+
    geom_point(data = country_df, aes(x = ymd(mid_week_date_recent), y = mid_week_prev_recent), size = 0.5) +
    geom_linerange(data = country_df, aes(x = ymd(mid_week_date_recent), ymin = mid_week_prev_recent -mid_week_prev_se_recent, 
                                          ymax = mid_week_prev_recent + mid_week_prev_se_recent)) +
    facet_wrap(~reference_date, ncol = 5)+
    theme_bw()+ 
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=.6, vjust = 0.5, 
                                     size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    scale_color_manual(values = lineage_colors)+
    scale_fill_manual(values = lineage_colors)+
    scale_x_date(date_labels = "%B %d")+
    coord_cartesian(ylim = c(0,1))+
    labs(x = '',
         y = 'BA.5 prevalence',
         fill = '',
         tag = 'B')+
    ggtitle(paste0(this_country, ' BA.5 MLE model'))
fig4b

cleaned_metrics_MC <- country_metrics %>% 
    select(country, reference_date, multicountry_period,
           BS_med_MC, BS_mean_MC, BS_lb_MC, BS_ub_MC) %>% 
    rename(BS_med = BS_med_MC,
           BS_mean = BS_mean_MC,
           BS_lb = BS_lb_MC,
           BS_ub = BS_ub_MC) %>% 
    mutate( model = 'multicountry')
cleaned_metrics_MLE <- country_metrics %>% 
    select(country, reference_date, multicountry_period,
           BS_med_MLE, BS_mean_MLE, BS_lb_MLE, BS_ub_MLE) %>% 
    rename(BS_med = BS_med_MLE,
           BS_mean = BS_mean_MLE,
           BS_lb = BS_lb_MLE,
           BS_ub = BS_ub_MLE) %>% 
    mutate( model = 'MLE')
cleaned_metrics <- bind_rows(cleaned_metrics_MC, cleaned_metrics_MLE)

dodge<- position_dodge(width = 2)
fig4d<- cleaned_metrics %>% filter(country == this_country) %>% 
    ggplot() + geom_point(aes(x = reference_date, y = BS_mean, color = model), position = dodge) +
    geom_errorbar(aes(x = reference_date, ymin = BS_lb, ymax = BS_ub, color = model), position = dodge, width = 0.8) +
    theme_bw()+ 
    scale_y_sqrt()+
    theme(axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2, size = 0.7*axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    scale_color_manual(values = model_colors)+
    scale_x_date(date_labels = "%B %d")+
    labs(x = 'Reference date',
         y = 'Brier score',
         fill = '',
         color = 'Model',
         tag = 'D')
fig4d


fig4c <- r_comb %>% 
    filter(country == this_country, lineage == 'BA.5') %>% 
    ggplot(aes(x = transmission_advantage, y = reference_date, 
               group = ref_date_num)) + 
    coord_flip()+
    geom_density_ridges(scale = 1, fill = model_colors[1], alpha = 0.5, show.legend = F)+
    geom_point(aes(x =trans_adv_median, y = reference_date), color = model_colors[1], show.legend = F)+
    #geom_errorbarh(aes(xmin = trans_adv_lb, xmax = trans_adv_ub, y = reference_date), width = 0.3, color = model_colors[1], show.legend = F)+
    geom_point(aes(x =trans_adv_MLE_mean, y = reference_date), color = model_colors[2], show.legend = F)+
    geom_errorbarh(aes(xmin = trans_adv_MLE_lb, xmax = trans_adv_MLE_ub, y = reference_date), width = 0.5, color = model_colors[2], show.legend = F)+
    theme_bw()+
    scale_y_date(date_labels = "%B %d")+
    scale_x_log10()+
    scale_color_manual(values = model_colors)+
    scale_fill_manual(values = model_colors)+
    theme(axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2,size = 0.7*axis_text_size),
          axis.title.y = element_text(size = 0.7*axis_title_size),
          axis.title.x = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    labs(x = paste0('Weekly fitness advantage of BA.5 relative to BA.2'),
         y = 'Reference date',
         fill = '',
         tag = 'C')
fig4c


# Put together panels -----------------------------------------------------

design <- "
  11113
  22224
"

fig <- fig4a + fig4b + fig4c + fig4d + plot_layout(design = design)
fig
# Save fig ----------------------------------------------------------------

ggsave('../data/output/figures/figure_4_Switzerland_BA5.png', 
       plot = fig,
       width = 15,
       height = 8)

fig4x <- r_comb %>% 
    filter(country == this_country, lineage == 'BA.5') %>% 
    ggplot(aes(x = transmission_advantage, y = reference_date, 
               group = ref_date_num)) + 
    coord_flip()+
    geom_density_ridges(scale = 1, fill = model_colors[1], alpha = 0.5, show.legend = F)+
    geom_point(aes(x =trans_adv_median, y = reference_date), color = model_colors[1], show.legend = F)+
    #geom_errorbarh(aes(xmin = trans_adv_lb, xmax = trans_adv_ub, y = reference_date), width = 0.3, color = model_colors[1], show.legend = F)+
    geom_point(aes(x =trans_adv_MLE_mean, y = reference_date), color = model_colors[2], show.legend = F)+
    geom_errorbarh(aes(xmin = trans_adv_MLE_lb, xmax = trans_adv_MLE_ub, y = reference_date), width = 0.5, color = model_colors[2], show.legend = F)+
    theme_bw()+
    scale_y_date(date_labels = "%B %d")+
    scale_x_log10()+
    scale_color_manual(values = model_colors)+
    scale_fill_manual(values = model_colors)+
    theme(axis.text.x = element_text(angle = -45, hjust=-.5, vjust = 2,size = 0.7*axis_text_size),
          axis.title.y = element_text(size = 0.7*axis_title_size),
          axis.title.x = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    labs(x = paste0('Weekly fitness advantage of BA.5 relative to BA.2'),
         y = 'Reference date',
         fill = '',
         tag = 'C')
fig4c



second_row <- ggplot()+
    geom_point(data = SA_calibration, aes(x = p_lineage_recent, y = p_hat_mean), size = 1, color = 'cornflowerblue')+
    geom_line(data = SA_all, aes(x =p_lineage_recent , y = p_lineage_recent), color = 'black')+
    geom_point(data = SA_forecast, aes(x = p_lineage_recent, y = p_hat_mean), size = 0.5, color = 'cornflowerblue')+
    facet_wrap(~reference_date, ncol = nrow(reference_data_df), scales = 'free') +
    theme_bw() + xlab('Observed') + ylab('Multicountry predicted')
second_row
third_row <- ggplot() + geom_point(data = SA_all, aes(x = ymd(mid_week_date_recent), y = mid_week_prev_recent), size = 0.5) +
    geom_linerange(data = SA_all, aes(x = ymd(mid_week_date_recent), ymin = mid_week_prev_recent -mid_week_prev_se_recent, 
                                      ymax = mid_week_prev_recent + mid_week_prev_se_recent)) +
    geom_line(data = SA_calibration, aes(x = ymd(collection_date), y = p_hat_MLE_mean), color = 'purple', size =1) +
    geom_ribbon(data = SA_calibration, aes(x = ymd(collection_date), ymin = p_hat_MLE_lb, ymax = p_hat_MLE_ub),fill = 'purple', alpha = 0.6)+
    geom_line(data = SA_forecast, aes(x = ymd(collection_date), y = p_hat_MLE_mean), color = 'purple', size =0.5) +
    geom_ribbon(data = SA_forecast, aes(x = ymd(collection_date), ymin = p_hat_MLE_lb, ymax = p_hat_MLE_ub),fill = 'purple', alpha = 0.2)+
    facet_wrap(~reference_date, ncol = nrow(reference_data_df))+
    theme_bw()+ labs(x = '', y = '') + coord_cartesian(ylim = c(0,1)) +
    ggtitle('Spain BA.5 MLE')
third_row
fourth_row <- ggplot()+
    geom_point(data = SA_calibration, aes(x = p_lineage_recent, y = p_hat_MLE_mean), size = 1, color = 'purple')+
    geom_line(data = SA_all, aes(x =p_lineage_recent , y = p_lineage_recent), color = 'black')+
    geom_point(data = SA_forecast, aes(x = p_lineage_recent, y = p_hat_MLE_mean), size = 0.5, color = 'purple')+
    facet_wrap(~reference_date, ncol = nrow(reference_data_df), scales = 'free') +
    theme_bw() + xlab('Observed') + ylab('MLE predicted')
fourth_row
fig4_draft<-grid.arrange(first_row, second_row, third_row, fourth_row, nrow = 4)




p1 <- mu_distrib %>% 
    filter(lineage == 'BA.5') %>% 
    ggplot(aes(x = global_transmission_advantage, y = reference_date, 
               group = ref_date_num, fill = lineage), show.legend = F) + 
    scale_y_continuous(label=function(x) strftime(chron(x), "%B %d"))+
    facet_wrap(~lineage)+
    scale_color_manual(values = lineage_colors)+
    scale_fill_manual(values = lineage_colors)+
    geom_density_ridges()+
    theme_bw()+
    theme(axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    labs(x = 'Global weekly fitness advantage relative to BA.2',
         y = element_blank(),
         fill = '',
         tag = 'B')

p1

p2 <- mu_distrib %>% 
    filter(lineage == 'BA.5') %>% 
    ggplot(aes(x = reference_date, y = global_transmission_advantage, 
               group = ref_date_num, 
               fill = lineage), show.legend = F) + 
    facet_wrap(~lineage)+
    scale_color_manual(values = lineage_colors)+
    scale_fill_manual(values = lineage_colors)+
    geom_vridgeline(state = "ydensity", trim = FALSE)+
    theme_bw()+
    theme(axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size),
          plot.tag = element_text(size = tag_size))+
    labs(y = 'Global weekly fitness advantage relative to BA.2',
         x = element_blank(),
         fill = '',
         tag = 'B')

p2


# Figure 5 ----------------------------------------------------------------

pal_hemi <- dutchmasters_pal(palette = 'staalmeesters')(7)

hemi_colors = c('Northern Hemisphere' = pal_capacity[5],
                'Southern Hemisphere' = pal_capacity[7])

pal_clade <- dutchmasters_pal(palette = "view_of_Delft")(12)

clade_colors <- c("A / H1N1" = pal_clade[1],
                  "Other" = pal_clade[12],
                  "B" = pal_clade[6],
                  "A / H3N2" = pal_clade[8],
                  "A / H3" = pal_clade[11])

southern_draws <- read_csv('data/output/flu/southern/draws_for_fig.csv')
northern_draws <- read_csv('data/output/flu/northern/draws_for_fig.csv')
r_hat <- read_csv('data/output/flu/r_hats.csv')
mu_hat <- read_csv('data/output/flu/mu_hats.csv')
combined <- read_csv('data/flu/combined.csv')

most_obs_subtype <- combined %>% 
  filter(country %in% c('Ghana', 'India', 
                        'United States', 'Panama',
                        'Australia', 'Brazil',
                        'Congo, the Democatic Republic of', 'South Africa')) %>% 
  ungroup() %>% 
  filter(n > 0, 
         !is.na(n)) %>% 
  mutate(collection_date = round_date(collection_date, unit = 'week')) %>% 
  group_by(country, collection_date, subtype) %>% 
  summarize(n = sum(n)) %>% 
  group_by(country, collection_date) %>% 
  filter(n == max(n)) %>% 
  mutate(n_remaining = n()) %>% 
  mutate(subtype = if_else(n_remaining == 1,
                           subtype,
                           NA_character_)) %>% 
  group_by(country) %>% 
  arrange(collection_date) %>% 
  tidyr::fill(subtype, .direction = 'updown') %>% 
  ungroup()  %>% 
  select(collection_date, country, subtype) %>% 
  rename(most_obs_subtype = subtype,
         date = collection_date) %>% 
  mutate(most_obs_subtype= if_else(most_obs_subtype %in% c('A / H3',
                                                           'A / H1N1',
                                                           'B',
                                                           'A / H3N2'),
                                    most_obs_subtype,
                                    'Other'))

p1 <- southern_draws %>% 
  ggplot()+
  geom_line(aes(date, value, group = interaction(draw, subtype), color = subtype), show.legend = F, alpha = 0.1)+
  #geom_point(aes(x = date, y = p, color = subtype))+
  facet_wrap(~country)+
  scale_color_manual(values = clade_colors)+
  scale_fill_manual(values = clade_colors)+
  theme_bw()+
  labs(x = element_blank(),
       y = 'Estimated subtype proportion',
       tag = 'A')+
  geom_xsidecol(aes(x = date, y = log10(n), fill = most_obs_subtype),
                size = 2,
                data = combined %>% 
                  filter(country %in% c( 'Australia', 'Brazil',
                                         'Congo, the Democatic Republic of', 'South Africa')) %>% 
                  mutate(collection_date = round_date(collection_date, unit = 'week')) %>% 
                  group_by(country, collection_date) %>% 
                  summarize(n = sum(n)) %>% 
                  ungroup() %>% 
                  left_join(most_obs_subtype %>% rename(collection_date = date)) %>% 
                  ungroup() %>% 
                  rename(date = collection_date),
                show.legend = F,
                alpha = 1)+
  scale_y_continuous(n.breaks = 3)+
  scale_x_date(date_labels = "%b\n%Y")+
  scale_xsidey_continuous(position = 'right',
                          minor_breaks = NULL,
                          n.breaks = 2,
                          guide = guide_axis(TeX(r"(Collected sequences (log$_{10}$))")),
                          limits = c(0, 3))+
  theme(legend.position = 'bottom',
        strip.background = element_blank(), # element_rect(colour = "black", fill = NA)
        strip.text = element_text(size = strip_text_size))+
  theme(
    axis.text = element_text(size = axis_text_size),
    ggside.panel.scale = .35,
    axis.title = element_text(size = 13),
    plot.tag = element_text(size = tag_size),
    strip.text = element_text(size = strip_text_size),
    legend.position= "bottom",
    strip.background = element_blank(), #  element_rect(colour = "black", fill = NA)
    legend.text = element_text(size=15),
    panel.grid.minor = element_blank()
  )

p2 <- northern_draws %>% 
  ggplot()+
  geom_line(aes(date, value, group = interaction(draw, subtype), color = subtype), alpha = 0.1)+
  geom_xsidecol(aes(x = date, y = log10(n), fill = most_obs_subtype), size = 1.5,
                data = combined %>% 
                  filter(country %in% c( 'United States', 'India',
                                         'Ghana', 'Panama')) %>% 
                  mutate(collection_date = round_date(collection_date, unit = 'week')) %>% 
                  group_by(country, collection_date) %>% 
                  summarize(n = sum(n)) %>% 
                  ungroup() %>% 
                  left_join(most_obs_subtype %>% rename(collection_date = date)) %>% 
                  ungroup() %>% 
                  mutate(most_obs_subtype = factor(most_obs_subtype)) %>% 
                  rename(date = collection_date),
                show.legend = F,
                alpha = 1)+
  scale_y_continuous(n.breaks = 3)+
  scale_xsidey_continuous(position = 'right',
                          minor_breaks = NULL,
                          breaks = c(0, 3),
                          guide = guide_axis(TeX(r"(Collected sequences (log$_{10}$))")),
                          limits = c(0, 3))+
  scale_color_manual(values = clade_colors)+
  scale_fill_manual(values = clade_colors)+
  facet_wrap(~country)+
  guides(colour = guide_legend(override.aes = list(alpha = 1,
                                                   size = 4),
                               nrow = 1))+
  theme_bw()+
  labs(x = element_blank(),
       y = 'Estimated subtype proportion',
       tag = 'B',
       color = element_blank())+
  theme(
    axis.text = element_text(size = axis_text_size),
    axis.title = element_text(size = 13),
    ggside.panel.scale = .35,
    plot.tag = element_text(size = tag_size),
    strip.text = element_text(size = strip_text_size),
    legend.position= "bottom",
    strip.background = element_blank(), #  element_rect(colour = "black", fill = NA)
    legend.text = element_text(size=15),
    panel.grid.minor = element_blank())+
  scale_x_date(date_labels = "%b\n%Y")
  

p3 <- r_hat %>% 
  filter(subtype %in% c('A / H1N1',
                        'B',
                        'A / H3')) %>% 
  ggplot()+
  geom_boxplot(aes(subtype,
                   mean,
                   color = type), 
               outlier.shape = NA,
               position = position_dodge(width = 1), 
               show.legend = F)+
  geom_point(aes(subtype, mean, color = type), 
             position = position_jitterdodge(jitter.height = 0, 
                                             jitter.width = 0.3), 
             alpha = 0.2,
             show.legend = F)+
  geom_hline(aes(yintercept = 0),
             alpha = y_intercept_alpha)+
  theme_bw()+
  scale_color_manual(values = hemi_colors)+
  labs(x = element_blank(),
       y = 'Country-specific subtype fitness advantage',
       color = element_blank(),
       tag = 'C')+
  theme(legend.position = 'bottom')+
  coord_flip()+
  theme(
    axis.text = element_text(size = axis_text_size),
    axis.title = element_text(size = axis_title_size),
    plot.tag = element_text(size = tag_size),
    strip.text = element_text(size = strip_text_size),
    legend.position= "bottom",
    strip.background = element_blank(), #  element_rect(colour = "black", fill = NA)
    legend.text = element_text(size=14))

p4 <- mu_hat %>% 
  filter(subtype %in% c('A / H1N1',
                        'B',
                        'A / H3')) %>% 
  ggplot()+
  geom_point(aes(subtype,  
                 mean, 
                 color = type, 
                 group = type),
             position = position_dodge(width = 0.5), 
             size = 4)+
  geom_linerange(aes(subtype, 
                     ymin = q5,
                     ymax = q95, 
                     color = type, 
                     group = type),
                 size = 1.5,
                 position = position_dodge(width = 0.5))+
  geom_hline(aes(yintercept = 0), 
             alpha = y_intercept_alpha)+
  theme_bw()+
  scale_color_manual(values = hemi_colors)+
  labs(x = element_blank(),
       y = 'Expected subtype fitness advantage',
       color = element_blank(),
       tag = 'D')+
  theme(legend.position = 'bottom')+
  coord_flip()+
  ylim(-.4, .4)+
  theme(
    axis.text = element_text(size = axis_text_size),
    axis.title = element_text(size = 14),
    plot.tag = element_text(size = tag_size),
    strip.text = element_text(size = strip_text_size),
    legend.position= "bottom",
    strip.background = element_blank(), #  element_rect(colour = "black", fill = NA)
    legend.text = element_text(size=13))

design <- "
  113
  224
"

fig <- p1 + p2 + p3 + p4 + plot_layout(design = design)

ggsave('data/output/figures/figure_5_new.png', 
       plot = fig,
       width = 15,
       height = 8)

ggsave('data/output/figures/figure_5A.png', 
       plot = p1,
       width = 15,
       height = 8)
ggsave('data/output/figures/figure_5B.png', 
       plot = p2,
       width = 15,
       height = 8)
ggsave('data/output/figures/figure_5C.png', 
       plot = p3,
       width = 15,
       height = 8)
ggsave('data/output/figures/figure_5D.png', 
       plot = p4,
       width = 15,
       height = 8)



