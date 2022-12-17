
library(tidyverse)

df <- read_csv('../data/output/validation/r_summary.csv')
clean_global_df <- read_csv('../data/output/validation/clean_global_df.csv')
n_seq_df <- clean_global_df %>% distinct() %>% group_by(reference_date, country, lineage) %>% summarise( n_seq = sum(n, na.rm = T)) %>% 
    filter(lineage %in% c('BA.1', 'BA.4', 'BA.5', 'BA.2.12.1', 'BA.2'))
df <- df %>% left_join(n_seq_df, by = c('reference_date', 'country', 'lineage'))

lineage_colors <- c('BA.1' = pal[10],
                    'BA.2.12.1' = pal[4],
                    'BA.4' = pal[13],
                    'BA.5' = pal[12])

df %>% 
  filter(lineage %in% c('BA.1', 'BA.4', 'BA.5', 'BA.2.12.1', 'BA.2'),
         model_convergence == T) %>% 
  mutate(dif_trans_adv_MLE = trans_adv_MLE_mean - trans_adv_median) %>% 
  filter(n_seq > 10) %>%
  mutate(z = (dif_trans_adv_MLE - mean(dif_trans_adv_MLE)) / sd(dif_trans_adv_MLE)) %>% 
  ggplot()+
  geom_point(aes(n_seq, dif_trans_adv_MLE, color = lineage), alpha = .5)+
  facet_wrap(~lineage)+
  geom_smooth(aes(n_seq, dif_trans_adv_MLE, color = lineage),
              method = 'lm', se = F)+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw()+
  scale_color_manual(values = lineage_colors)+
  theme(axis.text = element_text(size = axis_text_size),
        axis.title = element_text(size = axis_title_size),
        axis.title.x = element_text(vjust=-3),
        plot.tag = element_text(size = tag_size),
        strip.text = element_text(size = strip_text_size),
        strip.background = element_blank(),
        legend.position= "bottom",
        legend.text = element_text(size=15))+
  labs(y = TeX(r"($\hat{\beta}^{Single country}_{1ijt} - \hat{\beta}^{Multicountry}_{1ijt}$  (log$_{10}$))"),
       color = element_blank(),
       x = TeX(r"(Cumulative country-variant-day sequences ($log_{10}$))"),
       fill = element_blank())

ggsave('../data/output/figures/SFIG4_new.png',
       width = 10,
       height = 6)



df %>% 
  filter(lineage %in% c('BA.1', 'BA.4', 'BA.5', 'BA.2.12.1', 'BA.2'),
         model_convergence == T) %>% 
  select(country, lineage, n_seq, trans_adv_MLE_mean, trans_adv_median) %>% 
  rename(`Single country` = trans_adv_MLE_mean,
         Multicountry = trans_adv_median) %>% 
  pivot_longer(cols = c(`Single country`, Multicountry)) %>% 
  filter(value < 10) %>% 
  group_by(lineage, name) %>% 
  mutate(mean = mean(if_else(n_seq > quantile(n_seq, .9), value, NA_real_), na.rm = T)) %>% 
  ggplot(aes(n_seq, value, color = lineage))+
  geom_hline(aes(yintercept = mean))+
  geom_point(alpha = .1)+
  facet_grid(name ~ lineage)+
  scale_x_log10()+
  geom_smooth(method = 'lm', se = F)+
  theme_bw()+
  scale_color_manual(values = lineage_colors)+
  theme(axis.text = element_text(size = axis_text_size),
        axis.title = element_text(size = axis_title_size),
        axis.title.x = element_text(vjust=-3),
        plot.tag = element_text(size = tag_size),
        strip.text = element_text(size = strip_text_size),
        strip.background = element_blank(),
        legend.position= "bottom",
        legend.text = element_text(size=15))+
  labs(y = TeX(r"($\hat{\beta}_{1ijt}$$)"),
       color = element_blank(),
       x = TeX(r"(Cumulative country-variant-day sequences ($log_{10}$))"),
       fill = element_blank())
  
ggsave('../data/output/figures/SFIG5.png',
       width = 10,
       height = 6)

mu <- read_csv('../data/output/validation/mu_hat.csv')

df %>% 
  select(country, n_seq, reference_date, lineage) %>% 
  group_by(reference_date, lineage) %>% 
  summarize(n_seq = sum(n_seq)) %>% 
  left_join(mu)
  

mu %>% 
  ggplot(aes(last_collection_date, mu_hat_median))+
  geom_point()+
  geom_linerange(aes(x = last_collection_date, ymin = mu_hat_lb, ymax = mu_hat_ub))+
  facet_wrap(~lineage)+
  theme_bw()+
  theme(axis.text = element_text(size = axis_text_size),
        axis.text.x= element_text(angle = -90),
        axis.title = element_text(size = axis_title_size),
        axis.title.x = element_text(vjust=-3),
        plot.tag = element_text(size = tag_size),
        strip.text = element_text(size = strip_text_size),
        strip.background = element_blank(),
        legend.position= "bottom",
        legend.text = element_text(size=15))+
  labs(x = element_blank(),
       y = 'Estimated global mean fitness advantage\nfrom multicountry model')

ggsave('../data/output/figures/SFIG6.png',
       width = 10,
       height = 6)

