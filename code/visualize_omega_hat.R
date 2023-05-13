library(readr)
library(ggplot2)
library(dplyr)
library(patchwork)

df <- read_csv('data/processed/Omega_hat_2022-07-01.csv')

p <- df |>
  filter(comparison_id_1 > comparison_id_2) |>
  ggplot(aes(lineage_1,
             lineage_2,
             fill = Omega_hat_mean
             )
  ) +
    geom_raster() +
    theme_bw() +
    labs(x = '', y = '', fill = 'Expected  correlation') +
    scale_fill_viridis_c()

p1 <- df |>
    filter(comparison_id_1 > comparison_id_2) |>
    mutate(lineages_compared = as.factor(paste0(lineage_1,', ', lineage_2))) |>
    mutate(lineages_compared = forcats::fct_reorder(lineages_compared, Omega_hat_mean, max)) |>
    ggplot() +
      geom_hline(aes(yintercept = 0), alpha = 0.1) +
      geom_point(aes(x = lineages_compared,
                     y = Omega_hat_mean,
                     color = Omega_hat_mean
                     ),
                     show.legend = FALSE
                ) +
      geom_linerange(aes(x = lineages_compared,
                            ymin = Omega_hat_lb,
                            ymax = Omega_hat_ub,
                            color = Omega_hat_mean
                            ),
                     show.legend = FALSE
                       ) +
    scale_color_viridis_c() +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,
                                     vjust = 0.5,
                                     hjust = 1)
    ) +
    labs(x = element_blank(),
         y = "Expected pairwise correlation between lineages (95% CI)")

fig <- p1 + p
ggsave("figs/omega_hat.png",
       plot = fig,
       width = 15,
       height = 8)

ggsave("figs/omega_hat_heatmap.png",
       plot = p,
       )

ggsave("figs/omega_hat_point_interval.png",
       plot = p1,
       )






