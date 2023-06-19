library(duckdb)
library(dplyr)
library(ggplot2)

con <- dbConnect(duckdb())

df <- dbGetQuery(con, "

SELECT
    mu_hat_ub,
    mu_hat_lb,
    mu_hat_median,
    lineage,
    fit_date,
    dominant_lineage
FROM 'data/processed/mu_hat.parquet'
WHERE lineage NOT ILIKE 'other';

")

df |>
    group_by(lineage) |>
    mutate(n = n()) |>
    ungroup() |>
    filter(n > 35) |>
    mutate(dominant_lineage = if_else(fit_date < lubridate::ymd('2023-01-05'), 'BA.5', dominant_lineage)) |>
    ggplot(aes(fit_date, mu_hat_median, color = lineage)) +
    geom_hline(aes(yintercept = 0), alpha = 0.2) +
    geom_line(show.legend=T) +
    facet_wrap(~dominant_lineage, scale = 'free_x') +
    theme_bw()


#  -----------------------------------------------------------------------------

df <- dbGetQuery(con, "

SELECT
    mu_hat_ub,
    mu_hat_lb,
    mu_hat_median,
    lineage,
    fit_date,
    dominant_lineage
FROM 'data/processed/mu_hat.parquet'
WHERE lineage IN ('CH.1.1', 'XBB.1.5', 'BA.2.75', 'BF.7.14');

")


p <- df |>
    ggplot() +
    geom_hline(aes(yintercept = 0), alpha = 0.2) +
    geom_ribbon(aes(
                x = fit_date,
                ymin = mu_hat_lb,
                ymax = mu_hat_ub,
                fill = dominant_lineage,
                group = dominant_lineage
                ),
                alpha = 0.33,
                show.legend = FALSE)+
    geom_line(aes(
                x = fit_date,
                y = mu_hat_median,
                color = dominant_lineage,
                group = dominant_lineage
                )) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
    labs(
         x = 'Reference date for model fit',
         y = 'Fitted relative growth rate (95% CI)',
         color = 'Reference lineage'
         ) +
    facet_wrap(~lineage, scale = 'free_x') +
    scale_color_viridis_d() +
    scale_fill_viridis_d()

ggsave("figs/CH.1.1.jpeg", p)

#  -----------------------------------------------------------------------------

dbDisconnect(con, shutdown = TRUE)
