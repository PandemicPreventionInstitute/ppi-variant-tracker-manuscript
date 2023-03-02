library(tidyverse)
library(arrow)
library(DBI)
library(duckdb)
library(latex2exp)


source("code/plot_theme.R")

con <- dbConnect(duckdb::duckdb(),
  dbdir = ":memory:"
)

# Setup
dbExecute(con,
  statement = read_file("~/.duckdbrc")
)

path_prefix <- "s3://ppi-dev/validation/spline_fits/"
predicted <- paste0(path_prefix, "predicted/*.parquet")
contrast <- paste0(path_prefix, "contrast/*.parquet")
brier_score <- paste0(path_prefix, "brier_score/*.parquet")
baseline <- "s3://ppi-dev/validation/lineage_t_for_comp_2022-07-01.csv"

# Create a grid of all reference date-collection date-lineage combinations to join onto
grid <- dbGetQuery(
  con,
  paste0("
SELECT DISTINCT collection_date,
lineage,
regexp_extract(filename, '[0-9]*-[0-9]*-[0-9]*', 0) :: DATE AS reference_date
FROM read_parquet('", predicted, "', filename=TRUE)
WHERE country = 'Portugal'
AND lineage in ('BA.2', 'BA.4', 'BA.5', 'BA.2.12.1', 'BA.1', 'other')
")
) %>%
  expand(collection_date, lineage, reference_date)

# Make available for  joining
dbWriteTable(con, "grid", grid)

d <- dbGetQuery(
  con,
  paste0("
SELECT
pred.lineage,
pred.reference_date,
pred.country,
pred.collection_date,
pred.p_hat,
dat.p_lineage,
dat.p_lineage_se
FROM grid g
FULL JOIN
(SELECT
        lineage,
        country,
        collection_date,
        p_lineage,
        p_lineage_se
    FROM '", baseline, "'
    WHERE lineage IN ('BA.2', 'BA.4', 'BA.5', 'BA.2.12.1', 'BA.1', 'other')) dat
ON g.lineage = dat.lineage
AND g.collection_date = dat.collection_date
FULL JOIN (SELECT *,
      regexp_extract(filename, '[0-9]*-[0-9]*-[0-9]*', 0) :: DATE AS reference_date
      FROM read_parquet('", predicted, "', filename = true)) pred
ON pred.lineage = dat.lineage
AND pred.country = dat.country
AND pred.collection_date = dat.collection_date
AND pred.reference_date = g.reference_date
WHERE pred.country in ('Portugal', 'Brazil', 'Spain', 'India', 'Czechia')
AND pred.lineage IN ('BA.2', 'BA.4', 'BA.5', 'BA.2.12.1', 'BA.1', 'other')
                        ")
)

d %>%
  group_by(reference_date) %>%
  summarize(max_collection_date = max(collection_date)) %>%
  print()

for (lineage_to_show in c(
  "BA.2",
  "BA.4",
  "BA.5",
  "BA.2.12.1",
  "BA.1",
  "other"
)) {
  p1 <- d %>%
    filter(lineage == lineage_to_show) %>%
    ggplot() +
    geom_point(
      aes(collection_date,
        as.double(p_lineage),
        alpha = if_else(collection_date < reference_date, 0.25, 1)
      ),
      show.legend = FALSE
    ) +
    geom_line(aes(collection_date, p_hat)) +
    geom_vline(aes(xintercept = reference_date),
      alpha = 0.2
    ) +
    facet_grid(country ~ reference_date) +
    supp_theme(date = T) +
    labs(y = paste0('Prevalance of lineage: ', lineage_to_show))

  ggsave(paste0("figs/supp/spline_supp_", lineage_to_show, ".png"), p1)
}

p1 <- d %>%
  filter(country == "India") %>%
  ggplot() +
  geom_point(
    aes(collection_date,
      as.double(p_lineage),
      alpha = if_else(collection_date < reference_date, 0.25, 1)
    ),
    show.legend = FALSE
  ) +
  geom_line(aes(collection_date, p_hat, color = lineage)) +
  geom_vline(aes(xintercept = reference_date),
    alpha = 0.2
  ) +
  facet_grid(lineage ~ reference_date) +
  theme_bw()

ggsave(paste0("figs/supp/spline_supp_", "india", ".png"), p1)



p2 <- dbGetQuery(
  con,
  paste0("
SELECT *,
regexp_extract(filename, '[0-9]*-[0-9]*-[0-9]*', 0) :: DATE AS reference_date
FROM read_parquet('", contrast, "', filename  = TRUE)
                  ")
) %>%
  ggplot(aes(reference_date, estimate)) +
  geom_line() +
  geom_ribbon(
    aes(
      x = reference_date,
      ymin = estimate - 2 * SE,
      ymax = estimate + 2 * SE
    ),
    alpha = 0.2
  ) +
  facet_wrap(~contrast) +
  supp_theme(date = T) +
  labs(y = "Difference in estimated marginal means")

ggsave("figs/supp/contrasts.png", p2)

d <- dbGetQuery(
  con,
  paste0("
SELECT
    s.country,
    s.reference_date,
    TRY_CAST(s.spline_bs AS FLOAT) AS 'Spline-based model',
    TRY_CAST(m.BS_mean_mc AS FLOAT) AS 'Multicountry model (main text)'
FROM
(SELECT brier_score AS spline_bs,
country,
regexp_extract(filename, '[0-9]*-[0-9]*-[0-9]*', 0) :: DATE AS reference_date
FROM read_parquet('", brier_score, "', filename  = TRUE)) s
JOIN read_csv_auto('s3://ppi-dev/validation/country_metrics.csv') m
ON s.country = m.country
AND s.reference_date = m.reference_date
WHERE m.country in ('Portugal', 'Brazil', 'Bangladesh', 'India', 'Czechia', 'Russia',
                    'Morocco', 'Indonesia', 'Nigeria')
                  ")
) %>%
  pivot_longer(c("Spline-based model", "Multicountry model (main text)"))

p3 <- ggplot(d, aes(
  x = reference_date,
  y = value,
  group = name,
  color = name
)) +
  geom_point(position = position_dodge(width = 5)) +
  supp_theme(date=T,
             legend_color_size = 5) +
  facet_wrap(~country) +
  labs(x = 'Brier score (based on final reference dataset)',
       color = element_blank())

ggsave("figs/supp/brier_score.png", p3)


dbDisconnect(con,
  shutdown = TRUE
)
