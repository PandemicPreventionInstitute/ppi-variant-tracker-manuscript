library(readr)
library(dplyr)
library(countrycode)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(gtsummary)

lineages <- c(
  "BA.1",
  "BA.2.12.1",
  "BA.4",
  "BA.5"
)

df <- read_csv("data/processed/r_summary_2022-07-01.csv") |>
  filter(lineage %in% lineages)

df$region <- countrycode(
  sourcevar = df[[8]],
  origin = "country.name",
  destination = "un.region.name"
)

df$iso3c <- countrycode(
  sourcevar = df[[8]],
  origin = "country.name",
  destination = "iso3c"
)

####################
# Boxplot
####################

p <- df |>
  filter(!is.na(region)) |>
  ggplot() +
  geom_boxplot(
    aes(
      x = region,
      y = r_median
    ),
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(
      x = region,
      y = r_median,
      color = r_median
    ),
    width = 0.25,
    height = 0,
    alpha = 0.3,
    show.legend = FALSE
  ) +
  theme_bw() +
  facet_wrap(~lineage) +
  labs(
    x = element_blank(),
    y = "Median estimated country-specific fitness advantage"
  ) +
  scale_color_viridis_c()

ggsave("figs/region_trend.png", p)

####################
# World map
####################

for (l in lineages) {
  world <- ne_countries(scale = "medium", returnclass = "sf")

  p <- df |>
    filter(lineage == l) |>
    full_join(world %>%
      select(adm0_a3, name, geometry) %>%
      rename(iso3c = adm0_a3)) %>%
    filter(name != "Antarctica") %>%
    ggplot(aes(
      geometry = geometry,
      fill = r_median
    )) +
    geom_sf() +
    theme_classic() +
    scale_fill_viridis_c()

  ggsave(paste0("figs/", l, "_region_map.png"), p)
}

####################
# Regression
####################

# If it exists, remove the outfile. Needs to be done manually b/c loop appends
try(

  system("rm ./regression_tables.txt"),
  outFile = stdout()
)

lineage_outcome <- tibble(
  lineage = character(),
  t_stat = double()
)

for (l in lineages) {
  df_subset <- df |>
    filter(lineage == l) |>
    mutate(region = as.factor(region)) |>
    mutate(region = relevel(region, ref = "Europe"))

  fit <- lm(r_median ~ 1 + region,
    data = df_subset
  )

  tbl <- gtsummary::tbl_regression(fit,
    intercept = TRUE
  ) |>
    add_glance_source_note(include = c("r.squared", "p.value", "logLik")) |>
    as_gt() |>
    gt::tab_header(
      title = gt::md(paste("Model for variant:", l)), # nolint
      subtitle = gt::md("Linear model of
                                 country-specific fitness by region")
    ) |>
    gt::as_latex()

  cat(paste0(tbl, "\n"),
    file = "regression_tables.txt",
    append = TRUE
  )
}
