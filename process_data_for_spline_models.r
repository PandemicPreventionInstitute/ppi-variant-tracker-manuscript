# Databricks notebook source
# MAGIC %fs mkdirs dbfs:/mnt/ppi-test/validation/spline_fits/brier_score/

# COMMAND ----------

# Add code to install arrow package 
# fast install from https://arrow.apache.org/docs/r/articles/install.html#method-1a---binary-r-package-containing-libarrow-binary-via-rspmconda
options(
  HTTPUserAgent =
    sprintf(
      "R/%s R (%s)",
      getRversion(),
      paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])
    )
)
 
install.packages("arrow", repos = "https://packagemanager.rstudio.com/all/__linux__/focal/latest")

# COMMAND ----------

reference_dates <- c('2022-04-30', '2022-05-16', '2022-05-27', '2022-06-04', '2022-06-27', '2022-07-01')

# COMMAND ----------

# MAGIC %py
# MAGIC 
# MAGIC from  pyspark.sql.functions import input_file_name
# MAGIC 
# MAGIC sdf = spark.read.format('csv').options(header='True', inferSchema='True').load("dbfs:/mnt/ppi-test/validation/raw_references_dates/").withColumn("filename", input_file_name())
# MAGIC sdf.createOrReplaceTempView('raw_metadata')

# COMMAND ----------

# MAGIC %sql
# MAGIC 
# MAGIC CREATE OR REPLACE TEMP VIEW region_map AS SELECT DISTINCT TRIM(UPPER(regexp_extract(location, '([A-Za-z]* ([A-Za-z]*( )?)?)', 0))) AS region, country FROM raw_metadata

# COMMAND ----------

install.packages(c('RcppArmadillo', 'tidyverse', 'countrycode', 'emmeans'))

# COMMAND ----------

library(tidyverse)
library(lubridate)
library(countrycode)
library(emmeans)
library(arrow)

# COMMAND ----------

region_map <- as_tibble(SparkR::as.data.frame(SparkR::sql('SELECT * FROM region_map;')))

# COMMAND ----------

region_map$code <- countrycode(region_map$country, origin = 'country.name', destination = 'iso3c')

# COMMAND ----------

region_map <- as_tibble(region_map) %>% filter(!is.na(code)) %>% select(-country)

# COMMAND ----------

region_map

# COMMAND ----------

for (reference_date in reference_dates){
  
  print(paste('Fitting model for:', reference_date))
  
  # Read in lineage file
  df_full <- read_csv(paste0('/dbfs/mnt/ppi-test/validation/validation_data/lineage_t_', reference_date, '.csv')) %>%
      left_join(region_map, by = 'code') %>%
      filter(!is.na(region))
  
  df <- df_full %>%
    filter(period == 'calibration')

  max_t_calibration <- df %>% 
                        filter(t == max(t)) %>% 
                        pull(t) %>% unique()
  
  df$lineage <- relevel(as.factor(df$lineage), ref= 'BA.2')

  head(df)
  
  fit <- nnet::multinom(lineage ~ 1 + ns(t, df=2)+ns(t, df=2):region+country, 
                                       weights=n_seq, 
                                       data= df, 
                                       maxit=10000, MaxNWts=100000)
  
  # Save model
  saveRDS(fit, file=paste0("/dbfs/mnt/ppi-test/validation/spline_fits/fitted_model/", reference_date, ".rds"))
  
  # pass into emmeans
  rg_nnet <- ref_grid(fit, at = list(t = max_t_calibration))
  em_nnet <- emmeans(rg_nnet,
                   specs = pairwise ~ t|lineage)
  
  # regrid to get coefficients as logit
  em_nnet_logit <- regrid(em_nnet[[1]],
                        transform = "logit")
  
  # Get contrasts
  contrasts <- as_tibble(pairs(em_nnet_logit, simple = 'each', combine = TRUE)) %>% 
                filter(grepl("BA.2 - ", contrast, fixed = TRUE))
  
  write_csv(contrasts, paste0('/dbfs/mnt/ppi-test/validation/spline_fits/contrast/', reference_date, '.csv'))
  
  # Use full to include fitting period + forecast period
  newd <- df_full %>% 
    select(t, country, region) %>% 
    distinct()
  
  # Get brier score
  brier_score_df <- bind_cols(predict(fit, newd, type = 'p'), newd) %>% 
    pivot_longer(cols = !c(region, country, t), names_to = 'lineage', values_to = 'p') %>% 
    left_join(df_full) %>%
    mutate(fitted_count = tot_seq * p,
         brier_score = (fitted_count - n_seq)**2 / tot_seq) %>%
    select(region, country, collection_date, t, p, fitted_count, brier_score, tot_seq, n_seq, p_lineage_week)
  
  write_csv(brier_score_df, paste0('/dbfs/mnt/ppi-test/validation/spline_fits/brier_score/', reference_date, '.csv'))
  
}
