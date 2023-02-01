# Databricks notebook source
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

# MAGIC %py
# MAGIC 
# MAGIC from  pyspark.sql.functions import input_file_name
# MAGIC 
# MAGIC sdf = spark.read.format('csv').options(header='True', inferSchema='True').load("dbfs:/mnt/ppi-test/validation/raw_references_dates/").withColumn("filename", input_file_name())
# MAGIC sdf.createOrReplaceTempView('raw_metadata')

# COMMAND ----------

# MAGIC %sql
# MAGIC -- Pull in all the all the distinct region-country combinations from the metadata
# MAGIC CREATE OR REPLACE TEMP VIEW region_map AS SELECT DISTINCT TRIM(UPPER(regexp_extract(location, '([A-Za-z]* ([A-Za-z]*( )?)?)', 0))) AS region, country FROM raw_metadata

# COMMAND ----------

install.packages(c('tidyverse', 'countrycode', 'emmeans'))

# COMMAND ----------

library(tidyverse)
library(lubridate)
library(countrycode)
library(emmeans)
library(arrow)
library(splines)

# COMMAND ----------

# Pull the region-country combos into memory and map to iso3 code, dropping NAs
region_map <- as_tibble(SparkR::as.data.frame(SparkR::sql('SELECT * FROM region_map;')))
region_map$code <- countrycode(region_map$country, origin = 'country.name', destination = 'iso3c')
region_map <- as_tibble(region_map) %>% filter(!is.na(code)) %>% select(-country)

# COMMAND ----------

reference_dates <- c('2022-04-30', '2022-05-16', '2022-05-27', '2022-06-04', '2022-06-27', '2022-07-01')

# COMMAND ----------

for (reference_date in reference_dates){
  
  print(paste('Fitting model for:', reference_date))
  
  # If rerunning for some reason, skip already computed results
  if (!file.exists(paste0("/dbfs/mnt/ppi-test/validation/spline_fits/fitted_model/", reference_date, ".rds"))){
    
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
    
  # pass into emmeans to get marginal slopes at t = t_max
  # This is directly the logic on line 497 of https://github.com/tomwenseleers/LineageExplorer/blob/main/global%20analysis.R
  rg_nnet <- ref_grid(fit, 
                      at = list(t = max_t_calibration, 
                                # We only want to calculate main lineages of interest here
                                lineage = c('BA.2', 'BA.1', 'BA.2.12.1', 'BA.4', 'BA.5')))
  em_nnet <- emmeans(rg_nnet,
                   specs = ~ t|lineage,
                    mode = 'latent')
    
  # Get contrasts
  contrasts <- as_tibble(pairs(em_nnet, simple = 'each', combine = TRUE)) %>% 
                filter(grepl("BA.2 - ", contrast, fixed = TRUE)) %>%
                select(-lineage)

  arrow::write_parquet(contrasts, paste0('/dbfs/mnt/ppi-test/validation/spline_fits/contrast/', reference_date, '.parquet'))
    
  # Use full to include fitting period + forecast period
  newd <- df_full %>% 
    select(t, country, region) %>% 
    distinct()
  
  # Get brier score and predicted values
  brier_score_df <- bind_cols(predict(fit, newd, type = 'p'), newd) %>% 
    pivot_longer(cols = !c(region, country, t), names_to = 'lineage', values_to = 'p') %>% 
    left_join(df_full) %>%
    mutate(fitted_count = tot_seq * p,
         brier_score = (fitted_count - n_seq)**2 / tot_seq) %>%
    select(region, country, collection_date, t, period, p, fitted_count, brier_score, tot_seq, n_seq, p_lineage_week)
  
  arrow::write_parquet(brier_score_df, paste0('/dbfs/mnt/ppi-test/validation/spline_fits/brier_score/', reference_date, '.parquet'))
  

  }
}

# COMMAND ----------


