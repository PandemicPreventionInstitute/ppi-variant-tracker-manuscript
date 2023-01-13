# Databricks notebook source
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

install.packages(c('RcppArmadillo', 'marginaleffects'))

# COMMAND ----------

library(Rcpp)
library(RcppArmadillo)
library(tidyverse)
library(lubridate)
library(countrycode)
library(splines)
library(arrow)
library(marginaleffects)

# COMMAND ----------

sourceCpp(".//calc_infmatrix_arma.cpp")

# COMMAND ----------

# From T. Wenseleers: https://github.com/tomwenseleers/LineageExplorer
fastmultinomHess <- function(object, Z = model.matrix(object)) {
  
  probs <- object$fitted # predicted probabilities, avoid napredict from fitted.default
  
  coefs <- coef(object)
  if (is.vector(coefs)){ # ie there are only 2 response categories
    coefs <- t(as.matrix(coefs))
    probs <- cbind(1 - probs, probs)
  }
  coefdim <- dim(coefs)
  p <- coefdim[2L] # nr of parameters
  k <- coefdim[1L] # nr out outcome categories-1
  ncoefs <- k * p # nr of coefficients
  n <- dim(Z)[1L] # nr of observations
  
  #  Now compute the Hessian = the observed (= expected, in this case) 
  #  Fisher information matrix info
  
  # info <- calc_infmatrix(probs = probs[, -1, drop=F], # pure R function
  #                        Z = Z, 
  #                        row_totals = object$weights) 
  
  info <- calc_infmatrix_RcppArma(probs = probs[, -1, drop=F], # using faster RcppArmadillo function
                                  Z = Z, 
                                  row_totals = object$weights) 
  # note: this could still be parallelized either within Rcpp code with parallelReduce or
  # on R side
  
  Names <- dimnames(coefs)
  if (is.null(Names[[1L]])) Names <- Names[[2L]] else Names <- as.vector(outer(Names[[2L]], Names[[1L]],
                                                                               function(name2, name1)
                                                                                 paste(name1, name2, sep = ":")))
  dimnames(info) <- list(Names, Names)
  
  return(info)
}

# COMMAND ----------

region_map <- as_tibble(SparkR::as.data.frame(SparkR::sql('SELECT * FROM region_map;'))) %>%
  mutate(code = countrycode(country, origin = 'country.name', destination = 'iso3c')) %>% 
  filter(!is.na(code)) %>%
  select(-country)

# COMMAND ----------

reference_date = reference_dates[1]

# COMMAND ----------

df <- read_csv(paste0('/dbfs/mnt/ppi-test/validation/validation_data/lineage_t_', reference_date, '.csv')) %>%
  filter(period == 'calibration') %>%
  left_join(region_map, by = 'code') %>%
  filter(!is.na(region))

max_t_calibration <- df %>% filter(t == max(t)) %>% pull(t) %>% unique()

head(df)

# COMMAND ----------

df$lineage <- relevel(as.factor(df$lineage), ref= 'BA.2')

# COMMAND ----------

library(mgcv)

# COMMAND ----------

fit <- nnet::multinom(lineage ~ 1 + ns(t, df=2)+ns(t, df=2):region+country, 
                                       weights=n_seq, 
                                       data=df, 
                                       maxit=10000, MaxNWts=100000)

# COMMAND ----------

library(mgcv)

# COMMAND ----------

K = df %>% pull(lineage) %>% unique() %>% length()

# COMMAND ----------

K

# COMMAND ----------

fit <- mgcv::gam(list(n_seq ~ 1 + s(t, k = 3, bs = 'cs'), ~ ), family = multinom(K = 2), data = df %>% filter(n_seq > 0, ))

# COMMAND ----------

fit <- mgcv::bam()

# COMMAND ----------

# MAGIC %fs mkdirs dbfs:/mnt/ppi-test/validation/spline_fits/fitted_model

# COMMAND ----------

fit <- readRDS('/dbfs/mnt/ppi-test/validation/spline_fits/fitted_model/2022-04-30.rds')

# COMMAND ----------

fit$Hessian <- fastmultinomHess(fit, model.matrix(fit))
fit$vcov <- vcov(fit)
saveRDS(fit, file=paste0("/dbfs/mnt/ppi-test/validation/spline_fits/fitted_model/", reference_date, ".rds"))

# COMMAND ----------

# MAGIC %fs ls dbfs:/mnt/ppi-test/validation/spline_fits/fitted_model/

# COMMAND ----------

datagrid(newdata = data.frame(t = max_t_calibration))

# COMMAND ----------

df %>% filter(t == max(t)) %>% select(country, t, lineage)

# COMMAND ----------

marginaleffects::slopes(model = fit,
       newdata = df %>% filter(t == max(t)) %>% select(country, t, lineage))

# COMMAND ----------

meffects_byregion <- marginaleffects(fit, 
                                               type = "latent", # = additive log-ratio = growth rate advantage relative to BQ.1*
                                               variables = c("t"),
                                               by = c("group", "region"),
                                               vcov = fit$vcov,
                                               newdata = datagrid(t = 89,
                                                                  region = unique(df$region)))

# COMMAND ----------

# with new faster marginaleffects code
meffects <- marginaleffects(fit, 
                                               type = "latent", # = additive log-ratio = growth rate advantage relative to BQ.1*
                                               variables = c("t"),
                                               by = c("group"),
                                               vcov = fit$vcov,
                                               newdata = datagrid(t = max_t_calibration
                                               ))


# COMMAND ----------

# growth rate advantage compared to reference level BA.5
meffects_byregion <- marginaleffects(fit_best, 
                                               type = "link", # = additive log-ratio = growth rate advantage relative to BQ.1*
                                               variables = c("t"),
                                               by = c("group", "region"),
                                               vcov = fit_best$vcov,
                                               newdata = datagrid(date_num = today_num,
                                                                  region = unique(data_agbyweekcountry1$region)
                                               ))) # 30s


# COMMAND ----------

for (reference_date in reference_dates){
  
  print("Reading in data f")
  df <- read_csv(paste0('/dbfs/mnt/ppi-test/validation/validation_data/lineage_t_', reference_date, '.csv')) %>%
    filter(period == 'calibration') %>%
    left_join(region_map, by = 'code') %>%
    filter(!is.na(region))
  
  
  
  
  
}
