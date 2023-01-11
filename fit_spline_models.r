# Databricks notebook source
library(tidyverse)
library(lubridate)
library(mgcv)

# COMMAND ----------

# MAGIC %fs cp dbfs:/mnt/ppi-test/GISAID_metadata/2022-07-01_metadata.csv dbfs:/mnt/ppi-test/validation/raw_reference_dates/

# COMMAND ----------

# MAGIC %fs ls dbfs:/mnt/ppi-test/validation/validation

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
# MAGIC SELECT *, regexp_extract(filename, '[0-9]*-[0-9]*-[0-9]*', 0) :: DATE AS reference_date FROM raw_metadata

# COMMAND ----------

# MAGIC %sql
# MAGIC 
# MAGIC SELECT DISTINCT reference_date FROM clean_global LIMIT 10;

# COMMAND ----------

# MAGIC %fs ls dbfs:/mnt/ppi-test/GISAID_metadata

# COMMAND ----------


