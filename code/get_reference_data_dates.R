# Author: Kaitlyn Johnson

# Date initiated: 06-21-2022

# This script reads in the datasets saved in Domino and generates a list of:
# 1. The reference dates (tagged onto the end)
# 2. The last collection dates 
# 3. The last submission dates

USE_CASE = Sys.getenv("USE_CASE")
if(USE_CASE == ""){
    USE_CASE<-'local'
}

if (USE_CASE== 'domino'){
    install.packages("tidyverse", dependencies = TRUE, repos = 'http://cran.us.r-project.org')
    install.packages("janitor", dependencies = TRUE, repos = 'http://cran.us.r-project.org')
    install.packages("tibble", dependencies = TRUE, repos = 'http://cran.us.r-project.org')
    install.packages("countrycode", dependencies = TRUE, repos = 'http://cran.us.r-project.org')
    install.packages("lubridate", dependencies=TRUE, repos='http://cran.us.r-project.org')
    install.packages("R.utils", dependencies=TRUE, repos='http://cran.us.r-project.org')
    install.packages("stringr", dependencies=TRUE, repos='http://cran.us.r-project.org')
    install.packages("dplyr", dependencies=TRUE, repos='http://cran.us.r-project.org')
    install.packages("readr", dependencies=TRUE, repos='http://cran.us.r-project.org')
}

library(tidyverse) # data wrangling
library(tibble) # data wrangling
library(janitor) # column naming
library(countrycode) # country codes
library(lubridate) # date times
library(R.utils) # R utilities
library(stringr) # to parse strings in R
library(dplyr) # data wrangling
library(readr) # read_csv

# Load in reference data dataframe --------------------------------------------
if (USE_CASE == 'local'){
    reference_data_df<- read_csv("../data/processed/reference_data_df.csv")
}
if (USE_CASE == 'domino'){
    reference_data_df <- read_csv("/mnt/data/processed/reference_data_df.csv")
}

reference_date_list <- today()





for (i in 1:length(reference_date_list)){
    THIS_REF_DATE<- as.character(reference_date_list[i])
    if (USE_CASE == 'local'){
        GISAID_METADATA_PATH <- paste0('../data/raw/', THIS_REF_DATE, '_metadata.csv')
    }
    if (USE_CASE == 'domino'){
        GISAID_METADATA_PATH <- paste0('/mnt/data/raw/', THIS_REF_DATE, '_metadata.csv')
    }
    
    this_metadata <- read.csv(GISAID_METADATA_PATH)
    
    this_metadata$collection_date<-ymd(this_metadata$collection_date)
    this_metadata$submission_date<-ymd(this_metadata$submission_date)
    
    reference_date <- ymd(THIS_REF_DATE)
    last_submission_date <- ymd(max(this_metadata$submission_date, na.rm = T))
    last_collection_date <- ymd(max(this_metadata$collection_date, na.rm = T))
    n_pango_lineages <- length(unique(this_metadata$pango_lineage))
    
    this_df <- data.frame(reference_date, last_submission_date, last_collection_date, n_pango_lineages)
    reference_data_df <- bind_rows(reference_data_df, this_df)
    
}

reference_data_df <- unique(reference_data_df)
reference_data_df<- reference_data_df[order(reference_data_df$reference_date),]

# Unit test to check that GISAID feed is up to date
stopifnot('Last submission date is more than 3 days ago'= 
          as.numeric(max(reference_data_df$reference_date) - max(reference_data_df$last_submission_date))<=3)
if (USE_CASE == 'local'){
    write.csv(reference_data_df, '../data/processed/reference_data_df.csv', row.names = F)
}

if (USE_CASE == 'domino'){
    write.csv(reference_data_df, '/mnt/data/processed/reference_data_df.csv', row.names = F)
}