# Author: Kaitlyn Johnson


# This script loads in the gisaid metadata from 2022-07-01
# 1. The number of sequences for every country-pango_lineage-timepoint 
# 2. The number of sequences for every country-lineage-timepoint, where we use 
# script to specify the list of lineage buckets to aggregate to
# 3. The number of sequences for every country-variant-were-tracking+all-others-
# timepoint 
# 4. The number of cases per country-timepoint reported from OWID


rm(list = ls())

# Set global parameters-------------------------------------------------------

USE_CASE = Sys.getenv("USE_CASE")
if(USE_CASE == ""){
  USE_CASE<-'local'
}

# Install Libraries-------------------------------------------------------------

if (USE_CASE== 'domino'){
    install.packages("tidyverse", dependencies = TRUE, repos = 'http://cran.us.r-project.org')
    install.packages("janitor", dependencies = TRUE, repos = 'http://cran.us.r-project.org')
    install.packages("tibble", dependencies = TRUE, repos = 'http://cran.us.r-project.org')
    install.packages("countrycode", dependencies = TRUE, repos = 'http://cran.us.r-project.org')
    install.packages("lubridate", dependencies=TRUE, repos='http://cran.us.r-project.org')
    install.packages("readxl", dependencies=TRUE, repos='http://cran.us.r-project.org')
    install.packages("zoo", dependencies=TRUE, repos='http://cran.us.r-project.org')
    install.packages("R.utils", dependencies=TRUE, repos='http://cran.us.r-project.org')
    install.packages("stringr", dependencies=TRUE, repos='http://cran.us.r-project.org')
    install.packages("tsoutliers", dependencies=TRUE, repos='http://cran.us.r-project.org')
    install.packages("dplyr", dependencies=TRUE, repos='http://cran.us.r-project.org')
    install.packages("scales", dependencies=TRUE, repos='http://cran.us.r-project.org')
    install.packages("readr", dependencies=TRUE, repos='http://cran.us.r-project.org')
    install.packages("EpiEstim", dependencies=TRUE, repos='http://cran.us.r-project.org')
    install.packages("lme4", dependencies = TRUE, repos = 'http://cran.us.r-project.org')
}

library(tidyverse) # data wrangling
library(tibble) # data wrangling
library(janitor) # column naming
library(countrycode) # country codes
library(lubridate) # date times
library(readxl) # excel import
library(zoo) # calculate rolling averages
library(R.utils) # R utilities
library(stringr) # to parse strings in R
library(dplyr) # data wrangling
library(readr) # read_csv
library(lme4) # logistic regression

COMBINE_BA.4_BA.5 <- FALSE
# Set the variants we're tracking
if (COMBINE_BA.4_BA.5 == TRUE){
    variants_list<- c('BA.2', 'BA.1', 'BA.2.12.1', 'BA.4/BA.5') 
}
if (COMBINE_BA.4_BA.5 == FALSE){
    variants_list<- c('BA.2', 'BA.1', 'BA.2.12.1', 'BA.4', 'BA.5') 
}


variants_tracking<-variants_list
VOCs<- c('B.1.1.7', 'P.1', 'AY', 'C.37', 'B.1.621', 'B.1.1.529', 'B.1351', 'B.1.429', 'BA.1')
# set time window
duration<-90# 120



#local
if (USE_CASE == 'local'){
  REFERENCE_DATA_PATH <- '../data/processed/reference_data_df.csv'
  GISAID_METADATA_PATH<-"../data/raw/2022-07-01_metadata.csv" # from extracted datastream
  OWID_PATH<-url('https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv')
 
}
#Domino
if (USE_CASE == 'domino'){
  GISAID_METADATA_PATH<-"/mnt/data/raw/2022-07-01_metadata.csv" # from extracted datastream
  OWID_PATH<-url('https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv')

}

#Download metadata
metadata<-read.csv(GISAID_METADATA_PATH)
reference_data_df<-read_csv(REFERENCE_DATA_PATH)


TODAY_DATE<-reference_data_df$last_submission_date[reference_data_df$reference_date== ymd("2022-07-01")]
FIRST_DATE<-as.character(ymd(TODAY_DATE)- days(duration))# earliest date we want sequences for
LAST_DATE<-TODAY_DATE # want flexibility to change this if we want to pretend we were looking from a different time period
OBSERVATION_THRESHOLD<-50 # total global number of sequences needed for it to be fit to the model (go into lineage_t)

owid_raw<-read_csv(OWID_PATH)
print('OWID and GISAID data loaded successfully')





# Functions to process data -----------------------------------

clean_metadata<- function(metadata, FIRST_DATE, TODAY_DATE, START_DATE, duration, country_list){
    # input: raw metadata (with sub pango variants we're tracking added)
    # output: cleaned metadata subsetted to the time period we are looking at, 
    # and, if we include a country_list, the set of countries we want to focus on
    
    metadata<- metadata %>%
        separate(location,
                 into = c("continent", "country", "division", "location"),
                 sep = " / | /|/ |/")%>%
        filter(collection_date >= ymd(FIRST_DATE))
    #mutate(lag = ymd(submission_date) - ymd(collection_date))
    
    # 3. Fix random things in location names
    metadata$country[metadata$country == "USA"] <- "United States"
    # replace Usa acronym with United States
    metadata$country[metadata$country == "Usa"] <- "United States"
    # replace DC with District of Columbia
    metadata$division[metadata$division == "DC"] <- "District of Columbia"
    # capitalize first letter of country
    metadata$country <- capitalize(metadata$country)
    # correct mispelling
    metadata$country[metadata$country == "Cote dIvoire"] <- "Cote d'Ivoire"
    # correct mispelling
    metadata$country[metadata$country == "Niogeria"] <- "Nigeria"
    # correct mispelling
    metadata$country[metadata$country == "Republic of the Congo"] <- "Congo"
    # correct mispelling
    metadata$country[metadata$country == "Czech republic"] <- "Czech Republic"
    # correct misentry of Lithuania
    metadata$country[metadata$country == "Jonavos apskritis"] <- "Lithuania"
    # correct misreading
    metadata$country[metadata$country == "M?xico"] <- "Mexico"
    
    # 4. Deal with date issues
    # Assign submissions with missing date to the 15th of the month
    #metadata$collection_date <- ifelse(nchar(as.character(metadata$collection_date)) == 7, paste(as.character(metadata$collection_date), "15", sep = "-"), as.character(metadata$collection_date))
    # format as dates
    metadata$collection_date<-ymd(metadata$collection_date)
    metadata$submission_date<-ymd(metadata$submission_date)
    # exclude submissions earlier than 2019-12-01
    metadata<- metadata[metadata$collection_date >= as.Date(FIRST_DATE, format = "%Y-%m-%d"),]
    # exclude submissions dated to the future
    metadata <- metadata[metadata$collection_date <= as.Date(TODAY_DATE, format = "%Y-%m-%d") - days(1),]
    
    
    # 5. Generate country codes from GISAID country names
    metadata$code <- countrycode(metadata$country, origin = 'country.name', destination = 'iso3c')
    
    
    # inserts missing country codes
    metadata$code[metadata$country == "Micronesia (country)"] <- "FSM"
    metadata$code[metadata$country == "Timor"] <- "TLS"
    metadata$code[metadata$country == "Turks and Caicos Islands"] <- "TCA"
    metadata$code[metadata$country == "Nauru"] <- "NRU"
    metadata$code[metadata$country == "Kosovo"] <- "XKX"
    metadata$code[metadata$country == "Guernsey"] <- "GGY"
    metadata$code[metadata$country == "Falkland Islands"] <- "FLK"
    
    # Remove any rows that don't contain country codes (i.e. not a valid country)
    metadata<-metadata[!is.na(metadata$code),]
    # Remove any rows that don't contain valid collection or submission dates
    metadata<-metadata[!is.na(metadata$collection_date),]
    metadata<-metadata[!is.na(metadata$submission_date),]
    
    #6. Filter dataset so collection dates are from the past 120 days
    if (length(country_list)==0){
        metadata_subset<-metadata
    }
    
    if(length(country_list)>0){
        # Add a way to handle both codes and countries
        if (mean(nchar(country_list[1])) ==3){
            metadata_subset<-metadata %>% filter(code %in% country_list)
        }
        if (mean(nchar(country_list)) >3){
            metadata_subset<-metadata %>% filter(country %in% country_list)
        } # will eventually delete but leave for now
    }

     metadata_subset<-metadata_subset%>%
        #filter(country %in% country_list)%>%
        filter(submission_date <= (ymd(TODAY_DATE) - days(1)), # Exclude sequences submitted today 
               collection_date >= (ymd(START_DATE) - days(duration)))# exclude sequences more than 120 days ago
     print(paste0(length(unique(metadata_subset$code)), ' countries in dataset'))

    
    # name all missing pango_lineages 
    metadata_subset$pango_lineage[is.na(metadata_subset$pango_lineage) | (metadata_subset$pango_lineage == "") | (metadata_subset$pango_lineage == " null")]<-"Unassigned" # set those with NA to None as well
    most_recent_submission_date<-max(metadata_subset$submission_date, na.rm = T)
    print(paste0('Most recent submission from ', most_recent_submission_date))
    stopifnot('Most recent submission is more than 5 days ago' = as.integer(ymd(TODAY_DATE) - most_recent_submission_date ) >0)
    return(metadata_subset)
}


# Data aggregation
get_lineage_country_day<-function(metadata_subset, VARIANT_TYPE, LAST_DATE, OBSERVATION_THRESHOLD ){
    # input: metadata (line list) subsetted to dates we want
    # VARIANT_TYPE: level you want to aggregate to (pango_lineage, bucketed lineages, variants we're tracking)
    # LAST_DATE: date we want to cast out to, for mapping to model prediction. Could be the last day we observe data(max(collection_date)), or could be today's date, or a week into the future, etc.
    # output: dataset with every possible combination of pango lineage, country,
    # and date for time period of interest
    
    if (VARIANT_TYPE == 'pango_lineage'){
        gisaid_t<-metadata_subset %>%
            group_by(code, collection_date, pango_lineage) %>%
            summarise(n_seq = n()) %>%
            rename(lin = pango_lineage) # lin is generic term for VARIANT_TYPE
        n_lins<-length(unique(metadata_subset$pango_lineage))
    }
    
    if (VARIANT_TYPE == 'lineage'){
     # Find lineages that are below observation threshold
        lineages_rarely_observed <- metadata_subset %>% 
            group_by(lineage) %>% 
            summarize(n = n()) %>% 
            filter(n<OBSERVATION_THRESHOLD) %>% 
            pull(lineage)
    # Combine these into others    
        gisaid_t<- metadata_subset %>%
            mutate(lineage = case_when
                   (lineage %in% lineages_rarely_observed ~'other',
                    lineage == "Unassigned" ~ 'other',
                     TRUE ~ lineage)) %>% 
            group_by(code, collection_date, lineage) %>%
            summarise(n_seq = n()) %>%
            rename(lin = lineage) # lin is generic term for VARIANT_TYPE
        n_lins<-length(unique(gisaid_t$lin))
    }
    if (VARIANT_TYPE == 'variant'){
        gisaid_t<-metadata_subset %>%
            mutate(variant = case_when(
                variant == 'Unassigned' ~'other',
                TRUE ~ variant)) %>% 
            group_by(code, collection_date, variant) %>%
            summarise(n_seq = n()) %>%
            rename(lin = variant) # lin is generic term for VARIANT_TYPE
        n_lins<-length(unique(metadata_subset$variant))
    }
            
    
    gisaid_t<- gisaid_t %>%select(code, collection_date, n_seq, lin) %>%
            ungroup()
    
    # Complete variant-country-days
    unique_codes <- unique(gisaid_t$code)
    n_codes <- length(unique_codes)
    collection_dates<-seq(from = min(gisaid_t$collection_date), to = ymd(LAST_DATE)-days(1), by = 1)
    t <- seq(from = 1, to = length(collection_dates), by = 1)
    n_dates<-length(collection_dates)
    date_t<-data.frame(t,collection_dates)

    gisaid_orig<-gisaid_t
    gisaid_t<-gisaid_t %>%complete(collection_date =seq(min(gisaid_t$collection_date), ymd(LAST_DATE)-days(1), by = "days"))
    list_of_variant_country_days<-gisaid_t%>%
        complete(lin, code, collection_date)%>%
        select(lin, code, collection_date) %>%
        filter(!is.na(lin)) %>% 
        filter(!is.na(code)) %>% 
        mutate(period = case_when(
            collection_date <= max(gisaid_orig$collection_date) ~ 'calibration',
            collection_date > max(gisaid_orig$collection_date) ~ 'nowcast'
        ))


    gisaid_t<-left_join(list_of_variant_country_days, gisaid_t, by = c("collection_date", "code", "lin"))
    gisaid_t<-left_join(gisaid_t, date_t, by = c("collection_date" = "collection_dates"))
    # check that the length of the gisaid_t dataframe is all unique combos of dates, countries, lineages
    stopifnot('Number of rows of dataframe not equal to number of unique variant-country-days'= nrow(gisaid_t) == n_codes*n_lins*n_dates)
   
    
    # Replace n_seq with 0 if na
    gisaid_t$n_seq[is.na(gisaid_t$n_seq)]<-0
    # Adds in country and eweek, and total number of sequences
    gisaid_t<-gisaid_t%>%
        mutate(
            country = countrycode(code, origin = 'iso3c', destination = 'country.name'),
            eweek = epiweek(collection_date)) %>%
        group_by(code, collection_date) %>%
        mutate(tot_seq = sum(n_seq),
               p_lineage = n_seq/tot_seq,
               p_lineage_se = sqrt(p_lineage*(1-p_lineage)/tot_seq)) %>%
        ungroup() 
    
    gisaid_t$tot_seq[is.na(gisaid_t$tot_seq)]<-0
    
    
    if (VARIANT_TYPE == "pango_lineage"){
        gisaid_t<- gisaid_t %>% 
            rename(pango_lineage = lin)
        print('Dataset with every country-pango-lineage-timepoint is complete')
    }
    if (VARIANT_TYPE == "lineage"){
        gisaid_t<- gisaid_t %>% 
            rename(lineage = lin)
        print('Dataset with every country-lineage-timepoint is complete')
    }
    if (VARIANT_TYPE == "variant"){
        gisaid_t<- gisaid_t %>% 
            rename(variant = lin)
        print('Dataset with every country-variant-timepoint is complete')
    }
 
    return(gisaid_t)
    
}





add_bucketed_lineages<-function(gisaid_t, variants_tracking, COMBINE_BA.4_BA.5){
    # input: variant-country-day dataset and variants we're tracking
    # output: variant-country-dat dataset with additional columns for lineage buckets
    # and variants we're tracking vs all others
            
    gisaid_t <- gisaid_t %>% mutate(
            lineage = case_when( 
                pango_lineage == 'Unassigned' ~ 'Unassigned',
                pango_lineage %in% VOCs ~ pango_lineage,
                (COMBINE_BA.4_BA.5 == TRUE & grepl("BA.4", pango_lineage, fixed = TRUE)) ~ "BA.4/BA.5", # group BA.4 and BA.5 together
                (COMBINE_BA.4_BA.5 == TRUE & grepl("BA.5", pango_lineage, fixed = TRUE)) ~ "BA.4/BA.5",
                pango_lineage %in% variants_tracking ~ pango_lineage,
                grepl("AY", pango_lineage, fixed = TRUE) ~ "B.1.617.2",
                TRUE ~ str_match(pango_lineage, "[a-zA-Z]{1,2}[\\.]?[0-9]*")
            )
        )
    
    if (COMBINE_BA.4_BA.5 == TRUE){
        stopifnot("BA.4 and BA.5 aggregation not adding up" = sum(gisaid_t$pango_lineage== "BA.4") + 
                      sum(gisaid_t$pango_lineage == "BA.5") == sum(gisaid_t$lineage == "BA.4/BA.5"))
    }
    gisaid_t$lineage[is.na(gisaid_t$lineage)]<-'other'
    unique_lineages<-unique(gisaid_t$lineage)
    
    print(paste0('Bucketing lineages takes us from ', length(unique(gisaid_t$pango_lineage)), ' to ', length(unique_lineages)))
    
    # Add another column that is variants were tracking vs all others + the dominant lineage early on and later
    # get code and dominant lineage df for last week
    gisaid_prev_last_wk<-gisaid_t%>% group_by(code) %>%
        filter(collection_date >= (max(collection_date) - days(14))) %>%
        mutate(tot_seq = n())  %>%
        group_by(code, lineage)%>%
        summarise(n_seq = n(),
                  tot_seq = max(tot_seq),
                  p = n_seq/tot_seq) %>% 
        ungroup() %>% 
        group_by(code) %>% 
        mutate(
                  rank_lins = rank(-p)) %>%
        ungroup() %>%
        group_by(code) %>%
        summarise(
            dom_lin_last_wk = lineage[which.max(p)]
            #second_dom_lin_last_wk = pango_lineage[rank_lins ==2]
        )
    # same for the first week 
    gisaid_prev_first_wk<-gisaid_t%>% group_by(code) %>%
        filter(collection_date <= (min(collection_date) + days(7))) %>%
        mutate(tot_seq = n())  %>%
        group_by(code, lineage)%>%
        summarise(n_seq = n(),
                  tot_seq = max(tot_seq),
                  p = n_seq/tot_seq) %>%
        ungroup() %>%
        group_by(code) %>%
        summarise(
            dom_lin_first_wk = lineage[which.max(p)]
        )
        
        # Join to larger dataset
        gisaid_t<-left_join(gisaid_t, gisaid_prev_first_wk, by = "code")
        gisaid_t<-left_join(gisaid_t, gisaid_prev_last_wk, by = "code")
    
        gisaid_t<-gisaid_t %>% group_by(code, pango_lineage) %>% 
            mutate(
                variant = case_when(
                lineage %in% variants_tracking ~ lineage,
                #lineage == dom_lin_first_wk ~ lineage,
                #lineage == dom_lin_last_wk ~ lineage,
                !lineage %in% variants_tracking ~ 'other'
            )
        ) %>% ungroup()
    
    return(gisaid_t)
}



get_7d_avgs<-function(gisaid_t, VARIANT_TYPE){
    if (VARIANT_TYPE == 'pango_lineage'){
        gisaid_t<-gisaid_t%>%group_by(code, pango_lineage) %>% 
            rename(lin = pango_lineage)
    }
    if (VARIANT_TYPE == 'lineage'){
        gisaid_t<-gisaid_t%>%group_by(code, lineage) %>% 
            rename(lin = lineage)
    }
    if (VARIANT_TYPE == 'variant'){
        gisaid_t<-gisaid_t%>%group_by(code, variant) %>% 
            rename(lin = variant)
    }
    
    
    gisaid_t_calib<-gisaid_t %>% filter(period == "calibration") %>%
        group_by(code, lin) %>% 
        mutate(
            week = floor(t/7) +1,
            seq_7_d= zoo::rollapply(n_seq, 7,sum, align = "center", fill = NA),
            tot_seq_7_d = zoo::rollapply(tot_seq, 7,sum, align = "center", fill = NA),
            p_lineage_week = seq_7_d/tot_seq_7_d,
            p_lineage_week_se = sqrt(p_lineage_week*(1-p_lineage_week)/tot_seq_7_d)) %>% 
        ungroup () %>% 
        group_by(code, lin, week) %>% 
        mutate(
            mid_week_date = ymd(median(collection_date)),
            mid_week_p_lineage = ifelse(
                mid_week_date == collection_date, p_lineage_week, NA),
            mid_week_p_lineage_se = ifelse(
                mid_week_date == collection_date, p_lineage_week_se, NA))%>% 
                ungroup() %>% 
        select(collection_date, code, lin, seq_7_d, tot_seq_7_d, p_lineage_week, p_lineage_week_se, week, mid_week_date, mid_week_p_lineage, mid_week_p_lineage_se)


    gisaid_t<- left_join(gisaid_t, gisaid_t_calib, by = c("collection_date", "code", "lin"))
    
    
     if (VARIANT_TYPE == 'pango_lineage'){
         gisaid_t <- gisaid_t %>% 
             rename(pango_lineage = lin)
     }
    if (VARIANT_TYPE == 'lineage'){
        gisaid_t <- gisaid_t %>% 
            rename(lineage = lin)
    }
    if (VARIANT_TYPE == 'variant'){
        gisaid_t <- gisaid_t %>% 
            rename(variant = lin)
    }
    
    
    return(gisaid_t)
}

#  Clean & merge OWID data 
clean_owid<-function(owid_raw, TODAY_DATE, duration ){
    owid<-owid_raw %>% filter(
        date >= (ymd(TODAY_DATE) - days(duration) -1)) %>%
        mutate(date = ymd(date))
    
    #select column names
    owid<-owid%>%select(date,iso_code, new_cases, population,
                        people_fully_vaccinated) # option to add more here!
    
    
    # Fill in the missing country-days just in case there are some
    iso_code <-unique(owid$iso_code)
    date<-seq(from = (ymd(TODAY_DATE) - days(duration)), to = ymd(TODAY_DATE), by = 1)
    date_codes_owid<-expand_grid(date, iso_code)
    owid<-left_join(date_codes_owid,owid, by = c("iso_code", "date"))
    
    # Make NAs 0s for cases
    owid$new_cases[is.na(owid$new_cases)]<-0
    
    # add cases per 100k and smoothed cases
    owid<- owid %>% group_by(iso_code) %>%
        mutate(
            new_cases_7d_avg = zoo::rollmean(new_cases, fill = NA, 7, align = "right"),
            cases_per_100k = 1e5*new_cases_7d_avg/population,
            raw_cases = new_cases)
    
    

    return(owid)
    
}


# Identify a country for each variant
get_country_variant_summary_stats <- function(variant_t){
    country_variant_df<-variant_t %>% 
        filter( period == 'calibration') %>%
        filter(variant %in% variants_tracking) %>% 
        group_by(code, variant, country) %>%
        summarise(
            tot_seq_of_variant = sum(n_seq, na.rm = T),
            last_var_prev = mean(tail(p_lineage_week,10), na.rm = T),
            mean_lag_time_variant = first(tail(mean_lag_time_variant,4)),
            seq_last_week = first(tail(seq_7_d,4))
        ) %>% ungroup()%>% 
        group_by(variant) %>%
        arrange(desc(last_var_prev), by_group = TRUE) %>% 
        mutate(
            rank_var_prev = rank(-last_var_prev),
            rank_tot_seq_of_variant = rank(-tot_seq_of_variant),
            rank_seq_last_week = rank(-seq_last_week),
            sum_ranks = rank_var_prev + rank_tot_seq_of_variant + rank_seq_last_week
        ) %>% ungroup()
    return(country_variant_df)
}

pick_country_for_each_variant<-function(country_variant_df){
    # Quick and dirty way to find the best country for each variant (needs revising)
    variant_df <- country_variant_df %>%
        group_by(variant) %>%
        filter( sum_ranks == min(sum_ranks)) %>%
        mutate( code_lineage = as.factor(paste0(code,'_', variant)))
    return(variant_df)
}



# Run everything from start to finish-----------------------------------------

# Give the model an empty list of country codes, so it includes all countries with 
# sequences within 90 days of the reference date 
code_list <-c()
metadata_subset<-clean_metadata(metadata, FIRST_DATE, TODAY_DATE, LAST_DATE, duration, code_list)
print('Cleaned metadata')

# Add columns that collapse pango-lineages into "lineages" and "variants"
metadata_subset<-add_bucketed_lineages(metadata_subset, variants_tracking, COMBINE_BA.4_BA.5)
print('Added lineage buckets and variants were tracking to linelist data')


# Generate datasets with every country-pangolineage-timepoint combo
print('Generate dataset with every country-pango-lineage-timepoint combo')
gisaid_t<-get_lineage_country_day(metadata_subset, VARIANT_TYPE = 'pango_lineage', LAST_DATE = TODAY_DATE, OBSERVATION_THRESHOLD)

print('Generate dataset with every country-lineage-timepoint combo')
lineage_t<-get_lineage_country_day(metadata_subset, VARIANT_TYPE = 'lineage', LAST_DATE = TODAY_DATE, OBSERVATION_THRESHOLD)
unique_lineages<-unique(lineage_t$lineage)
print(paste0('Bucketed ', length(unique_lineages), ' lineages that exceed the observation threshold'))

print('Generate dataset with variants tracking + dominant lineages for each country-timepoint')
variant_t<-get_lineage_country_day(metadata_subset, VARIANT_TYPE = 'variant', LAST_DATE = TODAY_DATE, OBSERVATION_THRESHOLD)

unique_variants<-unique(variant_t$variant)
print(paste0('Variants were tracking: ', unique_variants))

# quick plot to check
variant_t %>% filter(country == "United States") %>% ggplot() + geom_area(aes(x = collection_date, y = p_lineage, fill = variant))+
    facet_wrap(~country) + theme_bw() + xlab('Date') + ylab('Variant prevalence')


# Get 7 day rolling averages for each of the 3 aggregated datasets
print('Get 7 day rolling averages of lineage prevalence and total sequences')
#gisaid_t<-get_7d_avgs(gisaid_t, VARIANT_TYPE = 'pango_lineage')
#lineage_t<-get_7d_avgs(lineage_t,  VARIANT_TYPE = 'lineage')
variant_t<-get_7d_avgs(variant_t, VARIANT_TYPE = 'variant')


print('Save files')
if (USE_CASE == "local"){
  write.csv(gisaid_t, '../data/processed/gisaid_t_2022-07-01.csv', row.names = F)
  write.csv(lineage_t, '../data/processed/lineage_t_2022-07-01.csv', row.names = F)
  write.csv(variant_t, '../data/processed/variant_t_2022-07-01.csv', row.names = F)
}


