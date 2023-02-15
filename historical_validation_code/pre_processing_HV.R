# Author: Kaitlyn Johnson

# Date initiated: 06-21-2022


# This script uses the list of reference dates and corresponding data sets and 
# generates the lineage_t dataset needed to run the multicountry_model


# This requires first running all of the datasets and getting the unique set of
# lineages for each, and then generating datasets aggregated to that lineage level

# Note: To regenerate the data as is presented in the manuscript, one must have the 
# GISAID metadata from those days, rather than just filtering for submission dates
# before the last submission date in each reference date. This is because the 
# pangolineage assignments get changed as the pangolineage model (which assigns 
# pangolineages to sequences) updates. 



rm(list = ls())

# Set global parameters-------------------------------------------------------

variants_list<- c('BA.2', 'BA.1','BA.2.12.1', 'BA.4', 'BA.5') 
sub_pango_variants<-c()
variants_tracking<-c(variants_list, sub_pango_variants)
VOCs<- c('B.1.1.7', 'P.1', 'AY', 'C.37', 'B.1.621', 'B.1.1.529', 'B.1351', 'B.1.429')
duration<-90 # time window to run model
OBSERVATION_THRESHOLD<-50 # total global number of sequences needed for a lineage to be fit to the model (go into lineage_t)


# Install Libraries-------------------------------------------------------------


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



# Load in reference date data & set country list -----------------------
REFERENCE_DATA_PATH <- '../data/processed/reference_data_df.csv'
reference_data_df <- read_csv(REFERENCE_DATA_PATH)

# Filter the dates we want to actually use
dates_list <-ymd(c("2022-04-30", "2022-05-16", "2022-05-27",
               "2022-06-04", "2022-06-27", "2022-07-01"))
reference_data_df <- reference_data_df %>% filter(reference_date %in% dates_list)

# leave this empty since we will run for all countries
country_list <-c()

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
    
    #6. Filter dataset so collection dates are from the past "duration" days
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
    
    # Make sure data isn't too far out of date
    most_recent_submission_date<-max(metadata_subset$submission_date, na.rm = T)
    print(paste0('Most recent submission from ', most_recent_submission_date))
    stopifnot('Most recent submission is more than 5 days ago' = as.integer(ymd(TODAY_DATE) - most_recent_submission_date ) >0)
    return(metadata_subset)
}


# Data aggregation
get_lineage_country_day<-function(metadata_subset, VARIANT_TYPE, LAST_DATE, OBSERVATION_THRESHOLD ){
    # input: metadata (line list) subsetted to dates we want
    # VARIANT_TYPE: level you want to agrgegate to (pango_lineage, bucketed lineages, variants we're tracking)
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
    
    
    if (VARIANT_TYPE != 'variant'){ # We don't want every country to have every other country's dominant lineage 
        # Makes sure every lineage in every country has a value
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
        # add a column for the period here 
        
        
        gisaid_t<-left_join(list_of_variant_country_days, gisaid_t, by = c("collection_date", "code", "lin"))
        gisaid_t<-left_join(gisaid_t, date_t, by = c("collection_date" = "collection_dates"))
        # check that the length of the gisaid_t dataframe is all unique combos of dates, countries, lineages
        stopifnot('Number of rows of dataframe not equal to number of unique variant-country-days'= nrow(gisaid_t) == n_codes*n_lins*n_dates)
    }
    # For the variant bucketing, we want to make sure that within a country every variant has a value on each date. 
    if (VARIANT_TYPE == 'variant'){ # We don't want every country to have every other country's dominant lineage 
        # Makes sure every lineage in every country has a value
        variant_t<-c()
        for (i in 1:n_codes){
            gisaid_i<-gisaid_t%>%filter(code == unique_codes[i])
            gisaid_orig<-gisaid_i
            n_country_lins <- length(unique(gisaid_i$lin))
            # Add all dates we want, here up til the max collection date in GISAID (but could also add until current date for nowcasting)
            gisaid_i<-gisaid_i %>%complete(collection_date =seq(min(gisaid_i$collection_date), ymd(LAST_DATE)-days(1), by = "days"))
            n_dates<-length(unique(gisaid_i$collection_date))
            list_of_variant_days<-gisaid_i%>%
                complete(lin, collection_date)%>%
                select(lin, collection_date) %>% 
                filter(!is.na(lin)) %>% 
                mutate(period = case_when(
                    collection_date <= max(gisaid_orig$collection_date) ~ 'calibration',
                    collection_date > max(gisaid_orig$collection_date) ~ 'nowcast'
                ))
            gisaid_i<-left_join(list_of_variant_days, gisaid_i, by = c("collection_date", "lin"))
            gisaid_i<-left_join(gisaid_i, date_t, by = c("collection_date" = "collection_dates"))
            gisaid_i$code[is.na(gisaid_i$code)]<-unique_codes[i]
            stopifnot('Number of rows of country-level dataframe not equal to number of unique variant-days'= nrow(gisaid_i) == n_country_lins*n_dates)
            variant_t <- bind_rows(variant_t, gisaid_i)
        }
        gisaid_t<-variant_t
    }
    
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




get_lineage_country_day_for_validation<-function(metadata_subset, VARIANT_TYPE = 'lineage', LAST_DATE, lineage_list){
    # input: metadata (line list) subsetted to dates we want
    # LAST_DATE: date we want to cast out to, for mapping to model prediction. Could be the last day we observe data(max(collection_date)), or could be today's date, or a week into the future, etc.
    # lineage_list: list of lineages we want to include 
    # output: dataset with every possible combination of pango lineage, country,
    # and date for time period of interest
    
    
    if (VARIANT_TYPE == 'lineage'){
        # Find lineages that are below observation threshold
    
        # Combine these into others    
        gisaid_t<- metadata_subset %>%
            mutate(lineage = case_when(
                       (!lineage %in% lineage_list ) ~'other',
                       lineage == "Unassigned" ~ 'other',
                       TRUE ~ lineage)) %>% 
            group_by(code, collection_date, lineage) %>%
            summarise(n_seq = n()) %>%
            rename(lin = lineage) # lin is generic term for VARIANT_TYPE
        # Join the lineage list to GISAID
        lin_df <- data.frame(lineage_list)
        gisaid_t<- gisaid_t %>% full_join(lin_df, by = c("lin" = "lineage_list"))
        # Arbitrarily add 0 sequences to a country and date so that it gets "completed" in the next section 
        if (any(is.na(gisaid_t$code) & is.na(gisaid_t$collection_date))){
            gisaid_orig<- gisaid_t
            gisaid_t$n_seq[is.na(gisaid_orig$n_seq) & is.na(gisaid_orig$collection_date)]<-0
            gisaid_t$code[is.na(gisaid_orig$n_seq) & is.na(gisaid_orig$collection_date)]<-gisaid_orig$code[1]
            gisaid_t$collection_date[is.na(gisaid_orig$n_seq) & is.na(gisaid_orig$collection_date)] <-gisaid_orig$collection_date[1]
        }
        n_lins<-length(unique(gisaid_t$lin))
    }
    
    
    gisaid_t<- gisaid_t %>%select(code, collection_date, n_seq, lin) %>%
        ungroup()
    
    # Complete variant-country-days
    unique_codes <- unique(gisaid_t$code)
    n_codes <- length(unique_codes)
    collection_dates<-seq(from = ymd(min(gisaid_t$collection_date, na.rm = T)), to = ymd(LAST_DATE)-days(1), by = "days")
    t <- seq(from = 1, to = length(collection_dates), by = 1)
    n_dates<-length(collection_dates)
    date_t<-data.frame(t,collection_dates)
    
    
    if (VARIANT_TYPE == 'lineage'){ # We don't want every country to have every other country's dominant lineage 
        # Makes sure every lineage in every country has a value
        # Option have this go to current date if we want to nowcast 
        gisaid_orig<-gisaid_t
        gisaid_t<-gisaid_t %>%complete(collection_date =seq(min(gisaid_t$collection_date, na.rm = T), ymd(LAST_DATE)-days(1), by = "days"))
        list_of_variant_country_days<-gisaid_t%>%
            complete(lin, code, collection_date)%>%
            select(lin, code, collection_date) %>%
            filter(!is.na(lin)) %>% 
            filter(!is.na(code)) %>% 
            mutate(period = case_when(
                collection_date <= max(gisaid_orig$collection_date) ~ 'calibration',
                collection_date > max(gisaid_orig$collection_date) ~ 'nowcast'
            ))
        # add a column for the period here 
        
        
        gisaid_t<-left_join(list_of_variant_country_days, gisaid_t, by = c("collection_date", "code", "lin"))
        gisaid_t<-left_join(gisaid_t, date_t, by = c("collection_date" = "collection_dates"))
        # check that the length of the gisaid_t dataframe is all unique combos of dates, countries, lineages
        stopifnot('Number of rows of dataframe not equal to number of unique variant-country-days'= nrow(gisaid_t) == n_codes*n_lins*n_dates)
    }
    
    
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
    
    
    
    
    if (VARIANT_TYPE == "lineage"){
        gisaid_t<- gisaid_t %>% 
            rename(lineage = lin)
        print('Dataset with every country-lineage-timepoint is complete')
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
   
    
    # If the pango lineage doesn't fall into any of these buckets, then we label it as other
    # NOTE: When we get to the point where we are identify variants algorithmically (not setting them)
    # we will want to keep the pango lineages of the variants we havent already bucketed.
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
                lineage == dom_lin_first_wk ~ lineage,
                lineage == dom_lin_last_wk ~ lineage,
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


# Run loop through reference dates to find consensus lineage list------------
# We must do this since we need the same set of lineages in all of the model runs for comparison

lineage_list<- c()
for (i in 1:nrow(reference_data_df)){
    THIS_REF_DATE<- as.character(reference_data_df$reference_date[i])
    this_metadata <- read.csv(paste0('../data/raw/', THIS_REF_DATE, '_metadata.csv'))
    LAST_DATE <- as.character(max(reference_data_df$last_submission_date[i]))
    FIRST_DATE <- as.character(ymd(LAST_DATE)- days(duration))
    
    metadata_subset<- clean_metadata(this_metadata, FIRST_DATE, LAST_DATE, LAST_DATE, duration, country_list)
    
    metadata_subset<-add_bucketed_lineages(metadata_subset, variants_tracking, COMBINE_BA.4_BA.5)

    
    lineage_t<-get_lineage_country_day(metadata_subset, VARIANT_TYPE = 'lineage', LAST_DATE, OBSERVATION_THRESHOLD)
    unique_lineages<-unique(lineage_t$lineage)
    print(paste0('There are ', length(unique_lineages), ' lineages from data from ', THIS_REF_DATE))
    
    lineage_list <- c(lineage_list, unique_lineages)
    
}
# List of lineages that need to be a part of every model fit 
lineage_list <- unique(lineage_list)
print(paste0('There are ', length(lineage_list), ' lineages from data across all reference datasets'))


# Run again, this time, making sure we have all of the lineages that are in any of the datasets
for (i in 1:nrow(reference_data_df)){
    THIS_REF_DATE<- as.character(reference_data_df$reference_date[i])
    this_metadata <- read.csv(paste0('../data/raw/', THIS_REF_DATE, '_metadata.csv'))
    LAST_DATE <- as.character(max(reference_data_df$last_submission_date[i]))
    FIRST_DATE <- as.character(ymd(LAST_DATE)- days(duration))
    
    metadata_subset<- clean_metadata(this_metadata, FIRST_DATE, LAST_DATE, LAST_DATE, duration, country_list)
    
    metadata_subset<-add_bucketed_lineages(metadata_subset, variants_tracking, COMBINE_BA.4_BA.5)
    print('Added lineage buckets and variants were tracking to linelist data')
    
    print('Generate dataset with every country-lineage-timepoint combo')
    lineage_t<-get_lineage_country_day_for_validation(metadata_subset, VARIANT_TYPE = 'lineage', LAST_DATE, lineage_list)
    unique_lineages<-unique(lineage_t$lineage)
    stopifnot('lineages do not all match across datasets'= length(unique_lineages) == length(lineage_list))
    print(paste0('There are ', length(unique_lineages), ' lineages from data from ', THIS_REF_DATE, 'now that we have set a consensus list of lineages'))
    
    # Get 7d averages
    lineage_t<-get_7d_avgs(lineage_t,  VARIANT_TYPE = 'lineage')
    
    # Write to a file
    write.csv(lineage_t, paste0('../data/processed/validation_data/lineage_t_', THIS_REF_DATE, '.csv'))
    
    
}

# Write the reference data list to be used in all the subsequent scripts
write.csv(reference_data_df, '../data/processed/validation_data/reference_data_used.csv')


# Make a lineage_t for the most recent data we are using with no cut off time----------

RECENT_DATA_DATE<- "2022-07-01"
this_metadata <- read.csv(paste0('../data/raw/', RECENT_DATA_DATE, '_metadata.csv'))
LAST_DATE <- as.character(max(reference_data_df$last_submission_date[reference_data_df$reference_date==RECENT_DATA_DATE]))
FIRST_DATE <- as.character(ymd(LAST_DATE) - days(365))
duration <- as.numeric(ymd(LAST_DATE) - ymd(FIRST_DATE))
metadata_subset<- clean_metadata(this_metadata, FIRST_DATE, LAST_DATE, LAST_DATE, duration, country_list)

metadata_subset<-add_bucketed_lineages(metadata_subset, variants_tracking, COMBINE_BA.4_BA.5)
print('Added lineage buckets and variants were tracking to linelist data')

print('Generate dataset with every country-lineage-timepoint combo')
lineage_t<-get_lineage_country_day_for_validation(metadata_subset, VARIANT_TYPE = 'lineage', LAST_DATE, lineage_list)
unique_lineages_for_comp<-unique(lineage_t$lineage)
stopifnot('Lineages must be the same for this data as for the historical dataset'= length(unique_lineages_for_comp) == length(lineage_list))
print(paste0('There are ', length(unique_lineages_for_comp), ' lineages from data from ', RECENT_DATA_DATE))


# This is the data from 2022-07-01 that we will use as a comparison
write.csv(lineage_t, paste0('../data/processed/validation_data/lineage_t_for_comp_', RECENT_DATA_DATE,'.csv'))
