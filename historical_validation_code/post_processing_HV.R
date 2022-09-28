# Author: Kaitlyn Johnson
# Variant tracker post-processing

# This script loads in the output from the hierarchical multinomial model (either
# multi-country or single country) + the outputs from the tracker_gisaid_processing_script.

# We start by identifying, for each variant we're tracking, which country to display it in
# Then we make visuals for each variant-country combo and package them for the Rmd

rm(list = ls())
USE_CASE = Sys.getenv("USE_CASE")
if(USE_CASE == ""){
    USE_CASE<-'local'
}

# Global variables-----------------------------------------

FINAL_RUN<-FALSE # set to TRUE if this is the run that will generate the weekly report, else set as false to avoid continuously adding to the dataset
FIRST_REPORT<-TRUE # if this is the first report, then we need to make the historical dataset new, else load in historical data from previous reports
Rt_METHOD<-'both' #EpiEstim' # or 'Ana method' or 'both' (EpiEstim uses built in package, Ana method uses MCMC)
ACCOUNT_FOR_DELAY<-FALSE # set to true if we want to shift back Rt to correspond to when those individuals likely got infected,
MULTICOUNTRY_MODEL<-TRUE # set to true if the output is coming from the multicountry model, to false if coming fromt the single country model
RT_INPUT <- 'model observed cases' # 'model observed cases' # 'model mean cases'
MEAN_GI <- 5.8 # need to update
STD_GI <- 2.3
n_draws<- 100 # number of draws from posterior 
nsim <- 1000 # number of MCMC draws 


# Install Libraries--------------------------------------

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
library(EpiEstim) # Estimate Rt
library(lme4) # logistic regression


# Load in data -----------------------------------------------------------------

if (USE_CASE == 'local'){
    VARIANT_T_PATH <- "../data/processed/variant_t.csv" # processed data by variants we're tracking
    VARIANT_COUNTRY_PATH <- "../data/processed/variant_df.csv" # list of variants we're tracking
    MODEL_OUTPUT_PATH <- "../data/processed/model_output.rds" # hierarchical multinomial model output as a stanfit object
    MODEL_PRED_P_PATH <- "../data/processed/model_pred_p.csv" # large data table of variant proportion distributions for S. Africa currently
    #MODEL_PRED_P_PATH <- "../data/processed/model_pred_p_USA.csv" # large data table of variant proportion distributions for US 
    MODEL_R_EST_PATH <- "../data/processed/r_distrib.csv" # large data table of relative growth estimates for S.Africa currently
    #MODEL_R_EST_PATH <- "../data/processed/r_distrib_USA.csv" # large data table of relative growth estimates for US
    HIST_VAR_COUNTRY_PATH <- "../data/output/variant_df_history.csv"
    HIST_GROWTH_ADV_PATH <- "../data/output/growth_adv_df_history.csv"
    HIST_VARIANT_T_ALL_PATH <- "../data/output/variant_t_all_history.csv"
    HIST_VAR_T_LONG_PATH <- "../data/output/var_t_long_history.csv"
    HIST_LAG_TIME_LONG_PATH <- "../data/output/lag_time_long_history.csv"
    MODEL_PATH <- "stancode/Rt_estimation.stan"
    if (MULTICOUNTRY_MODEL == TRUE){
        MODEL_PRED_P_PATH <- "../data/processed/multicountry_model_pred_p.csv" # all draws, all variants, all countries, all time points 
        MODEL_R_EST_PATH <- "../data/processed/multicountry_r_distrib.csv" # all draws, all variants, all countries, growth estimates
        # Other multicountry outputs are over all countries, and can directly be loaded into report
    }

}
#Domino
if (USE_CASE == 'domino'){
    VARIANT_T_PATH<-"/mnt/data/processed/variant_t.csv" # processed data by variants we're tracking
    VARIANTS_COUNTRY_PATH <- "/mnt/data/processed/variant_df.csv" # list of variants we're tracking
    MODEL_OUTPUT_PATH<-"/mnt/data/processed/model_output.rds" # hierarchical multinomial model output
    MODEL_PRED_P_PATH <- "/mnt/data/processed/model_pred_p.csv" # large data table of variant proportion distributions
    MODEL_R_EST_PATH <- "/mnt/data/processed/r_distrib.csv" # large data table of relative growth estimates
    HIST_VAR_COUNTRY_PATH <- "/mnt/data/output/variant_df_history.csv"
    HIST_GROWTH_ADV_PATH <- "/mnt/data/output/growth_adv_df_history.csv"
    HIST_VARIANT_T_ALL_PATH <- "/mnt/data/output/variant_t_all_history.csv"
    HIST_VAR_T_LONG_PATH <- "/mnt/data/output/var_t_long_history.csv"
    HIST_LAG_TIME_LONG_PATH <- "/mnt/data/output/lag_time_long_history.csv"
}



variant_t<-read_csv(VARIANT_T_PATH)
model_pred_p<-read_csv(MODEL_PRED_P_PATH)
r_distrib<-read_csv(MODEL_R_EST_PATH)
variant_df<-read_csv(VARIANT_COUNTRY_PATH)
variants_tracking<-variant_df$variant
print(variants_tracking)

# For now, remove BA.2
variant_df <- variant_df %>% filter(variant != "BA.2", variant != "BA.1")
variant_df$reference_lineage <- NA # to be filled in

if (FIRST_REPORT == FALSE){
    variant_df_history <- read.csv(HIST_VAR_COUNTRY_PATH)
    growth_adv_df_history <- read.csv(HIST_GROWTH_ADV_PATH)
    variant_t_all_history <- read.csv(HIST_VARIANT_T_ALL_PATH)
    var_t_long_history <- read.csv(HIST_VARIANT_T_LONG_PATH)
    lag_time_long_history <- read.csv(HIST_LAG_TIME_LONG_PATH)
}
if (FIRST_REPORT == TRUE){
    variant_df_history <- c()
    growth_adv_df_history <- c()
    variant_t_all_history <- c()
    var_t_long_history <- c()
    lag_time_long_history <- c()
}


# Get countries and variants for viz -------------------------------------------
codes<-variant_df$code
countries<-variant_df$country
lineages<-variant_df$variant
code_lineages<-variant_df$code_lineage



# Loop through each variant-country we are visualizing-------------------

growth_adv_df<-c()
variant_t_all<-c()
var_t_long<-c()
lag_time_long<-c()
coverage_df<-c()


#i<-3
#for (i in 1:length(code_lineages)){
for(i in 1:length(code_lineages)){
    # Designate the charactersistics for this country & variant
    this_code <- codes[i]
    this_country <- countries[i]
    this_variant <- lineages[i]
    this_variant_t<-variant_t %>% filter(code == this_code)
    these_variants_tracking <- unique(this_variant_t$variant)
    this_model_pred_p <- model_pred_p %>% filter( country == this_country)
    these_variants_est<- unique(this_model_pred_p$lineage)
    this_r_distrib <- r_distrib %>% filter (country == this_country & lineage == this_variant)
    reference_lineage <- this_r_distrib$dom_lin[1]
    variant_df$reference_lineage[i]<- reference_lineage
 # Growth advantage of each variant compared to dominant variant over the generation interval (5 days)         
    growth_adv_w_CIs<- round(100*(exp(MEAN_GI*(quantile(this_r_distrib$r, probs = c(0.025, 0.5, 0.975))))-1),2)
    ub <- growth_adv_w_CIs[[3]]
    lb <- growth_adv_w_CIs[[1]]
    median <- growth_adv_w_CIs[[2]]
    this_growth_adv_df<-data.frame(this_variant, this_code, this_country, median, ub, lb)
    this_growth_adv_df <- this_growth_adv_df %>% rename(variant = this_variant,
                                                        code = this_code,
                                                        country = this_country)
    
    print(paste0('Estimated growth advantage of ', this_variant, ' compared to ', reference_lineage, ' in ', this_country, ' is: ', growth_adv_w_CIs[2], '% [95% PI: ', growth_adv_w_CIs[1] ,'% , ', growth_adv_w_CIs[3], '%]'))
    stopifnot('Fitness advantage either missing or may be relevant to the wrong variant, check prior code for dominant lineage' = (growth_adv_w_CIs[2]>0 | 1>is.na(growth_adv_w_CIs)))
    
    

# Tidy data from model output -------------------------------------------

    
    # Tidy model output by aggregating lineages across a draw
    # so that we get the distribution of estimated prevalence
    # of the variant we're tracking vs all others at each time point
    
    p_distrib_t <- this_model_pred_p %>% 
        mutate( # collapses multinomial to only variants we're tracking vs all others
                # we may want to either 1. keep more lineages in this or 
                # collapse from lineage to variant here!
                variant = case_when(
                    lineage %in% these_variants_tracking ~ lineage,
                    !lineage %in% these_variants_tracking ~ 'other')
            ) %>% 
        group_by(t, draw, variant) %>%
        summarize(
               p_est = sum(p_hat), # p_hat for that draw (model predicted true mean)
               Y_tilde = sum(Y_tilde), # model predicted sequence numbers for that draw
               pred_obs_prev = sum(pred_obs_prev), # model predicted observed prevalence for that draw
               pred_7d_obs_prev = sum(pred_7d_obs_prev), #  uses 7 day averages of model predicted observed prevalences
               obs_7d_prev = sum(obs_7d_prev),
               date = max(date),
               code = this_code,
               ) %>%
        ungroup() %>% 
        mutate(var_draw = paste0(variant, '_', draw),
                             date = ymd(date))
    
    this_countrys_lineages_from_model<-unique(p_distrib_t$variant)
    n_draws<-length(unique(p_distrib_t$draw))
    
    stopifnot('Long model output is not the right size, check dates'= nrow(p_distrib_t) == n_draws*length(this_countrys_lineages_from_model)*length(unique(this_variant_t$t)))
    
    
    # check
    p_test <-p_distrib_t %>%group_by(t, draw) %>%
        # Want these to 
        summarise(p_test = sum(p_est), # this is the sum of the p_hats
                  model_obs_p = sum(pred_obs_prev),
                  pred_7d_obs_p = sum(pred_7d_obs_prev))
    stopifnot('All prevalences dont add to 1' = (round(max(p_test$p_test, na.rm = T),4) ==1 & round(min(p_test$p_test, na.rm = T),4) ==1))
    stopifnot('All prevalences dont add to 1' = (round(max(p_test$model_obs_p, na.rm = T),4) ==1 & round(min(p_test$model_obs_p, na.rm = T),4) ==1))
    stopifnot('All prevalences dont add to 1' = (round(max(p_test$pred_7d_obs_p, na.rm = T),4) ==1 & round(min(p_test$pred_7d_obs_p, na.rm = T),4) ==1))
    


    # join your long dataframe with all draws to the variant stats for that country-variant
    # leaving as a left join for now, but use inner join if you want to eliminate the variants we're tracking that are not in the countries multinomial output
    this_var_t_long <- p_distrib_t %>% left_join(this_variant_t, by = c("date" = "collection_date", "code", "variant", "t"))
    
    # Plot to check that prevalences look right 
    p_distrib_t %>% ggplot() +  geom_line(aes(x = date, y = p_est, group = var_draw, color = variant), alpha = 0.1)+
        xlab('Collection Date') + ylab('Estimated prevalence') + theme_bw() 
    
    # check that no p_ests are NAs
    stopifnot('NAs for estimated prevalence' = any(!is.na(p_distrib_t$p_est)))
    
    
    
    # Put summary (mean and 95% PI) into variants_t_subset
    p_summary_t <- p_distrib_t %>%
        group_by(variant, date, t) %>%
        summarise(
            obs_7d_prev = max(obs_7d_prev), # should all be the same
            obs_prev = max(obs_prev),
            median_prev = quantile(p_est, 0.5, na.rm = T),
            upper_prev = quantile(p_est, 0.975, na.rm = T),
            lower_prev = quantile(p_est, 0.025, na.rm = T),
            median_obs_prev = quantile(pred_7d_obs_prev, 0.5, na.rm = T),
            upper_obs_prev = quantile(pred_7d_obs_prev, 0.975, na.rm = T),
            lower_obs_prev = quantile(pred_7d_obs_prev, 0.025, na.rm = T),
            pred_obs_7d_prev_upper_90 = quantile(pred_7d_obs_prev, 0.95, na.rm = T),
            pred_obs_7d_prev_lower_90 = quantile(pred_7d_obs_prev, 0.05, na.rm = T),
            pred_obs_prev_upper_90 = quantile(pred_obs_prev, 0.95, na.rm = T),
            pred_obs_prev_lower_90 = quantile(pred_obs_prev, 0.05, na.rm = T),
            date = max(date, na.rm = T)
        ) %>% 
        ungroup() %>% 
        group_by(variant, date, t) %>% 
        mutate(
            code = this_code,
            is_within_bounds_7d_avg = ifelse( # will use this to calculate the percent coverage (over the whole duration)
                (obs_7d_prev<= pred_obs_7d_prev_upper_90 & obs_7d_prev >= pred_obs_7d_prev_lower_90), 1, 0), 
            is_within_bounds_daily = ifelse( # will use this to calculate the percent coverage (over the whole duration)
                (obs_prev<= pred_obs_prev_upper_90 & obs_prev >= pred_obs_prev_lower_90), 1, 0)
        )
    # Quick coverage estimates just out of curiosity
    this_coverage_df<- p_summary_t %>% group_by(variant) %>% 
        summarise(
            coverage_7d = sum(is_within_bounds_7d_avg, na.rm = T)/sum(!is.na(is_within_bounds_7d_avg)),
            coverage_daily = sum(is_within_bounds_daily, na.rm = T)/sum(!is.na(is_within_bounds_daily)),
        ) %>% ungroup() %>% 
        mutate(
            code = this_code,
            country = this_country,
            method = ifelse(MULTICOUNTRY_MODEL == TRUE, 'multicountry model', 'single country model')
        )
    
    # quick plot
    p_summary_t %>% ggplot() + geom_point(aes(x = date, y = obs_7d_prev, color = variant), size = 0.5) + 
        geom_ribbon(aes(x = date, ymin = pred_obs_prev_lower_90, ymax = pred_obs_prev_upper_90, fill = variant), alpha = 0.1) +
        theme_bw()
    
    
    
    # add summary stats on distributions to variant_t 
    this_variant_t<-left_join(this_variant_t, p_summary_t,  by = c("collection_date" = "date", "code", "variant", "t"))
    stopifnot(' NAs after joining to estimated prevalences to variant summary stats' = any(!is.na(this_variant_t$median_prev)))
    
    # test plot 
    calib <- this_variant_t %>% filter (period == 'calibration')
    nowcast <- this_variant_t %>% filter(period == "nowcast")
    # mean prevalence
    plot <- ggplot() + geom_line(data = calib, aes(x = collection_date, y = median_prev, group = variant, color = variant), size = 1) +
        geom_ribbon(data = calib, aes(x = collection_date, ymin = lower_prev, ymax = upper_prev, group = variant, fill = variant), alpha = 0.2) +
        geom_line(data = nowcast, aes(x = collection_date, y = median_prev, group = variant, color = variant), size = 0.5, alpha = 0.5) +
        geom_ribbon(data = nowcast, aes(x = collection_date, ymin = lower_prev, ymax = upper_prev, group = variant, fill = variant), alpha = 0.05) +
        geom_point(data = calib, aes(x = ymd(mid_week_date), y = mid_week_p_lineage, color = variant ))+
        geom_linerange(data = calib, aes(x = ymd(mid_week_date), ymin = mid_week_p_lineage - mid_week_p_lineage_se, ymax =mid_week_p_lineage + mid_week_p_lineage_se, color = variant)) +
        theme_bw()+ 
        scale_x_date(date_labels = "%m-%d") + theme(axis.text.x = element_text(angle=45, size = 10, hjust =1)) + 
        theme(text = element_text(size = 14))+
        xlab('Date') +
        ylab('Estimated and observed variant prevalence')
    plot
    # observed prevalence
    plot2 <- ggplot() + geom_line(data = calib, aes(x = collection_date, y = median_obs_prev, group = variant, color = variant), size = 1) +
        geom_ribbon(data = calib, aes(x = collection_date, ymin = lower_obs_prev, ymax = upper_obs_prev, group = variant, fill = variant), alpha = 0.2) +
        geom_line(data = nowcast, aes(x = collection_date, y = median_obs_prev, group = variant, color = variant), size = 0.5, alpha = 0.5) +
        geom_ribbon(data = nowcast, aes(x = collection_date, ymin = lower_obs_prev, ymax = upper_obs_prev, group = variant, fill = variant), alpha = 0.05) +
        geom_point(data = calib, aes(x = ymd(mid_week_date), y = mid_week_p_lineage, color = variant ))+
        geom_linerange(data = calib, aes(x = ymd(mid_week_date), ymin = mid_week_p_lineage - mid_week_p_lineage_se, ymax = mid_week_p_lineage + mid_week_p_lineage_se, color = variant)) +
        theme_bw()+ 
        scale_x_date(date_labels = "%m-%d") + theme(axis.text.x = element_text(angle=45, size = 10, hjust =1)) + 
        theme(text = element_text(size = 14))+
        xlab('Date') +
        ylab('7 day observed prevalence')
    plot2
    
    

# Infer cases with a  variant---------------------------- 
    
    # Use distribution of estimate prevalence * observed cases to
    # get case trajectories, 
    
    this_var_t_long<- this_var_t_long %>% group_by(variant, date, draw) %>%
        mutate(
            cases_per_100k_by_variant = cases_per_100k*p_est,
            cases_by_variant = round(raw_cases*p_est), # model predicted mean cases
            cases_per_100k_by_variant_obs = cases_per_100k*p_lineage_week,
            cases_by_variant_obs = round(raw_cases*p_lineage_week),
            cases_per_100k_by_variant_mod_pred = cases_per_100k*pred_7d_obs_prev,
            cases_by_variant_mod_pred = round(raw_cases*pred_obs_prev))# model predicted observed cases
    
    # Put summary (mean and 95% PI) into variants_t_subset 
    these_cases_summary<-this_var_t_long %>% group_by(date, variant) %>%
        summarize(cases_per_100k_by_var_median = quantile(cases_per_100k_by_variant, 0.5, na.rm = T),
                  cases_per_100k_by_var_upper = quantile(cases_per_100k_by_variant, 0.975, na.rm = T),
                  cases_per_100k_by_var_lower = quantile(cases_per_100k_by_variant, 0.025, na.rm = T),
                  # model predicted mean cases
                  raw_cases_med = quantile(cases_by_variant, 0.5, na.rm = T),
                  raw_cases_lb = quantile(cases_by_variant, 0.025, na.rm = T),
                  raw_cases_ub = quantile(cases_by_variant, 0.975, na.rm = T),
                  # model predicted observed cases (based onnumber of sequences collected in that day)
                  raw_cases_mod_pred_med = quantile(cases_by_variant_mod_pred, 0.5, na.rm = T),
                  raw_cases_mod_pred_lb = quantile(cases_by_variant_mod_pred, 0.025, na.rm = T),
                  raw_cases_mod_pred_ub = quantile(cases_by_variant_mod_pred, 0.975, na.rm = T),
                  period = max(period)
        )
    
    this_variant_t<-left_join(this_variant_t, these_cases_summary, by = c("collection_date" = "date", "variant", "period"))
    

    
    cases_calib<- this_variant_t %>% filter(period == 'calibration')
    cases_nowcast<- this_variant_t %>% filter(period == 'nowcast')
    plot<- ggplot() + geom_line(data = cases_calib, aes( x = collection_date, y = raw_cases_mod_pred_med, group = variant, color = variant), size = 1) + 
        geom_ribbon(data = cases_calib, aes(x = collection_date, group = variant, ymin = raw_cases_mod_pred_lb, ymax = raw_cases_mod_pred_ub, group = variant, fill = variant), alpha = 0.2)+
        #geom_line(data = cases_calib, aes(x = collection_date, y = raw_cases), color = "black")+
        #geom_line(data = cases_nowcast, aes( x = collection_date, y = cases_per_100k_by_var_median, group = variant, color = variant), size = 0.5, alpha = 0.5) + 
        #geom_ribbon(data = cases_nowcast, aes(x = collection_date, group = variant, ymin = cases_per_100k_by_var_lower, ymax = cases_per_100k_by_var_upper, group = variant, fill = variant), alpha = 0.05)+
        xlab('Date') + ylab('Inferred cases from model predicted sequences') + theme_bw() + 
        scale_x_date(date_labels = "%m-%d") + theme(axis.text.x = element_text(angle=45, size = 10, hjust =1))+
        theme(text = element_text(size = 14))
        
    plot
    plot2<- ggplot() + geom_line(data = cases_calib, aes( x = collection_date, y = cases_per_100k_by_var_median, group = variant, color = variant), size = 1) + 
        geom_ribbon(data = cases_calib, aes(x = collection_date, group = variant, ymin = cases_per_100k_by_var_lower, ymax = cases_per_100k_by_var_upper, group = variant, fill = variant), alpha = 0.2)+
        geom_line(data = cases_nowcast, aes( x = collection_date, y = cases_per_100k_by_var_median, group = variant, color = variant), size = 0.5, alpha = 0.5) + 
        geom_ribbon(data = cases_nowcast, aes(x = collection_date, group = variant, ymin = cases_per_100k_by_var_lower, ymax = cases_per_100k_by_var_upper, group = variant, fill = variant), alpha = 0.05)+
        xlab('Date') + ylab('Inferred cases per 100k') + theme_bw() + 
        scale_x_date(date_labels = "%m-%d") + theme(axis.text.x = element_text(angle=45, size = 10, hjust =1))+
        theme(text = element_text(size = 14))
    
    plot2
    
    plot2<- ggplot() + geom_area(data = cases_calib, aes( x = collection_date, y = cases_per_100k_by_var_median, group = variant, fill = variant), alpha = 0.8) + 
        geom_area(data = cases_nowcast, aes( x = collection_date, y = cases_per_100k_by_var_median, group = variant, fill = variant), alpha = 0.2) + 
        xlab('Date') + ylab('Inferred cases per 100k') + theme_bw() +scale_x_date(date_labels = "%m-%d") + theme(axis.text.x = element_text(angle=45, size = 10, hjust =1))+
        theme(text = element_text(size = 14))
    plot2
    
    

# Estimate Rt from case trajectories------------------------------------------- 

    if (Rt_METHOD == 'EpiEstim' | Rt_METHOD == 'both'){
        unique_variants<-unique(this_var_t_long$variant)
        draws <- unique(this_var_t_long$draw)
        for (v in 1:length(unique_variants)){
            for (m in 1:length(draws)){
                # this_var_t_long is the sampled stanfit object with other things added to it, 
                # so it contains all of the p hat trajectories for each variant, for this country
                if(RT_INPUT == 'model mean cases'){
                    cases_var_t_date <- this_var_t_long %>% filter(draw == draws[m], variant == unique_variants[v], period == 'calibration') %>% select(t, date, cases_by_variant, population) %>% 
                        rename(cases = cases_by_variant)
                }
                if(RT_INPUT == 'model observed cases'){
                    cases_var_t_date <- this_var_t_long %>% filter(draw == draws[m], variant == unique_variants[v], period == 'calibration') %>% select(t, date, cases_by_variant_mod_pred, population) %>% 
                        rename(cases = cases_by_variant_mod_pred)
                }
                config <- make_config(list(mean_si = MEAN_GI, std_si = STD_GI))
                Rt <- estimate_R(cases_var_t_date$cases, method = "parametric_si", config = config)
                
                # Rename columns
                colnames(Rt$R)[3]  <- "mean"
                colnames(Rt$R)[4]  <- "sd"
                
                # Make dataframe for means and stdev: each col is a draw, rows are days
                if (m == 1) {
                    Rt_means <- data.frame(Rt$R$t_end)
                    Rt_stds <- data.frame(Rt$R$t_end)
                    colnames(Rt_means) <- "time"
                    colnames(Rt_stds)   <- "time"
                }
                colname <- toString(m)
                Rt_means[[colname]] <- Rt$R$mean
                Rt_stds[[colname]]   <- Rt$R$sd
                
            }
            
            Rt_times <- Rt_means$time
            Rt_medians <- rep(0, length(Rt_times))
            Rt_upperbounds  <- rep(0, length(Rt_times))
            Rt_lowerbounds <- rep(0, length(Rt_times))
        
            
            # Summarize Rt estimates
            for (k in 1:nrow(Rt_means)) {
                Rt_aggregated_samples <- c()
                # For each Rt sample, use mean and std to 
                #   sample gamma, summarize all gamma draws
                #   for each day
                
                for (j in 2:length(draws)+1) { # can't use first column because its time!
                    mean <- Rt_means[k, j]
                    std  <- Rt_stds[k, j]
                    var  <- std ** 2
                    b <- mean / var
                    a <- mean * b
                    Rt_aggregated_samples <- append(
                        Rt_aggregated_samples,
                        rgamma(100, shape = a, rate = b)
                    ) # Supplement: https://rockfound.box.com/s/qjyt916q4t2qaq4qm32pnyp8tl9uecsu
                } # main text: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3816335/
                Rt_medians[k] <- quantile(Rt_aggregated_samples, 0.5)
                Rt_upperbounds[k] <- quantile(Rt_aggregated_samples, 0.975)
                Rt_lowerbounds[k] <- quantile(Rt_aggregated_samples, 0.025)
            }
            
            # Put together into an Rt summary
            mean_delay = MEAN_GI + 2 # time to develop symptoms + time to seek testing (2 days)
            time_window = 3.5
            if (ACCOUNT_FOR_DELAY == TRUE){
                dates<-seq(as.Date(min(cases_var_t_date$date) + time_window - mean_delay), as.Date(max(cases_var_t_date$date) - time_window - mean_delay), "days")
            }
            if (ACCOUNT_FOR_DELAY == FALSE){
                dates<-seq(as.Date(min(cases_var_t_date$date) + time_window ), as.Date(max(cases_var_t_date$date) - time_window ), "days")
            }
            dates<- dates[8:length(Rt_medians)] # remove initial values because unreliable
            Rt_summary_v<-data.frame(dates)
            # Name the two Rt method outputs differently so that we can compare
            Rt_summary_v['Rt_medians']<-Rt_medians[8:length(Rt_medians)]
            Rt_summary_v['Rt_lowerbounds']<-Rt_lowerbounds[8:length(Rt_medians)]
            Rt_summary_v['Rt_upperbounds']<-Rt_upperbounds[8:length(Rt_medians)]
            Rt_summary_v['variant']<-unique_variants[v]
            if (v ==1){
                Rt_summary<-Rt_summary_v
            }
            else{
                Rt_summary<-rbind(Rt_summary, Rt_summary_v)
            }
        } # end loop over each variant 
        Rt_summary$dates<-ymd(Rt_summary$dates)
        this_variant_t<- left_join(this_variant_t, Rt_summary, by = c("collection_date" = "dates", "variant"))
        
        plot<- Rt_summary %>%ggplot() + geom_line(aes(x = dates, y = Rt_medians, color = variant))+
            geom_ribbon(aes(x = dates, ymin = Rt_lowerbounds, ymax = Rt_upperbounds, fill = variant), alpha = 0.2) + 
            geom_hline(aes(yintercept = 1), linetype = "dashed") + coord_cartesian(ylim = c(0,5))+
            theme_bw() + xlab('Collection date') + ylab('Rt estimate')
        plot
        stopifnot('Rts didnt join' = any(this_variant_t$Rt_medians >0.5 & this_variant_t$Rt_medians <2))
    } # end if statement for Rt_METHOD == 'EpiEstim'
    
    if (Rt_METHOD == 'Ana method' | Rt_METHOD == 'both'){
        # from Sichuan Rt paper https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008467
        unique_variants<-unique(this_var_t_long$variant)
        draws <- unique(this_var_t_long$draw)
        #stanmodel <- stan_model(MODEL_PATH)  
        
        
       
        
        
        for (v in 1:length(unique_variants)){
            Rt_df <- data.frame()
            for (m in 1:length(draws)){
                
                # Draw the trajectory of cases (C_v(t)) from the posterior from multinomial model ---------------------------------
                
                if(RT_INPUT == 'model mean cases'){ # Pulls from the model estimated mean (independent of number of sequences on a day)
                    cases_var_t_date <- this_var_t_long %>% filter(draw == draws[m], variant == unique_variants[v], period == 'calibration') %>% 
                        rename(cases = cases_by_variant)
                }
                # Need to discuss/figure out which is best here to use phat but have higher uncertainty at low variant prevalence 
                if(RT_INPUT == 'model observed cases'){ # Pulls from the model predicted observed sequences (more variation on days with fewer sequences)
                    cases_var_t_date <- this_var_t_long %>% filter(draw == draws[m], variant == unique_variants[v], period == 'calibration') %>% 
                        rename(cases = cases_by_variant_mod_pred)
                }
                

                t<-cases_var_t_date$t
                cases<-cases_var_t_date$cases
                
                # Setting the generation time distribution ----------------------
                omega <- dgamma(t, shape = MEAN_GI, rate = 1/STD_GI) # this gets passed into log likelihood function 
                n<- nrow(cases_var_t_date)
                
                R<-matrix(nrow = n, ncol = nsim)
                # Define the loglikelihood -----------------------------------------
                loglikelihood <- function(t_int, R_int, cases, omega, n){
                    sumomega <- 0
                    omega2 <-rep(0, n)
                    # Gets the cumulative infectiousness as a function of time since infection onset
                    for (s in 1:t_int){
                        sumomega <- sumomega + omega[1+t_int-s] # indexing starting at 1 in R 
                    }
                    # normalizes it 
                    for (s in 1:t_int){
                        omega2[s] <- omega[s]/sumomega # ** This is inconsistent with C++ code :omega2[s] = omega[t-s][s]/sumomega; but rows of omega are all the same?
                    }
                    # find the sum of the probability of the previous cases weighted by generation time pdf,
                    # which we will multiply by R(t) to get the expected local cases at time t (shouldn't it be t-s though???) 
                    
                    lambda = 0 
                    for (s in 1:t_int){
                        lambda <-  lambda + cases[t_int-s+1]*omega2[s]
                    }
                    # cumulative sum of the infectiousness of previous cases (previous cases weighted by generation time distribution)
    
                    # multiply by Rt to get the poisson rate, then find the density of the observed cases at time t given that poisson rate
                    
                    loglike <- dpois(cases[t_int], lambda = R_int*lambda, log = T)
                    return(loglike)
                }
                
                
                # How I read the math: log likelihood over all time points, returns the best Rt vector 
                loglikelihood2 <- function(R_vec, cases, omega, n){
                   
                    lambda_vec<-rep(0, n)
                    for (t in 1:n){
                        
                        # Calculate lambda at each time point
                        lambda<-0
                        for (s in 1:t){ # add up phi(s)*C(t-s) for each s up to t
                            lambda<-lambda + omega[s]*cases[t-s+1]
                        }
                        lambda_vec[t]<-lambda
                    }
                    
                    loglike<- sum(dpois(cases, lambda = R_vec*lambda_vec, log = T))
                    return(loglike)
                    
                }
                # Altneratively, find the loglikelihood for just a single time point
                loglikelihood3<- function(t_int, R, cases, omega, n){
                    lambda_s<-0
                    for (s in 1:t_int){ # add up phi(s)*C(t-s) for each s up to t
                        lambda_s<-lambda_s + omega[s]*cases[t_int-s+1]
                    }
                    loglike<- dpois(cases[t_int], lambda = R*lambda_s, log = T)
                    return(loglike)
                }
                    
                   
                # loop over the time series
                # for (j in 1:n){
                #     
                #     data <- list(t_int = j,
                #                  cases = cases[1:j],
                #                  omega = omega[1:j],
                #                  prior_max = 1000)
                    
                    
                    # Load fitted model -------------------------------------------------------
                    
                    #stanfit <- sampling(stanmodel, data = data, cores = 4, chains = 4, iter = 10000)   
                    
                    
                    
                    
                    # set the starting seed for the first MCMC iteration
                for (j in 1:n) {
                    if (cases[j] >0){
                        R[j,1] <- 1 # set first MCMC iteration as 1
                        loglike <- loglikelihood3(j, R[j,1], cases, omega, n)
                        prior <- log(dunif(R[j], min = 0, max = 5)) # assumes a uniform prior (0, 1000) 
                        logpcurrent <- loglike # + prior?
                    }
                    else{
                        R[j,1] <- 0 # Assumes that if cases are 0 on that day, the first iteration value is 0, but will be updated
                    }
                    
                    # loop over MCMC iterations
                    for (e in 2:nsim){
                        R[j,e] = R[j, e-1] # assumes that if cases are 0 on that day, set Rt as previous value 
                        if (cases[j] >0){ # if cases are greater than 0, try a move for Rt 
                            Rtmp <- -1
                            while (Rtmp<0){ # why not sample R on the log scale? Makes sure R isn't negative 
                                Rtmp <- R[j, e-1] + rnorm(1, 0, 0.05) # fixed step-size
                            }
                            tmploglike <- loglikelihood3(j, Rtmp, cases, omega, n)
                            tmpprior <- log(dunif(Rtmp, min = 0, max = 5))
                            logptmp <- tmploglike # + prior?
                            # standard metropolis hasting accept/reject 
                            if (runif(1) < exp(logptmp - logpcurrent)){ # accept with some probability (p_accept = p_tmp/p_current)
                                logpcurrent <- logptmp
                                R[j,e] = Rtmp
                            } # end temperature annealing
                        } # end update Rt
                    } # end loop through MCMC iterations
                }# end loop through each time point
                
                # Look at the trace plots (across nsim in R matrix, for a random time point to check for convergence) 
               # plot(1:nsim, R[50,],type = "l")
                
                # Sample from the MCMC simulations? 
                samples<-sample(100:ncol(R), n_draws) # exlude the first 100 because of burn-in
                R_sampled <- R[, samples]
                R_df <- as.data.frame(R_sampled)
                #R_df <- as.data.frame(R)
                R_df['t']<- t
                Rt_df_m <- R_df %>% pivot_longer(!t, names_to = "samples", values_to = "Rt")
                Rt_df_m <- Rt_df_m %>% 
                    group_by(samples) %>% 
                    mutate(Rt_7d = zoo::rollmean(Rt, 7, fill = NA, align = "center"),
                           sample_draw = paste0(samples, m)) %>% 
                    ungroup()
                Rt_df <- bind_rows(Rt_df, Rt_df_m)
            }
            # assign time to dates
            mean_delay = MEAN_GI + 2 # time to develop symptoms + time to seek testing (2 days)
            if (ACCOUNT_FOR_DELAY == TRUE){
                dates<-seq(as.Date(min(cases_var_t_date$date) - mean_delay), as.Date(max(cases_var_t_date$date) - mean_delay), "days")
            }
            if (ACCOUNT_FOR_DELAY == FALSE){
                dates<-seq(as.Date(min(cases_var_t_date$date)), as.Date(max(cases_var_t_date$date)), "days")
            }
            
        Rt_summary_v <- Rt_df %>% group_by(t) %>% 
            summarise( Rt_median = quantile(Rt, probs = 0.5, na.rm = T),
                       Rt_lb = quantile(Rt, probs = 0.025, na.rm = T),
                       Rt_ub = quantile(Rt, probs = 0.975, na.rm = T),
                       Rt_median_7d = quantile(Rt_7d, probs = 0.5, na.rm = T),
                       Rt_lb_7d = quantile(Rt_7d, probs = 0.025, na.rm = T),
                       Rt_ub_7d = quantile(Rt_7d, probs = 0.975, na.rm = T))
        Rt_summary_v['variant']<-unique_variants[v]
        Rt_summary_v['dates']<-dates
        
        # Put together into an Rt summary
        if (v ==1){
            Rt_summary<-Rt_summary_v
        }
        else{
            Rt_summary<-rbind(Rt_summary, Rt_summary_v)
        }
        } # loop through variants

        Rt_summary$dates<-ymd(Rt_summary$dates)
        this_variant_t<- left_join(this_variant_t, Rt_summary, by = c("collection_date" = "dates", "variant", "t"))
        
       plot<- Rt_summary %>%ggplot() + geom_line(aes(x = dates, y = Rt_median, color = variant))+
            geom_ribbon(aes(x = dates, ymin = Rt_lb, ymax = Rt_ub, fill = variant), alpha = 0.2) + 
            geom_hline(aes(yintercept = 1), linetype = "dashed") +
            theme_bw() + xlab('Collection date') + ylab('Rt estimate')
       plot
       plot<- Rt_summary %>%ggplot() + geom_line(aes(x = dates, y = Rt_median_7d, color = variant))+
           geom_ribbon(aes(x = dates, ymin = Rt_lb_7d, ymax = Rt_ub_7d, fill = variant), alpha = 0.2) + 
           geom_hline(aes(yintercept = 1), linetype = "dashed") +
           theme_bw() + xlab('Collection date') + ylab('Rt estimate')
       plot
        stopifnot('Rts didnt join' = any(this_variant_t$Rt_median >0.5 & this_variant_t$Rt_median <2))
    } # end if statement for Ana method
    

 this_variant_t %>% ggplot() +  geom_line(aes(x = collection_date, y = Rt_medians), color = "green") +
     geom_ribbon(aes(x = collection_date, ymin = Rt_lowerbounds, ymax = Rt_upperbounds),fill = 'green', alpha = 0.2) +
     geom_line(aes(x = collection_date, y = Rt_median_7d), color = "blue") +
     geom_ribbon(aes(x = collection_date, ymin = Rt_lb_7d, ymax = Rt_ub_7d),fill = 'blue', alpha = 0.2) +
     geom_hline(aes(yintercept =1), linetype = "dashed")+ facet_wrap(~variant)+
     theme_bw() + xlab('Collection date') + coord_cartesian(ylim = c(0,5))
 
    multiplier <-max(this_variant_t$raw_cases_med[this_variant_t$variant == this_variant])/7
    if(RT_INPUT == 'model observed cases' & Rt_METHOD == 'both'){
    this_variant_t %>% filter(variant == this_variant) %>% ggplot() + 
        geom_bar(aes(x = collection_date, y = raw_cases_mod_pred_med/multiplier),alpha = 0.3, stat = "identity")+
        geom_errorbar(aes(x = collection_date, ymin = raw_cases_mod_pred_lb/multiplier, ymax = raw_cases_mod_pred_ub/multiplier)) +
        geom_line(aes(x = collection_date, y = Rt_medians), color = "green") +
        geom_ribbon(aes(x = collection_date, ymin = Rt_lowerbounds, ymax = Rt_upperbounds),fill = 'green', alpha = 0.2) +
        geom_line(aes(x = collection_date, y = Rt_median_7d), color = "blue") +
        geom_ribbon(aes(x = collection_date, ymin = Rt_lb_7d, ymax = Rt_ub_7d),fill = 'blue', alpha = 0.2) +
        geom_hline(aes(yintercept =1), linetype = "dashed")+
        scale_y_continuous(name = "R(t)", sec.axis = sec_axis( trans=~.*400, name=paste0("Inferred cases with ", this_variant)))+
        theme_bw() + xlab('Collection date')
    }
    
    if(RT_INPUT == 'model mean cases' & Rt_METHOD == 'both'){
        this_variant_t %>% filter(variant == this_variant) %>% ggplot() + 
            geom_bar(aes(x = collection_date, y = raw_cases_med/multiplier), alpha = 0.3, stat = "identity")+
            geom_errorbar(aes(x = collection_date, ymin = raw_cases_lb/multiplier, ymax = raw_cases_ub/multiplier)) +
            geom_line(aes(x = collection_date, y = Rt_medians), color = "green") +
            geom_ribbon(aes(x = collection_date, ymin = Rt_lowerbounds, ymax = Rt_upperbounds),fill = 'green', alpha = 0.2) +
            geom_hline(aes(yintercept =1), linetype = "dashed")+
            geom_line(aes(x = collection_date, y = Rt_median_7d), color = "blue") +
            geom_ribbon(aes(x = collection_date, ymin = Rt_lb_7d, ymax = Rt_ub_7d),fill = 'blue', alpha = 0.2) +
            scale_y_continuous(name = "R(t)", sec.axis = sec_axis( trans=~.*400, name=paste0("Inferred cases with ", this_variant)))+
            theme_bw() + xlab('Collection date')
    }
    
    
        
    

# Lag time calculations ------------------------------------------------------

    
    # Compare lag time in submission of VOC to all other variants
    mean_lag_time_df<- this_variant_t %>% filter(variant == this_variant) %>% group_by(eweek) %>%
        mutate(eweek_date = ymd(first(collection_date))) %>% 
        ungroup() %>% 
        select( eweek_date, mean_lag_time_variant, mean_lag_time_others) %>%
        pivot_longer( 
                      cols = starts_with("mean_lag_time_"),
                      names_to = "variant", 
                      names_prefix = "mean_lag_time_",
                      values_to = "mean_lag_time" ) %>% distinct()
    se_lag_time_df<- this_variant_t %>% filter(variant == this_variant) %>% group_by(eweek) %>%
        mutate(eweek_date = ymd(first(collection_date))) %>% 
        ungroup() %>% 
        select( eweek_date, lag_time_se_variant, lag_time_se_others) %>% 
        pivot_longer( 
            cols = starts_with("lag_time_se_"),
            names_to = "variant", 
            names_prefix = "lag_time_se_",
            values_to = "lag_time_se" ) %>%distinct()
    this_lag_time_long<-left_join(mean_lag_time_df, se_lag_time_df, by = c("eweek_date" = "eweek_date", "variant" = "variant"))
    this_lag_time_long$variant[this_lag_time_long$variant == 'variant']<- this_variant
    this_lag_time_long['country']<-this_country
    
    
    this_lag_time_long %>%ggplot() + 
        geom_point(aes(x = eweek_date, y = mean_lag_time, color = variant), size = 4) + 
        geom_linerange(aes(x = eweek_date, ymin = mean_lag_time - lag_time_se, ymax = mean_lag_time +lag_time_se, color = variant)) + 
        scale_color_manual(values = c("purple", "gray")) + 
        xlab('Date') + ylab('Lag time') + theme_bw() + 
        scale_x_date(date_labels = "%m-%d") + theme(axis.text.x = element_text(angle=45, size = 10, hjust =1))+
        theme(text = element_text(size = 14))
    
    
    

# Concatenate data for visuals in report ------------------------------------
    
    growth_adv_df<-bind_rows(growth_adv_df, this_growth_adv_df)
    variant_t_all<-bind_rows(variant_t_all, this_variant_t)
    var_t_long<-bind_rows(var_t_long, this_var_t_long)
    lag_time_long<-bind_rows(lag_time_long, this_lag_time_long)
    converage_df <- bind_rows(coverage_df, this_coverage_df)
}






# Save files for load into Rmd -------------------------------------------------

# add dates to all of these
UPDATE_DATE<- max(nowcast$collection_date, na.rm = T)+1 # edit 
growth_adv_df$date_data_run <- UPDATE_DATE
variant_df$date_data_run <- UPDATE_DATE
variant_t_all$date_data_run <- UPDATE_DATE
var_t_long$date_data_run <- UPDATE_DATE
lag_time_long$date_data_run <- UPDATE_DATE
coverage_df$date_data_run <- UPDATE_DATE
if (USE_CASE == 'local'){
    write.csv(variant_df, '../data/output/variant_df.csv', row.names = F)
    write.csv(growth_adv_df, '../data/output/growth_adv_df.csv',row.names = F)
    write.csv(variant_t_all, '../data/output/variant_t_all.csv', row.names = F)
    write.csv(lag_time_long, '../data/output/lag_time_long.csv', row.names = F)
    write.csv(coverage_df, '../data/output/coverage_df.csv', row.names = F)
    # Appended datasets
}

if (USE_CASE == 'domino'){
    write.csv(variant_df, '/mnt/data/output/variant_df.csv', row.names = F)
    write.csv(growth_adv_df, '/mnt/data/output/growth_adv_df.csv',row.names = F)
    write.csv(variant_t_all, '/mnt/data/output/variant_t_all.csv',row.names = F)
    write.csv(lag_time_long, '/mnt/data/output/lag_time_long.csv',row.names = F)
}


# Save files for retrospective analysis

if (FINAL_RUN == TRUE){
    growth_adv_history <- bind_rows(growth_adv_history, growth_adv_df) # this will be most critical 
    variant_df_history <- bind_rows(variant_df_history, variant_df)
    variant_t_all_history <- bind_rows(variant_t_all_history, variant_t_all) # also will want to save this
    var_t_long_history <- bind_rows(var_t_long_history, var_t_long) # all 1000 draws from model output, will want to save this for model validation
    lag_time_long_history <- bind_rows(lag_time_long_history, lag_time_long)

    if (USE_CASE == 'local'){
        write.csv(growth_adv_history, '../data/output/growth_adv_history.csv', row.names = F)
        write.csv(variant_df_history, '../data/output/variant_df_history.csv', row.names = F)
        write.csv(variant_t_all_history, '../data/output/variant_t_all_history.csv', row.names = F)
        write.csv(var_t_long_history, '../data/output/var_t_long_history.csv', row.names = F)
        write.csv(lag_time_history, '../data/output/lag_time_long_history.csv', row.names = F)
    }
    if (USE_CASE == 'domino'){
        write.csv(growth_adv_history, '/mnt/data/output/growth_adv_history.csv', row.names = F)
        write.csv(variant_df_history, '/mnt/data/output/variant_df_history.csv', row.names = F)
        write.csv(variant_t_all_history, '/mnt/data/output/variant_t_all_history.csv', row.names = F)
        write.csv(var_t_long_history, '/mnt/data/output/var_t_long_history.csv', row.names = F)
        write.csv(lag_time_history, '/mnt/data/output/lag_time_long_history.csv', row.names = F)
    }
}
