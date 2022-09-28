# Author: Zach Susswein, modified by Kaitlyn Johnson
# Date initiated: 07-8-2022

# This script runs the MLE estimation of the multinomial model individually on 
# each of the countries specificed. This will be used for comparison to the 
# hierarchical multicountry multinomial model

# It also calculates the Brier score for each country for the calibration and 
# forecast periods of the MLE 


rm(list = ls())
USE_CASE = Sys.getenv("USE_CASE")
if(USE_CASE == ""){
    USE_CASE<-'local'
}
# Libraries ---------------------------------------------------------------

library(tidyverse)
library(nnet)
library(lubridate)
library(DescTools)

# Constants --------------------------------------------------------------
mle_countries <- read.csv('../data/processed/validation_data/countries_fit_to_MLE.csv')
COUNTRIES<-mle_countries$x
reference_data_df <- read_csv('../data/processed/validation_data/reference_data_used.csv')
variant_t <- read_csv('../data/processed/variant_t_2022-07-01.csv') 
variants_were_tracking <- unique(variant_t$variant)
N_BOOT = 100 #5e3
TRANS_ADV_MULTIPLIER <- 7 # corresponds to a weekly fitness advantage 
NUM_DAYS_FORECAST <- 21 # number of days beyond the last collection date to forecast


df_last.orig<-read_csv('../data/processed/validation_data/lineage_t_for_comp_2022-07-01.csv')

LAST_COLLECTION_DATE <- as.character(max(df_last.orig$collection_date[df_last.orig$tot_seq>0]))


# MLE_COUNTRIES <-c('Brazil','Denmark', 'India', 'South Africa', 
#                   'United Kingdom', 'United States', 'Portugal')
#COUNTRIES<-MLE_COUNTRIES # set this as a smaller country list for now 

df_last <- df_last.orig %>% 
    filter(collection_date <= ymd(LAST_COLLECTION_DATE)) %>% 
    select(lineage, country, collection_date, n_seq, t) %>% 
    filter(country %in% COUNTRIES) %>%  
    group_by(lineage, country, collection_date, t) %>% 
    summarize(n = sum(n_seq)) %>% 
    ungroup() %>% 
    droplevels() %>% 
    rename(n_last = n,
           t_last = t) %>% 
    select(lineage, country, collection_date, n_last, t_last)






# Loop through each reference date -----------------------------------------
MLE_t<-c()
MLE_regressors <-c()
country_metrics<- c()
country_lineage_metrics<-c()
for (j in 1:nrow(reference_data_df)){
    THIS_REF_DATE<- as.character(reference_data_df$reference_date[j])
    dat <- read_csv(paste0('../data/processed/validation_data/data_for_nnet_', THIS_REF_DATE, '.csv'))
    LAST_COLLECTION_DATE<- as.character(reference_data_df$last_collection_date[j])
    
    
    # Set up tibble to hold summarized results
    full_results <- tibble(lineage = character(),
                           t = integer(),
                           country = character(),
                           p_hat_MLE_mean = double(),
                           p_hat_MLE_ub = double(),
                           p_hat_MLE_lb = double(),
                           reference_date = character(),# initially set as character then convert
                           model_convergence = logical()) # initially set as character then convert
    estimate_results <- tibble(country = character(),
                          lineage = character(),
                          r_MLE_mean = double(),
                          r_MLE_lb = double(),
                          r_MLE_ub = double(),
                          trans_adv_MLE_mean = double(),
                          trans_adv_MLE_lb = double(),
                          trans_adv_MLE_ub = double(),
                          intercept_MLE_mean = double(),
                          intercept_MLE_lb = double(),
                          intercept_MLE_ub = double(),
                          reference_date = character(),
                          model_convergence = logical())
    
    # Iterate through the countries to fit the multinomial MLE model to as well as 
    # bootstrap confidence intervals
    for(COUNTRY in COUNTRIES){
        
        # fit model ---------------------------------------------------------------
        
        # Y values are the matrix of observed counts
        Y <- dat %>% 
            filter(country == COUNTRY) %>% 
            select(-lineage, -country) %>% 
            as.matrix()
        
        colnames(Y) <- NULL
        
        X <- as.integer(colnames(dat)[3:ncol(dat)])
        
        # Check that there's nothing weird with the col names
        stopifnot(typeof(X) == 'integer')
        stopifnot(min(X) == 1)
        stopifnot(max(X) == length(X))
        stopifnot(!any(is.na(X)))
        
        # Set X as a data frame
        X <- data.frame("t" = X)
        
        # nnet::multinom requires that there be at least one observed sequence on a 
        # given day and behaves poorly if a lineage is never observed. Correct for 
        # this by removing lineages and days that are never observed.
        
        # First, save the days and lineages that are being fit to
        vars_analyzed <- unique(dat$lineage)[rowSums(Y) > 0] # use this to collapse lineages from df_last 
        t_analyzed = colSums(Y) > 0
        
        # Remove days and lineages with no observations
        Y <- Y[rowSums(Y) > 0,]
        X <- as.matrix(X) 
        X <- X[colSums(Y) > 0,]
        Y <- Y[,colSums(Y) > 0]
        
        # Fit the multinomial regression to the observed data
        fit <- nnet::multinom(t(Y) ~ 1 + X, maxit = 1000)
        
        # Check that model has converged
        if (fit$convergence!=0){
            MODEL_CONVERGENCE = FALSE
        }else{
            MODEL_CONVERGENCE = TRUE
        }
        # stopifnot(fit$convergence == 0)
        
        # Set up the array to hold the predictions
        preds <- tibble(cat = integer(),
                        t = integer(),
                        rep = integer(),
                        pred_p = double())
        
        # Set up the array to hold the intercepts and slopes
        regressor_estimates <- tibble(intercept = double(),
                                      r = double(),
                                      trans_adv = double())
        
        # Set up stuff that will be re-used a bunch of times in the bootstrap fits
        base_case <-  tibble(intercept = 0, 
                             r = 0,
                             trans_adv = 0)
        
        lineage_time_grid <- expand_grid(cat = 1:nrow(Y),
                                         t = 1:(length(t_analyzed)+NUM_DAYS_FORECAST)) %>% 
            left_join(tibble(lineage = vars_analyzed) %>% 
                          mutate(cat = row_number()),
                      by = 'cat')
        
        
        # Fit the bootstrap replicates
        for (i in 1:N_BOOT){
            
            # Sample with replacement the indices for the bootstrap fit
            boot_indices <- sample(1:length(X), size = length(X), replace = T)
            
            # Pull out the X and Y values for the bootstrap fit
            X_boot <- X[boot_indices]
            Y_boot <- Y[,boot_indices]
            
            # Fit the model, silencing printed output so that it doesn't overwhelm
            # the console. Only works on *nix & will throw an error on windows
            capture.output(fit_boot <- nnet::multinom(t(Y_boot) ~ 1 + X_boot), 
                           file = '/dev/null')
            
            # Save the output of the bootstrap replicate
            mean.rep <- base_case %>% 
                rbind(tibble(intercept = summary(fit_boot)$coefficients[,1],
                             r = summary(fit_boot)$coefficients[,2],
                             trans_adv = exp(r*TRANS_ADV_MULTIPLIER)-1)) %>% 
                mutate(cat = row_number()) %>% 
                full_join(lineage_time_grid,
                          by = "cat") %>% 
                mutate(pred_linear = intercept + r * t) %>% 
                group_by(t) %>% 
                mutate(pred_p = exp(pred_linear) / sum(exp(pred_linear))) %>% 
                ungroup() %>% 
                select(cat, t, pred_p) %>% 
                mutate(rep = i)
            
            # Save the regressor estimates output of the bootstrap replicate
            regressor.rep <- base_case %>% 
                rbind(tibble(intercept = summary(fit_boot)$coefficients[,1],
                    r = summary(fit_boot)$coefficients[,2],
                    trans_adv = exp(r*TRANS_ADV_MULTIPLIER)-1)) %>% 
                mutate(cat = row_number()) %>%  
                left_join(tibble(lineage = vars_analyzed) %>% 
                              mutate(cat = row_number()),
                          by = 'cat') %>% 
                mutate(rep = i)
            
            # Append the replicate's predictions to the tibble of predictions
            preds <- rbind(preds, mean.rep)
            
            # Append the replicate's regressor restimates to the tibble of estimates
            regressor_estimates <- rbind(regressor_estimates, regressor.rep)
            
            # Print the percentage completed if a multiple of 100
            if (i %% 100 == 0){
                
                print((i / N_BOOT) * 100)
            }
            
        } # end iterate through bootstrap replicates
        
        stopifnot('Variants not properly aligned'= vars_analyzed == unique(regressor_estimates$lineage))
        
        
        
        
        # Generate the full model fit predictions
        base_fit <- base_case %>% 
            rbind(tibble(intercept = summary(fit)$coefficients[,1],
                         r = summary(fit)$coefficients[,2],
                         trans_adv = exp(r*TRANS_ADV_MULTIPLIER)-1)) %>% 
            mutate(cat = row_number()) %>% 
            full_join(lineage_time_grid,
                      by = c("cat")) %>% 
            mutate(pred_linear = intercept + r * t) %>% 
            group_by(t) %>% 
            mutate(p_hat_MLE_mean = exp(pred_linear) / sum(exp(pred_linear))) %>% 
            ungroup() %>% 
            select(cat, t, p_hat_MLE_mean)
        
        # Generate the transmission advantage estimates
        base_estimates <- base_case %>% 
            rbind(tibble(intercept = summary(fit)$coefficients[,1],
                         r = summary(fit)$coefficients[,2],
                         trans_adv = exp(r*TRANS_ADV_MULTIPLIER)-1)) %>% 
            mutate(cat = row_number()) %>% 
            left_join(tibble(lineage = vars_analyzed) %>% 
                          mutate(cat = row_number())) %>% 
            rename(intercept_MLE_mean = intercept,
                   r_MLE_mean =r,
                   trans_adv_MLE_mean = trans_adv) %>% 
            select(-lineage)
        
        
        # Need to evaluate each draw 
        t <- seq(from = min(preds$t), to = max(preds$t), by = 1)
        date <- seq(from = ymd(LAST_COLLECTION_DATE) +days(NUM_DAYS_FORECAST) -(length(t)-1), 
                    to = ymd(LAST_COLLECTION_DATE) + days(NUM_DAYS_FORECAST), by = "days")
        dates <- tibble(date = date,
                        t = as.integer(t))
        p_hats<- preds %>% left_join(tibble(lineage = vars_analyzed, 
                                            cat = 1:length(vars_analyzed)),
                                     by = 'cat') %>% 
            left_join(base_fit, by = c('cat', 't')) %>%
            select(-cat) %>% 
            left_join(dates, by = "t") 
        
        # Need to merge real data, collapsing lineages for comparison
        this_country_future_data <- df_last %>% filter(country == COUNTRY) %>% 
            mutate(
                lineage_collapsed = ifelse(lineage %in% vars_analyzed, lineage, 'other') 
            ) %>% 
            group_by(country, collection_date, lineage_collapsed) %>% 
            summarise(
                n_lineages = sum(n_last, na.rm = T)
            )
            
        lin_counts<- p_hats %>% left_join(this_country_future_data,
                                         by = c("date" = "collection_date",
                                                "lineage" = "lineage_collapsed")) %>% 
            group_by(t, rep) %>% 
            mutate(N_tot_seq_day = sum(n_lineages, na.rm = T),
                   sum_pred_p = sum(pred_p, na.rm = T)) %>% 
            ungroup()
        # Run a loop through each rep and time point
        counts<-c()
        for(m in 1:N_BOOT){
            for(k in 1:length(unique(lin_counts$t))){
                # get the subsetted data
                slice <- lin_counts %>% filter(rep == m, t == k)
                if (any(is.na(slice$pred_p))){
                  slice$pred_p<-1
                }
                  Y_tildes<- rmultinom(n=1, size = slice$N_tot_seq_day, prob = slice$pred_p)
                slice$Y_tilde<-Y_tildes[,1]
                counts<-bind_rows(counts, slice)
            }
        }
        
        
        lin_counts<- counts %>% 
            mutate(lineage_collapsed = ifelse( lineage %in% variants_were_tracking, lineage, 'other')) %>% 
            group_by(t, date, rep, country, lineage_collapsed) %>% 
            summarise(
                pred_p = sum(pred_p),
                Y_tilde = sum(Y_tilde), 
                n_lineages = sum(n_lineages)) %>% 
            ungroup() %>% 
            mutate(
                multicountry_period =case_when( 
                    date <= ymd(LAST_COLLECTION_DATE) ~ 'multicountry calibration',
                    date > ymd(LAST_COLLECTION_DATE) ~ 'multicountry forecast')) %>% 
            group_by(rep) %>% 
            mutate(N_tot_seq = sum(n_lineages, na.rm = T)) %>% 
            ungroup() %>% 
            group_by(multicountry_period, rep) %>% 
            mutate(N_tot_seq_period = sum(n_lineages, na.rm = T)) %>% 
            ungroup() %>% 
            group_by(t, rep) %>% 
            mutate(N_tot_seq_day = sum(n_lineages, na.rm =T)) %>% 
            ungroup() %>% 
            rename(lineage = lineage_collapsed)
        lin_counts$n_lineages[is.na(lin_counts$n_lineages)]<-0
        
        lin_counts<- lin_counts %>% left_join(base_fit %>% left_join(tibble(lineage = vars_analyzed, 
                                                                        cat = 1:length(vars_analyzed)),
                                                                 by = 'cat'),
                                              by = c('lineage', 't'))
        
        # CCC metrics-------------------------------------------------
        
        
        # CCC by country 
        CCC_all <- lin_counts %>% 
            group_by(rep) %>% 
            summarise(
                CCC_by_country = CCC(n_lineages, Y_tilde)$rho.c$est
            ) %>% 
            ungroup() %>% 
            summarise(
                CCC_country_med_MLE = quantile(CCC_by_country, probs = 0.5, na.rm = T),
                CCC_country_mean_MLE = mean(CCC_by_country, na.rm = T),
                CCC_country_lb_MLE = quantile(CCC_by_country, probs = 0.025, na.rm = T),
                CCC_country_ub_MLE = quantile(CCC_by_country, probs = 0.975, na.rm = T)
            ) %>% 
            mutate(country = COUNTRY)
        
        # CCC by country and variant 
        CCC_country_lin <- lin_counts %>% 
            group_by(rep, lineage) %>% 
            summarise(
                CCC_by_country_lineage = CCC(n_lineages, Y_tilde)$rho.c$est
            ) %>% 
            group_by(lineage) %>% 
            summarise(
                CCC_country_lin_med_MLE = quantile(CCC_by_country_lineage, probs = 0.5, na.rm = T),
                CCC_country_lin_mean_MLE = mean(CCC_by_country_lineage, na.rm = T),
                CCC_country_lin_lb_MLE = quantile(CCC_by_country_lineage, probs = 0.025, na.rm = T),
                CCC_country_lin_ub_MLE = quantile(CCC_by_country_lineage, probs = 0.975, na.rm = T)
            ) %>% 
            mutate(country = COUNTRY)
        
        # CCC by country, period, and variant
        CCC_lin_period <- lin_counts %>% 
            group_by(rep, lineage, multicountry_period) %>% 
            summarise(
                CCC_by_country_lineage_period = CCC(n_lineages, Y_tilde)$rho.c$est
            ) %>% 
            group_by(lineage, multicountry_period) %>% 
            summarise(
                CCC_country_lin_period_med_MLE = quantile(CCC_by_country_lineage_period, probs = 0.5, na.rm = T),
                CCC_country_lin_period_mean_MLE = mean(CCC_by_country_lineage_period, na.rm = T),
                CCC_country_lin_period_lb_MLE = quantile(CCC_by_country_lineage_period, probs = 0.025, na.rm = T),
                CCC_country_lin_period_ub_MLE = quantile(CCC_by_country_lineage_period, probs = 0.975, na.rm =T)
            )%>% 
            mutate(country = COUNTRY)
       
         # CCC by country and period 
        CCC_period <- lin_counts %>% 
            group_by(rep, multicountry_period) %>% 
            summarise(
                CCC_by_country_period = CCC(n_lineages, Y_tilde)$rho.c$est
            ) %>% 
            group_by(multicountry_period) %>% 
            summarise(
                CCC_country_period_med_MLE = quantile(CCC_by_country_period, probs = 0.5, na.rm = T),
                CCC_country_period_mean_MLE = mean(CCC_by_country_period, na.rm = T),
                CCC_country_period_lb_MLE = quantile(CCC_by_country_period, probs = 0.025, na.rm = T),
                CCC_country_period_ub_MLE = quantile(CCC_by_country_period, probs = 0.975, na.rm = T)
            ) %>% 
            mutate(country = COUNTRY)
        
        #Combine the country variant level CCC metrics
        CCC_lineages <- CCC_lin_period %>% 
            left_join(CCC_country_lin, by = c('country', 'lineage')) %>% 
            mutate(reference_date = THIS_REF_DATE,
                    country = COUNTRY) %>% 
            relocate(country, .before = everything())
        country_lineage_metrics <- bind_rows(country_lineage_metrics, CCC_lineages)
    
        
            
        # Brier score overall
        BS_all <- lin_counts %>% 
            group_by(rep, lineage) %>% 
            summarise(
                sum_over_t = sum( n_lineages - 2*pred_p*n_lineages + N_tot_seq_day*(pred_p^2)), # using sufficient stat for for (p_hat - Y)^2
                N_tot_seq = max(N_tot_seq),
                sum_over_t_mean = sum( n_lineages - 2*p_hat_MLE_mean*n_lineages + N_tot_seq_day*(p_hat_MLE_mean^2))
            ) %>% ungroup() %>% 
            group_by(rep) %>% 
            summarise( Brier_score = 1/max(N_tot_seq)*sum(sum_over_t),
                       Brier_score_mean = 1/max(N_tot_seq)*sum(sum_over_t_mean)) %>% 
            summarise(BS_med_MLE = quantile(Brier_score, probs = 0.5, na.rm = T),
                      BS_lb_MLE = quantile(Brier_score, probs = 0.025, na.rm = T),
                      BS_ub_MLE = quantile(Brier_score, probs = 0.975, na.rm = T),
                      BS_mean_MLE = max(Brier_score_mean))
        
        
        # Brier score by period 
        BS_by_period <- lin_counts %>% 
            group_by(rep, lineage, multicountry_period) %>% 
            summarise(
                sum_over_t = sum( n_lineages - 2*pred_p*n_lineages + N_tot_seq_day*(pred_p^2), na.rm = T),
                N_tot_seq_period = max(N_tot_seq_period),
                sum_over_t_mean = sum( n_lineages - 2*p_hat_MLE_mean*n_lineages + N_tot_seq_day*(p_hat_MLE_mean^2))
            ) %>% ungroup() %>% 
            group_by(rep, multicountry_period) %>% 
            summarise(Brier_score_period = 1/max(N_tot_seq_period)*sum(sum_over_t, na.rm = T),
                      Brier_score_mean = 1/max(N_tot_seq_period)*sum(sum_over_t_mean)) %>% 
            group_by(multicountry_period) %>% 
            summarise(BS_med_period_MLE = quantile(Brier_score_period, probs = 0.5, na.rm = T),
                      BS_lb_period_MLE = quantile(Brier_score_period, probs = 0.025, na.rm = T),
                      BS_ub_period_MLE = quantile(Brier_score_period, probs = 0.975, na.rm = T),
                      BS_mean_period_MLE = max(Brier_score_mean))
        
        # combine BS and CCC 
        this_country_metrics<- BS_by_period %>% cbind(BS_all) %>% 
            mutate(country = COUNTRY,
                   reference_date = THIS_REF_DATE) %>%
            left_join(CCC_period, by = c('country', 'multicountry_period')) %>% 
            left_join(CCC_all, by = 'country') %>% 
            mutate(country = COUNTRY) %>% 
            relocate(country, .before = everything())
       
         # add this country to the collection of BSs
        country_metrics <- bind_rows(country_metrics, this_country_metrics)
        
        # Save the bootstrap CIs appending the mean of the full data fit
        boot_se <- preds %>% 
            group_by(cat, t) %>% 
            summarize(p_hat_MLE_lb = quantile(pred_p, 0.025, na.rm = T),
                      p_hat_MLE_ub = quantile(pred_p, 0.975, na.rm = T),
                      .groups = 'drop') %>% 
            left_join(tibble(lineage = vars_analyzed,
                             cat = 1:length(vars_analyzed)),
                      by = 'cat') %>% 
            left_join(base_fit,
                      by = c('cat', 't')) %>% 
            mutate(country = COUNTRY) %>% 
            # add collection_date to the saved output 
            select(-cat) %>% left_join(dat %>% filter(country == COUNTRY)
                                       %>% select(lineage)) %>% 
            mutate(reference_date = THIS_REF_DATE,
                   model_convergence = MODEL_CONVERGENCE)# maps to all ts!
        
        full_results <- rbind(full_results, boot_se)
        
        
        boot_regressor_se <- regressor_estimates %>% 
            group_by(cat) %>% 
            summarise(r_MLE_lb = quantile(r, 0.025, na.rm = T),
                      r_MLE_ub = quantile(r, 0.975, na.rm = T),
                      intercept_MLE_lb = quantile(intercept, 0.025, na.rm = T),
                      intercept_MLE_ub = quantile(intercept, 0.975, na.rm = T),
                      trans_adv_MLE_lb = quantile(trans_adv, 0.025, na.rm = T),
                      trans_adv_MLE_ub = quantile(trans_adv, 0.975, na.rm = T),
                      .groups = 'drop') %>% 
            left_join(tibble(lineage = vars_analyzed,
                             cat = 1:length(vars_analyzed)),
                      by = 'cat') %>% 
            left_join(base_estimates, by = 'cat') %>% 
            mutate(country = COUNTRY) %>% 
            select(country, lineage, r_MLE_mean, r_MLE_lb, r_MLE_ub, 
                   trans_adv_MLE_mean, trans_adv_MLE_lb, trans_adv_MLE_ub, 
                   intercept_MLE_mean, intercept_MLE_lb, intercept_MLE_ub) %>% 
            mutate(reference_date = THIS_REF_DATE,
                   model_convergence = MODEL_CONVERGENCE)
        
        estimate_results <- bind_rows(estimate_results, boot_regressor_se)
            
    } # end loop over countries
    
    # Combine for all dates the prevalence estimates over time
    full_results$reference_date<- THIS_REF_DATE
    MLE_t <- bind_rows(MLE_t, full_results)
    
    # Combine for all reference dates the fitness advantage estimates over time
    estimate_results$reference_date <- THIS_REF_DATE
    MLE_regressors <- bind_rows(MLE_regressors, estimate_results)
    
    
    
}# end loop over all reference dates

country_metrics <- country_metrics %>% 
  relocate(reference_date, .before = multicountry_period)
country_lineage_metrics <- country_lineage_metrics %>% 
  relocate(reference_date, .before = multicountry_period)



MLE_t$reference_date <- ymd(MLE_t$reference_date)
MLE_regressors$reference_date <- ymd(MLE_regressors$reference_date)
country_metrics$reference_date <- ymd(country_metrics$reference_date)
country_lineage_metrics$reference_date <- ymd(country_lineage_metrics$reference_date)
   
# Quick fix to append last reference date----------------------------------------
# MLE_regressors_first<- read_csv('../data/processed/validation_data/MLE_regressors.csv')
# MLE_t_first <- read_csv('../data/processed/validation_data/MLE_t.csv')
# country_metrics_MLE<- read_csv('../data/processed/validation_data/country_metrics_MLE.csv')
# country_lineage_metrics_MLE <- read_csv('../data/processed/validation_data/country_lineage_metrics_MLE.csv')
# 
# MLE_t<- bind_rows(MLE_t_first, MLE_t)
# MLE_regressors <- bind_rows(MLE_regressors_first, MLE_regressors)
# country_metrics <- bind_rows(country_metrics_MLE, country_metrics)
# country_lineage_metrics <- bind_rows(country_lineage_metrics_MLE, country_lineage_metrics)



# Save the concatenated output for merging into the evaluation-----------------


# p_hat with upper and lower bounds for each lineage-country-time point
write.csv(MLE_t, '../data/processed/validation_data/MLE_t.csv', row.names = F)
# fitness advantage estimates with upper and lower bounds for each lineage-country
write.csv(MLE_regressors, '../data/processed/validation_data/MLE_regressors.csv', row.names = F)
write.csv(country_metrics, '../data/processed/validation_data/country_metrics_MLE.csv', row.names = F)
write.csv(country_lineage_metrics, '../data/processed/validation_data/country_lineage_metrics_MLE.csv', row.names = F)