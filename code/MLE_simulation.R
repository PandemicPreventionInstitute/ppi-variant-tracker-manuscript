# Author: Kaitlyn Johnson
# Date initiated: 2022-07-29

# The purpose of this script is to use simulation to identify how sequencing 
# quantity affects the accuracy of the MLE estimates of the intercept and growth 
# advantage.

#rm(list = ls())
#USE_CASE = Sys.getenv("USE_CASE")
#if(USE_CASE == ""){
#  USE_CASE<-'local'
#}
# Libraries ---------------------------------------------------------------

library(tidyverse)
library(nnet)
library(lubridate)
library(DescTools)


THIS_COUNTRY <- 'Portugal'
DURATION <- 90
TRANS_ADV_MULTIPLIER <- 7
if (USE_CASE == 'databricks'){
  DATA_PATH <- 'dbfs/FileStore/tables/ppi-variant-tracker/data/processed/validation_data/lineage_t_for_comp_2022-07-01.csv'
}else{
    DATA_PATH <- paste0('../data/processed/validation_data/lineage_t_for_comp_2022-07-01.csv')
}

df.orig <- read_csv(DATA_PATH) %>% filter(country == THIS_COUNTRY)

# Going to model starting this many days in the past
REF_DATE <- "2022-07-01"
FIRST_DATE <- ymd(REF_DATE) - days(DURATION)
df.orig <- df.orig %>% filter(collection_date >= FIRST_DATE)




# Wrangle data ------------------------------------------------------------

df <- df.orig %>% 
  select(lineage, country, collection_date, n_seq, t) %>% 
  group_by(lineage, collection_date, t) %>% 
  summarize(n = sum(n_seq)) %>% 
  ungroup() %>% 
  group_by(collection_date, t) %>% 
  mutate(N= sum(n)) %>% 
  ungroup() %>% 
  droplevels() %>% 
  mutate(t = t - min(t) + 1)


df %>% ggplot() + geom_line(aes(x = collection_date, y = N))
most_prevalent_lineage <- 'BA.2'

N_t_baseline <-df %>% 
    select(t, N) %>% 
    distinct()
avg_seq_per_day = sum(N_t_baseline$N)/length(N_t_baseline$t)

# UK about 1000, Portugal 62, South Africa 50, Bangladesh 0.5
seq_per_day_vec <-c(1000, 500, 100, 50, 25, 10, 5, 1, 0.5)

# get scale factor vec to be ued to modify the number of samples sequenced,
# keeping the same dynamics of the sequences as observed in the chosen country. 
scale_factor_vec <- seq_per_day_vec/avg_seq_per_day







# Functions to do the fitting and generation of synthetic data--------------------


fit_data<- function(df, most_prevalent_lineage){
    n_lins_global <- length(unique(df$lineage))
    
    lineage_map <- tibble(lineage = df %>% 
                              filter(lineage != most_prevalent_lineage) %>%
                              pull(lineage) %>%
                              unique(),
                          lin_id = 2:n_lins_global) %>% 
        bind_rows(tibble(lineage = most_prevalent_lineage,
                         lin_id = 1))
    
    # Set up the data for nnet------------------------------------
    data_nnet <- df %>% 
        full_join(lineage_map) %>%
        arrange(lin_id) %>% 
        mutate(n = as.integer(n)) %>% 
        group_by(lineage) %>% 
        arrange(collection_date) %>% 
        ungroup() %>% 
        select(-collection_date) %>% 
        pivot_wider(id_cols = c(lineage),
                    names_from = t,
                    values_from = n) %>% 
        left_join(lineage_map) %>% 
        relocate(lin_id, .before = everything()) %>% 
        select(-lin_id)
    
    Y <- data_nnet %>% 
        select(-lineage) %>% 
        as.matrix()
    
    colnames(Y) <- NULL
    
    X <- as.integer(colnames(data_nnet)[2:ncol(data_nnet)])
    
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
    vars_analyzed <- unique(data_nnet$lineage)[rowSums(Y) > 0] # use this to collapse lineages from df_last 
    t_analyzed = colSums(Y) > 0
    
    if(length(vars_analyzed) ==1){
        base_estimates <- data.frame( intercept_MLE_mean = NA, r_MLE_mean = NA, lineage = NA, MODEL_CONVERGENCE = FALSE)
        base_fit <- data.frame (t = 1:DURATION, p_hat_MLE_mean = NA, lineage = NA, MODEL_CONVERGENCE = FALSE)
        output_list <- list(base_estimates, base_fit)
    }else{
    
    # Remove days and lineages with no observations
    Y <- Y[rowSums(Y) > 0,]
    X <- as.matrix(X) 
    X <- X[colSums(Y) > 0,]
    Y <- Y[,colSums(Y) > 0]
    
    # Fit the multinomial regression to the observed data
    fit <- nnet::multinom(t(Y) ~ 1 + X, maxit = 1000)
    
    if (fit$convergence!=0){
        MODEL_CONVERGENCE = FALSE
    }else{
        MODEL_CONVERGENCE = TRUE
    }
    
    lineage_time_grid <- expand_grid(cat = 1:nrow(Y),
                                     t = 1:length(t_analyzed)) %>% 
        left_join(tibble(lineage = vars_analyzed) %>% 
                      mutate(cat = row_number()))
    
    
    base_case <-  tibble(intercept = 0, 
                         r = 0)
    
    
    
    base_fit <- base_case %>% 
        rbind(tibble(intercept = summary(fit)$coefficients[,1],
                     r = summary(fit)$coefficients[,2])) %>% 
        mutate(cat = row_number()) %>% 
        full_join(lineage_time_grid,
                  by = c("cat")) %>% 
        mutate(pred_linear = intercept + r * t) %>% 
        group_by(t) %>% 
        mutate(p_hat_MLE_mean = exp(pred_linear) / sum(exp(pred_linear))) %>% 
        ungroup() %>% 
        select(cat, t, p_hat_MLE_mean) %>% 
        left_join(tibble(lineage = vars_analyzed, 
                         cat = 1:length(vars_analyzed)),
                  by = 'cat') %>% 
        select(-cat) %>% 
        mutate(MODEL_CONVERGENCE = MODEL_CONVERGENCE)
    stopifnot('Variants not properly aligned'= vars_analyzed == unique(base_fit$lineage))
    
    
    base_estimates <- base_case %>% 
        rbind(tibble(intercept = summary(fit)$coefficients[,1],
                     r = summary(fit)$coefficients[,2])) %>% 
        mutate(cat = row_number()) %>% 
        left_join(tibble(lineage = vars_analyzed) %>% 
                      mutate(cat = row_number())) %>% 
        rename(intercept_MLE_mean = intercept,
               r_MLE_mean =r) %>% 
        select(-cat) %>% 
        mutate(MODEL_CONVERGENCE = MODEL_CONVERGENCE)
    
    output_list <- list(base_estimates, base_fit)
    }
    return(output_list)
}




# Generate simulated data from those sets of p_hats ---------------------------

generate_sim_df<- function(base_fit, N_t){
    sim_df <- c()
    df_sim<- c()
    for (k in 1:DURATION){
        N_seq_day <- N_t$N[k]
        Y_tildes <- rmultinom(n=1, size = N_seq_day, prob = base_fit$p_true[base_fit$t==k])
        n<- Y_tildes
        df_sim_k <- data.frame(n)
        df_sim_k <- df_sim_k %>% 
            mutate(cat = row_number()) %>% 
            left_join(tibble(lineage = unique(base_fit$lineage)) %>% 
                          mutate(cat = row_number()), 
                      by = 'cat') %>% 
            mutate(
                t = k,
                N = N_seq_day, 
                collection_date = ymd(FIRST_DATE) + days(k-1)) %>% 
            select(-cat) %>% 
            select(lineage, collection_date, t, n, N) 
        df_sim <- bind_rows(df_sim, df_sim_k)
    }
    
    df_sim <- df_sim %>% arrange(lineage)
    
    return(df_sim)
}


# Get true lineage proportions over time from the parameter values
get_base_fit_truth <- function(base_estimates){
    
    lineage_time_grid <- expand_grid(cat = 1:nrow(base_estimates),
                                     t = 1:DURATION) %>% 
        left_join(tibble(lineage =base_estimates$lineage) %>% 
                      mutate(cat = row_number()))
    
    base_fit <- tibble(intercept = base_estimates$intercept,
                             r = base_estimates$r, 
                             trans_adv = exp(r*TRANS_ADV_MULTIPLIER)-1) %>% 
        mutate(cat = row_number()) %>% 
        full_join(lineage_time_grid,
                  by = c("cat")) %>% 
        mutate(pred_linear = intercept + r * t) %>% 
        group_by(t) %>% 
        mutate(p_true = exp(pred_linear) / sum(exp(pred_linear))) %>% 
        ungroup() %>% 
        select(cat, t, p_true) %>% 
        left_join(tibble(lineage = base_estimates$lineage, 
                         cat = 1:length(base_estimates$lineage)),
                  by = 'cat') %>% 
        select(-cat)
    return(base_fit)
}



# Test things out-----------------------------
fit_sim_list <- fit_data( df, most_prevalent_lineage)
# Parameter estimates
pars_fit_test <- fit_sim_list[[1]]
# Estimated proportions
fit_p_t_test<- fit_sim_list[[2]]







# Write out code for doing this 
lineage_list <- unique(df$lineage[df$n>0])
SAMPLE_MULTIPLIERS <- scale_factor_vec
N_PARAMETER_DRAWS <-100
pars_comp_df <-c()
p_t_comp<-c()
for (k in 1:length(SAMPLE_MULTIPLIERS)){
    # Sets the number of sequences per day
    N_t<- round(N_t_baseline * SAMPLE_MULTIPLIERS[k])
    avg_seq<- SAMPLE_MULTIPLIERS[k] * avg_seq_per_day
    for (i in 1: N_PARAMETER_DRAWS){
        # Generate the parameter sets for the X many lineages in the multinomial
        r <- rnorm( n= length(lineage_list)-1, mean = -1, sd =0.5) # same prior as in multinomial model
        # Generates a student t of student_t( nu = 3, mu = -5, sigma = 5) https://mc-stan.org/docs/stan-users-guide/reparameterization.html
        tau <- rgamma(n = length(lineage_list) -1, shape = 3/2, rate = 3/2)
        alpha <- rnorm(n=length(lineage_list)-1, mean = -5, sd = 5)
        intercept <- alpha/sqrt(tau)
        #intercept <- rt(n=length(lineage_list)-1, df = length(lineage_list)-1, ncp = -5) # same prior as in multinomial model ??? Need to fix this 
        lineage<- lineage_list[lineage_list != most_prevalent_lineage]
        # Parameters
        other_pars_df <- data.frame(lineage, r, intercept)
        pars_df<- tibble(lineage = most_prevalent_lineage,
                          r = 0, 
                          intercept =0) %>% bind_rows(other_pars_df)
        # True proportions from parameters
        p_t<- get_base_fit_truth(pars_df)
        # Simulated proportions/sequences from the true proportions and number of sequences per day
        df_sim <- generate_sim_df(p_t, N_t)
        # Fit the simulated data with the MLE
        fit_sim_list <- fit_data( df_sim, most_prevalent_lineage)
        # Parameter estimates
        pars_fit_df <- fit_sim_list[[1]]
        # Estimated proportions
        fit_p_t<- fit_sim_list[[2]]
        
        # Combine the true and the fit parameters
        pars_comp_df_i <- left_join(pars_df, pars_fit_df, by = 'lineage') 
        pars_comp_df_i <- pars_comp_df_i %>% 
            mutate( draw = i,
                    avg_seq = avg_seq) %>% 
            relocate(r_MLE_mean, .before = intercept)
        
        # Combine the true and the fit variant proportions
        p_t_comp_i <- left_join(p_t, fit_p_t) 
        p_t_comp_i <- p_t_comp_i %>%
            mutate(draw = i,
                   avg_seq = avg_seq) %>% 
            relocate(p_hat_MLE_mean, .before = lineage)
        
        
        p_t_comp_i %>% ggplot() + geom_line(aes(x = t, y = p_true, color = lineage)) +
            geom_point(aes(x = t, y = p_hat_MLE_mean, color = lineage))
        
        # Aggregate with other draws and other sample frequencies
        pars_comp_df <- bind_rows(pars_comp_df, pars_comp_df_i)
        p_t_comp <- bind_rows(p_t_comp, p_t_comp_i)
    }
}


# Get residuals and some summary stats

pars_comp_df <- pars_comp_df %>% 
    mutate( r_residual = r_MLE_mean - r,
            intercept_residual = intercept_MLE_mean - intercept)

pct_converge <- pars_comp_df %>% group_by(avg_seq, draw) %>% 
    summarise(MODEL_CONVERGENCE_AT_ALL = any(MODEL_CONVERGENCE==TRUE))%>% 
    ungroup() %>% 
    group_by(avg_seq) %>% 
    summarise(
    pct_converged = 100*sum(MODEL_CONVERGENCE_AT_ALL==TRUE, na.rm = T)/n())

pct_converge %>% ggplot() + geom_line(aes(x = avg_seq, y = pct_converged)) + scale_x_log10() + theme_bw()+
    xlab('Average number of sequences per day') + ylab('Percent of model runs that converged')
    
pars_comp_df %>% ggplot() + geom_point(aes(x = avg_seq, y = r_residual)) + scale_x_log10()
pars_comp_df %>% ggplot(aes(x = (r_residual), y = avg_seq, group = avg_seq)) + geom_density_ridges() + scale_y_log10()+
    xlim(c(-0.01,0.01)) + ylab('Avg number of sequences per day') + xlab('Estimated r - true r') + theme_bw()
pars_comp_df %>% ggplot() + geom_point(aes(x = avg_seq, y = intercept_residual))+ scale_x_log10()
pars_comp_df %>% ggplot(aes(x = (intercept_residual), y = avg_seq, group = avg_seq)) + geom_density_ridges() + scale_y_log10()+
    xlim(c(-0.1,0.1)) + ylab('Avg number of sequences per day') + xlab('Estimated intercept - true intercept') + theme_bw()

pars_comp_df %>% ggplot() + geom_point(aes(x = r, y = avg_seq, fill = r_residual))+
    scale_y_log10() +
    scale_fill_gradient(low = 'blue', high = 'red', na.value = 'gray',
                        guide = "colourbar", aesthetics = "fill",
                        breaks = c(-0.5,0, 1, 2), labels= c(-0.5,0,1, 2), limits = c(-0.6,2))

pars_comp_df %>% ggplot() + geom_point(aes(x = r, y = avg_seq, fill = intercept_residual))+scale_y_log10() +
    scale_fill_gradient(low = 'blue', high = 'red', na.value = 'gray',
                        guide = "colourbar", aesthetics = "fill",
                        breaks = c(-0.5,0, 1, 2), labels= c(-0.5,0,1, 2), limits = c(-0.6,2))






if(USE_CASE == 'databricks'){
    PARS_PATH <- 'dbfs/FileStore/tables/ppi-variant-tracker/data/processed/MLE_sim_parameter_comparison.csv'
    P_T_PATH <- 'dbfs/FileStore/tables/ppi-variant-tracker/data/processed/MLE_sim_proportion_t_comparison.csv'
}else{
    PARS_PATH <- '../data/processed/MLE_sim_parameter_comparison.csv'
    P_T_PATH <- '../data/processed/MLE_sim_proportion_t_comparison.csv'
}


write.csv(pars_comp_df, PARS_PATH, row.names = F)
write.csv(p_t_comp, P_T_PATH, row.names = F)
