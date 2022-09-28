The purpose of this folder is to mimic the code folder (which will serve
as the main pipeline for the variant tracker reports) except for multiple
"reference" datasets, which refer to timestamped downloads from GISAID.

The pipeline works as follows:

pre_processing_HV.R 
inputs: reference_data_df, timestamped metadata
outputs: reference_data_used, timestamped lineage_t

prepare_data_for_cmdstan_HV.R
inputs: reference_data_used, timestamped lineage_t, stan model
output: timestamped data_for_cmdstan, countries to fit MLE,
    timetsamped data_for_nnet

run_model_using_shell_scripts.R
inputs: reference_data_used
dependencies: run_cmstan_HV.sh, process_cmdstan_HV.sh
outputs: time stamped stan objects (4 chains for each), 
    time stamped p_hat summaries
    
process_results_from_stansummary_HV.R
inputs: reference_data_used, time stamped lineage_t, time stamped p_hat summaries
outputs: time stamped clean_global_df

process_results_from_cmdstan_HV.R
inputs: reference_data_used, time stamped lineage_t,
    time stamped stan objects
outputs: time stamped mu_hat, time stamped r_summary

run_MLE_fit_HV.R
inputs: reference_data_used, countries to fit MLE, timestamped data for nnet
outputs: MLE_t, MLE_regressors (p_hat and r estimates for MLE multinomial for
    each lineage and country in countries fit to MLE)

evaluate_output_from_HV.R
inputs: reference_data_used, time stamped r_summary, time stamped mu_hat,
time stamped clean_global_df, MLE_t, MLE_regressors
outputs: figures and metrics for comparison