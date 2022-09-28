#!/usr/bin/env bash
# This file runs the AL1 model, checks it for convergence, and saves the 
# output as a csv

# Set paths
output_path="data/output/multicountry_output"
cmdstan_path="$HOME/cmdstan/bin"

# Run sampler
for i in {1..4}
do
      ./code/stancode/multivariate_variant_multinomial_ncp \
       method=sample adapt delta=.8 algorithm=hmc engine=nuts \
       max_depth=12 num_warmup=2500 num_samples=500 init=0  data \
             file=data/processed/AL1_data_for_cmdstan.json \
      output file="$output_path"/AL1_output_${i}.csv &
done

wait

echo "Model run completed and output saved"

# Check model diagnostics
"$cmdstan_path"/diagnose "$output_path"/AL1_output_*.csv 

wait

# Save output and suppress printing to stdout
if [ -f "data/processed/AL1_summary.csv" ]; then
    rm "data/processed/AL1_summary.csv"
fi

"$cmdstan_path"/stansummary -c "data/processed/AL1_summary.csv" -s 8 "$output_path"/AL1_output_*.csv >/dev/null
echo "Output processed. Saved in data/processed/AL1_summary.csv"

exit 0
