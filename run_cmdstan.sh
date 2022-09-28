#!/usr/bin/env bash
# This script assumes that the cmdstan directory is located in the home 
# directory

# Where to store the output
prefix='data/output/multicountry_output'

# Run the model
for i in {1..4}
do
      ./code/stancode/multivariate_variant_multinomial_ncp \
       method=sample adapt delta=.8 algorithm=hmc engine=nuts \
       max_depth=12 num_warmup=2500 num_samples=500 init=0  data \
             file=data/processed/data_for_cmdstan.json \
      output file="$prefix"/output_${i}.csv &
done
wait
echo "Model run completed and output saved"


# Process the output
# Use the cmdstan directory (assumed to be in $HOME)
~/cmdstan/bin/stansummary -c "$prefix"/processed_output.csv "$prefix"/output_"${1..4}".csv &
wait
echo "Model run processed and saved"


exit 0

