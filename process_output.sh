#!/usr/bin/env bash
# This script assumes that the cmdstan directory is located in the home 
# directory


echo "Processing output from current model"


# Where to store the output
prefix='data/output/multicountry_output'

# Process the output
# Use the cmdstan directory (assumed to be in $HOME)

~/cmdstan/bin/stansummary -c "$prefix"/processed_output.csv "$prefix"/output_2.csv "$prefix"/output_3.csv "$prefix"/output_4.csv &

wait
echo "Model run processed and saved"

exit 0

