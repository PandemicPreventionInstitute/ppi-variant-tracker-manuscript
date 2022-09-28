#!/usr/bin/env bash
# This script assumes that the cmdstan directory is located in the home 
# directory

echo "Processing model output from from: $1"

# Where to store the output
prefix='data/output/multicountry_output/validation'

# Process the output
# Use the cmdstan directory (assumed to be in $HOME)
~/cmdstan/bin/stansummary -p 2.5,25,50,75,97.5 -c "$prefix"/processed_output_$1.csv "$prefix"/output_$1_1.csv "$prefix"/output_$1_2.csv "$prefix"/output_$1_3.csv "$prefix"/output_$1_4.csv &
wait
echo "Model run processed and saved"

exit 0

