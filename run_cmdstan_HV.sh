#!/usr/bin/env bash

echo "Running model for data from: $1"

for i in {1..4}
do
      ./code/stancode/multivariate_variant_multinomial_ncp_HV \
       method=sample adapt delta=.8 algorithm=hmc engine=nuts \
       max_depth=12 num_warmup=2500 num_samples=500 init=0  data \
             file=data/processed/validation_data/data_for_cmdstan_$1.json \
      output file=data/output/multicountry_output/validation/output_$1_${i}.csv &
done

wait

echo "Model run completed and output saved"