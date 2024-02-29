#!/bin/bash

rm -rf output*
mkdir -p output
../build/solve_maxcut \
  --id_path data/chr6_variant_id.txt \
  --graph_path data/chr6_contacts.csv \
  --output_dir output/ \
  --core_iterations 200 \
  --sample_size 30 \
  --n_rounds 5 \
  --threads 4

