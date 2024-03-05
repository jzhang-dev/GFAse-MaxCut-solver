#!/usr/bin/python3
# -*- coding: utf-8 -*-

INPUT_FOLDER = "data/batch_1"
OUTPUT_FOLDER = "output/batch_1"

import os
import subprocess

for config in range(8):
    contacts_path = os.path.join(INPUT_FOLDER, f"config_{config}_chr6_contacts.csv")
    variant_path = os.path.join(INPUT_FOLDER, f"config_{config}_chr6_variant_id.txt")
    output_path = os.path.join(OUTPUT_FOLDER, f"config_{config}")
    os.makedirs(output_path, exist_ok=True)

    command = ['/workspace/GFAse-MaxCut-solver/build/solve_maxcut',
    "--id_path", variant_path, "--graph_path", contacts_path, 
    "--output_dir", output_path, "--core_iterations", "200", "--sample_size", "30", 
    "--n_rounds", "5", "--threads", "4"]

    print(" ".join(command))
    subprocess.run(command)

    os.rename(os.path.join(output_path, "components_final.csv"),
              os.path.join(OUTPUT_FOLDER, f"config_{config}_components_final.csv"))











