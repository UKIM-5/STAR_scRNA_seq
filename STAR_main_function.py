#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 17:00:41 2024

@author: manueltrebo

#Note: This script allows to generate the adata objects as well as the annotations as seen in our STAR Protocols Scheiber et al 2024.

"""

import os
import argparse
from adata_setup.adata_processing import  adata_processing
from adata_setup.adata_annotation import cell_annotation

def main(sample_type):
    
    setup_dir = os.getcwd()
    raw_data_dir = os.path.join(f"{setup_dir}/raw_data/")
    results_dir = os.path.join(f"{setup_dir}/results/")
    
    #setup result directory:
    if not os.path.exists(results_dir):
      os.mkdir(results_dir)
      print("Folder %s Created!" % results_dir)
    else:
      print("Folder %s already exists" % results_dir)
    
    file_base_analysis = os.path.join(f'{setup_dir}/results/{sample_type}')
    
    adata_processing(sample_type, file_base_analysis, setup_dir)
    cell_annotation(sample_type, file_base_analysis, setup_dir)

    print("FINISHED")
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scheiber et al 2024 STAR protocols")
    parser.add_argument("--sample_type", type=str, help="Choose sample type (possibilities: lung, liver, prostate)")
    
    args = parser.parse_args()
    main(args.sample_type)
