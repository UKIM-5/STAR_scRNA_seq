#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 08:52:35 2024

@author: manueltrebo
"""

import argparse
import os
import scanpy as sc
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import leidenalg as ld

def adata_processing(sample_type: str, file_base_analysis: str, setup_dir: str):
    
    raw_data_dir = os.path.join(f"{setup_dir}/raw_data/")

    sc.settings.figdir = os.path.join(file_base_analysis, "figures")
    sc.settings.set_figure_params(figsize=(6,6), frameon=False, dpi = 1200)
    
    #setup major directories:
    if not os.path.exists(file_base_analysis):
      os.mkdir(file_base_analysis)
      print("Folder %s Created!" % file_base_analysis)
    else:
      print("Folder %s already exists" % file_base_analysis)
      
    os.chdir(file_base_analysis) #set wd
    
    dir_structure = ["tables","figures","saves"]
    
    #create directories if not existing
    for dirName in dir_structure:
        dir_loc = file_base_analysis + "/" + dirName
        
        if not os.path.exists(dir_loc):
            os.makedirs(dir_loc)
            print("Directory " , dir_loc ,  " Created!")
        else:    
            print("Directory " , dir_loc ,  " already exists!")
    
    #load data
    file_path = f'{raw_data_dir}{sample_type}_raw.h5ad'

    adata = sc.read_h5ad(file_path)
    
    print(adata)
    
    #filter needed in case sample was multiplexed - keep only necessary columns and rows for this analysis
    try:
        adata.obs["Sample_Name"]
        adata.obs = adata.obs.iloc[:,[1,2]]
        print("Filter out so only GG TU remains")
        print(f"    Before: {adata.shape[0]}")
        adata = adata[(adata.obs["Sample_Name"] == "GG")|(adata.obs["Sample_Name"] == "TU")].copy()
        print(f"    After: {adata.shape[0]}")
        adata.obs.rename(columns={"Sample_Name":'sample_origin'}, inplace=True)
        adata.obs['sample_origin'] = adata.obs['sample_origin'].replace({'GG': 'benign', 'TU': 'tumor'})
    except:
        print("This sample does not have Sample-Tag information added - adding sample column now.")
        if sample_type == "liver":
            adata.obs["sample_origin"] = "benign"

    #calculate QC_metrics:
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], log1p=False, inplace=True, percent_top=None
        )
    
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        groupby="sample_origin",
        size = 0.25,
        jitter=0.2,
        multi_panel=True,
        save=f"_plot_{sample_type}_preqc.png",
        show=False,
        xlabel = sample_type
        )
    
    #Quality Control via filtering - values are from each respective publication (see README)
    print(f"    Before: {adata.shape[1]}")
    sc.pp.filter_genes(adata, min_counts=3)
    print(f"    After: {adata.shape[1]}")
    
    filter_values = {
    "lung": {"max_total_cts": 60000, "min_total_cts": 1000, "max_genes": 8000, "min_genes": 200, "pct_counts_mt": 30, "var_genes": 6000},
    "liver": {"max_total_cts": 100000, "min_total_cts": 1000, "max_genes": 8000, "min_genes": 250, "pct_counts_mt": 30, "var_genes": 4000},
    "prostate": {"max_total_cts": 50000, "min_total_cts": 1000, "max_genes": 8000, "min_genes": 200, "pct_counts_mt": 30,"var_genes": 2000}
    }
    

    #filter adata object on counts, genes and percentage of mitochondrial genes
    print("Filter out based on thresholds")
    print(f"    Before: {adata.shape[0]}")
    
    sc.pp.filter_cells(adata, min_counts=filter_values[sample_type]["min_total_cts"])
    sc.pp.filter_cells(adata, max_counts=filter_values[sample_type]["max_total_cts"])
    sc.pp.filter_cells(adata, min_genes=filter_values[sample_type]["min_genes"])
    sc.pp.filter_cells(adata, max_genes=filter_values[sample_type]["max_genes"])
    
    adata = adata[adata.obs["pct_counts_mt"] < filter_values[sample_type]["pct_counts_mt"]].copy()
    
    print(f"    After: {adata.shape[0]}")
    
    #plot for QC post filtering
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        groupby="sample_origin",
        size = 0.25,
        jitter=0.2,
        multi_panel=True,
        save=f"_plot_{sample_type}_postqc.png",
        show=False,
        xlabel = sample_type
    )
    
    #normalization and transformation
    adata.layers["counts"] = adata.X #optional - saves raw counts
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    
    #variable gene selection:
    
    sc.pp.highly_variable_genes(
        adata,
        layer = "counts",
        n_top_genes=filter_values[sample_type]["var_genes"],
        subset=True,
        flavor="seurat_v3",
    )
    
    sc.tl.pca(adata, mask_var="highly_variable")
    
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)
    
    sc.pl.umap(adata, color=["leiden"], 
               cmap="viridis", 
               save=f"_plot_{sample_type}.png",
               show=False,
               legend_loc="on data",
               )
    
    sc.pl.umap(adata, color=["pct_counts_mt","sample_origin"], 
               cmap="viridis", 
               save=f"_plot_{sample_type}.png",
               show=False
               )
    
    sc.pl.umap(adata, 
               color="total_counts", 
               vmax=20000, cmap="Reds", 
               size=15, 
               save=f"_plot_{sample_type}_total_counts.png",
               show=False
        )
    
    #write the file
    adata.write(f"{file_base_analysis}/saves/{sample_type}_adata_processed.h5ad")