#!/usr/bin/env python
# coding: utf-8

# The following script reproduces the plots and adata objects from Scheiber et al 2024. 
#We hereby include a manual annotation based on genes regularly associated with certain cell types. For a more detailled cell type annotation please refer to the originally published papers (REFERENZ).

#import dependencies needed for this section:
import os
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import leidenalg as ld
import csv

def cell_annotation(sample_type: str, file_base_analysis: str, setup_dir: str):

    sc.settings.figdir = os.path.join(file_base_analysis, "figures")
    sc.settings.set_figure_params(figsize=(4,4), frameon=False, dpi = 1200)
        
    #load data
    file_path = f'{file_base_analysis}/saves/{sample_type}_adata_processed.h5ad'
    
    adata = sc.read_h5ad(file_path)
    
    # Initialize an empty dictionary
    marker_genes = {}

    genes_csv = setup_dir + "/tables/marker_genes.csv"

    # Read the CSV file to get the marker genes used in this STAR_Protocols
    # Load the CSV file with the correct delimiter
    df = pd.read_csv(genes_csv, delimiter=';')
    
    # Convert the DataFrame to a dictionary with the desired structure
    marker_genes = {}
    for _, row in df.iterrows():
        cell_type = row['Cell Type']
        gene = row['Gene']
        if cell_type not in marker_genes:
            marker_genes[cell_type] = []
        marker_genes[cell_type].append(gene)
    
    #subset list only to markers found in data
    marker_genes_in_data = dict() 
    for ct, markers in marker_genes.items():
        markers_found = list()
        for marker in markers:
            if marker in adata.var.index:
                markers_found.append(marker)
        marker_genes_in_data[ct] = markers_found
    
    cell_list = list(marker_genes.keys())
    
    #plot the genes on DotPlot:
    sc.pl.dotplot(
        adata,
        groupby="leiden",
        var_names=marker_genes_in_data,
        standard_scale="var",  # standard scale: normalize each gene to range from 0 to 1
        show=False,
        save=f"_plot_{sample_type}_marker_genes.png"
    )
    
    #plot the genes un UMAP plots:
    for ct in cell_list:
        print(f"{ct.upper()}:") 
        sc.pl.embedding(
            adata,
            basis="X_umap",
            color=marker_genes_in_data[ct],
            vmin=0,
            vmax="p99",
            sort_order=False,  
            frameon=False,
            size=4,
            cmap="Reds",
            save=f"_plot_{sample_type}_{ct}.png",
            show=False
        )
        print("\n\n\n")  # print white space for legibility
        
    #plot the genes based on highest rank:
    sc.tl.rank_genes_groups(adata, groupby="leiden")
    sc.pl.rank_genes_groups_dotplot(adata, save=f"_plot_{sample_type}_genes_ranked.png", show=False)
    
    #after coarse manual identification assign cluster ID
    ct_map = {
        "lung": {
        "B cells": [13,14,17],
        "Epithelial cells": [1,7,10,12],
        "Endothelial cells": [20],
        "Stromal cells": [16],
        "Mast cells": [15],
        "Neutrophils": [0,4],
        "Monocytes_Macrophages": [3,5,9,18],
        "pDCs": [19],
        "T cells": [2,8,11],
        "NK cells": [6,8,11,21],
        "Progenitors": [],
        },
        "liver": {
        "B cells": [11],
        "Epithelial cells": [10, 13],
        "Endothelial cells": [7],
        "Stromal cells": [],
        "Mast cells": [],
        "Neutrophils": [1,5],
        "Monocytes_Macrophages": [4,12],
        "pDCs": [19],
        "T cells": [3],
        "NK cells": [0,2,6],
        "Progenitors": [5,8]
        },
        "prostate": {
        "B cells": [12],
        "Epithelial cells": [0,2,5,6,9,10],
        "Endothelial cells": [4,13],
        "Stromal cells": [8,14],
        "Mast cells": [7],
        "Neutrophils": [],
        "Monocytes_Macrophages": [3],
        "pDCs": [19],
        "T cells": [1],
        "NK cells": [11],
        "Progenitors": []
        }
    }
    
    cell_type_colors = {
        "B cells": "#1f77b4",  # blue
        "Epithelial cells": "#ff7f0e",  # orange
        "Endothelial cells": "#2ca02c",  # green
        "Stromal cells": "#8B8000",  # yellowish
        "Mast cells": "#9467bd",  # purple
        "Neutrophils": "#8c564b",  # brown
        "Monocytes_Macrophages": "#e377c2",  # pink
        "pDCs": "#7f7f7f",  # gray
        "T cells": "#cc0000",  # reddish
        "NK cells": "#99cc00",  # light green
        "Progenitors": "#17becf", # strong cyan
        "other": "#ffffcc" #light yellow
    }
    
    
    #assigns cell type to adata
    def annotate_cell_types(
            adata,
            cell_type_map,
            *,
            key_added="cell_type",
            default_cell_type="other",
            column="leiden",
        ):
            """Generate a cell-type column from a Mapping cell_type -> [list of clusters]"""
            res = np.full(adata.shape[0], default_cell_type, dtype=object)
            for ct, clusters in cell_type_map.items():
                clusters = [str(x) for x in clusters]
                res[adata.obs[column].isin(clusters)] = ct
            adata.obs[key_added] = res
    
    annotate_cell_types(adata, ct_map[sample_type])
    
    #assigns color to cell type
    adata.obs['color'] = adata.obs['cell_type'].map(cell_type_colors)
    
    #plot final plots:
    sc.pl.umap(
        adata, 
        color="cell_type", 
        palette=cell_type_colors, 
        size=10,
        save=f"_plot_{sample_type}_cell_annotation",
        show=False
    )
    
    #swap to backend
    #plt.switch_backend('Agg')

    
    # Plot each variable with the specified figure size and color scheme
    for obs, figsize in zip(["cell_type"], [(6,5), (7,5)]):
        with plt.rc_context({"figure.figsize": figsize}):
            for var, label in zip(["n_genes_by_counts", "pct_counts_mt","total_counts"], ["n_genes", "% mitochondrial reads","Transcript counts"]):
                fig, ax = plt.subplots()
                sc.pl.violin(adata, var, groupby=obs, rotation=90, ylabel=label, ax=ax, show=False)
                fig.savefig(f"{file_base_analysis}/figures/{sample_type}_{var}_violin_plot.png", bbox_inches="tight")
                plt.close(fig) 
        
    adata.write(f"{file_base_analysis}/saves/{sample_type}_adata_annotated.h5ad")

