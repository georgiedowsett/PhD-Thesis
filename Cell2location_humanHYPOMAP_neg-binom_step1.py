#this file will cread the NucSeq gene models for each cluster level of the neuronal subset of the human_hypomap nucseq data 

import os
import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
os.environ["THEANO_FLAGS"] = 'device=0,floatX=float32,force_device=True'
import cell2location
from cell2location.models import RegressionModel
import scvi
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rcParams
from cell2location.utils import select_slide
from scipy import sparse
from cell2location.utils.filtering import filter_genes
from cell2location.plt import plot_spatial

clusters = ['C1','C2','C3','C4','C5','C6']

for cluster in clusters:
    print(cluster)
    
    results_folder = f"./results/230312/{cluster}/"
    sp_data_folder = './data/st-test/'
    nucseq_folder = './data/'

    # create paths and names to results folders for reference regression and cell2location models
    ref_run_name = f'{results_folder}/reference_signatures'
    run_name = f'{results_folder}/cell2location_map'

    #loading scRNAseq ref data
    # Read data
    adata_ref = anndata.read_h5ad(nucseq_folder + "human_hypo_combined.h5ad")


    #in the trial data they work on the raw count data and not the normalised data. 
    #move raw counts to main datafeame
    adata_ref = adata_ref.raw.to_adata()



    #adding in gene names and ensid
    adata_ref.var = adata_ref.var.rename(columns={"_index": "SYMBOL"})
    adata_ref.var.index = adata_ref.var['SYMBOL']
    #adata_ref.var
    #adata_snrna.var['gene_ids'] = t.var_names
    #adata_snrna.var.set_index('gene_ids', drop=True, inplace=True)
    del adata_ref.raw

    #filtering genes
    selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.08, nonz_mean_cutoff=1.4)
    adata_ref = adata_ref[:, selected].copy()

    # prepare anndata for the regression model
    cell2location.models.RegressionModel.setup_anndata(adata=adata_ref, batch_key='Sample_ID', labels_key=cluster, categorical_covariate_keys=['Donor_ID', 'sex', 'Dataset'])

                            
	
    mod = RegressionModel(adata_ref)
    print(mod.view_anndata_setup())

    mod.train(max_epochs=250, use_gpu=True)

    plt.clf()
    mod.plot_history(20)
    plt.savefig(results_folder + 'elbo-plot-training.png')

    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_ref = mod.export_posterior(adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True})

    # Save model
    mod.save(f"{ref_run_name}", overwrite=True) 
    # Save anndata object with results
    adata_file = f"{ref_run_name}/sc.h5ad"
    adata_ref.write(adata_file)
    adata_file

    #figure out how to save both plots instead of just one
    plt.clf()
    mod.plot_QC()
    plt.savefig(results_folder + 'model-qc-plots.png')
